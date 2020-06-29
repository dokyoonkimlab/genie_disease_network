library(MendelianRandomization)
library(ieugwasr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

#input_dirs = "/cbica/projects/kimlab_hpeace/projects/genie_mendelian_randomization/genie/convert_genie_output/results/out"
input_dirs_iv = args[1]
input_dirs_exp = args[2]
input_dirs_out = args[3]
output_dir = args[4]
PVAL_CUTOFF = as.numeric(args[5])
ADJUST_CORRELATION = as.logical(args[6])
RUN_MAX_LIK = as.logical(args[7])
SGE_INDEX  = as.integer(args[8])
RUN_ONLY_ONE_PAIR_FOR_EACH_SGE_INDEX = as.logical(args[9])
# PVAL_CUTOFF = 10^-4
#ADJUST_CORRELATION = T

dir.create(output_dir)

input_dirs_iv = strsplit(input_dirs_iv, split = ',')[[1]]
input_dirs_exp = strsplit(input_dirs_exp, split = ',')[[1]]
input_dirs_out = strsplit(input_dirs_out, split = ',')[[1]]

phe_pairs = read.csv(file = "/cbica/projects/kimlab_hpeace/projects/genie_mendelian_randomization/genie/cross_phenotype.csv", stringsAsFactors = F)
if (RUN_ONLY_ONE_PAIR_FOR_EACH_SGE_INDEX) {
  phe_pairs = phe_pairs[SGE_INDEX, ]
} else {
  exposure_phe = unique(phe_pairs[,1])[SGE_INDEX]
  phe_pairs = phe_pairs[phe_pairs[,1] == exposure_phe, ]
}


find_file_path = function(x, input_dirs) {
  for (input_dir in input_dirs) {
    f = paste(input_dir, paste("result_", x, sep = ''), "output.txt", sep = .Platform$file.sep)
    if (file.exists(f)) {
      return(f)
    }
  }
  return(NA)
}

file_pairs = phe_pairs
file_pairs[,3] = sapply(file_pairs[,1], find_file_path, input_dirs = input_dirs_iv)
file_pairs[,1] = sapply(file_pairs[,1], find_file_path, input_dirs = input_dirs_exp)
file_pairs[,2] = sapply(file_pairs[,2], find_file_path, input_dirs = input_dirs_out)

phe = NA
if (ADJUST_CORRELATION) {
  discrete_phe = read.table(file = "/cbica/projects/kimlab_hpeace/projects/genie_phewas/QC/imputed/association1/discrete.phe", header = T, stringsAsFactors = F)
  continous_phe = read.table(file = "/cbica/projects/kimlab_hpeace/projects/genie_phewas/QC/imputed/association1/continous.phe", header = T, stringsAsFactors = F)
  phe = merge(discrete_phe, continous_phe)
}

ivs = fread(file = file_pairs[1, 3], sep = '\t')
ivs = ivs[ivs$Var1_Pval < PVAL_CUTOFF,]
d = data.frame(rsid=ivs$Var1_ID, pval=ivs$Var1_Pval, id=ivs$Outcome)
exp_clumped = ld_clump(d, clump_kb=10000, clump_r2=0.001, clump_p=1, 
                      bfile="/cbica/projects/kimlab_hpeace/projects/genie_phewas/QC/imputed/5_hwe/merged_chr1-23_maf0.01_rsq0.7_addrmsnps", 
                      plink_bin="/cbica/projects/kimlab/software_dir/modules/common/GeneralSoftwares/plink")

exposure = fread(file = file_pairs[1, 1], sep = '\t')
exposure = exposure[exposure$Var1_ID %in% exp_clumped$rsid, ]

calc_MR = function(exp_phe, out_phe, exposure, outcome_file) {
  correlation = 0
  if(ADJUST_CORRELATION) {
    correlation = cor(phe[,exp_phe], phe[,out_phe], method = "pearson", use = "pairwise.complete.obs")
    cat(paste("Correlation used", exp_phe, out_phe, correlation, '\n'))
  }
  outcome = fread(file = outcome_file, sep = '\t')
  # outcome = outcome[outcome$Var1_ID %in% exposure$Var1_ID, ]
  merged_exp_out = merge(exposure, outcome, by = "Var1_ID")
  merged_exp_out = merged_exp_out[merged_exp_out$Var1_Pval.y >= 0.01,]
    
  MRInputObject = mr_input(merged_exp_out$Var1_beta.x, bxse = merged_exp_out$Var1_SE.x, by = merged_exp_out$Var1_beta.y, byse = merged_exp_out$Var1_SE.y)

  if(nrow(merged_exp_out) > 2) { 
    ivw = mr_ivw(MRInputObject, psi = correlation)

    result = c(exp_phe, out_phe, "IVW", ivw@Estimate, ivw@StdError, ivw@CILower, ivw@CIUpper, ivw@Pvalue)
    egger = mr_egger(MRInputObject)
    result = rbind(result, c(exp_phe, out_phe, "Egger", egger@Estimate, egger@StdError.Est, egger@CILower.Est, egger@CIUpper.Est, egger@Pvalue.Est))
    median = mr_median(MRInputObject)
    result = rbind(result, c(exp_phe, out_phe, "median", median@Estimate, median@StdError, median@CILower, median@CIUpper, median@Pvalue))
    if (RUN_MAX_LIK) {
      r = tryCatch({
        cat("start mr_maxlik\n")
        flush.console()
        maxlike = mr_maxlik(MRInputObject, psi = correlation)
        cat("done mr_maxlik\n")
        flush.console()
        return(data.frame(exp_phe, out_phe, "MaxLik", maxlike@Estimate, maxlike@StdError, maxlike@CILower, maxlike@CIUpper, maxlike@Pvalue))
      }, error = function(e) {
        cat(e)
        cat("\n")
        return(data.frame(exp_phe, out_phe, "MaxLik", -1, -1, -1, -1, -1))
      })
      result = rbind(result, r)
    }
    result = data.frame(result, stringsAsFactors = F)
    colnames(result) = c("Exposure", "Outcome", "Method", "Estimate", "StdError", "CILower", "CIUpper", "Pvalue")
    return(result)
  }
  return(data.frame())
}

res = -1
for (i in 1:nrow(file_pairs)) {
  calc_res = calc_MR(phe_pairs[i, 1], phe_pairs[i, 2], exposure, file_pairs[i, 2])
  if (!is.data.frame(res)) {
    res = calc_res
  } else {
    res = rbind(res, calc_res)
  }
}

if (SGE_INDEX == 1) {
  colnames(res) = c("Exposure", "Outcome", "Method", "Estimate", "StdError", "CILower", "CIUpper", "Pvalue")
  write.table(res, file = paste(output_dir, SGE_INDEX, sep = .Platform$file.sep), row.names = F, quote = F, sep = '\t')
} else {
  write.table(res, file = paste(output_dir, SGE_INDEX, sep = .Platform$file.sep), row.names = F, col.names = F, quote = F, sep = '\t')
}
