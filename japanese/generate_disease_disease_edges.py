import sys
import os
from pathlib import Path

input_jap_gwas_results_dir = "../../genie_download_ukbb_jap_results/jap_results"
input_jap_manifest_file = "../../genie_download_ukbb_jap_results/manifest_files/jap_riken_manifest.txt"
out_file = "jap_snpid_disease_disease.csv"
out_file_list = "jap_disease_disease_list.csv"

PVAL_CUTOFF = 0.0001

# Read japanese manifest file and create a map with Phenotype desc and filename on disk
# The result files were downloaded using the ID number as dir name
jap_phe_to_file_map = {}
with open(input_jap_manifest_file) as inf:
	for line in inf:
		line = line.strip()
		row = line.split('\t')
		row[0] = row[0].strip()
		row[1] = row[1].strip('\"')
		row[1] = row[1].strip()
		in_dir = os.path.join(input_jap_gwas_results_dir, row[0])
		if row[8] == "Yes":
			in_dir = os.path.join(input_jap_gwas_results_dir, "QTL_" + row[0])
		input_files = Path(in_dir).rglob('*.txt')
		for input_file in input_files:
			filename = os.path.basename(input_file)
			if filename != "README.txt" and "chrX" not in filename and "chrx" not in filename:
				if row[1] in jap_phe_to_file_map:
					print("ERROR: Multiple files found for " + row[1] + " in directory " + in_dir + ". Cannot resolve, delete unnecessary files from directory")
					sys.exit(1)
				jap_phe_to_file_map[row[1]] = input_file


pval_colnames = ["P","p.value","PVALUE","P_BOLT"]
snpid_phe_map = {}
for phe in jap_phe_to_file_map:
	input_file = jap_phe_to_file_map[phe]
	print("Processing:" + str(input_file))
	with open(input_file) as inf:
		header = next(inf)
		header = header.strip()
		header = header.split()

		chr_index = header.index("CHR")
		pos_index = None
		if "POS" in header:
			pos_index = header.index("POS")
		else:
			pos_index = header.index("BP")
		pval_index = None
		for colname in pval_colnames:
			if colname in header:
				pval_index = header.index(colname)
				break

		for line in inf:
			line = line.strip()
			row = line.split()

			if float(row[pval_index]) < PVAL_CUTOFF:
				chrpos = row[chr_index] + ':' + row[pos_index]
				if chrpos not in snpid_phe_map:
					# set is not required, just in case if duplicate lines are present!
					snpid_phe_map[chrpos] = set()

				snpid_phe_map[chrpos].add(phe)

with open(out_file, 'w') as of:
	for snpid in snpid_phe_map:
		phe_list = list(snpid_phe_map[snpid])
		for i in range(0, len(phe_list)):
			for j in range(i + 1, len(phe_list)):
				of.write(snpid + ',' + phe_list[i] + ',' + phe_list[j] + '\n')


phe_phes_map = {}
for snpid in snpid_phe_map:
	phe_set = snpid_phe_map[snpid]
	for phe in phe_set:
		if phe not in phe_phes_map:
			phe_phes_map[phe] = set()
		phe_phes_map[phe] = phe_phes_map[phe].union(phe_set)

max_len = None
data_arr = []
for phe, val in sorted(phe_phes_map.items(), key=lambda item: len(item[1]), reverse = True):
	phe_phes_map[phe].remove(phe)
	if not max_len:
		max_len = len(phe_phes_map[phe]) + 2
	data_list = [phe, str(len(phe_phes_map[phe]))]
	data_list.extend(list(phe_phes_map[phe]))

	data_arr.append(data_list)

with open(out_file_list, 'w') as of:
	for i in range(0, max_len):
		for j in range(0, len(data_arr)):
			if (i < len(data_arr[j])):
				of.write(data_arr[j][i] + ',')
			else:
				of.write(',')
		of.write('\n')
