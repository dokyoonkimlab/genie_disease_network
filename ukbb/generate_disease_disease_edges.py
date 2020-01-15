import sys
import os
from pathlib import Path

input_ukbb_gwas_results_dir = "../../genie_download_ukbb_jap_results/ukbb_results"
input_ukbb_manifest_file = "../../genie_download_ukbb_jap_results/manifest_files/UKBB GWAS Imputed v3 - File Manifest Release 20180731 - Manifest 201807.txt"
out_file = "ukbb_snpid_disease_disease.csv"
out_file_list = "ukbb_disease_disease_list.csv"

PVAL_CUTOFF = 0.0001

# Read ukbb manifest file and create a map with Phenotype desc and filename on disk
# The result files were downloaded using the ID number as dir name
ukbb_phe_to_file_map = {}
with open(input_ukbb_manifest_file) as inf:
	for i in range(0, 26):
		next(inf)
	for line in inf:
		line = line.strip()
		row = line.split('\t')
		row[0] = row[0].strip()
		row[1] = row[1].strip('\"')
		row[1] = row[1].strip()
		if row[3] == "both_sexes" and not row[0].endswith("irnt"):
			in_dir = os.path.join(input_ukbb_gwas_results_dir, row[0])
			if os.path.exists(in_dir):
				input_files = []
				for in_f in Path(in_dir).rglob('*.tsv'):
					input_files.append(in_f)
				if len(input_files) > 1:
					print("ERROR: Multiple files found for " + row[1] + " in directory " + in_dir + ". Cannot resolve, delete unnecessary files from directory")
					sys.exit(1)
				ukbb_phe_to_file_map[row[1]] = input_files[0]

snpid_phe_map = {}
for phe in ukbb_phe_to_file_map:
	input_file = ukbb_phe_to_file_map[phe]
	print("Processing:" + str(input_file))
	with open(input_file) as inf:
		header = next(inf)
		header = header.strip()
		header = header.split()
		snpid_index = header.index("variant")
		pval_index = header.index("pval")
		low_conf_index = header.index("low_confidence_variant")
		for line in inf:
			line = line.strip()
			row = line.split()
			if row[low_conf_index] == "false" and float(row[pval_index]) < PVAL_CUTOFF:
				if row[snpid_index] not in snpid_phe_map:
					# set is not required, just in case if duplicate lines are present!
					snpid_phe_map[row[snpid_index]] = set()

				snpid_phe_map[row[snpid_index]].add(phe)

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
