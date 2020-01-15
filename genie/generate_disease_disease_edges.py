import sys
import os

input_file = "../../genie_compare_ukbb_jap/data/imputedNARD_results_pval10-4.tsv"
# input_file = "/cbica/projects/kimlab_hpeace/projects/genie_phewas/QC/imputed/association/plato/results/imputedNARD_continous_results_pval10-4.csv"
out_file = "genie_snpid_disease_disease.csv"
out_file_list = "genie_disease_disease_list.csv"

snpid_phe_map = {}
with open(input_file) as inf:
	header = next(inf)
	header = header.strip()
	header = header.split(',')
	phe_index = header.index("Outcome")
	snp_id = header.index("Var1_ID")
	for line in inf:
		line = line.strip()
		row = line.split(',')
		if row[1] not in snpid_phe_map:
			# set is not required, just in case if duplicate lines are present!
			snpid_phe_map[row[1]] = set()

		snpid_phe_map[row[1]].add(row[0])

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
