import sys
import os
import csv
from pathlib import Path

input_match_manifest_file = "../../genie_compare_ukbb_jap/manifest_files/Matching_phe_genie_ukbb_jap.txt"
input_genie_ddn_file = "../genie/genie_disease_disease.csv"
input_jap_ddn_file = "../japanese/jap_disease_disease_uniq.csv"
input_ukbb_ddn_file = "../ukbb/ukbb_disease_disease_uniq.csv"

PRINT_WARNING = False

class DDConnection:

	def __init__(self, d1, d2, convert_map = None):
		if convert_map:
			self.d1 = convert_map[d1]
			self.d2 = convert_map[d2]
		else:
			self.d1 = d1
			self.d2 = d2

	@staticmethod
	def create_ddcon_instances(d1, d2, convert_map):
		if not convert_map:
			return [DDConnection(d1, d2)]

		ddcon_objs = []
		for od1 in convert_map[d1]:
			for od2 in convert_map[d2]:
				ddcon_objs.append(DDConnection(od1, od2))
		return ddcon_objs

	def convert_disease_names(self, convert_map):
		self.d1 = convert_map[self.d1]
		self.d2 = convert_map[self.d2]

	def is_self_loop(self):
		return self.d1 == self.d2

	def get_disease_set(self):
		return(frozenset([self.d1, self.d2]))

	def __hash__(self):
		return hash(self.get_disease_set())

	def __eq__(self, other):
		return self.get_disease_set() == other.get_disease_set()

	def __repr__(self):
		return self.d1 + "<->" + self.d2

	def __str__(self):
		return self.d1 + "<->" + self.d2


def parse_ddns(input_file, phe_restrict_set = None, phe_convert_map = None):
	ddn_set = set()
	with open(input_file) as inf:
		for line in inf:
			line = line.strip()
			row = line.split(',')
			row = sorted(row)
			if phe_restrict_set:
				if row[0] in phe_restrict_set and row[1] in phe_restrict_set:
					ddc_objs = DDConnection.create_ddcon_instances(row[0], row[1], phe_convert_map)
					for ddc in ddc_objs:
						if not ddc.is_self_loop():
							ddn_set.add(ddc)
				else:
					if PRINT_WARNING:
						print("WARNING: removed pair:" + line)
			else:
				ddc_objs = DDConnection.create_ddcon_instances(row[0], row[1], phe_convert_map)
				for ddc in ddc_objs:
					if not ddc.is_self_loop():
						ddn_set.add(ddc)
	return ddn_set

def get_values_from_map(phe_map):
	all_set = set()
	for phe in phe_map:
		all_set.update(phe_map[phe])
	return all_set

# Read the matching jap phenptypes for genie phenotypes manually curated to a map
jap_genie_phe_map = {}
with open(input_match_manifest_file) as inf:
	reader = csv.reader(inf, delimiter = '\t')
	next(reader)
	for parts in reader:
		missing = True
		for i in range(3,8):
			if parts[i]:
				missing = False
		if not missing:
			for i in range(1,3):
				parts[i] = parts[i].strip('\"')
				parts[i] = parts[i].strip()
				if parts[i] and "chrX" not in parts[i]:
					if parts[i] not in jap_genie_phe_map:
						jap_genie_phe_map[parts[i]] = set()	
					jap_genie_phe_map[parts[i]].add(parts[0])

# Read the matching ukbb phenptypes for genie phenotypes manually curated to a map
ukbb_genie_phe_map = {}
with open(input_match_manifest_file) as inf:
	reader = csv.reader(inf, delimiter = '\t')
	next(reader)
	for parts in reader:
		missing = True
		for i in range(1,3):
			if parts[i] and "chrX" not in parts[i]:
				missing = False
		if not missing:
			for i in range(3,8):
				parts[i] = parts[i].strip()
				if parts[i]:
					if parts[i] not in ukbb_genie_phe_map:
						ukbb_genie_phe_map[parts[i]] = set()
					ukbb_genie_phe_map[parts[i]].add(parts[0])

print(len(get_values_from_map(jap_genie_phe_map)))
genie_ddn = parse_ddns(input_genie_ddn_file, get_values_from_map(jap_genie_phe_map))
jap_ddn = parse_ddns(input_jap_ddn_file,  jap_genie_phe_map.keys(), jap_genie_phe_map)
ukbb_ddn = parse_ddns(input_ukbb_ddn_file, ukbb_genie_phe_map.keys(), ukbb_genie_phe_map)

# Compare jap vs genie
print(str(len(genie_ddn)) + " " +  str(len(jap_ddn)) + " " + str(len(ukbb_ddn)))
print(str(len(genie_ddn & jap_ddn)) + " " + str(len(genie_ddn & ukbb_ddn)) + " " + str(len(jap_ddn & ukbb_ddn)))

print(str(len(genie_ddn & jap_ddn & ukbb_ddn)))

print(str(len(genie_ddn - jap_ddn)))
print(len(genie_ddn - (jap_ddn | ukbb_ddn)))
print(genie_ddn & jap_ddn & ukbb_ddn)

#genie_ddn1 = parse_ddns(input_genie_ddn_file, get_values_from_map(ukbb_genie_phe_map), genie_ukbb_phe_convert_map)

#print(str(len(genie_ddn1)) + " " +  str(len(ukbb_ddn)) + " " + str(len(genie_ddn1 & ukbb_ddn)) + " " + str(len(genie_ddn1 - ukbb_ddn)))

#print(genie_ddn - (jap_ddn | ukbb_ddn))


