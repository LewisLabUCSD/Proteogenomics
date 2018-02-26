import os
import sys
import math
import re
import pickle

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r

def inserts(original, new, pos):
	return original[:pos] + new + original[pos:]


def	CleanPeptideString(pep_str):
	#print pep_str
	pep_str = pep_str.replace("0","")
	pep_str = pep_str.replace("1","")
	pep_str = pep_str.replace("2","")
	pep_str = pep_str.replace("3","")
	pep_str = pep_str.replace("4","")
	pep_str = pep_str.replace("5","")
	pep_str = pep_str.replace("6","")
	pep_str = pep_str.replace("7","")
	pep_str = pep_str.replace("8","")
	pep_str = pep_str.replace("9","")
	pep_str = pep_str.replace("+","")
	pep_str = pep_str.replace(".","")
	pep_str = pep_str.replace("?","_")
	pep_str = pep_str.replace("_","")
	#pep_str = inserts(pep_str,".",1)
	#pep_str = inserts(pep_str,".",-1)
	#print pep_str
	return pep_str


input_pickle_file_name = ""

RefSeq_file_name = sys.argv[1]
PSM_file_name = sys.argv[2]
Filtered_out_file_name = sys.argv[3]

decoy_str = "XXX"


file_name_col = 0	##SpecFile	
spec_id_col = 1	#SpecID	
index_col = 2	#ScanNum	
frag_meth_col = 3	#FragMethod	
precursor_col = 4	#Precursor	
iso_error_col = 5	#IsotopeError	
pre_error_col = 6	#PrecursorError(ppm)	
charge_col = 7	#Charge	
pep_seq_col = 8	#Peptide	
protein_db_col = 9	#Protein	
denov_score_col = 10	#DeNovoScore	
msgf_score_col = 11	#MSGFScore	
spec_e_val_score_col = 12	#SpecEValue	
e_value_col = 13	#EValue

#assign score
score_col = spec_e_val_score_col
score_higher_better = 0

SourceFile = open(RefSeq_file_name, "r")
PSMSourceFile = open(PSM_file_name, "r")
FilteredOutFile = open(Filtered_out_file_name, "w")

combined_tvs_dic = {}
for line in SourceFile:
	if line == "":
		continue
	if line.startswith("#"):
		caption = line
		continue
	data = re.split('\t',line)
	if len(data) <= score_col:
		continue
	try:
	    float(data[score_col])
	except:
	    continue
	score = float(data[score_col])
	file_name_ind = data[file_name_col]
	spec_ind = data[index_col]
	key_ind = file_name_ind + spec_ind
	#print key_ind
	if combined_tvs_dic.has_key(key_ind):
		print key_ind, combined_tvs_dic[key_ind], "already in the file"
	else:
		combined_tvs_dic[key_ind] = 1


for line in PSMSourceFile:
	if line == "":
		continue
	if line.startswith("#"):
		FilteredOutFile.write(line)
		continue
	data = re.split('\t',line)
	if len(data) <= score_col:
		continue
	try:
	    float(data[score_col])
	except:
	    continue
	score = float(data[score_col])
	file_name_ind = data[file_name_col]
	spec_ind = data[index_col]
	key_ind = file_name_ind + spec_ind
	if combined_tvs_dic.has_key(key_ind):
		continue
	else:
		FilteredOutFile.write(line)


SourceFile.close()
PSMSourceFile.close()
FilteredOutFile.close()






