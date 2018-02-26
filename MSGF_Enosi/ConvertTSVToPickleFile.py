import os
import sys
import math
import re
import cPickle as pickle

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r

def inserts(original, new, pos):
	return original[:pos] + new + original[pos:]


def	CleanPeptideString(pep_str):
	#print pep_str
    if pep_str[1] == '.':
        pep_str = pep_str[2:-2]
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
	pep_str = pep_str.replace("-","")
	pep_str = pep_str.replace("*","")
    pep_str = pep_str.replace("(","")
    pep_str = pep_str.replace(")","")
    pep_str = pep_str.replace("[","")
    pep_str = pep_str.replace("]","")
    pep_str = pep_str.replace(",","")
    return pep_str


input_pickle_file_name = ""

file_name = sys.argv[1]
output_pickle_file_name = sys.argv[2]


decoy_str = "XXX"

# Combine MSGFPlus tvs results
print "\nCombine MSGFPlus tvs results"

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


score_col = spec_e_val_score_col

count_lines = 0
count_peptides = 0
pepdic = {}
SourceFile = open(file_name, "r")
for line in SourceFile:
	if line == "":
		continue
	if line.startswith("#"):
		caption = line
        data = line.strip().split('\t')
        for ind,item in enumerate(data):
            if item == "Peptide" or item == "peptide":
                pep_seq_col = ind
		continue
	data = re.split('\t',line)
	if len(data) < 12:
		continue
	try:
	    float(data[spec_e_val_score_col])
	except:
	    continue
	count_lines += 1

	clean_pep_seq = CleanPeptideString(data[pep_seq_col])

	if pepdic.has_key(clean_pep_seq):
		pepdic[clean_pep_seq][0].append(line)
		pepdic[clean_pep_seq][3] = pepdic[clean_pep_seq][3] + 1
		if float(data[15]) < pepdic[clean_pep_seq][4]:
			pepdic[clean_pep_seq][4] = float(data[15])
	else:
		pepdic[clean_pep_seq] = [[],{},[],1,float(data[15]),0] # [[spectra file information],[location information],[extra needed],peptides count,fdr,location count] 
		pepdic[clean_pep_seq][0].append(line)
		count_peptides += 1

SourceFile.close()

pickle.dump(pepdic,open(output_pickle_file_name,"wb"))

print "Count Peptides:", count_peptides
print "Count  PSMs   :", count_lines












