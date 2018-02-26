'''
Created on Mar 3, 2015

@author: s3cha
'''

import os
import sys

parameter_file = sys.argv[1]
spectrum_folder = sys.argv[2]
tsv_file = sys.argv[3]
output_file_name = sys.argv[4]

def CleanPepSeq(peptide):
    peptide = peptide.replace('+','')
    peptide = peptide.replace('0','')
    peptide = peptide.replace('1','')
    peptide = peptide.replace('2','')
    peptide = peptide.replace('3','')
    peptide = peptide.replace('4','')
    peptide = peptide.replace('5','')
    peptide = peptide.replace('6','')
    peptide = peptide.replace('7','')
    peptide = peptide.replace('8','')
    peptide = peptide.replace('9','')
    peptide = peptide.replace(':','')
    peptide = peptide.replace('[','')
    peptide = peptide.replace(']','')
    peptide = peptide.replace('(','')
    peptide = peptide.replace(')','')
    peptide = peptide.replace(',','')
    if peptide[1] == '.':
        peptide = peptide[2:-2]
    peptide = peptide.replace('.','')
    peptide = peptide.replace('-','')
    return peptide


f = open(parameter_file,'r')

file_map_dictionary = {}
for line in f:
    if line.find('upload_file_mapping')>-1:
        data = line.split('>')[1].split('<')[0].split('|')
        file_name = os.path.split(data[1])[1]
        if not file_map_dictionary.has_key(file_name):
            file_map_dictionary[file_name] = data[0]
        else:
            print 'Same file name exist multiple times: ',file_name
f.close()

tsv_file_iterator = open(tsv_file,'r')
s = open(output_file_name,'w')

SpecFile_column = 0
Peptide_column = 0
for line in tsv_file_iterator:
    if line.startswith('#'):
        line = line.replace("internalFilename","nameFileInternal")
        data = line.strip().split('\t')
        for index,item in enumerate(data):
            if item.startswith('#SpecFile') or item.startswith('SpecFile') or item.startswith('SpectrumFile'): #add any type of spectrum file header identifier here
                SpecFile_column = index
            elif item.startswith('Peptide'):
                Peptide_column = index
        s.write(line.strip()+'\tpep\n')
        continue
    data = line.strip().split('\t')
    if file_map_dictionary.has_key(data[0]):
        filename = file_map_dictionary.get(data[0])
    else:
        filename = data[0]
    pep = data[Peptide_column]
    data.pop(0)
    temp = '\t'.join(data)
    s.write(filename+'\t'+temp+'\t'+CleanPepSeq(pep)+'\n')

print 'END:'
