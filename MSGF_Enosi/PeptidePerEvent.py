'''
Created on Mar 3, 2015

@author: s3cha
'''

import sys
import os

input_event_file = sys.argv[1]
output_file = sys.argv[2]

f = open(input_event_file,'r')
s = open(output_file,'w')
s.write('#Num\tEvent\tPeptide\tLocation\n')
event_column = 0
pep_column = 0
group_column = 0
location_column = 0 
def CleanPepSeq(string):
    string = string.replace(':','')
    if string[1] == '.':
        string = string[2:-2]
    return string

for line in f:
    if line.startswith('#'):
        line = line.split('\t')
        for index,item in enumerate(line):
            if item == 'Event':
                event_column = index
            elif item.startswith('Peptide'):
                pep_column = index
            elif item.startswith('GroupInfo'):
                group_column = index
            elif item.startswith('Location'):
                location_column = index 
        continue
    line = line.strip().split('\t')
    peptide = []
    peptide.append([CleanPepSeq(line[pep_column]),line[location_column]])
    if line[group_column] != '':
        group_info = line[group_column].strip('|').split('|')
        for group_part in group_info:
            pep = CleanPepSeq(group_part.split('/')[0])
            loc = group_part.split('/')[1]
            peptide.append([pep,loc])
    for pep in peptide:
        s.write(line[0]+'\t'+line[1]+'\t'+pep[0]+'\t'+pep[1]+'\n')    

print 'END:'