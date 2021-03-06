'''
Created on Dec 7, 2015

@author: s3cha
Input: Event ( or any other type of tap seperated file including peptide ), p-file ( including PSM information ), spectrum
Output: Event to spectrum match on spectrum viewer.

'''

import os
import sys
import cPickle as pickle

if len(sys.argv)>=2:
    input_event_filename = sys.argv[1]
    event_folder = sys.argv[2]
    p_filename = sys.argv[3]
    parameter_file = sys.argv[4]
    spectrum_folder = sys.argv[5]
else:
    input_event_filename = ''
    event_folder = ''
    p_filename = ''
    parameter_file = ''
    spectrum_folder = ''

os.system('mkdir '+event_folder)
os.system('mv '+input_event_filename+' '+event_folder)    
event_filename = event_folder +'/'+ os.path.split(input_event_filename)[1]

event_file = open(event_filename,'r')
p_folder = os.listdir(p_filename)
pepdic = {}
input_folder = os.path.split(sys.argv[0])
system_folder = input_folder[0]
if input_folder[0] == '':
    system_folder = '.'

def CleanPepSeq(pep):
    return pep.replace(':','').replace('.','').replace('-','').replace('?','')

def Replace_mod(pep_string):
    while pep_string.find('+')>-1:
        index = pep_string.find('+')
        pep_string = '{0}({1},{2}.{3}){4}'.format(pep_string[:index-1],pep_string[index-1],pep_string[index+1:].split('.')[0],pep_string[index+1:].split('.')[1][:3],'.'.join([pep_string[index+1:].split('.')[1][3:],'.'.join(pep_string[index+1:].split('.')[2:])]))
    return pep_string

# Replace_mod('K.STS+125.312GGTAALGC+57.021LVK+144.132.D')

pep_column = None
peptide = []
for line in event_file:
    if line.startswith('#'):
        line = line.strip().split('\t')
        for index,header in enumerate(line):
            if header == 'Peptide':
                pep_column = index
                break
        continue
    line = line.strip().split('\t')
    peptide.append(CleanPepSeq(line[pep_column]))

for file in p_folder:
    temp_pepdic = pickle.load(open(p_filename+'/'+file,'rb'))
    for key in temp_pepdic:
#         if key not in peptide:
#             continue
        if key not in pepdic:
            pepdic[key] = temp_pepdic[key]
        else:
            for temp in temp_pepdic[key][0]:
                pepdic[key][0].append(temp)

s = open('./temp.txt','w')
#s.write('#SpecFile\tSpecID\tScanNum\tFragMethod\tPrecursor\tIsotopeError\tPrecursorError(ppm)\tCharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\tSpecEValue\tEValue\tQValue\tPepQValue\tFDR\tPepFDR\n')
header = '#SpecFile\tSpecID\tScanNum\tFragMethod\tPrecursor\tIsotopeError\tPrecursorError(ppm)\tCharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\tSpecEValue\tEValue\tQValue\tPepQValue\tFDR\tPepFDR'.split('\t')

mismatch = 1

error_check = True
for key in pepdic:
    for item in pepdic[key][0]:
        temp_item = item.split('\t')
        if error_check:
            error_check= False
            size_vec = len(temp_item)
	    s.write('\t'.join(header[:size_vec])+'\n')
        if len(temp_item) != size_vec:
            continue
        if temp_item[0].split('.')[1].lower() == 'mgf':
            if temp_item[1].startswith('index='):
                temp_item[1] = 'index='+str(int(temp_item[1].replace('index=',''))+mismatch)
            elif temp_item[1].startswith('scan='):
                temp_item[1] = 'scan='+str(int(temp_item[1].replace('scan=',''))+mismatch)
            temp_item[2] = str(int(temp_item[2])+mismatch)
        
#         temp = temp_item[8].split('.')
#         temp_item[8] = '.'.join(temp)
        temp_item[8] = Replace_mod(temp_item[8])
        item = '\t'.join(temp_item)
        s.write(item)
s.close()


if os.path.isfile(system_folder+'/PeptidePerEvent.py'):
    os.system('python '+system_folder+'/PeptidePerEvent.py '+event_filename+' '+os.path.splitext(event_filename)[0]+'_pepinfo.txt')
if parameter_file != '' and spectrum_folder != '' and os.path.isfile(system_folder+'/MingleFilename.py'):
    os.system('python '+system_folder+'/MingleFilename.py '+parameter_file+' '+spectrum_folder+' temp.txt '+os.path.splitext(event_filename)[0]+'_specinfo.txt')
os.system('rm temp.txt')








