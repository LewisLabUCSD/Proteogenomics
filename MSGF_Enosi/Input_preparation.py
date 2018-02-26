'''
Created on Dec 2, 2015

@author: s3cha
'''

import sys
import os


'''sys.argv: download folder, gff folder, known db folder, variant folder, splice folder'''
print sys.argv
print os.listdir('./')
file_list = os.listdir(sys.argv[1])
print file_list

for file in file_list:
    folder = sys.argv[1]
    if file.find('.gff')>-1:
        if not os.listdir(sys.argv[2]):
            os.system('mv '+folder+'/'+file+' '+sys.argv[2])
    elif (file.find('GRC')>-1 or file.find('Known')>-1) and file.find('.gff')<0:
        if not os.listdir(sys.argv[3]):
            os.system('mv '+folder+'/'+file+' '+sys.argv[3])
    elif file.find('splice')>-1 or file.find('Splice')>-1:
        os.system('mv '+folder+'/'+file+' '+sys.argv[4])
    elif file.find('VGraph')>-1 or file.find('Vgraph')>-1:
        os.system('mv '+folder+'/'+file+' '+sys.argv[5])

# python /home/s3cha/data/ProteoSaFe/ProteoSAFe-1.2.6_beta-linux64/tools/Enosi/2014.1208/Input_preparation.py dbset fasta1 sequence1 sequence2 sequence3