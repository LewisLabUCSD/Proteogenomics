'''
Created on Oct 29, 2014

@author: s3cha
'''

import pdb
import commands
import sys
import time
import os
from natsort import natsorted

betime = time.clock()
print 'BEGIN: ',time.ctime()

even_check = 0
identifier = ''
tsv_result = []
database = []
ref_seq = ''
dna_folder = ''
output_file_name = ''
parameter_file = ''
spectrum_folder = ''
input_folder = os.path.split(sys.argv[0])
for argument in sys.argv:
    even_check += 1   
    if even_check % 2 == 0:
        identifier = argument
    else:
        if identifier   == '-i': # -i main input of enosi, output of PSM
            tsv_result.append(argument)
        elif identifier == '-db': # -db database. could have multiple db, could be form of directory 
            database.append(argument)
        elif identifier == '-r': # -r reference. should be form of gff
            ref_seq = argument
        elif identifier == '-a': # -a dna folder. optional 
            dna_folder = argument
        elif identifier == '-o':
            output_file_name = argument
        elif identifier == '-s':
            spectrum_folder = argument
        elif identifier == '-param':
            parameter_file = argument
        elif even_check != 1 or identifier == '-h':
            print '***Enosi version 1.0 / input options***\n-h\thelp\n-i\toutput of PSM(if multiple files, use \'-i\' for every file name. could be directory)\n-db\tdatabase. could have multiple db, could be directory(only .fa or .fasta will be readed)',
            print '\n-r\treference. should be form of .gff\n-a\tdna folder. optional'
            sys.exit()
    if identifier == '-h':
        print '***Enosi version 1.0 / input options***\n-h\thelp\n-i\toutput of PSM(if multiple files, use \'-i\' for every file name. could be directory)\n-db\tdatabase. could have multiple db, could be directory(only .fa or .fasta will be readed)',
        print '\n-r\treference. should be form of .gff\n-a\tdna folder. optional'
        sys.exit()

os.system('mkdir '+output_file_name)
s = open(output_file_name+'/log2.ini','w')
s.write('Code begins\n')        
output_file_name = output_file_name + '/event.txt'
output_folder = os.path.split(output_file_name)
output_name = os.path.splitext(output_file_name)
string_question = ['Enter the PSM result file(if multiple files, seperate by space. could be directory), Q/q to exit: '
                   ,'Enter the Database file(if multiple files, seperate by space. could be directory), Q/q to exit: '
                   ,'Enter the Reference genome file(.gff), Q/q to exit: ']

s.write('argument check\n')#
def CheckReadableFileFormat(user_input,index,extension):
    input_file_list = []
    my_list = user_input
    while(1):
        if user_input == []:
            my_list = raw_input(string_question[index])
        else:
            if my_list.__class__ == [].__class__:
                my_list = ' '.join(my_list)
        for file in my_list.strip().split(' '):
            if os.path.isdir(file) or os.path.isdir('./'+file):
                file_list = os.listdir(file)
                directory = file.strip('/').strip('\\')
                for current_file in file_list:
                    current_file = './'+directory+'/'+current_file
                    if not os.path.isfile(current_file):
                        continue
                    elif os.path.splitext(current_file)[1] in extension:
                        input_file_list.append(current_file.replace('\\','/'))
            elif os.path.isfile(file):
                ext = os.path.splitext(file)[1]
                if ext in extension:
                    input_file_list.append(file)
                else:
                    print 'Unreadable extensiton: ',ext
        if input_file_list != []:
            input_file_list = natsorted(list(set(input_file_list)))
            break
        elif my_list == 'q' or my_list == 'Q':
            sys.exit()
        else:
            print 'Error: unreadable format of input or file not exist: ',my_list
            my_list = raw_input(string_question[index])            
    return input_file_list

if tsv_result == ['Result']:
    input_file_list = ['Result/Novel.txt']
else:
    print tsv_result
    input_file_list = CheckReadableFileFormat(tsv_result,0,['.txt','.tsv'])
print '[Input files] '
for i in input_file_list:
    print i

print '[Databases] '
db_list = CheckReadableFileFormat(database,1,['.fa','.fasta'])
for i in db_list:
    print i

print '[Reference file] '
ref_seq_file = CheckReadableFileFormat(ref_seq,2,['.gff3'])[0]
print ref_seq_file

if dna_folder != '':
    print '[DNA]'
    if not os.path.isdir(dna_folder):
        print 'Unable to read DNA folder: ',dna_folder 
        dna_folder = ''

input_file_string = ' '.join(input_file_list)
temp_filename = output_folder[0]+'/temp.txt'
os.system('mkdir '+output_folder[0])
os.system('cat '+input_file_string+' > '+temp_filename)
s.write('P file generation\n')

system_folder = input_folder[0]
if input_folder[0] == '':
    system_folder = '.'
print '***Generating P file***'
os.system('python '+system_folder+'/ConvertTSVToPickleFile.py '+temp_filename+' '+output_name[0]+'_pepdic.p')
print
s.write('Location file generation\n')
print '***Generating Location file***'
print 'python '+system_folder+'/Location05242013.py '+','.join(db_list)+' '+output_name[0]+'_location.txt '+output_name[0]+'_pepdic.p'
os.system('python '+system_folder+'/Location05242013.py '+','.join(db_list)+' '+output_name[0]+'_location.txt '+output_name[0]+'_pepdic.p')

print
s.write('Event file generation\n')
print '***Generating Event file***'
print 'python '+system_folder+'/newEventCaller3.py '+output_name[0]+'_location.txt '+system_folder+'/Location.txt '+ref_seq_file+' '+output_file_name+' '+dna_folder
os.system('python '+system_folder+'/newEventCaller3.py '+output_name[0]+'_location.txt '+system_folder+'/Location.txt '+ref_seq_file+' '+output_file_name+' '+dna_folder)

if os.path.isfile(system_folder+'/PeptidePerEvent.py'):
    os.system('python '+system_folder+'/PeptidePerEvent.py '+output_file_name+' '+os.path.splitext(output_file_name)[0]+'_pepinfo.txt')
if parameter_file != '' and spectrum_folder != '' and os.path.isfile(system_folder+'/MingleFilename.py'):
    os.system('python '+system_folder+'/MingleFilename.py '+parameter_file+' '+spectrum_folder+' '+temp_filename+' '+os.path.splitext(output_file_name)[0]+'_specinfo.txt')
os.system('rm '+output_folder[0]+'/temp.txt')

if len(sys.argv) > 1:
    unknown_location_filename = sys.argv[1]
    known_location_filename = sys.argv[2]
    ref_seq_filename = sys.argv[3]
    output_filename = sys.argv[4]
    if len(sys.argv)>5:
        DNA_folder = sys.argv[5]
else:
    #VU run
    unknown_location_filename = '/home/s3cha/data/VU/Cosmic/Location_VU_SPECFDR01_Cosmicoverlap2.txt'#'/home/s3cha/data/VU/Cosmic/Location_VU_SPECFDR01_Cosmicoverlap2.txt'#'/home/s3cha/data/VU/Somatic/Location_VU_RefSeq_PEPFDR01_rna.txt'#'/home/s3cha/data/PNNL/ov_somatic_mutations/Locations_PNNL_20140202_removed_dbSNP_somatic.txt'#'/home/s3cha/data/VU/Somatic/Location_VU_SPECFDR01_splice_vcf_dbSNP_somatic.txt'#'/home/s3cha/data/JHU/Location_JHU_SV_Result_SVrna_dbSNP_somatic.txt'#'/home/s3cha/data/PNNL/Location_PNNL_Single_recalculated.txt'#'C:/Users/Akard3/Dropbox/Workspace/SpliceGraph/src/Enosi/Locations_PNNL_20140202_removed_dbSNP.txt'#Location_PNNL_Single_recalculated.txt'#sample.txt'#
    known_location_filename = 'Location_Combine123_RefSeq_recalculated_Correction.txt'#'/home/s3cha/data/VU/Location_VU_RefSeq_PEPFDR01.txt'#
    ref_seq_filename = 'Homo_sapiens.GRCh37.70_formatted.gff'
    output_filename = '/home/s3cha/data/VU/Cosmic/Event3_VU_somatic_cosmic.txt'#'/home/s3cha/data/PNNL/ov_somatic_mutations/EventLocations_PNNL_20140202_removed_dbSNP_somatic.txt'#'/home/s3cha/data/VU/Somatic/Event_VU_SPECFDR01_splice_vcf_dbSNP_somatic.txt'#'/home/s3cha/data/PNNL/Event_PNNL_Single_recalculated.txt'#'/home/s3cha/Dropbox/Workspace/SpliceGraph/src/Enosi/sample_event2.txt'#Event_PNNL_20140202_Concat.txt'
    unknown_location_filename = '/home/s3cha/data/VU/Normal/Location_Normal_colon_SV_merged_somatic.txt'#'/home/s3cha/data/VU/Cosmic/Location_VU_SPECFDR01_Cosmicoverlap2.txt'#'/home/s3cha/data/VU/Somatic/Location_VU_RefSeq_PEPFDR01_rna.txt'#'/home/s3cha/data/PNNL/ov_somatic_mutations/Locations_PNNL_20140202_removed_dbSNP_somatic.txt'#'/home/s3cha/data/VU/Somatic/Location_VU_SPECFDR01_splice_vcf_dbSNP_somatic.txt'#'/home/s3cha/data/JHU/Location_JHU_SV_Result_SVrna_dbSNP_somatic.txt'#'/home/s3cha/data/PNNL/Location_PNNL_Single_recalculated.txt'#'C:/Users/Akard3/Dropbox/Workspace/SpliceGraph/src/Enosi/Locations_PNNL_20140202_removed_dbSNP.txt'#Location_PNNL_Single_recalculated.txt'#sample.txt'#
    known_location_filename = 'Location_Combine123_RefSeq_recalculated_Correction.txt'#'/home/s3cha/data/VU/Location_VU_RefSeq_PEPFDR01.txt'#
    ref_seq_filename = 'Homo_sapiens.GRCh37.70_formatted.gff'
    output_filename = '/home/s3cha/data/VU/Normal/Event_Normal_colon_SV_merged_somatic.txt'#'/home/s3cha/data/PNNL/ov_somatic_mutations/EventLocations_PNNL_20140202_removed_dbSNP_somatic.txt'#'/home/s3cha/data/VU/Somatic/Event_VU_SPECFDR01_splice_vcf_dbSNP_somatic.txt'#'/home/s3cha/data/PNNL/Event_PNNL_Single_recalculated.txt'#'/home/s3cha/Dropbox/Workspace/SpliceGraph/src/Enosi/sample_event2.txt'#Event_PNNL_20140202_Concat.txt'
    
    DNA_folder = '/home/s3cha/s3cha/data/dna/70'



print 'END: ',time.ctime()
print 'TIME:',time.clock()-betime,' sec'


                         
