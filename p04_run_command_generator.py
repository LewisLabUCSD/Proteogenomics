import os
import argparse

path = os.getcwd()

parser = argparse.ArgumentParser(description='Generate command file to transfer mzid to tsv')
parser.add_argument('-p','--path',action='store',dest='path',help='path to the main foler of the pipeline',default=path)
parser.add_argument('-c','--code',action='store',dest='code',help='path to python code of mzidtotsv',
                             default='/data/s3cha/CHO_ENOSI_JOB/tool/2014.1208/MzidToTsv_withmod.py')
parser.add_argument('-s',action='store',dest='split',help='split folder that has all mzid files',default='MSGFsplit')
parser.add_argument('-o',action='store',dest='out',help='output folder that store all tsv files',default='out_tsv')
parser.add_argument('-m',action='store',dest='cmd',help='file that store all command to run',default='tsv_cmd.txt')
args = parser.parse_args()

path = args.path
python_code = args.code
split_folder = args.split
output_folder = args.out
output_command = args.cmd

if not os.path.exists(output_folder): os.mkdir(output_folder)

sequence_list = os.listdir(split_folder)

s = open(output_folder+'/'+output_command,'w')

for sequence in sequence_list:
    os.system('mkdir '+output_folder+'/'+sequence)
    folder_list = os.listdir(split_folder+'/'+sequence)
    for folder in folder_list:
        s.write('python '+python_code+' '+path+'/'+split_folder.split('/')[-1]+'/'+sequence+'/'+folder+' '+
                path+'/'+output_folder.split('/')[-1]+'/'+sequence+'/'+sequence+'_'+folder+'.tsv\n')


