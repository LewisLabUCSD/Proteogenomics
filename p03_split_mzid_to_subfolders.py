'''
Created on Jun 13, 2017

@author: s3cha
'''

import sys
import os
import argparse


parser = argparse.ArgumentParser(description='Converts spectrum to indexed spectrum file')
parser.add_argument('-m','--mzid',action='store',dest='mzid',help='MSGF result mzid folder',default='MSGFresult')
parser.add_argument('-o','--out',action='store',dest='out',help='output folder,mzid correspond to one spectrum file will goes into the same folder',
					default='MSGFsplit')
args = parser.parse_args()


mzid_folder = args.mzid
out_folder = args.out
if not os.path.exists(out_folder):os.mkdir(out_folder)

dbs = os.listdir(mzid_folder)
for db in dbs:
	out = out_folder+'/'+db
	if not os.path.exists(out):os.mkdir(out)
	f_list = os.listdir(mzid_folder+'/'+db)
	file_dict = {}

	for file in f_list:
		spec_name = file.split('_')[0]
		if spec_name not in file_dict:
			file_dict[spec_name] = []
		file_dict[spec_name].append(file)
	for spec_name in file_dict:
		os.system('mkdir '+out+'/'+spec_name)
		for file in file_dict[spec_name]:
			os.system('mv '+mzid_folder+'/'+db+'/'+file+' '+out+'/'+spec_name)

