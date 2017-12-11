import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Filter Novel results from the Known results')
parser.add_argument('-p',action='store',dest='tsv',help='path to tsv folder and cmd file',default='out_tsv')

args = parser.parse_args()

# path = '/data/s3cha/CHO_ENOSI_JOB/e22076b9e7ab40ef8fc2c63eddf67eba/group1/out_tsv'
path = args.tsv

novel = path +'/Novel.txt'
known = path +'/Known.txt'

known_df = pd.read_csv(known,sep='\t',header=0,low_memory=False,na_filter=False)
novel_df = pd.read_csv(novel,sep='\t',header=0,low_memory=False,na_filter=False)
known_df['id'] = known_df['#SpecFile'] + known_df['SpecID']
novel_df['id'] = novel_df['#SpecFile'] + novel_df['SpecID']

known_ids = known_df['id'].tolist()
novel_df = novel_df[~novel_df['id'].isin(known_ids)]
del novel_df['id']
novel_df.to_csv(novel[:-3]+'pre.txt',sep='\t',index=False)

# get the known peptides
def overlap(pep,known_peps):
	'''check if novel peptides in knwon peptides
	'''
	keep = True
	if pep.startswith('-'):
		pep = ''.join([p for p in pep if p.isalpha()])
		if pep in known_peps:
			keep = False
		# for peptide in known_peps:
		# 	if pep in peptide:
		# 		keep=False
				# return keep
	return keep

known_peps = ''
for pep in known_df['Peptide']:
	peptide = ''.join([p for p in pep if p.isalpha()])
	known_peps += peptide + 'X'

novel_df = novel_df[novel_df['Peptide'].map(lambda x:overlap(x,known_peps))]
novel_df.to_csv(novel+'1',sep='\t',index=False)
