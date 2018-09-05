'''this code filters the novel results identified using MSGF+.
for a peptide starts or ends with -, which means the peptide is the start or end
of the reference protein. MSGF+ considers them the same as R or K. but for spliceDB
or snpDB, many protein seqeunces are not the real start or end of a protein. eg:
for a identified peptide  -.ALK.R, in known database, a protein includes ALKR. Since 
in the known protein, the AA before A is not K or R, so MSGF+ would not consider it
when mapping to known database, however, when mapping to novel database, A is the 
beginning of the peptide, so MSGF considers it as mapping, so what we are doing here
is to map the identified peptides to the known database, if the peptide maps perfectly to 
known database, it means the peptide is actually known peptide instead of novel.
'''
import pandas as pd
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Filter Novel results from the Known results')
parser.add_argument('-p',action='store',dest='tsv',help='path to tsv folder and cmd file',default='out_tsv')
parser.add_argument('-k',action='store',dest='known',help='full path to known fa database')

args = parser.parse_args()

# path = '/data/s3cha/CHO_ENOSI_JOB/e22076b9e7ab40ef8fc2c63eddf67eba/group1/out_tsv'
path = args.tsv
known_fa_fn = args.known

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
	if pep.startswith('-') or pep.endswith('-'):
		pep = ''.join([p for p in pep if p.isalpha()])
		if pep in known_peps:
			keep = False
	return keep

known_peps = ''
for record in SeqIO.parse(known_fa_fn,'fasta'):
    known_peps += str(record.seq) + 'X'
cri = novel_df['Peptide'].map(lambda x:overlap(x,known_peps))
novel_pst_df = novel_df[cri]
novel_pst_df.to_csv(novel[:-4]+'_post.txt',sep='\t',index=False)

novel_known_df = novel_df[~cri]
novel_known_df.to_csv(novel[:-4]+'_novel2known.txt',sep='\t',index=False)
# known_peps = ''
# for pep in known_df['Peptide']:
#     peptide = ''.join([p for p in pep if p.isalpha()])
#     known_peps += peptide + 'X'
# 