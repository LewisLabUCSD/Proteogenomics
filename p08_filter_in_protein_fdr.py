import pandas as pd

group = 1
path = '/data/s3cha/CHO_ENOSI_JOB/e22076b9e7ab40ef8fc2c63eddf67eba/group'+ str(group) +\
        '/out_tsv/'
raw_fn = path + '/sequence1/sequence1.tsv'
pep_fdr_fn = path + '/Known_pep_fdr.txt'

raw_df = pd.read_csv(raw_fn, sep='\t')
# add a column to test if it maps to decoy or not
raw_df['isdecoy'] = raw_df['Protein'].map(lambda x: 'T' if x.startswith('XXX') else 'F')
# sort by evalue
raw_df = raw_df.sort_values(by=['EValue'])
raw_df = raw_df.reset_index(drop=True)

# calculate fdr
visited_target = set()
visited_decoy = set()
fdr = []
for isdecoy, pr in zip(raw_df['isdecoy'], raw_df['Protein']):
    if pr.startswith('XXX'):
        visited_decoy.add(pr)
    else:
        visited_target.add(pr)
    fdr.append(len(visited_decoy)/float(len(visited_target)))

raw_df['proFDR'] = pd.Series(fdr)

pep_fdr_df = pd.read_csv(pep_fdr_fn, sep='\t')

pro_fdr_df = pd.merge(pep_fdr_df, raw_df[['#SpecFile','SpecID','proFDR']],
         on=['#SpecFile','SpecID'])

pro_fdr_df = pro_fdr_df.query('proFDR < 0.01')
pro_fdr_df.to_csv(path+'/known_pro_fdr.txt', sep='\t', index=False)