import pandas as pd
import argparse

def filter_known_psm_pep_level(known_res, known_pep_res):
    '''this function filter peptide in 1% FDR level'''
    df = pd.read_csv(known_res,sep='\t')
    print 'there are', df.shape[0],'PSM in 1% FDR'
    # extract peptides
    df['pep'] = df['Peptide'].map(lambda x: ''.join([p for p in x if p.isalpha()][:-1]))
    df = df[df['pep'].map(lambda x: len(x) >= 9)]
    # filter in peptide level at 1& FDR
    pep_df = df.query('PepFDR < 0.01')
    print 'there are', pep_df['pep'].unique().shape[0], 'peptide in 1% FDR'
    del pep_df['pep']
    pep_df.to_csv(known_pep_res, sep='\t', index=False)


parser = argparse.ArgumentParser(description='filter PSM in peptide-level')
parser.add_argument('-i','--in',action='store',dest='input',help='PSM result in 1% PSM FDR')
parser.add_argument('-o','--out',action='store',dest='output',help='PSM result in 1% peptide level FDR')

args = parser.parse_args()

in_fn = args.input
out_fn = args.output

filter_known_psm_pep_level(in_fn, out_fn)

# # filter peptide level for known peptide in draft only
# fn = '/data/shangzhong/Proteogenomics/event_results/known_res.txt'
# fn1 = '/data/shangzhong/Proteogenomics/event_results/known_pep_res.txt'
# filter_known_psm_pep_level(fn, fn1)

# # filter peptide level for known peptide in draft and refseq 
# fn = '/data/s3cha/CHO_ENOSI_JOB/dft_ref/known_res.txt'
# fn1 = '/data/s3cha/CHO_ENOSI_JOB/dft_ref/known_pep_res.txt'
# filter_known_psm_pep_level(fn, fn1)


