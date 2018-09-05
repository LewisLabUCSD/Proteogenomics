import pandas as pd
import re

gff = '/data/genome/hamster/ncbi_refseq/hamster.gff'
#===============================================================================
#             1. some statistics of annotation file
#===============================================================================
stats_dic = {}
# total gene number
gff_df = pd.read_csv(gff,sep='\t',header=None,comment='#',names=['chr','s',
                    'feature','start','end','dot','strand','none','anno'])
gff_df = gff_df.drop(['s','dot','none'],axis=1)
gene_df = gff_df[gff_df['feature'].values=='gene']
gene_number = gene_df.shape[0]
stats_dic['total_gene_number'] = gene_number
# add geneid and gene symbol column
gene_df = gene_df.reset_index(drop=True)
gene_df.loc[:,'gene'] = gene_df['anno'].apply(lambda x: re.search('(?<=gene=).+?(?=;|$)',x).group())
gene_df.loc[:,'geneid'] = gene_df['anno'].apply(lambda x: re.search('(?<=GeneID:).+?(?=[;,]|$)',x).group())


