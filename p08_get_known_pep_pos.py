'''
Created on Jul 7, 2017

@author: s3cha
'''
import re,sys
from Bio import SeqIO
import pandas as pd
import argparse


class TagSearch(object):
    tag_table = {}
    tag_length = None
    def __init__(self):
        pass
    
    def makeTagTable(self,fasta,tag_length):
        self.tag_length = tag_length
        for index in range(len(fasta)-tag_length+1):
            tag = fasta[index:index+tag_length]
            if tag not in self.tag_table:
                self.tag_table[tag] = []
            self.tag_table[tag].append(index)
        pass
    
    def getAllLocation(self,peptide,fasta):
        if len(peptide) < self.tag_length:
            return [m.start() for m in re.finditer(peptide,fasta)]
        else:
            result = []
            tag = peptide[:self.tag_length]
            candidates = self.tag_table.get(tag)
            if candidates == None:
                return result
            for candi in candidates:
                if fasta[candi:candi+len(peptide)] == peptide:
                    result.append(candi)
            return result
            
def binary_search(array,key,imin,imax):
    if (imax < imin):
        return imax
    else:
        imid = (imin+imax)/2
        if array[imid] > key:
            return binary_search(array,key,imin,imid-1)
        elif array[imid] < key:
            return binary_search(array,key,imid+1,imax)
        else:
            return imid

def find_pos(ref,pep):
    ref = str(ref).replace('I','L')
    pep = pep.replace('I','L')
    try:
        return ref.index(pep)
    except:
        print ref,pep
        return None


#============== parameters =====================
parser = argparse.ArgumentParser(description='Get peptie position in proteins')
parser.add_argument('-k','--known',action='store',dest='known',help='fa file of known protein database')
parser.add_argument('-r','--res',action='store',dest='res',help='proteomics result file')
parser.add_argument('-o','--out',action='store',dest='out',help='outfile store the position of peptide')
args = parser.parse_args()
known_DB = sys.argv[1] # reference fa file
known_res = sys.argv[2] # merged known peptide results from proteogenomics
out_fn = sys.argv[3] # have rna, pep, index
known_DB = args.known
known_res = args.res
out_fn = args.out

# known_DB = '/data/shangzhong/Picr_assembly/Annotation/PASA/pasa_stringtie/03_pasa_stringtie_pr.fa'
# known_res = '/data/shangzhong/Proteogenomics/event_results/known_res.txt'
# out_fn = '/data/shangzhong/Proteogenomics/event_results/known_pos.txt'


x = TagSearch()
fasta_database = open(known_DB,'r')
fasta_seq = ''
fasta_info = []
index_info = [0]
total_len = 0
for record in SeqIO.parse(fasta_database,'fasta'):
    fasta_info.append(record.id)
    seq = str(record.seq)
    fasta_seq += seq.replace('I','L') + 'X'
    total_len += len(seq) + 1
    index_info.append(total_len)
fasta_seq = fasta_seq[:-1]

tagS = TagSearch()
tagS.makeTagTable(fasta_seq,4)



peptides = []
with open(known_res) as f:
    for line in f:
        if line.startswith('#'): continue
        item = line.split('\t')
        pep = ''.join([p for p in item[8] if p.isalpha()][:-1])
        peptides.append(pep)
peptides = list(set(peptides))

with open(out_fn,'w') as f:
    for pep in peptides:
        locis = tagS.getAllLocation(pep.replace('I','L'),fasta_seq)
        for loci in locis:
            idx = binary_search(index_info,loci,0,len(index_info)-1)
            header = fasta_info[idx]
            rna = header.split(' ')[0]
            f.write(rna+'\t'+pep+'\n')


pr_index = SeqIO.index(known_DB,'fasta')
df = pd.read_csv(out_fn,sep='\t',header=None)
df = df.drop_duplicates()
df['pos'] = df.apply(lambda x:find_pos(pr_index[x[0]].seq,x[1]),axis=1)
df.to_csv(out_fn,sep='\t',header=None,index=False)


known_res_df = pd.read_csv(known_res,sep='\t',header=0,low_memory=False,usecols=[0,1,8])
known_res_df['Peptide'] = known_res_df['Peptide'].map(lambda x: ''.join([p for p in x if p.isalpha()][:-1]))
pep_spec_num_dic = known_res_df['Peptide'].value_counts()
out_df = pd.read_csv(out_fn,sep='\t',header=None,names=['rna','pep','pos'])
out_df['spec_num'] = out_df['pep'].map(lambda x:pep_spec_num_dic[x])

out_df = out_df.query('pos == pos')
out_df = out_df.reset_index(drop=True)
out_df['pos'] = out_df['pos'].astype('int')
out_df.to_csv(out_fn,sep='\t',index=False)