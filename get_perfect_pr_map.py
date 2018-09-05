from Bio import SeqIO
import re,sys
import time
start = time.time()

# query_fn = '/data/shangzhong/Picr_assembly/Annotation/PASA/pasa_stringtie/03_pasa_stringtie_pr.fa'
# ref_fn = '/data/shangzhong/Picr_assembly/Annotation/PASA/pasa_stringtie/hamster_pr.fa'
# out = '/data/shangzhong/Picr_assembly/Annotation/PASA/pasa_stringtie/12_perfect_draft_prs2refseq1.txt'
query_fn = sys.argv[1]
ref_fn = sys.argv[2]
out = sys.argv[3]


pr_seqs = ''
pr_pos = [0]
prs = []
for record in SeqIO.parse(ref_fn, 'fasta'):
    pr_seqs += str(record.seq) + '#'
    pr_pos.append(len(pr_seqs))
    prs.append(record.id)
pr_seqs = '#' + pr_seqs


with open(out,'w') as f:
	for record in SeqIO.parse(query_fn,'fasta'):
		seq = '#' + str(record.seq) + '#'
		try:
			idxes = [m.start() for m in re.finditer(seq,pr_seqs)]
			for idx in idxes:
				new_idx = pr_pos.index(idx)		
				f.write(record.id + '\t' + str(prs[new_idx]) + '\n')
		except:
			continue

print time.time() - start