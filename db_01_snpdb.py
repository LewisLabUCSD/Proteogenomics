from Bio import SeqIO
import sarge,os,glob
from natsort import natsorted
from multiprocessing import Pool
import pandas as pd
fa = '/data/genome/hamster/picr_old/picr.fa'
gff = '/data/genome/hamster/picr_old/pasa_stringtie.gff3'
snpDB_path = '/data/shangzhong/Proteogenomics/Database/SNP_Database_Generate'
snpDB_code_path = '/data/shangzhong/Proteogenomics/Database/SNP_Database_Generate'
chr_num = 1830
thread = 24

os.chdir(snpDB_path)
graph = snpDB_code_path + '/bin/ConstructSpliceGraphVariant.jar'
dna_folder = snpDB_path + '/dna'
gff_folder = snpDB_path + '/gff'
vcf_folder = snpDB_path + '/vcf'
spl_folder = snpDB_path + '/spl'
msdb_folder= snpDB_path + '/msdb'
dbfa_folder = snpDB_path + '/dbfa'
if not os.path.exists(dna_folder): os.mkdir(dna_folder)
if not os.path.exists(gff_folder): os.mkdir(gff_folder)
if not os.path.exists(spl_folder): os.mkdir(spl_folder)
if not os.path.exists(msdb_folder): os.mkdir(msdb_folder)
if not os.path.exists(dbfa_folder): os.mkdir(dbfa_folder)
#----- 1. convert vcf to spl
out_spl = snpDB_path+'/spl.out'
cmd = ('python {spl_code} temp {out_spl} {vcf_folder}').format(spl_code=snpDB_code_path+'/bin/SPLgenerator.py',out_spl=out_spl,vcf_folder=vcf_folder)
print(cmd)
sarge.run(cmd)
# # 1) split spl
all_picr_vcfs = [[] for n in range(chr_num)]
with open(out_spl) as fin:
    for line in fin:
        if line.startswith('#'):
            for n in range(chr_num):
                all_picr_vcfs[n].append(line)
        else:
            scaff = int(line.split('\t')[0][5:])
            all_picr_vcfs[scaff].append(line)
for n in range(chr_num):
    scaff = 'picr_'+str(n)
    spl = spl_folder + '/'+scaff+'.spl'
    with open(spl,'w') as f:
        for line in all_picr_vcfs[n]:
            f.write(line)
# ----- 2. creates ms2db
# 1) split fa files
seqh = SeqIO.index(fa,'fasta')
for record in seqh:
    out_fa = dna_folder+'/'+record+'.fa'
    SeqIO.write(seqh[record],out_fa,'fasta')
# 2) split gff files
def write_gff(df):
    chrom = df.iloc[0,0]
    outfn = gff_folder + '/' + chrom + '.gff'
    df.to_csv(outfn,sep='\t',header=None,index=False)
      
gff_df = pd.read_csv(gff,sep='\t',header=None)
gff_df.groupby(0).apply(write_gff)
  
# 3) ms2db
def run_cmd(cmd):
    print(cmd)
    sarge.run(cmd)
gffs = natsorted([f for f in os.listdir(gff_folder) if f.endswith('.gff')])
p = Pool(processes=thread)
for gff in gffs:
    fa = dna_folder + '/' + gff[:-4] + '.fa'
    msdb = msdb_folder + '/'+gff[:-4] + '.ms2db'
    spl = spl_folder + '/' + gff[:-4] + '.spl'
    cmd = ('java -jar {graph} -p {spl} -w {msdb} -r {gff} -s {fa}').format(graph=graph,spl=spl,msdb=msdb,gff=gff_folder+'/'+gff,fa=fa)
    p.apply_async(run_cmd,args=(cmd,))
p.close()
p.join()
#--------- 3. ms2db to fasta
code = snpDB_code_path + '/bin/ACGT03102013.py'
msdbs = natsorted(glob.glob(msdb_folder+'/*.ms2db'))
p = Pool(processes=thread)
for msdb in msdbs:
    out_fa = dbfa_folder + '/' + os.path.basename(msdb).split('.')[0] + '.fa'
    cmd = ('python {code} {msdb} {out_fa} 30').format(code=code,msdb=msdb,out_fa=out_fa)
    p.apply_async(run_cmd,args=(cmd,))
p.close()
p.join()
#--------- 4. reduce snp db
fas = natsorted(glob.glob(dbfa_folder+'/*.fa'))
merge_fa = snpDB_path + '/snpdb.fa'
cmd = ('cat {fas} > {merge}').format(fas=' '.join(fas),merge=merge_fa)
sarge.run(cmd)
reduce_code = snpDB_code_path + '/bin/SNPreduce.py'
reduce_fa = snpDB_path+'/snpdb_reduce.fa'
cmd = ('python {c} {i} {o}').format(c=reduce_code,i=merge_fa,o=reduce_fa)
sarge.run(cmd)
# #-------- 5. remove duplicated sequence
# sed -e '/^>/s/$/*/' -e 's/^>/#/' snpdb_reduce.fa | tr -d '\n' | tr "#" "\n" | tr "*" "\t" | sort -u -t$'\t' -f -k 2,2  | sed -e 's/^/>/' -e 's/\t/\n/'