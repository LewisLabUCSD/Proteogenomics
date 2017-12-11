import os,glob,sarge
import multiprocessing as mp
import shutil
from natsort import natsorted
from Bio import SeqIO




def bam2chrbam(raw_path,target_path,chroms):
    '''this function copies bam file and split it based on chromosome
    '''
    bams = natsorted(glob.glob(raw_path+'/*.bam'))
    for bam in bams[3:]:
        # copy file
        out = target_path+'/' + os.path.split(bam)[1]
        if not os.path.exists(out):
            cmd = ('cp {i} {o}').format(i=bam,o=out)
            print cmd
            sarge.run(cmd)
        if not os.path.exists(out+'.bai'):
            cmd = ('samtools index {b}').format(b=out)
            print cmd
            sarge.run(cmd)
        # split bam file
        bam_folder = target_path+'/'+os.path.split(bam)[1][:-4]
        if not os.path.exists(bam_folder):
            os.mkdir(bam_folder)
        for scaff in chroms:
            cmd = ('samtools view -bh {i} {c} > {out}').format(i=out,c=scaff,out=bam_folder+'/'+scaff+'.bam')
            print cmd
            sarge.run(cmd)
        os.remove(out)
            




def run_spliceDB(scaff,raw_path,target_path,code_path):
    '''
    * scaff: scaffold 
    * target_path: folder that have all results
    * code_path: splice db codes path
    '''
    # create folder
    bam_path = target_path + '/bam'+scaff
    if os.path.exists(bam_path): 
        shutil.rmtree(bam_path)
    os.mkdir(bam_path)
    spl_path = target_path + '/spl'+scaff
    if os.path.exists(spl_path): 
        shutil.rmtree(spl_path)
    os.mkdir(spl_path)
    dna_path = target_path + '/dna'+scaff
    if os.path.exists(dna_path): 
        shutil.rmtree(dna_path)
    os.mkdir(dna_path)
    fa_res_path = target_path+'/prfa'
    if not os.path.exists(fa_res_path): os.mkdir(fa_res_path)
    bam2chrbam(raw_path,target_path,[scaff])
    # copy bam file and extract chromosome specific bam file
    bams = glob.glob(raw_path+'/*/'+scaff+'.bam')
    for bam in bams:
        sarge.run('cp {before} {after}'.format(before=bam,after=bam_path+'/'+bam.split('/')[-2]+'.bam'))
    # create fa file
    global ref_dic
    fa = dna_path+'/'+scaff+'.fa'
    if not os.path.exists(fa):
        with open(dna_path+'/'+scaff+'.fa','w') as handle:
            SeqIO.write(ref_dic[scaff], handle,'fasta')
    # run splice db
    buildsp = code_path + '/buildSpliceGraph.py'
    sam = "/usr/local/bin" # samtools path
     
    cmd = 'python {build} -a {bam} -p {spl} -t {sam} -f {fa} -c {c}'.format(build=buildsp,bam=bam_path,spl=spl_path,sam=sam,fa=dna_path,c=code_path)
    sarge.run(cmd)
    # move fa file to prfa folder
    fa = glob.glob(spl_path+'/*.fa')
    for f in fa:
        shutil.move(f,fa_res_path)
    # remove those folders
    shutil.rmtree(bam_path)
    shutil.rmtree(spl_path)
    shutil.rmtree(dna_path)


import time
start = time.time()

# bam2chrbam(raw_path,target_path,list(ref_dic.keys()))

# raw_path = '/data/shangzhong/GoogleDrive/Research_Project/CHO_Project/Bam_files/Proteogenomics'
raw_path = '/media/lewislab/Dropbox (UCSD SBRG)/users/shangzhong/proteogenomics_bam'
target_path = '/data/shangzhong/Proteogenomics/Database/Splice_Db'
ref_fa = '/data/genome/hamster/picr_old/picr.fa'
ref_dic = SeqIO.index(ref_fa,'fasta')

scaffs = ['picr_'+str(n) for n in range(1829)]
code_path = '/home/shangzhong/Codes/Proteogenomics/SpliceDB/bin'
thread = 9

pool = mp.Pool(processes=thread)
for scaff in scaffs:
    pool.apply_async(run_spliceDB,args=(scaff,raw_path,target_path,code_path,))
pool.close()
pool.join()


print time.time() - start