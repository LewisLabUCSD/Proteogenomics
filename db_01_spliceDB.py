import os,glob,sarge,subprocess
import multiprocessing as mp
import shutil
from natsort import natsorted
from Bio import SeqIO



'''There are two strategies. 1. copy all bam files to a foler and then split them into 
each scaffold, thus in each scaffold folder, it has all bam files. then run the pipeline'''
def bam2splitBam(raw_path,target_path,scaffs):
    bams = natsorted(glob.glob(raw_path+'/*.bam'))
    for bam in bams:
        out = target_path+'/' + os.path.split(bam)[1]
        if not os.path.exists(out):
            cmd = ('cp \'{i}\' {o}').format(i=bam,o=out)
            print cmd
            sarge.run(cmd)
            cmd = ('cp \'{i}.bai\' {o}.bai').format(i=bam,o=out)
            print cmd
            subprocess.call(cmd,shell=True)
        for scaff in scaffs:
            bam_path = target_path + '/bam'+scaff
            if not os.path.exists(bam_path): os.mkdir(bam_path)
            out_fn = bam_path+'/'+bam.split('/')[-1]
            if not os.path.exists(out_fn):
                cmd = ('samtools view -bh \'{i}\' {c} > {out}').format(i=bam,c=scaff,out=out_fn)
                print cmd
                subprocess.call(cmd,shell=True)
        os.remove(out)
        os.remove(out+'.bai')   

def run_spliceDB_old(scaff,raw_path,target_path,code_path):
    '''
    * scaff: scaffold 
    * target_path: folder that have all results
    * code_path: splice db codes path
    '''
    global ref_dic
    # create folder
    bam_path = target_path + '/bam'+scaff
    if not os.path.exists(bam_path): 
        assert False,'file not exitts'
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
    
    fa = dna_path+'/'+scaff+'.fa'
    if not os.path.exists(fa):
        with open(dna_path+'/'+scaff+'.fa','w') as handle:
            SeqIO.write(ref_dic[scaff], handle,'fasta')
    # run splice db
    buildsp = code_path + '/buildSpliceGraph.py'
    sam = "/opt/miniconda2/bin" # samtools path
    cmd = 'python {build} -a {bam} -p {spl} -t {sam} -f {fa} -c {c}'.format(build=buildsp,bam=bam_path,spl=spl_path,sam=sam,fa=dna_path,c=code_path)
    subprocess.call(cmd,shell=True)
    # move fa file to prfa folder
    fa = glob.glob(spl_path+'/*.fa')
    for f in fa:
        shutil.move(f,fa_res_path)
    # remove those folders
    shutil.rmtree(bam_path)
    # shutil.rmtree(spl_path)
    shutil.rmtree(dna_path)


# import time
# start = time.time()


# # raw_path = '/data/shangzhong/GoogleDrive/Research_Project/CHO_Project/Bam_files/Proteogenomics'
# raw_path = '/media/lewislab/Dropbox (UCSD SBRG)/users/shangzhong/proteogenomics_bam'
# target_path = '/data/shangzhong/Proteogenomics/Database/Splice_Db'
# ref_fa = '/data/genome/hamster/picr_old/picr.fa'
# ref_dic = SeqIO.index(ref_fa,'fasta')


# scaffs = ['picr_'+str(n) for n in [0,1]]
# code_path = '/home/shangzhong/Codes/Proteogenomics/SpliceDB/bin'
# thread = 9

# # bam2splitBam(raw_path,target_path,list(ref_dic.keys()))

# pool = mp.Pool(processes=thread)
# for scaff in scaffs:
#     pool.apply_async(run_spliceDB,args=(scaff,raw_path,target_path,code_path,))
# pool.close()
# pool.join()


# print time.time() - start

# before this line is code to copy bam to your server first, the following is directly extract chr bams.



'''The second method is directly extract bam file for each scaffold and then run the pipeline.
'''
def bam2splitBam(raw_path,target_path,scaff):
    bams = natsorted(glob.glob(raw_path+'/*.bam'))
    for bam in bams:
        out_fn = target_path+'/'+bam.split('/')[-1]
        if not os.path.exists(out_fn):
            cmd = ('samtools view -bh \'{i}\' {c} > {out}').format(i=bam,c=scaff,out=out_fn)
            subprocess.call(cmd,shell=True)
            print cmd



def wirte_fa(target_path,ref_dic,scaffs):
    ''''''
    for scaff in scaffs:
        dna_path = target_path + '/dna_'+scaff
        if os.path.exists(dna_path): 
            shutil.rmtree(dna_path)
        os.mkdir(dna_path)
        fa = dna_path+'/'+scaff+'.fa'
        if not os.path.exists(fa):
            with open(dna_path+'/'+scaff+'.fa','w') as handle:
                SeqIO.write(ref_dic[scaff], handle,'fasta')

def run_spliceDB(scaff,raw_path,target_path,code_path):
    '''
    * scaff: scaffold 
    * target_path: folder that have all results
    * code_path: splice db codes path
    '''
    # create folder
    bam_path = target_path + '/bam_'+scaff
    if os.path.exists(bam_path): 
        shutil.rmtree(bam_path)
    os.mkdir(bam_path)
    spl_path = target_path + '/spl_'+scaff
    if os.path.exists(spl_path): 
        shutil.rmtree(spl_path)
    os.mkdir(spl_path)
    dna_path = target_path + '/dna_'+scaff
    # if os.path.exists(dna_path): 
    #     shutil.rmtree(dna_path)
    # os.mkdir(dna_path)
    fa_res_path = target_path+'/prfa'
    if not os.path.exists(fa_res_path): os.mkdir(fa_res_path)
    # create fa file
    # global ref_dic
    fa = dna_path+'/'+scaff+'.fa'
    if not os.path.exists(fa):
        with open(dna_path+'/'+scaff+'.fa','w') as handle:
            SeqIO.write(ref_dic[scaff], handle,'fasta')
    bam2splitBam(raw_path,bam_path,scaff)
    # run splice db
    buildsp = code_path + '/buildSpliceGraph.py'
    sam = "/usr/bin" # samtools path
    cmd = 'python {build} -a {bam} -p {spl} -t {sam} -f {fa} -c {c}'.format(
        build=buildsp,bam=bam_path,spl=spl_path,sam=sam,fa=dna_path,c=code_path)
    subprocess.call(cmd,shell=True)
    print cmd
    # move fa file to prfa folder
    fa = glob.glob(spl_path+'/*.fa')
    for f in fa:
        shutil.move(f,fa_res_path+'/'+f.split('/')[-1])
    # remove those folders
    shutil.rmtree(bam_path)
    # shutil.rmtree(spl_path)
    shutil.rmtree(dna_path)

import time
start = time.time()

# bam2chrbam(raw_path,target_path,list(ref_dic.keys()))


raw_path = '/media/lewislab/Dropbox (UCSD SBRG)/users/shangzhong/proteogenomics_bam'
target_path = '/data/shangzhong/Proteogenomics/Database/Splice_Db_old_samtools'
if not os.path.exists(target_path): os.mkdir(target_path)
ref_fa = '/data/genome/hamster/picr_old/picr.fa'
ref_dic = SeqIO.index(ref_fa,'fasta')


chrom = range(1829)
scaffs = ['picr_'+str(n) for n in chrom]
code_path = '/home/shangzhong/Codes/Proteogenomics/SpliceDB'
thread = 1

wirte_fa(target_path,ref_dic,scaffs)

pool = mp.Pool(processes=thread)
for scaff in scaffs:
    pool.apply_async(run_spliceDB,args=(scaff,raw_path,target_path,code_path))
pool.close()
pool.join()


print time.time() - start