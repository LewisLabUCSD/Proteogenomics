import pandas as pd
import os
import sarge
from Bio import SeqIO,Seq
from natsort import natsorted
import shutil,yaml
import glob
import sys
import multiprocessing as mp
"""This file prepare the bam file for spliceDB
   and vcf file for the pipeline"""


def chunk(l,n):
    n = max(1,n)
    result = [l[i:i+n] for i in range(0,len(l),n)]
    return result

def copy_fq_files(path,target_path,index=''):
    '''
    path: path of raw fastq files
    target_path: where to copy fastq files to
    index: should be a list of integers, '' means it copy all samples
    '''
    os.chdir(path)
    folders = [f for f in os.listdir(path) if os.path.isdir(f)]
    folders = natsorted(folders)
    if index == '':
        index = range(0,len(folders))
    
    for folder in [folders[i] for i in index]:
        fq_path = path + '/' + folder
        os.chdir(fq_path)
        fqFiles = [f for f in os.listdir(fq_path) if f.endswith('.fastq.gz')]
        fqFiles = natsorted(fqFiles)
        fst = [f for f in fqFiles if 'R1' in f]
        snd = [f for f in fqFiles if 'R2' in f]
        
        cmd = ('cat {input} > {tar}/{folder}_1.fq.gz').format(
                    input=' '.join(fst),folder=folder,tar=target_path)
        print(cmd)
        sarge.run(cmd)
        
        cmd = ('cat {input} > {tar}/{folder}_2.fq.gz').format(
                    input=' '.join(snd),folder=folder,tar=target_path)
        print(cmd)
        sarge.run(cmd)

#===============================================================================
#                 1. run RNAseq count pipeline for copied RNAseq and extract alignment based on chromosome
#===============================================================================
def run_star(indexes,rnaseq_path,target_path,rna_pipeline_fn,parameter_fn,splice_path,ps):
    '''indexes is a list, and each item in it is also a list
    '''
    for index in indexes:
        # 1. copy files
        copy_fq_files(rnaseq_path,target_path,index)
        # 2. map using STAR
        os.chdir(target_path)
        cmd = ('python {pipe} {pipe_param}').format(pipe=rna_pipeline_fn,pipe_param=parameter_fn)
        sarge.run(cmd)
        # 3. move files to target files
        if ps == 1:
            res_fns = glob.glob('sortBam/*.tab')
        elif ps == 2:
            res_fns = glob.glob('sortBam/*.bam')
        for f in res_fns:
            shutil.move(f,splice_path)
        shutil.rmtree(target_path);os.mkdir(target_path)

def rna_get_bam(rnaseq_path,target_path,final_path,batch,rna_pipeline_fn,parameter_fn,pass_num,indexes=[],chrom=''):
    """ This functions does the followings: 1. map reads to genome without annotation. 2. merge all .tab file. as annotation
    3. run 2nd round of STAR using annotation generated in step 2.
    * batch: int. How many RNAseq samples to run at each batch.
    * chrom: if not sepcified, will run for all chromosomes.
    * target_path: path for running variant calling for fq files.
    * final_path: path to store final results. should not inside target_path.
    * pass_num: indicates if this run is 1st pass or 2nd pass of the STAR run.
    """
    if os.path.exists(target_path): 
        shutil.rmtree(target_path)
    os.mkdir(target_path)
    if not os.path.exists(final_path): os.mkdir(final_path)
    splice_path = final_path+'/splice'
    if not os.path.exists(splice_path): os.mkdir(splice_path)
    os.chdir(rnaseq_path)
    folders = natsorted([f for f in os.listdir(rnaseq_path) if os.path.isdir(f)])
    if indexes != []:
        index = chunk(indexes,batch)
    else:
        index = chunk(range(len(folders)),batch)
    tabs = glob.glob(splice_path+'/*.tab')
    if pass_num == (1 or 'both'):
        # 1. run the first round to generate the splice junction files
        run_star(index,rnaseq_path,target_path,rna_pipeline_fn,parameter_fn,splice_path,1)
    elif (pass_num == (2 or 'both')) and (len(tabs) == len(folders)):
        # 2. run the second pass of star
        # change files
        with open(parameter_fn) as f:
            doc = yaml.load(f)
        doc['other_params'] = ['--sjdbFileChrStartEnd ' + ' '.join(tabs),'--limitSjdbInsertNsj 1999999']
        with open(parameter_fn,'w') as f:
            yaml.dump(doc,f)
        # run 2nd pass using STAR
        bam_path = final_path + '/bam'
        if not os.path.exists(bam_path): os.mkdir(bam_path)
        run_star(index,rnaseq_path,target_path,rna_pipeline_fn,parameter_fn,bam_path,2)
        # extract alignment by chromosome
        if chrom != '':
            os.chdir(final_path)
            if not os.path.exists(final_path+'/chrom'): os.mkdir(final_path+'/chrom')
            bams = [f for f in os.listdir(final_path) if f.endswith('.bam')]
            sam_cmd = ''
            rm_cmd = ''
            for bam in bams:
                out = bam[5:]
                sam_cmd = sam_cmd + ('samtools view -bh {input} {chr} > chrom/{out} & ').format(input=bam,chr=chrom,out=out)
                rm_cmd = rm_cmd + ('rm {fn} & rm {fn1} && ').format(fn=bam,fn1=bam+'.bai')
            sarge.run(sam_cmd[:-3])
            sarge.run(rm_cmd[:-3])

# path = '/media/lewislab/Dropbox (UCSD SBRG)/LewisPub/Data/DataForProteogenomics/DBH_CHO_mass_spec_data/YS_transfered_fastq_files'
# target_path = '/data/shangzhong/Proteogenomics/bam_fq'
# final_path = '/data/shangzhong/Proteogenomics/01_bam'
# batch= 8
# pass_num = 2
# rna_pipeline_file = '/home/shangzhong/Codes/NGS-Pipeline/STAR_get_bam.py'
# rna_pipeline_param = '/data/shangzhong/Proteogenomics/STAR_get_bam.yaml'
# indexes = range(36,68)
# 
# rna_get_bam(path,target_path,final_path,batch,rna_pipeline_file,rna_pipeline_param,pass_num,indexes,chrom='')

#===============================================================================
#                 2. run RNAseq variant calling pipeline
#===============================================================================
def rna_vari_call(rnaseq_path,target_path,final_path,batch,rna_vcf_pipeline_fn,parameter_fn,indexes=[],chrom =''):
    '''run RNAseq variant calling for all files
    * rnaseq_path: path of raw fq files
    * target_path: target path to copy fq files to
    * indexes: decide which fq files to include
    
    '''
    if os.path.exists(target_path): 
        shutil.rmtree(target_path)
    os.mkdir(target_path)
    if not os.path.exists(final_path): os.mkdir(final_path)
    os.chdir(rnaseq_path)
    folders = natsorted([f for f in os.listdir(rnaseq_path) if os.path.isdir(f)])
    if indexes != []:
        index = chunk(indexes,batch)
    else:
        index = chunk(range(len(folders)),batch)
    for i in index:
        # 1. copy files
        try:
            copy_fq_files(rnaseq_path,target_path,index=i)
        except:
            print('copy files failed')
        # 2. change read group in  paramter file
        os.chdir(target_path)
        with open(parameter_fn) as f:
            doc = yaml.load(f)
        
        sub_folders = [folders[j] for j in i]
        rg = ['@RG\\tID:'+folder+'\\tSM:'+folder for folder in sub_folders]
        print 'readgroup is:',rg
        doc['read_groups'] = rg
        with open(parameter_fn,'w') as f:
            yaml.dump(doc,f)
        # 3. run variant call pipeline
        cmd = ('python {pipe} {pipe_param}').format(pipe=rna_vcf_pipeline_fn,pipe_param=parameter_fn)
        sarge.run(cmd)
        # 4. move results to anohter place
        res_fns = natsorted(glob.glob('f12_FinalVcf/*'))
        for f in res_fns:
            shutil.move(f,final_path)
        shutil.rmtree(target_path);os.mkdir(target_path)
        # 5. extract alignment by chromosome
        if chrom != '':
            chr_path = final_path + '/chrom'
            if not os.path.exists(chr_path): os.mkdir(chr_path)
            os.chdir(final_path)
            vcfs = glob.glob('*.filter.vcf')
            for vcf in vcfs:
                out = 'chrom/'+vcf
                outhandle = open(out,'w')
                for line in open(vcf):
                    if line.startswith('#') or line.starswith(chrom):
                        outhandle.write('line')
                    outhandle.close()
                   
    
# rnaseq_path = '/media/lewislab/Dropbox (UCSD SBRG)/LewisPub/Data/DataForProteogenomics/DBH_CHO_mass_spec_data/YS_transfered_fastq_files'
# target_path = '/data/shangzhong/Proteogenomics/fq'
# final_path = '/data/shangzhong/Proteogenomics/vcf'
#   
# batch = 6
# indexes= range(12,68)
# rna_vcf_pipeline_fn = '/home/shangzhong/Codes/NGS-Pipeline/GATK_RNA_CHO.py'
# parameter_fn = '/data/shangzhong/Proteogenomics/vcf/GATK_RNA_CHO.yaml'
#   
# rna_vari_call(rnaseq_path,target_path,final_path,batch,rna_vcf_pipeline_fn,parameter_fn,indexes,chrom ='')

#===============================================================================
#                     3. run stringtie for RNAseq
#===============================================================================
'''
rnaseq_path = '/media/lewislab/Dropbox (UCSD SBRG)/LewisPub/Data/DataForProteogenomics/DBH_CHO_mass_spec_data/YS_transfered_fastq_files'
target_path = '/data/shangzhong/Proteogenomics'

stringtie_ppl_fn = '/home/shangzhong/Codes/NGS-Pipeline/StringTie_quant.py'
stringtie_par_fn = '/home/shangzhong/Codes/NGS-Pipeline/Parameters/StringTie_quant.yaml'

thread = 6
indexes = range(3,68)
batch = 11
index_batch = chunk(indexes,batch)
def run_stringtie_ppl(indexes,rna_path,target_path,ppl_fn,par_fn):
    fq_path = target_path +'/fq'    
    quant_path = target_path+'/quant'
    if not os.path.exists(quant_path):os.mkdir(quant_path)
    transcript_path = target_path+'/rna_assemble'
    if not os.path.exists(transcript_path):os.mkdir(transcript_path)
    map_trpts_path = target_path +'/map_rna_assembl'
    if not os.path.exists(map_trpts_path):os.mkdir(map_trpts_path)
    
    for index in indexes:
        if os.path.exists(fq_path):shutil.rmtree(fq_path)
        os.mkdir(fq_path)
        # 1. copy files
        copy_fq_files(rnaseq_path,fq_path,index)
        # 2. run stringtie function
        cmd = ('python {ppl} {par}').format(ppl=ppl_fn,par=par_fn)
        sarge.run(cmd)
        # 3. move file 
        quant_fns = glob.glob(fq_path+'/stringtie/*.tab')
        for f in quant_fns:
            shutil.move(f,quant_path)
        # 4. move mapped transcripts
        map_assemble_fns = glob.glob(fq_path+'/stringtie/*cov_ref.gtf')
        for f in map_assemble_fns:
            shutil.move(f,map_trpts_path)
        # 5. move assembled transcripts
        assemble_fns = glob.glob(fq_path+'/stringtie/*.gtf')
        for f in assemble_fns:
            shutil.move(f,transcript_path)
        
# run_stringtie_ppl(index_batch,rnaseq_path,target_path,stringtie_ppl_fn,stringtie_par_fn)  
'''

