'''
Created on Jan 21, 2015

@author: s3cha
'''

import os
import time
import sys

indicator = -1
database = []
spectra = []
temp = ''
arguments = ''
dbout = ''
fdr = 0.01
for argument in sys.argv:
    if argument.startswith('-'):
        if argument == '-d':
            indicator = 0 #database
        elif argument == '-s':
            indicator = 1 #spectra
        elif argument == '-o':
            indicator = 2
        elif argument == '-fdr':
            indicator = 3
        elif argument == '-db':
            indicator = 4
        else:
            indicator = 5
            temp = argument
    else:
        if indicator == 0:
            database.append(argument)
        elif indicator == 1:
            spectra += [argument+'/'+x for x in os.listdir(argument) if os.path.splitext(x)[1] != '.gz']
        elif indicator == 2:
            main_folder = argument
        elif indicator == 3:
            fdr = argument
        elif indicator == 4:
            dbout = argument
        elif indicator == 5:
            arguments = arguments + temp + ' '+argument+' '
            
    
        

# database = ['/home/s3cha/data/demo_del/demo','/home/s3cha/data/demo_del/PNNL']
# database = ['demo','PNNL']
# database = ['demo','abc','abcd']
# spectra = ['/home/s3cha/data/ProteoSaFe/ProteoSAFe-1.2.6_beta-linux64/data/uploads/default/test/foo.mgf','/home/s3cha/data/VU/Spectrum/mgf/TCGA-AA-A017-01A-22_W_VU_20120817_A0218_2B_R_FR02.mgf']
# arguments = '-t 30ppm -tda 1 -m 0 -inst 0 -e 1 -nnet 1 -thread 2 -showFDR 0 -replicate 1'
# arguments = '-t 30ppm -tda 1 -m 0 -e 1 -thread 2 -ntt 1'
# fdr = 0.01
# main_folder = 'multistageFDR'
print sys.argv
print arguments
print main_folder
print database
print spectra



max_db_size = 100000000

def SplitBySize(size,myfile,folder):
    filesize = os.path.getsize(os.path.join(folder+'/sum/',myfile))
    if filesize < size:
        os.system('cp '+str(folder.rstrip('/')+'/sum/'+myfile)+' '+folder.rstrip('/')+'/split')
        os.system('rm '+str(folder.rstrip('/')+'/sum/*'))
        return
    new_size = filesize/(filesize/size+1)
    f = open(os.path.join(folder+'/sum/',myfile),'r')
    index = 1
    out_folder = folder+'/split/'
    filename = os.path.splitext(myfile)[0]
    out_filename = out_folder+filename+'_'+str(index)+'.fa'
    s = open(out_filename,'w')
    count = 1
    for line in f:
        if line.startswith('>'):
            if os.path.getsize(out_filename)>new_size:
                index += 1
                s.close()
                out_filename = out_folder+filename+'_'+str(index)+'.fa'
                s = open(out_filename,'w')
            s.write('>'+str(count)+'\n')
            count += 1
        else:
            s.write(line)
    os.system('rm '+str(folder.rstrip('/')+'/sum/*'))
    return

def BuildSA(db_piece,db):
#     print ':::::::mkdir '+db+'/prep/'+os.path.splitext(db_piece)[0]
#     print ':::::::mv '+db+'/split/'+db_piece+' '+db+'/prep/'+os.path.splitext(db_piece)[0]
#     print ':::::::java -Xmx3500M -cp '+system_folder+'/MSGFPlus.20130403/MSGFPlus.jar edu.ucsd.msjava.msdbsearch.BuildSA -d '+db+'/prep/'+os.path.splitext(db_piece)[0]+'/'+db_piece+' -tda 2'
    os.system('mkdir '+db+'/prep/'+os.path.splitext(db_piece)[0])
    os.system('mv '+db+'/split/'+db_piece+' '+db+'/prep/'+os.path.splitext(db_piece)[0])
    os.system('java -Xmx10000M -cp '+system_folder+'/MSGFPlus.20130403/MSGFPlus.jar edu.ucsd.msjava.msdbsearch.BuildSA -d '+db+'/prep/'+os.path.splitext(db_piece)[0]+'/'+db_piece+' -tda 2')
    return

def MultiStageFDR(db_order,folder):
    os.system('mkdir '+folder)
    my_list = []
    for db in db_order:
        os.system('mkdir '+folder+'/'+db)
        os.system('python '+system_folder+'/MzidToTsv_withmod.py '+db+'/result '+folder+'/'+db+'/'+db+'.tsv')
        last_filename = folder+'/'+db+'/'+db+'.tsv'
        num = 'f'
#         print 4444
#         print my_list
        for item in my_list:
            temp = os.path.splitext(os.path.split(item)[1])[0]
#             print 5555
            new_filename = os.path.splitext(last_filename)[0]+'_'+temp+'_filtered'+'.tsv'
            os.system('python '+system_folder+'/FilterOutRefSeqIdentifications.py '+item+' '+last_filename+' '+new_filename)
            last_filename = new_filename
#         print 33333,last_filename
#         print 'java -Xmx3500M -cp '+system_folder+'/../../MSGFDB/2012.0607/MSGFDB.jar fdr.ComputeFDR -f '+last_filename+' 9 XXX -i 0 -n 2 -p 8 -s 12 0 -fdr '+str(fdr)+' -o '+folder+'/'+db+'.txt'
        os.system('java -Xmx1000M -cp '+system_folder+'/../../MSGFDB/2012.0607/MSGFDB.jar fdr.ComputeFDR -f '+last_filename+' 9 XXX -i 0 -n 2 -p 8 -s 12 0 -fdr '+str(fdr)+' -o '+folder+'/'+db+'.txt')
        my_list.append(folder+'/'+db+'.txt')
    return my_list

system_folder = os.path.split(sys.argv[0])[0]


index = -1
db_library = []
db_order = []
database = database[::-1]
os.system('mkdir '+dbout)
if not os.path.exists(main_folder):
    os.mkdir(main_folder)
for db in database[::-1]:
    filelist = [f for f in os.listdir(db) if os.path.isfile(os.path.join(db,f))]
    if filelist == []:
        database.remove(db)
        continue
    target_db = dbout+'/'+db
    db_order.append(target_db)
    index += 1
    db_library.append([])
    os.system('mkdir '+target_db)
    os.system('mkdir '+target_db+'/sum')
    os.system('mkdir '+target_db+'/split')
    os.system('mkdir '+target_db+'/prep')
    os.system('mkdir '+target_db+'/result')
    filename = target_db.split('/')[-1]+'.fa'
    if len(filelist)>1:
        os.system('java -cp '+system_folder+'/../../scripts/2012.0815/CCMSWorkflowUtils.jar edu.ucsd.workflow.Merge -input '+db+' -output '+target_db+'/sum/'+filename)
    else:
        os.system('cp '+db+'/'+filelist[0]+' '+target_db+'/sum/'+filename)
    SplitBySize(max_db_size,os.listdir(target_db+'/sum')[0],target_db)
    my_list = os.listdir(target_db+'/split')
    for db_piece in my_list:
        db_library[index].append(target_db+'/prep/'+os.path.splitext(db_piece)[0]+'/'+db_piece)
        BuildSA(db_piece,target_db)
        print db_piece,target_db



execution_string = []
for index,db in enumerate(db_library):
    for db_piece in db:
        for spectrum in spectra:
            #output = db_order[index]+'/result/'+os.path.splitext(os.path.split(spectrum)[1])[0]+'_'+os.path.splitext(os.path.split(db_piece)[1])[0]+'.mzid'
            output = 'MSGFresult/'+db_order[index].split('/')[-1]+'/'+os.path.splitext(os.path.split(spectrum)[1])[0]+'_'+os.path.splitext(os.path.split(db_piece)[1])[0]+'.mzid'
            execution_string.append('java -Xmx10G -jar '+system_folder+'/MSGFPlus.20130403/MSGFPlus.jar -s '+spectrum+' -d '+db_piece+' -o '+output+' '+arguments)
#             s_exe = open(main_folder+'/'+os.path.splitext(os.path.split(output)[1])[0]+'.txt','w')
#             s_exe.write('java -Xmx3500M -jar '+system_folder+'/MSGFPlus.20130403/MSGFPlus.jar -s '+spectrum+' -d '+db_piece+' -o '+output+' '+arguments)
#             s_exe.write('java -Xmx3500M -jar '+system_folder+'/../../MSGFPlus/2012.0830/MSGFPlus.jar -s '+spectrum+' -d '+db_piece+' -o '+output+' '+arguments)
#            s_exe.write('java -Xmx3500M -jar '+system_folder+'/MSGFPlus.20130403/MSGFPlus.jar -s '+spectrum+' -d '+db_piece+' -o '+output+' '+arguments)
#             os.system('java -Xmx3500M -jar '+system_folder+'/MSGFPlus.20130403/MSGFPlus.jar -s '+spectrum+' -d '+db_piece+' -o '+output+' '+arguments)
#             s_exe.close()
num_pal = 40
if len(execution_string) < num_pal:
    for index,string in enumerate(execution_string):
        s = open(main_folder+'/'+'LaunchMSGF_'+str(index)+'.txt','w')
        s.write(string)
else:
    for index in range(num_pal):
        s = open(main_folder+'/'+'LaunchMSGF_'+str(index)+'.txt','w')
        for coeffi in range(len(execution_string)/num_pal):
            s.write(execution_string[coeffi*num_pal+index]+'\n')
        if index < len(execution_string)%num_pal:
            s.write(execution_string[len(execution_string)/num_pal*num_pal+index]+'\n')






