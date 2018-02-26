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
mzid_folder = ''
fdr = 0.01
for argument in sys.argv:
    if argument.startswith('-'):
        if argument == '-d':
            indicator = 0 #database
        elif argument == '-s':
            indicator = 1 #spectra
        elif argument == '-o':
            indicator = 2
        elif argument == '-mz':
            indicator = 3
        elif argument == '-fdr':
            indicator = 4
        else:
            indicator = 5
            temp = argument
    else:
        if indicator == 0:
            database.append(argument)
        elif indicator == 1:
            spectra = [argument+'/'+x for x in os.listdir(argument)]
        elif indicator == 2:
            main_folder = argument
        elif indicator == 3:
            mzid_folder = argument            
        elif indicator == 4:
            fdr = argument
        elif indicator == 5:
            arguments = arguments + temp + ' '+argument+' '
        



def MultiStageFDR(db_order,folder,mzid_folder):
    os.system('mkdir '+folder)
    my_list = []
    for db in db_order:
        os.system('mkdir '+folder+'/'+db)
        os.system('python '+system_folder+'/MzidToTsv_withmod.py '+mzid_folder+'/'+db+' '+folder+'/'+db+'/'+db+'.tsv')
        last_filename = folder+'/'+db+'/'+db+'.tsv'
#         num = 'f'
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
        os.system('java -Xmx3500M -cp '+system_folder+'/../../MSGFDB/2012.0607/MSGFDB.jar fdr.ComputeFDR -f '+last_filename+' 9 XXX -i 0 -n 2 -p 8 -s 12 0 -fdr '+str(fdr)+' -o '+folder+'/'+db+'.txt')
        my_list.append(folder+'/'+db+'.txt')
    return my_list

system_folder = os.path.split(sys.argv[0])[0]


index = -1
db_library = []
db_order = []
database = database[::-1]
for db in database[::-1]:
    filelist = [f for f in os.listdir(db) if os.path.isfile(os.path.join(db,f))]
    if filelist == []:
        database.remove(db)
        continue
    print "DB: ",db
    db_order.append(db)


result_list = MultiStageFDR(db_order,main_folder,mzid_folder)
s = open(os.path.split(result_list[0])[0]+'/Novel.txt','w')
for result in result_list:
    check_double = []
    if result == result_list[0]:
        os.system('mv '+result+' '+os.path.split(result)[0]+'/Known.txt')
    else:
        f = open(result,'r')
        for line in f:
            if line.startswith('#') and result == result_list[1]:
                s.write(line)
                continue
            elif line.startswith('#'):
                continue
            data = line.split('\t')
            identity = data[0]+data[1]
            if identity in check_double:
                continue
            check_double.append(identity)
            s.write(line)
            
print db_library
    