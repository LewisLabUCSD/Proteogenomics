'''
Created on 2013. 3. 19.

@author: Akard3
'''
import time
import re
import cPickle as pickle
import os
import sys

# print 'Begins: ',time.ctime()

case = 0

if case == 0:
    fasta_database_name = sys.argv[1]
    output_file_name = sys.argv[2]
    pickle_file_name = sys.argv[3]
    ntt = sys.argv[4]
#if len(sys.argv) > 3:
#        tmp_file_name = sys.argv[3]
else:
    #spectrum_file_name = 'CPTAC_BRCA_VGraph_TCGA_VCF_750_A0201_WashU-Pool1_Aliquot7_Process1_20120203_r_klc_x_A1_t1_fr01_VCF_chr1_30_2.tsv'
#     fasta_database_name = 'RNAseqMatch/v1-23_chr1.fa'
#     output_file_name = 'RNAseqMatch/output_v1-23_chr1.txt'
#     tmp_file_name = 'RNAseqMatch/v1-23_chr1.tmp'

    fasta_database_name = '/data/s3cha/CHO_ENOSI_JOB/e22076b9e7ab40ef8fc2c63eddf67eba/sequence3/splice.fa'
    output_file_name = '/data/s3cha/CHO_ENOSI_JOB/e22076b9e7ab40ef8fc2c63eddf67eba/group1/event_out/event_location.txt'
    pickle_file_name = '/data/s3cha/CHO_ENOSI_JOB/e22076b9e7ab40ef8fc2c63eddf67eba/group1/event_out/event_pepdic.p'
    ntt = '2'


#spectrum_file = open(spectrum_file_name,'r')
# print fasta_database_name
#fasta_database = open(fasta_database_name,'r')
output = open(output_file_name,'w')
pickle_file = pickle_file_name
# print fasta_database_name
pepdic = pickle.load( open(pickle_file,"rb") )
db_list = fasta_database_name.split(',')
#############################################################
#functions
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
        

def    CleanPeptideString(pep_str):
#     if pep_str[1] == '.':
#         pep_str = pep_str[2:-2]
    pep_str = pep_str.replace("0","")
    pep_str = pep_str.replace("1","")
    pep_str = pep_str.replace("2","")
    pep_str = pep_str.replace("3","")
    pep_str = pep_str.replace("4","")
    pep_str = pep_str.replace("5","")
    pep_str = pep_str.replace("6","")
    pep_str = pep_str.replace("7","")
    pep_str = pep_str.replace("8","")
    pep_str = pep_str.replace("9","")
    pep_str = pep_str.replace("+","")
    pep_str = pep_str.replace(".","")
    pep_str = pep_str.replace("?","_")
    pep_str = pep_str.replace("_","")
    pep_str = pep_str.replace("(","")
    pep_str = pep_str.replace(")","")
    pep_str = pep_str.replace("[","")
    pep_str = pep_str.replace("]","")
    pep_str = pep_str.replace(",","")
    return pep_str

def get_location(start_seq,end_seq,start,end,length,strand,chrNum):
    i = -1
    j = -1
    while start_seq >= 0:
        i += 1
        start_seq -= length[i]
    while end_seq > 0:
        j += 1
        end_seq -= length[j]
    if strand == 0:
        tmp_start = start[i] - start_seq
        tmp_end = start[j] - end_seq
    else:
        tmp_start = end[i] + start_seq
        tmp_end = end[j] + end_seq
    #print ' start,end,i,j ',tmp_start,tmp_end,i,j,
    location = [[],[],[]]
    if strand == 0:
        location.append(tmp_end)
        location.append(tmp_start)
        if j-i <= 0:
            location[0].append(tmp_end)
            location[1].append(tmp_start)
            location[2].append(tmp_start-tmp_end)
        else:
            location[0].append(start[i])
            location[1].append(tmp_start)
            location[2].append(tmp_start-start[i])
            for k in range(j-i-1):
                location[0].append(start[k+i+1])
                location[1].append(end[k+i+1])
                location[2].append(end[k+i+1]-start[k+i+1])
            location[0].append(tmp_end)
            location[1].append(end[j])
            location[2].append(end[j]-tmp_end)
    else:
        location.append(tmp_start)
        location.append(tmp_end)
        if j-i <= 0:
            location[0].append(tmp_start)
            location[1].append(tmp_end)
            location[2].append(tmp_end-tmp_start)
        else:
            location[0].append(tmp_start)
            location[1].append(end[i])
            location[2].append(end[i]-tmp_start)
            for k in range(j-i-1):
                location[0].append(start[k+i+1])
                location[1].append(end[k+i+1])
                location[2].append(end[k+i+1]-start[k+i+1])
            location[0].append(start[j])
            location[1].append(tmp_end)
            location[2].append(tmp_end-start[j])
    '''
    location = ''
    if strand == 0:
        if j-i <= 0:
            location = str(tmp_end) + '-' + str(tmp_start)
        else:
            location += str(start[i])+'-'+str(tmp_start)+';'
            for k in range(j-i-1):
                location += str(start[k+i]) +'-' + str(end[k+i]) +';'
            location += str(tmp_end) + '-'+str(end[j])
    else:
        if j-i <= 0:
            location = str(tmp_start) + '-' + str(tmp_end)
        else:
            location = str(tmp_start)+'-'+str(end[i])+';'
            for k in range(j-i-1):
                location += str(start[k+i]) +'-' + str(end[k+i]) +';'
            location += str(start[j]) + '-'+str(tmp_end)
    '''
    location.append(strand)
    location.append(chrNum)
    return location

# replace dictioanry key
pep_list = pepdic.keys()
for pep in pep_list:
    key = pep
    pep = pepdic[pep][0][0].split('\t')[8]
    pep = CleanPeptideString(pep)
    if ntt == '2':
        if pep[0] == '-':
            pep = pep[1:-1]
        else:
            pep = pep[:-1]
    elif ntt == '1':
        if pep[-2] in ['R','K']:
            pep = pep[1:-1]
        elif pep[0] in ['R','K']:
            pep = pep[:-1]
        else:
            raise ValueError('protein is not tryptic peptides')
    pepdic[pep] = pepdic.pop(key)
##data base read######################################################
for db in db_list:
    print 'Reading database: ',db
    fasta_database = open(db,'r')
    dummy, fileExtension = os.path.splitext(fasta_database_name)
    
    if fileExtension.strip() == '.txt':
        file_case = 0
    elif fileExtension.strip() == '.fa':
        file_case = 1
    elif fileExtension.strip() == '.fasta':
        file_case = 1
        
    if file_case == 1:
        fasta = fasta_database.readlines()
        fasta_info = []
        index_info = [0]
        fasta_seq = ''
        #sequence_info = []
        
        count = 0
        for i in range(len(fasta)):
            if count % 2 == 0:
                fasta_info.append(fasta[i])
            else:
                #sequence_info.append(fasta[i].strip())
                fasta_seq += fasta[i].strip() + 'X'
                index_info.append(len(fasta_seq))
            count += 1
            
        del fasta
    else:
        fasta_info = []
        index_info = [0]
        fasta_seq = ''
        for file in fasta_database:
            input = open(file.strip(),'r')
            fasta = input.readlines()
            count = 0
            for i in range(len(fasta)):
                if count % 2 == 0:
                    fasta_info.append(fasta[i])
                else:
                    #sequence_info.append(fasta[i].strip())
                    fasta_seq += fasta[i].strip() + 'X'
                    index_info.append(len(fasta_seq))
                count += 1
                
            del fasta
    #>Splice@chr1@176@15444@0;20082126-20082213;20073655-20073753;20072948-20073092;20072024-20072144;20070131-20070219;  example of splice_info
    
    pep_list = pepdic.keys()
    for pep in pep_list:
        try:
            pep = CleanPeptideString(pep)
        except:
            print pep
            assert False,'is wrong'
        location_index = [m.start() for m in re.finditer(pep,fasta_seq)]
        pep_location = []
        for location in location_index:
            index = binary_search(index_info,location,0,len(index_info)-1)
            splice_info = fasta_info[index].split('@')
            #full_seq = sequence_info[index]
            #start_seq = full_seq.find(pep)
            start_seq = location - index_info[index]
            end_seq = start_seq + len(pep)
            start_seq = start_seq *3
            end_seq = end_seq *3
            if splice_info[0].find('Splice')>-1 or splice_info[0].find('Varient')>-1 or splice_info[0].find('rs')>-1:
                chrNum = splice_info[1]
                splice_in = splice_info[4].split(';')
                #print start_seq,end_seq, sequence, full_seq,
                start = []
                end = []
                length = []
                for i,splice in enumerate(splice_in):
                    if i == 0:
                        strand = int(splice)
                    elif i == len(splice_in)-1:
                        continue
                    else:
                        if splice.find('/')>-1:
                            splice = splice.split('/')
                        else:
                            splice = splice.split('-')
                        start.append(int(splice[0]))
                        end.append(int(splice[1]))
                        length.append(int(splice[1])-int(splice[0]))
            else:
                if len(splice_info)<4:
                    continue
                strand = int(splice_info[3])
                start = [int(splice_info[1])]
                end = [int(splice_info[2])]
                length = [end[0]-start[0]]
                chrNum = splice_info[0][1:]
                        
            pep_location.append(get_location(start_seq,end_seq,start,end,length,strand,chrNum))
        
        prev = 0
        i = 0
        while i < len(pep_location):  # delete same location seq
            if prev == pep_location[i][3]:
                del pep_location[i]
            else:
                prev = pep_location[i][3]
                i += 1
        
        ### save pep location in pepdic with checking the occurence in previous pepdic
        for loc in pep_location:
            chrNum = loc[6]#[-1]
            if not pepdic[pep][1].has_key(chrNum):
                pepdic[pep][1][chrNum] = [loc]
            else:
                list_in_pepdic = pepdic[pep][1].get(chrNum)
                is_in_list = 0
                for i in list_in_pepdic:
                    if i[3] == loc[3]:
                        is_in_list = 1
                if is_in_list == 0:
                    pepdic[pep][1][chrNum].append(loc)
                    
        

    '''
    for i in pep_location:
        output.write(pep+'\t')
        for j in range(len(i[0])):
            output.write(str(i[0][j])+'/'+str(i[1][j])+';')
        output.write('\n')


print pepdic[pepdic.keys()[0]]
print pepdic['RRRRRKRK']

for i in pepdic:
    if len(pepdic[i][1]) == 0:
        continue
    for j in pepdic[i][1].values():
        for k in j:
            output2.write(i+'\t')
            for l in range(len(k[0])):
                output2.write(str(k[0][l])+'/'+str(k[1][l])+';')
            output2.write('\n')
'''
del fasta_info
del index_info
del fasta_seq


####################################################################### grouping

#location count fill in
for pep in pepdic:
    count = 0
    for chr in pepdic[pep][1]:
        count += len(pepdic[pep][1][chr])
    if len(pepdic[pep]) == 5:
        pepdic[pep].append(count)
    else:
        pepdic[pep][5] = count


#copy list
coordi_list = {}
for pep in pepdic:
    for chrNum in pepdic[pep][1]:
        if not chrNum in coordi_list.keys():
            coordi_list[chrNum] = {}
        for list in pepdic[pep][1][chrNum]:
            coordi_list[chrNum][list[3]] = pep
#coordi_list = {'chr1':{100000:'PEPTIDE',1000100:'PEPDIDE'}}
#list_in_chr = [1000000,10000100,...] coordinates of peptide in each chromosome
#grouping_inf = [0,1,3,6,10,...] index of beginning groups
output.write('#chr\tstart-end\tPEP\tspec_count\tlocation_count\tFDR\tSprob\n')


for chrNum in coordi_list:
    list_in_chr = coordi_list[chrNum].keys()
    list_in_chr.sort()
    grouping_info = []
    gstart = -100
    for index in range(len(list_in_chr)):
        if gstart + 100 < list_in_chr[index]:
            grouping_info.append(index)
        gstart = list_in_chr[index]

    Sprob = []
    for group in range(len(grouping_info)-1):
        prob = 1
        for index in range(grouping_info[group],grouping_info[group+1]):
            prob *= (1- (1-pepdic[coordi_list[chrNum][list_in_chr[index]]][4]) / pepdic[coordi_list[chrNum][list_in_chr[index]]][5])
        prob = 1- prob
        Sprob.append(prob)    
    prob = 1
    for index in range(grouping_info[len(grouping_info)-1],len(coordi_list[chrNum])):
        prob *= (1- (1-pepdic[coordi_list[chrNum][list_in_chr[index]]][4])/ pepdic[coordi_list[chrNum][list_in_chr[index]]][5])
    prob = 1- prob    
    Sprob.append(prob)


    for coord in list_in_chr:
        Sp_index = binary_search(grouping_info,list_in_chr.index(coord),0,len(grouping_info)-1)
        pep = coordi_list[chrNum][coord]
        output.write(chrNum+'\t')
        for i in pepdic[pep][1].get(chrNum):
            if i[3] != coord:
                continue
            for j in range(len(i[0])):
                output.write(str(i[0][j])+'/'+str(i[1][j])+';')
            output.write('\t'+pep+'\t'+str(pepdic[pep][3])+'\t'+str(pepdic[pep][5])+'\t')
            output.write(str(pepdic[pep][4])+'\t')
            #output.write(str(Sprob[Sp_index])+'\n')
            output.write(str(Sprob[Sp_index])+'\t'+str(i[-2])+'\n')






#######################################################################
# print 'End: ',time.ctime()
