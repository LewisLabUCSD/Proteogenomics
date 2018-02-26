'''
Created on Nov 8, 2013

@author: s3cha
'''
import time
import sys
import os
# need 3 input. unknown locations, known locations, ref-seq gff

btime = time.clock()
# print 'Begins: ', time.ctime()
DNA_folder = ''
if len(sys.argv) > 1:
    unknown_location_filename = sys.argv[1]
    known_location_filename = sys.argv[2]
    ref_seq_filename = sys.argv[3]
    output_filename = sys.argv[4]
    if len(sys.argv)>5:
        DNA_folder = sys.argv[5]
else:
    #VU run
    unknown_location_filename = '/home/s3cha/data/VU/Cosmic/Location_VU_SPECFDR01_Cosmicoverlap2.txt'#'/home/s3cha/data/VU/Cosmic/Location_VU_SPECFDR01_Cosmicoverlap2.txt'#'/home/s3cha/data/VU/Somatic/Location_VU_RefSeq_PEPFDR01_rna.txt'#'/home/s3cha/data/PNNL/ov_somatic_mutations/Locations_PNNL_20140202_removed_dbSNP_somatic.txt'#'/home/s3cha/data/VU/Somatic/Location_VU_SPECFDR01_splice_vcf_dbSNP_somatic.txt'#'/home/s3cha/data/JHU/Location_JHU_SV_Result_SVrna_dbSNP_somatic.txt'#'/home/s3cha/data/PNNL/Location_PNNL_Single_recalculated.txt'#'C:/Users/Akard3/Dropbox/Workspace/SpliceGraph/src/Enosi/Locations_PNNL_20140202_removed_dbSNP.txt'#Location_PNNL_Single_recalculated.txt'#sample.txt'#
    known_location_filename = 'Location_Combine123_RefSeq_recalculated_Correction.txt'#'/home/s3cha/data/VU/Location_VU_RefSeq_PEPFDR01.txt'#
    ref_seq_filename = 'Homo_sapiens.GRCh37.70_formatted.gff'
    output_filename = '/home/s3cha/data/VU/Cosmic/Event3_VU_somatic_cosmic.txt'#'/home/s3cha/data/PNNL/ov_somatic_mutations/EventLocations_PNNL_20140202_removed_dbSNP_somatic.txt'#'/home/s3cha/data/VU/Somatic/Event_VU_SPECFDR01_splice_vcf_dbSNP_somatic.txt'#'/home/s3cha/data/PNNL/Event_PNNL_Single_recalculated.txt'#'/home/s3cha/Dropbox/Workspace/SpliceGraph/src/Enosi/sample_event2.txt'#Event_PNNL_20140202_Concat.txt'
    unknown_location_filename = '/home/s3cha/data/VU/Normal/Location_Normal_colon_SV_merged_somatic.txt'#'/home/s3cha/data/VU/Cosmic/Location_VU_SPECFDR01_Cosmicoverlap2.txt'#'/home/s3cha/data/VU/Somatic/Location_VU_RefSeq_PEPFDR01_rna.txt'#'/home/s3cha/data/PNNL/ov_somatic_mutations/Locations_PNNL_20140202_removed_dbSNP_somatic.txt'#'/home/s3cha/data/VU/Somatic/Location_VU_SPECFDR01_splice_vcf_dbSNP_somatic.txt'#'/home/s3cha/data/JHU/Location_JHU_SV_Result_SVrna_dbSNP_somatic.txt'#'/home/s3cha/data/PNNL/Location_PNNL_Single_recalculated.txt'#'C:/Users/Akard3/Dropbox/Workspace/SpliceGraph/src/Enosi/Locations_PNNL_20140202_removed_dbSNP.txt'#Location_PNNL_Single_recalculated.txt'#sample.txt'#
    known_location_filename = 'Location_Combine123_RefSeq_recalculated_Correction.txt'#'/home/s3cha/data/VU/Location_VU_RefSeq_PEPFDR01.txt'#
    ref_seq_filename = 'Homo_sapiens.GRCh37.70_formatted.gff'
    output_filename = '/home/s3cha/data/VU/Normal/Event_Normal_colon_SV_merged_somatic.txt'#'/home/s3cha/data/PNNL/ov_somatic_mutations/EventLocations_PNNL_20140202_removed_dbSNP_somatic.txt'#'/home/s3cha/data/VU/Somatic/Event_VU_SPECFDR01_splice_vcf_dbSNP_somatic.txt'#'/home/s3cha/data/PNNL/Event_PNNL_Single_recalculated.txt'#'/home/s3cha/Dropbox/Workspace/SpliceGraph/src/Enosi/sample_event2.txt'#Event_PNNL_20140202_Concat.txt'
    
    DNA_folder = '/home/s3cha/s3cha/data/dna/70'
#     DNA_folder = 'D:/workspace/DNA/GRCh37.70'
    
    
#     unknown_location_filename = '/home/s3cha/data/VU/Comet/temp/Location_VU_DT_Comet.txt'#'/home/s3cha/data/VU/Somatic/Location_VU_RefSeq_PEPFDR01_rna.txt'#'/home/s3cha/data/PNNL/ov_somatic_mutations/Locations_PNNL_20140202_removed_dbSNP_somatic.txt'#'/home/s3cha/data/VU/Somatic/Location_VU_SPECFDR01_splice_vcf_dbSNP_somatic.txt'#'/home/s3cha/data/JHU/Location_JHU_SV_Result_SVrna_dbSNP_somatic.txt'#'/home/s3cha/data/PNNL/Location_PNNL_Single_recalculated.txt'#'C:/Users/Akard3/Dropbox/Workspace/SpliceGraph/src/Enosi/Locations_PNNL_20140202_removed_dbSNP.txt'#Location_PNNL_Single_recalculated.txt'#sample.txt'#
#     known_location_filename = 'Location_Combine123_RefSeq_recalculated_Correction.txt'#'/home/s3cha/data/VU/Location_VU_RefSeq_PEPFDR01.txt'#
#     ref_seq_filename = 'Homo_sapiens.GRCh37.70_formatted.gff'
#     output_filename = '/home/s3cha/data/VU/Comet/temp/Event3_VU_DT_Comet.txt'#'/home/s3cha/data/PNNL/ov_somatic_mutations/EventLocations_PNNL_20140202_removed_dbSNP_somatic.txt'#'/home/s3cha/data/VU/Somatic/Event_VU_SPECFDR01_splice_vcf_dbSNP_somatic.txt'#'/home/s3cha/data/PNNL/Event_PNNL_Single_recalculated.txt'#'/home/s3cha/Dropbox/Workspace/SpliceGraph/src/Enosi/sample_event2.txt'#Event_PNNL_20140202_Concat.txt'
    
    
    
    
    #PNNL run
    #unknown_location_filename = '/home/s3cha/data/PNNL/Result/Location_Combine_PNNL_JHU.txt'
    #known_location_filename = 'Location_Combine123_RefSeq_recalculated_Correction.txt'#'/home/s3cha/data/VU/Location_VU_RefSeq_PEPFDR01.txt'#
    #ref_seq_filename = 'Homo_sapiens.GRCh37.70_formatted.gff'
    #output_filename = '/home/s3cha/data/PNNL/Result/Event3_Combine_PNNL_JHU.txt'#'/home/s3cha/data/PNNL/ov_somatic_mutations/EventLocations_PNNL_20140202_removed_dbSNP_somatic.txt'#'/home/s3cha/data/VU/Somatic/Event_VU_SPECFDR01_splice_vcf_dbSNP_somatic.txt'#'/home/s3cha/data/PNNL/Event_PNNL_Single_recalculated.txt'#'/home/s3cha/Dropbox/Workspace/SpliceGraph/src/Enosi/sample_event2.txt'#Event_PNNL_20140202_Concat.txt'
    
    #unknown_location_filename = '/home/s3cha/Dropbox/Workspace/SpliceGraph/src/Enosi/Locations_PNNL_20140202_removed_dbSNP.txt'#sample.txt'#Locations_PNNL_Concat_FDR_20140210_dbSNP.txt'#
    #known_location_filename = 'Location_Combine123_RefSeq_recalculated.txt'
    #ref_seq_filename = 'Homo_sapiens.GRCh37.70_formatted.gff'
    
DNA_folder = DNA_folder.strip()
filter = True
# filter = False
if DNA_folder != '':
    dna_file_list = os.listdir(DNA_folder)
    for file in dna_file_list:
    
        if os.path.splitext(file)[1] == '.trie' and file.find('chr10')>-1:
            file = file.split('chr10')
            dna_string1 = file[0]
            dna_string2 = file[1]
            break

'''
if os.path.exists(output_filename):
    usrin = raw_input('Output file is already exist. Do you want to overwrite? (y/n): ')
    if usrin == 'y' or usrin == 'Y':
        s = open(output_filename,'w')
else:
    s = open(output_filename,'w')
'''  
s = open(output_filename, 'w')
header = '#Num\tEvent\tPeptide\tChromosome\tLocation\tNumOfNovel\tNumOfKnown\tspec_count\tlocation_count\tFDR\tSprob\tStrand\tGbrowserLink\tGroupInfo\tRelatedGene\n'
s.write(header)
s2 = open(os.path.split(output_filename)[0]+'/event_num.txt','w')
s2.write('#Event\tCount\n')
log = open(os.path.splitdrive(output_filename)[0] + 'log.ini', 'w')

db_str = 'hg19'
usa_gffile_str = 'Locations_PNNL_20140210_UCSC.gff'
usa_gffile_str = 'Locations_VU_20140428_UCSC_SPEC.gff'
usa_gffile_str = 'Locations_JHU_20140523_UCSC_SPEC.gff'
#usa_gffile_str = 'Locations_VU_20140428_UCSC_PEP.gff'
                  
ForwardCode = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
               "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
               "TAT":"Y", "TAC":"Y", "TAA":"X", "TAG":"X",
               "TGT":"C", "TGC":"C", "TGA":"X", "TGG":"W",
               "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
               "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
               "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
               "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
               "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
               "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
               "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
               "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
               "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
               "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
               "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
               "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
               }
ReverseCode =   {"AAA":"F", "AAG":"F", "AAT":"L", "AAC":"L",
               "AGA":"S", "AGG":"S", "AGT":"S", "AGC":"S",
               "ATA":"Y", "ATG":"Y", "ATT":"X", "ATC":"X",
               "ACA":"C", "ACG":"C", "ACT":"X", "ACC":"W",
               "GAA":"L", "GAG":"L", "GAT":"L", "GAC":"L",
               "GGA":"P", "GGG":"P", "GGT":"P", "GGC":"P",
               "GTA":"H", "GTG":"H", "GTT":"Q", "GTC":"Q",
               "GCA":"R", "GCG":"R", "GCT":"R", "GCC":"R",
               "TAA":"I", "TAG":"I", "TAT":"I", "TAC":"M",
               "TGA":"T", "TGG":"T", "TGT":"T", "TGC":"T",
               "TTA":"N", "TTG":"N", "TTT":"K", "TTC":"K",
               "TCA":"S", "TCG":"S", "TCT":"R", "TCC":"R",
               "CAA":"V", "CAG":"V", "CAT":"V", "CAC":"V",
               "CGA":"A", "CGG":"A", "CGT":"A", "CGC":"A",
               "CTA":"D", "CTG":"D", "CTT":"E", "CTC":"E",
               "CCA":"G", "CCG":"G", "CCT":"G", "CCC":"G",
               }

def CompareLocation(FS, FE, SS, SE):
    case = 0
    if FE < SS or FS > SE:
        case = 10
    elif FS > SS:
        if FE > SE:
            case = 1
        elif FE < SE:
            case = 2
        elif FE == SE:
            case = 3
    elif FS < SS:
        if FE > SE:
            case = 4
        elif FE < SE:
            case = 5
        elif FE == SE:
            case = 6
    elif FS == SS:
        if FE > SE:
            case = 7
        elif FE < SE:
            case = 8
        elif FE == SE:
            case = 9
    
    return case

def FindEventOld(event, each_gene, loc_info):
    case = CompareLocation(each_gene[0], each_gene[1], loc_info[3], loc_info[4])
    if case == 10:
        print 'error'
    if (each_gene[3][6] == "+" and loc_info[2] == 0) or (each_gene[3][6] == "1" and loc_info[2] == 0):
        event[12] = 1
    elif (each_gene[3][6] == "-" and loc_info[2] == 1) or (each_gene[3][6] == "0" and loc_info[2] == 1):
        event[12] = 1
    elif case in [1, 5]:  # case intersection
        if (case == 1 and (event[4] == 2 or event[5] == 2)) or (case == 5 and (event[4] == 1 or event[5] == 1)):
            event[6] = 1
        else:
            for index in range(len(loc_info[0])):
                if each_gene[5] != []:
                    UTR_case = CompareLocation(each_gene[5][3], each_gene[5][4], loc_info[0][index], loc_info[1][index])
                    if UTR_case != 10:
                        event[5] = 1
                if each_gene[6] != []:
                    UTR_case = CompareLocation(each_gene[6][3], each_gene[6][4], loc_info[0][index], loc_info[1][index])
                    if UTR_case != 10:
                        event[5] = 2  # translated UTR
            if event[5] == 0:  # gene boundary
                if case == 1:
                    event[4] = 1
                else:
                    event[4] = 2
    elif case == 4:
        for piece in range(len(loc_info[0])):  # piece is piece of peptide
            check_novelexon = 0
            for exon in each_gene[4]:
                exon_case = CompareLocation(exon[3], exon[4], loc_info[0][piece], loc_info[1][piece])
                if exon_case == 10:
                    continue
                check_novelexon = 1
                if exon_case in [1, 2, 3, 5, 8]:
                    event[7] = 1  # 'exon boundary'
                elif len(loc_info[0]) > 1:
                    # splicing case
                    if (piece == 0 and exon_case == 6) or (piece == len(loc_info[0]) - 1 and exon_case == 7) or ((piece != 0 and piece != len(loc_info[0]) - 1) and exon_case in [3, 6, 7, 8]):
                        event[8] = 1  # 'alternative'
                else:
                    # 6-frame case
                    if exon_case in [4, 6, 7]:
                        event[9] = 1  # 'frame shift'
            if check_novelexon == 0:
                event[10] = 1  # 'novel exon'
        if len(loc_info[0]) > 1 and event[7] + event[8] + event[10] == 0:
            event[11] = 1  # 'novel splice'
    else:
        event[-1] = 1
    # loc_info = [start,end,strand,pep_start,pep_end]
    # event = [0,0,0,0,0,0,0,0,0,0,0,0,0] 
    # 0:mutation , 1:insertion, 2:deletion, 3:novel gene, 4:gene boundary, 5:tranlated UTR, 
    # 6:fusion gene, 7: exon boundary, 8: alternative splice, 9: frame shift, 10: novel exon,
    # 11: novel splice, 12: reverse strand
    return event

def CheckStopCodon(dna_string,strand):
    if strand == 0:
        dna_string = dna_string[::-1]
    for i in range(len(dna_string)/3):
        if strand == 1:
            if ForwardCode.get(dna_string[i*3:(i+1)*3]) == 'X':
                
                return False
        else:
            if ReverseCode.get(dna_string[i*3:(i+1)*3]) == 'X':
                return False
    return True


def compare_start(item1, item2):
    if int(item1[0]) < int(item2[0]):
        return -1
    elif int(item1[0]) > int(item2[0]):
        return 1
    else:
        return 0

def compare_pseudo(item1,item2):
    if item1[0][0] < item2[0][0]:
        return -1
    elif item1[0][0] > item2[0][0]:
        return 1
    else:
        return 0

def sort_location(item1, item2):
    if item1[-1] == 1:
        position1 = item1[0][0]
    else:
        position1 = item1[0][-1]
    if item2[-1] == 1:
        position2 = item2[0][0]
    else:
        position2 = item2[0][-1]
    
    if position1 < position2:
        return -1
    elif position1 > position2:
        return 1
    else:
        return 0

def sort_by_chromosome(item1, item2):
    position1 = item1[1]
    position2 = item2[1]
    
    if position1 < position2:
        return -1
    elif position1 > position2:
        return 1
    else:
        return 0

def WriteLocation(location):
    beg = location[0][0]
    end = location[0][1]
    string = ''
    for i in range(len(beg)):
        if beg[i] < 0:
            if beg[i] < -1000:
                string += ('M' + str(end[i] - beg[i]))
            else:
                string += ('I' + str(end[i] - beg[i]))
        else:
            string += (str(beg[i]) + '-' + str(end[i]))
        if i != len(beg) - 1:
            string += (';')
    return string


def CheckRelation(start, end, gframe, ref_seq):
    exon_case = 10
    for cds in ref_seq[4]:
        if end < cds[3]:
            continue
        elif start > cds[4]:
            continue
        exon_case = CompareLocation(start, end, cds[3], cds[4])
        if cds[7] == '.':
            cds_gframe = 0
        elif ref_seq[3][6] == '+':
            cds_gframe = (cds[3] % 3 + int(cds[7]))%3
        else:
            cds_gframe = (cds[3] % 3 + (cds[4] - cds[3] - int(cds[7])) % 3)%3
        if cds_gframe != gframe:
            exon_case = exon_case - 10
    return exon_case

def Bigger(x, y):
    if x > y:
        return x
    else:
        return y


def Recruit(current_location, unknown_nonsplice, known,unknown_splice,unknown_indel,indicator,strand): 
    known = []
    group = [[current_location], []]
    NearConstant = 1000
    start = current_location[0][0][0]
    end = current_location[0][1][-1]
    start_gframe = start % 3
    end_gframe = (current_location[0][0][-1] - (sum(current_location[0][1][:-1]) - sum(current_location[0][0][:-1]) % 3)) % 3
    
    #if start < 0 or end < 0:
    #    print 'negative start or end error',
    #    print group,start,end
    #    return group
    
    if start < 0:
        start = current_location[0][0][1]
        start_gframe = (start - (current_location[0][1][0] - current_location[0][0][0]))%3
    elif end < 0:
        end = current_location[0][1][-2]
        end_gframe = (current_location[0][0][-1] - (sum(current_location[0][1][:-2]) - sum(current_location[0][0][:-2]) % 3)) % 3
    # grouping the splice event ( if the two different peptide share it's junction, combine them together )
    if indicator == 1 and unknown_splice != None:
        index = 0
        start_c = start
        end_c = end
        while index < len(unknown_splice):
            candidate = unknown_splice[index]
            if candidate[5] != strand:
                index += 1
                continue
            if start_c > candidate[0][1][-1] + NearConstant:
                index += 1
                continue
            if end_c < candidate[0][0][0] - NearConstant:
                break
            isOverlap = True
            check = False
            for i in range(len(current_location[0][0])-1 if len(current_location[0][0])<len(candidate[0][0]) else len(candidate[0][0])-1):
                if current_location[0][1][i] != candidate[0][1][i] or current_location[0][0][i+1] != candidate[0][0][i+1]:
                    isOverlap = False
                else:
                    check = True
            if isOverlap and check and current_location != candidate and start_gframe == candidate[0][0][0] % 3:
                group[0].append(candidate)
                unknown_splice.remove(candidate)
                index -= 1
            else:
                index += 1
    
    elif indicator == 2 and unknown_indel != None:
        index = 0
        start_c = start
        end_c = end
        while index < len(unknown_indel):
            candidate = unknown_indel[index]
            if candidate[5] != strand:
                index += 1
                continue
            if start_c > candidate[0][1][-1] + NearConstant:
                index += 1
                continue
            if end_c < candidate[0][0][0] - NearConstant:
                break
            isOverlap = True
            check = False
            for i in range(len(current_location[0][0])-1 if len(current_location[0][0])<len(candidate[0][0]) else len(candidate[0][0])-1):
                if current_location[0][1][i] != candidate[0][1][i] or current_location[0][0][i+1] != candidate[0][0][i+1]:
                    isOverlap = False
                else:
                    check = True
            if isOverlap and check and current_location != candidate and start_gframe == candidate[0][0][0] % 3:
                group[0].append(candidate)
                unknown_indel.remove(candidate)
                index -= 1
            else:
                index += 1
        
    
    # grouping the event
    # find relative unknown sequence
    index = 0
    start_c = start
    end_c = end
    if unknown_nonsplice != None:
        while index < len(unknown_nonsplice):
            candidate = unknown_nonsplice[index]
            if candidate[5] != strand:
                index += 1
                continue
            candidate_end_gframe = (candidate[0][0][-1] - (sum(candidate[0][1][:-1]) - sum(candidate[0][0][:-1]) % 3)) % 3
            if start_c > candidate[0][1][-1] + NearConstant:
                index += 1
                continue
            if end_c < candidate[0][0][0] - NearConstant:
                break
            if candidate[0][1][-1] < current_location[0][1][0] and start_gframe == candidate[0][0][0] % 3:
                group[0].append(candidate)
                unknown_nonsplice.remove(candidate)
                start_c = candidate[0][0][0]
                index -= 1
            elif candidate[0][0][0] > current_location[0][0][-1] and end_gframe == candidate_end_gframe:
                group[0].append(candidate)
                unknown_nonsplice.remove(candidate)
                end_c = candidate[0][1][-1]
            else:
                index += 1
    
    index = 0
    start_c = start
    end_c = end
    while index < len(known):
        candidate = known[index]
        candidate_end_gframe = (candidate[0][0][-1] - (sum(candidate[0][1][:-1]) - sum(candidate[0][0][:-1]) % 3)) % 3
        if start_c > candidate[0][1][-1] + NearConstant:
            index += 1
            continue
        if end_c < candidate[0][0][0] - NearConstant:
            break
        if candidate[0][1][-1] < current_location[0][1][0] and start_gframe == candidate[0][0][0] % 3:
            group[1].append(candidate)
            # known.remove(candidate)
            # start_c = candidate[0][0][0]
        elif candidate[0][0][0] > current_location[0][0][-1] and end_gframe == candidate_end_gframe:
            group[1].append(candidate)
            # known.remove(candidate)
            # end_c = candidate[0][1][-1]
        index += 1
    
    '''
    itr = iter(unknown_nonsplice)
    while 1:
        try:
            candidate = next(itr)
            candidate_end_gframe = (candidate[0][0][-1]-(sum(candidate[0][1][:-1])-sum(candidate[0][0][:-1])%3))%3
            if start > candidate[0][1][-1] + NearConstant:
                continue
            if end < candidate[0][0][0] - NearConstant:
                break
            if candidate[0][1][-1] < current_location[0][1][0] and start_gframe == candidate[0][0][0]%3:
                'recrute current'
                'delete current from unknown list'
                
            elif candidate[0][0][0] > current_location[0][0][-1] and end_gframe == candidate_end_gframe:
                'recrute current'
                'delete current from unknown list'
        except StopIteration:
            break
    #find relative known sequence
    itr = iter(known)
    while 1:
        try:
            candidate = next(itr)
            candidate_end_gframe = (candidate[0][0][-1]-(sum(candidate[0][1][:-1])-sum(candidate[0][0][:-1])%3))%3
            if start > candidate[0][1][-1] + NearConstant:
                continue
            if end < candidate[0][0][0] - NearConstant:
                break
            if candidate[0][1][-1] < current_location[0][1][0] and start_gframe == candidate[0][0][0]%3:
                'recrute current'
            elif candidate[0][0][0] > current_location[0][0][-1] and end_gframe == candidate_end_gframe:
                'recrute current'
        except StopIteration:
            break
    '''
    return group



def FindSpliceEvent(current_location, ref_seq, event):
    if ref_seq[3][6] == '+':
        ref_strand = 1
    else:
        ref_strand = 0
    if ref_strand != int(current_location[5]):
        event[-1] = Bigger(event[-1], 1)
        # return event########## reverse
    case = CompareLocation(current_location[0][0][0], current_location[0][1][-1], ref_seq[0], ref_seq[1])
    if case == 10:
        return event
    if case == 1 and event[1] != -1: ## check fusion
        gframe = current_location[0][0][0] % 3
        exon_case = CheckRelation(current_location[0][0][0], current_location[0][1][0], gframe, ref_seq)
        if not (exon_case == 10 or exon_case == 0):
            if event[1] == 2 or event[1] == 3:
                event[1] = 3
            else:
                event[1] = 1 #fusion indicator. 1 : left side overlap 2: right side overlap, 3: both side
    elif case == 5 and event[1] != -1:
        gframe = (current_location[0][0][-1] - (sum(current_location[0][1][:-1]) - sum(current_location[0][0][:-1]) % 3)) % 3
        exon_case = CheckRelation(current_location[0][0][-1], current_location[0][1][-1], gframe, ref_seq)
        if not (exon_case == 10 or exon_case == 0):
            if event[1] == 1 or event[1] == 3:
                event[1] = 3
            else:
                event[1] = 2
    if case == 2:
        gframe = current_location[0][0][0] % 3
        exon_case1 = CheckRelation(current_location[0][0][0], current_location[0][1][0], gframe, ref_seq)
        exon_case2 = CheckRelation(current_location[0][0][-1], current_location[0][1][-1], gframe, ref_seq)
        if exon_case1 not in [0,10] or exon_case2 not in [0,10]:
            event[1] = -1
    ref_seq_junction = ref_seq[7]
    UTR_junction = []
    UTR_area = []
    if len(ref_seq[6])>0:
        UTR_area.append([ref_seq[6][0][3],ref_seq[6][0][4]])
    if len(ref_seq[5])>0:
        UTR_area.append([ref_seq[5][0][3],ref_seq[5][0][4]])
    # [start, end, gene_name_str, mRNA, CDS, five_prime_UTR, three_prime_UTR, splice junctions, gene_readable_name_str]
    for index_three in range(len(ref_seq[6])-1):
        UTR_area.append([ref_seq[6][index_three+1][3],ref_seq[6][index_three+1][4]])
        UTR_junction.append([ref_seq[6][index_three][4],ref_seq[6][index_three+1][3]])
    for index_five in range(len(ref_seq[5])-1):
        UTR_area.append([ref_seq[5][index_five+1][3],ref_seq[5][index_five+1][4]])
        UTR_junction.append([ref_seq[5][index_five][4],ref_seq[5][index_five+1][3]])
    
    current_junction = []
    for index in range(len(current_location[0][0])-1):
        current_junction.append([current_location[0][1][index],current_location[0][0][index+1]])
    
    for junction in UTR_junction:
        for c_j in current_junction:
            if junction == c_j:
                event[7] = 2
    if len(current_junction) < 2:
        for index,junction in enumerate(ref_seq_junction):
            if current_junction[0] == junction: # possible case : different frame, exon boundary 
                if ref_strand == 1:
                    if current_location[0][0][0]%3 != (ref_seq[4][index][3]+int(ref_seq[4][index][7]))%3 and event[10] > -1: #check frame shift
                        event[10] = 2
                    elif current_location[0][0][0] < ref_seq[4][index][3] or current_location[0][1][-1] > ref_seq[4][index+1][4]: #check exon boundary
                        for i in UTR_area:
                            if CompareLocation(current_location[0][0][0],current_location[0][1][0],i[0],i[1]) != 10 or CompareLocation(current_location[0][0][-1],current_location[0][1][-1],i[0],i[1]) != 10:
                                event[7] = 2
                            else:
                                event[8] = 2
                        event[10] = -1
                    else:
                        #print 'Splice two boundary error: ',current_location,ref_seq
                        event[10]=0                      
                else:
                    if current_location[0][1][-1]%3 != (ref_seq[4][index][4]-int(ref_seq[4][index][7]))%3 and event[10] > -1:
                        event[10] = 2
                    elif current_location[0][0][0] < ref_seq[4][index+1][3] or current_location[0][1][-1] > ref_seq[4][index][4]: #check exon boundary
                        for i in UTR_area:
                            if CompareLocation(current_location[0][0][0],current_location[0][1][0],i[0],i[1]) != 10 or CompareLocation(current_location[0][0][-1],current_location[0][1][-1],i[0],i[1]) != 10:
                                event[7] = 2
                            else:
                                event[8] = 2
                        event[10] = -1
                    else:
                        #print 'Splice two boundary error: ',current_location,ref_seq
                        event[10]=0
            elif current_junction[0][0] == junction[0]:
                if event[0] == 4 or event[0] == 5:
                    event[0] = 5 #alternative splicing
                else:
                    event[0] = 3 # left side overlap
            elif current_junction[0][1] == junction[1]:
                if event[0] ==3 or event[0] == 5:
                    event[0] = 5 #alternative splicing
                else:
                    event[0] = 4 # right side overlap
            elif event[0] == 0:
                event[0] = 2
                continue
        if len(ref_seq_junction)<1 and event[1] > 0:
            event[0] = 2        
    else:
        event[0] = 2
    return event

'''
def FindSpliceEvent(current_location, ref_seq, event):
    
    if ref_seq[3][6] == '+':
        ref_strand = 1
    else:
        ref_strand = 0
    if ref_strand != int(current_location[5]):
        event[-1] = Bigger(event[-1], 1)
        # return event########## reverse
    case = CompareLocation(current_location[0][0][0], current_location[0][1][-1], ref_seq[0], ref_seq[1])
    if case == 10 or case == 4:
        return event
    elif case == 1 or case == 5:  #### fusion specification. 1: different strand 2: non overlap with cds, 3: overlap with cds but out of frame, 4: overlap with cds within frame match 5: overlap junction
        if sum(event[1:]) > 1:
            'if event contain some other event, ignore fusion case'
            return event
        elif case == 1:
            gframe = current_location[0][0][0] % 3
            exon_case = CheckRelation(current_location[0][0][0], current_location[0][1][0], gframe, ref_seq)
            if exon_case > 0:
                if exon_case in [3, 6, 9]:
                    event[0] = 5
                elif exon_case == 10:
                    event[0] = Bigger(event[0], 2)
                else:
                    event[0] = Bigger(event[0], 4)
            else:
                event[0] = Bigger(event[0], 3)
            return event
            'first half save it possible fusion or gene boundary'
        else:
            gframe = (current_location[0][0][-1] - (sum(current_location[0][1][:-1]) - sum(current_location[0][0][:-1]) % 3)) % 3 
            exon_case = CheckRelation(current_location[0][0][-1], current_location[0][1][-1], gframe, ref_seq)
            if exon_case == 10 or exon_case == 0:
                return event
            elif event[0] > 0 :
                if (exon_case in [7, 8, 9, -3, -2, -1] and event[0] > 2) or (event[0] > 4):
                    event[0] = 10  # determine fusion
                elif event[0] > 2:
                    event[0] = Bigger(event[0], 9)  # weak fusion
                else:
                    event[0] = Bigger(event[0], 8)  # very weak fusion
            else:
                if exon_case in [7, 8, 9, -3, -2, -1]:
                    event[1] = 1
                else:
                    event[2] = 1
            return event
            'second half if event contain prev case 1 fusion case then save it for a possible fusion if not save it for a gene boundary'
    else:
        exon = []
        for estart in range(len(current_location[0][0])):
            gframe = (current_location[0][0][estart] - (sum(current_location[0][1][:estart]) - sum(current_location[0][0][:estart]) % 3)) % 3
            exon.append(CheckRelation(current_location[0][0][estart], current_location[0][1][estart], gframe, ref_seq))
        if any([item in exon[:-1] for item in [3, 6, 9, -7, -4, -1]]) or any([item in exon[1:] for item in [7, 8, 9, -3, -2, -1]]):
            event[1] = 2
        else:
            event[2] = 2
        'for every part in current location, check the overlap status first'
        'if one junction is matched with splice side, record alternative splice (possibly be a exon boundary or novel exon needed to be check later)'
        'else if none junction is matched and have both side overlap record strong novel splice'
        'else if one side or half side overlap record week novel splice'
        'else if none side overlap record very week novel splice'
        
    
    return event
'''
def FindNonSpliceEvent(current_location, ref_seq, event):
    start = current_location[0][0][0]
    end = current_location[0][1][-1]
    gframe = start % 3
    if ref_seq[3][6] == '+':
        ref_strand = 1
    else:
        ref_strand = 0
    if ref_strand != int(current_location[5]):
        event[-1] = Bigger(event[-1], 1)
        return event  ########## reverse
    case = CompareLocation(current_location[0][0][0], current_location[0][1][-1], ref_seq[0], ref_seq[1])
    exon_case = CheckRelation(start, end, gframe, ref_seq)

    if case == 1 or case == 5:
        event[6] = 1
    if len(ref_seq[6]) > 0 :
        for utr in ref_seq[6]:
            #if utr[3] < start < utr[4] or utr[3] < end < utr[4]:
            if utr[3] <= start and end <= utr[4]:
                event[7] = 1
    if len(ref_seq[5]) > 0: 
        for utr in ref_seq[5]:
            if utr[3] <= start and end <= utr[4]:
            #if utr[3] < start < utr[4] or utr[3] < end < utr[4]:
                event[7] = 1
    if exon_case in [0, 10]:
        event[9] = 1
    elif exon_case in [1, 4, 5, 6, 7]:
        utr = ref_seq[:]
        utr[4] = ref_seq[6]
        if CheckRelation(start,end,gframe,utr) not in [0,10]:
            event[7] = 1
        utr[4] = ref_seq[5]
        if CheckRelation(start,end,gframe,utr) not in [0,10]:
            event[7] = 1
        event[8] = 1
    elif exon_case < 0:
        event[10] = 2
    else:
        #print 'exon_case error', exon_case, current_location, ref_seq, event
        event[10]=0
    #print event, case,exon_case, ref_seq,current_location
    return event

def PepString(peptide, location):
    beg = location[0][0]
    end = location[0][1]
    index = 0
    pivot = 2
    if location[-1] == 1:
        new_string = peptide[0:2]
        for i in range(len(beg) - 1):
            index += (end[i] - beg[i])
            if index < 3 or end[i] < 0:
                continue
            if peptide[0] != '-' and i == 0:
                index = index - 3
            new_string = new_string + peptide[pivot:index / 3 + 2] + ':'
            pivot = index / 3 + 2
        new_string += peptide[pivot:]
    else:
        new_string = peptide[-2:]
        for i in range(len(beg) - 1):
            index += (end[i] - beg[i])
            if index < 3 or end[i] < 0:
                continue
            if peptide[-1] != '-' and i == 0:
                index = index - 3
            new_string = ':' + peptide[-(index / 3 + 2):-pivot] + new_string
            pivot = index / 3 + 2
        new_string = peptide[:-pivot] + new_string
            
    return new_string

novelgene = 0
fusiongene = 0
alternativesplice = 0
novelsplice = 0
insertion = 0
mutation = 0
deletion = 0
translatedutr = 0
exonboundary = 0
novelexon = 0
geneboundary = 0
frameshift = 0
reversestrand = 0
pseudo = 0
IG = 0
indelcount = 0
na = 0


def CheckPseudo(event_group,chr):
    ret_val = [False,'','']
    if pseudo_gff_dic.get(chr) == None:
        return ret_val
    pseudo_ref = pseudo_gff_dic.get(chr)
    for pseudo in pseudo_ref:
        if pseudo[1][-1] < event_group[0][0][0][0][0]:
            continue
        elif pseudo[0][0] > event_group[0][0][0][1][-1]:
            break
        for index in range(len(pseudo[0])):
            for splice_index in range(len(event_group[0][0][0][0])):
                #print index,pseudo,pseudo[0][index],pseudo[1][index],event_group[0][0][0][0][splice_index],event_group[0][0][0][1][splice_index]
                pseudo_event = CompareLocation(event_group[0][0][0][0][splice_index],event_group[0][0][0][1][splice_index],pseudo[0][index],pseudo[1][index])
                if pseudo_event != 10:
                    ret_val[0] = True
                    ret_val[1] = pseudo[2]
                    ret_val[2] = pseudo[3]
    
    return ret_val


def ChooseEvent(event,event_group,chr,IG_check):
    global novelgene
    global fusiongene
    global alternativesplice
    global novelsplice
    global insertion
    global mutation
    global deletion
    global translatedutr
    global exonboundary
    global novelexon
    global geneboundary
    global frameshift
    global reversestrand
    global pseudo
    global IG
    global indelcount
    global na
 
    if IG_check:
        sevent = 'IG gene'
        IG += 1
    elif sum(event[:]) == 0 or (sum(event[:]) - event[1]) == 0:  # and event[-1] != 2:
        sevent = 'novel gene'
        pseudogene = CheckPseudo(event_group,chr)
        if pseudogene[0]:
            sevent = pseudogene[2]
            sevent = 'pseudogene'
            pseudo += 1
        else:
            novelgene += 1
    elif event[1] == 3:#event[0] in [9, 10]:# and event[1] == 0 and event[2] == 0:
        sevent = 'fusion gene'
        fusiongene += 1
    elif event[7] != 0 and event[8] != 0:
        sevent = 'translated UTR'
        translatedutr += 1
    elif event[8] != 0:
        sevent = 'exon boundary'
        exonboundary += 1 
    elif event[10] != 0:
        sevent = 'frame shift'
        frameshift += 1
    elif event[0] == 5:
        sevent = 'alternative splice'
        alternativesplice += 1
 
    elif event[7] != 0:
        sevent = 'translated UTR'
        translatedutr += 1    
    elif (event[0] > 0)and event[5] == 0:
        sevent = 'novel splice'
        novelsplice += 1
    elif event[3] != 0:
        sevent = 'insertion'
        indelcount += 1
        insertion += 1
    elif event[4] != 0:
        sevent = 'mutation'
        indelcount += 1
        mutation += 1
    elif event[5] != 0:
        sevent = 'deletion'
        indelcount += 1
        deletion += 1

    elif event[9] != 0:
        sevent = 'novel exon'
        pseudogene = CheckPseudo(event_group,chr)
        if pseudogene[0]:
            sevent = pseudogene[2]
            sevent = 'pseudogene'
            pseudo += 1
        else:
            novelexon += 1
    elif event[6] != 0:
        sevent = 'gene boundary'
        geneboundary += 1
    #elif event[10] != 0:
    #    sevent = 'frame shift'
    #    frameshift += 1
    elif event[11] != 0:
        sevent = 'reverse strand'
        pseudogene = CheckPseudo(event_group,chr)
        if pseudogene[0]:
            sevent = pseudogene[2]
            sevent = 'pseudogene'
            pseudo += 1
        else:
            reversestrand += 1
    else:
        sevent = 'NA'
        na += 1
    return sevent


def GbrowserLink(tmpGbrowserCoor):
    tmpGbrowserLink = "http://genome.ucsc.edu/cgi-bin/hgTracks?db="
    tmpGbrowserLink += db_str
    tmpGbrowserLink += "&position=" 
    tmpGbrowserLink += tmpGbrowserCoor
    tmpGbrowserLink += "&hgt.customText=http://bix.ucsd.edu/tmp/"
    tmpGbrowserLink += usa_gffile_str
    return tmpGbrowserLink



# read unknown location file
f = open(unknown_location_filename, 'r')
unknown_splice = {}
unknown_nonsplice = {}
unknown_indel = {}

indelCount = 0
for line in f:
    if line.startswith('#'):
        continue
    if line.strip() == '':
        continue
    line = line.strip()
    line = line.split('\t')
    chromosome = line[0]
    location_temp = line[1].split(';')
    pep = line[2]
    spec_count = line[3]
    location_count = line[4]
    fdr = line[5]
    strand = int(line[7].strip())
    original_pep = pep
    remaining = ''
    if len(line) > 8:
        if line[8][1] == '.':
            original_pep = line[8].strip()
            if len(line) > 9:
                remaining = line[9].strip()
                if len(line) > 10:
                    remaining = remaining + '\t' + line[10].strip()
                    if len(line) > 11:
                        remaining = remaining + '\t' + line[11].strip()
                        if len(line) > 12:
                            remaining = remaining + '\t' + line[12].strip()
        else:
            remaining = line[8].strip()
            
        
    
    location = [[], []]
    check = True
    for i in range(len(location_temp) - 1):
        if strand == 1:
            temp = location_temp[i].split('/')
        else:
            temp = location_temp[len(location_temp) - 2 - i].split('/')
        location[0].append(int(temp[0]))
        location[1].append(int(temp[1]))
        if int(temp[0]) < 0:
            check = False
    pep = PepString(original_pep, [location, strand])
          
    if len(location[0]) > 1 and check:
        if unknown_splice.has_key(chromosome):
            unknown_splice[chromosome].append([location, pep, spec_count, location_count, fdr, strand, remaining])
        else:
            unknown_splice[chromosome] = [[location, pep, spec_count, location_count, fdr, strand, remaining]]
    elif check == False:
        #indelCount += 1
        if unknown_indel.has_key(chromosome):
            unknown_indel[chromosome].append([location, pep, spec_count, location_count, fdr, strand, remaining])
        else:
            unknown_indel[chromosome] = [[location, pep, spec_count, location_count, fdr, strand, remaining]]
    else:
        if unknown_nonsplice.has_key(chromosome):
            unknown_nonsplice[chromosome].append([location, pep, spec_count, location_count, fdr, strand, remaining])
        else:
            unknown_nonsplice[chromosome] = [[location, pep, spec_count, location_count, fdr, strand, remaining]]

for i in unknown_nonsplice:
    unknown_nonsplice[i].sort(sort_location)
for i in unknown_indel:
    unknown_indel[i].sort(sort_location)
for i in unknown_splice:
    unknown_splice[i].sort(sort_location)

f.close()
log.write('Unknown file read complete\n')

# read known location file
f = open(known_location_filename, 'r')
known = {}

for line in f:
    if line.startswith('#'):
        continue
    line = line.split('\t')
    chromosome = line[0]
    location_temp = line[1].split(';')
    pep = line[2]
    spec_count = line[3]
    location_count = line[4]
    fdr = line[5]
    strand = int(line[7].strip())
    
    # original_pep   = line[8].strip()
    
    location = [[], []]
    for i in range(len(location_temp) - 1):
        temp = location_temp[i].split('/')
        location[0].append(int(temp[0]))
        location[1].append(int(temp[1]))  
    
    if known.has_key(chromosome):
        known[chromosome].append([location, pep, spec_count, location_count, fdr, strand])
    else:
        known[chromosome] = [[location, pep, spec_count, location_count, fdr, strand]]
        
for i in known:
    known[i].sort(sort_location)
        
f.close()
log.write('Known file read complete\n')

# read ref_seq file 
f = open(ref_seq_filename, 'r')
ref_seq_gff_dic = {}
pseudo_gff_dic = {}
gene_gff_dic = {}

prev_CDS_parent_ID = ""
prev_CDS_start = -1
prev_CDS_end = -1
prev_pseudo_ID = ''

for line in f:
    curr_part = line.strip()
    curr_part = curr_part.split('\t')
    if len(curr_part) != 9:
        continue
    chr_str = curr_part[0]
    curr_part[3] = int(curr_part[3]) - 1  # 0 base inclusive transfer
    curr_part[4] = int(curr_part[4])

    if curr_part[2] == "mRNA":
        data = curr_part[-1].split(";")
        for i in range(0, len(data)):
            if data[i].startswith("ID="):
                gene_name_str = data[i].replace("ID=", "")
            if data[i].startswith("Name="):
                gene_readable_name_str = data[i].replace("Name=", "")        
                break
        if ref_seq_gff_dic.has_key(chr_str):
            ref_seq_gff_dic[chr_str].append([curr_part[3], curr_part[4], gene_name_str , curr_part, [], [], [], [], gene_readable_name_str])
            curr_ind = len(ref_seq_gff_dic[chr_str]) - 1
        else:
            ref_seq_gff_dic[chr_str] = []
            ref_seq_gff_dic[chr_str].append([curr_part[3], curr_part[4], gene_name_str , curr_part, [], [], [], [], gene_readable_name_str])
            curr_ind = 0
            # [start, end, gene_name_str, mRNA, CDS, five_prime_UTR, three_prime_UTR, splice junctions, gene_readable_name_str]
    else:
        data = curr_part[-1].split(";")
        parent_name_str = ""
        for i in range(0, len(data)):
            if data[i].startswith("Parent="):
                parent_name_str = data[i].replace("Parent=", "")
                break

    if curr_part[2] == "CDS":
        if parent_name_str != ref_seq_gff_dic[chr_str][curr_ind][2]:
            print "No Parent", ref_seq_gff_dic[chr_str][curr_ind][2]
            break
        ref_seq_gff_dic[chr_str][curr_ind][4].append(curr_part)

        if prev_CDS_parent_ID == parent_name_str:
            splice_info = []
            if prev_CDS_start < int(curr_part[3]):
                j_start = prev_CDS_end
                j_end = int(curr_part[3])
            else:
                j_start = int(curr_part[4])
                j_end = prev_CDS_start
            splice_info.append(j_start)
            splice_info.append(j_end)
            ref_seq_gff_dic[chr_str][curr_ind][7].append(splice_info)

        prev_CDS_parent_ID = parent_name_str
        prev_CDS_start = int(curr_part[3])
        prev_CDS_end = int(curr_part[4])

    elif curr_part[2] == "five_prime_UTR":
        if parent_name_str != ref_seq_gff_dic[chr_str][curr_ind][2]:
            print "No Parent", ref_seq_gff_dic[chr_str][curr_ind][2]
            break
        ref_seq_gff_dic[chr_str][curr_ind][5].append(curr_part)

    elif curr_part[2] == "three_prime_UTR":
        if parent_name_str != ref_seq_gff_dic[chr_str][curr_ind][2]:
            print "No Parent", ref_seq_gff_dic[chr_str][curr_ind][2]
            break
        ref_seq_gff_dic[chr_str][curr_ind][6].append(curr_part)
    
    elif curr_part[1].find('pseudogene')>-1 and curr_part[2] == 'exon':
        if curr_part[8].find('Parent=') < 0:
            continue
        if not pseudo_gff_dic.has_key(curr_part[0]):
            pseudo_gff_dic[curr_part[0]] = [] #{chr:[[start1,start2],[end1,end2],'parent_name','pseudo name']
        curr_parent_ID = curr_part[8].split('Parent=')[1]
        if curr_parent_ID != prev_pseudo_ID:
            prev_pseudo_ID = curr_parent_ID
            pseudo_gff_dic[curr_part[0]].append([[int(curr_part[3])],[int(curr_part[4])],curr_parent_ID,curr_part[1]])
        else:
            pseudo_gff_dic[curr_part[0]][-1][0].append(int(curr_part[3]))
            pseudo_gff_dic[curr_part[0]][-1][1].append(int(curr_part[4]))
    
    elif curr_part[2] == "gene":
        if not gene_gff_dic.has_key(curr_part[0]):
            gene_gff_dic[curr_part[0]] = [] #{chr:[[start],[end],'parent_name}
        curr_part_ID = curr_part[8].split('Name=')[1]
        gene_gff_dic[curr_part[0]].append([[int(curr_part[3])],[int(curr_part[4])],curr_part_ID])
        
            
            

        

    else:
        continue

for i in ref_seq_gff_dic:
    ref_seq_gff_dic[i].sort(compare_start)

for i in pseudo_gff_dic:
    pseudo_gff_dic[i].sort(compare_pseudo)

for i in gene_gff_dic:
    gene_gff_dic[i].sort(compare_pseudo)

f.close()
log.write('Ref seq file read complete\n')
# Find event

event_set = []
event_group = []  # event_group[index - group of event indicated][0 - unknown, 1 - known][unknown_splice info / known_seq_info]
chromosome = []
related_gene_set = []
for chr in unknown_splice:
    for current_location in unknown_splice[chr]:
        strand = current_location[5]
        # current location : [location,pep,spec_count,location_count,fdr,strand], location : [[start1,start2],[end1,end2]]
        # record the start and end position of current peptide location
        beg_position = current_location[0][0][0]
        end_position = current_location[0][1][-1]
        event = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # [0:possible fusion, 1:alternative splice, 2:novel splice, 3:in, 4:mu, 5:del] 6:geneboundary, 7:translated UTR, 8:exon boundary, 9:novel exon, 10:frame shift 11: reverse_strand
        related_gene = []
        deletion = False
        for i in range(len(current_location[0][0]) - 1):
            if current_location[0][0][i + 1] - current_location[0][1][i] < 10:
                deletion = True
        if deletion == True:
            event[5] = 1
        if ref_seq_gff_dic.has_key(chr):  
            for ref_seq in ref_seq_gff_dic[chr]:
                if ref_seq[0] > end_position:  # pass any of the ref-seq if theres no overlap
                    break
                if ref_seq[1] < beg_position:
                    continue
                event = FindSpliceEvent(current_location, ref_seq, event)
                related_gene.append(ref_seq[8].split('-')[0])
        # print event, chr, current_location
        event_set.append(event)
        event_group.append(Recruit(current_location, unknown_nonsplice.get(chr), known.get(chr),unknown_splice.get(chr),unknown_indel.get(chr),1,strand))
        chromosome.append(chr)
        related_gene_set.append(related_gene)

for chr in unknown_indel:
    for current_location in unknown_indel[chr]:
        strand = current_location[5]
        if current_location[0][0][0] > 0:
            beg_position = current_location[0][0][0]
        else:
            beg_position = current_location[0][0][1]
        if current_location[0][1][-1] > 0:
            end_position = current_location[0][1][-1]
        else:
            end_position = current_location[0][1][-2]
        event = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        related_gene = []
        for i in current_location[0][1]:
            if i < 0:
                if i < -1000:
                    event[4] = 1
                else:
                    event[3] = 1
        if ref_seq_gff_dic.has_key(chr):                
            for ref_seq in ref_seq_gff_dic[chr]:
                if ref_seq[0] > end_position:
                    break
                if ref_seq[1] < beg_position:
                    continue
                related_gene.append(ref_seq[8].split('-')[0])
        event_set.append(event)
        event_group.append(Recruit(current_location, unknown_nonsplice.get(chr), known[chr],unknown_splice.get(chr),unknown_indel.get(chr),2,strand))
        chromosome.append(chr)
        related_gene_set.append(related_gene)

''' temporary stop the code. Some event are not detected in the novel databases (so they don't include any splice or indel). need to find out why this happen.

for chr in unknown_nonsplice:
    index = 0
    while index < len(unknown_nonsplice[chr]):   
    # for current_location in unknown_nonsplice[chr]:
        strand = current_location[5]
        current_location = unknown_nonsplice[chr][index]
        beg_position = current_location[0][0][0]
        end_position = current_location[0][1][-1]
        event = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        related_gene = []
        if ref_seq_gff_dic.has_key(chr):  
            for ref_seq in ref_seq_gff_dic[chr]:
                if ref_seq[0] > end_position:
                    break
                if ref_seq[1] < beg_position:
                    continue
                event = FindNonSpliceEvent(current_location, ref_seq, event)
                related_gene.append(ref_seq[8].split('-')[0])
        event_set.append(event)
        event_group.append(Recruit(current_location, unknown_nonsplice[chr], known[chr],unknown_splice.get(chr),unknown_indel.get(chr),3,strand))
        chromosome.append(chr)
        related_gene_set.append(related_gene)
        del unknown_nonsplice[chr][index]
'''
        # index += 1
        

# Write output file



novel_gene_index = []
novel_gene_list = []
for index in range(len(event_set)):
    if sum(event_set[index]) == 0:
        novel_gene_index.append([index,chromosome[index]])
        novel_gene_list.append(index)
        
novel_gene_index.sort(sort_by_chromosome)


dbCount = 0
Total_event = {}
for index in range(len(event_set)):
    string_out = ''
    if index in novel_gene_list:
        continue
    gene = set(related_gene_set[index])
    if gene == ['']:
        gene = ['NA']
    IG_check = False
    for i in gene:
        if i.find('IG')>-1:
            IG_check = True
    
    # peptide = MakePeptide(event_group[index][0][0][1],event_group[index][0][0][0])
    
    sprob = 1
    for i in range(len(event_group[index][0])):
        sprob = sprob * (1 - (1 - float(event_group[index][0][i][4])) / float(event_group[index][0][i][3]))
    sprob = 1 - sprob

    
    if filter:
        if sprob <= 0.5:
            continue
        false_tryptic = 0
        peptide_sequence = event_group[index][0][0][1]
        for i in range(len(peptide_sequence)-1):
            if ((peptide_sequence[i] == 'K' or peptide_sequence[i] == 'R') and peptide_sequence[i+1] != 'P'):
                false_tryptic += 1
        if false_tryptic > 3:
            continue

    current_event = ChooseEvent(event_set[index],event_group[index],chromosome[index],IG_check)
    string_out += (current_event + '\t' + event_group[index][0][0][1] + '\t' + chromosome[index] + '\t')
    string_out += WriteLocation(event_group[index][0][0])
    string_out += ('\t' + str(len(event_group[index][0])) + '\t' + str(len(event_group[index][1])) + '\t')
    spec_count = 0
    for i in range(len(event_group[index][0])):
        spec_count += int(event_group[index][0][i][2])
    string_out += (str(spec_count) + '\t' + str(event_group[index][0][0][3]) + '\t' + str(event_group[index][0][0][4]) + '\t')
    if event_group[index][0][0][5] == 1:
        strand = '+'
    else:
        strand = '-'
    string_out += (str(sprob) + '\t' + strand + '\t')
    scoord = event_group[index][0][0][0][0][0]
    ecoord = event_group[index][0][0][0][1][-1]
    if scoord < 0 :
        scoord = event_group[index][0][0][0][0][1]
    if ecoord < 0 :
        ecoord = event_group[index][0][0][0][1][-2]
    tmpGbrowserCoor = chromosome[index] + ':' + str(scoord - 100) + '-' + str(ecoord + 100)
    # print tmpGbrowserCoor
    string_out += (GbrowserLink(tmpGbrowserCoor) + '\t')
    
    strand_error_check = 0
    for i in range(len(event_group[index][0])):
        # strand_error_check += event_group[index][0][i][-1]
        strand_error_check += event_group[index][0][i][5]  # changed
    if strand_error_check != 0 and strand_error_check != len(event_group[index][0]):
        print 'Err: different strand come into same group.'
        print event_group[index][0], strand_error_check, len(event_group[index][0])
    
    for i in range(len(event_group[index][0]) - 1):
        string_out += (event_group[index][0][i + 1][1] + '/')
        string_out += WriteLocation(event_group[index][0][i + 1])
        string_out += ('|')
    string_out += ('\t')
    
    for i in gene:
        string_out += (i + ';')
    if len(event_group[index][0][0]) > 6:
        #if event_group[index][0][0][-1] != '':
        if event_group[index][0][0][-1].find('\t')> -1:# != '':
            dbCount += 1
        string_out += ('\t' + event_group[index][0][0][-1])
    string_out += ('\t\n')
    if Total_event.has_key(current_event):
        Total_event[current_event].append(string_out)
    else:
        Total_event[current_event] = [string_out]
        

novelgene_ws = 0
novelgene = 0
transcriptgene = 0
old_chr = ''
for item in novel_gene_index:
    string_out = ''
    index = item[0]
    chr = item[1]
    pseudogene = CheckPseudo(event_group[index],chr)
    related_gene_name = ''  
    if pseudogene[0]:
        #event = pseudogene[2]+' w/o'
        event = 'pseudogene'
        related_gene_set[index].append(pseudogene[1])
        pseudo += 1
    else:
        novel_gene_check = True
        for gene_gff in gene_gff_dic[chr]:
            if event_group[index][0][0][0][0][0] > gene_gff[1][0]:
                continue
            elif event_group[index][0][0][0][1][-1] < gene_gff[0][0]:
                break
            novel_gene_check = False
            event = 'transcript gene(non CDS)'
            related_gene_name = gene_gff[2]
            transcriptgene += 1
        
        if novel_gene_check:
            if DNA_folder != '':
                if old_chr != chr:
                    dna_file = open(DNA_folder+'/'+dna_string1+chr+dna_string2,'r')
                    dna = dna_file.readline()
                    old_chr = chr
                    #print DNA_folder+'/'+dna_string1+chr+dna_string2
                if not CheckStopCodon(dna[event_group[index][0][0][0][0][0]-90:event_group[index][0][0][0][0][0]],event_group[index][0][0][5]):
                    event = 'novel gene w/s'
                    novelgene_ws += 1
                elif not CheckStopCodon(dna[event_group[index][0][0][0][1][-1]:event_group[index][0][0][0][1][-1]+90],event_group[index][0][0][5]):
                    event = 'novel gene w/s'
                    novelgene_ws += 1
                else:
                    event = 'novel gene'
                    novelgene += 1
            else:
                event = 'novel gene'
                novelgene += 1
            
    
    sprob = 1
    for i in range(len(event_group[index][0])):
        sprob = sprob * (1 - (1 - float(event_group[index][0][i][4])) / float(event_group[index][0][i][3]))
    sprob = 1 - sprob
    if filter:
        if sprob <= 0.5:
            continue
        false_tryptic = 0
        peptide_sequence = event_group[index][0][0][1]
        for i in range(len(peptide_sequence)-1):
            if ((peptide_sequence[i] == 'K' or peptide_sequence[i] == 'R') and peptide_sequence[i+1] != 'P'):
                false_tryptic += 1
        if false_tryptic > 3:
            continue
    current_event = event
    if current_event == 'novel gene w/s':
        current_event = 'novel gene'
    string_out += (event + '\t' + event_group[index][0][0][1] + '\t' + chromosome[index] + '\t')
    string_out += WriteLocation(event_group[index][0][0])
    string_out += ('\t' + str(len(event_group[index][0])) + '\t' + str(len(event_group[index][1])) + '\t')
    spec_count = 0
    for i in range(len(event_group[index][0])):
        spec_count += int(event_group[index][0][i][2])
    string_out += (str(spec_count) + '\t' + str(event_group[index][0][0][3]) + '\t' + str(event_group[index][0][0][4]) + '\t')
    if event_group[index][0][0][5] == 1:
        strand = '+'
    else:
        strand = '-'
    string_out += (str(sprob) + '\t' + strand + '\t')
    scoord = event_group[index][0][0][0][0][0]
    ecoord = event_group[index][0][0][0][1][-1]
    if scoord < 0 :
        scoord = event_group[index][0][0][0][0][1]
    if ecoord < 0 :
        ecoord = event_group[index][0][0][0][1][-2]
    tmpGbrowserCoor = chromosome[index] + ':' + str(scoord - 100) + '-' + str(ecoord + 100)
    # print tmpGbrowserCoor
    string_out += (GbrowserLink(tmpGbrowserCoor) + '\t')
    
    strand_error_check = 0
    for i in range(len(event_group[index][0])):
        # strand_error_check += event_group[index][0][i][-1]
        strand_error_check += event_group[index][0][i][5]  # changed
    if strand_error_check != 0 and strand_error_check != len(event_group[index][0]):
        print 'Err: different strand come into same group.'
        print event_group[index][0], strand_error_check, len(event_group[index][0])
    
    for i in range(len(event_group[index][0]) - 1):
        string_out += (event_group[index][0][i + 1][1] + '/')
        string_out += WriteLocation(event_group[index][0][i + 1])
        string_out += ('|')
    string_out += ('\t')
    gene = set(related_gene_set[index])
    string_out += (related_gene_name)
    if gene == set([]):
        gene = (['NA'])
    for i in gene:
        string_out += (i + ';')
    if len(event_group[index][0][0]) > 6:
        if event_group[index][0][0][-1] != '':
            dbCount += 1
        string_out += ('\t' + event_group[index][0][0][-1])
    string_out += ('\t\n')
    if Total_event.has_key(current_event):
        Total_event[current_event].append(string_out)
    else:
        Total_event[current_event] = [string_out]

def MergeEvent(list_event):
    string = ''
    index = 0
    for event in list_event:
        if index == 0:
            noCDS = event.split('\t')
            start = int(noCDS[3].split('-')[0]) -300
            end = int(noCDS[3].split('-')[-1]) +300
            Sprob = [float(noCDS[9])]
        else:
            current = event.split('\t')
            noCDS[12] += current[1]+'/'+current[3]+'|'+current[12]

            temp_start = int(current[3].split('-')[0])
            temp_end = int(current[3].split('-')[-1])
            if start > temp_start:
                start = temp_start-300
            if end < temp_end:
                end = temp_end+300
            noCDS[4] = int(noCDS[4]) + int(current[4])
            noCDS[5] = int(noCDS[5]) + int(current[5])
            noCDS[6] = int(noCDS[6]) + int(current[6])
            Sprob.append(float(noCDS[9]))
        index += 1
    prob = 1
    for s in Sprob:
        prob = prob * (1-s)
    prob = 1 - prob
    noCDS[9] = str(prob)
    temp = noCDS[11].split('&')
    noCDS[11] = temp[0]+'&'+temp[1].split(':')[0]+':'+str(start)+'-'+str(end)+'&'+temp[2]
    for i in range(len(noCDS)):
        noCDS[i] = str(noCDS[i])
    string = '\t'.join(noCDS)
    return string


for event in Total_event:
    if event.startswith('transcript') or event.startswith('pseudo') or event.startswith('fusion') or event.startswith('translated'):
        sub_event = Total_event[event]
        noCDS_event = {}
        for each_event in sub_event:
            data = each_event.split('\t')
            gene = data[13]
            if noCDS_event.has_key(data[13]):
                noCDS_event[data[13]].append(each_event)
            else:
                noCDS_event[data[13]] = [each_event] 
        
        substitution_event = []
        for gene in noCDS_event:
            new_event_string = MergeEvent(noCDS_event[gene])
            substitution_event.append(new_event_string)
            
        Total_event[event] = substitution_event
    if event.startswith('novel gene'):
        sub_event = Total_event[event]
        novelgene_event = {}
        for each_event in sub_event:
            data = each_event.split('\t')
            chr = data[2]
            coord = data[3].split('-')
            coord_key = (int(coord[0]) + int(coord[-1]))/2
            is_written = False
            for keys in novelgene_event:
                if (keys[1] + coord_key) /2 < 6000 and keys[0] == chr:
                    novelgene_event[keys].append(each_event)
                    is_written = True
            if not is_written:
                novelgene_event[(chr,coord_key)] = [each_event]
        substitution_event = []
        for keys in novelgene_event:
            new_event_string = MergeEvent(novelgene_event[keys])
            substitution_event.append(new_event_string)
        Total_event[event] = substitution_event
                    
                    


ind = 0
header_length = len(header.split('\t'))
for event in Total_event:
    for each in Total_event[event]:
        if len(each.strip().split('\t')) != header_length-1:
            print 'ERR: header is not appropriately matched with the contents'
            continue
        else:
            each = each.strip()+'\n'
        ind += 1
        s.write(str(ind)+'\t'+each)
    

    
    
    #print event_group[index][0][0]
#[[[26679959, 26680266], [26680045, 26680288]], 'R.ALYAIPG:LDYVSHEDILPYTSTDQVPIQHELFER.F', '35', '1', '0.0', 0, '']


s.close()



'''
s = open(output_filename,'r')

bcount = 0
for line in s:
    if line.startswith('#') or line == "":
        continue
    
    line = line.split('\t')
    chr = line[3]
    if line[1] != 'novel gene':
        continue
    if line[4].find('--')>-1:
        continue
    location = line[4].split(';')
    start = []
    end = []
    for i in range(len(location)):
        temp = location[i].split('-')
        start.append(int(temp[0]))
        end.append(int(temp[1]))
    #print line, start, end
    for ref_seq in ref_seq_gff_dic[chr]:
        
        if ref_seq[0] > end[-1]: # pass any of the ref-seq if theres no overlap
            break
        if ref_seq[1] < start[0]:
            continue
        bcount += 1
        
        print start[0],end[-1],ref_seq[0],ref_seq[1], line
        case = CompareLocation(start[0],end[-1],ref_seq[0],ref_seq[1])
        if case != 10:
            print case, ref_seq[-1],line, line[2]
'''
        

'''
for j,i in enumerate(event_group):
    if len(i[0])>1:#+ len(i[1])>1:
        print event_set[j],i
'''
indel_count = 0
total_count = 0
for event in Total_event:
    total_count += len(Total_event.get(event))
    print event ,': ',len(Total_event.get(event))
    s2.write(event+'\t'+str(len(Total_event.get(event)))+'\n')
    if event in ['insertion','deletion','mutation']:
        indel_count += len(Total_event.get(event))
print 'total number of event group: ', total_count, ', indelCount: ', indel_count, ', dbCount: ', dbCount
s2.write('Sum\t'+str(total_count))
'''
print 'total number of event group: ', len(event_group), ', indelCount: ', indelCount, ', dbCount: ', dbCount

print    'novel gene: ',novelgene
print    'novel gene with stopcodon: ',novelgene_ws
print    'transcript gene(non CDS): ',transcriptgene
print    'fusiton gene: ',fusiongene
print    'alternative splice: ',alternativesplice
print    'novel splice: ',novelsplice
print    'insertion: ',insertion
print    'mutation: ',mutation
print    'deletion: ',deletion
print    'translatedUTR: ',translatedutr
print    'exon boundary: ',exonboundary
print    'novel exon: ',novelexon
print    'gene boundary: ', geneboundary
print    'frame shift: ',frameshift
print    'reverse strand: ',reversestrand
print    'pseudo gene: ',pseudo
print    'IG gene: ',IG
print    'na: ',na
'''

# print 'END: ', time.ctime()
# print 'TIME: ', time.clock() - btime, ' sec'
log.write('TIME: ' + str(time.clock() - btime) + ' sec')

