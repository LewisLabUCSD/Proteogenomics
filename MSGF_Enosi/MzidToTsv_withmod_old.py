'''
Created on Apr 16, 2013

@author: Akard3
'''

import sys
import os




if len(sys.argv)>1:
    Mzid_file = sys.argv[1]
    Tsv_file = sys.argv[2]
else:
    Mzid_file = '/home/s3cha/s3cha/data/MzidToTsv'
    Tsv_file = '/home/s3cha/s3cha/data/MzidToTsv/sample.tsv'
    

###############################################
#[*SpecFile]    [*SpecID   = ScanNum]    FragMethod    [*Precursor]    IsotopeError    
#PrecursorError(ppm)    [*Charge]    [*Peptide]    [*Protein]    *DeNovoScore    *MSGFScore    
#**SpecEValue    *EValue    QValue    PepQValue



def readMzid(Mzid_file):
    input_file = open(Mzid_file,'r')
    Protein = []
    Peptide = []
    FragMethod = []
    SpecID = []
    Precursor = []
    IsotopeError = []
    Charge = []
    DeNovoScore = []
    MSGFScore = []
    SpecEValue = []
    EValue = []
    QValue = []
    PepQValue = []
    SpecFile = []
    
    hash_key = {}
    
    PepID = {}
    ProtID = {}
    PepEvID = {}
    
    index_num = 0
    for line in input_file:
        #if count > 10:
        #    break
        if line.find('<DBSequence')>-1:
            line = line.split('"')
            ProtID[line[7]] = line[1]
            continue
        if line.find('<Peptide id=')>-1:
            pepid = line.split('"')[1]
        if line.find('<PeptideSequence>')>-1:
            line = line.split('>')[1].split('<')[0]
            peplength = len(line)
            PepID[pepid] = line
            Pep = line
            continue
        
        
        if line.find('<Modification monoisotopicMassDelta')>-1:
            line = line.split('"')
            mass = line[1].split('.')[0]+'.'+line[1].split('.')[1][:3]
            location = peplength-int(line[3])
            #print location,Pep,mass
            if not mass.startswith('-'):
                mass = '+'+mass
            if location == 0:
                Pep = Pep+mass
            else:
                Pep = Pep[:-location] + mass + Pep[-location:]  
            
            #Peptide[-1] = Peptide[-1][:-location]+mass+Peptide[-1][-location:]
            
            PepID[pepid] = Pep
            continue


        if line.find('<PeptideEvidence isDe')>-1:
            line = line.split('"')
            PepEvID[line[15]] = [line[5]+'.'+PepID[line[11]]+'.'+line[3],ProtID[line[13]]]
            continue
        if line.find('<SpectraData location=')>-1:
            line = line.split('"')
            SpecFiletemp = line[3]
            continue
        if line.find('<SpectrumIdentificationResult')>-1:
            FragMethod.append('')
            SpecID.append('')
            Peptide.append('')
            Protein.append('')
            SpecFile.append(SpecFiletemp)
            Precursor.append(0)
            IsotopeError.append(0)
            Charge.append(0)
            DeNovoScore.append(0)
            MSGFScore.append(0)
            SpecEValue.append(10000)
            EValue.append(0)
            QValue.append(0)
            PepQValue.append(0)
            line = line.split('"')
            SpecID[-1] = line[3]
            hash_key[SpecFile[-1]+SpecID[-1]] = index_num
            index_num += 1
            continue
        line = line.split('"')
        if ' experimentalMassToCharge=' in line:
            Precursor[-1] = round(float(line[9]),4)
            Charge[-1] = int(line[11])
            #pep = line[5]
            continue
        if line[0].find('<PeptideEvidenceRef peptideEvidence_ref=')>-1:
            Peptide[-1] = PepEvID[line[1]][0]
            Protein[-1] = PepEvID[line[1]][1]
        if 'MS-GF:RawScore' in line:
            MSGFScore[-1] = int(line[5])
            continue
        if 'MS-GF:DeNovoScore' in line:
            DeNovoScore[-1] = int(line[5])
            continue
        if 'MS-GF:SpecEValue' in line:
            SpecEValue[-1] = float(line[5])
            continue
        if 'MS-GF:EValue' in line:
            EValue[-1] = float(line[5])
            continue
        if 'MS-GF:QValue' in line:
            QValue[-1] = float(line[5])
            continue
        if 'MS-GF:PepQValue' in line:
            PepQValue[-1] = float(line[5])
            continue
        if 'IsotopeError' in line:
            IsotopeError[-1] = int(line[1])
            continue
        if 'AssumedDissociationMethod' in line:
            FragMethod[-1] = line[1]
            continue
    return [SpecFile,SpecID,FragMethod,Precursor,IsotopeError,Charge,Peptide,Protein,
            DeNovoScore,MSGFScore,SpecEValue,EValue,QValue,PepQValue,hash_key]
        
def writeTsv(series,Tsv_file):
    s = open(Tsv_file,'w')
    s.write('#SpecFile\tSpecID\tScanNum\tFragMethod\tPrecursor\tIsotopeError\tPrecursorError(ppm)\tCharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\tSpecEValue\tEValue\tQValue\tPepQValue\n')
    #\tSpecID\tScanNum\tFragMethod\tPrecursor\tIsotopeError\t
    #PrecursorError(ppm)\tCharge\tPeptide\tProtein\tDeNovoScore\tMSGFScore\t
    #SpecEValue\tEValue\tQValue\tPepQValue
    for index in range(len(series[1])):
        if series[10][index] == 10000:
            continue
        s.write(str(series[0][index])+'\t') #SpecFile
        s.write(str(series[1][index])+'\t') #SpecID
        for i in range(len(series[1][index].split('='))):
            if series[1][index].split('=')[i].find("scan")>-1 or series[1][index].split('=')[i].find("index")>-1:
                s.write(str(series[1][index].split('=')[i+1])+'\t')
        #s.write(str(series[1][index].split('=')[1])+'\t') #ScanNum
        s.write(str(series[2][index])+'\t') #FragMethod
        s.write(str(series[3][index])+'\t') #Precursor
        s.write(str(series[4][index])+'\t') #Isotope
        s.write('NA'+'\t') #PrecursorError
        s.write(str(series[5][index])+'\t') #charge
        s.write(str(series[6][index])+'\t') #Peptide
        s.write(str(series[7][index])+'\t') #Protein
        s.write(str(series[8][index])+'\t') #DeNovoScore
        s.write(str(series[9][index])+'\t') #MSGFScore
        s.write(str(series[10][index])+'\t') #SpecEvalue
        s.write(str(series[11][index])+'\t') #EValue
        s.write(str(series[12][index])+'\t') #QValue
        s.write(str(series[13][index])+'\n') #PepQvalue
       
def getFileList(dirname):
    current_path = dirname
    return_list = []
    file_list = os.listdir(dirname)
    for filename in file_list:
        filename = current_path+'/'+filename
        if os.path.isdir(filename):
            return_list += getFileList(filename)
        else:
            return_list.append(filename)
    return return_list

def addSeries(stack_series,series):
    for key in series[-1].keys():
        if stack_series[-1].has_key(key):
            stack_index = stack_series[-1][key]
            series_index = series[-1][key]
            #if len(stack_series[10]) <= stack_index:
		#print "Invalid line ", stack_series[10]
		#continue
            if stack_series[10][stack_index] < series[10][series_index]:
                continue
            else:
                stack_series[2][stack_index] = series[2][series_index]
                stack_series[3][stack_index] = series[3][series_index]
                stack_series[4][stack_index] = series[4][series_index]
                stack_series[5][stack_index] = series[5][series_index]
                stack_series[6][stack_index] = series[6][series_index]
                stack_series[7][stack_index] = series[7][series_index]
                stack_series[8][stack_index] = series[8][series_index]
                stack_series[9][stack_index] = series[9][series_index]
                stack_series[10][stack_index] = series[10][series_index]
                stack_series[11][stack_index] = series[11][series_index]
                stack_series[12][stack_index] = series[12][series_index]
                stack_series[13][stack_index] = series[13][series_index]
        else:
            series_index = series[-1][key]
            stack_series[0].append(series[0][series_index])
            stack_series[1].append(series[1][series_index])
            stack_series[2].append(series[2][series_index])
            stack_series[3].append(series[3][series_index])
            stack_series[4].append(series[4][series_index])
            stack_series[5].append(series[5][series_index])
            stack_series[6].append(series[6][series_index])
            stack_series[7].append(series[7][series_index])
            stack_series[8].append(series[8][series_index])
            stack_series[9].append(series[9][series_index])
            stack_series[10].append(series[10][series_index])
            stack_series[11].append(series[11][series_index])
            stack_series[12].append(series[12][series_index])
            stack_series[13].append(series[13][series_index])
            stack_series[14][key] = len(stack_series[1]) - 1
            continue
    return stack_series


    
stack_series = []
if os.path.isdir(Mzid_file):
    file_list = getFileList(Mzid_file)
    for filename in file_list:
        if os.path.splitext(filename)[1] != '.mzid':
            continue
        print 'Read ' + filename
        series = readMzid(filename)
        if stack_series == []:
            stack_series = series
        else:
            stack_series = addSeries(stack_series,series)
    writeTsv(stack_series,Tsv_file)
        
else:
    print 'Read ' + Mzid_file
    series = readMzid(Mzid_file)
    writeTsv(series,Tsv_file)

    



print 'END: '
  
