'''
this file converts spectrum to indexed spectrum file
'''
from natsort import natsorted
import os,subprocess
import argparse

path = os.getcwd()

parser = argparse.ArgumentParser(description='Converts spectrum to indexed spectrum file')
parser.add_argument('-p','--path',action='store',dest='path',help='path to the main foler of the enosi',default=path)
parser.add_argument('-s','--snum',action='store',dest='start_index',type=int,help='start index of the spectrum file')
parser.add_argument('-n','--npath',action='store',dest='new_path',help='new result path',default='convert_spectrum')
args = parser.parse_args()


n = args.start_index
path = args.path
new_spec_path = args.new_path


os.chdir(path)

mgfs = natsorted(subprocess.check_output('find spectrum -name *.mgf',shell=True).split('\n'))
 
with open('mgf_fn_match.txt','a') as f:
    for mgf in mgfs:
        if mgf == '':
            continue
        new_mgf = 'spectrum-'+str(n).zfill(5)+'.mgf'
        if not os.path.exists(new_spec_path): os.mkdir(new_spec_path)
        os.rename(mgf,new_spec_path+'/'+new_mgf)
        f.write(new_mgf+'\t'+mgf+'\n')
        n += 1
