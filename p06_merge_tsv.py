import os,glob
import argparse
import subprocess
from natsort import natsorted

parser = argparse.ArgumentParser(description='Generate command file to transfer mzid to tsv')
parser.add_argument('-p',action='store',dest='tsv',help='path to tsv folder and cmd file',default='out_tsv')


args = parser.parse_args()

folder = args.tsv

# cat files
dbs = natsorted([f for f in os.listdir(folder) if os.path.isdir(folder+'/'+f)])
for db in dbs:
    files = glob.glob(folder+'/'+db+'/*.tsv')
    if files == []: raise ValueError('empty folder')
    out_fn = folder+'/'+db+'/'+db+'.tsv'
    with open(out_fn,'w') as out_f:
        n = 0
        for f in files:
            with open(f) as in_f:
                for line in in_f:
                    if not line.startswith('#'):
                        out_f.write(line)
                    else:
                        if n == 0:
                            out_f.write(line)
                        else:
                            pass
            n += 1
    for f in files:
        os.remove(f)