import os,glob
import argparse
import multiprocessing as mp
import subprocess
from natsort import natsorted

def run_cmd(cmd):
    subprocess.call(cmd.strip(),shell=True)

path = os.getcwd()

parser = argparse.ArgumentParser(description='Generate command file to transfer mzid to tsv')
parser.add_argument('-p',action='store',dest='tsv',help='path to tsv folder and cmd file',default='out_tsv')
parser.add_argument('-t',action='store',dest='thread',type=int,help='number of threads')

args = parser.parse_args()

folder = args.tsv
thread = args.thread

cmd_fn = glob.glob(folder+'/*.txt')[0]
handle = open(cmd_fn)
cmds = natsorted(handle.readlines())

pool = mp.Pool(processes=thread)
for cmd in cmds:
    pool.apply_async(run_cmd,args=(cmd,))
pool.close()
pool.join()



