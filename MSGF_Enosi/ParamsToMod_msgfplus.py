'''
Created on Apr 21, 2015

@author: s3cha
'''

import os
import sys

param_filename = sys.argv[1]
mod_filename = sys.argv[2]

f = open(param_filename,'r')

MOD = {"ptm.DEAMIDATION":"H-1N-1O1,NQ,opt,any,Deamidation",
       "ptm.LYSINE_METHYLATION":"+14.015650,K,opt,any,Methylation",
       "ptm.NTERM_ACETYLATION":"C2H2O,*,opt,Prot-N-term,Acetylation",
       "ptm.NTERM_CARBAMYLATION":"C2H3NO,*,opt,N-term,Carbamidomethylation",
       "ptm.OXIDATION":"O1,M,opt,any,Oxidation",
       "ptm.PHOSPHORYLATION":"HO3P,STY,opt,any,Phosphorylation",
       "ptm.PYROGLUTAMATE_FORMATION":"H-3N-1,Q,opt,N-term,Pyro-glu",
       "c57":"C2H3N1O1,C,fix,any,Carbamidomethyl",
       "c58":"C2H2O2,C,fix,any,Carboxymethyl",
       "c99":"99,C,fix,any,NIPIA",
       "None":"###None"}

OPTION = {"fix_nterm":"fix,N-term,CUSTOM","fix":"fix,any,CUSTOM","opt_nterm":"opt,N-term,CUSTOM","opt":"opt,any,CUSTOM"}

string = ""
for line in f:
    if line.find("ptm.") == -1 and line.find("cysteine_protease.cysteine") == -1:
        continue
    index = line.find("ptm.")
    data = line.strip()[index:].split('"')[0]
    if MOD.has_key(data):
        string += MOD.get(data)+'\n'
    else:
        if line.find("ptm.custom_PTM")>-1:
            data = line.split('>')[1].split('<')[0].strip()
            string += ','.join(data.split(',')[:-1])+','+OPTION.get(data.split(',')[-1])+'\n'
        elif line.find("ptm.mods")>-1:
            data = line.split('>')[1].split('<')[0].strip()
            string = "NumMods="+data+'\n'+string
	elif line.find("cysteine_protease.cysteine")>-1:
	    data = line.split('>')[1].split('<')[0].strip()
	    if data in MOD:
		string += MOD.get(data)+'\n'
string = string.strip()
  

s = open(mod_filename,'w')
s.write(string)

