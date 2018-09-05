import pandas as pd
class event:
    def __init__(self,fn):
        self.fn = fn
    
    def count_uniq_peptide(self):
        peptide = []
        with open(self.fn) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                item = line.split('\t')
                peptide.append(item[2])
        return len(set(peptide))
    
    def event_stats(self):
        event_dic = {}
        peptide = []
        with open(self.fn) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                item = line.split('\t')
                event = item[1]
                pep = item[2]
                if pep not in peptide:
                    peptide.append(pep)
                    if event in event_dic:
                        event_dic[event] += 1
                    else:
                        event_dic[event] = 1
        return event_dic