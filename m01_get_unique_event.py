import sys
import pandas as pd
event_fn = sys.argv[1]
out_fn = sys.argv[2]

df = pd.read_csv(event_fn,sep='\t',header=0,usecols=[1,2,3,4,11])
new_df = df.drop_duplicates()
new_df.to_csv(out_fn,sep='\t',index=False)