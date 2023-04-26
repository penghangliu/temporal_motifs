import numpy as np
import pandas as pd
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

def rpb(st):
    if st[0]==st[2]:
        if st[1]==st[3]:
            return "r"
        else:
            return "o"
    elif st[1]==st[2]:
        if st[0]==st[3]:
            return "p"
        else:
            return "c"
    elif st[0]==st[3]:
        return "w"
    elif st[1]==st[3]:
        return "i"
    else:
        return "_"

def label(st):
    first = st[0:4]
    second = st[2:6]
    third = st[4:8]
    return rpb(first) + rpb(second)


a = ["r","p","o","i","c","w"]
A = []
for i in range(6):
    for j in range(6):
        A.append(a[i]+a[j])

def count(data):
    # data = "calls.txt_1000_4000_4_4"
    M = np.loadtxt(data, dtype=str, usecols=(0,1))
    dfR = pd.DataFrame(M)
    dfR.columns = ['motif', 'count']
    df = dfR[3:]
    df = df.reset_index(drop=True)
    df['count'] = df['count'].astype(int)
    
    dfL = pd.DataFrame({'label':A})
    dfL['count'] = np.zeros((dfL.shape[0]), dtype=int)
    
    b = 0
    for i in range(df.shape[0]):
        l = label(df['motif'][i])
        b += df['count'][i]
        dfL['count'][dfL['label']==l] += df['count'][i]

    return dfL['count']

fdf = pd.DataFrame({'label':A})
fdf['count'] = count(infile)

fdf.to_csv(outfile, index=None)
