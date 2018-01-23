#!/usr/bin/python3
from sys import argv
from os import path
import numpy as np
import pandas as pd
import glob
from scipy.interpolate import griddata
from matplotlib import pyplot as plt

try: 
    dirname = argv[1]
except:
    dirname = '.'
    
try: 
    res_frame = int(argv[2])
except:
    res_frame = -1

# MESH INFO
f = open(path.join(dirname,'reentry.plt'))

nxy_data= f.readlines()
f.close()
start = [ i for i,t in enumerate(nxy_data) if 'coordinates' in t ][0]+1
end = [ i for i,t in enumerate(nxy_data) if 'boundary sides' in t ][0]
nxy_data = '\n'.join(nxy_data[start:end])
nxy_data = np.fromstring(nxy_data,sep=' ').reshape((-1,5))[:,0:3]
xny_df = pd.DataFrame(data=nxy_data[:,1:],index = nxy_data[:,0].astype(np.int),columns=['x','y'])
#PARTITIONING INFO
part_data = np.loadtxt(glob.glob(path.join(dirname,'reentry.con.npart.*'))[0])
pd_df = pd.DataFrame(data=part_data,index = np.arange(1,part_data.size+1), columns = ['P'])

part_dfs = [ pd_df.loc[ pd_df['P'] == i ] for i in range(int(pd_df.values.min()),int(pd_df.values.max()+1))]
parts = [ part_df.join(xny_df)[['x','y']].values for part_df in part_dfs ]

# RESULTS 
def parse_results(file_name):
    f = open(file_name)
    r_data = f.readlines()
    f.close()
    idxs = [i for i,t in enumerate(r_data) if 'RESULTS' in t] + [len(r_data)]
    starts = idxs[:-1]
    ends = idxs[1:]
    results = []
    for s,e in zip(starts,ends):
        data = '\n'.join(r_data[s+1:e])
        data = np.fromstring(data,sep=' ').reshape((-1,4))
        results.append(data)
    return results

r1 = parse_results(path.join(dirname,'RESULTS1.RES'))
r2 = parse_results(path.join(dirname,'RESULTS2.RES'))

nd = r1[-1][:,1]
vx = r1[-1][:,2]
vy = r1[-1][:,3]
rho  = r2[-1][:,1]
ps   = r2[-1][:,2]
temp = r2[-1][:,3]

r1_dfs = [ pd.DataFrame(data=r1_data[:,1:],
                        index = r1_data[:,0].astype(np.int),
                        columns=['ND','U','V']) for r1_data in r1 ] 

r2_dfs = [ pd.DataFrame(data=r2_data[:,1:],
                        index = r2_data[:,0].astype(np.int),
                        columns=['RHO','PS','TEMP']) for r2_data in r2 ] 
all_data = [xny_df.join(r1_df) for r1_df in r1_dfs ]
all_data = [d.join(r2_df) for d,r2_df in zip(all_data,r2_dfs) ]

x = all_data[0][['x']].values.flatten()
y = all_data[0][['y']].values.flatten()
#nds = [ data[['ND']].values for data in all_data ]
#Us = [ data[['U']].values for data in all_data ]
#Vs = [ data[['V']].values for data in all_data ]
#RHOs = [ data[['RHO']].values for data in all_data ]
#PSs = [ data[['PS']].values for data in all_data ]
TEMPs = [ data[['TEMP']].values for data in all_data ]

xi = np.linspace(x.min(),x.max(),800)
yi = np.linspace(y.min(),y.max(),800)

#for qname in ['ND','U','V','RHO','PS','TEMP']:
for qname in ['TEMP']:
    qs = [ data[[qname]].values for data in all_data ] 
    plt.figure()
    plt.title(qname)
    for part in parts:
        plt.plot(part[:,0],part[:,1],linestyle='None',marker='+')
    zi = griddata((x,y),qs[res_frame].flatten(), (xi[None,:], yi[:,None]), method='nearest')
    plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
    plt.colorbar()

plt.show()
