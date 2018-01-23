#!/usr/bin/python3
import numpy as np
from sys import argv
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
times = np.loadtxt(argv[1])

#speedup vs total number of working ranks

# grouping per p-space partitions
plt.figure()
t_per_p= [ (times[times[:,0]==c ][:,[1,3]],c) for c in set(times[:,0]) ]
for arr,c in t_per_p:
    plt.plot(arr[:,0],arr[0,1]/arr[:,1]/arr[:,0],label = "Pp="+str(c))
      
plt.xlabel("V-partitions")
plt.ylim([0,2])
plt.ylabel("Deviation from exp.scaling")
plt.legend()

plt.figure()
for arr,c in t_per_p:
    nranks = arr[:,0] * c
    speedup = times[0,3]/arr[:,1]
    plt.plot(nranks,speedup,label = "Pp="+str(c), linestyle = "None", marker = '+')
      
plt.xlabel("no of working ranks")
plt.ylabel("speedup")
plt.legend()


# grouping per v-space partitions
plt.figure()
t_per_v= [ (times[times[:,1]==c ][:,[0,3]],c) for c in set(times[:,1])]
for arr,c in t_per_v:
    plt.plot(arr[:,0],arr[0,1]/arr[:,1]/arr[:,0],label = "Pv="+str(c))
    
plt.xlabel("P-partitions")
plt.ylabel("Deviation from exp.scaling")
plt.ylim([0,2])
plt.legend()

plt.figure()
for arr,c in t_per_v:
    nranks = arr[:,0] * c
    speedup = times[0,3]/arr[:,1]
    plt.plot(nranks,speedup,label = "Pv="+str(c), linestyle = "None", marker = '+')

plt.xlabel("no of working ranks")
plt.ylabel("speedup")
plt.legend()



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
z = times[:,0]*times[:,1]*times[:,3]
ax.set_xlabel('Pp')
ax.set_ylabel('Pv')
ax.set_zlabel('Deviation from exp.scaling')
ax.plot_trisurf(times[:,0], times[:,1],z[0]/z)

plt.show()
