import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import os

#print("run from 'Research' directory")

#files=glob.glob('MISHAPS_F1_N*_r/MN*RAdec.txt')
files=glob.glob('MISHAPS_F1_*_r/MN*_1RAdec.txt')

#print(os.getcwd())
#print(files)

#plot asteroids	
for i in files:
	#print(i)							#use this line to find out which file is causing the code to break
	d=pd.read_csv(i, sep = '\s+', header = None)
	meanX=np.mean(d.iloc[:,0])
	meanY=np.mean(d.iloc[:,1])
	asteroid=os.path.basename(i)
	plt.text(meanX, meanY, asteroid, fontsize=6, color='cadetblue')
	plt.scatter(d.iloc[:,0], d.iloc[:,1], s=1, c='darkblue', marker='o')

#plot vertices
vertices=pd.read_csv("MISHAPS_F1.vertices", header=None, sep='\s+', skip_blank_lines=False)
plt.plot(vertices.iloc[:,0], vertices.iloc[:,1], 'k-')
plt.gca().set_aspect(aspect=1.0/np.cos(28.25*np.pi/180.0))

#plot centers
centers=pd.read_csv("MISHAPS_F1.centers", header=None, sep='\s+', skip_blank_lines=False)
for i in centers.iterrows():
	plt.text(i[1].iloc[1], i[1].iloc[2], i[1].iloc[0][11:], ha='center', va='center')

plt.show()

