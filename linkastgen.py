import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import sys
import os
'''for linking an unspecified number of asteroids'''

#manager = plt.get_current_fig_manager()
#manager.resize(*manager.window.maxsize())

#if len(sys.argv)!=3:
	#print("usage: %s <asteroid1name> <asteroid2name>" % (sys.argv[0]))
	#sys.exit()

#number of asteroids
n=len(sys.argv)-1

#appends arbitrary number of arguments into a list
asteroids=[]	
for i in sys.argv:
	if i == sys.argv[0]:
		continue
	asteroids.append(i)

#replaces arbitrary number of arguments with asteroid names
astnames=[]
for i in range(0, n):
	astnames.append(os.path.basename(asteroids[i]))

#loads data files 
data=[]
RAdec=[]
for i in range(0, n):
	data.append(pd.read_csv(asteroids[i]+".photast", sep = '\s+', header = None))
for i in range(0,n):
	RAdec.append(pd.read_csv(asteroids[i]+"RAdec.txt", sep = '\s+', header = None))

#merge asteroid data
jointData=[]
filenames=[]
for i in range(0,n):
	filenames.append(asteroids[i]+'.photast')
with open('jointData.photast', 'w') as outfile:
	for x in filenames:
		with open(x) as infile:
			for line in infile:
				outfile.write(line)
				
filenames=[]
for i in range(0,n):
	filenames.append(asteroids[i]+'RAdec.txt')
with open('jointDataRAdec.txt', 'w') as outfile:
	for x in filenames:
		with open(x) as infile:
			for line in infile:
				outfile.write(line)
				
with open('jointData.txt', 'w') as file3:
    with open('jointData.photast', 'r') as file1:
        with open('jointDataRAdec.txt', 'r') as file2:
            for line1, line2 in zip(file1, file2):
                #print >>file3, line1.strip(), line2.strip()		#python 2 syntax
                print (line1.strip(), line2.strip(), file=file3)	#python 3 syntax
                
#joint asteroid data
jointData=pd.read_csv("jointData.txt", sep = '\s+', header = None)
#jointData=jointData[jointData.iloc[:,3]<2458666.1]
jointData.iloc[:,-3]-=2458660
                
#quadratic fitting function
def quad(x, a, b, c):
	return a*x*x + b*x + c

#linear fitting function
def lin(x, a, b):
	return a*x+b
	
#sinusoidal fitting function
def sin(x, a, b,c):
	per=86164.0905/24/3600
	#return a*np.sin((t*2*np.pi)/c+b)
	return a*np.sin(x*2*np.pi/c)+b*np.cos(x*2*np.pi/c)
	
def quad_sin(x, a, b, c, d, e, f):
	return quad(x, a, b, c) + sin(x, d, e, f)

#exclude outliers
jointData=jointData[(jointData.iloc[:,6]-jointData.iloc[:,-9])**2+(jointData.iloc[:,7]-jointData.iloc[:,-8])**2<jointData.iloc[:,-7]**2]

#data points for joint asteroids
jointT=jointData.iloc[:,-3]
jointX=jointData.iloc[:,-6]
jointY=jointData.iloc[:,-5]


########### RA PLOT 1 ###########

fig = plt.figure(figsize=(18, 10))
#fig = plt.figure()
plt.subplot(2,2,1)

#plot points for joint asteroids
plt.plot(jointT, jointX, 'ko', markersize = 4, label='_nolegend_')

#initial guess for quadratic fit
averageX=np.mean(jointT)
slope=(jointData.iloc[-1,-6]-jointData.iloc[0,-6])/(jointData.iloc[-1,-3]-jointData.iloc[0,-3])
aValue=0.1*slope
p0=[aValue, slope, averageX]

#quadratic fit
popt, pcov = curve_fit(quad, jointT, jointX, p0=p0)
print(popt)
plt.plot(jointT, quad(jointT, *popt), c='dodgerblue', linestyle='-', label='Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt))

#linear fit
poptl, pcovl = curve_fit(lin, jointT, jointX) 
print(poptl)
plt.plot(jointT, lin(jointT, *poptl), c='orange', linestyle='--', label='Linear Fit: a=%5.4f, b=%5.4f' % tuple(poptl))

#initial guess for quad_sin
p0=[popt[0]+0, popt[1]+0, popt[2]+0, 0.0003, 0.0004, 1.0]
#p0=[popt[0], popt[1], popt[2], 0.0007, -108.0, 1.00005]
#p0=[popt[0], popt[1], popt[2], 0.0005, 1.06, 0.07]

#time values evenly spaced for plotting
tplot=np.linspace(np.min(jointT), np.max(jointT), 1000)

#quad_sin fit
popts, pcovs = curve_fit(quad_sin, jointT, jointX, p0)
print(popts)
#plt.plot(jointT, quad_sin(jointT, *popts), c='red', linestyle=':')
plt.plot(tplot, quad_sin(tplot, *popts), c='red', linestyle=':', label='Sinusoidal Fit: a=%5.4f, b=%5.4f, c=%5.4f, d=%5.4f, e=%5.4f, f=%5.4f' % tuple(popts))

#legend
#plt.legend(('Data', 'Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt), 'Linear Fit: a=%5.4f, b=%5.4f' % tuple(poptl), 'Sinusoidal Fit: a=%5.4f, b=%5.4f, c=%5.4f, d=%5.4f, e=%5.4f, f=%5.4f' % tuple(popts)))
plt.legend()
plt.title('RA Fit')
plt.xlabel('Time')
plt.ylabel('Right Ascension')

midtime=np.median(jointT)
print(midtime)
print(popts[1:3])
#print("Parameters and Uncertainties:")
#print(popts)
#print(pcovs)
#print(np.sqrt(np.diagonal(pcovs)))


########### RA PLOT 2 ###########

plt.subplot(2,2,3)

#linear component
lplot=lin(tplot, *poptl)
#lplot=lin(tplot, *(popts[1:3]))
l=lin(jointT, *poptl)
#l=lin(jointT, *(popts[1:3]))

#plots data points minus linear 
plt.plot(jointT, jointX-l, 'ko', markersize = 4, label='_nolegend_')

#plots quadratic function minus linear
#plt.plot(jointT, quad(jointT, *popt)-l, label='Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt))	
plt.plot(tplot, quad(tplot, *popt)-lplot, label='Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt))		

#plots quad_sin minus linear component
plt.plot(tplot, quad_sin(tplot, *popts)-lplot, c='red', linestyle='-', label='Sinusoidal Fit: a=%5.4f, b=%5.4f, c=%5.4f, d=%5.4f, e=%5.4f, f=%5.4f' % tuple(popts))

#plots initial guess for quad_sin
#plt.plot(tplot, quad_sin(tplot, *p0)-lplot, c='g', linestyle=':')

#legend
#plt.legend(('Data', 'Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt), 'Sinusoidal Fit: a=%5.4f, b=%5.4f, c=%5.4f, d=%5.4f, e=%5.4f, f=%5.4f' % tuple(popts)))
plt.legend()
plt.title('Linear Component Subtracted (RA)')
plt.xlabel('Time')
plt.ylabel('RA minus linear')


########### DEC PLOT 1 ###########

plt.subplot(2,2,2)

#plot points for joint asteroids
plt.plot(jointT, jointY, 'ko', markersize = 4, label='_nolegend_')

#initial guess for quadratic fit
averageX=np.mean(jointT)
slope=(jointData.iloc[-1,-5]-jointData.iloc[0,-5])/(jointData.iloc[-1,-3]-jointData.iloc[0,-3])
aValue=0.1*slope
p0=[aValue, slope, averageX]

#quadratic fit
popt, pcov = curve_fit(quad, jointT, jointY, p0=p0)
print(popt)
plt.plot(jointT, quad(jointT, *popt), c='dodgerblue', linestyle='-', label='Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt))

#linear fit
poptl, pcovl = curve_fit(lin, jointT, jointY) 
print(poptl)
plt.plot(jointT, lin(jointT, *poptl), c='orange', linestyle='--', label='Linear Fit: a=%5.4f, b=%5.4f' % tuple(poptl))

#initial guess for quad_sin
p0=[popt[0]+0, popt[1]+0, popt[2]+0, 0.0003, 0.0004, 1.0]
#p0=[popt[0], popt[1], popt[2], 0.0007, -108.0, 1.00005]
#p0=[popt[0], popt[1], popt[2], 0.0005, 1.06, 0.07]

#time values evenly spaced for plotting
tplot=np.linspace(np.min(jointT), np.max(jointT), 1000)

#quad_sin fit
popts, pcovs = curve_fit(quad_sin, jointT, jointY, p0)
print(popts)
#plt.plot(jointT, quad_sin(jointT, *popts), c='red', linestyle=':')
plt.plot(tplot, quad_sin(tplot, *popts), c='red', linestyle=':', label='Sinusoidal Fit: a=%5.4f, b=%5.4f, c=%5.4f, d=%5.4f, e=%5.4f, f=%5.4f' % tuple(popts))

#legend
#plt.legend(('none', 'Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt), 'Linear Fit: a=%5.4f, b=%5.4f' % tuple(poptl), 'Sinusoidal Fit: a=%5.4f, b=%5.4f, c=%5.4f, d=%5.4f, e=%5.4f, f=%5.4f' % tuple(popts)))
plt.legend()
plt.title('Dec Fit')
plt.xlabel('Time')
plt.ylabel('Declination')


midtime=np.median(jointT)
print(midtime)
print(popts[1:3])
#print("Parameters and Uncertainties:")
#print(popts)
#print(pcovs)
#print(np.sqrt(np.diagonal(pcovs)))


########### DEC PLOT 2 ###########

plt.subplot(2,2,4)

#linear component
lplot=lin(tplot, *poptl)
#lplot=lin(tplot, *(popts[1:3]))
l=lin(jointT, *poptl)
#l=lin(jointT, *(popts[1:3]))

#plots data points minus linear 
plt.plot(jointT, jointY-l, 'ko', markersize = 4, label='_nolegend_')

#plots quadratic function minus linear
#plt.plot(jointT, quad(jointT, *popt)-l, label='Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt))	
plt.plot(tplot, quad(tplot, *popt)-lplot, label='Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt))		

#plots quad_sin minus linear component
plt.plot(tplot, quad_sin(tplot, *popts)-lplot, c='red', linestyle='-', label='Sinusoidal Fit: a=%5.4f, b=%5.4f, c=%5.4f, d=%5.4f, e=%5.4f, f=%5.4f' % tuple(popts))

#plots initial guess for quad_sin
#plt.plot(tplot, quad_sin(tplot, *p0)-lplot, c='g', linestyle=':')

#convert linear parameters to strings
#poptl =["%.8f" % i for i in poptl]

#prints linear parameters into a text file
#f=open("/mnt/c/Users/marzt/Documents/Research/linear_params_Dec.txt", "a+")
#f.write(ast1+'_'+ast2name+'\t'+poptl[0]+'\t'+poptl[1]+'\n')
#f.close()

#legend
#plt.legend(('none', 'Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt), 'Sinusoidal Fit: a=%5.4f, b=%5.4f, c=%5.4f, d=%5.4f, e=%5.4f, f=%5.4f' % tuple(popts)))
plt.legend()
plt.title('Linear Component Subtracted (Dec)')
plt.xlabel('Time')
plt.ylabel('Dec minus linear')


#save/show plot
fig.suptitle(astnames[0]+'_'+astnames[1]+'_'+astnames[2]+' RA and Dec Plots')
plt.savefig(astnames[0]+'_'+astnames[1]+'_'+astnames[2]+'_linkast.pdf')
plt.show()

