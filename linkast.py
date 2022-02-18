import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import sys
import os
import math
import scipy.stats as stats

#manager = plt.get_current_fig_manager()
#manager.resize(*manager.window.maxsize())

if len(sys.argv)!=3:
	print("usage: %s <asteroid1name> <asteroid2name>" % (sys.argv[0]))
	sys.exit()

ast1=sys.argv[1]
ast2=sys.argv[2]
ast2name=os.path.basename(ast2)

astdir1=os.getcwd()
astdir1 = os.path.basename(astdir1) 
astdir2=os.getcwd()
astdir2 = os.path.basename(astdir2)

calib = pd.read_csv('/mnt/c/Users/marzt/Documents/Research/MISHAPS_F1_bothcal.txt', sep = '\s+', header = None)

data = pd.read_csv(ast1+".photast", sep = '\s+', header = None)
data2 = pd.read_csv(ast2+".photast", sep = '\s+', header = None)

RAdec1 = pd.read_csv(ast1+"RAdec.txt", sep = '\s+', header = None)
RAdec2 = pd.read_csv(ast2+"RAdec.txt", sep = '\s+', header = None)

#merge asteroid data
filenames=[ast1+'.photast', ast2+'.photast']
with open('jointData.photast', 'w') as outfile:
	for x in filenames:
		with open(x) as infile:
			for line in infile:
				outfile.write(line)
				
filenames=[ast1+'RAdec.txt', ast2+'RAdec.txt']
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
                
                
########### RA PLOT 1 ###########

fig = plt.figure(figsize=(18, 10))
plt.subplot(3,2,1)

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
	
#function that searches for a string in a file
def search_string_in_file(file_name, string_to_search):
	line_number = 0
	line_with_data=[]
	with open(file_name, 'r') as read_obj:
		for line in read_obj:
			line_number += 1
			if string_to_search in line:
				line_with_data.append(line_number)
	return line_with_data

#exclude outliers
jointData=jointData[(jointData.iloc[:,6]-jointData.iloc[:,-9])**2+(jointData.iloc[:,7]-jointData.iloc[:,-8])**2<jointData.iloc[:,-7]**2]

#data points for joint asteroids
jointT=jointData.iloc[:,-3]
jointX=jointData.iloc[:,-6]
jointY=jointData.iloc[:,-5]

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

plt.subplot(3,2,3)

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

#sinusoidal fitting function
def sin2(x, a, b,c):
	per=86164.0905/24/3600
	return a*np.sin((x*2*np.pi)/c+b)
	#return a*np.sin(x*2*np.pi/c)+b*np.cos(x*2*np.pi/c)
	
def quad_sin2(x, a, b, c, d, e, f):
	return quad(x, a, b, c) + sin2(x, d, e, f)

#performs quad_sin2 fit	
popts2, pcovs2 = curve_fit(quad_sin2, jointT, jointX, p0)

#plt.plot(tplot, quad_sin2(tplot, *popts2)-lplot, c='green', linestyle='-') 

#maximum parallax amplitude
amplitude=popts2[3]

#calculates sinusoidal parameter error bars
errors=np.sqrt(np.diagonal(pcovs2))

#calculate distance to asteroid
unc_parallax=errors[3]*3600
radius=6378000
#baseline=radius/149600000000
baseline=0.0000426342967
unc_distance=baseline/unc_parallax
#latitude=30.169*(np.pi/180)
#parallax=amplitude*(np.pi/180)
parallax=amplitude*3600
distance=(baseline/parallax)*206265
#distance=(radius*2/parallax)/149600000000
#distance=(radius*np.cos(latitude)/parallax)/149600000000

#converts variables to strings
asteroid1 = str(ast1)
asteroid2 = str(ast2name)
parallax = str(parallax)
distance = str(distance)
unc_parallax = str(unc_parallax)
unc_distance = str(unc_distance)

#creates a text file for each asteroid 
f=open("/mnt/c/Users/marzt/Documents/Research/"+asteroid1+"_"+asteroid2+"distance_RA.txt", "w+")
f.write(ast1+'\t'+ast2name+'\t'+parallax+'\t'+distance+'\t'+unc_parallax+'\t'+unc_distance+'\n')
f.close()

#legend
#plt.legend(('Data', 'Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt), 'Sinusoidal Fit: a=%5.4f, b=%5.4f, c=%5.4f, d=%5.4f, e=%5.4f, f=%5.4f' % tuple(popts)))
plt.legend()
plt.title('Linear Component Subtracted (RA)')
plt.xlabel('Time')
plt.ylabel('RA minus linear')


########### PHOT PLOT 1 ###########

plt.subplot(3,2,5)

photdata1 = pd.read_csv(ast1+"_2.photast", sep = '\s+', header = None)
	
#defines flux data
fluxData=photdata1.iloc[:,2]
tData=photdata1.iloc[:,0]

#passes directory through the search_string_in_file function
matched_line = search_string_in_file('/mnt/c/Users/marzt/Documents/Research/MISHAPS_F1_bothcal.txt', astdir1)
matched_line=int(matched_line[0])

#exclude outliers
zscore=stats.zscore(fluxData)

photdata1=photdata1[zscore<1.25]

fluxData=photdata1.iloc[:,2]
tData=photdata1.iloc[:,0]

zscore=stats.zscore(fluxData)

photdata1=photdata1[zscore>-1.25]

#defines data excluding outliers
fluxData=photdata1.iloc[:,2]
tData=photdata1.iloc[:,0]

#converts flux to magnitude
row=matched_line-1
c=calib.iloc[row,1]
d=calib.iloc[row,3]
#photData=[]
#for i in fluxData:
	#photData.append(c-d-2.5*math.log10(-i))
photData=(c-d-2.5*np.log10(-fluxData))

#defines data excluding outliers
#photData=photdata1.iloc[:,2]
#tData=photdata1.iloc[:,0]

#plots photometry data
#plt.errorbar(tData, photData, fmt='ko', yerr=photdata1.iloc[:,3])
plt.errorbar(tData, photData, fmt='ko')

#initial guess
#averageX=np.mean(tData)
#slope=(photdata1.iloc[-1,2]-photdata1.iloc[0,2])/(photdata1.iloc[-1,0]-photdata1.iloc[0,0])
#aValue=0.1*slope
#p0=[aValue, slope, averageX]

#quadratic fit
#popt, pcov = curve_fit(func, tData, photData, p0)
#print(popt)

#initial guess for sin_func fit
#p0=[popt[0], popt[1], popt[2], 200, 0.07, 0.0005]

#sin_func fit
#popts, pcovs = curve_fit(sin_func, tData, photData, p0)
#print(popts)
#tplot=np.linspace(np.min(tData), np.max(tData), 1000)
#plt.plot(tplot, sin_func(tplot, *popts), c='red', linestyle='-')

#plots sin_func initial guess
#plt.plot(tplot, sin_func(tplot, *p0), c='g', linestyle=':')

plt.title('Asteroid 1 Photometry')
plt.xlabel('Time')
plt.ylabel('Magnitude')
#plt.ylim([20, 22])


########### DEC PLOT 1 ###########

plt.subplot(3,2,2)

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
#plt.legend(('Data', 'Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt), 'Linear Fit: a=%5.4f, b=%5.4f' % tuple(poptl), 'Sinusoidal Fit: a=%5.4f, b=%5.4f, c=%5.4f, d=%5.4f, e=%5.4f, f=%5.4f' % tuple(popts)))
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

plt.subplot(3,2,4)

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
#plt.legend(('Data', 'Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt), 'Sinusoidal Fit: a=%5.4f, b=%5.4f, c=%5.4f, d=%5.4f, e=%5.4f, f=%5.4f' % tuple(popts)))
plt.legend()
plt.title('Linear Component Subtracted (Dec)')
plt.xlabel('Time')
plt.ylabel('Dec minus linear')


########### PHOT PLOT 2 ###########

plt.subplot(3,2,6)

photdata2 = pd.read_csv(ast2+"_2.photast", sep = '\s+', header = None)
	
#defines data
fluxData=photdata2.iloc[:,2]
tData=photdata2.iloc[:,0]

#passes directory through the search_string_in_file function
matched_line = search_string_in_file('/mnt/c/Users/marzt/Documents/Research/MISHAPS_F1_bothcal.txt', astdir2)
matched_line=int(matched_line[0])

#exclude outliers
zscore=stats.zscore(fluxData)

photdata2=photdata2[zscore<1.25]

fluxData=photdata2.iloc[:,2]
tData=photdata2.iloc[:,0]

zscore=stats.zscore(fluxData)

photdata2=photdata2[zscore>-1.25]

#defines data excluding outliers
fluxData=photdata2.iloc[:,2]
tData=photdata2.iloc[:,0]

#converts flux to magnitude
row=matched_line-1
c=calib.iloc[row,1]
d=calib.iloc[row,3]
#photData=[]
#for i in fluxData:
	#photData.append(c-d-2.5*math.log10(abs(i)))
photData=(c-d-2.5*np.log10(-fluxData))

#defines data excluding outliers
#photData=photdata2.iloc[:,2]
#tData=photdata2.iloc[:,0]

#plots photometry data
#plt.errorbar(tData, photData, fmt='ko', yerr=photdata2.iloc[:,3])
plt.errorbar(tData, photData, fmt='ko')

#initial guess
#averageX=np.mean(tData)
#slope=(photdata2.iloc[-1,2]-photdata2.iloc[0,2])/(photdata2.iloc[-1,0]-photdata2.iloc[0,0])
#aValue=0.1*slope
#p0=[aValue, slope, averageX]

#quadratic fit
#popt, pcov = curve_fit(func, tData, photData, p0)
#print(popt)

#initial guess for sin_func fit
#p0=[popt[0], popt[1], popt[2], 200, 0.07, 0.0005]

#sin_func fit
#popts, pcovs = curve_fit(sin_func, tData, photData, p0)
#print(popts)
#tplot=np.linspace(np.min(tData), np.max(tData), 1000)
#plt.plot(tplot, sin_func(tplot, *popts), c='red', linestyle='-')

#plots sin_func initial guess
#plt.plot(tplot, sin_func(tplot, *p0), c='g', linestyle=':')

plt.title('Asteroid 2 Photometry')
plt.xlabel('Time')
plt.ylabel('Magnitude')
#plt.ylim([20, 22])



#save/show plot
fig.suptitle(ast1+'_'+ast2name+' RA and Dec Plots')
plt.savefig(ast1+'_'+ast2name+'_linkast.pdf')
#plt.ylim([mean-2000, mean+2000])
plt.show()

