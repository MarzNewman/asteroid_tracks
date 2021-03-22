import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import sys
import os
import math
import scipy.stats as stats
from matplotlib.backends.backend_pdf import PdfPages

if len(sys.argv)!=3:
	print("usage: %s <asteroidname> <iteration>" % (sys.argv[0]))
	sys.exit()


astname=sys.argv[1]
iteration=sys.argv[2]

#############################
########## OPTIONS ##########
#############################

sigthresh=2

#############################


astdir=os.getcwd()
print(astdir)
astdir = os.path.basename(astdir) 
#curdir=str(curdir)
print(astdir)

if iteration=="1":
	data = pd.read_csv(astname+"_1RAdec.photast", sep = '\s+', header = None)
elif iteration=="2":	
	data = pd.read_csv(astname+"_2RAdec.photast", sep = '\s+', header = None)
	
calib = pd.read_csv('/mnt/c/Users/marzt/Documents/Research/MISHAPS_F1_bothcal.txt', sep = '\s+', header = None)

#exclude outliers
data=data[(data.iloc[:,6]-data.iloc[:,-5])**2+(data.iloc[:,7]-data.iloc[:,-4])**2<data.iloc[:,-3]**2]


###############################
########## FUNCTIONS ##########
###############################

#quadratic fitting function
def quad(x, a, b, c):
	midtime=np.median(data.iloc[:,0])
	t=x-midtime
	return a*t*t + b*t + c
	
#linear fitting function
def lin(x, a, b):
	#midtime=np.median(x)
	#t=x-midtime
	midtime=np.median(data.iloc[:,0])
	t=x-midtime
	#t=x-2458660
	return a*t+b
	
#sinusoidal fitting function
def sin(x, a, b, c):
	per=86164.0905/24/3600
	#return a*np.sin((t*2*np.pi)/c+b)
	midtime=np.median(data.iloc[:,0])
	t=x-midtime
	return a*np.sin(t*2*np.pi/c)+b*np.cos(t*2*np.pi/c)
	
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


################################
########## X TIME FIT ##########
################################

fig = plt.figure(figsize=(18, 10))
plt.figure()
ax=plt.subplot(2,1,1)

#x and y data points
tdata = data.iloc[:, 0]
xdata = data.iloc[:, 6]

#plot points
plt.plot(tdata, xdata, 'ko', label = 'data')

#initial guess
averageX=np.mean(xdata)
slope=(data.iloc[-1,6]-data.iloc[0,6])/(data.iloc[-1,0]-data.iloc[0,0])
aValue=0.1*slope
p0=[aValue, slope, averageX]

#perform quadratic curve-fit
popt, pcov = curve_fit(quad, tdata, xdata, p0)
print(popt)

#plot the fitted quadratic function
plt.plot(tdata, quad(tdata, *popt), 'r-', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

#perform linear curve-fit
poptl, pcovl = curve_fit(lin, tdata, xdata)
print(poptl)

#calculates linear parameter error bars
errorl=np.sqrt(np.diagonal(pcovl))

#plots points minus linear component
l=lin(tdata, *poptl)
#plt.plot(tdata, xdata-l, 'bo')

#plots quadratic function minus linear
#plt.plot(tdata, quad(tdata, *popt)-l, label='Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt))

#creates a list of optimized positions
opt_position=[]
n=len(tdata)
if iteration=="1":
	for i in tdata:
		#opt_position.append(np.interp(i, tdata, xdata))
		#opt_position.append(popt[0]*i*i + popt[1]*i + popt[2])
		opt_position.append(quad(i, *popt))				#for quadratic fit without subtracting linear component
	#for i in range(0,n):	
		#opt_position.append(quad(tdata[i], *popt)-l[i])		#for quadratic fit minus linear component
	
#converts optimized positions to strings
if iteration=="1":
	opt_position = ["%.8f" % i for i in opt_position]

#creates a text file with optimized positions using the model
if iteration=="1":
	f=open(astname+'_optx.txt','w')
	for element in opt_position:
	     f.write(element)
	     f.write('\n')
	f.close()

#initial guess for sin fit
#p0=[popt[0]+0, popt[1]+0, popt[2]+0, 0.0003, 0.0004, 1.0]
p0=[popt[0]+0, popt[1]+0, popt[2]+0, 0.0005, 0.0005, 1.0]

#perform quad_sin fitting function
#popts, pcovs = curve_fit(quad_sin, tdata, xdata, p0)
#popts, pcovs = curve_fit(quad_sin, tdata, xdata)

#plot the fitted quad_sin function
#plt.plot(tdata, quad_sin(tdata, *popts), 'b:', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

#calculates sinusoidal parameter error bars
#errors=np.sqrt(np.diagonal(pcovs))
#errors=np.sqrt(abs(np.diagonal(pcovs)))

#converts linear fit pixels to arcseconds
poptl[0]=poptl[0]*0.2637
poptl[1]=poptl[1]*0.2637
errorl[0]=errorl[0]*0.2637
errorl[1]=errorl[1]*0.2637

#converts linear parameters and error bars to strings
poptl = ["%.8f" % i for i in poptl]
errorl = ["%.8f" % i for i in errorl]

#creates a text file for each asteroid in the linear parameters table
f=open("/mnt/c/Users/marzt/Documents/Research/"+astname+"linparams_"+iteration+"x.txt", "w+")
f.write(astname+'\t'+poptl[0]+'\t'+poptl[1]+'\t'+errorl[0]+'\t'+errorl[1]+'\n')
f.close()

#converts sin parameters and error bars to strings
#popts = ["%.8f" % i for i in popts]
#errors = ["%.8f" % i for i in errors]

#creates a text file for each asteroid in the sin parameters table
#f=open("/mnt/c/Users/marzt/Documents/Research/"+astname+"sinparams_"+iteration+"x.txt", "w+")
#f.write(astname+'\t'+popts[0]+'\t'+popts[1]+'\t'+popts[2]+'\t'+popts[3]+'\t'+popts[4]+'\t'+popts[5]+'\t'+errors[0]+'\t'+errors[1]+'\t'+errors[2]+'\t'+errors[3]+'\t'+errors[4]+'\t'+errors[5]+'\n')
#f.close()

#prints linear parameters into a text file
#f=open("/mnt/c/Users/marzt/Documents/Research/linear_params_x.txt", "a+")
#f.write(astname+'\t'+poptl[0]+'\t'+poptl[1]+'\n')
#f.close()

#creates an array of residuals
#n=len(xdata)
#residuals=[]
#for i in range(0,n):
#	residuals.append(xdata[i]-quad(tdata[i], *popt))
residuals=xdata-quad(tdata, *popt)

#passes residuals through median absolute deviation function
series = pd.Series(residuals)
mad=series.mad()
badx=np.abs(residuals)>sigthresh*mad
outdata=data.copy()
outdata.insert(len(outdata.columns), 'badx'+iteration, badx)
print(mad)

plt.xlabel('time')
plt.ylabel('x')
plt.legend()
#plt.savefig(astname+'timefitx.pdf')
#plt.show()


################################
########## Y TIME FIT ##########
################################

ay=plt.subplot(2,1,2)

#x and y data points
tdata = data.iloc[:, 0]
ydata = data.iloc[:, 7]

#plot points
#with PdfPages(astname+'_'+iteration+'all_fits.pdf') as pdf:

plt.plot(tdata, ydata, 'ko', label = 'data')
timefity=plt.plot(tdata, ydata, 'ko', label = 'data')
#pdf.savefig()
#plt.close()

#initial guess
averageX=np.mean(ydata)
slope=(data.iloc[-1,7]-data.iloc[0,7])/(data.iloc[-1,0]-data.iloc[0,0])
aValue=0.1*slope
p0=[aValue, slope, averageX]

#perform quadratic curve-fit
popty, pcovy = curve_fit(quad, tdata, ydata, p0)
print(popty)

#plot the fitted quadratic function
plt.plot(tdata, quad(tdata, *popty), 'r-', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popty))

#perform linear curve-fit
poptly, pcovly = curve_fit(lin, tdata, ydata)
print(poptly)

#calculates linear parameter error bars
errorly=np.sqrt(np.diagonal(pcovly))

#plots points minus linear component
l=lin(tdata, *poptly)
#plt.plot(tdata, ydata-l, 'bo')

#plots quadratic function minus linear
#plt.plot(tdata, quad(tdata, *popt)-l, label='Quadratic Fit: a=%5.4f, b=%5.4f, c=%5.4f' % tuple(popt))

#creates a list of optimized positions
opt_position=[]
n=len(tdata)
if iteration=="1":
	for i in tdata:
		#opt_position.append(np.interp(i, tdata, ydata))
		#opt_position.append(popt[0]*i*i + popt[1]*i + popt[2])
		opt_position.append(quad(i, *popty))				#for quadratic fit without subtracting linear component
	#for i in range(0,n):	
		#opt_position.append(quad(tdata[i], *popt)-l[i])		#for quadratic fit minus linear component

#converts optimized positions to strings
if iteration=="1":
	opt_position = ["%.8f" % i for i in opt_position]
	
#print('optimization')
#print(opt_position)

#creates a text file with optimized positions using the model
if iteration=="1":
	f=open(astname+'_opty.txt','w')
	for element in opt_position:
	     f.write(element)
	     f.write('\n')
	f.close()

#initial guess for sin fit
p0=[popty[0]+0, popty[1]+0, popty[2]+0, 0.0003, 0.0004, 1.0]

#perform quad_sin fitting function
#popts, pcovs = curve_fit(quad_sin, tdata, ydata, p0)

#plot the fitted quad_sin function
#plt.plot(tdata, quad_sin(tdata, *popts), 'b:', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

#calculates sinusoidal parameter error bars
#errors=np.sqrt(np.diagonal(pcovs))
#errorsy=np.sqrt(abs(np.diagonal(pcovsy)))

#converts linear fit pixels to arcseconds
poptly[0]=poptly[0]*0.2637
poptly[1]=poptly[1]*0.2637
errorly[0]=errorly[0]*0.2637
errorly[1]=errorly[1]*0.2637

#converts linear parameters and error bars to strings
poptly = ["%.8f" % i for i in poptly]
errorly = ["%.8f" % i for i in errorly]

#creates a text file for each asteroid in the linear parameters table
f=open("/mnt/c/Users/marzt/Documents/Research/"+astname+"linparams_"+iteration+"y.txt", "w+")
f.write(astname+'\t'+poptly[0]+'\t'+poptly[1]+'\t'+errorly[0]+'\t'+errorly[1]+'\n')
f.close()

#converts sin parameters and error bars to strings
#popts = ["%.8f" % i for i in popts]
#errors = ["%.8f" % i for i in errors]

#creates a text file for each asteroid in the sin parameters table
#f=open("/mnt/c/Users/marzt/Documents/Research/"+astname+"sinparams_"+iteration+"y.txt", "w+")
#f.write(astname+'\t'+popts[0]+'\t'+popts[1]+'\t'+popts[2]+'\t'+popts[3]+'\t'+popts[4]+'\t'+popts[5]+'\t'+errors[0]+'\t'+errors[1]+'\t'+errors[2]+'\t'+errors[3]+'\t'+errors[4]+'\t'+errors[5]+'\n')
#f.close()

#prints linear parameters into a text file
#f=open("/mnt/c/Users/marzt/Documents/Research/linear_params_y.txt", "a+")
#f.write(astname+'\t'+poptl[0]+'\t'+poptl[1]+'\n')
#f.close()

#creates an array of residuals
#n=len(ydata)
#residuals=[]
#for i in range(0,n):
#	residuals.append(ydata[i]-quad(tdata[i], *popt))
residualsy=ydata-quad(tdata, *popty)

#print(residuals)
	
#passes residuals through median absolute deviation function
#dev=[]
#for i in residuals:
	#dev.append(stats.median_absolute_deviation(i))

series = pd.Series(residualsy)
mad=series.mad()
bady=np.abs(residualsy)>sigthresh*mad

notoutliers=np.logical_not(badx | bady)
outdata.insert(len(outdata.columns), 'bady'+iteration, bady)
outdata.to_csv(astname+'_'+iteration+'outliers.txt', sep=' ', header=False)
	
plt.plot(tdata[notoutliers], ydata[notoutliers], 'bo', label = 'outliers')
ax.plot(tdata[notoutliers], xdata[notoutliers], 'bo', label = 'outliers')
plt.close

print(badx, bady)	
	
#print(mad)

plt.xlabel('time')
plt.ylabel('y')
plt.legend()
#plt.savefig(astname+'timefity.pdf')
#plt.show()

plt.savefig(astname+'_'+iteration+'xytimefits.pdf')


###############################################
########## MERGE OPTIMIZED POSITIONS ##########
###############################################

if iteration=="1":
	with open(astname+'_deleteme.txt', 'w') as file3:
		with open(astname+'_optx.txt', 'r') as file1:
        		with open(astname+'_opty.txt', 'r') as file2:
            			for line1, line2 in zip(file1, file2):
                			print >> file3, line1.strip(), line2.strip()		#python 3 syntax
                			#print (line1.strip(), line2.strip(), file=file3)	#python 2 syntax
   

################################
########## PHOTOMETRY ##########
################################

plt.figure()

#defines flux data
fluxData=data.iloc[:,2]
tData=data.iloc[:,0]

#passes directory through the search_string_in_file function
matched_line = search_string_in_file('/mnt/c/Users/marzt/Documents/Research/MISHAPS_F1_bothcal.txt', astdir)
matched_line=int(matched_line[0])

#exclude outliers
#zscore=stats.zscore(fluxData)
#data=data[zscore<1.25]
#fluxData=data.iloc[:,2]
#zscore=stats.zscore(fluxData)
#data=data[zscore>-1.25]
#fluxData=data.iloc[:,2]
#tData=data.iloc[:,0]p

#converts flux to magnitude
row=matched_line-1
c=calib.iloc[row,1]
d=calib.iloc[row,3]
photData=(c-d-2.5*np.log10(-fluxData))
#photData=(c-d-2.5*np.log10(abs(fluxData)))

#converts uncertainties to magnitude
num=len(fluxData)
#unc=[]
#for i in range(0,num):
#	unc.append(data.iloc[i,3]*-2.5/(np.log(10)*abs(data.iloc[i,2])))
unc=(data.iloc[:,3]*2.5/(np.log(10)*abs(data.iloc[:,2])))
	
#plots photometry data 	
plt.errorbar(tData, photData, fmt='ko', yerr=unc)
plt.errorbar(tData[notoutliers], photData[notoutliers], fmt="bo", yerr=unc[notoutliers])
medianmag=np.nanmedian(photData[notoutliers])
plt.ylim([medianmag+.5,medianmag-.5])				#creates a range for the y-axis of the magnitude plot
					

#calculates mean magnitude and standard deviation of the magnitude
mean_mag=str(np.mean(photData))
std_mag=str(np.std(photData))

#appends mean magnitude and standard deviation to linparams_x file
f=open("/mnt/c/Users/marzt/Documents/Research/"+astname+"linparams_"+iteration+"x.txt", "a")
f.write(mean_mag+'\t'+std_mag+'\n')
f.close()

#appends mean magnitude and standard deviation to linparams_y file
f=open("/mnt/c/Users/marzt/Documents/Research/"+astname+"linparams_"+iteration+"y.txt", "a")
f.write(mean_mag+'\t'+std_mag+'\n')
f.close()

plt.savefig(astname+'_'+iteration+'photom.pdf')
plt.close()
#plt.show()

#creates multi-page pdf with all plots
with PdfPages(astname+'_'+iteration+'all_fits.pdf') as pdf:
	#time fits for x and y
	plt.figure('Position-Time')
	plt.subplot(2,1,1)
	plt.plot(tdata, xdata, 'bo')
	plt.plot(tdata[notoutliers], xdata[notoutliers], 'ko')
	plt.title('Position-Time ('+astname+'_'+iteration+')')
	plt.subplot(2,1,2)
	plt.plot(tdata, ydata, 'bo')
	plt.plot(tdata[notoutliers], ydata[notoutliers], 'ko')
	pdf.savefig()
	plt.close()
	
	plt.figure('Residuals')
	plt.subplot(2,1,1)
	plt.plot(tdata, residuals, 'bo')
	plt.plot(tdata[notoutliers], residuals[notoutliers], 'ko', label='%6.3f'%(np.std(residuals[notoutliers])))
	plt.title('Residuals ('+astname+'_'+iteration+')')
	plt.legend()
	plt.subplot(2,1,2)
	plt.plot(tdata, residualsy, 'bo')
	plt.plot(tdata[notoutliers], residualsy[notoutliers], 'ko', label='%6.3f'%(np.std(residualsy[notoutliers])))
	plt.legend()
	pdf.savefig()
	plt.close()
	
	#photometry
	plt.figure('Magnitude Photometry')
	plt.errorbar(tData, photData, fmt='bo', yerr=unc)
	plt.errorbar(tData[notoutliers], photData[notoutliers], fmt='ko', yerr=unc[notoutliers])
	plt.title('Magnitude Photometry ('+astname+'_'+iteration+')')
	pdf.savefig()
	plt.close()
	

