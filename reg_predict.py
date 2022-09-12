import numpy as np
import pandas as pd
import re 
import sys
import os
import glob
from os.path import exists
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from all_fits.py import residuals

print("\n\n\nRunning reg_predict.py")
print("To run this program successfully, you need a predict.reg regfile with the first three regions of a night")

if len(sys.argv)!=4:
	print("usage: %s <astname> <numNewRegions> <print new region to MNpredict.reg? [True or False]>" % (sys.argv[0]))
	sys.exit()

#rms_residuals = np.sqrt(np.mean(residuals**2))
#print("rms_residuals:", rms_residuals)

astname=sys.argv[1]
#print(astname)
numNewRegions=int(sys.argv[2])	#defines the number of new regions to produce in the script
#print(numNewRegions)
#print("numNewRegions:", numNewRegions)
print2reg=sys.argv[3]
#print(print2reg)
regfile="/project/marz746/astdata/regfiles/" + astname + ".reg"
#print(regfile)
#predictFile=astname + ".reg"
#imsep=5
headfile="/project/marz746/astdata/regfiles/MN_asteroid_tracking.txt"
#print(headfile)

#load asteroid_tracking header
header = pd.read_csv(headfile, sep = '\s+', header = 0)
#print(header)

#defines residuals
#resfile=open("/project/marz746/astdata/output/"+astname+"_2xresiduals.txt", "r")
########resfile="/project/marz746/astdata/output/"+astname+"_2xresiduals.txt"
#print(resfile)
#residuals=resfile.readlines()
##residuals=resfile.read()
##residuals=residuals.split("\n")
##residuals=residuals[:-1]
############residuals=np.loadtxt(resfile)
#print(residuals)
#residuals=np.asarray(residuals)
#print("residuals:", residuals)

#calculates rms of residuals
###############rms_res=np.sqrt(np.mean(residuals**2))
#print(rms_res)
#print("rms of residuals:", rms_res)

#defines photometry (used to find out where the algorithm went wrong)
##photfile="/project/marz746/astdata/output/"+astname+"_1phot.txt"
##print(photfile)
##photdata=np.loadtxt(photfile)
##print(photdata)
#print("phot data:", photdata)

#find line with asteroid data in header
headData = header[header.loc[:, "Asteroid"]==astname]
print(headData)
headData.reset_index(drop=True, inplace=True)
print(headData)
night=headData["Night"]
#print(night[0])
image=headData["First_image"]
#print(image)

imsubstring= "c4d_"+str(night[0])
#image = 'c4d_210706_001008_ooi_r_v1.fits.fz'
holdings = "/project/marz746/astdata/MISHAPS_F1_r.good.holdings"

#checks if astrometry output exists
if os.path.exists("/project/marz746/astdata/output/"+astname+"_predInput.txt") == False:
	print(astname+"_predInput.txt does not exist. Run astrometry on regfile and re-run")
	sys.exit()

#loads in astrometry output file
optoutput = pd.read_csv("/project/marz746/astdata/output/"+astname+"_predInput.txt", sep = '\s+', header = None)
optTimes = optoutput.iloc[:, 2]
optx = optoutput.iloc[:, 0]
#print(optx)
opty = optoutput.iloc[:, 1]

#loads good.holdings file for the timestamp data
timedata = pd.read_csv(holdings, sep = '\s+', header = None)

#finds the number of images taken in a night
#numimages = len(glob.glob('conv*'))

#exits if opt files don't exist
#if (os.path.exists("MNpredict_optx.txt") or os.path.exists("MNpredict_opty.txt")) == False:
#	print("MNpredict_optx.txt and/or MNpredict_opty.txt do(es) not exist. Run astrometry on MNpredict.reg and re-run")
#	sys.exit()

'''
#defines data for optimized positions
optx = np.loadtxt("MNpredict_optx.txt", dtype='float', delimiter='\n')
#print(optx)
opty = np.loadtxt("MNpredict_opty.txt", dtype='float', delimiter='\n')
#print(opty)
'''

'''
#define time index for the first image from the good.holdings 
##timeindex = timedata[timedata.iloc[:,0] == image]
##timeindex.reset_index(drop=True, inplace=True)
##print(timeindex)
timeindex = []
for i, line in enumerate(open(holdings, 'r').readlines()): 
	if image in line:
		#extract line index for lines that contain string
		timeindex.append(i)
'''

#define time indexes in the good.holdings file
holdingsData = timedata[timedata.iloc[:,0].str.contains(imsubstring)]
holdingsData.reset_index(drop=True, inplace=True)
#print(holdingsData)
include=np.zeros(shape=holdingsData.shape[0])+1
#print(holdingsData.iloc[0,0])

#all Times of the holdings data
allTimes=holdingsData.iloc[:,1]

#boolean indexing if the first image for the asteroid is not the same as the first image of the night
#if image[0] != holdingsData.iloc[0,0]:
#print(str(holdingsData.iloc[0,0]).contains(str(image[0])[5:22]))
if str(image[0])[5:22] not in holdingsData.iloc[0,0]:
	print('\nAsteroid is not observed on first image of the night; entering Boolean indexing loop')
	#print(str(image[0])[5:22])
	#matchIdx_tmp = holdingsData[holdingsData.iloc[:,0].str.contains(str(image[0])[5:22])]
	matchIdx = holdingsData.iloc[:,0].str.contains(str(image[0])[5:22])
	#print('fffffff')
	#print(matchIdx_tmp)
	#matchIdx = matchIdx_tmp.index[0]
	#print(matchIdx)
	#boolVector=[]
	started=False
	for i,row in holdingsData.iterrows():
		#print(row)
		#print(row[0])
		if matchIdx[i]==True:
			started=True
		matchIdx[i]=started
		#else:
		#	boolVector.append(True)
	#print(matchIdx_tmp)
	#holdingsData=pd.DataFrame(holdingsData, index=boolVector)
	holdingsData=holdingsData[matchIdx]
	#print(holdingsData)
	#print(holdingsData.loc[True])

holdingsData.reset_index(drop=True, inplace=True)	
times = holdingsData.iloc[:,1]
#times.reset_index(drop=True, inplace=True)
#print(times)
#print(times)

#print(timeindex)
#numTimes = len(timeindex)
#print(numTimes)
#timeindex = timeindex[0]
#print(timeindex)
#timeEnd = times[-1]
#print(timeEnd)

'''
#defines a list of times for the night
times = timedata.iloc[timeindex:timeindex+numimages-2,1]
times.reset_index(drop=True, inplace=True)
#print(times)
print(times[0])
'''

'''
#checks if region file already exists
if os.path.exists(astname + ".reg") == True:
	#print(astname + ".reg already exists; delete and try again")
	#sys.exit()
	print(astname + ".reg already exists; deleting")
	os.remove(astname + ".reg")
'''

#defines data from predict.reg file
regdata = np.loadtxt(regfile, skiprows=3, dtype='str')

'''
#function that searches for a string in a file
def search_string_in_file(file_name, string_to_search):
	line_number = 0
	line_with_data=[]
	with open(file_name, 'r') as read_obj:
		for line in read_obj:
			line_number += 1
			if string_to_search in line:
				line_with_data.append(line_number)
	return 
'''

#defines median time to subtract (makes the initial guesses work better for the fitting)
midtime=np.median(optTimes)

#defines a linear function for the motion of the asteroid
def lin(t, m, b):
	t=t-midtime
	return m*t+b

#sinusoidal function to use for 9 new regions and above
def sin(t, a, b, c):
	return a*np.sin(t*2*np.pi/c)+b*np.cos(t*2*np.pi/c)
	
#quadratic function
def quad(t, a, b, c):
	t=t-midtime
	return a*t*t + b*t + c
	
#function that finds the slope between two positions
def find_slope(x1, x2, y1, y2):
	return (float(y2)-float(y1))/(float(x2)-float(x1))
	
#function that finds y-intercept
def find_yintercept(x, y, slope):
	return float(y)-slope*(float(x))


'''
#defines times to use from the good.holdings file
timeIndex = search_string_in_file("MISHAPS_F1_r.good.holdings", "2458664.51183237")
print("time index:", timeIndex)
'''

#finds initial guesses
#reads in data from regfile
xvals=[]
yvals=[]
for line in regdata:
	#print(line)
	line=line.split(",")
	#print(line)
	x=line[0][7:15]		#defines x position
	xvals.append(x)		#adds x position to list
	#print("x position:", x)
	y=line[1]			#defines y position
	#print("y position:", y)
	yvals.append(y)		#adds x position to list
	regsize=line[2][0:3]		#defines regsize

print("regsize:", regsize)

'''
##makes the MN****.reg file with optimized positions##
#checks if region file already exists
if os.path.exists(astname + ".reg") == True:
	#print(astname + ".reg already exists; delete and try again")
	#sys.exit()
	print(astname + ".reg already exists; deleting")
	os.remove(astname + ".reg")
	
#load MN****.reg file
predictFile = astname + ".reg"
f=open(predictFile, 'a+')
f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nphysical\n')

#write the optimized strings to MN****.reg file
for i in range(0, len(optx)):
	f.write("circle(" + str(optx[i]) + "," + str(opty[i]) + "," + str(regsize) + ")\n")
	
f.close()	
'''

'''
#numRegs = len(xvals)
optRegData = np.loadtxt(astname+".reg", skiprows=3, dtype='str')
optxvals=[]
for line in optRegData:
	line=line.split(",")
	x=line[0][7:15]		#defines x position
	#print(x)
	optxvals.append(x)		#adds x position to list
numRegs = len(optxvals)
'''

numRegs=len(optx)
print("\nNumber of current regions:", numRegs)
print('\n')	
#print("xpositions:", xvals)
#print("ypositions:", yvals)

#this section is to address the indexing issue that happens because of outliers being removed during astrometry
#print(numRegs)
#print(len(xvals))
if numRegs != len(xvals):
	print("WARNING: Length of optimized values does not match length of unoptimized values;")
	outlierDiff = len(xvals)-numRegs	#defines how many outliers were removed in astrometry
else:
	outlierDiff = 0
print(outlierDiff, "outlier(s) were removed during astrometry\n")

#finds the slope between the first two positions
slopex=find_slope(times[0], times[1], optx[1], optx[2])
print("x slope:", slopex)
slopey=find_slope(times[0], times[1], opty[1], opty[2])
print("y slope:", slopey)

#finds the y-intercept
y0x=find_yintercept(times[0], optx[0], slopex)
print("\nx y-intercept:", y0x)
y0y=find_yintercept(times[0], opty[0], slopey)
print("y y-intercept:", y0y)

#linear fit from optimized positions for both x and y
if numNewRegions < 7:
        print("performing linear fit")
        p0x=[slopex, y0x]
        p0y=[slopey, y0y]

        plt.plot(optTimes[0:numRegs], optx, color='black', linestyle="", marker='o')
        poptx, pcovx = curve_fit(lin, optTimes[0:numRegs], optx, p0x)
        #poptx, pcovx = curve_fit(lin, times[0:numRegs], optx)
        plt.plot(times, lin(times, *poptx), 'r-')
        print("\nx parameters:", poptx)
        #plt.show()

        plt.plot(optTimes[0:numRegs], opty, color='black', linestyle="", marker='o')
        popty, pcovy = curve_fit(lin, optTimes[0:numRegs], opty, p0y)
        plt.plot(times, lin(times, *popty), 'r-')
        print("y parameters:", popty)
        #plt.show()
        plt.close()

        #finds next regions using the fit
        newx=[]
        newy=[]
        print(numNewRegions+numRegs, len(times))
        if numNewRegions+numRegs >= (len(times)):
                print("\nEnd of night reached; predictions have been performed on all images")
                numNewRegions=len(times)-numRegs-outlierDiff
                f=open("/project/marz746/astdata/output/endOfNight.txt", "w")
                f.close()
                print("updated number of new regions:", numNewRegions)

        for i in range(0, numNewRegions):
                newx.append(lin(times[i + numRegs + outlierDiff], *poptx))
                newy.append(lin(times[i + numRegs + outlierDiff], *popty))
                #print(newx, newy)
		
        #plots new region to see if it's good
        plt.plot(optTimes[0:numRegs], optx, color='black', linestyle="", marker='o')
        #print("reg_predict.py Check 2")
        for i in range(0, numNewRegions):
                plt.plot(times[i + numRegs + outlierDiff], newx[i], 'bo')
        #poptx, pcovx = curve_fit(lin, times[0:3], optx)
        plt.plot(times, lin(times, *poptx), 'r-')
        plt.show()
        #print("reg_predict.py Check 3")
        
        plt.plot(optTimes[0:numRegs], opty, color='black', linestyle="", marker='o')
        for i in range(0, numNewRegions):
                plt.plot(times[i + numRegs + outlierDiff], newy[i], 'bo')
        #print("reg_predict.py Check 4")
        #popty, pcovy = curve_fit(lin, times[0:3], opty)
        plt.plot(times, lin(times, *popty), 'r-')
        plt.show()
        #print("reg_predict.py Check 5")
else:
        print("\nperforming quadratic fit\n")
        p0x = [0.01*slopex, 0.01*slopex, 1.5]
        p0y = [0.01*slopey, 0.01*slopey, 1.5]
        plt.plot(optTimes[0:numRegs], optx, color='black', linestyle="", marker='o')
        poptx, pcovx = curve_fit(quad, optTimes[0:numRegs], optx, p0x)
        #poptx, pcovx = curve_fit(quad, times[0:numRegs], optx)
        plt.plot(times, quad(times, *poptx), 'r-')
        print("x parameters:", poptx)
        #plt.show()
        
        plt.plot(optTimes[0:numRegs], opty, color='black', linestyle="", marker='o')
        popty, pcovy = curve_fit(quad, optTimes[0:numRegs], opty, p0y)
        plt.plot(times, quad(times, *popty), 'r-')
        print("y parameters:", popty)
        #plt.show()
        plt.close()
        
        #finds next regions using the fit
        newx=[]
        newy=[]
        print("numNewRegions+number of current regions:", numNewRegions+numRegs)
        #print(times)
        print("number of images in a night:", len(times))
        if numNewRegions+numRegs >= (len(times)):
                print("\nEnd of night reached; predictions have been performed on all images\n")
                numNewRegions=len(times)-numRegs-outlierDiff
                f=open("/project/marz746/astdata/output/endOfNight.txt", "w")
                f.close()
                print("updated number of new regions:", numNewRegions)
        for i in range(0, numNewRegions):
                #print("i=", i)
                #print("key:", i + numRegs + outlierDiff)
                print("time of new region:", times[i + numRegs + outlierDiff])
                newx.append(quad(times[i + numRegs + outlierDiff] , *poptx))
                newy.append(quad(times[i + numRegs + outlierDiff], *popty))
                #print(newx, newy)
        
        #plots new region to see if it's good
        plt.plot(optTimes[0:numRegs], optx, color='black', linestyle="", marker='o')
        for i in range(0, numNewRegions):
                plt.plot(times[i + numRegs + outlierDiff], newx[i], 'bo')
        #poptx, pcovx = curve_fit(quad, times[0:3], optx)
        plt.plot(times, quad(times, *poptx), 'r-')
        plt.show()
        
        plt.plot(optTimes[0:numRegs], opty, color='black', linestyle="", marker='o')
        for i in range(0, numNewRegions):
                plt.plot(times[i + numRegs + outlierDiff], newy[i], 'bo')
        #popty, pcovy = curve_fit(quad, times[0:3], opty)
        plt.plot(times, quad(times, *popty), 'r-')
        plt.show()

#writes new region to the predict.reg file as input for astrometry
print("AAAAAAAAAAAAAAAAAAAA")
if print2reg == "True":
        print("printing new regions to regfile")
        f=open(regfile, 'a+')
        for i in range(0, numNewRegions):
                f.write("circle(" + str(newx[i]) + "," + str(newy[i]) + "," + str(regsize) + "\n")#")\n")

'''
#reads in data from regfile
xvals=[]
yvals=[]
for line in regdata:
	print(line)
	line=line.split(",")
	print(line)
	x=line[0][7:15]		#defines x position
	xvals.append(x)		#adds x position to list
	print("x position:", x)
	y=line[1]			#defines y position
	print("y position:", y)
	yvals.append(y)		#adds x position to list
	regsize=line[2][0:1]		#defines regsize
	
print('\n')	
print("xpositions:", xvals)
print("ypositions:", yvals)
'''

'''
#finds the average slope between the three positions
#slope1=find_slope(xvals[0], xvals[1], yvals[0], yvals[1])
#slope2=find_slope(xvals[1], xvals[2], yvals[1], yvals[2]) 
#slope=(slope1+slope2)/2
slope=find_slope(xvals[0], xvals[2], yvals[0], yvals[2])
print("slope:", slope)  

#finds the y-intercept
#y01=find_yintercept(xvals[1], yvals[1], slope1)
#y02=find_yintercept(xvals[2], yvals[2], slope2)
#y0=(y01+y02)/2
y0=find_yintercept(xvals[2], yvals[2], slope)

#average x distance between two points
sumTerms=0
avdelx=0
for i in range(0, len(xvals)-1):
	#print("i=", i)
	#print(float(xvals[i+1])-float(xvals[i]))
	sumTerms+=float(xvals[i+1])-float(xvals[i])
	avdelx=sumTerms/(i+1)/float(imsep)
#print(avdelx)

#extrapolates the next region
newxval=float(xvals[0])+float(avdelx)
print("new x value:", newxval)
print("y intercept:", y0)
newyval=lin(slope, newxval, y0)
print("new y value:", newyval)

#appends first region to regfile
f=open(predictFile, 'a+')
f.write('# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nphysical\n')
f.write("circle(" + str(xvals[0]) + "," + str(yvals[0]) + "," + str(regsize) + ")\n")
f.write("circle(" + str(newxval) + "," + str(newyval) + "," + str(regsize) + ")\n")

#iterates over 45 for no particular reason and prints the predictions to the regfile
for i in range (0, 45):
	if (i+3)%11 == 0:
		newxval+=float(avdelx)
		continue
	else:
		newxval+=float(avdelx)
		newyval=lin(slope, newxval, y0)
		f.write("circle(" + str(newxval) + "," + str(newyval) + "," + str(regsize) + ")\n")
'''		
