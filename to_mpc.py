import pandas as pd
import numpy as np
import sys
import os
import datetime
from astropy.time import Time
from time import gmtime, strftime
import astropy
#print(astropy.__version__)
from astropy.time import TimeISO
from astropy.coordinates import Angle

# iteration: an integer value of 1 or 2. An iteration value of 1 is for un-optimized asteroid data and an iteration value of 2 is for optimized data.
if len(sys.argv)!=3:
	print("usage: %s <asteroidname> <iteration>" % (sys.argv[0]))	
	sys.exit()

class TimeYearDayTimeCustom(TimeISO):
	name='mpc_format'
	subfmts=(('data_hms', 
		'%Y-%m-%d',
		'{year:d} {mon:02d} {day:02d}'),
                #'{year:d}-{yday:02d} {hour:02d} {min:02d} {sec:02d}'),
               ('date_hm',
                '%Y-%jT%H:%M',
                '{year:d}-{yday:03d}T{hour:02d}:{min:02d}'),
               ('date',
                '%Y-%j',
                '{year:d}-{yday:03d}'))
                
class TimeYearDayTimeCustom(TimeISO):
	name='mpc_hms'
	subfmts=(('data_hms', 
		'%H:%M:%S',
		'{hour:02d}:{min:02d}:{sec:05.2f}'),
               ('date_hm',
                '%Y-%jT%H:%M',
                '{year:d}-{yday:03d}T{hour:02d}:{min:02d}'),
               ('date',
                '%Y-%j',
                '{year:d}-{yday:03d}'))
                
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

astname=sys.argv[1]
iteration=sys.argv[2]

tracking=pd.read_csv("asteroid_tracking.txt",sep='\s+',header=None)
#photast=pd.read_csv(astname+".photast",sep='\s+',header=None)
RAdec=pd.read_csv(astname+"_"+iteration+"RAdec_deg.txt",sep='\s+',header=None)
calib = pd.read_csv('/mnt/c/Users/marzt/Documents/Research/MISHAPS_F1_bothcal.txt', sep = '\s+', header = None)
outliers = pd.read_csv(astname+"_"+iteration+"outliers.txt",sep='\s+',header=None)
#linparams=pd.read_csv("/mnt/c/Users/marzt/Documents/Research/"+astname+"linparams_y.txt", sep='\s+', header=None)

not_outliers=outliers[np.logical_not(outliers.iloc[:,-1]) & np.logical_not(outliers.iloc[:,-2])]
for i,row in not_outliers.iterrows():
	#good_lines=search_string_in_file(astname+"outliers.txt", "False False")
	print(row.iloc[0], RAdec.iloc[i,0], RAdec.iloc[i,1], RAdec.iloc[i,1][0:3], RAdec.iloc[i,1][4:6], str(RAdec.iloc[i,1][7:11]))
	
#print(good_lines)


	

######## HEADER ########

#finds the number of lines of data
tfile=not_outliers.iloc[:,1]
num_of_lines=len(tfile)
n=str(num_of_lines)

#gets current UTC time
curtime=strftime("%Y-%m-%d %H:%M:%S", gmtime())

#writes header
f=open("/mnt/c/Users/marzt/Documents/Research/"+astname+"_"+iteration+"mpc.txt", "w+")
f.write('COD W84'+'\n')
f.write('CON M. T. Penny, Louisiana State University'+'\n')
f.write('CON [penny1@lsu.edu]'+'\n')
f.write('OBS M. T. Penny, T. G. Beatty'+'\n')
f.write('MEA M. L. Newman'+'\n')
f.write('TEL 4.0-m CTIO Blanco reflector'+'\n')
f.write('NET Gaia-DR2'+'\n')
f.write('BND r'+'\n')
f.write('COM MISHAPS survey'+'\n')
f.write('NUM '+n+'\n')
f.write('ACK '+astname+' '+curtime+'\n')
f.write('AC2 penny1@lsu.edu,mnewm25@lsu.edu'+'\n')



######## MPC INFO ########

#finds decimal date of observation, columns 16-32
#t=[]
#for i in range(0, num_of_lines):
#	t.append(Time(tfile[i], format="jd", scale="utc"))
#	t[i].format='iso'

#print(t)
#print(t[0].mpc_format)
#print(("%7.5f" % (Angle(t[0].mpc_hms, unit="hourangle").degree/360))[1:])
#f.write(t.mpc_format+("%7.5f" % (Angle(t.mpc_hms, unit="hourangle").degree/360))[1:]+" ")

#right ascension, columns 33-44 and declination, columns 45-56
RAfile=RAdec.iloc[:,0]
RA=[]
RA_hm=[]
RA_sec=[]
decfile=RAdec.iloc[:,1]
dec=[]
dec_dm=[]
dec_sec=[]
for i in range(0, num_of_lines):
	RA.append(RAfile[i].replace(':',' '))
	RA_hm.append(RA[i][:5])
	RA_sec.append(RA[i][6:])
	RA_sec[i]=float(RA_sec[i])
	dec.append(decfile[i].replace(':',' '))
	dec_dm.append(dec[i][:6])
	dec_sec.append(dec[i][7:])
	dec_sec[i]=float(dec_sec[i])
	
#f.write(RA_hm+' '+'{0:.2f}'.format(RA_sec)+' ')

#declination, columns 45-56
#decfile=RAdec.iloc[:,1]
#dec=[]
#dec_dm=[]
#dec_sec=[]
#for i in range (0, num_of_lines):
#	dec.append(decfile[i].replace(':',' '))
#	dec_dm.append(dec[i][:6])
#	dec_sec.append(dec[i][7:])
#	dec_sec[i]=float(dec_sec[i])

#f.write(dec_dm+' '+'{0:.1f}'.format(dec_sec)+' ')

#columns 57-65 are blank
#f.write('         ')



#columns 66-71, observed magnitude and band #
#assigns current working directory to variable
astdir=os.getcwd()
astdir = os.path.basename(astdir)

#finds line with magnitude calibration data from MISHAPS_F1_bothcal.txt
matched_line = search_string_in_file('/mnt/c/Users/marzt/Documents/Research/MISHAPS_F1_bothcal.txt', astdir)
matched_line=int(matched_line[0])

#converts flux data to magnitudes
row=matched_line-1
c=calib.iloc[row,1]
d=calib.iloc[row,3]
fluxData=not_outliers.iloc[:,3]
#magData=(c-d-2.5*np.log10(-fluxData))
magData=(c-d-2.5*np.log10(abs(fluxData)))

#print(magData[0])
#magfile=linparams.iloc[0,5]
#f.write('{0:.2f}'.format(avg_mag)+'r')

#for i in range(0, num_of_lines):
#	f.write('     '+astname+' '+' '+' '+' '+t[i].mpc_format+("%7.5f" % (Angle(t[i].mpc_hms, unit="hourangle").degree/360))[1:]+' '+RA_hm[i]+' '+'{:05.2f}'.format(RA_sec[i])+' '+dec_dm[i]+' '+'{:04.1f}'.format(dec_sec[i])+' '+'         '+'{:04.1f}'.format(magData[i])+' r'+'      '+'W84'+'\n')
	
for i,row in not_outliers.iterrows():
	t=Time(tfile[i], format="jd", scale="utc")
	t.format='iso'
	#f.write(str(RAdec.iloc[i,1][0:3])+str(RAdec.iloc[i,1][4:6]))
	f.write('     '+astname+' '+' '+' '+' '+t.mpc_format+("%7.5f" % (Angle(t.mpc_hms, unit="hourangle").degree/360))[1:]+' '+str(RAdec.iloc[i,0][0:2])+' '+str(RAdec.iloc[i,0][3:5])+' '+str(RAdec.iloc[i,0][6:11])+' '+str(RAdec.iloc[i,1][0:3])+' '+str(RAdec.iloc[i,1][4:6])+' '+str(RAdec.iloc[i,1][7:11])+' '+'         '+'{:4.1f}'.format(magData[i])+' r'+'      '+'W84'+'\n')
	
f.close()
