#Code for plotting lightcurve of each peak in each quadrant for veto. 
#Written by Priyanka Shahane 
#Date: 18-june-2018 Time 11:35 AM
###################################################################################################################
from astropy.table import Table
from numpy import arange
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astroML.plotting import hist
import argparse

## Input Informations:---------------------------------------------------------------------------------------------
print ('\n \n Required Input files: \n (1)	*.mkf \n ')

#----------------------------------User inputs------------------------------#

parser	=	argparse.ArgumentParser()
parser.add_argument("inputfile_mkf",type=str,help='Enter the path of quad clean double evt file')
parser.add_argument("outfile_png", nargs="?", type=str, help="Stem to be prepended to all output files")
parser.add_argument("--tmark",type=float,help='trigger time')
args	=	parser.parse_args()

tmark=args.tmark
#----------------------------------User inputs------------------------------#
event_filename = args.inputfile_quad_clean_event
outputfits =args.outfile_png
pha_list = fits.open(event_filename, memmap=True)
#---------------------------------------------------------------------------#
#Get info of the fits file
pha_list.info()
     

pha_data1 = pha_list[1].data

Veto_Q1=pha_data1['Q1_VetoCounter']
Veto_Q2=pha_data1['Q2_VetoCounter']
Veto_Q3=pha_data1['Q3_VetoCounter']
Veto_Q4=pha_data1['Q4_VetoCounter']

#Time Q1
Time1 = pha_data1['TIME'] 



#Taking time 50s before peak time
less=  float(tmark) -500.0
    #Taking time 50s after peak time
greater= float(tmark) +500.0

#Array of time 50s before and after peak time for Q1. 
op1= Time1[ (Time1 >= less) & (Time1 <= greater)]

#Array of Veto counts 50s before and after peak time for Q1. 
ct1= Veto_Q1[ (Time1 >= less) & (Time1 <= greater)]
ct2= Veto_Q2[ (Time1 >= less) & (Time1 <= greater)]
ct3= Veto_Q3[ (Time1 >= less) & (Time1 <= greater)]
ct4= Veto_Q4[ (Time1 >= less) & (Time1 <= greater)]



#plotting peak time in each quadrant
fig=plt.figure(1)
plt.suptitle('Peak Time='+str(tmark))
plt.subplot(221)
plt.axvline([tmark], color='red', linestyle='dashed', label='Marked times')
plt.xlabel('Time(s)',labelpad=20)
plt.ylabel('Counts/Bin')

plt.hist(op1,len(ct1), weights=ct1,histtype='step',color="green")
plt.title('Q0')
plt.subplots_adjust(hspace=0.5,wspace=0.5)
plt.subplot(222)
plt.axvline([tmark], color='red', linestyle='dashed', label='Marked times')
plt.hist(op1,len(ct2), weights=ct2,histtype='step',color="green")
plt.xlabel('Time(s)',labelpad=20)
plt.ylabel('Counts/Bin')

plt.title('Q1')
plt.subplots_adjust(hspace=0.5,wspace=0.5)

plt.subplot(223)
plt.axvline([tmark], color='red', linestyle='dashed', label='Marked times')
plt.hist(op1,len(ct3), weights=ct3,histtype='step',color="green")
plt.xlabel('Time(s)',labelpad=20)
plt.ylabel('Counts/Bin')
plt.title('Q2')
plt.subplots_adjust(hspace=1.0,wspace=0.5)

plt.subplot(224)
plt.axvline([tmark], color='red', linestyle='dashed', label='Marked times')
plt.hist(op1,len(ct4), weights=ct4,histtype='step',color="green")
plt.xlabel('Time(s)',labelpad=20)
plt.ylabel('Counts/Bin')
plt.title('Q3')
plt.subplots_adjust(hspace=1.0,wspace=0.5)

fig.savefig(outputfits) 
plt.show()
