#Code for plotting lightcurve of each peak in each quadrant. 
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
print ('\n \n Required Input files: \n (1)	*_quad_clean.evt \n ')

#----------------------------------User inputs------------------------------#

parser	=	argparse.ArgumentParser()
parser.add_argument("inputfile_quad_clean_event",type=str,help='Enter the path of quad clean double evt file')
parser.add_argument("outfile_png", nargs="?", type=str, help="Stem to be prepended to all output files")
parser.add_argument("--tmark",type=float,help='trigger time')
args	=	parser.parse_args()

tmark=args.tmark
#----------------------------------User inputs------------------------------#
event_filename = args.inputfile_quad_clean_event
#q=event_filename[14:]
#print q[:16]
outputfits =args.outfile_png
pha_list = fits.open(event_filename, memmap=True)
#pha_list1 = fits.open(outputfits, memmap=True)
#---------------------------------------------------------------------------#
#Get info of the fits file
pha_list.info()
 

pha_data1 = pha_list[1].data
pha_data2 = pha_list[2].data
pha_data3 = pha_list[3].data
pha_data4 = pha_list[4].data


#Time Q1
Time1 = pha_data1['TIME'] 
#assumes the gti, bti and livetime corrected inpute file Time array
#Time Q2
Time2 = pha_data2['TIME']
#Time Q3
Time3 = pha_data3['TIME']
#Time Q4
Time4 = pha_data4['TIME']


  #------------------- Loop over the quadrants Peak_time --------------------#

    #Taking time 50s before peak time
less= float( tmark) -100.0
    #Taking time 50s after peak time
greater= float(tmark) +100.0
    #Array of time 50s before and after peak time for Q1. 
op1= Time1[ (Time1 >= less) & (Time1 <= greater)]
op2= Time2[ (Time2 >= less) & (Time2 <= greater)]
op3= Time3[ (Time3 >= less) & (Time3 <= greater)]
op4= Time4[ (Time4 >= less) & (Time4 <= greater)]

#-----------------------Plotting peak time in each quadrant----------------#

fig=plt.figure(1)
plt.suptitle('Peak Time='+str(tmark))
plt.subplot(221)
plt.axvline([tmark], color='red', linestyle='dashed', label='Marked times')
plt.hist(op1,bins=(int)(op1[-1]-op1[0]) ,histtype='step',color="green")
plt.xlabel('Time(s)',labelpad=20)
plt.ylabel('Counts/Bin')
plt.title('Q0')
plt.subplots_adjust(hspace=0.5,wspace=0.5)

plt.subplot(222)
plt.axvline([tmark], color='red', linestyle='dashed', label='Marked times')
plt.hist(op2,bins=(int)(op2[-1]-op2[0]),histtype='step',color="green")
plt.xlabel('Time(s)',labelpad=20)
plt.ylabel('Counts/Bin')
plt.title('Q1')
plt.subplots_adjust(hspace=0.5,wspace=0.5)

plt.subplot(223)
plt.axvline([tmark], color='red', linestyle='dashed', label='Marked times')
plt.hist(op3,bins=(int)(op3[-1]-op3[0]),histtype='step',color="green")
plt.xlabel('Time(s)',labelpad=20)
plt.ylabel('Counts/Bin')
plt.title('Q2')
plt.subplots_adjust(hspace=1.0,wspace=0.5)

plt.subplot(224)
plt.axvline([tmark], color='red', linestyle='dashed', label='Marked times')
plt.hist(op4,bins=(int)(op4[-1]-op4[0]),histtype='step',color="green")
plt.xlabel('Time(s)',labelpad=20)
plt.ylabel('Counts/Bin')
plt.title('Q3')
plt.subplots_adjust(hspace=1.0,wspace=0.5)
fig.savefig(outputfits) 
plt.show()
