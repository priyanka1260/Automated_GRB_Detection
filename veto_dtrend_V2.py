#Automatic code for identifying non-Poissonian peaks in veto. 
#Version 1.2
#Written by Shabnam Iyyani 
#Date: 25-march-2018
#Co-authors in the project:Vidushi Sharma and Priyanka Shahane
#Last Edited by Priyanka        Date&Time 25/03/2018 03:05 pm
#----------------------------------------------------------------------------------------------------------#
import numpy as np
import timeit
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import collections
from astropy.io import ascii
import matplotlib.patches as mpatches
from astropy.io import fits
from astropy.table import Table
from numpy import arange
#from astroML.plotting import hist
#%matplotlib inline
#%matplotlib notebook
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.misc import factorial
from scipy.special import gamma 
from astropy.io.fits import getheader
from astropy.io.fits import update
from scipy.stats import poisson
import itertools
import argparse

## Input Informations:---------------------------------------------------------------------------------------------
print ('\n \n Required Input files: \n (1)	*_quad_clean.evt \n ')

#----------------------------------User inputs------------------------------#

parser	=	argparse.ArgumentParser()
parser.add_argument("inputfile_quad_clean_event",type=str,help='Enter the path of quad clean double evt file')
parser.add_argument("outfile_fits", nargs="?", type=str, help="Stem to be prepended to all output files")
#parser.add_argument("--tbin",type=float,help='Binning time for  lightcurve, default=1.0',default=1.0)
parser.add_argument("--sigma",type=int,help='Sigma clipping (3,5,7), default=3',default=3)
args	=	parser.parse_args()

#Input file name 
event_filename=args.inputfile_quad_clean_event 
#event_filename ='/home/cztipoc/AS1A02_174T04_9000000938_06826cztM0_level2_common_clean.evt'
#livetime_filename='/home/cztipoc/AS1G06_083T01_9000000984_07208cztM0_level2_quad_livetime.fits'
#Output filename
output_filename = args.outfile_fits
#Threshold hold significance cut of S sigma
S=args.sigma
#Binning timescales 
binwidth = [1.0]#[0.001, 0.01, 0.1, 1.0, 10.0]
#---------------------------------------------------------------------------#

#---------------Poisson function definition------------------#
# poisson function, parameter lamb is the fit parameter
def poisson(k, lamb):
    np.seterr(divide='ignore', invalid='ignore')
    return ((lamb**k/factorial(k)) * np.exp(-lamb))
def trend(time, a, b, c):
    # For given time series and some parameters, return trend value
    return a* time**2 + b* time + c

#Open the fits file 
veto_file = fits.open(event_filename, memmap=True)
#pha_list_livetime = fits.open(livetime_filename, memmap=True)

#Get info of the fits file
veto_file.info()

# load the data in separate variable
#pha_data1 = pha_list[1].data
#pha_data2 = pha_list[2].data
#pha_data3 = pha_list[3].data
#pha_data4 = pha_list[4].data

#pha_list_livetime1 = pha_list_livetime[2].data
#pha_list_livetime2 = pha_list_livetime[2].data
#p=len(pha_list_livetime1)
#print p
hdu_data=veto_file[1].data
Veto_Q1=hdu_data['Q1_VetoCounter']
#Veto_Q2=hdu_data['Q2_VetoCounter']
#Veto_Q3=hdu_data['Q3_VetoCounter']
#Veto_Q4=hdu_data['Q4_VetoCounter']
Time =hdu_data['Time']


max_time = max([max(Time)])
min_time = min([min(Time)])
#print max_time
#print min_time
#--------------------------------------#
#Dictionary declaration 
QM_tb={}
for width in binwidth:
 #---------------------------------#
    #Nested dictionary declaration for a different time binning
    QM_tb[width]={}
 #---------------------------------##Dictionary declaration
    d_E_Ch ={}

    for k in range(1,5):
        Vcnt="'Veto_Q%d'%k'" 
        Vcnt=hdu_data['Q%d_VetoCounter' %k]
                #Vcnt=hdu_data['Q1_VetoCounter']
        #print "Vcnt"
        #print Vcnt
     #-------------------------------------------------------------------------------#
        #Nested Dictionary declaration for different quadrant pks:cts 
        d_E_Ch['Q%1d' %k] ={}
        dchunk = 500* 1
        #No: of chunks 
        N_chunks = int((max(Time) -min(Time))/dchunk)-1
 
        EM_cts={}
        #Declaration of mean array
        V_m=[]
        for n in range(N_chunks+1):#N_chunks+1):
            #print 'Chunk'
            Tstart =min(Time) + n*dchunk 
            Tstop = min(Time) + (n+1)*dchunk
            d_C=[]    #chunk dictionary declaration
            d_C=Time[ (Time >= Tstart) & (Time <= Tstop)]
            #len(d_C)
            #d_C=Time[Tstart+i < Time <=Tstop]
            #i =0.00001
           # print "Q"
            #print Tstart
            #print Tstop
           
            fig1 = plt.figure()
            Vcnt1=Vcnt[(Time >= Tstart) & (Time <= Tstop)]
            d_CEd = {Vcnt1[i]:d_C[i] for i in range(len(Vcnt1))}
            Cts1=[]
            Tbin_mid1=[]
            for a in Vcnt1:
                if a>0:
                    Cts1.append(a)
                    Tbin_mid1.append(d_CEd[a])
                        #if Cts1 ==[]:
                         #   continue 
            #print Cts1
     
                        # calculate binmiddles
                        #Tbin_mid = 0.5*(Edges[1:] + Edges[:-1])
                        #print Tbin_mid
                        #Create a dictionary for counts: bin_edges
            d_CEd = {Tbin_mid1[i]:Cts1[i] for i in range(len(Tbin_mid1))}
            #print "counts################"
            #print d_CEd
           # print len(Cts1)
            #print len(Tbin_mid1)
            #entries, bin_edges, patches = plt.hist(Vcnt1, bins=len(Vcnt1), normed=True)
            CTS = np.array(Cts1)
            Tbin_Mid= np.array(Tbin_mid1)
            #print len(CTS)
            #print len(Tbin_Mid)
            params, cov = curve_fit(trend, Tbin_Mid, Cts1)
            polyfit_counts=trend(Tbin_Mid, *params)
            #print "**********************************************************"
           # print polyfit_counts
            FCts = CTS - polyfit_counts
            d_CEd1 = {Tbin_Mid[i]:FCts[i] for i in range(len(FCts))}
            #print d_CEd1
            #print Tbin_Mid[0]
            #print FCts
            #print FCts
            d_CEd3 = {key:d_CEd[key] for key in d_CEd1 if key in d_CEd}
            #print d_CEd3
            Ct=[]
            Tbin_m=[]
            for key in d_CEd3:
                    Ct.append(d_CEd3[key])
                    Tbin_m.append(key)
            
            z=np.polyfit(Tbin_m, Ct, 2)
            #print z
            dof=len(Tbin_m)-len(z)
            #print dof
            #chi_squared = np.sum((np.polyval(z, Tbin_m) - Ct) ** 2)/np.polyval(z, Tbin_m)
            p=np.polyval(z, Tbin_m)
            val=[x1 - x2 for (x1, x2) in zip(p, Ct)]
            m_val=[b1 ** 2 for b1 in val]
            #p1=np.sum(p)
            p2=np.sum(m_val)
            op=[p2 / a2 for a2 in p]
            #p=chi/dof
           # print op
            p3=np.sum(op)/dof
            print p3
            #x_plot =np.linspace(294.815, 500, 1000)
            #plt.plot(x_plot, poisson(x_plot, *z), 'r-', lw=2)
            #plt.show()
            #op=[a1 / a2 for (a1, a2) in zip(m_val, p)]
            #chi_squared=np.sum(op)
            #print op
            
           # print chi_squared
            #r=chi_squared/dof
            #print "reduced degrre 0"
            #print r 
            
