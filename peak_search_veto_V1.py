#Automatic code for identifying non-Poissonian peaks in veto. 
#Version 1.2
#Written by Shabnam Iyyani 
#Date: 25-march-2018
#Priyanka Shahane       Date&Time 25/03/2018 03:05 pm
#Automated Veto_GRB_search code for finding GRB by using veto counts.
#input-MKF file
#Output- Fits file
#eg. Python2.7 peak_search_veto_V1.py AS1G06_029T01_9000000964_07040czt_level2.mkf 3 
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
print ('\n \n Required Input files: \n (1)	*.mkf \n ')

#----------------------------------User inputs------------------------------#

parser	=	argparse.ArgumentParser()
parser.add_argument("inputfile_mkf",type=str,help='Enter the path of quad clean double evt file')
parser.add_argument("outfile_fits", nargs="?", type=str, help="Stem to be prepended to all output files")
#parser.add_argument("--tbin",type=float,help='Binning time for  lightcurve, default=1.0',default=1.0)
parser.add_argument("--sigma",type=int,help='Sigma clipping (3,5,7), default=3',default=3)
args	=	parser.parse_args()

#Input file name 
event_filename=args.inputfile_mkf
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
#def poisson(k, lamb):
   # np.seterr(divide='ignore', invalid='ignore')
    #return ((lamb**k/factorial(k)) * np.exp(-lamb))


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
print max_time
print min_time
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
        print "Vcnt"
        print Vcnt
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
            #entries, bin_edges, patches = plt.hist(Vcnt1, bins=len(Vcnt1), normed=True)
            #bin_middles = 0.5*(bin_edges[1:] + bin_edges[:-1])
            a, bins, patches = plt.hist(Vcnt1, 100, normed = 1)
            M=np.mean(Vcnt1)
            V_m.append(M)
            X = np.arange( 0, 500 )
            plt.plot( X, poisson.pmf(X,M), 'r-' )
            #plt.show()
            Th_Cts = int(round(M)) + S*int(round(math.sqrt(M)))
            print "th_cnts"
            print Th_Cts
            #print n
            EM_cts['chunk%02d' % n]={}
            for k1 in Vcnt1:
                if k1 >= Th_Cts:
                    EM_cts['chunk%02d' % n][d_CEd[k1]]=k1
                    print EM_cts
                else:
                     print 'no'
        Ch_E_pks=[]
        Ch_E_sig=[]
        for key in EM_cts:
            for k2 in EM_cts[key]:
                Ch_E_pks.append(k2)
                Ch_E_sig.append(EM_cts[key][k2])
                print "###################################"
                print EM_cts[key][k2]
                print Ch_E_sig
       #Nested dictionary declaration
        d_E_Ch['Q%1d' %k] ={Ch_E_pks[z]: Ch_E_sig[z] for z in range(len(Ch_E_pks))}
    #-------------------------------Matching the peaks in different quadrants----------------------------#
    print d_E_Ch 
    QD=[]
    for key, q in map(None,d_E_Ch,range(4)):
        QD.append([])
        for s in d_E_Ch[key]:
            QD[q].append(s)
    #--------------------Matching the peak time arrays of diff quadrants-------------------------------#
    A1=list(itertools.combinations(QD, r=2))  #Generating all combinations of the elements in array QD
    z1=[]
    for i in range(len(A1)):
        z1.append(np.intersect1d(A1[i][0],A1[i][1])) 
    L1 =list(set([val for sublist in z1 for val in sublist]))
    #Creating dictionary of matched peaks and its significance for diff quadrants
    QM_tb[width]={}
    QM_tb[width]={l:(d_E_Ch['Q1'].setdefault(l,0)+d_E_Ch['Q2'].setdefault(l,0)+d_E_Ch['Q3'].setdefault(l,0)+d_E_Ch['Q4'].setdefault(l,0)-np.mean(V_m))/np.mean(V_m) for l in L1}
print QM_tb 
        #--------------Extracting the peak array found for different time binning into a single array-------------#
tb=[]
Pks=[]
Sig_all=[]
for key,q in map(None,QM_tb,range(len(binwidth))):
    tb.append([])
    Pks.append([])
    Sig_all.append([])
    for h in QM_tb[key]:
        tb[q].append(key)
        Pks[q].append(h)
        Sig_all[q].append(QM_tb[key][h])
merged_Pks=[]
merged_tbs=[]
merged_sigall=[]

for g in range(len(binwidth)):
    merged_Pks+=Pks[g]
    merged_tbs+=tb[g]
    merged_sigall+=Sig_all[g]
print merged_Pks,merged_sigall
#--------------------------------#  
#QM_tb5={}
#QM_tb5={merged_Pks[z]:[merged_sigall[z],merged_sn_ratio[z]] for z in range(len(merged_Pks))}

#zipped = zip(terminal, [0] * len(terminal))
#table = collections.OrderedDict([(lhs, collections.OrderedDict(zipped)) for lhs in left])
#orderedDict = collections.OrderedDict(sorted(QM_tb5.iteritems(), key=lambda (k,v):(v,k)))
#print orderedDict
#list1=list(orderedDict.items())
#print list1
#PKs_1=[]
#val=[]
#Sigall_1=[]
#Sn_1=[]
#PKs_1 = [item[0] for item in list1]
#val = [item[1] for item in list1]
#Sigall_1=[item[0] for item in val]
#Sn_1=[item[1] for item in val]

#print Sigall_1,Sn_1




n=range(0,50)                             
hdu = fits.PrimaryHDU(n)
a1 = merged_Pks
a2 = merged_tbs
a3 = merged_sigall


col1 = fits.Column(name='Peak time(s)', format='20A', array=a1)
col2 = fits.Column(name='Binning scale', format='E', array=a2)
col3 = fits.Column(name='Significance', format='E', array=a3)

#col4 = fits.Column(name='Lightcurve data', format='E', array=a3)                             
cols = fits.ColDefs([col1, col2, col3])
tbhdu = fits.BinTableHDU.from_columns(cols)
tbhdu.name='All peaks'
tbhdu.header['Orbit filename'] =""
tbhdu.header.comments['Orbit filename']=event_filename                             
tbhdu.header['Threshold significance (sigma)']=S   
tbhdu.header.comments['Threshold significance (sigma)']='/Threshold significance for which peaks are identified'
tbhdu.header['Time binnings']= binwidth[0] 
tbhdu.header.comments['Time binnings'] = '/The time binnings for which the analysis is done'
tbhdu.header['Energy selections']="<50.0, 50.0 -100.0, >100.0"
tbhdu.header.comments['Energy selections']='/The time binnings for which the analysis is done'
tbhdu.header['SN_ratio']="Signle/noise Ratio"
tbhdu.header.comments['SN_ratio']='/for analyasis purpose'
#tbhdu.header['Vmean']=""
#tbhdu.header.comments['SN_ratio']='/for analyasis purpose'
#tbhdu.header['SN_ratio']="Signle/noise Ratio"
#tbhdu.header.comments['SN_ratio']='/for analyasis purpose'
#tbhdu.header['SN_ratio']="Signle/noise Ratio"
#tbhdu.header.comments['SN_ratio']='/for analyasis purpose'
#tbhdu.header['SN_ratio']="Signle/noise Ratio"
#tbhdu.header.comments['SN_ratio']='/for analyasis purpose'

                             
                             
                             
#c1 = L2
#c2 = maxsig_tb 
#c3 = maxSig 

#col21 = fits.Column(name='Peak time(s)', format='20A', array=c1)
#col22 = fits.Column(name='Binning scale', format='E', array=c2)
#col23 = fits.Column(name='Significance', format='E', array=c3)
#col24 = fits.Column(name='Lightcurve data', format='E', array=c4) 
#cols2 = fits.ColDefs([col21, col22, col23])
#tbhdu2 = fits.BinTableHDU.from_columns(cols2)
#tbhdu2.name='Highest Significant peaks'
#tbhdu2.header['Orbit filename'] =""
#tbhdu2.header.comments['Orbit filename']=event_filename                             
#tbhdu2.header['Threshold significance (sigma)']=S   
#tbhdu2.header.comments['Threshold significance (sigma)']='/Threshold significance for which peaks are identified'
#bhdu2.header['Time binnings']= binwidth[0] 
#tbhdu2.header.comments['Time binnings'] = '/The time binnings for which the analysis is done'
#tbhdu2.header['Energy selections']="<50.0, 50.0 -100.0, >100.0"
#tbhdu2.header.comments['Energy selections']='/The time binnings for which the analysis is done'

hdulist=fits.HDUList([hdu, tbhdu])
hdulist.writeto(output_filename)
#stop = timeit.default_timer()
