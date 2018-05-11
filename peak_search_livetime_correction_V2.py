#Automatic code for identifying non-Poissonian peaks in the CZTI orbit data in different time binning life time correction Added. 
#Version 1.2
#Written by Shabnam Iyyani 
#Date: 7 Sept 2017
#Co-authors in the project:Vidushi Sharma and Priyanka Shahane
#Last Edited by Priyanka        Date&Time 23/04/18  2:34 pm
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
from astroML.plotting import hist
#%matplotlib inline
#%matplotlib notebook
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.misc import factorial
from scipy.special import gamma 
from astropy.io.fits import getheader
from astropy.io.fits import update
import itertools
from numpy import inf
import argparse

## Input Informations:---------------------------------------------------------------------------------------------
print ('\n \n Required Input files: \n (1)	*_quad_clean.evt \n (2) *_quad_livetime.fits \n ')

#----------------------------------User inputs------------------------------#

parser	=	argparse.ArgumentParser()
parser.add_argument("inputfile_quad_clean_event",type=str,help='Enter the path of quad clean evt file')
parser.add_argument("inputfile_quad_clean_livetime_event",type=str,help='Enter the path of quad clean  evt file')
parser.add_argument("outfile_fits", nargs="?", type=str, help="Stem to be prepended to all output files")
#parser.add_argument("--tbin",type=float,help='Binning time for  lightcurve, default=1.0',default=1.0)
parser.add_argument("--sigma",type=int,help='Sigma clipping (3,5,7), default=3',default=3)
args	=	parser.parse_args()

#Input file name 
event_filename=args.inputfile_quad_clean_event 
#event_filename ='/home/cztipoc/AS1A02_174T04_9000000938_06826cztM0_level2_common_clean.evt'
livetime_filename=args.inputfile_quad_clean_livetime_event
#Output filename
output_filename = args.outfile_fits
#Threshold hold significance cut of S sigma
S=args.sigma
#Binning timescales 
#binwidth=[]
binwidth =[1.0] #[0.001, 0.01, 0.1, 1.0, 10.0]
#---------------------------------------------------------------------------#

#---------------Poisson function definition------------------#
# poisson function, parameter lamb is the fit parameter
def poisson(k, lamb):
    np.seterr(divide='ignore', invalid='ignore')
    return ((lamb**k/factorial(k)) * np.exp(-lamb))


#Open the fits file 
pha_list = fits.open(event_filename, memmap=True)
pha_list_livetime = fits.open(livetime_filename, memmap=True)

#Get info of the fits file
pha_list.info()

# load the data in separate variable
pha_data1 = pha_list[1].data
pha_data2 = pha_list[2].data
pha_data3 = pha_list[3].data
pha_data4 = pha_list[4].data


pha_data_live1 = pha_list_livetime[1].data
pha_data_live2 = pha_list_livetime[2].data
pha_data_live3 = pha_list_livetime[3].data
pha_data_live4 = pha_list_livetime[4].data

#pha_list_livetime1 = pha_list_livetime[2].data
#pha_list_livetime2 = pha_list_livetime[2].data
#p=len(pha_list_livetime1)
#print p
#Time Q1
Time1 = pha_data1['TIME']   #assumes the gti, bti and livetime corrected inpute file Time array
#Time Q2
Time2 = pha_data2['TIME']
#Time Q3
Time3 = pha_data3['TIME']
#Time Q4
Time4 = pha_data4['TIME']



max_time = max([max(Time1),max(Time2),max(Time3),max(Time4)])
min_time = min([min(Time1),min(Time2),min(Time3),min(Time4)])


#--------------------------------------#
#Dictionary declaration 
QM_tb={}
QM_tb1={}
QM_tb2={}
#------------------------------Choose a given time binning--------------------#
for width in binwidth:
    #print 'bin'

    #---------------Creating a uniform binning timescale----------------#
    #Bin width
    #binwidth=1    #sec
    #Bins info 
    #bins =arange(int(min(T_E)),int(max(T_E)) + binwidth, binwidth)
    n= int((max_time - min_time)/width)
    print n
    bins=[]
    for j in range(0,n):
        #print 'YY'
        bins.append(min_time+(j*width))
    
    #---------------------------------#
    #Nested dictionary declaration for a different time binning
    QM_tb[width]={}
    #---------------------------------#
    #Dictionary declaration
    d_E_Ch ={}
    d_E_Ch1 ={}
    d_E_Ch2 ={}
    #------------------- Loop over the quadrants Q0 --1 , Q1 --2, Q2 --3, Q3 --4 --------------------#
    for k in range(1,5):
        #print 'Quadrant'
  
        
        # load the data in separate variable
       
        pha_data1 = pha_list[k].data
        pha_data_live1=pha_list_livetime[k].data
        #Time
        Time = pha_data1['TIME']   #assumes the gti, bti and livetime corrected inpute file Time array
        #Time from livetime file
        Time_l1 = pha_data_live1['TIME'] 
        #fractional Exposure From Livetime file
        fx1 = pha_data_live1['FRACEXP']
        #Energy 
        E =pha_data1['Energy']
        
        #Creating a dictionary of time and energy for the entire input file
        d_ET = {Time[n]: E[n] for n in range(len(E))}
        
        #-------------------------------------------------------------------------------#
        #Nested Dictionary declaration for different quadrant pks:cts 
        d_E_Ch['Q%1d' %k] ={}
        d_E_Ch1['Q%1d' %k] ={}
        d_E_Ch2['Q%1d' %k] ={}
        #---------------------------------Creating the Chunks------------------------# 
        #Chunk duration
        dchunk = 500* width
        #No: of chunks 
        N_chunks = int((max(Time) -min(Time))/dchunk)-1
        #print "nchunk"
        #print N_chunks 
        #Choosing the chunk interval 
        i=0
        #-------------------------------------#
        #Dictionary declaration
        EM_cts={}
        SN_cts={}
        VM_cts={}
        #---------------Starts looping over different chunks-------------------#
        for n in range(N_chunks+1):#N_chunks+1):
            print 'Chunk'
            print n
            Tstart =min(Time) + n*dchunk 
            Tstop = min(Time) + (n+1)*dchunk
            d_C={}    #chunk dictionary declaration
            
            for m in Time:
                if Tstart+i<= m <=Tstop:
                    #Create a new dictionary of Time and energy for the chunk 
                    d_C[m] =d_ET[m]
            i =0.00001   #least decimal place where the Time elements change
            
            #Sorting the times in chunk according to the energy selections
            T_E1=[]
            T_E2=[]
            T_E3=[]
            for key in d_C:
                if d_C[key] <50.0:
                    T_E1.append(key)
                elif d_C[key] >=50 and d_C[key] <=100:
                    T_E2.append(key)
                else:
                    T_E3.append(key)
            
            #Dictionary declaration for pk:cts for different chunks             
            EM_cts['chunk%02d' % n]={}
            SN_cts['chunk%02d' % n]={}
            VM_cts['chunk%02d' % n]={}
            #Dictionary declaration
            d_E ={}
            #--------------------analysing each time array for each energy selection---------#
            #Declaration of mean array
            V_m=[]
            t=[]
            for T_E,x in map(None,[T_E1,T_E2,T_E3],[1,2,3]):
                #print 'Energy'
                #count=1
                fig1 = plt.figure()
                #Histogram of event times
                Cts, Edges, patches=hist(T_E,bins, histtype='step')
                # calculate binmiddles
                Tbin_mid = 0.5*(Edges[1:] + Edges[:-1])
                #Taking time same as in livetime file time 
                Tbin_mid_new = Time_l1[1:-1]
                fx_new = fx1[1:-1]
                #Dividing counts with fractional Exposure
                Cts_new=[a/b for a, b in zip(Cts, fx_new)]
               
                d_CEd = {Cts_new[i]:Tbin_mid_new[i] for i in range(len(Cts_new))}
                #print d_CEd
                Cts1=[]
                Tbin_mid1=[]
                for z in Cts_new:
                    if z>0 and z != inf:
                        Cts1.append(z)
                        Tbin_mid1.append(d_CEd[z])
                if Cts1 ==[]:
                    continue
                d_CEd = {Cts1[i]:Tbin_mid1[i] for i in range(len(Cts1))}
               #----------------------Poisson fit to the counts hist-----------------------------#
               #Plotting
                fig2 = plt.figure()
               # the bins should be of integer width, because poisson is an integer distribution
                entries, bin_edges, patches = plt.hist(Cts1, bins=100, normed=True)
                #print entries
                # calculate binmiddles
                bin_middles = 0.5*(bin_edges[1:] + bin_edges[:-1])
                # fit with curve_fit
                M, cov_matrix = curve_fit(poisson, bin_middles, entries, p0=np.mean(Cts1),maxfev=1000)
                V_m.append(M)
                #print V_m
                #plot poisson-deviation with fitted parameter
                x_plot =np.linspace(0, 200, 1000)
                plt.plot(x_plot, poisson(x_plot, *M), 'r-', lw=2)
                #plt.xlim([0,40])
                #plt.show()
                #fig1.savefig("/home/cztipoc/plot_970_07092"+str(n)+"_"+str(x)+".png")
                #fig2.savefig("/home/cztipoc/fit_964_07028_chunk_"+str(n)+"_Q"+str(k)+"_energy"+str(x)+".png")
                #-------------------------Estimating the threshold----------------------------------#
                Th_Cts = int(round(M)) + S*int(round(math.sqrt(M)))
                print Th_Cts

                #-----------------------------------Identifying the peaks---------------------------#
                #Identify the bins where counts are greater than Th_Cts from the dictionary d_CEd-- peak time:Cts
                d_E['E%02d' % x] ={}
                for k1 in d_CEd:
                    if k1 >= Th_Cts:
                        #print 'yes'
                        d_E['E%02d' % x][d_CEd[k1]]=k1
                    else:
                        print 'no'
                        #count += 1 
                        #print count
            #---------------------------Match the energy arrays for this chunk-----------------------#
            #-----Extracting the peak Times of different energy cuts into diff array ----------------#
            #print d_E
            if d_E=={}:
                continue
            T=[[],[],[]]
            for key,c in map(None,d_E,[0,1,2]):
                    #if key!=None:
                for l in d_E[key]:
                    T[c].append(l)
            #---------------Matching the peak time arrays for diff energy cuts------------------------#                 
            Q =list(set([val for sublist in [np.intersect1d(T[0],T[1]),np.intersect1d(T[0],T[2]),np.intersect1d(T[1],T[2])] for val in sublist]))
                #print Q
            #Creating the dictionary of matched peaks and its significance for energy cuts-- 
            #matched peak times: significance (added counts in matched energies/average M of three energy cuts) for this chunk 
            #Nested dictionary declaration
            EM_cts['chunk%02d' % n]={}

            EM_cts['chunk%02d' % n]= {m:(d_E['E01'].setdefault(m, 0)+d_E['E02'].setdefault(m, 0)+d_E['E03'].setdefault(m, 0))/np.mean(V_m) for m in Q}
            #For extracting all the peak times in different chunks into a single dictionary
            #matched peak times: signle to noise ratios = (added counts in matched energies - average M of three energy cuts)/ sqrt(added counts in matched energies) 
            #Nested dictionary declaration
            SN_cts['chunk%02d' % n]={}
            SN_cts['chunk%02d' % n]= {m:(d_E['E01'].setdefault(m, 0)+d_E['E02'].setdefault(m, 0)+d_E['E03'].setdefault(m, 0)-np.mean(V_m))/math.sqrt(d_E['E01'].setdefault(m, 0)+d_E['E02'].setdefault(m, 0)+d_E['E03'].setdefault(m, 0)) for m in Q}
            VM_cts['chunk%02d' % n]={}

            VM_cts['chunk%02d' % n]= {m:sum(V_m) for m in Q}
            #print "Vm_cnts"
             #   print VM_cts
                #print VM_cts
        #For extracting all the peak times in different chunks into a single dictionary
        Ch_E_pks=[]
        Ch_E_sig=[]
        Ch_E_sn=[]
        Ch_E_vm=[]
        for key in EM_cts:
            for k2 in EM_cts[key]:
                Ch_E_pks.append(k2)
                Ch_E_sig.append(EM_cts[key][k2])
                #print Ch_E_pks 
        for key in SN_cts:
            for k2 in SN_cts[key]:
                Ch_E_sn.append(SN_cts[key][k2])

        for key in VM_cts:
            for k2 in VM_cts[key]:
                Ch_E_vm.append(VM_cts[key][k2])
        print Ch_E_vm      
        #Nested dictionary declaration
        d_E_Ch['Q%1d' %k] ={Ch_E_pks[z]:Ch_E_sig[z] for z in range(len(Ch_E_pks))}
        d_E_Ch1['Q%1d' %k] ={Ch_E_pks[z]:Ch_E_sn[z] for z in range(len(Ch_E_pks))}
        d_E_Ch2['Q%1d' %k] ={Ch_E_pks[z]:Ch_E_vm[z] for z in range(len(Ch_E_pks))}
    #-------------------------------Matching the peaks in different quadrants----------------------------#
    #print d_E_Ch
    #print d_E_Ch1
    #print d_E_Ch2
    #-------------------------------Matching the peaks in different quadrants----------------------------#
    #print d_E_Ch1
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
            #print z1
    L1 =list(set([val for sublist in z1 for val in sublist]))
            #print L1
    QM_tb[width]={}
    QM_tb1[width]={}
    QM_tb2[width]={}
    QM_tb[width]={l:d_E_Ch['Q1'].setdefault(l,0)+d_E_Ch['Q2'].setdefault(l,0)+d_E_Ch['Q3'].setdefault(l,0)+d_E_Ch['Q4'].setdefault(l,0) for l in L1}
    QM_tb1[width]={l:d_E_Ch1['Q1'].setdefault(l,0)+d_E_Ch1['Q2'].setdefault(l,0)+d_E_Ch1['Q3'].setdefault(l,0)+d_E_Ch1['Q4'].setdefault(l,0)/float(q) for l in L1}
    QM_tb2[width]={l:[d_E_Ch2['Q1'].setdefault(l,0),d_E_Ch2['Q2'].setdefault(l,0),d_E_Ch2['Q3'].setdefault(l,0),d_E_Ch2['Q4'].setdefault(l,0)] for l in L1}
            #QM_tb2[width]={l:[d_E_Ch2['Q1'],d_E_Ch2['Q2'],d_E_Ch2['Q3'],d_E_Ch2['Q4']]for l in L1}
            #print QM_tb
            #print QM_tb2
            #print d_E_Ch
            #print QM_tb1
        #print QM_tb 
#--------------Extracting the peak array found for different time binning into a single array-------------#
tb=[]
Pks=[]
Sig_all=[]
Sn_ratio=[]
Vm_mean=[]
for key,q in map(None,QM_tb,range(len(binwidth))):
    tb.append([])
    Pks.append([])
    Sig_all.append([])
    for h in QM_tb[key]:
        tb[q].append(key)
        Pks[q].append(h)
        Sig_all[q].append(QM_tb[key][h])
for key,q in map(None,QM_tb1,range(len(binwidth))):
    Sn_ratio.append([])
    for k2 in QM_tb1[key]:
        Sn_ratio[q].append(QM_tb1[key][k2])

for key,q in map(None,QM_tb2,range(len(binwidth))):
    Vm_mean.append([])
    for k2 in QM_tb2[key]:
        Vm_mean[q].append(QM_tb2[key][k2])
        #print Vm_mean
Vm_Q1 = []
Vm_Q2 = []
Vm_Q3 = []
Vm_Q4 = []
Vm_Q1 = [ y[0] for x in Vm_mean for y in x ]
Vm_Q2 = [ y[1] for x in Vm_mean for y in x ]
Vm_Q3 = [ y[2] for x in Vm_mean for y in x ]
Vm_Q4 = [ y[3] for x in Vm_mean for y in x ]

        #print Vm_Q1
        #print Vm_Q2
        #print Vm_Q3
        #print Vm_Q4
#------------------------------Writing into the output fits file----------------------------#
#--------ALL PEAKS INPUT--------#
merged_Pks=[]
merged_tbs=[]
merged_sigall=[]
merged_sn_ratio=[]
#merged_Vm_Q1=[]
#merged_Vm_Q2=[]
#merged_Vm_Q3=[]
#merged_Vm_Q4=[]
for g in range(len(binwidth)):
    merged_Pks+=Pks[g]
    merged_tbs+=tb[g]
    merged_sigall+=Sig_all[g]
    merged_sn_ratio+=Sn_ratio[g] 
    merged_Vm_Q1+=Vm_Q1[g]
            #merged_Vm_Q2+=Vm_Q2[g]
            #merged_Vm_Q3+=Vm_Q3[g]
            #merged_Vm_Q4+=Vm_Q4[g]    
#print merged_Pks,merged_sigall
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
a4 = merged_sn_ratio
a5=Vm_Q1
a6=Vm_Q2
a7=Vm_Q3
a8=Vm_Q4
col1 = fits.Column(name='Peak time(s)', format='20A', array=a1)
col2 = fits.Column(name='Binning scale', format='E', array=a2)
col3 = fits.Column(name='Significance', format='E', array=a3)
col4 = fits.Column(name='Signle/noise Ratio', format='E', array=a4)
col5 = fits.Column(name='Vm_Q1', format='E', array=a5)
col6 = fits.Column(name='Vm_Q2', format='E', array=a6)
col7 = fits.Column(name='Vm_Q3', format='E', array=a7)
col8 = fits.Column(name='Vm_Q4', format='E', array=a8)
        #col4 = fits.Column(name='Lightcurve data', format='E', array=a3)                             
cols = fits.ColDefs([col1, col2, col3, col4,col5,col6,col7,col8])
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
     

hdulist=fits.HDUList([hdu, tbhdu])
hdulist.writeto(output_filename)
        #stop = timeit.default_timer()

        #print stop - start



