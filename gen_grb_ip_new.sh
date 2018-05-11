#!/bin/bash
#gen_grb_ip_new.sh
#Priyanka Shahane,
#09 Aug 2017
#Script takes bunch_clean.evt and bunch_clean_livetime.fits and genarate quad_clean.evt for grb code  
#eg. ./gen_grb_ip_new.sh AS1A04_199T01_9000002040_13792cztM0_level2_bc.evt AS1A04_199T01_9000002040_13792cztM0_level2_bc_livetime.fits /data2/czti/testarea/priyanka/test_orbit5

############################################################################################################################################

special=$3
bcevt=$1 
live_fits=$2
#tct=$3
quad_bc_ds_evt=`echo $bcevt|sed 's/bc.evt/quad_bc_ds.evt/'`
echo $quad_bc_ds_evt
quad_bc_ds_evt=$special/$(basename $quad_bc_ds_evt)


cztdatasel infile=$bcevt gtifile=$bcevt gtitype="QUAD" outfile=$quad_bc_ds_evt clobber="y" history="y"


quad_pc_evt=`echo $quad_bc_ds_evt |sed 's/quad_bc_ds.evt/quad_bc_ds_pc.evt/'`
quad_pc_evt=$special/$(basename $quad_pc_evt)
echo $quad_pc_evt

quad_pc_double_evt=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_bc_ds_pc.dblevt/'`
quad_pc_double_evt=$special/$(basename $quad_pc_double_evt)
echo $quad_pc_double_evt

quad_livetime=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_livetime.fits/'`
quad_livetime=$special/$(basename $quad_livetime)
echo $quad_livetime

quad_badpix=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_badpix.fits/'`
quad_badpix=$special/$(basename $quad_badpix)
echo $quad_badpix

cztpixclean par_infile=$quad_bc_ds_evt par_inlivetimefile=$live_fits par_outfile1=$quad_pc_evt par_outlivetimefile=$quad_livetime  par_badpixfile=$quad_badpix par_outfile2=$quad_pc_double_evt par_writedblevt="y" par_nsigma=5 par_det_tbinsize=1 par_pix_tbinsize=1 par_det_count_thresh=1000 par_pix_count_thresh=100


quad_clean_evt=`echo $quad_bc_ds_evt|sed 's/quad_bc_ds.evt/quad_clean.evt/'`
quad_clean_evt=$special/$(basename $quad_clean_evt)
cztevtclean infile=$quad_pc_evt outfile=$quad_clean_evt alphaval="0" vetorange="0" clobber="y" isdoubleEvent="n" history="y"





rm $quad_pc_evt  $quad_badpix $quad_bc_ds_evt $quad_pc_double_evt $bcevt $bunch_clean_livetime 
