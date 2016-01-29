#!/bin/sh
#   Copyright (C) 2012 University of Oxford
#
#   Part of FSL - FMRIB's Software Library
#   http://www.fmrib.ox.ac.uk/fsl
#   fsl@fmrib.ox.ac.uk
#   
#   Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
#   Imaging of the Brain), Department of Clinical Neurology, Oxford
#   University, Oxford, UK
#   
#   
#   LICENCE
#   
#   FMRIB Software Library, Release 5.0 (c) 2012, The University of
#   Oxford (the "Software")
#   
#   The Software remains the property of the University of Oxford ("the
#   University").
#   
#   The Software is distributed "AS IS" under this Licence solely for
#   non-commercial use in the hope that it will be useful, but in order
#   that the University as a charitable foundation protects its assets for
#   the benefit of its educational and research purposes, the University
#   makes clear that no condition is made or to be implied, nor is any
#   warranty given or to be implied, as to the accuracy of the Software,
#   or that it will be suitable for any particular purpose or for use
#   under any specific conditions. Furthermore, the University disclaims
#   all responsibility for the use which is made of the Software. It
#   further disclaims any liability for the outcomes arising from using
#   the Software.
#   
#   The Licensee agrees to indemnify the University and hold the
#   University harmless from and against any and all claims, damages and
#   liabilities asserted by third parties (including claims for
#   negligence) which arise directly or indirectly from the use of the
#   Software or the sale of any products based on the Software.
#   
#   No part of the Software may be reproduced, modified, transmitted or
#   transferred in any form or by any means, electronic or mechanical,
#   without the express permission of the University. The permission of
#   the University is not required if the said reproduction, modification,
#   transmission or transference is done without financial return, the
#   conditions of this Licence are imposed upon the receiver of the
#   product, and all original and amended source code is included in any
#   transmitted product. You may be held legally responsible for any
#   copyright infringement that is caused or encouraged by your failure to
#   abide by these terms and conditions.
#   
#   You are not permitted under this Licence to use this Software
#   commercially. Use for which any financial return is received shall be
#   defined as commercial use, and includes (1) integration of all or part
#   of the source code or the Software into a product for sale or license
#   by or on behalf of Licensee to third parties or (2) use of the
#   Software or any derivative of it for research with the final aim of
#   developing software products for sale or license to a third party or
#   (3) use of the Software or any derivative of it for research with the
#   final aim of developing non-software products for sale or license to a
#   third party, or (4) use of the Software to provide any service to an
#   external organisation for which payment is received. If you are
#   interested in using the Software commercially, please contact Isis
#   Innovation Limited ("Isis"), the technology transfer company of the
#   University, to negotiate a licence. Contact details are:
#   innovation@isis.ox.ac.uk quoting reference DE/9564.
export LC_ALL=C
subjdir=$1

numfib=`${FSLDIR}/bin/imglob ${subjdir}.bedpostX/diff_slices/data_slice_0000/f*samples* | wc -w | awk '{print $1}'`

if [ `${FSLDIR}/bin/imtest ${subjdir}.bedpostX/diff_slices/data_slice_0000/f0samples` -eq 1 ];then
    numfib=$(($numfib - 1))
fi

fib=1
while [ $fib -le $numfib ]
do
    ${FSLDIR}/bin/fslmerge -z ${subjdir}.bedpostX/merged_th${fib}samples `${FSLDIR}/bin/imglob ${subjdir}.bedpostX/diff_slices/data_slice_*/th${fib}samples*`
    ${FSLDIR}/bin/fslmerge -z ${subjdir}.bedpostX/merged_ph${fib}samples `${FSLDIR}/bin/imglob ${subjdir}.bedpostX/diff_slices/data_slice_*/ph${fib}samples*`
    ${FSLDIR}/bin/fslmerge -z ${subjdir}.bedpostX/merged_f${fib}samples  `${FSLDIR}/bin/imglob ${subjdir}.bedpostX/diff_slices/data_slice_*/f${fib}samples*`
    ${FSLDIR}/bin/fslmaths ${subjdir}.bedpostX/merged_th${fib}samples -Tmean ${subjdir}.bedpostX/mean_th${fib}samples
    ${FSLDIR}/bin/fslmaths ${subjdir}.bedpostX/merged_ph${fib}samples -Tmean ${subjdir}.bedpostX/mean_ph${fib}samples
    ${FSLDIR}/bin/fslmaths ${subjdir}.bedpostX/merged_f${fib}samples -Tmean ${subjdir}.bedpostX/mean_f${fib}samples

    ${FSLDIR}/bin/make_dyadic_vectors ${subjdir}.bedpostX/merged_th${fib}samples ${subjdir}.bedpostX/merged_ph${fib}samples ${subjdir}.bedpostX/nodif_brain_mask ${subjdir}.bedpostX/dyads${fib}
    if [ $fib -ge 2 ];then
	${FSLDIR}/bin/maskdyads ${subjdir}.bedpostX/dyads${fib} ${subjdir}.bedpostX/mean_f${fib}samples
	${FSLDIR}/bin/fslmaths ${subjdir}.bedpostX/mean_f${fib}samples -div ${subjdir}.bedpostX/mean_f1samples ${subjdir}.bedpostX/mean_f${fib}_f1samples
	${FSLDIR}/bin/fslmaths ${subjdir}.bedpostX/dyads${fib}_thr0.05 -mul ${subjdir}.bedpostX/mean_f${fib}_f1samples ${subjdir}.bedpostX/dyads${fib}_thr0.05_modf${fib}
	${FSLDIR}/bin/imrm ${subjdir}.bedpostX/mean_f${fib}_f1samples
    fi

    fib=$(($fib + 1))

done

if [ `${FSLDIR}/bin/imtest ${subjdir}.bedpostX/mean_f1samples` -eq 1 ];then
    ${FSLDIR}/bin/fslmaths ${subjdir}.bedpostX/mean_f1samples -mul 0 ${subjdir}.bedpostX/mean_fsumsamples
    fib=1
    while [ $fib -le $numfib ]
    do
	${FSLDIR}/bin/fslmaths ${subjdir}.bedpostX/mean_fsumsamples -add ${subjdir}.bedpostX/mean_f${fib}samples ${subjdir}.bedpostX/mean_fsumsamples
	fib=$(($fib + 1))
    done	
fi


if [ `${FSLDIR}/bin/imtest ${subjdir}.bedpostX/diff_slices/data_slice_0000/mean_dsamples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${subjdir}.bedpostX/mean_dsamples  `${FSLDIR}/bin/imglob ${subjdir}.bedpostX/diff_slices/data_slice_*/mean_dsamples*`
fi

if [ `${FSLDIR}/bin/imtest ${subjdir}.bedpostX/diff_slices/data_slice_0000/mean_d_stdsamples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${subjdir}.bedpostX/mean_d_stdsamples  `${FSLDIR}/bin/imglob ${subjdir}.bedpostX/diff_slices/data_slice_*/mean_d_stdsamples*`
fi

if [ `${FSLDIR}/bin/imtest ${subjdir}.bedpostX/diff_slices/data_slice_0000/mean_Rsamples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${subjdir}.bedpostX/mean_Rsamples  `${FSLDIR}/bin/imglob ${subjdir}.bedpostX/diff_slices/data_slice_*/mean_Rsamples*`
fi

if [ `${FSLDIR}/bin/imtest ${subjdir}.bedpostX/diff_slices/data_slice_0000/mean_f0samples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${subjdir}.bedpostX/mean_f0samples  `${FSLDIR}/bin/imglob ${subjdir}.bedpostX/diff_slices/data_slice_*/mean_f0samples*`
fi

if [ `${FSLDIR}/bin/imtest ${subjdir}.bedpostX/diff_slices/data_slice_0000/mean_S0samples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${subjdir}.bedpostX/mean_S0samples  `${FSLDIR}/bin/imglob ${subjdir}.bedpostX/diff_slices/data_slice_*/mean_S0samples*`
fi

if [ `${FSLDIR}/bin/imtest ${subjdir}.bedpostX/diff_slices/data_slice_0000/mean_tausamples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${subjdir}.bedpostX/mean_tausamples  `${FSLDIR}/bin/imglob ${subjdir}.bedpostX/diff_slices/data_slice_*/mean_tausamples*`
fi


echo Removing intermediate files

if [ `${FSLDIR}/bin/imtest ${subjdir}.bedpostX/merged_th1samples` -eq 1 ];then
  if [ `${FSLDIR}/bin/imtest ${subjdir}.bedpostX/merged_ph1samples` -eq 1 ];then
    if [ `${FSLDIR}/bin/imtest ${subjdir}.bedpostX/merged_f1samples` -eq 1 ];then
      rm -rf ${subjdir}.bedpostX/diff_slices
      rm -f ${subjdir}/data_slice_*
      rm -f ${subjdir}/nodif_brain_mask_slice_*
      if [ `${FSLDIR}/bin/imtest ${subjdir}/grad_dev_slice_0000` -eq 1 ];then
	  rm -f ${subjdir}/grad_dev_slice_*
      fi	  
    fi
  fi
fi

echo Creating identity xfm

xfmdir=${subjdir}.bedpostX/xfms
echo 1 0 0 0 > ${xfmdir}/eye.mat
echo 0 1 0 0 >> ${xfmdir}/eye.mat
echo 0 0 1 0 >> ${xfmdir}/eye.mat
echo 0 0 0 1 >> ${xfmdir}/eye.mat

echo Done
