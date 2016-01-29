#!/bin/sh

# Mark Jenkinson, FMRIB Image Analysis Group
#
# Copyright (C) 2004 University of Oxford
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

if [ $# -lt 2 ] ; then
  echo "Usage: $0 <input phase image> <output phase image>"
  exit 1;
fi

arg1=`$FSLDIR/bin/remove_ext $1`

dtype=`$FSLDIR/bin/fslval $arg1 datatype`;
if [ $dtype -eq 32 ] ; then
    phaseimt=${arg1}_realphase
    absim=${arg1}_realabs
    $FSLDIR/bin/fslcomplex -realpolar $arg1 $absim $phaseimt 
else
    phaseimt=${arg1}_phase;
    $FSLDIR/bin/fslmaths ${arg1} $phaseimt
fi

phaseout=${arg1}_phaseout

tt=0;
tmax=`fslval $arg1 dim4`;

while [ $tt -lt $tmax ] ; do

  phaseim=${phaseimt}_tmp
  $FSLDIR/bin/fslroi $phaseimt $phaseim $tt 1

  if [ $tt -eq 0 ] ; then
    xdim=`$FSLDIR/bin/fslval $phaseim pixdim1`;
    ydim=`$FSLDIR/bin/fslval $phaseim pixdim2`;
    zdim=`$FSLDIR/bin/fslval $phaseim pixdim3`;
    
    pixdim=`echo "scale=10; $xdim / 3.14159265 " | bc`;
    piydim=`echo "scale=10; $ydim / 3.14159265 " | bc`;
    pizdim=`echo "scale=10; $zdim / 3.14159265 " | bc`;
    
    echo "dimension / PI = $pixdim , $piydim , $pizdim"
    
    tmpfile=${phaseim}_$$
    
    echo "$pixdim 0 0 0" > $tmpfile
    echo "0 $piydim 0 0" >> $tmpfile
    echo "0 0 $pizdim 0" >> $tmpfile
    echo "0 0 0 1" >> $tmpfile
    
    echo "Creating xyzramp image"
    
    $FSLDIR/bin/convertwarp -m $tmpfile -o xyzramp -r $phaseim
    
    echo "Splitting into separate ramps"
    
    $FSLDIR/bin/fslroi xyzramp xramp 0 1
    $FSLDIR/bin/fslroi xyzramp yramp 1 1
    $FSLDIR/bin/fslroi xyzramp zramp 2 1
  fi

  echo "Adding pi ramps to original phase image"
  
  $FSLDIR/bin/fslmaths $phaseim -add xramp -add yramp -add zramp $2 -odt float
  
  echo "Wrapping phase back to +/- pi range"
  
  # Implements: ph_wrapped = ph - 2*pi*round(ph/(2*pi))
  #  where round(x) is implemented as int(x + 0.5) - (x<-0.5)
  $FSLDIR/bin/fslmaths $2 -div 6.2831853 -add 0.5 ${2}_tmp -odt float 
  $FSLDIR/bin/fslmaths ${2}_tmp ${2}_tmp -odt int
  $FSLDIR/bin/fslmaths ${2} -uthr -0.5 -abs -bin -mul -1 -add ${2}_tmp ${2}_tmp -odt float
  $FSLDIR/bin/fslmaths ${2}_tmp -mul -6.2831853 -add ${2} ${2} -odt float

  if [ $tt -ge 1 ] ; then
    $FSLDIR/bin/fslmerge -t $phaseout $phaseout ${2}
  else
    $FSLDIR/bin/fslmaths ${2} ${phaseout}
  fi

  tt=`echo $tt + 1 | bc`;
  # end while loop  
done

# Recombine phase with abs if complex arg1
if [ $dtype -eq 32 ] ; then
    $FSLDIR/bin/fslmaths $phaseout $phaseim
    $FSLDIR/bin/fslcomplex -complexpolar $absim $phaseim ${2}
else
    $FSLDIR/bin/fslmaths $phaseout ${2}
fi


echo "Cleaning up files"

/bin/rm -f xyzramp.* xramp.* yramp.* zramp.* $tmpfile ${2}_tmp.* ${phaseim}.* ${phaseout}.* ${phaseimt}.*
if [ $dtype -eq 32 ] ; then
    /bin/rm -f ${absim}.*
fi

