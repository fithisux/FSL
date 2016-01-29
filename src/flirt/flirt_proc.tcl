#

# FLIRT - FMRIB's Linear Image Registration Tool
#
# Mark Jenkinson and Stephen Smith, FMRIB Image Analysis Group
#
# Copyright (C) 1999-2006 University of Oxford
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


proc flirt:proc { regmode refname testname testname2 nstats statslist output dof doftwo bins searchrxmin searchrxmax searchrymin searchrymax searchrzmin searchrzmax disablesearch_yn cost interp sincwidth sincwindow refweight inweight inweight2 popups } {

global PXHOME FSLDIR USER HOME

# setup options

set flirtoptions "-bins $bins -cost $cost"
set dofoptions ""
set doftwooptions ""

if { $dof == "2D" } {
    set dofoptions "-2D -dof 12"
} else {
  if { $dof == "TRANS" } {
      set dofoptions "-dof 6 -schedule ${FSLDIR}/etc/flirtsch/sch3Dtrans_3dof"
  } else {
    set dofoptions "-dof $dof"
  }
}

if { $doftwo == "2D" } {
    set doftwooptions "-2D -dof 12"
} else {
  if { $doftwo == "TRANS" } {
      set doftwooptions "-dof 6 -schedule ${FSLDIR}/etc/flirtsch/sch3Dtrans_3dof"
  } else {
    set doftwooptions "-dof $doftwo"
  }
}


if { $disablesearch_yn } {
    set flirtoptions "$flirtoptions -searchrx 0 0 -searchry 0 0 -searchrz 0 0"
} else {
    set flirtoptions "$flirtoptions -searchrx $searchrxmin $searchrxmax -searchry $searchrymin $searchrymax -searchrz $searchrzmin $searchrzmax"
}

set flirtweights1 ""
set flirtweights2 ""
if { $refweight != "" } {
    set flirtweights1 "$flirtweights1 -refweight $refweight"
}
if { $inweight != "" } {
    set flirtweights1 "$flirtweights1 -inweight $inweight"
    set flirtweights2 "$flirtweights2 -refweight $inweight"
}
if { $inweight2 != "" } {
    set flirtweights2 "$flirtweights2 -inweight $inweight2"
}

set flirtinterp "-interp $interp"
if { $interp == "sinc" } {
    set flirtinterp "$flirtinterp -sincwidth $sincwidth -sincwindow $sincwindow"
}

set outroot [ remove_ext $output ]

# tell fsl:exec to pass to batch system if setup, and set job durations to 10 minutes (potentially used by batch scheduler)
if { $regmode == 1 } {
    fsl:exec "${FSLDIR}/bin/flirt -in $testname -ref $refname -out $output -omat ${outroot}.mat $flirtoptions $dofoptions $flirtweights1 $flirtinterp" -t 10
} else {
    fsl:exec "${FSLDIR}/bin/flirt -in $testname -ref $refname -omat ${outroot}1.mat $flirtoptions $dofoptions $flirtweights1" -t 10
    fsl:exec "${FSLDIR}/bin/flirt -in $testname2 -ref $testname -omat ${outroot}2.mat $flirtoptions $doftwooptions $flirtweights2" -t 10
    fsl:exec "${FSLDIR}/bin/convert_xfm -concat ${outroot}1.mat -omat ${outroot}.mat ${outroot}2.mat"
    fsl:exec "${FSLDIR}/bin/flirt -in $testname2 -ref $refname -out $output -applyxfm -init ${outroot}.mat $flirtinterp"
}

for { set i 1 } { $i <= $nstats } { incr i 1 } {
    set statsname [ lindex $statslist [ expr $i - 1 ] ]
    fsl:exec "${FSLDIR}/bin/flirt -in $statsname -ref $refname -out [ file rootname ${output} ]_shadowreg_[ file tail [ file rootname $statsname ] ] -applyxfm -init ${outroot}.mat $flirtinterp"
}

puts Finished

set returnval 0

    return $returnval
}

