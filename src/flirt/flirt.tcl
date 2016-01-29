#

# FLIRT - FMRIB's Linear Image Registration Tool
#
# Mark Jenkinson, Stephen Smith and Matthew Webster, FMRIB Image Analysis Group
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

#{{{ setup

source [ file dirname [ info script ] ]/fslstart.tcl

catch { exec sh -c "${FSLDIR}/bin/flirt 2>&1 | grep -i 'flirt version' | awk '{ print \$3 }' -" } VERSION

#}}}
#{{{ flirt

proc flirt { w } {

    #{{{ setup main window etc

    global reg entries USER FSLDIR argc argv PWD PADY IGS VERSION gui_ext

    set reg($w,maxnstats) 20

    # ---- Set up Frames ----
    toplevel $w
    wm title $w "FLIRT - FMRIB's Linear Image Registration Tool - v${VERSION}"
    wm iconname $w "FLIRT"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
    frame $w.f

    TitleFrame $w.f.basic -relief groove 
    set lfbasic [ $w.f.basic getframe ]

	set PADY 3
	set IG  "image"
	set IGS "images"

#}}}
    #{{{ number of secondary images

TitleFrame $w.f.stats -relief groove 
set lfstats [ $w.f.stats getframe ]

set reg($w,nstats) 0
LabelSpinBox $w.f.nstats -label "Number of secondary $IGS to apply transform to " -textvariable reg($w,nstats) -range " 0 $reg($w,maxnstats) 1"  -command " $w.f.nstats.spin.e validate; flirt:updatestats $w" -modifycmd  " flirt:updatestats $w"
pack $w.f.nstats -in $lfstats -side top -anchor w -pady 3 -padx 5

#}}}
    #{{{ standard image

set entries($w,1) ${FSLDIR}/data/standard/MNI152_T1_2mm_brain

FileEntry  $w.f.ref -textvariable entries($w,1) -label "Reference image   " -title "Select" -width 50 -filedialog directory  -filetypes IMAGE

#}}}

    #{{{ DOF
  LabelFrame $w.f.dof -text "  Model/DOF (input to ref)"
  optionMenu2 $w.f.dof.menu reg($w,dof) 2Dmenu  "   2D to 2D registration" 2D "Rigid Body (3 parameter model)" 3Dmenu "   3D to 3D registration" TRANS "Translation Only (3 parameter model)" 6 "Rigid Body (6 parameter model)" 7 "Global Rescale (7 parameter model)" 9 "Traditional (9 parameter model)" 12 "Affine (12 parameter model)"
  pack $w.f.dof.menu 

  $w.f.dof.menu.menu entryconfigure 0 -state disabled -background black
  $w.f.dof.menu.menu entryconfigure 2 -state disabled -background black

    set reg($w,dof) 12

#}}}

    #{{{ input/high res image

FileEntry  $w.f.test -textvariable entries($w,2) -label "Input image          " -title "Select" -width 50 -filedialog directory  -filetypes IMAGE

#}}}
    #{{{ low res image

FileEntry  $w.f.test2 -textvariable entries($w,3) -label "Low res image     " -title "Select" -width 50 -filedialog directory  -filetypes IMAGE

# DOF 2
LabelFrame $w.f.doftwo -text "  Model/DOF (lowres to highres)"
optionMenu2 $w.f.doftwo.menu reg($w,doftwo) 2Dmenu  "   2D to 2D registration" 2D "Rigid Body (3 parameter model)" 3Dmenu "   3D to 3D registration" TRANS "Translation Only (3 parameter model)" 6 "Rigid Body (6 parameter model)" 7 "Global Rescale (7 parameter model)" 9 "Traditional (9 parameter model)" 12 "Affine (12 parameter model)"
pack $w.f.doftwo.menu 

$w.f.doftwo.menu.menu entryconfigure 0 -state disabled -background black
$w.f.doftwo.menu.menu entryconfigure 2 -state disabled -background black
set reg($w,doftwo) 12

#}}}
    #{{{ mode
set reg($w,mode) 1
LabelFrame $w.f.mode -text "Mode " 
optionMenu2 $w.f.mode.menu  reg($w,mode)  -command "flirt:updatemode $w" 1 "Input $IG -> Reference image" 2 "Low res $IG -> High res $IG -> Reference image"
pack $w.f.mode.menu 
#}}}
    #{{{ output image

FileEntry  $w.f.output -textvariable entries($w,4) -label "Output image       " -title "Select" -width 50 -filedialog directory  -filetypes IMAGE

#}}}

    pack $w.f.mode $w.f.ref $w.f.dof $w.f.test $w.f.output -in $lfbasic -side top -anchor w -pady $PADY -padx 5

    #{{{ secondary image(s)

set i 1
while { $i <= $reg($w,maxnstats) } {
    FileEntry  $w.f.second$i -textvariable entries($w,[ expr $i + 4 ]) -label "Secondary image $i" -title "Select" -width 50 -filedialog directory  -filetypes IMAGE
    incr i 1
}

#}}}
    pack $w.f.basic $w.f.stats -in $w.f -side top -anchor w -pady 0 -padx 5
    #{{{ advanced options

    # ---- Optional stuff ----

    collapsible frame $w.f.opts -title "Advanced Options"    


   NoteBook $w.nb -side top -bd 2 -tabpady {5 10} -arcradius 3
   $w.nb insert 0 search  -text "Search"
   $w.nb insert 1 cost    -text "Cost Function"
   $w.nb insert 2 interp  -text "Interpolation"
   $w.nb insert 3 weights -text "Weighting Volumes"

    #{{{ Search

    set lf [$w.nb getframe search]
    set reg($w,searchrxmin) -90
    set reg($w,searchrxmax) 90
    set reg($w,searchrymin) -90
    set reg($w,searchrymax) 90
    set reg($w,searchrzmin) -90
    set reg($w,searchrzmax) 90
    frame  $w.searchf

    label $w.searchf.angleslabel -text "Search Angles"
    LabelSpinBox  $w.searchf.rxmin -label "X-axis (degrees): min " -textvariable reg($w,searchrxmin) -range {-180.0 180 1 } 
    LabelSpinBox  $w.searchf.rxmax -label  "  max" -textvariable reg($w,searchrxmax) -range {-180.0 180 1 } 
    LabelSpinBox  $w.searchf.rymin -label "Y-axis (degrees): min " -textvariable reg($w,searchrymin) -range {-180.0 180 1 } 
    LabelSpinBox  $w.searchf.rymax -label  "  max" -textvariable reg($w,searchrymax) -range {-180.0 180 1 }   
    LabelSpinBox  $w.searchf.rzmin -label "Z-axis (degrees): min " -textvariable reg($w,searchrzmin) -range {-180.0 180 1 } 
    LabelSpinBox  $w.searchf.rzmax -label  "  max" -textvariable reg($w,searchrzmax) -range {-180.0 180 1 }  
    
    grid  $w.searchf.angleslabel -row 1 -column 1 -padx 3 -pady 3 -sticky w
    grid  $w.searchf.rxmin -row 2 -column 1 -padx 3 -pady 3
    grid  $w.searchf.rxmax -row 2 -column 2 -padx 3 -pady 3
    grid  $w.searchf.rymin -row 3 -column 1 -padx 3 -pady 3
    grid  $w.searchf.rymax -row 3 -column 2 -padx 3 -pady 3
    grid  $w.searchf.rzmin -row 4 -column 1 -padx 3 -pady 3
    grid  $w.searchf.rzmax -row 4 -column 2 -padx 3 -pady 3

   
    LabelFrame  $w.searchrange -text "Images" 
    optionMenu2 $w.searchrange.menu reg($w,search) -command "flirt:updatesearch $w $lf" 0 "Already virtually aligned (no search)" 1 "Not aligned, but same orientation" 2 "Incorrectly oriented"
    pack  $w.searchrange.menu 
#-in  [$w.searchrange getframe ]

    set reg($w,search) 1
    set reg($w,disablesearch_yn) 0

    pack $w.searchrange $w.searchf -in $lf -side top -anchor w -padx 3 -pady 3
$w.nb raise search


#}}}
    #{{{ Cost Function

    set costlf [$w.nb getframe cost]

    radiobutton $w.corratio -text "Correlation Ratio" \
	    -variable reg($w,cost) -value corratio -anchor w -command "flirt:updatecost $w $costlf"
    radiobutton $w.mutualinfo -text "Mutual Information" \
	    -variable reg($w,cost) -value mutualinfo -anchor w -command "flirt:updatecost $w $costlf"
    radiobutton $w.nmi -text "Normalised Mutual Information" \
	    -variable reg($w,cost) -value normmi -anchor w -command "flirt:updatecost $w $costlf"
    radiobutton $w.normcorr -text "Normalised Correlation (intra-modal)" \
	    -variable reg($w,cost) -value normcorr -anchor w -command "flirt:updatecost $w $costlf"
    radiobutton $w.leastsq -text "Least Squares (intra-modal)" \
	    -variable reg($w,cost) -value leastsq -anchor w -command "flirt:updatecost $w $costlf"
    set reg($w,bins) 256
    LabelSpinBox  $w.bins -label "Number of Histogram Bins " -textvariable reg($w,bins) -range {1 5000 1 } 
    # ---- pack ----
    pack $w.corratio $w.mutualinfo $w.nmi $w.normcorr $w.leastsq $w.bins -in $costlf -side top -anchor w -padx 3
    set reg($w,cost) corratio

#}}}
    #{{{ Interpolation

    set interplf [$w.nb getframe interp]

    label $w.interpbanner -text "Final Interpolation Method (Reslice Only)"
    radiobutton $w.trilinear -text "Tri-Linear" \
	    -variable reg($w,interp) -value trilinear -anchor w -command "flirt:updateinterp $w $interplf"
    radiobutton $w.nearestneighbour -text "Nearest Neighbour" \
	    -variable reg($w,interp) -value nearestneighbour -anchor w -command "flirt:updateinterp $w $interplf"
    radiobutton $w.spline -text "Spline" \
	    -variable reg($w,interp) -value spline -anchor w -command "flirt:updateinterp $w $interplf"
    radiobutton $w.sinc -text "Sinc" \
	    -variable reg($w,interp) -value sinc -anchor w -command "flirt:updateinterp $w $interplf"
    set reg($w,sincwidth) 7
    LabelSpinBox $w.sincwidth -label " Width of Sinc Window (full width - voxels)" -textvariable reg($w,sincwidth) -range {1 5000 1 } 
    frame $w.swinopt
    label $w.swinbanner -text "Sinc Window Options"
    radiobutton $w.rectangular -text "Rectangular" \
	    -variable reg($w,sincwindow) -value rectangular -anchor w
    radiobutton $w.hanning -text "Hanning" \
	    -variable reg($w,sincwindow) -value hanning -anchor w
    radiobutton $w.blackman -text "Blackman" \
	    -variable reg($w,sincwindow) -value blackman -anchor w
    set reg($w,sincwindow) hanning
    
    # ---- pack ----
    pack $w.interpbanner $w.trilinear -in $interplf -side top -anchor w -padx 3
    pack $w.nearestneighbour $w.spline $w.sinc -in $interplf -side top -anchor w -padx 3
    set reg($w,interp) trilinear

    pack $w.swinbanner -in $w.swinopt -side top -anchor w -padx 3
    pack $w.rectangular $w.hanning $w.blackman -in $w.swinopt -side left -anchor w -padx 3

#}}}
    #{{{ Weightings

set weightlf [$w.nb getframe weights]

set entries($w,35) ""

FileEntry  $w.wgt -textvariable entries($w,35) -label "Reference weighting  " -title "Select" -width 50 -filedialog directory  -filetypes IMAGE

set entries($w,36) ""

FileEntry  $w.iwgt -textvariable entries($w,36) -label "Input weighting         " -title "Select" -width 50 -filedialog directory  -filetypes IMAGE

set entries($w,37) ""

FileEntry  $w.iwgt2 -textvariable entries($w,37) -label  "Low res weighting     " -title "Select" -width 50 -filedialog directory  -filetypes IMAGE

pack $w.wgt $w.iwgt -in $weightlf -side top -anchor w -padx 3 -pady $PADY

#}}}

    frame $w.f.advopts
    pack $w.nb -in $w.f.advopts -side top
    pack $w.f.advopts -in $w.f.opts.b -side left -padx 8 -pady 6 -expand yes -fill both
    pack $w.f.opts -in $w.f -side left -padx 5 -pady 5

#}}}
    #{{{ buttons

    # ---- Button Frame ----

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.apply     -command "flirt:apply $w keep" \
	    -text "Go" -width 5
    bind $w.apply <Return> {
	[winfo toplevel %W].apply invoke
    }

    button $w.cancel    -command "flirt:destroy $w" \
	    -text "Exit" -width 5
    bind $w.cancel <Return> {
	[winfo toplevel %W].cancel invoke
    }

    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/redirects/flirt.html" \
	    -text "Help" -width 5

    #{{{ Utils

menubutton $w.utils -text "Utils" -menu $w.utils.menu -relief raised -bd 2

menu $w.utils.menu

$w.utils.menu add command -label "Apply FLIRT transform" -command { exec sh -c "${FSLDIR}/bin/ApplyXFM$gui_ext" & }
$w.utils.menu add command -label "Concat FLIRT transforms" -command { exec sh -c "${FSLDIR}/bin/ConcatXFM$gui_ext" & }
$w.utils.menu add command -label "Invert FLIRT transform" -command { exec sh -c "${FSLDIR}/bin/InvertXFM$gui_ext" & }

#}}}

    pack $w.btns.b -side bottom -fill x
    pack $w.apply $w.cancel $w.help $w.utils -in $w.btns.b \
	    -side left -expand yes -padx 3 -pady 10 -fill y
    
    pack $w.f $w.btns -expand yes -fill both

#}}}
}

#}}}
#{{{ flirt:apply

proc flirt:apply { w dialog } {
    global reg entries

    if { $reg($w,nstats) == 0 } {
	set statslist "grot"
    } else {
	set statslist ""
    }
    set i 1
    while { $i <= $reg($w,nstats) } {
	lappend statslist $entries($w,[ expr $i + 4 ])
	incr i 1
    }
    set status [ flirt:proc $reg($w,mode) $entries($w,1) $entries($w,2) $entries($w,3) $reg($w,nstats) $statslist $entries($w,4) $reg($w,dof) $reg($w,doftwo) $reg($w,bins) $reg($w,searchrxmin) $reg($w,searchrxmax) $reg($w,searchrymin) $reg($w,searchrymax) $reg($w,searchrzmin) $reg($w,searchrzmax) $reg($w,disablesearch_yn) $reg($w,cost) $reg($w,interp) $reg($w,sincwidth) $reg($w,sincwindow) $entries($w,35) $entries($w,36) $entries($w,37) 1 ]

    update idletasks
    
    # destroy if the OK button was used AND a normal exit occurred
    if { $dialog == "destroy" && $status == 0 } {
	flirt:destroy $w
    }
}

#}}}
#{{{ flirt:destroy

proc flirt:destroy { w } {
    destroy $w
}

#}}}
#{{{ flirt:updatemode

proc flirt:updatemode { w } {
    global reg PADY IGS

    if { $reg($w,mode) == 1 } {
	$w.f.dof configure -text "  Model/DOF (input to ref)"
	$w.f.test configure -label "Input image          "
	$w.iwgt configure -label "Input weighting         "
	$w.f.nstats configure -label "Number of secondary $IGS to apply transform to "
	pack forget $w.f.test2
	pack forget $w.f.doftwo
	pack forget $w.iwgt2
    } else {
	$w.f.dof configure -text "  Model/DOF (highres to ref)"
	$w.f.test configure -label "High res image     "
	$w.iwgt configure -label "High res weighting    "
	$w.f.nstats configure -label "Number of secondary $IGS to apply combined transform to "
	pack $w.f.doftwo $w.f.test2 -in $w.f -side top -anchor w -pady $PADY -padx 5 -after $w.f.test
	pack $w.iwgt2 -in [$w.nb getframe weights] -side top -anchor w -padx 3 -pady $PADY
    }
}

#}}}
#{{{ flirt:updatestats

proc flirt:updatestats { w } {
    global reg PADY

    set i 1
    while { $i <= $reg($w,maxnstats) } {
	pack forget $w.f.second$i
	incr i 1
    }

    if { $reg($w,nstats) > 0 } {
	pack $w.f.second1 -in $w.f -side top -anchor w -pady $PADY -padx 5 -after $w.f.nstats
	$w.f.second1 configure -label "Secondary image "
    }
    set i 2
    while { $i <= $reg($w,nstats) } {
	$w.f.second1 configure -label "Secondary image 1"
	pack $w.f.second$i -in $w.f -side top -anchor w -pady $PADY -padx 5 -after $w.f.second[ expr $i - 1 ]
	incr i 1
    }
}

#}}}
#{{{ flirt:updatesearch

proc flirt:updatesearch { w lf } {
    global reg

    if { $reg($w,search) == 0 } {
	set reg($w,disablesearch_yn) 1
    } else {
	set reg($w,disablesearch_yn) 0

	if { $reg($w,search) == 1 } {
	    set reg($w,searchrxmin) -90
	    set reg($w,searchrxmax) 90
	    set reg($w,searchrymin) -90
	    set reg($w,searchrymax) 90
	    set reg($w,searchrzmin) -90
	    set reg($w,searchrzmax) 90
	} else {
	    set reg($w,searchrxmin) -180
	    set reg($w,searchrxmax) 180
	    set reg($w,searchrymin) -180
	    set reg($w,searchrymax) 180
	    set reg($w,searchrzmin) -180
	    set reg($w,searchrzmax) 180
	}
    }

    if { $reg($w,disablesearch_yn) } {
	pack forget $w.searchf
    } else {
	pack $w.searchf -in $lf -side top -anchor w -padx 3 -pady 3
    }
}

#}}}
#{{{ flirt:updatecost

proc flirt:updatecost { w costlf } {
    global reg

    if { [ string match $reg($w,cost) "normcorr" ] == 1  || 
         [ string match $reg($w,cost) "leastsq" ] == 1 } {
	pack forget $w.bins
    } else {
	pack $w.bins -in $costlf -side top -anchor w -padx 3
    }
}

#}}}
#{{{ flirt:updateinterp

proc flirt:updateinterp { w interplf } {
    global reg

    if { [ string match $reg($w,interp) "sinc" ] == 1 } {
	pack $w.swinopt -in $interplf -side top -anchor w -padx 40
	pack $w.sincwidth -in $interplf -side top -anchor w -padx 40
    } else {
	pack forget $w.swinopt
	pack forget $w.sincwidth
    }
    $w.nb compute_size
}

#}}}

#{{{ start TK window

wm withdraw .
flirt .rename
tkwait window .rename

#}}}

