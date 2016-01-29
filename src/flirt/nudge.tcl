#{{{ copyright and setup

#   nudge - interactive (if slow) alignment correction
#
#   Stephen Smith and Matthew Webster FMRIB Image Analysis Group
#
#   Copyright (C) 2005-2007 University of Oxford
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

source [ file dirname [ info script ] ]/fslstart.tcl

set VARS(history) {}

#}}}
#{{{ nudge

proc nudge { w } {
    global FSLDIR VARS PWD nvars TMP

    toplevel $w
    wm title      $w "Nudge"
    wm iconname   $w "Nudge"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

    set TMP [ exec sh -c "${FSLDIR}/bin/tmpnam /tmp/nudge" ]

    #{{{ select images and initial transform

frame $w.inputs -relief raised -borderwidth 1

    FileEntry  $w.inputs.input -textvariable nvars(input) -label "Input image" -title "Select the input image" -width 40 -filedialog directory  -filetypes IMAGE -command "nudge_inputrange"

    FileEntry  $w.inputs.reference  -textvariable nvars(reference) -label "Reference image" -title "Select the reference image" -width 40 -filedialog directory  -filetypes IMAGE 

    FileEntry  $w.inputs.initial_xfm -textvariable nvars(initial_xfm) -label "Initial transformation (optional)" -title "Select the initial transformation matrix" -width 40 -filedialog directory  -filetypes *.mat

pack $w.inputs.input $w.inputs.reference $w.inputs.initial_xfm -in $w.inputs -padx 5 -pady 5 -side top -anchor w -expand yes

#}}}
    #{{{ nudge

frame $w.nudge -relief raised -borderwidth 1

#{{{ rot

frame $w.nudge.rot

label $w.nudge.rot.label -text "Rotation (degrees):"

set nvars(xrot) 0
set nvars(yrot) 0
set nvars(zrot) 0
set nvars(rotinc) 5

LabelSpinBox $w.nudge.rot.x -label " X " -textvariable nvars(xrot) -range "-10000.0 10000 $nvars(rotinc)" -width 5
LabelSpinBox $w.nudge.rot.y -label " Y " -textvariable nvars(yrot) -range "-10000.0 10000 $nvars(rotinc)" -width 5
LabelSpinBox $w.nudge.rot.z -label " Z " -textvariable nvars(zrot) -range "-10000.0 10000 $nvars(rotinc)" -width 5
LabelSpinBox $w.nudge.rot.inc -label "increment " -textvariable nvars(rotinc) -range "-10000.0 10000 5" -width 5 -command "$w.nudge.rot.inc.spin.e validate; nudge_update $w" -modifycmd "nudge_update $w"

pack $w.nudge.rot.label $w.nudge.rot.x $w.nudge.rot.y $w.nudge.rot.z $w.nudge.rot.inc -in $w.nudge.rot -padx 5 -side left -expand yes

#}}}
#{{{ trans

frame $w.nudge.trans

label $w.nudge.trans.label -text "Translation (mm):"

set nvars(xtrans) 0
set nvars(ytrans) 0
set nvars(ztrans) 0
set nvars(transinc) 5

LabelSpinBox $w.nudge.trans.x -label " X " -textvariable nvars(xtrans) -range "-10000.0 10000 $nvars(transinc)" -width 5
LabelSpinBox $w.nudge.trans.y -label " Y " -textvariable nvars(ytrans) -range "-10000.0 10000 $nvars(transinc)" -width 5
LabelSpinBox $w.nudge.trans.z -label " Z " -textvariable nvars(ztrans) -range "-10000.0 10000 $nvars(transinc)" -width 5
LabelSpinBox $w.nudge.trans.inc -label "increment " -textvariable nvars(transinc) -range "-10000.0 10000 5" -width 5 -command "$w.nudge.rot.inc.spin.e validate; nudge_update $w" -modifycmd "nudge_update $w"

pack $w.nudge.trans.label $w.nudge.trans.x $w.nudge.trans.y $w.nudge.trans.z $w.nudge.trans.inc -in $w.nudge.trans -padx 5 -side left -expand yes

#}}}
#{{{ scale

frame $w.nudge.scale

label $w.nudge.scale.label -text "Scaling (X):"

set nvars(xscale) 1
set nvars(yscale) 1
set nvars(zscale) 1
set nvars(scaleinc) .01

LabelSpinBox $w.nudge.scale.x -label " X " -textvariable nvars(xscale) -range "-10000.0 10000 $nvars(scaleinc)" -width 5
LabelSpinBox $w.nudge.scale.y -label " Y " -textvariable nvars(yscale) -range "-10000.0 10000 $nvars(scaleinc)" -width 5
LabelSpinBox $w.nudge.scale.z -label " Z " -textvariable nvars(zscale) -range "-10000.0 10000 $nvars(scaleinc)" -width 5
LabelSpinBox $w.nudge.scale.inc -label "increment " -textvariable nvars(scaleinc) -range "-10000.0 10000 0.01" -width 5 -command "$w.nudge.rot.inc.spin.e validate; nudge_update $w" -modifycmd "nudge_update $w"

pack $w.nudge.scale.label $w.nudge.scale.x $w.nudge.scale.y $w.nudge.scale.z $w.nudge.scale.inc -in $w.nudge.scale -padx 5 -side left -expand yes

#}}}

pack $w.nudge.rot $w.nudge.trans $w.nudge.scale -in $w.nudge -padx 5 -pady 5 -side top -anchor w

#}}}
    #{{{ view options

frame $w.fslview -relief raised -borderwidth 1

label $w.fslview.label -text "FSLView options for reference and input"

entry $w.fslview.reference -textvariable nvars(fslview_reference) -width 50

set nvars(fslview_input) "-l Red-Yellow -b 3,10"
entry $w.fslview.input -textvariable nvars(fslview_input) -width 50

pack $w.fslview.label $w.fslview.reference $w.fslview.input -in $w.fslview -padx 5 -pady 5 -side top -anchor w -expand yes

#}}}
    #{{{ bottom buttons

frame $w.btns -relief raised -borderwidth 1

button $w.btns.apply  -text "View" -command "nudge_run $w 0"

pack $w.btns.apply -in $w.btns -padx 5 -pady 5 -side left -expand yes

#}}}

    pack $w.inputs $w.nudge $w.fslview $w.btns -in $w -padx 5 -pady 5 -fill x 
}

#}}}
#{{{ nudge_update

proc nudge_update { w } {
    global FSLDIR nvars

    $w.nudge.rot.x configure   -range "-10000.0 10000 $nvars(rotinc)" 
    $w.nudge.rot.y configure   -range "-10000.0 10000 $nvars(rotinc)" 
    $w.nudge.rot.z configure   -range "-10000.0 10000 $nvars(rotinc)" 

    $w.nudge.trans.x configure -range "-10000.0 10000 $nvars(transinc)"
    $w.nudge.trans.y configure -range "-10000.0 10000 $nvars(transinc)"
    $w.nudge.trans.z configure -range "-10000.0 10000 $nvars(transinc)"

    $w.nudge.scale.x configure -range "-10000.0 10000 $nvars(scaleinc)"
    $w.nudge.scale.y configure -range "-10000.0 10000 $nvars(scaleinc)"
    $w.nudge.scale.z configure -range "-10000.0 10000 $nvars(scaleinc)"
}


#}}}
#{{{ nudge_inputrange

proc nudge_inputrange { { dummy "" } } {
    global FSLDIR nvars

    set minmax [ fsl:exec "${FSLDIR}/bin/fslstats $nvars(input) -l 0.1 -r" ]
    set nvars(fslview_input) "-l Red-Yellow -b [ lindex $minmax 0 ],[ lindex $minmax 1 ]"
}

#}}}
#{{{ nudge_run

proc nudge_run { w dummy } {
    global FSLDIR nvars fslviewpid TMP

    fsl:exec "${FSLDIR}/bin/makerot -o ${TMP}_tmp.xfm --cov=$nvars(reference) -a 1,0,0 -t $nvars(xrot)"
    if { [ file exists $nvars(initial_xfm) ] } {
	fsl:exec "${FSLDIR}/bin/convert_xfm -omat ${TMP}.xfm -concat ${TMP}_tmp.xfm $nvars(initial_xfm)"
    } else {
	fsl:exec "cp ${TMP}_tmp.xfm ${TMP}.xfm"
    }

    fsl:exec "${FSLDIR}/bin/makerot -o ${TMP}_tmp.xfm --cov=$nvars(reference) -a 0,1,0 -t $nvars(yrot)"
    fsl:exec "${FSLDIR}/bin/convert_xfm -omat ${TMP}.xfm -concat ${TMP}_tmp.xfm ${TMP}.xfm"

    fsl:exec "${FSLDIR}/bin/makerot -o ${TMP}_tmp.xfm --cov=$nvars(reference) -a 0,0,1 -t $nvars(zrot)"
    fsl:exec "${FSLDIR}/bin/convert_xfm -omat ${TMP}.xfm -concat ${TMP}_tmp.xfm ${TMP}.xfm"

    set nudgexfm [ open ${TMP}_tmp.xfm "w" ]
    puts $nudgexfm "$nvars(xscale) 0 0 $nvars(xtrans)
0 $nvars(yscale) 0 $nvars(ytrans)
0 0 $nvars(zscale) $nvars(ztrans)
0 0 0 1
"
    close $nudgexfm
    fsl:exec "${FSLDIR}/bin/convert_xfm -omat ${TMP}.xfm -concat ${TMP}_tmp.xfm ${TMP}.xfm"

    fsl:exec "${FSLDIR}/bin/flirt -ref $nvars(reference) -in $nvars(input) -out $TMP -applyxfm -init ${TMP}.xfm"

    if { [ info exists fslviewpid ] } {
	set fslviewpid [ expr $fslviewpid + 1 ]
	catch { exec sh -c "kill -9 $fslviewpid" } errmsg
    }

    set FSLVIEW /usr/local/bin/fslview
    if { [ file exists ${FSLDIR}/bin/fslview ] } {
	set FSLVIEW ${FSLDIR}/bin/fslview
    }

    set fslviewpid [ exec sh -c "$FSLVIEW $nvars(reference) $nvars(fslview_reference) $TMP $nvars(fslview_input)" & ]

    puts "final transform is in ${TMP}.xfm"
}

#}}}
#{{{ call GUI and wait

if { ! [ info exists INGUI ] } {
    wm withdraw .
    nudge .r
    tkwait window .r
}

#}}}

