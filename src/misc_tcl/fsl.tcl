#   fsl - meta-menu for fsl tools
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 1999-2003 University of Oxford
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

#{{{ setups

source [ file dirname [ info script ] ]/fslstart.tcl

set VARS(history) {}

#}}}

proc fsl { w } {

    #{{{ vars and setup

global FSLDIR argc argv PWD HOME gui_ext

set FSLVERSION [ exec sh -c "cat ${FSLDIR}/etc/fslversion" ]

toplevel $w -menu $w.f1.fmrib.menu
wm title $w "FSL $FSLVERSION"
wm iconname $w "FSL"
wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

set graphpic [image create photo -file ${FSLDIR}/doc/images/fsl-logo.gif ]

frame $w.f1

button $w.f1.label -image $graphpic -borderwidth 0 -relief raised -borderwidth 2 -command "FmribWebHelp file: ${FSLDIR}/doc/redirects/index.html"

pack $w.f1.label -in $w.f1

#}}}
    #{{{ main menu

if { [ file exists ${FSLDIR}/tcl/loadvarian.tcl ] } {
    button $w.f1.loadvarian -text "Load Varian" \
	    -command { exec sh -c "${FSLDIR}/bin/LoadVarian$gui_ext" & }
    pack $w.f1.loadvarian -in $w.f1 -fill x -padx 1 -pady 1
}

if { [ file exists ${FSLDIR}/tcl/bet.tcl ] } {
    button $w.f1.bet -text "BET brain extraction" \
	    -command { exec sh -c "${FSLDIR}/bin/Bet$gui_ext" & }
    pack $w.f1.bet -in $w.f1 -fill x -padx 1 -pady 1
}

if { [ file exists ${FSLDIR}/tcl/susan.tcl ] } {
    button $w.f1.susan -text "SUSAN noise reduction" \
	    -command { exec sh -c "${FSLDIR}/bin/Susan$gui_ext" & }
    pack $w.f1.susan -in $w.f1 -fill x -padx 1 -pady 1
}

if { [ file exists ${FSLDIR}/tcl/fast.tcl ] } {
    button $w.f1.fast -text "FAST Segmentation" \
	    -command { exec sh -c "${FSLDIR}/bin/Fast$gui_ext" & }
    pack $w.f1.fast -in $w.f1 -fill x -padx 1 -pady 1
}

if { [ file exists ${FSLDIR}/tcl/flirt.tcl ] } {
    button $w.f1.flirt -text "FLIRT linear registration" \
	    -command { exec sh -c "${FSLDIR}/bin/Flirt$gui_ext" & }
    pack $w.f1.flirt -in $w.f1 -fill x -padx 1 -pady 1
}

if { [ file exists ${FSLDIR}/tcl/feat.tcl ] } {
    button $w.f1.feat -text "FEAT FMRI analysis" \
	    -command { exec sh -c "${FSLDIR}/bin/Feat$gui_ext" & }
    pack $w.f1.feat -in $w.f1 -fill x -padx 1 -pady 1
}

if { [ file exists ${FSLDIR}/tcl/melodic.tcl ] } {
    button $w.f1.melodic -text "MELODIC ICA" \
	    -command { exec sh -c "${FSLDIR}/bin/Melodic$gui_ext" & }
    pack $w.f1.melodic -in $w.f1 -fill x -padx 1 -pady 1
}

if { [ file exists ${FSLDIR}/tcl/fdt.tcl ] } {
    button $w.f1.fdt -text "FDT diffusion" \
	    -command { exec sh -c "${FSLDIR}/bin/Fdt$gui_ext" & }
    pack $w.f1.fdt -in $w.f1 -fill x -padx 1 -pady 1
}

if { [ file exists ${FSLDIR}/tcl/possum.tcl ] } {
    button $w.f1.possum -text "POSSUM MRI simulator" \
	    -command { exec sh -c "${FSLDIR}/bin/Possum$gui_ext" & }
    pack $w.f1.possum -in $w.f1 -fill x -padx 1 -pady 1
}

if { [ file exists ${FSLDIR}/bin/fslview ] || [ file exists ${FSLDIR}/bin/fslview.exe ] } {
    button $w.f1.fslview -text "FSLView" \
	    -command { exec sh -c "${FSLDIR}/bin/fslview" & }
    pack $w.f1.fslview -in $w.f1 -fill x -padx 1 -pady 1
}

#}}}
    #{{{ Button Frame

    frame $w.btns

    #{{{ misc menu

menubutton $w.menub -text "Misc" -menu $w.menub.menu -relief raised -bd 2

menu $w.menub.menu

if { [ file exists ${FSLDIR}/tcl/applyxfm.tcl ] } {
    $w.menub.menu add command -label "Apply FLIRT transform" -command { exec sh -c "${FSLDIR}/bin/ApplyXFM$gui_ext" & }
}

if { [ file exists ${FSLDIR}/tcl/concatxfm.tcl ] } {
    $w.menub.menu add command -label "Concat FLIRT transforms" -command { exec sh -c "${FSLDIR}/bin/ConcatXFM$gui_ext" & }
}

if { [ file exists ${FSLDIR}/tcl/invertxfm.tcl ] } {
    $w.menub.menu add command -label "Invert FLIRT transform" -command { exec sh -c "${FSLDIR}/bin/InvertXFM$gui_ext" & }
}

if { [ file exists ${FSLDIR}/tcl/glm.tcl ] } {
    $w.menub.menu add command -label "GLM Setup" -command { exec sh -c "${FSLDIR}/bin/Glm$gui_ext" & }
}

if { [ file exists ${FSLDIR}/tcl/make_flobs.tcl ] } {
    $w.menub.menu add command -label "Make FLOBS" -command { exec sh -c "${FSLDIR}/bin/Make_flobs$gui_ext" & }
}

if { [ file exists ${FSLDIR}/tcl/featquery.tcl ] } {
    $w.menub.menu add command -label "Featquery" -command { exec sh -c "${FSLDIR}/bin/Featquery$gui_ext" & }
}

if { [ file exists ${FSLDIR}/tcl/renderhighres.tcl ] } {
    $w.menub.menu add command -label "Renderhighres" -command { exec sh -c "${FSLDIR}/bin/Renderhighres$gui_ext" & }
}

$w.menub.menu add command -label "Offline Documentation" -command "FmribWebHelp file: ${FSLDIR}/doc/wiki/index.html"
#}}}
 
    button $w.cancel -command "destroy $w" -text "Exit"
 
    button $w.help -command  "FmribWebHelp file: ${FSLDIR}/doc/redirects/index.html" -text "Help"

    pack $w.menub $w.cancel $w.help -in $w.btns -side left -anchor s -expand yes -padx 2 -pady 2 -fill x
 
    pack $w.f1 $w.btns -expand yes -fill x

#}}}
}

wm withdraw .
fsl .rename
tkwait window .rename
