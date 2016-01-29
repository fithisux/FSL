#{{{ copyright and setup 

#   GLM - setup GLM design files
#
#   Stephen Smith, Matthew Webster, FMRIB Image Analysis Group
#
#   Copyright (C) 2006 University of Oxford
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

#{{{ GLM GUI

proc glm { w } {
    global fmri FSLDIR USER feat_files VARS PWD HOME tempSpin
 
    #{{{ main window basic setup

feat5:setupdefaults

set fmri(analysis) 2
set fmri(filtering_yn) 0
set fmri(poststats_yn) 0
set fmri(reg_yn) 0
set fmri(wizard_type) 1
set fmri(r_count) 30
set fmri(a_count) 30
set fmri(infeat) 0
set fmri(level) 1

toplevel $w

wm title      $w "GLM Setup"
wm iconname   $w "GLM"
wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

#}}}
    #{{{ main window widgets


optionMenu2 $w.level fmri(level) -command "glm:updatelevel $w"  1 "Timeseries design" 2 "Higher-level / non-timeseries design"
set fmri(npts) 100

LabelSpinBox  $w.npts -textvariable fmri(npts) -label " # timepoints " -range {2 10000 1} -command "$w.npts.spin.e validate; glm:updatenpts $w" -modifycmd "glm:updatenpts $w"
LabelSpinBox  $w.tr -textvariable fmri(tr) -label " TR (s) " -range {0.001 10000 0.25} 
LabelSpinBox  $w.paradigm_hp -textvariable fmri(paradigm_hp) -label " High pass filter cutoff (s) " -range {1.0 10000 5} -width 5 

pack $w.level $w.npts -in $w -side top -anchor w -padx 3 -pady 3

#{{{ button Frame

frame $w.btns
    
button $w.btns.wizard -command "feat5:wizard $w" -text "Wizard"

button $w.btns.save -command "feat_file:setup_dialog $w a a a [namespace current] *.fsf {Save Feat setup} {feat5:write $w 1 0 0} {}" -text "Save"

button $w.btns.load -command "feat_file:setup_dialog $w a a a [namespace current] *.fsf {Load Feat setup} {glm:sourcefsf $w } {}" -text "Load"

button $w.btns.cancel -command "destroy $w" -text "Exit"

button $w.btns.help -command "FmribWebHelp file: ${FSLDIR}/doc/redirects/feat.html" -text "Help"

pack $w.btns.wizard $w.btns.save $w.btns.load $w.btns.cancel $w.btns.help -in $w.btns -side left -expand yes

pack $w.btns -in $w -side bottom -fill x -padx 10 -pady 10

#}}}

#}}}

    glm:updatelevel $w 
}

#}}}
#{{{ glm:sourcefsf

#This will correctly source a design fsf file (can't just put "source" in feat_file as need to have global variables in scope

proc glm:sourcefsf { w filename } {
    global fmri FSLDIR USER feat_files VARS PWD HOME tempSpin
    source $filename
    if { [ info exists fmri(w_model) ] } {
	if { [ winfo exists $fmri(w_model) ] } {
	    destroy $fmri(w_model)
	}
    }
    feat5:setup_model $w
}

#}}}
#{{{ glm:updatelevel

proc glm:updatelevel { w } {
    global fmri

    pack forget $w.tr $w.paradigm_hp
    if { $fmri(level)==1 } {
	$w.npts configure -label " # timepoints "
	pack $w.tr $w.paradigm_hp -in $w -after $w.npts -side top -anchor w -padx 3 -pady 3
    } else {
	set fmri(multiple) $fmri(npts)
	$w.npts configure -label " # inputs "
    }

    if { [ info exists fmri(w_model) ] } {
	if { [ winfo exists $fmri(w_model) ] } {
	    destroy $fmri(w_model)
	}
    }

    feat5:setup_model_vars_simple $w
    feat5:setup_model $w
}

#}}}
#{{{ glm:updatenpts

proc glm:updatenpts { w } {
    global fmri

    if { $fmri(level) == 2 && $fmri(multiple)!=$fmri(npts) } {
	glm:updatelevel $w 
    }
}

#}}}

#{{{ call GUI 

if { ! [ info exists INGUI ] } {
    wm withdraw .
    glm .r
    tkwait window .r
}

#}}}

