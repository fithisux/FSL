#
#   GUI for MELODIC
#
#   Stephen Smith & Christian Beckmann, FMRIB Analysis Group
#
#   Copyright (C) 1999-2007 University of Oxford
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

#{{{ melodic:updatelevel

proc melodic:updatelevel { w } {

    global fmri

    feat5:updatereg $w

    pack forget $fmri(poststatsf).models.subject_mat $fmri(poststatsf).models.subject_con
    if { $fmri(icaopt) > 1 } {
	pack $fmri(poststatsf).models.subject_mat $fmri(poststatsf).models.subject_con -in $fmri(poststatsf).models -after $fmri(poststatsf).models.ts_con -anchor w -side top
    }

    $w.nb compute_size
}

#}}}
#{{{ melodic:updatedim

proc melodic:updatedim { w } {

    global fmri
    set f $fmri(statsf)

    if { $fmri(dim_yn) == 1 } {
	pack forget $f.dim.n
    } else {
	pack $f.dim.n -in $f.dim -side left -padx 5 -after $f.dim.yn
    }
    $w.nb compute_size
}

#}}}
#{{{ melodic:updatethresh

proc melodic:updatethresh { w } {

    global fmri
    set f $fmri(poststatsf).thresh

    if { $fmri(thresh_yn) == 0 } {
	set fmri(mmthresh) 0
	pack forget $f.n
    } else {
	set fmri(mmthresh) 0.5
	pack $f.n -in $f -side left -padx 5 -after $f.yn
    }
}

#}}}
#{{{ melodic:apply

proc melodic:apply { w } {
    global fmri feat_files unwarp_files unwarp_files_mag initial_highres_files highres_files FSLDIR HOME

    if { [ feat5:write $w 0 1 1 $fmri(feat_filename) ] } {
	return 1
    }

    set FSFROOT [ file rootname $fmri(feat_filename) ]

    catch { exec sh -c "$FSLDIR/bin/feat $FSFROOT" & } junk

    set fmri(feat_filename) [ exec sh -c "${FSLDIR}/bin/tmpnam /tmp/feat" ].fsf

    update idletasks
}

#}}}

proc melodic { w } {
    global fmri FSLDIR USER feat_files unwarp_files unwarp_files_mag initial_highres_files highres_files VARS argc argv PWD gui_ext HOME tempSpin
 
    #{{{ main window

feat5:setupdefaults

set tempSpin -1

toplevel $w

catch { exec sh -c "${FSLDIR}/bin/melodic -V | grep '\[0-9\]'" } melodicversion
set fmri(version) [ lindex $melodicversion 2 ]

wm title      $w $melodicversion
wm iconname   $w $melodicversion
wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

set fmri(inmelodic) 1
set fmri(level) 1
set fmri(analysis) 7
set fmri(icaopt) 1
set fmri(regstandard_res) 4

#}}}
    #{{{ notebook

NoteBook $w.nb -side top -bd 2 -tabpady {5 10} -arcradius 3 
$w.nb insert 0 misc      -text "Misc"    
$w.nb insert 1 data      -text "Data"     
$w.nb insert 2 filtering -text "Pre-Stats"  
$w.nb insert 3 reg       -text "Registration"
$w.nb insert 4 stats     -text "Stats"
$w.nb insert 5 poststats -text "Post-Stats"

#{{{ Misc

set fmri(miscf) [ $w.nb getframe misc ]

feat5:misc_gui $w

#}}}
#{{{ Data

set fmri(dataf) [ $w.nb getframe data ]

feat5:data_gui $w

#}}}
#{{{ Pre-ICA

set fmri(filteringf) [ $w.nb getframe filtering ]

feat5:prestats_gui $w

#}}}
#{{{ ICA

set fmri(statsf) [ $w.nb getframe stats ]

set f $fmri(statsf)

#{{{ variance normalisation

set fmri(varnorm) 1
checkbutton $f.varnorm -variable fmri(varnorm) -text "Variance-normalise timecourses"
balloonhelp_for $f.varnorm "When switched on, Melodic will rescale each time series so 
that the estimation is more influenced by the voxel-wise 
temporal dynamics and less by a voxels' mean signal. "

#}}}
#{{{ output components

frame $f.dim

set fmri(dim_yn) 1
checkbutton $f.dim.yn -variable fmri(dim_yn) -text "Automatic dimensionality estimation" -command "melodic:updatedim $w"

pack $f.dim.yn -in $f.dim -side left

balloonhelp_for $f.dim "In order to avoid overfitting, Melodic will attempt to estimate the number of 
components from the data using Bayesian estimators for the model 
order and use PCA to reduce the data prior to the IC estimation."

set fmri(dim) 1
LabelSpinBox $f.dim.n -label "Output components" -textvariable fmri(dim) -range { 1 200000 1}

#}}}
#{{{ ICA level

optionMenu2 $f.icaopt fmri(icaopt) -command "melodic:updatelevel $w" 1 "Single-session ICA" 2 "Multi-session temporal concatenation" 3 "Multi-session Tensor-ICA"

balloonhelp_for $f.icaopt "

Here you can select the kind of ICA-based analysis you want to run.

If you want to process a single-session FMRI dataset (or several
datasets with the same setup, but independently of each other) then
select \"Single-session ICA\". This runs a standard (space X time) ICA
decomposition on each of the input datasets separately.

If you want to process several sessions or subjects simultaneously,
constraining the spatial maps to be the same across sessions/subjects,
but with no constraint on the timecourses being the same, then you
might want to choose the \"Multi-session temporal concatenation\"
option. After transforming all sessions' 4D datasets into standard
space, they are all temporally concatenated together, resulting in one
large 4D dataset. Then a standard (space X time) ICA decomposition is
run.

If you want to process several sessions or subjects simultaneously,
and you want to constrain both spatial maps and timecourses to be the
same across the session/subjects, then you should choose the
\"Multi-session tensor ICA\" option. After transforming all sessions'
4D datasets into standard space, they are combined together and tensor
(space X time X sessions/subjects) ICA is run.


For analysis of multi-session or multi-subject resting FMRI data, in
order to find resting-state networks, we generally recommend using the
temporal concatenation option. In all other cases in general we
recommend using the Tensor ICA option."

#}}}

pack $f.varnorm $f.dim $f.icaopt -in $f -anchor w -side top

#}}}
#{{{ Post-ICA

set fmri(poststatsf) [ $w.nb getframe poststats ]

set f $fmri(poststatsf)

#{{{ thresholding

frame $f.thresh

set fmri(thresh_yn) 1
checkbutton $f.thresh.yn -variable fmri(thresh_yn) -text "Threshold IC maps" -command "melodic:updatethresh $w"

balloonhelp_for $f.thresh "MELODIC uses a mixture model approach to assign significance to individual 
voxels within a spatial map. The mixture model of a single Gaussian 
distribution (for the background noise within the spatial maps) and 
2 Gamma distributions (which model the 'active' voxels contained in 
the tails of the Gaussian) is fitted to the intensity histogram of 
the Z-transformed IC maps using a restricted EM algorithm. 

From this mixture model fit, MELODIC calculates voxel-wise probabilities 
of 'activation' (as the ratio of a voxels' intensity being in the 
non-background class relative to probability of the intensity being 
background noise).
Voxels above a certain threshold level are overlayed on top of 
an example volume. The default level of 0.5 will report any voxel 
where the probability of belonging to the non-background mixtures 
exceeds the probability of the voxel belonging to the background 
noise Gaussian."

set fmri(mmthresh) 0.5
entry $f.thresh.n -textvariable fmri(mmthresh) -width 10

pack $f.thresh.yn $f.thresh.n -in $f.thresh -side left

#}}}
#{{{ background image for group stats

frame $f.bgimage

label $f.bgimage.label -text "Background image "
optionMenu2 $f.bgimage.menu fmri(bgimage) 1 "Mean highres" 3 "Mean functional" 5 "Standard space template"

pack $f.bgimage.label $f.bgimage.menu -in $f.bgimage -side top -side left

balloonhelp_for $f.bgimage "With \"Higher-level analysis\" you can select what image will be used
as the background image for the activation colour overlays. The
default of \"Mean highres\" is probably the best for relating
activation to underlying structure. For a sharper underlying image,
(but one which is not so representative of the group of subjects), you
can instead choose to use the highres image from the first selected
subject.

You can alternatively choose to use the original lowres functional
data for the overlays, or the standard-space template image."

#}}}
#{{{ output options

set fmri(ostats) 0
checkbutton $f.ostats -variable fmri(ostats) -text "Output full stats folder"

balloonhelp_for $f.ostats "
When switched on, Melodic will save the thresholded IC 
maps and the probability maps inside a folder \/stats. 
This will substantially increase the amount of space used, 
so only switch this on if you intend to use these maps."

#}}}
#{{{ model

frame $f.models

FileEntry $f.models.ts_mat      -textvariable fmri(ts_model_mat)      -label " Timeseries model                 " -title "Timeseries design.mat model" -width 35 -filedialog directory -filetypes *.mat
FileEntry $f.models.ts_con      -textvariable fmri(ts_model_con)      -label " Timeseries contrasts           " -title "Timeseries design.con contrasts" -width 35 -filedialog directory -filetypes *.con
FileEntry $f.models.subject_mat -textvariable fmri(subject_model_mat) -label " Session/subjects model       " -title "Session/subjects design.mat model" -width 35 -filedialog directory -filetypes *.mat
FileEntry $f.models.subject_con -textvariable fmri(subject_model_con) -label " Session/subjects contrasts " -title "Session/subjects design.con contrasts" -width 35 -filedialog directory -filetypes *.con

balloonhelp_for $f.models "If you select a timeseries \"design.mat\" model and \"design.con\"
contrast file (e.g. as used in a first-level FEAT model-based FMRI
analysis), these will be used by MELODIC in ordering the ICA
components, and in providing richer timeseries reporting information.

If you are doing multi-session/multi-subjects ICA, you can also
optionally select a subject \"design.mat\" model and \"design.con\"
contrast file (e.g. as used in a higher-level FEAT model-based FMRI
analysis). These will be used by MELODIC in ordering the ICA
components, and in providing richer reporting information about the
multiple sessions/subjects."

pack $f.models.ts_mat $f.models.ts_con -in $f.models -anchor w -side top

#}}}

pack $f.thresh $f.bgimage $f.ostats -in $f -anchor w -side top
pack $f.models -in $f -anchor w -side top -pady 20

#}}}
#{{{ Registration

set fmri(regf) [ $w.nb getframe reg ]

feat5:reg_gui $w

#}}}

set fmri(level) 1
set fmri(analysis) 7

set tmpval $fmri(paradigm_hp)
feat5:updatelevel $w 
set fmri(paradigm_hp) $tmpval

$w.nb raise data

#}}}
    #{{{ button Frame

frame $w.btns
    
button $w.btns.apply -command "melodic:apply $w" -text "Go"

button $w.btns.save -command "feat_file:setup_dialog $w a a a [namespace current] *.fsf {Save Feat setup} {feat5:write $w 0 1 0} {}" -text "Save"

button $w.btns.load -command "feat_file:setup_dialog $w a a a [namespace current] *.fsf {Load Feat setup} {feat5:load $w 1 } {}" -text "Load"

button $w.btns.cancel -command "destroy $w" -text "Exit"

button $w.btns.help -command "FmribWebHelp file: ${FSLDIR}/doc/redirects/melodic.html" -text "Help"

#{{{ Utils

menubutton $w.btns.utils -text "Utils" -menu $w.btns.utils.menu -relief raised -bd 2

menu $w.btns.utils.menu

$w.btns.utils.menu add command -label "High-res FEAT stats colour rendering" -command { exec sh -c "${FSLDIR}/bin/Renderhighres$gui_ext" & }

#}}}

#pack $w.btns.apply $w.btns.save $w.btns.load $w.btns.cancel $w.btns.help $w.btns.utils -in $w.btns -side left -expand yes
pack $w.btns.apply $w.btns.save $w.btns.load $w.btns.cancel $w.btns.help -in $w.btns -side left -expand yes

#}}}

    pack $w.nb -in $w -side top -anchor n -padx 10 -pady 10 
    pack $w.btns -in $w -side bottom -fill x -padx 10 -pady 10 

    #{{{ load fsf file

if { $argc > 0 } {

    set inputname [ lindex $argv 0 ]

    if { [ string first / $inputname ] != 0 && [ string first ~ $inputname ] != 0 } {
	set inputname ${PWD}/$inputname
    }

    if { [ string compare [ file extension $inputname ] .fsf ] == 0 } {

	if { [ file readable $inputname ] } {
	    puts "Loading FEAT setup file $inputname"
	    feat5:load $w 1 $inputname
	} else {
	    MxPause "setup file $inputname doesn't exist!"
	}
    }

}

#}}}
    #{{{ updates needed after the loading of settings

if { $fmri(perfsub_yn) } {
    pack $fmri(temp).tcmenu -in $fmri(temp) -after $fmri(temp).ps_yn -side top -side left -padx 5
}

#}}}

    melodic:updatelevel $w
    feat5:updateanalysis $w
    $w.nb compute_size
}

if { ! [ info exists INGUI ] } {
    wm withdraw .
    melodic .r
    tkwait window .r
}

