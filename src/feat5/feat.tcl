#
#   GUI for FEAT - FMRI Expert Analysis Tool
#
#   Stephen Smith & Matthew Webster, FMRIB Analysis Group
#
#   Copyright (C) 1999-2008 University of Oxford
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

proc feat5 { w } {
    global fmri PXHOME FSLDIR USER feat_files unwarp_files unwarp_files_mag initial_highres_files highres_files VARS argc argv PWD gui_ext HOME tempSpin
 
    # main window

feat5:setupdefaults

set tempSpin -1

toplevel $w

wm title      $w "FEAT - FMRI Expert Analysis Tool v$fmri(version)"
wm iconname   $w "FEAT [ expr int($fmri(version)) ]"
wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

set fmri(donemodel) 0

#
    # mode

frame $w.mode

optionMenu2 $w.mode.level fmri(level) -command "feat5:updatelevel $w" 1 "First-level analysis" 2 "Higher-level analysis" 
#bind  $w.mode.level.menu <Leave>  "feat5:updatelevel $w"
#trace variable fmri(level) w "feat5:updatelevel $w" (need to put 3 dummys in proc header)
#two other ways of running "feat5:updatelevel $w" when the menu is used, each with their own drawbacks
#NB The trace is probably the WORST method, and should be replaced for other optionmenus

balloonhelp_for $w.mode.level "Use \"First-level analysis\" for analysing each session's data 
-i.e. the time-series analysis of the raw 4D FMRI data.

Use \"Higher-level analysis\" for combining first-level analyses. You
can use this hierarchically - for example at second-level to analyse
across several sessions and then at third-level to analyse across
several subjects."                       

    optionMenu2 $w.mode.analysis fmri(analysis) -command "feat5:updateanalysis $w" 7 "Full analysis" 1 "Preprocessing" 2 "Statistics" 

balloonhelp_for $w.mode.analysis "You can run a full analysis - Preprocessing and Statistics
- or either of these of these options individually.

If you select \"Statistics\" you will need to select a FEAT directory 
(or directories) instead of starting with 4D image data; the results 
already produced in those FEAT directories will then be used as appropriate.

Note that if you want to run only \"Statistic\", you must select the
FEAT directory/directories before editing the contrasts or
thresholding parameters, as these will get reset on selection of the
FEAT directory/directories."

pack $w.mode.level $w.mode.analysis -in $w.mode -side left -anchor w

#
    # notebook

NoteBook $w.nb -side top -bd 2 -tabpady {5 10} -arcradius 3 
$w.nb insert 0 misc -text "Misc"    
$w.nb insert 1 data      -text "Data"     
$w.nb insert 2 filtering -text "Pre-stats"  
$w.nb insert 3 reg       -text "Registration"
$w.nb insert 4 stats     -text "Stats"
$w.nb insert 5 poststats -text "Post-stats"

# Misc

set fmri(miscf) [ $w.nb getframe misc ]

feat5:misc_gui $w

#
# Data

set fmri(dataf) [ $w.nb getframe data ]

feat5:data_gui $w

#
# Pre-statistics processing

set fmri(filteringf) [ $w.nb getframe filtering ]

feat5:prestats_gui $w

#
# Stats

set fmri(statsf) [ $w.nb getframe stats ]

checkbutton $w.prewhiten -text "Use FILM prewhitening" -variable fmri(prewhiten_yn)
balloonhelp_for $w.prewhiten "For normal first-level time series analysis you should use
prewhitening to make the statistics valid and maximally efficient. For
other data - for example, very long TR (>30s) FMRI data, PET data or
data with very few time points (<50) - this should be turned off."

optionMenu2 $w.motionevs fmri(motionevs) 0 "Don't Add Motion Parameters" 1 "Standard Motion Parameters" 2 "Standard + Extended Motion Parameters" 
balloonhelp_for $w.motionevs "You may want to include the head motion parameters (as estimated by
MCFLIRT motion correction in the Pre-stats processing) as confound EVs
in your model. This can sometimes help remove the residual effects of
motion that are still left in the data even after motion correction.

This is not strongly recommended as there is still much to learn about
residual motion effects; simply adding such confound EVs is quite a
simplistic solution. We would recommend instead turning on MELODIC in
the FEAT Pre-stats and using ICA-based denoising as a better
alternative to removing residual motion effects (see the FEAT manual
for more information on that). However, if you do wish to include
motion parameters in your model then select this option. If you do
this, then once the motion correction has been run, the translation
and rotation parameters are added as extra confound EVs in the model.

If you select Standard Motion Parameters then only the components of the main EVs
that are orthogonal to the motion confound EVs will be used in
determining significance of the effects of interest. These are also the 6 regressors used 
in the motion parameter orthogonalisation in older versions of FEAT.

Standard + Extended Motion Parameters will generate a further 18 regressors that have
been shown to reduce motion problems in some analyses

You cannot use this option unless you are carrying out both pre-stats
and stats in the same FEAT run."

FileEntry $w.motionevsbeta -textvariable fmri(motionevsbeta) -label "Voxelwise Confound List" -title "Select text file" -width 30 -filedialog directory
FileEntry $w.scriptevsbeta -textvariable fmri(scriptevsbeta) -label "BETA OPTION: Apply external script" -title "Select text file" -width 30 -filedialog directory


frame $w.confoundevs
checkbutton $w.confoundevs.yn -text "Add additional confound EVs" -variable fmri(confoundevs) -command "feat5:updatemotionevs $w"
button $w.confoundevs.enter -text "Select confound EVs text file(s)" -command "feat5:multiple_select $w 20 \"Select confound EVs text file(s)\" "
balloonhelp_for $w.confoundevs "If you want to add other confound EVs than motion parameters, that you
have precomputed for your data, then turn this option on and then
enter the filename of a raw text file (or, if you are setting up
multiple analyses, enter one text file for each input FMRI dataset to
be processed). The file can contain as many confound EVs as you like,
each in a separate column in the text file."

pack $w.confoundevs.yn -in $w.confoundevs -side left


checkbutton $w.robust -text "Use automatic outlier de-weighting" -variable fmri(robust_yn)
balloonhelp_for $w.robust "If you turn this on then FLAME will automatically detect outlier
datapoints (for each voxel, each subject's data is considered with
respect to the other subjects regarding whether it appears to be an
outlier). Outliers are then automatically de-weighted in the
multi-subject statistics.

Outlier de-weighting is only available for the mixed effects options
as it doesn't make sense in the context of a fixed effects model. It
inceases the computation time considerably.

The estimated outlier behaviour is stored in the <b>stats</b> directory in the higher-level FEAT directory. The prob_outlier1.nii.gz file is a 4D niftii file giving the probability that each subject has outlier data on a voxelwise basis. The global_prob_outlier1.nii.gz file is a 3D niftii file that indicates the size of the outlier population expressed as the proportion of subjects that are outliers. Note there are versions of these files for each variance group in the analysis."

button $w.wizard -width 20 -text "Model setup wizard" -command "feat5:wizard $w"
balloonhelp_for $w.wizard "This lets you easily setup simple common experimental designs.

At first level, the options are regular rest-A-rest-A... or
rest-A-rest-B-rest-A-rest-B... designs (block or single-event) for
normal BOLD FMRI, or a rest-A-rest-A... design for full modelling of
perfusion FMRI data.

At second level, the options are one-group t-test, two-group-unpaired
and two-group-paired t-tests.

If you need to further adjust the resulting setup, use \"Model setup
wizard\" first, then press the \"Full model setup\" button."

button $w.model -width 20 -text "Full model setup" -command "feat5:setup_model $w"
balloonhelp_for $w.model "This allows complete control of the model-based analysis to be used."

set fmri(w_model) 0

optionMenu2 $w.mixed fmri(mixed_yn) -command "feat5:updatestats $w 0; feat5:setMinimumInputs" 3 "Fixed effects" 0 "Mixed effects: Simple OLS" 2 "Mixed effects: FLAME 1" 1 "Mixed effects: FLAME 1+2"

balloonhelp_for $w.mixed "The main choice here is between fixed effects (FE) and mixed effects (ME) higher-level modelling. FE modelling is more \"sensitive\" to activation than ME, but is restricted in the inferences that can be made from its results; because FE ignores cross-session/subject variance, reported activation is with respect to the group of sessions or subjects present, and not representative of the wider population. ME does model the session/subject variability, and it therefore allows inference to be made about the wider population from which the sessions/subjects were drawn.

The FE option implements a standard weighted fixed effects model. No random effects variances are modelled or estimated. The FE error variances are the variances (varcopes) from the previous level. Weighting is introduced by allowing these variances to be unequal (heteroscedastic). Degrees-of-freedom are calculated by summing the effective degrees-of-freedom for each input from the previous level and subtracting the number of higher-level regressors.

We now discuss the different ME options.

OLS (ordinary least squares) is a fast estimation technique which ignores all lower-level variance estimation and applies a very simple higher-level model. This is the least accurate of the ME options.

For the most accurate estimation of higher-level activation you should use FLAME (FMRIB's Local Analysis of Mixed Effects) modelling and estimation. This is a sophisticated two-stage process using Bayesian modelling and estimation (for example it allows separate modelling of the variance in different subject groups, and forces random effects variance to be non-negative).

The first stage of FLAME is significantly more accurate than OLS, and nearly as fast. The second stage of FLAME increases accuracy slightly over the first stage, but is quite a lot slower (typically 45-200 minutes). It takes all voxels which FLAME stage 1 shows to be near threshold and carries out a full MCMC-based analysis at these points, to get the most accurate estimate of activation.

We generally recommend using \"FLAME 1\", as it is MUCH faster than running both stages, and nearly as accurate. The added value from running full \"FLAME 1+2\" is most significant in a highest-level analysis when you have a small number of subjects (say <10).


If you are carrying out a mid-level analysis (e.g., cross-sessions) and will be feeding this into an even higher-level analysis (e.g., cross-subjects), then you should not use the \"FLAME 1+2\" option, as it is not possible for FLAME to know in advance of the highest-level analysis what voxels will ultimately be near threshold. With respect the question of whether to use fixed-effects or mixed-effects for such mid-level analyses, it could be argued that a mixed-effects analysis should be done at the mid-level. A mixed-effects analysis would assume that the sessions are randomly sampled from a \"population\" of sessions that that subject could produce. This includes estimation of each subject's session-to-session variance. However, it is common for only a small number of sessions to be collected for each subject, making estimation of each subject's session-to-session variance impractical. One solution to this is to assume a common session-to-session variance for all subjects, thereby providing enough data for the session-to-session variance to be estimated. However, this has a downside in that you lose information about which subjects are good (i.e. low variance) and which subjects are bad (i.e. high variance). Hence, when only a small number of sessions has been collected for each subject (say, less than 10), it is recommended that you use a fixed-effects analysis at the mid-level. This in effect treats the multiple first-level sessions (for each subject) as if they were one long session. Although this does ignore the session-session variability, it is arguable that this is not of interest anyway (this is a somewhat philosophical debate). In short, fixed-effects is favoured as it avoids practical problems associated with esimating the session-to-session variance (when there are not many sessions per subject), at the same time as maintaining information about which subjects are good and bad.


If you do decide to run \"FLAME 1+2\" and the FEAT logs indicate a large difference between the stage 1 and stage 2 estimations (or, for example, the final thresholded zstat image looks \"speckled\"), this is an indication that your data is highly non-Gaussian (e.g., has one or more strong outlier subjects, or has two clearly different groups of subjects being modelled as a single group). In such a case, stage 1 estimation is quite inaccurate (OLS even more so), hence the larger-than-normal difference between stages 1 and 2. The only really good solution is to investigate in your data what is going on - for example, to find the bad outlier."

#
# Post-Stats

set fmri(poststatsf) [ $w.nb getframe poststats ]

# edit contrasts
# pre-thresholding masking

FileEntry $fmri(poststatsf).threshmask -textvariable fmri(threshmask) -label "Pre-threshold masking" -title "Select mask" -width 30 -filedialog directory  -filetypes IMAGE

balloonhelp_for $fmri(poststatsf).threshmask "If you choose a mask for \"Pre-threshold masking\" then all stats
images will be masked by the chosen mask before thresholding. There
are two reasons why you might want to do this. The first is that you
might want to constrain your search for activation to a particular
area. The second is that in doing so, you are reducing the number of
voxels tested and therefore will make any
multiple-comparison-correction in the thresholding less stringent.

The mask image chosen does not have to be a binary mask - for example,
it can be a thresholded stats image from a previous analysis (in the
same space as the data to be analysed here); only voxels containing
zero in the \"mask\" image will get zeroed in this masking process."

#
# thresholding

TitleFrame   $w.thresh -text "Thresholding" -relief groove 
set fmri(lfthresh) [ $w.thresh getframe ]

optionMenu2 $w.thresh.menu fmri(thresh) -command "feat5:updatepoststats $w" 0 "None" 1 "Uncorrected" 2 "Voxel" 3 "Cluster"

LabelSpinBox $w.prob_thresh -label "Cluster P threshold" -textvariable fmri(prob_thresh) -range {0.0 1 0.005 }  
LabelSpinBox $w.z_thresh -label "Z threshold" -textvariable fmri(z_thresh) -range {0.0 10000 0.1 } 

pack $w.thresh.menu -in $fmri(lfthresh) -side top -padx 5 -side left
balloonhelp_for $w.thresh "After carrying out the initial statistical test, the resulting Z
statistic image is then normally thresholded to show which voxels or
clusters of voxels are activated at a particular significance level.

If \"Cluster\" thresholding is selected, a Z statistic threshold is
used to define contiguous clusters.  Then each cluster's estimated
significance level (from GRF-theory) is compared with the cluster
probability threshold. Significant clusters are then used to mask the
original Z statistic image for later production of colour blobs. This
method of thresholding is an alternative to \"Voxel\"-based
correction, and is normally more sensitive to activation. You may well
swant to increase the cluster creation \"Z threshold\" if you have high
levels of activation.

If \"Voxel\" thresholding is selected, GRF-theory-based maximum height
thresholding is carried out, with thresholding at the level set, using
one-tailed testing. This test is less overly-conservative than
Bonferroni correction.

You can also choose to simply threshold the uncorrected Z statistic
values, or apply no thresholding at all."

#
# contrast masking

button $w.conmask -text "Contrast masking" -command "feat5:setup_conmask $w"

set fmri(conmask_help) "Setup the masking of contrasts by other contrasts; after thresholding
of all contrasts has taken place you can further \"threshold\" a given
Z statistic image by masking it with non-zeroed voxels from other
contrasts.

This means that of the voxels which passed thresholding in the
contrast (or F-test) of interest, only those which also survived
thresholding in the other contrasts (or F-tests) are kept.

As a further option, the generated masks can be derived from all
positive Z statistic voxels in the mask contrasts rather than all
voxels that survived thresholding."

balloonhelp_for $w.conmask $fmri(conmask_help) 

#
# rendering

TitleFrame  $w.render -text "Rendering" -relief groove 
set fmri(lfrendering) [ $w.render getframe ]

set fmri(lfrenderingtop) [ frame $fmri(lfrendering).top ]

# Z display min and max

set tmpvalzdisplay $fmri(zdisplay)

optionMenu2 $w.zmaxmenu fmri(zdisplay) -command "feat5:updatepoststats $w" 0 "Use actual Z min/max" 1 "Use preset Z min/max"

LabelSpinBox $w.zmin -label "Min" -textvariable fmri(zmin) -range {0.0 10000 1 } 
LabelSpinBox $w.zmax -label "Max" -textvariable fmri(zmax) -range {0.0 10000 1 } 
balloonhelp_for $w.zmaxmenu "The Z statistic range selected for rendering is automatically
calculated by default, to run from red (minimum Z statistic after
thresholding) to yellow (maximum Z statistic). If more than one colour
rendered image is to be produced (i.e., when multiple constrasts are
created) then the overall range of Z values is automatically found
from all of the Z statistic images, for consistent Z statistic
colour-coding.

If multiple analyses are to be carried out, \"Use preset Z min/max\"
should be chosen, and the min/max values set by hand. Again, this
ensures consistency of Z statistic colour-coding - if several
experiments are to be reported side-by-side, colours will refer to the
same Z statistic values in each picture. When using this option, you
should choose a conservatively wide range for the min and max (e.g.,
min=1, max=15), to make sure that you do not carry out unintentional
thresholding via colour rendering."

#
# render type

set tmpvalrendertype $fmri(rendertype)

optionMenu2 $w.rendertype fmri(rendertype) 0 "Solid blobs" 1 "Transparent blobs"
balloonhelp_for $w.rendertype "With \"Solid colours\" you don't see any sign of the background images
within the colour blobs; with \"Transparent colours\" you will see
through the colour blobs to the background intensity"

#

pack $w.zmaxmenu $w.rendertype -in $fmri(lfrenderingtop) -side left -anchor n
pack $fmri(lfrenderingtop) -in $fmri(lfrendering) -anchor w

#

pack $fmri(poststatsf).threshmask -in $fmri(poststatsf) -side top -anchor w -padx 5 -pady 5
pack $w.thresh $w.render          -in $fmri(poststatsf) -side top -anchor w



set fmri(zdisplay) $tmpvalzdisplay
set fmri(rendertype) $tmpvalrendertype

# background image for group stats

frame $w.bgimage

label $w.bgimage.label -text "Background image "
optionMenu2 $w.bgimage.menu fmri(bgimage) 1 "Mean highres" 2 "First highres" 3 "Mean functional" 4 "First functional" 5 "Standard space template"

pack $w.bgimage.label $w.bgimage.menu -in $w.bgimage -side top -side left


balloonhelp_for $w.bgimage "With \"Higher-level analysis\" you can select what image will be used
as the background image for the activation colour overlays. The
default of \"Mean highres\" is probably the best for relating
activation to underlying structure. For a sharper underlying image,
(but one which is not so representative of the group of subjects), you
can instead choose to use the highres image from the first selected
subject.

You can alternatively choose to use the original lowres functional
data for the overlays, or the standard-space template image."

#


checkbutton $fmri(poststatsf).tsplot_yn -text "Create time series plots" -variable fmri(tsplot_yn)
balloonhelp_for $fmri(poststatsf).tsplot_yn "If you do not wish to create the time series plots in the web page report, turn this option off."
pack $fmri(poststatsf).tsplot_yn -in $fmri(poststatsf) -side top -anchor w

#
# Registration

set fmri(regf) [ $w.nb getframe reg ]

feat5:reg_gui $w

#

set fmri(level) 1
set fmri(analysis) 7

set tmpval $fmri(paradigm_hp)
feat5:updatelevel $w 
set fmri(paradigm_hp) $tmpval

$w.nb raise data

#
    # button Frame

frame $w.btns
    
button $w.btns.apply -command "feat5:apply $w ; $w.btns.apply configure -state disabled ; after 30000 $w.btns.apply configure -state active" -text "Go"

button $w.btns.save -command "feat_file:setup_dialog $w a a a [namespace current] *.fsf {Save Feat setup} {feat5:write $w 1 1 0} {}" -text "Save"

button $w.btns.load -command "feat_file:setup_dialog $w a a a [namespace current] *.fsf {Load Feat setup} {feat5:load $w 1} {}" -text "Load"

button $w.btns.cancel -command "destroy $w" -text "Exit"

button $w.btns.help -command "FmribWebHelp file: ${FSLDIR}/doc/redirects/feat.html" -text "Help"

# Utils

menubutton $w.btns.utils -text "Utils" -menu $w.btns.utils.menu -relief raised -bd 2

menu $w.btns.utils.menu

$w.btns.utils.menu add command -label "Make_flobs - create optimal basis set (of HRF convolution kernels)" -command { exec sh -c "${FSLDIR}/bin/Make_flobs$gui_ext" & }

$w.btns.utils.menu add command -label "Featquery - get FEAT stats from ROI mask or co-ordinates" -command { exec sh -c "${FSLDIR}/bin/Featquery$gui_ext" & }

$w.btns.utils.menu add command -label "High-res FEAT stats colour rendering" -command { exec sh -c "${FSLDIR}/bin/Renderhighres$gui_ext" & }

#

pack $w.btns.apply $w.btns.save $w.btns.load $w.btns.cancel $w.btns.help $w.btns.utils -in $w.btns -side left -expand yes

#

    pack $w.mode $w.nb -in $w -side top -anchor n -padx 10 -pady 10 
    pack $w.btns -in $w -side bottom -fill x -padx 10 -pady 10 

    # load fsf file

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

#
    # updates needed after the loading of settings

if { $fmri(perfsub_yn) } {
    pack $fmri(temp).tcmenu -in $fmri(temp) -after $fmri(temp).ps_yn -side top -side left -padx 5
}

#
}

if { ! [ info exists INGUI ] } {
    wm withdraw .
    feat5 .r
    tkwait window .r
}

