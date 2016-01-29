#   FSL interface for FDT (BEDPOSTX and PROBTRACKX)
#
#   Timothy Behrens, Heidi Johansen-Berg, Dave Flitney, Matthew Webster, Saad Jbabdi and Stam Sotiropoulos, FMRIB Image Analysis Group
#
#   Copyright (C) 2007 University of Oxford
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

#TO DO replace -filetypes * with -filetypes { } for directory selectors
source [ file dirname [ info script ] ]/fslstart.tcl
option add *FileEntry*Entry*width 35
set TCLPATH [file dirname [ info script ] ]
regsub tcl $TCLPATH bin BINPATH
regsub tcl $TCLPATH doc/redirects HTMLPATH

set VERSION "3.0"

proc mm_to_voxels { X Y Z mask } {

    global FSLDIR

    upvar $X cX
    upvar $Y cY
    upvar $Z cZ

    set vcX [ exec sh -c "echo $cX $cY $cZ | $FSLDIR/bin/std2imgcoord -img $mask -std $mask -vox - 2>/dev/null | awk '{print \$1}'" ]    
    set vcY [ exec sh -c "echo $cX $cY $cZ | $FSLDIR/bin/std2imgcoord -img $mask -std $mask -vox - 2>/dev/null | awk '{print \$2}'" ] 
    set vcZ [ exec sh -c "echo $cX $cY $cZ | $FSLDIR/bin/std2imgcoord -img $mask -std $mask -vox - 2>/dev/null | awk '{print \$3}'" ] 
    set cX $vcX
    set cY $vcY
    set cZ $vcZ
}

proc fdt:dialog { w tclstartupfile } {

    global eddy bedpost registration dtifit probtrack HTMLPATH FSLDIR VERSION INMEDX VARS
    set probtrack(tool) "probtrackx"

    if [winfo exists $w] {
        wm deiconify $w
        raise $w
        return
    }

    toplevel $w
    wm title $w "FDT - FMRIB's Diffusion Toolbox $VERSION"
    wm iconname $w "FDT"
    wm iconbitmap $w @$FSLDIR/tcl/fmrib.xbm

    #-------- Stage and Mode Options -------- 

    frame $w.tool
    optionMenu2 $w.tool.menu probtrack(tool) -command "fdt:select_tool $w" eddy_current "Eddy current correction" bedpostx "BEDPOSTX Estimation of diffusion parameters"  registration "Registration" probtrackx "PROBTRACKX Probabilistic tracking" xutilssx "----------------------------------------------------" dtifit "DTIFIT Reconstruct diffusion tensors" 
    $w.tool.menu.menu entryconfigure 4 -state disabled -background black
    pack $w.tool.menu -side left -pady 3 -padx 6 -anchor nw

    #-------- Tool Options... -------- 

    frame $w.opts

    #------- Registration --------

    frame $w.registration

    proc registration_set_directory { w dirname } {
	global registration

	set struct [ file join $dirname struct_brain ]

	if { [ imtest $struct ] } {
	    set registration(struct_image) $struct
	} else {
	    set registration(struct_image) ""
	}
    }

    FileEntry $w.registration.directory -textvariable registration(directory) -label "BEDPOSTX directory:" -title "Choose directory" -filetypes * -command "registration_set_directory $w" 

    frame       $w.registration.struct
    checkbutton $w.registration.struct.yn -variable registration(struct_yn) -command "registration_packframe $w"
    label       $w.registration.struct.lb -text "Main structural image"
    TitleFrame  $w.registration.struct.tf -text "Main structural image" 
    frame       $w.registration.nonlin
    checkbutton $w.registration.nonlin.yn -variable registration(nonlin_yn) -text "Non-betted structural (for nonlinear reg)" -command "registration_packframe $w"
    FileEntry   $w.registration.nonlin.file -textvariable registration(nonlin_image) -filetypes IMAGE 


    optionMenu2 $w.registration.struct.tf.search registration(struct_search) 0 "No search" 90 "Normal search" 180 "Full search"
    optionMenu2 $w.registration.struct.tf.dof registration(struct_dof)   6 "6 DOF" 7 "7 DOF" 9 "9 DOF" 12 "12 DOF"  
    optionMenu2 $w.registration.struct.tf.costfn registration(struct_costfn) corratio "Correlation ratio" mutualinfo "Mutual information"
    FileEntry   $w.registration.struct.tf.file -textvariable registration(struct_image) -filetypes IMAGE -width 45
    pack $w.registration.struct.tf.file -side top -in [ $w.registration.struct.tf getframe ]
    pack $w.registration.struct.tf.search $w.registration.struct.tf.dof $w.registration.struct.tf.costfn -side left  -in [ $w.registration.struct.tf getframe ]
    set registration(struct_costfn) corratio
    set registration(struct_dof) 6
    set registration(struct_search) 90
    set registration(struct_yn) 0

    frame       $w.registration.standard
    checkbutton $w.registration.standard.yn -variable registration(standard_yn)  -command "registration_packframe $w"
    TitleFrame  $w.registration.standard.tf -text "Standard space"
    label       $w.registration.standard.lb -text "Standard space"
    optionMenu2 $w.registration.standard.tf.search registration(standard_search) 0 "No search" 90 "Normal search" 180 "Full search"
    optionMenu2 $w.registration.standard.tf.dof registration(standard_dof)   6 "6 DOF" 7 "7 DOF" 9 "9 DOF" 12 "12 DOF"  
    optionMenu2 $w.registration.standard.tf.costfn registration(standard_costfn) corratio "Correlation ratio" mutualinfo "Mutual information"
    FileEntry   $w.registration.standard.tf.file -textvariable registration(standard_image) -filetypes IMAGE -width 45 
    pack $w.registration.standard.tf.file -side top -in [ $w.registration.standard.tf getframe ]
    pack $w.registration.standard.tf.search $w.registration.standard.tf.dof $w.registration.standard.tf.costfn -side left -in [ $w.registration.standard.tf getframe ]
    set registration(standard_yn) 1
    set registration(standard_dof) 12
    set registration(standard_search) 90
    set registration(standard_image) [ file join ${FSLDIR} data standard MNI152_T1_2mm_brain ]

    pack $w.registration.directory $w.registration.struct $w.registration.nonlin $w.registration.standard -side top -padx 3 -pady 3 -anchor w



    proc registration_packframe { w } {
       global registration
       pack forget $w.registration.struct.yn $w.registration.struct.tf $w.registration.struct.yn $w.registration.struct.lb $w.registration.nonlin.yn $w.registration.nonlin.file
       pack forget $w.registration.standard.yn $w.registration.standard.tf $w.registration.standard.yn $w.registration.standard.lb
       if { $registration(struct_yn) } { pack $w.registration.struct.yn $w.registration.struct.tf -side left -anchor w } else { pack $w.registration.struct.yn  $w.registration.struct.lb -side left -anchor w}
       if { $registration(struct_yn) } { pack $w.registration.nonlin.yn -side left }
       if { $registration(struct_yn) && $registration(nonlin_yn) } { pack $w.registration.nonlin.file -side left }
       if { $registration(standard_yn) } { pack $w.registration.standard.yn $w.registration.standard.tf -side left -anchor w } else { pack $w.registration.standard.yn  $w.registration.standard.lb -side left -anchor w}
    }
    
    registration_packframe $w
    #------- ECC --------
    frame $w.ecc

    proc ecc_update_files { w filename } {
	global eddy
	set eddy(output) [ file join [file dirname $eddy(input)] data ]
    }    

    FileEntry $w.ecc.input -textvariable eddy(input) -label "Diffusion weighted data:" -title "Choose diffusion weighted image" -filetypes IMAGE -command "ecc_update_files $w"

    FileEntry $w.ecc.output -textvariable eddy(output) 	-label "Corrected output data:" -title  "Choose output image name" -filetypes IMAGE -command "ecc_update_files $w"

    set eddy(refnum) 0
    LabelSpinBox  $w.ecc.refnum -label "Reference volume"  -textvariable eddy(refnum) -range { 0 100 1 } -width 6 
    FileEntry $w.ecc.bvecdata -textvariable eddy(bVecData) 	-label "bvecs file:" -title  "Choose bvecs name" -filetypes IMAGE 

    pack $w.ecc.input $w.ecc.output $w.ecc.refnum -side top -padx 3 -pady 3 -expand yes -anchor w

   #------- DTIFit --------

    frame $w.dtifit

    FileEntry $w.dtifit.directory -textvariable dtifit(directory) -label  "Input directory:" -title "Choose directory" -command "set_working_directory dtifit(cwd)"

    proc dtifit_toggle_expert { w } {
	global dtifit

	if { $dtifit(expert_yn) } {
	    pack forget $w.dtifit.directory
	    pack $w.dtifit.expert -in $w.dtifit -after $w.dtifit.expert_yn
	} else {
	    pack forget $w.dtifit.expert
	    pack $w.dtifit.directory -in $w.dtifit -before $w.dtifit.expert_yn
	}
    }

    checkbutton $w.dtifit.expert_yn -text "Specify input files manually" \
	-variable dtifit(expert_yn) -command "dtifit_toggle_expert $w"

    frame $w.dtifit.expert

    proc set_working_directory { cwd filename } {
	global dtifit
	set dirname [file dirname $filename]
	puts "switching from $dtifit(cwd) to $dirname" 
	set dtifit(cwd) $dirname
    }

    proc dtifit_update_files { w filename } {
	global dtifit

	set dtifit(output) [ file join [file dirname $dtifit(input)] dti ]
	set_working_directory dtifit(cwd) $dtifit(input)
    }
    
    set dtifit(cwd) [ pwd ]

#All the below orignally had -directory $dtifit(cwd) 
    option add *dtifit.expert.FileEntry*labf*width 27
    FileEntry $w.dtifit.expert.input -textvariable dtifit(input) -label  "Diffusion weighted data:" -title "Choose diffusion weighted image" -filetypes IMAGE -command "dtifit_update_files $w" 
    FileEntry $w.dtifit.expert.mask -textvariable dtifit(mask) -label "BET binary brain mask:" -title "Choose BET brain mask file" -filetypes IMAGE -command "set_working_directory dtifit(cwd)"
    FileEntry $w.dtifit.expert.output -textvariable dtifit(output) -label "Output basename:" -title  "Choose output base name" -filetypes * -command "set_working_directory dtifit(cwd)"
    FileEntry $w.dtifit.expert.bvecs -textvariable dtifit(bvecs) -label "Gradient directions:" -title  "Choose bvecs file" -filetypes * -command "set_working_directory dtifit(cwd)"
    FileEntry $w.dtifit.expert.bvals -textvariable dtifit(bvals) -label  "b values:" -title  "Choose bvals file" -command "set_working_directory dtifit(cwd)"

    pack $w.dtifit.expert.input $w.dtifit.expert.mask $w.dtifit.expert.output \
	$w.dtifit.expert.bvecs $w.dtifit.expert.bvals \
	-side top -padx 3 -pady 3 -expand yes -anchor w

    pack $w.dtifit.directory $w.dtifit.expert_yn -side top -padx 3 -pady 3 -expand yes -anchor w

   collapsible frame $w.dtifit.advanced -title "Advanced Options"

    checkbutton $w.dtifit.advanced.w -text "Weighted Least-Squares" -variable dtifit(doWLS)  
    checkbutton $w.dtifit.advanced.sse -text "Save Sum-squared Error" -variable dtifit(doSaveSSE)  
    checkbutton $w.dtifit.advanced.tensor     -text "Save tensor" -variable dtifit(doSaveTensor)  

    pack $w.dtifit.advanced.w $w.dtifit.advanced.sse $w.dtifit.advanced.tensor -in $w.dtifit.advanced.b -anchor w
    pack $w.dtifit.directory $w.dtifit.advanced -side top -padx 3 -pady 3 -expand yes -anchor w


    #------- BEDPOST --------

    frame $w.bedpost

    FileEntry $w.bedpost.directory -textvariable bedpost(directory) -label "Input directory:" -title "Choose directory" -filetypes * -command "set_working_directory dtifit(cwd)"

   collapsible frame $w.bedpost.advanced -title "Advanced Options"
   set bedpost(nfibres) 2
   set bedpost(weight)  1
   set bedpost(burnin)  1000
   LabelSpinBox  $w.bedpost.advanced.nfibres -label "Fibres  "  -textvariable bedpost(nfibres) -range {1 1000000000 1 } 
   LabelSpinBox  $w.bedpost.advanced.weight  -label "Weight  "  -textvariable bedpost(weight) -range {0.0 100000000.0 1 } 
   LabelSpinBox  $w.bedpost.advanced.burnin  -label "Burn In "  -textvariable bedpost(burnin) -range {1 1000000000 1 } 

    checkbutton $w.bedpost.advanced.multishell -text "Multi-Shell model" -variable bedpost(useMultiShell)  
    checkbutton $w.bedpost.advanced.noisefloor -text "Model Noise Floor" -variable bedpost(useNoiseFloor)  
    checkbutton $w.bedpost.advanced.rician     -text "Rician Noise" -variable bedpost(useRician)  



    set bedpost(ecc_yn) 0
    pack $w.bedpost.advanced.nfibres $w.bedpost.advanced.weight $w.bedpost.advanced.burnin $w.bedpost.advanced.multishell $w.bedpost.advanced.noisefloor $w.bedpost.advanced.rician  -in $w.bedpost.advanced.b -anchor w
    pack $w.bedpost.directory $w.bedpost.advanced -side top -padx 3 -pady 3 -expand yes -anchor w

    #-------- ProbTrackX -------- 
    NoteBook $w.probtrack -bd 2 -tabpady {5 10} -arcradius 3
    $w.probtrack insert 0 data -text "Data"
    $w.probtrack insert 1 options -text "Options"
    #-------- Mode specific option --------
    frame $w.data
    FileEntry $w.data.directory -textvariable probtrack(bedpost_dir) -label "BEDPOSTX directory" -title "Choose BEDPOSTX directory" -filetypes { } -command "probtrack_update_files $w"

    TitleFrame  $w.data.seed -text "Seed Space"
 
    optionMenu2 $w.data.seed.menu probtrack(mode) -command "fdt:probtrack_mode $w" simple "Single voxel" seedmask "Single mask" network "Multiple masks"
    pack $w.data.seed.menu -in $w.data.seed.f -side top -anchor w -pady 2

    set probtrack(x) 0
    set probtrack(y) 0
    set probtrack(z) 0
    set probtrack(units) vox
    #Co-ordinate edit frame
    frame $w.data.seed.voxel
    LabelSpinBox $w.data.seed.voxel.x -label "X" -textvariable probtrack(x) -range {-1000000 1000000 1 } 
    LabelSpinBox $w.data.seed.voxel.y -label "Y" -textvariable probtrack(y) -range {-1000000 1000000 1 } 
    LabelSpinBox $w.data.seed.voxel.z -label "Z" -textvariable probtrack(z) -range {-1000000 1000000 1 } 
    radiobutton $w.data.seed.voxel.vox -text "vox" -value vox -variable probtrack(units)
    radiobutton $w.data.seed.voxel.mm  -text "mm"  -value mm  -variable probtrack(units)
    
   
 
    option add *seed*FileEntry*labf*width 24

    frame  $w.data.seed.ssf
    set probtrack(mode) simple
    checkbutton $w.data.seed.ssf.ssd -text "Seed space is not diffusion" -variable probtrack(usereference_yn)  -command " fdt:probtrack_mode $w "
    checkbutton $w.data.seed.ssf.nonlinear -text "nonlinear" -variable probtrack(useNonlinear)  -command " fdt:probtrack_mode $w "
    checkbutton $w.data.seed.ssf.useSurface -text "surface" -variable probtrack(useSurface)  -command " fdt:probtrack_mode $w "
    FileEntry $w.data.seed.ssf.xfm -textvariable probtrack(xfm)  -label "Select Seed to diff transform" -title "Select seed-space to DTI-space transformation matrix" -filetypes *
    FileEntry $w.data.seed.ssf.invxfm -textvariable probtrack(invxfm)  -label "Select diff to Seed transform" -title "Select seed-space to DTI-space transformation matrix" -filetypes *

    labelframe $w.data.seed.ssf.typelabel -text "Mesh convention" -labelanchor w -relief flat
    optionMenu2  $w.data.seed.ssf.typelabel.type probtrack(meshspace) caret "Caret" freesurfer "FreeSurfer" first "FIRST" vox "Voxel"
    FileEntry $w.data.seed.ssf.surfref -textvariable probtrack(surfref)  -label "Surface Reference Image" -title "Select surface reference volume" -filetypes IMAGE

    FileEntry $w.data.seed.ssf.reference -textvariable probtrack(reference) -label "Seed Image/Surface" -title "Choose Image/Surface" -filetypes "{*.nii,*.nii.gz,*.gii}"
    pack $w.data.seed.ssf.ssd -side top -anchor nw
    pack $w.data.seed.voxel.x $w.data.seed.voxel.y $w.data.seed.voxel.z $w.data.seed.voxel.vox $w.data.seed.voxel.mm -side left -padx 2
    pack $w.data.seed.voxel $w.data.seed.ssf -in $w.data.seed.f -side left -anchor w -pady 2
    


    TitleFrame $w.data.seed.target -text "Masks list"    
    listbox $w.data.seed.targets -height 6 -width 50 -yscrollcommand "$w.data.seed.sb set"
    scrollbar $w.data.seed.sb -command "$w.data.seed.targets yview " 
    frame $w.data.seed.tb
    button $w.data.seed.tb.add -text "Add Image"  -command "feat_file:setup_dialog $w a a a [namespace current] {{*.nii,*.nii.gz,*.gii}} {Select File} {fdt_add $w $w.data.seed.targets} {}" 
    button $w.data.seed.tb.del -text "Remove Image"  -command "fdt_sub $w $w.data.seed.targets" 
    button $w.data.seed.tb.imp -text "Load List" -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_imp $w $w.data.seed.targets} {}"
    button $w.data.seed.tb.exp -text "Save List" -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_exp $w $w.data.seed.targets} {}"
    pack $w.data.seed.tb.add $w.data.seed.tb.del $w.data.seed.tb.imp $w.data.seed.tb.exp -side left
    pack $w.data.seed.tb -in [$w.data.seed.target getframe ] -side bottom  -expand yes -fill x -anchor w -padx 3 -pady 3
    pack $w.data.seed.targets $w.data.seed.sb -in [$w.data.seed.target getframe ] -side left  -expand yes -fill y -anchor w -padx 3 -pady 3
    
    TitleFrame  $w.data.targets -text "Optional Targets"

    proc fdt_add { w listbox filename } {
    set filename [ fix_cygwin_filename $filename ]
    $listbox insert end $filename
    }

    proc fdt_sub { w listbox} {
    set count 0
    foreach file [ $listbox get 0 end ] {
	if { [ $listbox selection includes $count ] == 1 } {
	    $listbox delete $count
	    incr count -1
	}
	incr count
    } 
    }

    proc fdt_imp { w listbox filename } {
    if { ![ file exists $filename ] } {
	MxPause "Warning: Bad or missing file!"
	return
    }
    set fd [ open $filename ]
    $listbox  delete 0 end
    while { [ gets $fd file ] >= 0 } {
	$listbox insert end $file
    }
    close $fd
    }

    proc fdt_exp { w listbox filename } {
    set fd [ open $filename w ]
    foreach file [ $listbox get 0 end ] {
	puts $fd $file
    }
    close $fd
    }

    frame $w.data.targets.wf    
    checkbutton $w.data.targets.wf.sct -text "Waypoints masks" -variable probtrack(waypoint_yn)  -command " pack forget $w.data.targets.wf.tf ; if { \$probtrack(waypoint_yn) } { pack $w.data.targets.wf.tf } ; $w.probtrack compute_size"
    TitleFrame $w.data.targets.wf.tf -text "Waypoints list"    
    listbox $w.data.targets.wf.tf.targets -height 6 -width 50 -yscrollcommand "$w.data.targets.wf.tf.sb set"
    scrollbar $w.data.targets.wf.tf.sb -command "$w.data.targets.wf.tf.targets yview " 
    frame $w.data.targets.wf.tf.tb
    button $w.data.targets.wf.tf.tb.add -text "Add Image"  -command "feat_file:setup_dialog $w a a a [namespace current] {{*.nii,*.nii.gz,*.gii}} {Select File} {fdt_add $w $w.data.targets.wf.tf.targets} {}"
    button $w.data.targets.wf.tf.tb.del -text "Remove Image"  -command "fdt_sub $w $w.data.targets.wf.tf.targets" 
    button $w.data.targets.wf.tf.tb.imp -text "Load List" -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_imp $w $w.data.targets.wf.tf.targets} {}"
    button $w.data.targets.wf.tf.tb.exp -text "Save List" -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_exp $w $w.data.targets.wf.tf.targets} {}"
    pack $w.data.targets.wf.tf.tb.add $w.data.targets.wf.tf.tb.del $w.data.targets.wf.tf.tb.imp $w.data.targets.wf.tf.tb.exp -side left
    pack $w.data.targets.wf.tf.tb -in [$w.data.targets.wf.tf getframe ] -side bottom  -expand yes -fill x -anchor w -padx 3 -pady 3
    pack $w.data.targets.wf.tf.targets $w.data.targets.wf.tf.sb -in [$w.data.targets.wf.tf getframe ] -side left  -expand yes -fill y -anchor w -padx 3 -pady 3
    pack  $w.data.targets.wf.sct -side top -anchor nw
    pack $w.data.targets.wf 
  
    option add *targets*Checkbutton*width 18
    option add *targets*Checkbutton*anchor w
    frame  $w.data.targets.ef
    checkbutton $w.data.targets.ef.srt -text "Exclusion mask" -variable probtrack(exclude_yn)  -command " pack forget $w.data.targets.ef.rubbish ; if { \$probtrack(exclude_yn) } { pack $w.data.targets.ef.rubbish } ; $w.probtrack compute_size"
    FileEntry $w.data.targets.ef.rubbish -textvariable probtrack(exclude) -title "Select exclusion image" -filetypes "{*.nii,*.nii.gz,*.gii}"
    pack $w.data.targets.ef.srt -side left

    frame  $w.data.targets.sf
    checkbutton $w.data.targets.sf.sst -text "Termination mask" -variable probtrack(terminate_yn)  -command " pack forget $w.data.targets.sf.stop ; if { \$probtrack(terminate_yn) } { pack $w.data.targets.sf.stop } ; $w.probtrack compute_size"
    FileEntry $w.data.targets.sf.stop -textvariable probtrack(stop) -title "Select termination image" -filetypes "{*.nii,*.nii.gz,*.gii}"
    pack $w.data.targets.sf.sst -side left

    frame $w.data.targets.cf    
    checkbutton $w.data.targets.cf.sct -text "Classification targets" -variable probtrack(classify_yn)  -command " pack forget $w.data.targets.cf.tf ; if { \$probtrack(classify_yn) } { pack $w.data.targets.cf.tf } ; $w.probtrack compute_size"
    TitleFrame $w.data.targets.cf.tf -text "Targets list"    
    listbox $w.data.targets.cf.tf.targets -height 6 -width 50 -yscrollcommand "$w.data.targets.cf.tf.sb set"
    scrollbar $w.data.targets.cf.tf.sb -command "$w.data.targets.cf.tf.targets yview " 
    frame $w.data.targets.cf.tf.tb
    button $w.data.targets.cf.tf.tb.add -text "Add Image"  -command "feat_file:setup_dialog $w a a a [namespace current] {{*.nii,*.nii.gz,*.gii}} {Select File} {fdt_add $w $w.data.targets.cf.tf.targets} {}"
    button $w.data.targets.cf.tf.tb.del -text "Remove Image"  -command "fdt_sub $w $w.data.targets.cf.tf.targets" 
    button $w.data.targets.cf.tf.tb.imp -text "Load List" -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_imp $w $w.data.targets.cf.tf.targets} {}"
    button $w.data.targets.cf.tf.tb.exp -text "Save List" -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_exp $w $w.data.targets.cf.tf.targets} {}"
    pack $w.data.targets.cf.tf.tb.add $w.data.targets.cf.tf.tb.del $w.data.targets.cf.tf.tb.imp $w.data.targets.cf.tf.tb.exp -side left
    pack $w.data.targets.cf.tf.tb -in [$w.data.targets.cf.tf getframe ] -side bottom  -expand yes -fill x -anchor w -padx 3 -pady 3
    pack $w.data.targets.cf.tf.targets $w.data.targets.cf.tf.sb -in [$w.data.targets.cf.tf getframe ] -side left  -expand yes -fill y -anchor w -padx 3 -pady 3
    pack  $w.data.targets.cf.sct -side top -anchor nw
    pack $w.data.targets.cf 


    pack $w.data.targets.wf $w.data.targets.ef $w.data.targets.sf $w.data.targets.cf -in $w.data.targets.f -anchor w




    set probtrack(xfm) ""
    set probtrack(basename) "merged"
    set probtrack(mask) "nodif_brain_mask"

    proc probtrack_update_files { w filename } {
	global probtrack
	global FSLDIR

	if { ($probtrack(bedpost_dir) != "") && ($probtrack(reference) != "") } {
	    set probtrack(output) \
		[ file join $probtrack(bedpost_dir) [ file tail [ exec $FSLDIR/bin/remove_ext $probtrack(reference) ] ] ]
	}
    }

    FileEntry $w.data.dir -textvariable probtrack(output) -label  "Output directory:" -title  "Name the output directory" -filetypes { }

    pack $w.data.directory $w.data.seed $w.data.targets $w.data.dir -padx 3 -pady 3 -anchor nw

    pack $w.data -in  [$w.probtrack getframe data] -padx 3 -pady 3 -anchor nw -expand yes -fill both

    #-------- ...Options... --------
    TitleFrame $w.options -text "Basic Options"

    checkbutton $w.options.verbose -text "Verbose" -variable probtrack(verbose_yn)
    
    set probtrack(nparticles) 5000
    LabelSpinBox $w.options.nparticles -label  "Number of samples" -textvariable probtrack(nparticles) -range {1 1e24 100 } -width 6 

    set probtrack(curvature) 0.2
    LabelSpinBox $w.options.curvature -label "Curvature threshold" -textvariable probtrack(curvature) -range {0.0 1.0 0.01 }

    set probtrack(loopcheck_yn) 1
    checkbutton $w.options.loopcheck -text "Loopcheck" -variable probtrack(loopcheck_yn)

    collapsible frame $w.advanced -title "Advanced Options" -command "$w.probtrack compute_size; set dummy" 

    set probtrack(nsteps) 2000
    LabelSpinBox $w.advanced.nsteps -label "Maximum number of steps" -textvariable probtrack(nsteps) -range {2 1000000 10 } -width 6

    set probtrack(steplength) 0.5
    LabelSpinBox $w.advanced.steplength -label "Step length (mm)" -textvariable probtrack(steplength) -range {0.001 10000.0 0.1} 

    set probtrack(modeuler_yn) 0
    checkbutton $w.advanced.modeuler -text "Use modified Euler streamlining" -variable probtrack(modeuler_yn)

    set probtrack(pd) 0
    checkbutton $w.advanced.pd -text "Use Distance correction" -variable probtrack(pd)

    set probtrack(usef_yn) 0
    checkbutton $w.advanced.usef -text "Use anisotropy to constrain tracking" -variable probtrack(usef_yn)

    set probtrack(fibthresh) 0.01
    LabelSpinBox $w.advanced.fibthresh -label "Subsidary fibre volume fraction threshold" -textvariable probtrack(fibthresh) -range {0.00 1.0 0.01} 

    set probtrack(distthresh) 0.0
    LabelSpinBox $w.advanced.distthresh -label "Minimum length threshold (mm)" -textvariable probtrack(distthresh) -range {0.00 1000.0 1.0} 

    set probtrack(sampvox) 0.0
    LabelSpinBox $w.advanced.sampvox -label "Seed sphere sampling (mm)" -textvariable probtrack(sampvox) -range {0.00 100.0 1.0} 


    collapsible frame $w.wayadvanced -title "Waypoint Options" -command "$w.probtrack compute_size; set dummy"

    set probtrack(oneway_yn) 1
    checkbutton $w.wayadvanced.oneway_yn -text "Apply waypoint independently to both directions" -variable probtrack(oneway_yn)

    set probtrack(wayorder_yn) 0
    checkbutton $w.wayadvanced.wayorder_yn -text "Force waypoint crossing in listed order" -variable probtrack(wayorder_yn)

    set probtrack(waycond) AND

    labelframe $w.wayadvanced.waylabel -text "Waypoint condition" -labelanchor w -relief flat
    optionMenu2 $w.wayadvanced.waylabel.waycond probtrack(waycond) AND "AND" OR "OR"
    pack $w.wayadvanced.waylabel $w.wayadvanced.waylabel.waycond -side top -anchor w -pady 2

    collapsible frame $w.matadvanced -title "Matrix Options" -command "$w.probtrack compute_size; set dummy"

    # Matrix Stuff
    set probtrack(omatrix1_yn) 0
    checkbutton $w.matadvanced.omatrix1_yn -text "Matrix1: Seed x Seed Matrix" -variable probtrack(omatrix1_yn)

    set probtrack(omatrix2_yn) 0
    checkbutton $w.matadvanced.omatrix2_yn -text "Matrix2: Seed x Mask2 Matrix" -variable probtrack(omatrix2_yn) -command "fdt:matrix_mode $w"
    FileEntry $w.matadvanced.omatrix2_mask -textvariable probtrack(mask22) -label "Tract space mask" -title "Select tract space mask" -filetypes "{*.nii,*.nii.gz,*.gii}"

    set probtrack(omatrix3_yn) 0
    checkbutton $w.matadvanced.omatrix3_yn -text "Matrix3: Mask1 x Mask2 Matrix" -variable probtrack(omatrix3_yn) -command "fdt:matrix_mode $w"
    FileEntry $w.matadvanced.omatrix3_mask1 -textvariable probtrack(mask31) -label "Row space mask" -title "Select row space mask" -filetypes "{*.nii,*.nii.gz,*.gii}"
    FileEntry $w.matadvanced.omatrix3_mask2 -textvariable probtrack(mask32) -label "Column space mask" -title "Select column space mask" -filetypes "{*.nii,*.nii.gz,*.gii}"

    pack \
	$w.options.nparticles \
	$w.options.curvature \
	$w.options.verbose \
	$w.options.loopcheck \
	-in [$w.options getframe ] -side top -pady 3 -padx 6 -anchor nw


    pack \
	$w.advanced.modeuler \
	$w.advanced.nsteps \
	$w.advanced.steplength \
	$w.advanced.usef \
	$w.advanced.pd \
	$w.advanced.fibthresh \
	$w.advanced.distthresh \
	$w.advanced.sampvox \
	-in $w.advanced.b  -side top -pady 3 -padx 6 -anchor nw

    pack \
	$w.wayadvanced.oneway_yn \
	$w.wayadvanced.wayorder_yn \
	$w.wayadvanced.waylabel \
	-in $w.wayadvanced.b -side top -pady 3 -padx 6 -anchor nw


    pack \
	$w.matadvanced.omatrix1_yn \
	$w.matadvanced.omatrix2_yn \
	$w.matadvanced.omatrix3_yn \
	-in $w.matadvanced.b -side top -pady 3 -padx 6 -anchor nw

    pack \
	$w.options \
	$w.advanced \
	$w.wayadvanced \
	$w.matadvanced \
	-in [$w.probtrack getframe options] -side top -pady 3 -padx 6 -anchor nw -expand yes -fill both


    #-------- Buttons --------

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.apply     -command "fdt:apply $w keep" \
        -text "Go" -width 5
    bind $w.apply <Return> {
        [winfo toplevel %W].apply invoke}
 
    button $w.cancel    -command "fdt:destroy $w" \
        -text "Exit" -width 5
    bind $w.cancel <Return> {
        [winfo toplevel %W].cancel invoke}

    button $w.help -command "FmribWebHelp file: $HTMLPATH/fdt.html" \
	    -text "Help" -width 5
    bind $w.help <Return> {
	[winfo toplevel %W].help invoke}
 
    pack $w.btns.b -side bottom -fill x -padx 3 -pady 5
    pack $w.apply $w.cancel $w.help -in $w.btns.b \
        -side left -expand yes -padx 3 -pady 10 -fill y
 
    pack $w.tool $w.opts $w.btns -side top -expand yes -fill both
    
 

    $w.probtrack raise data 
    fdt:select_tool $w 
    
    fdt:probtrack_mode $w

    update idletasks
    if { $tclstartupfile != "" } {
	puts "Reading $tclstartupfile"
	source $tclstartupfile
	fdt:select_tool $w 
	fdt:probtrack_mode $w
    }
}

proc fdt:eddycorrect_mode { w } {
    global eddy 
    pack forget $w.ecc.bvecdata
}

proc fdt:matrix_mode { w } {
    global probtrack

    pack forget $w.matadvanced.omatrix2_mask $w.matadvanced.omatrix3_mask1 $w.matadvanced.omatrix3_mask2

    if { $probtrack(omatrix2_yn) } { pack $w.matadvanced.omatrix2_yn $w.matadvanced.omatrix2_mask $w.matadvanced.omatrix3_yn -in $w.matadvanced.b -anchor w -pady 2 }
    if { $probtrack(omatrix3_yn) } { pack $w.matadvanced.omatrix3_yn $w.matadvanced.omatrix3_mask1 $w.matadvanced.omatrix3_mask2 -in $w.matadvanced.b -anchor w -pady 2 }
    
    $w.probtrack compute_size
}

proc fdt:probtrack_mode { w } {
    global probtrack FSLDIR

    pack forget $w.data.seed.voxel $w.data.seed.ssf  $w.data.seed.ssf.xfm $w.data.seed.ssf.reference $w.data.seed.bcf $w.data.seed.target $w.data.targets.cf $w.data.seed.ssf.invxfm $w.data.seed.ssf.nonlinear $w.data.seed.ssf.meshspace $w.data.seed.ssf.typelabel $w.data.seed.ssf.useSurface  $w.data.seed.ssf.surfref
    $w.data.dir configure -label  "Output directory:" -title  "Name the output directory" -filetypes *
    pack $w.data.seed.ssf -in $w.data.seed.f -side bottom -anchor w -pady 2

    switch -- $probtrack(mode) {
  	simple {
	    pack $w.data.seed.voxel -in $w.data.seed.f -side bottom -anchor w -pady 2
	    if { $probtrack(usereference_yn) } { pack $w.data.seed.ssf.reference -side bottom -anchor w -pady 2 }
	    $w.data.seed.ssf.reference configure -label "Seed reference image:" -title "Choose reference image" 
	    $w.data.dir configure -label  "Output file:" -title  "Name the output file" -filetypes IMAGE
	    if { $probtrack(usereference_yn) } { pack $w.data.seed.ssf.nonlinear $w.data.seed.ssf.xfm -side top -anchor w -pady 2 }
	    if { $probtrack(useNonlinear) && $probtrack(usereference_yn) } { pack $w.data.seed.ssf.xfm $w.data.seed.ssf.invxfm -side top -anchor w -pady 2 }	    
    	}
	seedmask {
            pack forget $w.data.seed.ssf.ssd
	    pack $w.data.targets.cf -in $w.data.targets.f -anchor w
	    $w.data.seed.ssf.reference configure -label "Seed Image/Surface:" -title "Choose Image/Surface" 
            pack $w.data.seed.ssf.reference  $w.data.seed.ssf.ssd -side top -anchor w -pady 2

	    if { $probtrack(usereference_yn) } { pack $w.data.seed.ssf.nonlinear $w.data.seed.ssf.xfm -side top -anchor w -pady 2 }
	    if { $probtrack(useNonlinear) && $probtrack(usereference_yn) } { pack $w.data.seed.ssf.xfm $w.data.seed.ssf.invxfm -side top -anchor w -pady 2 }	    

	    if { $probtrack(usereference_yn) } { pack $w.data.seed.ssf.useSurface -side top -anchor w -pady 2 }
	    if { $probtrack(useSurface) && $probtrack(usereference_yn) } { pack $w.data.seed.ssf.typelabel $w.data.seed.ssf.typelabel.type -side top -anchor w -pady 2 }
	    if { $probtrack(useSurface) && $probtrack(usereference_yn) } { pack $w.data.seed.ssf.surfref -side top -anchor w -pady 2 }
  	}
	network {
	    if { $probtrack(usereference_yn) } { pack $w.data.seed.ssf.nonlinear $w.data.seed.ssf.xfm -side top -anchor w -pady 2 }
	    if { $probtrack(useNonlinear) && $probtrack(usereference_yn) } { pack $w.data.seed.ssf.xfm $w.data.seed.ssf.invxfm -side top -anchor w -pady 2 }	    

	    if { $probtrack(usereference_yn) } { pack $w.data.seed.ssf.useSurface -side top -anchor w -pady 2 }
	    if { $probtrack(useSurface) && $probtrack(usereference_yn) } { pack $w.data.seed.ssf.typelabel $w.data.seed.ssf.typelabel.type -side top -anchor w -pady 2 }
	    if { $probtrack(useSurface) && $probtrack(usereference_yn) } { pack $w.data.seed.ssf.surfref -side top -anchor w -pady 2 }
	    pack  $w.data.seed.target -in $w.data.seed.f -side bottom -anchor w -pady 2
	}
    }
    if { $probtrack(waypoint_yn) } { pack $w.data.targets.wf.tf } 
    if { $probtrack(classify_yn) } { pack $w.data.targets.cf.tf }

    $w.probtrack compute_size
}

proc fdt:select_tool { w } {
    global probtrack
    pack forget $w.ecc
    pack forget $w.probtrack
    pack forget $w.bedpost
    pack forget $w.registration
    pack forget $w.dtifit
    if {$probtrack(tool) == "bedpostx"} { 
	pack $w.bedpost -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    } elseif {$probtrack(tool) == "probtrackx"}  { 
	pack $w.probtrack -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    } elseif {$probtrack(tool) == "dtifit"}  { 
	pack $w.dtifit -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    } elseif {$probtrack(tool) == "eddy_current"}  { 
	pack $w.ecc -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    } elseif {$probtrack(tool) == "registration"} {
	pack $w.registration -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    }
}
proc fdt_monitor_short { w cmd } {
    global debugging OSFLAVOUR FSLPARALLEL

    puts "$cmd"

    if { $OSFLAVOUR != "cygwin" } {
	set oldcursor [ $w configure -cursor { watch red white } ]

	catch {
	    update idletasks
	    if { ! $debugging } {
		set fd [ open "|$cmd" r ]
#		set fd [ open "|qrsh -V -now n -q short.q $cmd" r ]
		while { ( [ gets $fd line ] >= 0 ) } {
		    update idletasks
		    puts $line
		}
		close $fd
	    }
	} junk

	$w configure -cursor $oldcursor

    } else {
	catch { exec sh -c $cmd } junk
    }

    if { $junk != "" } {
	MxPause "Errors: $junk"
    } 

    puts "Done!"
}

proc fdt_monitor { w cmd } {
    global debugging OSFLAVOUR

    puts "$cmd"

    if { $OSFLAVOUR != "cygwin" } {
	set oldcursor [ $w configure -cursor { watch red white } ]

	catch {
	    update idletasks
	    if { ! $debugging } {
		set fd [ open "|$cmd" r ]
#		set fd [ open "|qrsh -V -now n -q long.q $cmd" r ]
		while { ( [ gets $fd line ] >= 0 ) } {
		    update idletasks
		    puts $line
		}
		close $fd
	    }
	} junk

	$w configure -cursor $oldcursor

    } else {
	catch { exec sh -c $cmd } junk
    }

    if { $junk != "" } {
	MxPause "Errors: $junk"
    } 

    puts "Done!"
}

proc fdt:apply { w dialog } {
    global probtrack BINPATH FSLDIR FSLPARALLEL

    switch -- $probtrack(tool) {
	eddy_current {
	    global eddy

	    set errorStr ""
	    if { $eddy(input) == "" } { set errorStr "You need to specify the input image! " }
	    if { $eddy(output) == "" } { set errorStr "$errorStr You need to specify an output image!" }
	    if { $errorStr != "" } {
		MxPause $errorStr
		return
	    }

	    #	    check output!=input
	    set canwrite 1
	    if { $eddy(input) == $eddy(output) } {
		set canwrite [ YesNoWidget "Output and input images have the same name. Overwrite input?" Yes No ]
	    }
	    if { $canwrite } {
		fdt_monitor $w "${FSLDIR}/bin/eddy_correct $eddy(input) $eddy(output) $eddy(refnum) "
	    }
	}
	dtifit {
	    global dtifit

	    if { ! $dtifit(expert_yn) } {
		set dtifit(input)  [ file join $dtifit(directory) data ]
		set dtifit(output) [ file join $dtifit(directory) dti ]
		set dtifit(mask)   [ file join $dtifit(directory) nodif_brain_mask ]
		set dtifit(bvecs)  [ file join $dtifit(directory) bvecs ]
		set dtifit(bvals)  [ file join $dtifit(directory) bvals ]
	    }

	    set errorStr ""
	    if { $dtifit(directory) == "" && ! $dtifit(expert_yn) } { set errorStr "You must specify the input directory!" }
	    if { $dtifit(input) == "" } { set errorStr "You need to specify the diffusion weighted data image!" }
	    if { $dtifit(output) == "" } { set errorStr "$errorStr You need to specify the output basename!" }
	    if { $dtifit(mask) == "" } { set errorStr "$errorStr You need to specify a mask image!" }
	    if { $dtifit(bvecs) == "" } { set errorStr "$errorStr Please select a gradient directions file!" }
	    if { $dtifit(bvals) == "" } { set errorStr "$errorStr Please select a b values file!" }
	    if { $errorStr != "" } {
		MxPause $errorStr
		return
	    }

	    set canwrite 1
	    if { [file exists $dtifit(output) ] } {
		set canwrite [ YesNoWidget "Overwrite $dtifit(output)?" Yes No ]
	    }
	    if { $canwrite } {
		set flags "--data=$dtifit(input) --out=$dtifit(output) --mask=$dtifit(mask) --bvecs=$dtifit(bvecs) --bvals=$dtifit(bvals)"
		if { $dtifit(doWLS) } { set flags "$flags --wls" }
		if { $dtifit(doSaveSSE) } { set flags "$flags --sse" }
		if { $dtifit(doSaveTensor) } { set flags "$flags --save_tensor" }
		
		fdt_monitor_short $w "${FSLDIR}/bin/dtifit $flags"
	    }
	}
	bedpostx {
	    global bedpost

	    set errorStr ""
	    if { $bedpost(directory) == ""  } { set errorStr "You must specify the bedpost directory!" }
	    if { $errorStr != "" } {
		MxPause $errorStr
		return
	    }

	    set canwrite 1
	    if { [file exists ${bedpost(directory)}.bedpost ] } {
		set canwrite [ YesNoWidget "Overwrite ${bedpost(directory)}.bedpostX?" Yes No ]
		if { $canwrite } {
		    puts "rm -rf ${bedpost(directory)}.bedpostX"
		    catch { exec rm -rf ${bedpost(directory)}.bedpost } errmsg
		}
	    }
	    if { $canwrite } {
		set flags  "$bedpost(directory) --nf=$bedpost(nfibres) --fudge=$bedpost(weight)  --bi=$bedpost(burnin)"
		if { $bedpost(useMultiShell) } { set flags "$flags --model=2" }
		if { $bedpost(useNoiseFloor) } { set flags "$flags --f0 --ardf0" }
		if { $bedpost(useRician) } { set flags "$flags --rician" }
		
		puts "bedpostx $flags"
                
                set filebase $bedpost(directory)/bedpostcom
	        set logfile "${filebase}_log.tcl"
	        set log [open "$logfile" w]
	        puts $log "bedpostx $flags"
                close $log

		fdt_monitor $w "${FSLDIR}/bin/bedpostx $flags"
	    }
	}
	probtrackx {
	    global probtrack env
	    set errorStr ""
            set FSLPARALLEL 0
            if { [ info exists env(SGE_ROOT) ] && $env(SGE_ROOT) != "" } { set FSLPARALLEL 1 }
	    if { $probtrack(bedpost_dir) == ""  } { set errorStr "You must specify the bedpostX directory!" }
	    if { $probtrack(mode) == "simple" && $probtrack(usereference_yn) && $probtrack(reference) == "" } { set errorStr "$errorStr You must specify a reference image" } 
	    if { $probtrack(mode) == "seedmask" && $probtrack(reference) == "" } { set errorStr "$errorStr You must specify a mask image" } 
	    if { $probtrack(exclude_yn) && $probtrack(exclude) == "" } { set errorStr "$errorStr You must specify the exclusion mask!" }
	    if { $probtrack(useNonlinear) && $probtrack(usereference_yn) && $probtrack(invxfm) == "" } { set errorStr "$errorStr You must specify the inverse transform!" }
	    if { $probtrack(useSurface) && $probtrack(usereference_yn) && $probtrack(surfref) == "" } { set errorStr "$errorStr You must specify the surface reference!" }
            if { $probtrack(terminate_yn) && $probtrack(stop) == ""} { set errorStr "$errorStr You must specify the termination mask!" }
	    if { $probtrack(output) == ""  } { set errorStr "$errorStr You must specify the output basename!" }
	    set flags ""
	    if { $probtrack(verbose_yn) == 1 } { set flags "$flags -V 1" }
	    if { $probtrack(loopcheck_yn) == 1 } { set flags "$flags -l" }
	    if { $probtrack(usef_yn) == 1 } { set flags "$flags -f" }
	    if { $probtrack(modeuler_yn) == 1 } { set flags "$flags --modeuler" }
	    if { $probtrack(oneway_yn) == 1 } { set flags "$flags --onewaycondition" }
	    if { $probtrack(wayorder_yn) == 1 } { set flags "$flags --wayorder" }

	    if { $probtrack(omatrix1_yn) == 1 } { set flags "$flags --omatrix1" }
	    if { $probtrack(omatrix2_yn) && $probtrack(mask22) == "" } { set errorStr "$errorStr You must specify a tract space mask!" }
	    if { $probtrack(omatrix2_yn) == 1 } { set flags "$flags --omatrix2 --target2=$probtrack(mask22)" }
	    if { $probtrack(omatrix3_yn) && $probtrack(mask31) == "" } { set errorStr "$errorStr You must specify a row space mask!" }
	    if { $probtrack(omatrix3_yn) && $probtrack(mask32) == "" } { set errorStr "$errorStr You must specify a column space mask!" }
	    if { $probtrack(omatrix3_yn) == 1 } { set flags "$flags --omatrix3 --target3=$probtrack(mask31) --lrtarget3=$probtrack(mask32)" }

            if { $probtrack(pd) } { set flags "$flags --pd"  }
	    set flags "$flags -c $probtrack(curvature) -S $probtrack(nsteps) --steplength=$probtrack(steplength) -P $probtrack(nparticles)"
	    set flags "$flags --fibthresh=$probtrack(fibthresh)"
	    set flags "$flags --distthresh=$probtrack(distthresh)"
	    set flags "$flags --sampvox=$probtrack(sampvox)"

	    if { $errorStr != "" } {
       		MxPause $errorStr
       		return
      	    }
	    set canwrite 1
      	    if { [ file exists $probtrack(output) ] } {
      		set canwrite [  YesNoWidget "Overwrite $probtrack(output)?" Yes No ]
	    }

	    set deslashedFileName [ rmSlash  $probtrack(output) ]
            if {[string index $deslashedFileName 0] != "/"} {
		#remove ./ if name only typed in by tailing value
		set deslashedFileName [file tail $deslashedFileName ]
		set deslashedFileName [pwd]/$deslashedFileName
		set probtrack(output) [rmSlash $deslashedFileName]
	    }

       	    if { $canwrite } {
       		puts "rm -rf $probtrack(output)"
       		exec rm -rf $probtrack(output)
	        puts "mkdir -p $probtrack(output)"
		exec mkdir -p $probtrack(output)
       	    }

	    set filebase $probtrack(output)/fdt
	    set logfile "${filebase}_log.tcl"
	    set log [open "$logfile" w]
	    puts $log "set tool $probtrack(tool)"
	    set copylog ""

	    if { $probtrack(usereference_yn) } {
		set flags "$flags --xfm=$probtrack(xfm)"
      		puts $log "set probtrack(usereference_yn) $probtrack(usereference_yn)"
      		puts $log "set probtrack(xfm) $probtrack(xfm)"
		if { $probtrack(useNonlinear) } { 
		    set flags "$flags --invxfm=$probtrack(invxfm)" 
		    puts $log "set $probtrack(useNonlinear) $probtrack(useNonlinear)"
		    puts $log "set probtrack(invxfm) $probtrack(invxfm)"
		}
		if { $probtrack(useSurface) } {
		    set flags "$flags --meshspace=$probtrack(meshspace) --seedref=$probtrack(surfref)"
		    puts $log "set probtrack(useSurface) $probtrack(useSurface)"
		    puts $log "set probtrack(meshspace) $probtrack(meshspace)"
		    puts $log "set probtrack(surfref) $probtrack(surfref)"
		}
	    }	    

	    if { $probtrack(exclude_yn) == 1 } {
		set flags "$flags --avoid=$probtrack(exclude)"
		puts $log "set probtrack(exclude_yn) $probtrack(exclude_yn)"
		puts $log "set probtrack(exclude) $probtrack(exclude)"
	    }

	    if { $probtrack(terminate_yn) == 1 } {
		set flags "$flags --stop=$probtrack(stop)"
		puts $log "set probtrack(terminate_yn) $probtrack(terminate_yn)"
		puts $log "set probtrack(stop) $probtrack(stop)"
	    }

	    set flags "$flags --forcedir --opd -s $probtrack(bedpost_dir)/merged -m $probtrack(bedpost_dir)/nodif_brain_mask  --dir=$probtrack(output)" 
    	    foreach entry {bedpost_dir xfm mode exclude_yn usereference_yn verbose_yn loopcheck_yn modeuler_yn curvature nsteps steplength nparticles fibthresh distthresh sampvox oneway_yn wayorder_yn waycond} {
		puts $log "set probtrack($entry) $probtrack($entry)"
	    }
		    set singleFileName $probtrack(output)
            switch $probtrack(mode) {

	       simple { 
		    set singleFileName [ file tail $probtrack(output) ]
		    set fd [ open "${filebase}_coordinates.txt" w ]
		    set x $probtrack(x)
		    set y $probtrack(y)
		    set z $probtrack(z)
		   if { ! $probtrack(usereference_yn) } {
                       set probtrack(reference) [ file join $probtrack(bedpost_dir) nodif_brain_mask ]
		   }
		    if { $probtrack(units) == "mm" } {
			if { $probtrack(reference) != "" } {
			    mm_to_voxels x y z $probtrack(reference)
			} else {
			    mm_to_voxels x y z [ file join $probtrack(bedpost_dir) nodif_brain_mask ]
			}			    
			puts $fd "$x $y $z"
			puts "$probtrack(x) $probtrack(y) $probtrack(z) (mm) -> $x $y $z (voxels)"
		    } else {
			puts $fd "$probtrack(x) $probtrack(y) $probtrack(z)"
		    }
		    close $fd
 		    puts $log "set probtrack(x) $probtrack(x)"
		    puts $log "set probtrack(y) $probtrack(y)"
		    puts $log "set probtrack(z) $probtrack(z)"
		    puts $log "set probtrack(units) $probtrack(units)"
		   set flags "--simple --seedref=$probtrack(reference) -o ${singleFileName} -x ${filebase}_coordinates.txt $flags"
	       } 
               seedmask {
                   set flags " -x $probtrack(reference) $flags"  
	       }
	       network {
                   fdt_exp w $w.data.seed.targets $probtrack(output)/masks.txt
		   set flags "--network -x $probtrack(output)/masks.txt $flags"
		   puts $log  " $w.data.seed.targets insert end [  $w.data.seed.targets get 0 end ]"
	       }
	    }
	    puts $log "set probtrack(reference) $probtrack(reference)"
	    puts $log "set probtrack(output) $probtrack(output)"
       	    if { $canwrite } {
       		set copylog "$probtrack(output)/fdt.log"
	        if { $probtrack(waypoint_yn) == 1 } {
                    fdt_exp w $w.data.targets.wf.tf.targets $probtrack(output)/waypoints.txt
		    puts $log "set probtrack(waypoint_yn) $probtrack(waypoint_yn)"
                    puts $log " $w.data.targets.wf.tf.targets insert end [  $w.data.targets.wf.tf.targets get 0 end ]"
                    set flags "$flags --waypoints=$probtrack(output)/waypoints.txt "
		    set flags "$flags --waycond=$probtrack(waycond)"
	        } 
	        if { $probtrack(classify_yn) == 1 } {
                    fdt_exp w $w.data.targets.cf.tf.targets $probtrack(output)/targets.txt
		    puts $log "set probtrack(classify_yn) $probtrack(classify_yn)"
		    puts $log " $w.data.targets.cf.tf.targets insert end [  $w.data.targets.cf.tf.targets get 0 end ]"
                    set flags "$flags --targetmasks=$probtrack(output)/targets.txt --os2t "
                }
		close $log
		if { $FSLPARALLEL } {
                    set script [open "${filebase}_script.sh" w]
                    puts "${filebase}_script.sh"
                    exec chmod 777 ${filebase}_script.sh
                    puts $script "#!/bin/sh"
                    puts $script "cd $probtrack(output)"
                    puts $script "$FSLDIR/bin/probtrackx2 $flags"		    
                    if { $probtrack(classify_yn) == 1 } {
			puts $script "$FSLDIR/bin/find_the_biggest seeds_to_* $probtrack(output)/biggest >> fdt_seed_classification.txt"
		    }
                    #if { $probtrack(mode) == "simple" } {
                    #puts $script "rm ${filebase}_coordinates.txt"
		    #}
                    puts $script "mv $logfile $copylog"
                    puts $script "rm ${filebase}_script.sh"
		    close $script
		    exec $FSLDIR/bin/fsl_sub -q long.q ${filebase}_script.sh
		} else {

		    fdt_monitor_short $w "$FSLDIR/bin/probtrackx2 $flags"
		    if { $probtrack(classify_yn) == 1 } {
			exec sh -c "$FSLDIR/bin/find_the_biggest $probtrack(output)/seeds_to_* $probtrack(output)/biggest >> $probtrack(output)/fdt_seed_classification.txt"
		    }
		}
       	    }
            if { !$FSLPARALLEL } {
		#if { $probtrack(mode) == "simple" } {
		    #puts "rm ${filebase}_coordinates.txt"
		    #exec rm ${filebase}_coordinates.txt
		#}
		if { $copylog != "" } {
		    puts "mv $logfile $copylog"
		    exec mv $logfile $copylog
		} else {
		    puts "rm $logfile"
		    exec rm $logfile
		}
	    }
	}
	registration {
	    global registration

	    set errorStr ""
	    if { $registration(directory) == ""  } { set errorStr "You must specify the bedpostX directory!" }
	    if { $registration(struct_yn) && $registration(struct_image) == ""  } { set errorStr "$errorStr You must specify the structural image!" }
	    if { $registration(struct_yn) && $registration(nonlin_yn) && $registration(nonlin_image) == ""  } { set errorStr "$errorStr You must specify the non-betted structural image!" }
	    if { $registration(standard_yn) && $registration(standard_image) == ""  } { set errorStr "$errorStr You must specify the standard image!" }
	    if { $errorStr != "" } {
		MxPause $errorStr
		return
	    }

	    exec mkdir -p [ file join $registration(directory) xfms ]
	    set eyefd [ open [ file join $registration(directory) xfms eye.mat ] w ]
	    puts $eyefd "1 0 0 0"
	    puts $eyefd "0 1 0 0"
	    puts $eyefd "0 0 1 0"
	    puts $eyefd "0 0 0 1"
	    close $eyefd

	    set diff2str   [ file join $registration(directory) xfms diff2str.mat ]
	    set str2diff   [ file join $registration(directory) xfms str2diff.mat ]
	    set str2stand  [ file join $registration(directory) xfms str2standard.mat ]
	    set str2stand_warp  [ file join $registration(directory) xfms str2standard_warp ]
	    set stand2str  [ file join $registration(directory) xfms standard2str.mat ]
	    set stand2str_warp  [ file join $registration(directory) xfms standard2str_warp ]
	    set diff2stand [ file join $registration(directory) xfms diff2standard.mat ]
	    set diff2stand_warp [ file join $registration(directory) xfms diff2standard_warp ]
	    set stand2diff [ file join $registration(directory) xfms standard2diff.mat ]
	    set stand2diff_warp [ file join $registration(directory) xfms standard2diff_warp ]
	    set diff       [ file join $registration(directory) nodif_brain ]
	    if { $registration(struct_yn) } {
		set searchrx  "-searchrx -$registration(struct_search) $registration(struct_search)"
		set searchry  "-searchry -$registration(struct_search) $registration(struct_search)"
		set searchrz  "-searchrz -$registration(struct_search) $registration(struct_search)"
		set options   "$searchrx $searchry $searchrz -dof $registration(struct_dof)"
		fdt_monitor $w "${FSLDIR}/bin/flirt -in $diff -ref $registration(struct_image) -omat $diff2str $options -cost $registration(struct_costfn)"
		fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $str2diff -inverse $diff2str"
		if { $registration(standard_yn) } {
		    set searchrx  "-searchrx -$registration(standard_search) $registration(standard_search)"
		    set searchry  "-searchry -$registration(standard_search) $registration(standard_search)"
		    set searchrz  "-searchrz -$registration(standard_search) $registration(standard_search)"
		    set options   "$searchrx $searchry $searchrz -dof $registration(standard_dof)"
		    fdt_monitor $w "${FSLDIR}/bin/flirt -in $registration(struct_image) -ref $registration(standard_image) -omat $str2stand $options -cost $registration(standard_costfn)"
		    fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $stand2str -inverse $str2stand"
		    fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $diff2stand -concat $str2stand $diff2str"
		    fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $stand2diff -inverse $diff2stand"
		    if { $registration(nonlin_yn) } {
			fdt_monitor $w "${FSLDIR}/bin/fnirt --in=$registration(nonlin_image) --aff=$str2stand --cout=$str2stand_warp --config=T1_2_MNI152_2mm"
			fdt_monitor $w "${FSLDIR}/bin/invwarp -w $str2stand_warp -o $stand2str_warp -r $registration(struct_image)"
			fdt_monitor $w "${FSLDIR}/bin/convertwarp -o $diff2stand_warp -r ${FSLDIR}/data/standard/MNI152_T1_2mm -m $diff2str -w $str2stand_warp"
			fdt_monitor $w "${FSLDIR}/bin/convertwarp -o $stand2diff_warp -r ${diff}_mask -w $stand2str_warp --postmat=$str2diff"
		    }
		}
	    } elseif { $registration(standard_yn) } {
		set searchrx  "-searchrx -$registration(standard_search) $registration(standard_search)"
		set searchry  "-searchry -$registration(standard_search) $registration(standard_search)"
		set searchrz  "-searchrz -$registration(standard_search) $registration(standard_search)"
		set options   "$searchrx $searchry $searchrz -dof $registration(standard_dof)"
		fdt_monitor $w "${FSLDIR}/bin/flirt -in $diff -ref $registration(standard_image) -omat $diff2stand $options"
		fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $stand2diff -inverse $diff2stand"
	    }
	    puts "Done!"
	    # Fudge to make the logic work
	    set canwrite 1
	}
    }

    if { $canwrite } { 
	if { $FSLPARALLEL == 1  && [ string match *x* $probtrack(tool) ]} { MxPause " Job submitted to queue" } else  { MxPause "  Done!  " }
	update idletasks
    }

    if {$dialog == "destroy"} {
        fdt:destroy $w
    }
}


proc fdt:destroy { w } {
    destroy $w
}    

set debugging 0

while {[llength $argv] > 0 } {
    set flag [lindex $argv 0]
    switch -- $flag {
	"-debugging" {
	    set debugging 1
	    set argv [lrange $argv 1 end]
	    puts "Debug mode!"
	}
	default { break }
    }
}


wm withdraw .
if { [ info exists env(MRDATADIR) ] } {
    set MRDATADIR $env(MRDATADIR)
} else {
    set MRDATADIR ~/MRdata
}

fdt:dialog .fdt $argv
tkwait window .fdt
