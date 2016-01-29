# 
# ASL
# Michael Chappell  and Matthew Webster, FMRIB Image Analysis Group
#
# Copyright (C) 2010-2013 University of Oxford
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

set MYEXEC    [ info nameofexecutable ]
set MYSHELL   [ file tail $MYEXEC ]
if { [  string match -nocase *wish* $MYSHELL ] } {
option add *LabelEntry.e.background grey95
}

proc asl { w } {
    global FSLDIR Asl
    # ---- Set up Frames ----
    toplevel $w
    wm title $w "ASL"
    wm iconname $w "ASL"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
    NoteBook $w.nb -side top -bd 2 -tabpady {5 10} -arcradius 3 
    set Asl(NB) $w.nb
    $w.nb insert 0 data         -text "Data" 
    $w.nb insert 1 analysis     -text "Analysis"    
    $w.nb insert 2 registration -text "Registration"     
    $w.nb insert 3 calibration  -text "Calibration"  

    set Asl(data) [ $w.nb getframe data ]
    FileEntry $Asl(data).input -textvariable Asl(input) -label "Input Filename  " -title "Select the input image" -width 35 -filedialog directory -filetypes IMAGE
    LabelEntry $Asl(data).inversionTimes -label "Inversion Times" -textvariable Asl(inversionTimes) -width 35
    set Asl(bolusDuration) 1
    LabelSpinBox $Asl(data).bolus -label "Bolus duration "  -textvariable Asl(bolusDuration) -range { 0.1 10.0 .1 } -width 5 
    frame $Asl(data).lf
    label $Asl(data).lf.label -text "Labelling: "
    set Asl(labeling) 0
    optionMenu2 $Asl(data).lf.labeling Asl(labeling) -command "set Asl(inversionEfficiency) 0.98; if { \$Asl(labeling) == 1 } { set Asl(inversionEfficiency) 0.85 }" 0 "pASL" 1 "cASl/pcASL" 

    frame $Asl(data).tf
    label $Asl(data).tf.label -text "Data is tag-control pairs: "
    set Asl(tagControl) 1
    checkbutton $Asl(data).tf.button -variable Asl(tagControl) -anchor w -command "updateAslMode"

    frame $Asl(data).tff
    label $Asl(data).tff.label -text "Data has control first: "
    set Asl(controlFirst) 0
    checkbutton $Asl(data).tff.button -variable Asl(controlFirst) -anchor w 

    labelframe $Asl(data).useStructural -text "Structural image" -labelanchor w -relief flat
    checkbutton $Asl(data).useStructural.button -variable Asl(useStructural) -command "updateRegistration"
    FileEntry $Asl(data).useStructural.file -textvariable Asl(structural) -title "Select the structural image" -width 35 -filedialog directory -filetypes IMAGE

    set Asl(betStructural) 1
    labelframe $Asl(data).betStructural -text "Run BET on Structural image" -labelanchor w -relief flat
    checkbutton $Asl(data).betStructural.button -variable Asl(betStructural) -command "updateRegistration"

    frame $Asl(data).df
    label $Asl(data).df.label -text "Data order (grouped by): "
    optionMenu2 $Asl(data).df.order Asl(dataOrder)  0 "repeats" 1 "TIs"

    frame $Asl(data).sf
    set Asl(staticTissue) "normal"
    label $Asl(data).sf.label -text "Static tissue: "
    optionMenu2 $Asl(data).sf.type Asl(staticTissue) -command "updateAslMode" "normal" "normal" "suppressed" "background suppressed" "preSaturation" "pre saturation"

    pack $Asl(data).input -anchor w
    pack $Asl(data).inversionTimes -anchor w
    pack $Asl(data).bolus -anchor w
    pack $Asl(data).lf -anchor w
    pack $Asl(data).lf.label $Asl(data).lf.labeling -side left
    pack $Asl(data).tf -anchor w
    pack $Asl(data).tf.label $Asl(data).tf.button -side left
    pack $Asl(data).tff -anchor w
    pack $Asl(data).tff.label $Asl(data).tff.button -side left
    pack $Asl(data).df -anchor w
    pack $Asl(data).df.label $Asl(data).df.order -side left
    pack $Asl(data).sf -anchor w
    pack $Asl(data).sf.label $Asl(data).sf.type -side left
    pack $Asl(data).useStructural -anchor w
    pack $Asl(data).useStructural.button  -side left
    pack $Asl(data).betStructural.button -anchor w

    set Asl(analysis) [ $w.nb getframe analysis ]
    FileEntry $Asl(analysis).outputdir -textvariable Asl(outputdir) -label "Output directory" -title "Name the output directory" -width 35 -filedialog directory -filetypes { }
    FileEntry $Asl(analysis).brainMask -textvariable Asl(brainMask) -label "Optional Brain Mask" -title "Select Brain Mask" -width 35 -filedialog directory -filetypes IMAGE

    labelframe $Asl(analysis).outputVariance -text "Output parameter variance" -labelanchor w -relief flat
    checkbutton $Asl(analysis).outputVariance.button -variable Asl(outputVariance)
    set Asl(bolusArrival) 0.7
    LabelSpinBox $Asl(analysis).bolusArrival -label "Bolus arrival time  "  -textvariable Asl(bolusArrival) -range {0.1 100.0 .1 } -width 5 
    set Asl(t1) 1.3
    LabelSpinBox $Asl(analysis).t1 -label "T1   "  -textvariable Asl(t1) -range {0.1 10.0 .1 } -width 5 
    set Asl(t1b) 1.6
    LabelSpinBox $Asl(analysis).t1b -label "T1b "  -textvariable Asl(t1b) -range {0.1 10.0 .1 } -width 5
    set Asl(inversionEfficiency) 0.98
    LabelSpinBox $Asl(analysis).invE -label "Inversion efficiency "  -textvariable Asl(inversionEfficiency) -range {0.0 1.0 .01 } -width 5
    labelframe $Asl(analysis).adaptiveSmoothing -text "Use adaptive spatial smoothing on CBF" -labelanchor w -relief flat
    checkbutton $Asl(analysis).adaptiveSmoothing.button -variable Asl(adaptiveSmoothing) 
    labelframe $Asl(analysis).useUncertainty -text "Incorporate T1 value uncertainty" -labelanchor w -relief flat
    checkbutton $Asl(analysis).useUncertainty.button -variable Asl(useUncertainty) 
    labelframe $Asl(analysis).flowSuppression -text "Include macro vascular component" -labelanchor w -relief flat
    checkbutton $Asl(analysis).flowSuppression.button -variable Asl(flowSuppression)
    set Asl(fixBolus) 1
    labelframe $Asl(analysis).fixBolus -text "Fix bolus duration" -labelanchor w -relief flat
    checkbutton $Asl(analysis).fixBolus.button -variable Asl(fixBolus) 
    pack $Asl(analysis).outputdir -anchor w
    pack $Asl(analysis).brainMask -anchor w
    pack $Asl(analysis).outputVariance.button
    pack $Asl(analysis).outputVariance -anchor w
    pack $Asl(analysis).bolusArrival -anchor w
    pack $Asl(analysis).t1 -anchor w
    pack $Asl(analysis).t1b -anchor w
    pack $Asl(analysis).invE -anchor w
    pack $Asl(analysis).adaptiveSmoothing.button
    pack $Asl(analysis).adaptiveSmoothing -anchor w
    pack $Asl(analysis).useUncertainty.button
    pack $Asl(analysis).useUncertainty -anchor w
    pack $Asl(analysis).flowSuppression.button
    pack $Asl(analysis).flowSuppression -anchor w
    pack $Asl(analysis).fixBolus.button
    pack $Asl(analysis).fixBolus -anchor w

    set Asl(registration) [ $w.nb getframe registration ]

    label $Asl(registration).needStructural -text "Registration requires a structural image input ( via Data tab )" -font {Helvetica 12 bold}

    labelframe $Asl(registration).useTransform -text "Structural to standard space transform" -labelanchor w -relief flat
    checkbutton $Asl(registration).useTransform.button -variable Asl(useTransform) -command "pack forget $Asl(registration).useTransform.file; if { \$Asl(useTransform) } { pack $Asl(registration).useTransform.file -side left -anchor w}; $w.nb compute_size"
    FileEntry $Asl(registration).useTransform.file -textvariable Asl(transform) -title "Select the structural to standard space transformation" -width 35 -filedialog directory -filetypes {{*.mat,*.nii.gz,*.nii,*.img,*.img.gz}}

    labelframe $Asl(registration).useStandard -text "Alternate standard brain image" -labelanchor w -relief flat
    checkbutton $Asl(registration).useStandard.button -variable Asl(useStandard) -command "pack forget $Asl(registration).useStandard.file; if { \$Asl(useStandard) } { pack $Asl(registration).useStandard.file -side left -anchor w}; $w.nb compute_size "
    FileEntry $Asl(registration).useStandard.file -textvariable Asl(standard) -title "Select the alternate standard image" -width 35 -filedialog directory -filetypes IMAGE

    pack $Asl(registration).needStructural
    pack $Asl(registration).useTransform.button -side left
    pack $Asl(registration).useStandard.button -side left

    set Asl(calibration) [ $w.nb getframe calibration ]
    set Asl(doCalibration) 0
    labelframe  $Asl(calibration).useCalibration -text "Perform calibration"  -labelanchor w -relief flat
    checkbutton $Asl(calibration).useCalibration.button -variable Asl(doCalibration) -command "updateCalibration"
    labelframe $Asl(calibration).reference -text "Reference Tissue" -labelanchor n
    labelframe $Asl(calibration).reference.type -text "Reference Tissue Type: " -labelanchor w -relief flat
    labelframe  $Asl(calibration).reference.mask -text "Reference Tissue Mask" -labelanchor w -relief flat
    checkbutton $Asl(calibration).reference.mask.button -variable Asl(userefmask) -disabledforeground "yellow" -command "pack forget $Asl(calibration).reference.mask.file; if { \$Asl(userefmask) } { pack $Asl(calibration).reference.mask.file -side left -anchor w; $w.nb compute_size}"
    FileEntry $Asl(calibration).reference.mask.file -textvariable Asl(mask) -title "Select the mask image" -width 35 -filedialog directory -filetypes IMAGE
    pack  $Asl(calibration).reference.mask.button -side left -anchor w

    set Asl(tissueType) "csf"
    optionMenu2 $Asl(calibration).reference.type.control Asl(tissueType) -command "updateReference" "csf" "CSF" "wm" "White matter" "gm" "Grey matter" "none" "none"
    set Asl(ReferenceT1) 4.3
    LabelSpinBox $Asl(calibration).reference.t1 -textvariable Asl(ReferenceT1) -label "Reference T1(s)" -range {0.0 10.0 0.1 } 
    set Asl(ReferenceT2) 0.75
    LabelSpinBox $Asl(calibration).reference.t2 -textvariable Asl(ReferenceT2) -label "Reference T2(s)" -range {0.0 10.0 0.1 } 
    set Asl(BloodT2) 0.15
    LabelSpinBox $Asl(calibration).reference.t2b -textvariable Asl(BloodT2) -label "Blood T2(s)" -range {0.0 10.0 0.1 } 


    labelframe $Asl(calibration).sequence -text "Sequence parameters" -labelanchor n
    set Asl(SequenceTR) 3.2
    LabelSpinBox $Asl(calibration).sequence.tr -textvariable Asl(SequenceTR) -label "Sequence TR(s)" -range {0.0 30.0 0.1 } 
    set Asl(SequenceTE) 0.0
    LabelSpinBox $Asl(calibration).sequence.te -textvariable Asl(SequenceTE) -label "Sequence TE(s)" -range {0.0 10.0 0.1 } 


    frame $Asl(calibration).mode 
    frame $Asl(calibration).mode.mf
    label $Asl(calibration).mode.mf.label1 -text "Mode: "
    set Asl(Mode) 0
    optionMenu2 $Asl(calibration).mode.mf.control Asl(Mode) -command "updateAslMode" 0 "Long TR" 1 "Saturation Recovery"
    $Asl(calibration).mode.mf.control.menu entryconfigure 1 -state disabled	

    FileEntry $Asl(calibration).mode.structural -textvariable Asl(M0image) -label "M0 calibration image" -title "Select the calibration image" -width 35 -filedialog directory -filetypes IMAGE 
    labelframe $Asl(calibration).mode.useCoil -text "Use Coil sensitivity reference image" -labelanchor w -relief flat
    checkbutton $Asl(calibration).mode.useCoil.button -variable Asl(useCoil) -command "pack forget $Asl(calibration).mode.useCoil.file; if { \$Asl(useCoil) } { pack $Asl(calibration).mode.useCoil.file -side left -anchor w; $w.nb compute_size}"
    FileEntry $Asl(calibration).mode.useCoil.file -textvariable Asl(coilImage) -title "Select the mask image" -width 35 -filedialog directory -filetypes IMAGE

    set Asl(calibrationGain) 1.0
    LabelSpinBox $Asl(calibration).mode.gain -textvariable Asl(calibrationGain) -label "Calibration Gain" -range {0.0 30.0 0.1 } 

    pack $Asl(calibration).useCalibration.button -side left
    pack $Asl(calibration).useCalibration -anchor w
    #pack $Asl(calibration).mode -anchor w
    pack $Asl(calibration).mode.mf -anchor w
    pack $Asl(calibration).mode.mf.label1 $Asl(calibration).mode.mf.control -side left -anchor w
    pack $Asl(calibration).mode.useCoil.button -side left
    pack $Asl(calibration).mode.useCoil -anchor w
    pack $Asl(calibration).mode.gain -side left
    #pack $Asl(calibration).reference -anchor w
    pack $Asl(calibration).reference.type $Asl(calibration).reference.type.control -anchor w
    pack $Asl(calibration).reference.mask -anchor w
    pack $Asl(calibration).reference.t1 $Asl(calibration).reference.t2 $Asl(calibration).reference.t2b  -anchor w
    #pack $Asl(calibration).sequence -anchor w
    pack $Asl(calibration).sequence.tr $Asl(calibration).sequence.te -anchor w

    pack $w.nb
    $w.nb raise data
								    
    button $w.execute -command "aslLaunch" -text "Go"
    button $w.cancel -command "destroy $w" -text "Exit" 
    pack $w.execute $w.cancel -side left  -padx 120 -anchor center								    
}

proc updateCalibration { } {
    global FSLDIR Asl FSLDEVDIR 
    pack forget $Asl(calibration).mode $Asl(calibration).reference $Asl(calibration).sequence
    $Asl(calibration).reference.mask.button configure -state active

    if { $Asl(doCalibration) } { 
	pack  $Asl(calibration).mode $Asl(calibration).reference $Asl(calibration).sequence -anchor w
    }
    if { ! $Asl(useStructural) } {
	if { !$Asl(userefmask)} {
	    $Asl(calibration).reference.mask.button invoke
	}
        $Asl(calibration).reference.mask.button configure -state disabled
	pack $Asl(calibration).reference.mask.file -side left -anchor w
    }
    $Asl(NB) compute_size
}

proc updateRegistration { } {
    global FSLDIR Asl FSLDEVDIR 
    pack forget $Asl(data).useStructural.file $Asl(data).betStructural $Asl(registration).needStructural $Asl(registration).useTransform $Asl(registration).useStandard

    if { $Asl(useStructural) } {
	pack $Asl(data).useStructural.file -side left -anchor w
	pack $Asl(data).betStructural -side left -anchor w
	pack $Asl(registration).useTransform $Asl(registration).useStandard -anchor w
    } else {
	pack $Asl(registration).needStructural
    }
    updateCalibration
}

proc updateReference { } {
    global FSLDIR Asl FSLDEVDIR 
    pack forget  $Asl(calibration).reference.mask.file $Asl(calibration).reference.mask.button
    if { $Asl(tissueType) == "none" } {
	set Asl(userefmask) 1
    } else {
	pack  $Asl(calibration).reference.mask.button -side left -anchor w

    }
    if { $Asl(userefmask) } { pack $Asl(calibration).reference.mask.file -side left -anchor w}
}

proc updateAslMode { } {
    global FSLDIR Asl FSLDEVDIR 
    $Asl(calibration).mode.mf.control.menu entryconfigure 0 -state normal
    $Asl(calibration).mode.mf.control.menu entryconfigure 1 -state normal
    pack forget optionMenu2 $Asl(data).sf $Asl(data).useStructural $Asl(data).betStructural

    if  {  $Asl(tagControl) == 1 } { 
	pack $Asl(data).sf -anchor w
    }
    pack $Asl(data).useStructural -anchor w
    pack $Asl(data).betStructural -anchor w

    if {  $Asl(tagControl) == 1 && $Asl(staticTissue) == "normal" } { 
	set Asl(Mode) 0
	$Asl(calibration).mode.mf.control.menu entryconfigure 1 -state disabled	
    }
    if {  $Asl(tagControl) == 1 && $Asl(staticTissue) == "preSaturation" } { 
	set Asl(Mode) 1
	$Asl(calibration).mode.mf.control.menu entryconfigure 0 -state disabled	
    }
    pack forget $Asl(calibration).mode.structural $Asl(calibration).mode.useCoil $Asl(calibration).mode.gain
    if { $Asl(tagControl) == 1 && $Asl(staticTissue) == "suppressed" } {
	pack  $Asl(calibration).mode.structural -anchor w
    }	
    if { $Asl(tagControl) == 0 } {
	pack  $Asl(calibration).mode.structural -anchor w
    }	

    pack $Asl(calibration).mode.useCoil  $Asl(calibration).mode.gain -anchor w
    $Asl(NB) compute_size
}

proc aslLaunch {  } {
    global FSLDIR Asl FSLDEVDIR 


    #  if { ![ file exists $filename ] } {
    #     MxPause "Warning: Bad or missing file!"
    #      return
    #  }
    if { $Asl(input)=="" } {
	MxPause "You have not specified an input file!"
	return
    }
    if { $Asl(inversionTimes)=="" } {
	MxPause "You have not specified any inversion times!"
	return
    }

    fsl:exec "mkdir $Asl(outputdir)"
    fsl:exec "mkdir $Asl(outputdir)/native_space"

    set tisListLength [ llength [ split $Asl(inversionTimes) , ] ]
    puts "$tisListLength"
    #todo check ntis with MC
    set aslFileCommand "$FSLDIR/bin/asl_file --data=$Asl(input) --out=$Asl(outputdir)/native_space/diffData --obf=rpt --ntis=$tisListLength" 
    if { $Asl(dataOrder) == 0 } { 
	set aslFileCommand "$aslFileCommand --ibf=rpt" 
    } else {
	set aslFileCommand "$aslFileCommand --ibf=tis" 
    }

    set aslDataForm "--iaf=tc" 
    if { $Asl(controlFirst) == 1 } {
	set aslDataForm "--iaf=ct"
    }


    if { $Asl(tagControl) == 1 } {  #0 "none" 1 "pairwise"
	set aslFileCommand "$aslFileCommand $aslDataForm --diff" 
    } else {
	set aslFileCommand "$aslFileCommand --iaf=diff" 
    }

    puts $aslFileCommand
    fsl:exec "$aslFileCommand"

    set aslCommand "$FSLDIR/bin/oxford_asl -i $Asl(outputdir)/native_space/diffData -o $Asl(outputdir) --tis $Asl(inversionTimes)"
    if { $Asl(outputVariance) } {
	set aslCommand "$aslCommand --vars "
    }
    set aslCommand "$aslCommand --bolus $Asl(bolusDuration) --bat $Asl(bolusArrival) --t1 $Asl(t1) --t1b $Asl(t1b) --alpha $Asl(inversionEfficiency)"
    if { $Asl(adaptiveSmoothing) } {
	set aslCommand "$aslCommand --spatial "
    }
    if { $Asl(useUncertainty) } {
	set aslCommand "$aslCommand --infert1 "
    }
    if { !$Asl(flowSuppression) } {
	set aslCommand "$aslCommand --artoff "
    }
    if { $Asl(fixBolus) } {
	set aslCommand "$aslCommand --fixbolus "
    }
    if { $Asl(brainMask) != "" } {
	set aslCommand "$aslCommand -m $Asl(brainMask) "
    }
    if { $Asl(labeling) } {
	set aslCommand "$aslCommand  --casl "
    }
    if { $Asl(useStructural) } {
	if { $Asl(betStructural) } {
	    fsl:exec "${FSLDIR}/bin/bet $Asl(structural) $Asl(outputdir)/structural_brain"
	} else {
	    fsl:exec "$FSLDIR/bin/imcp $Asl(structural) $Asl(outputdir)/structural_brain"
	}

	set aslCommand "$aslCommand -s $Asl(outputdir)/structural_brain"
	if { $Asl(useStandard) } {
	    set aslCommand "$aslCommand  -t $Asl(transform) -S $Asl(standard) "
	}
    }



    set foundRegTarget 0

    if { $Asl(Mode) == 1 || ( $Asl(tagControl) == 1 && $Asl(staticTissue) != "suppressed" ) } {
	fsl:exec "${FSLDIR}/bin/asl_file --data=$Asl(input) --ntis=$tisListLength $aslDataForm --spairs --out=$Asl(outputdir)/asldata_mc"
	if { $Asl(Mode) == 1 } {
	    fsl:exec "${FSLDIR}/bin/bet $Asl(outputdir)/asldata_mc_odd $Asl(outputdir)/asldata_mc_odd_brain"
	    set aslCommand "$aslCommand --regfrom $Asl(outputdir)/asldata_mc_odd_brain "
	} else {
	    fsl:exec "${FSLDIR}/bin/bet $Asl(outputdir)/asldata_mc_even $Asl(outputdir)/asldata_mc_even_brain"
	    set aslCommand "$aslCommand --regfrom $Asl(outputdir)/asldata_mc_even_brain "
	}
	set foundRegTarget 1
    }

    if { $Asl(doCalibration) && !$foundRegTarget } {
	fsl:exec "${FSLDIR}/bin/bet $Asl(M0image) $Asl(outputdir)/MOimage_brain"
	set aslCommand "$aslCommand --regfrom  $Asl(outputdir)/MOimage_brain"
    }

    puts $aslCommand
    fsl:exec "$aslCommand"
    if { $Asl(doCalibration) } {
	set aslCalibCommand "$FSLDIR/bin/asl_calib -i $Asl(outputdir)/native_space/perfusion --tissref $Asl(tissueType) --t1r $Asl(ReferenceT1) --t2r $Asl(ReferenceT2) --t2b $Asl(BloodT2) --te $Asl(SequenceTE) -o $Asl(outputdir)/calibration"

	if { $Asl(userefmask) } {
	    set aslCalibCommand "$aslCalibCommand -m $Asl(mask)"
	} else {
	    set aslCalibCommand "$aslCalibCommand -s $Asl(outputdir)/structural_brain -t $Asl(outputdir)/native_space/asl2struct.mat"
	}

	if { $Asl(brainMask) != "" } {
	    set aslCalibCommand "$aslCalibCommand --bmask $Asl(brainMask)"
	}
    
	if { $Asl(Mode) == 0 } {  #long TR
	    set calibImage $Asl(M0image)
	    if  { $Asl(tagControl) == 1 && $Asl(staticTissue) == "normal" } {
		set calibImage $Asl(outputdir)/asldata_mc_even
	    }
	    set aslCalibCommand "$aslCalibCommand -c $calibImage --mode longtr --tr $Asl(SequenceTR) --cgain $Asl(calibrationGain)" 
	    if { $Asl(useCoil) } {
		set aslCalibCommand "$aslCalibCommand --cref $Asl(coilImage)"
	    }
	}
	if { $Asl(Mode) == 1 } { #Saturation Recovery - maybe change aslcalib -c option to even - for MC to think about!!!
	    set aslCalibCommand "$aslCalibCommand --tis $Asl(inversionTimes) -c $Asl(outputdir)/asldata_mc_odd"
	}
	puts $aslCalibCommand
	fsl:exec "$aslCalibCommand"
    }
    set aslDivCommand ""
    if { $Asl(useCoil) } {
	set aslDivCommand "-div $Asl(coilImage)"
    }	 

    if { [ imtest $Asl(outputdir)/perfusion ] } {
	set M0 [ fsl:exec "cat $Asl(outputdir)/calibration/M0.txt" ] 
	fsl:exec "${FSLDIR}/bin/fslmaths $Asl(outputdir)/perfusion -div $M0 $aslDivCommand -mul 6000 $Asl(outputdir)/perfusion_calib"
    }
    if { [ imtest $Asl(outputdir)/standard_space/perfusion ] } {
	set M0 [ fsl:exec "cat $Asl(outputdir)/calibration/M0.txt" ] 
	fsl:exec "${FSLDIR}/bin/fslmaths $Asl(outputdir)/standard_space/perfusion -div $M0 $aslDivCommand -mul 6000 $Asl(outputdir)/standard_space/perfusion_calib"
    }
    if { [ imtest $Asl(outputdir)/structural_space/perfusion ] } {
	set M0 [ fsl:exec "cat $Asl(outputdir)/calibration/M0.txt" ] 
	fsl:exec "${FSLDIR}/bin/fslmaths $Asl(outputdir)/structural_space/perfusion -div $M0 $aslDivCommand -mul 6000 $Asl(outputdir)/structural_space/perfusion_calib"
    }


    if { $Asl(useStandard) == 1 } {
	set regOptions "--prefix=$Asl(outputdir)/native_space/asl2struct.mat -r $Asl(standard) "
	if { [ imtest $Asl(transform) ] } {
	    set regOptions "$regOptions -w $Asl(transform) "
	} else {
	    set regOptions "$regOptions --postmat=$Asl(transform) "
	}
	if { [ imtest $Asl(outputdir)/native_space/perfusion_calib ] } {
	    fsl:exec "${FSLDIR}/bin/applywarp $regOptions -i $Asl(outputdir)/native_space/perfusion_calib -o $Asl(outputdir)/native_space/perfusion_calib_standard"
	} else {
	    fsl:exec "${FSLDIR}/bin/applywarp $regOptions -i $Asl(outputdir)/native_space/perfusion -o $Asl(outputdir)/native_space/perfusion_standard"
	}
	fsl:exec "${FSLDIR}/bin/applywarp $regOptions -i $Asl(outputdir)/native_space/arrival -o $Asl(outputdir)/native_space/arrival_standard"
    }
}

if { [  string match -nocase *wish* $MYSHELL ] } {
    wm withdraw .
    asl .rename
    tkwait window .rename
}
