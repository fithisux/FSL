#

# PNM - the GUI for calculating a Physiological Noise Model suitable for FEAT
#
# Mark Jenkinson, FMRIB Image Analysis Group
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

#{{{ pnm

proc pnm { w } {

    global entries FSLDIR PWD

    # ---- Set up Frames ----
    toplevel $w
    wm title $w "PNM"
    wm iconname $w "PNM"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
    frame $w.f

    TitleFrame $w.f.input -text "Input"  -relief groove 
    set lfinput [ $w.f.input getframe ]

    # Input physio data

    set entries($w,intrace) ""

    FileEntry  $w.f.intrace -textvariable entries($w,intrace) -label "Input Physiological Recordings  " -title "Select" -width 40 -filedialog directory  -filetypes *.txt

    # TimeSeries

    set entries($w,invols) ""

    FileEntry  $w.f.invols -textvariable entries($w,invols)  -label "Input TimeSeries (4D)                "  -title "Select" -width 40 -filedialog directory  -filetypes IMAGE

    # Column specifications
    frame $w.f.columns
    set entries($w,colcard) "4"
    set entries($w,colresp) "2"
    set entries($w,coltrig) "3"
    LabelSpinBox $w.f.columns.colcard -label "Column number of data: Cardiac " -textvariable entries($w,colcard) -range {0 100 1} 
    LabelSpinBox $w.f.columns.colresp -label "Respiratory " -textvariable entries($w,colresp) -range {0 100 1} 
    LabelSpinBox $w.f.columns.coltrig -label "Scanner triggers " -textvariable entries($w,coltrig) -range {0 100 1} 
    pack $w.f.columns.colcard $w.f.columns.colresp $w.f.columns.coltrig -in $w.f.columns -side left -padx 3 -pady 3

    

    # sampling etc.
    frame $w.f.triggers
    label $w.f.triggers.idlabel -text "Pulse Ox Triggers     "
    set entries($w,poxtriggers) 0
    checkbutton $w.f.triggers.idbutton -variable entries($w,poxtriggers)
    set entries($w,samp) "200"
    LabelSpinBox $w.f.triggers.sampling -label "    Sampling Rate (Hz) " -textvariable entries($w,samp) -range {0 10000 1} 
    set entries($w,tr) "3.0"
    LabelSpinBox $w.f.triggers.tr -label "         TR (sec) " -textvariable entries($w,tr) -range {0.0 100.0 0.1} 
    pack $w.f.triggers.idbutton $w.f.triggers.idlabel $w.f.triggers.sampling $w.f.triggers.tr -in $w.f.triggers -side left -padx 3 -pady 3

    # slice order
    frame $w.f.sliceorder

    label $w.f.sliceorder.label -text "Slice Order: "

    set entries($w,sliceorder) "up"
    radiobutton $w.f.sliceorder.up -text "Up" \
	    -variable entries($w,sliceorder) -value up -anchor w -command "pnm:slicefile $w"
    radiobutton $w.f.sliceorder.down -text "Down" \
	    -variable entries($w,sliceorder) -value down -anchor w -command "pnm:slicefile $w"
    radiobutton $w.f.sliceorder.intup -text "Interleaved Up" \
	    -variable entries($w,sliceorder) -value interleaved_up -anchor w -command "pnm:slicefile $w"
    radiobutton $w.f.sliceorder.intdown -text "Interleaved Down" \
	    -variable entries($w,sliceorder) -value interleaved_down -anchor w -command "pnm:slicefile $w"
    radiobutton $w.f.sliceorder.usrfile -text "User Specified (via file)" \
	    -variable entries($w,sliceorder) -value usrfile -anchor w -command "pnm:slicefile $w"

    pack $w.f.sliceorder.label $w.f.sliceorder.up $w.f.sliceorder.down $w.f.sliceorder.intup $w.f.sliceorder.intdown $w.f.sliceorder.usrfile -in $w.f.sliceorder -side left -padx 3 -pady 3

    frame $w.f.slicefileframe
    set entries($w,slicefile) ""
    FileEntry  $w.f.slicefileframe.slicefile -textvariable entries($w,slicefile) -label "    User slice timing file   " -title "Select" -width 40 -filedialog directory  -filetypes IMAGE -disabledbackground gray
    $w.f.slicefileframe.slicefile configure -state disabled

    pack $w.f.slicefileframe.slicefile -in $w.f.slicefileframe -side top -anchor w -padx 3 -pady 3

   # slice direction
    frame $w.f.slicedir

    label $w.f.slicedir.label -text "Scanner Slice Direction (wrt voxel axes of the current image): "

    set entries($w,slicedir) z
    radiobutton $w.f.slicedir.x -text "X" \
	    -variable entries($w,slicedir) -value x -anchor w
    radiobutton $w.f.slicedir.y -text "Y" \
	    -variable entries($w,slicedir) -value y -anchor w
    radiobutton $w.f.slicedir.z -text "Z" \
	    -variable entries($w,slicedir) -value z -anchor w

    pack $w.f.slicedir.label $w.f.slicedir.x $w.f.slicedir.y $w.f.slicedir.z -in $w.f.slicedir -side left -padx 3 -pady 3

    # general pack

    pack $w.f.intrace $w.f.invols $w.f.columns $w.f.triggers $w.f.sliceorder $w.f.slicefileframe $w.f.slicedir -in $lfinput -side top -anchor w -pady 3 -padx 5

    #output volume size

    # output volume
    
    set entries($w,outvol) ""

    TitleFrame $w.f.output -text "Output" -relief groove 
    set lfoutput [ $w.f.output getframe ]
   
    FileEntry  $w.f.outvol -textvariable entries($w,outvol) -label "Output Basename   " -title "Select" -width 40 -filedialog directory  -filetypes IMAGE

    pack $w.f.outvol -in $lfoutput -side top -anchor w -pady 3 -padx 5

    # Specify what EVs to generate
    TitleFrame $w.f.evs -text "EVs" -relief groove 
    set lfevs [ $w.f.evs getframe ]

    set entries($w,ocard) "4"
    LabelSpinBox $w.f.evs.card -label " Order for Cardiac EVs                           " -textvariable entries($w,ocard) -range {0 100 1} 
    set entries($w,oresp) "4"
    LabelSpinBox $w.f.evs.resp -label " Order for Respiratory EVs                     " -textvariable entries($w,oresp) -range {0 100 1} 
    set entries($w,omultc) "2"
    LabelSpinBox $w.f.evs.multc -label " Order for Cardiac Interaction EVs         " -textvariable entries($w,omultc) -range {0 100 1} 
    set entries($w,omultr) "2"
    LabelSpinBox $w.f.evs.multr -label " Order for Respiratory Interaction EVs   " -textvariable entries($w,omultr) -range {0 100 1} 

    pack $w.f.evs.card $w.f.evs.resp $w.f.evs.multc $w.f.evs.multr -in $lfevs -side top -anchor w -pady 3 -padx 5

    frame $w.f.evs.buttons
    label $w.f.evs.buttons.rvtlabel -text "RVT             "
    set entries($w,rvt) 0
    checkbutton $w.f.evs.buttons.rvtbutton -variable entries($w,rvt)
    label $w.f.evs.buttons.hrlabel -text "HeartRate"
    set entries($w,hr) 0
    checkbutton $w.f.evs.buttons.hrbutton -variable entries($w,hr)
    label $w.f.evs.buttons.csflabel -text "CSF"
    set entries($w,csf) 0
    checkbutton $w.f.evs.buttons.csfbutton -variable entries($w,csf) -command "pnm:csfmask $w"
    pack $w.f.evs.buttons.rvtbutton $w.f.evs.buttons.rvtlabel $w.f.evs.buttons.hrbutton $w.f.evs.buttons.hrlabel $w.f.evs.buttons.csfbutton $w.f.evs.buttons.csflabel -in $w.f.evs.buttons -side left -padx 3 -pady 3

    set entries($w,csfmask) ""
    FileEntry  $w.f.evs.csfmask -textvariable entries($w,csfmask) -label "CSF mask   " -title "Select" -width 40 -filedialog directory  -filetypes IMAGE -disabledbackground gray
    $w.f.evs.csfmask configure -state disabled
    pack $w.f.evs.buttons $w.f.evs.csfmask -in $lfevs -side top -anchor w -padx 3 -pady 3

    # overall pack

    pack $w.f.input $w.f.output $w.f.evs -in $w.f -side top -anchor w -pady 0 -padx 5


    # advanced options

    # ---- Optional stuff ----

    collapsible frame $w.f.opts -title "Advanced Options"    


    NoteBook $w.nb -side top -bd 2 -tabpady {5 10} -arcradius 3
    $w.nb insert 0 smoothing -text "Smoothing"
    $w.nb insert 1 misc -text  "Miscellaneous"

    $w.nb raise smoothing

    # Smoothing

    set smoothinglf [$w.nb getframe smoothing]

    set entries($w,smoothcard) 0.1    
    set entries($w,smoothresp) 0.1    
    set entries($w,smoothrvt) 10    
    set entries($w,smoothhr) 10    

    LabelSpinBox $w.smoothc -label " Cardiac Smoothing (sec)    " -textvariable entries($w,smoothcard) -range {0 5000 0.01 } 
    LabelSpinBox $w.smoothr -label " Respiratory Smoothing (sec)    " -textvariable entries($w,smoothresp) -range {0 5000 0.01 } 
    LabelSpinBox $w.smoothhr -label " Heart Rate Smoothing (sec)    " -textvariable entries($w,smoothhr) -range {0 5000 0.5 } 
    LabelSpinBox $w.smoothrvt -label " RVT Smoothing (sec)    " -textvariable entries($w,smoothrvt) -range {0 5000 0.5 } 
    
    # ---- pack ----
    pack $w.smoothc $w.smoothr $w.smoothhr $w.smoothrvt -in $smoothinglf -side top -anchor w -padx 3 -expand yes

    # Misc

    set misclf [$w.nb getframe misc]

    set entries($w,cleanup) 1
    checkbutton $w.cleanup  -text " Apply cleanup stages   " -variable entries($w,cleanup)
    set entries($w,invresp) 0
    checkbutton $w.invresp  -text " Invert Respiratory Trace" -variable entries($w,invresp)
    set entries($w,invcard) 0
    checkbutton $w.invcard  -text " Invert Cardiac Trace" -variable entries($w,invcard)

    # ---- pack ----

    pack $w.cleanup $w.invresp $w.invcard -in $misclf -side top -anchor w -padx 3 -expand yes

    set entries($w,cleanup) 1


    # ---- pack ----

    frame $w.f.advopts
    pack $w.nb -in $w.f.advopts -side top
    pack $w.f.advopts -in $w.f.opts.b -side left -padx 8 -pady 6 -expand yes -fill both
    pack $w.f.opts -in $w.f -side left -padx 5 -pady 5



    # ---- Button Frame ----

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.apply     -command "pnm:apply $w" \
	    -text "Go" -width 5

    button $w.cancel    -command "destroy $w" \
	    -text "Exit" -width 5

    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/redirects/feat.html" \
	    -text "Help" -width 5

    pack $w.btns.b -side bottom -fill x
    pack $w.apply $w.cancel $w.help -in $w.btns.b \
	    -side left -expand yes -padx 3 -pady 10 -fill y
    
    pack $w.f $w.btns -expand yes -fill both
}


proc pnm:csfmask { w } {
    global entries
    if { $entries($w,csf) == 1 } {
	$w.f.evs.csfmask configure -state normal
    } else {
	$w.f.evs.csfmask configure -state disabled
    }
    update idletasks
}

proc pnm:slicefile { w } {
    global entries
    if { $entries($w,sliceorder) == "usrfile" } {
	$w.f.slicefileframe.slicefile configure -state normal
    } else {
	$w.f.slicefileframe.slicefile configure -state disabled
    }
    update idletasks
}


proc pnm:apply { w } {
    global entries

    set status [ pnm:proc $entries($w,intrace) $entries($w,invols) $entries($w,colcard) $entries($w,colresp) $entries($w,coltrig) $entries($w,poxtriggers) $entries($w,samp) $entries($w,tr) $entries($w,sliceorder) $entries($w,slicedir) $entries($w,slicefile) $entries($w,outvol) $entries($w,ocard) $entries($w,oresp) $entries($w,omultc) $entries($w,omultr) $entries($w,rvt) $entries($w,hr) $entries($w,csf) $entries($w,csfmask) $entries($w,smoothcard) $entries($w,smoothresp) $entries($w,smoothhr) $entries($w,smoothrvt) $entries($w,cleanup) $entries($w,invresp) $entries($w,invcard) ]

    update idletasks
    puts "Done"
}

#}}}
#{{{ proc

proc pnm:proc { infile invol colcard colresp coltrig poxtrig samprate tr sliceorder slicedir slicefile outname oc or multc multr rvt hr csf csfmask smoothc smoothr smoothhr smoothrvt cleanup invr invc } {

    global FSLDIR
    
    set poppargs "-i ${outname}_input.txt -o ${outname} -s $samprate --tr=$tr --smoothcard=$smoothc --smoothresp=$smoothr --resp=$colresp --cardiac=$colcard --trigger=$coltrig"
    set pnmcommand "${FSLDIR}/bin/pnm_evs -i $invol -c ${outname}_card.txt -r ${outname}_resp.txt -o ${outname} --tr=$tr --oc=$oc --or=$or --multc=$multc --multr=$multr"    

    if { $poxtrig == 1 } { set poppargs "$poppargs --pulseox_trigger" }
    if { $rvt == 1 } { 
	set poppargs "$poppargs --rvt" 
	set pnmcommand "$pnmcommand --rvt=${outname}_rvt.txt --rvtsmooth=${smoothrvt}"
    }
    if { $hr == 1 } { 
	set poppargs "$poppargs --heartrate" 
	set pnmcommand "$pnmcommand --heartrate=${outname}_hr.txt --heartratesmooth=${smoothhr}"
    }
    if { $csf == 1 } { set pnmcommand "$pnmcommand --csfmask=$csfmask" }
    if { $cleanup == 0 } { set poppargs "$poppargs --noclean1 --noclean2" }
    if { $invr == 1 } { set poppargs "$poppargs --invertresp" }
    if { $invc == 1 } { set poppargs "$poppargs --invertcardiac" }
    set pnmsliceconf "--sliceorder=$sliceorder --slicedir=$slicedir"
    if { $sliceorder == "usrfile" } { set pnmsliceconf "--slicetiming=$slicefile" }
    set pnmcommand "$pnmcommand $pnmsliceconf"

    puts "$FSLDIR/bin/fslFixText $infile ${outname}_input.txt"
    puts "$FSLDIR/bin/pnm_stage1 $poppargs"
    #puts "$pnmcommand"

    set pnms3 [open "${outname}_pnm_stage3" "w"] 
    puts $pnms3 "#!/bin/sh"
    puts $pnms3 "$pnmcommand"
    puts $pnms3 "ls -1 `$FSLDIR/bin/imglob -extensions \$\{obase\}ev0*` > ${outname}_evlist.txt"
    close $pnms3
    fsl:exec "$FSLDIR/bin/fslFixText $infile ${outname}_input.txt"
    fsl:exec "$FSLDIR/bin/pnm_stage1 $poppargs"

    # MJ TODO - this doesn't seem to wait until commands are run!
    #if { [ imtest ${outname}ev001 ] == 0 } {
	#puts "No output saved!"
	#return 4
    #}

    return 0
}

#}}}

wm withdraw .
pnm .rename
tkwait window .rename

