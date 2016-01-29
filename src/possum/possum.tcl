# Possum - GUI for simulating FMRI
#
# Ivana Drobnjak and Mark Jenkinson, FMRIB Image Analysis Group
#
# Copyright (C) 2006-2007 University of Oxford
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


source $FSLDIR/tcl/fslstart.tcl
set VARS(history) {}

# RE:POSSUMDIR
# For all the users POSSUMDIR will be empty and therefore automatically become FSLDIR
# For me POSSUMDIR is my FSLDEVDIR directory. This allows me to run POSSUM on the cluster 
# without having to make it stable and wait for a day or find way to run my binaries. 
# It just seemed the easiest way to get around this. It is in all POSSUM scripts so please
# leave it that way if possible.
if [ info exists env(POSSUMDIR) ] {
    set POSSUMDIR $env(POSSUMDIR)
} else {
   set POSSUMDIR $FSLDIR
}
# The following two lines are just for me. Erase for the stable version. 
puts ""
puts "Possum is running from POSSUMDIR=$POSSUMDIR"
proc possum { w } {
    global entries guivars FSLDIR PWD HOME 
    # ---- Set up Frames ----
    toplevel $w
    wm title $w "POSSUM"
    wm iconname $w "Possum"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
    #    tixBalloon    $w.bhelp
    frame $w.f
    NoteBook $w.nb -side top -bd 2 -tabpady {5 5} -arcradius 3
    $w.nb insert 0 object -text "Object"
    $w.nb insert 1 pulse -text "Pulse Sequence"
    $w.nb insert 2 b0field -text "B0 Field"
    $w.nb insert 3 motion -text "Motion"
    $w.nb insert 4 activation -text "Activation"
    $w.nb insert 5 noise -text "Noise"
    $w.nb insert 6 output -text "Run POSSUM"
    $w.nb raise object
    
    #----- Object-------
    set objectlf [$w.nb getframe object]
    set entries($w,obvol) ${FSLDIR}/data/possum/brain.nii.gz
    possum:updateOBprop $w ${FSLDIR}/data/possum/brain.nii.gz
    frame $w.objim
    FileEntry $w.obvol \
	-textvariable entries($w,obvol) \
	-filetypes IMAGE \
	-label "Input Object      " \
	-title "Select" \
	-width 40 \
	-filedialog directory \
        -command  "possum:updateOBprop $w; possum:updatecomptime $w "
    pack $w.obvol -in $objectlf -anchor w -padx 3 -pady 3
    label $w.objim.label -image "" -text " "
#tejas-edit: change command to '' to get back original
    button $w.objim.preview -text "Preview Image" -command "possum:previewimagedialog $w"
#tejas-end
#$w.objim.label
    pack $w.objim.preview -in $w.objim -pady 10
    pack $w.objim -in  $objectlf -anchor w -padx 3 -pady 3

    #-------- Pulse sequence -----------
    set pulself [$w.nb getframe pulse]
#    LabelFrame $w.pul -text "EPI" -font {Helvetica 11 bold}
#    pack $w.pul -in $pulself -side top -anchor w -padx 3 -pady 3
    
    # set up default values
    set entries($w,te) 0.030
    set entries($w,tr) 3
    set entries($w,trslc) 0.12
    set entries($w,autotrslc) 1
    set entries($w,outsize_nx) 64
    set entries($w,outsize_ny) 64
    set entries($w,outsize_nz) 1   
    set entries($w,outsize_dx) 4.0
    set entries($w,outsize_dy) 4.0
    set entries($w,outsize_dz) 1.0
    set entries($w,numvol) 1
    set entries($w,gap) 0
    set entries($w,bw) 100000
    set entries($w,zstart) 70
    set entries($w,readgrad) x
    set entries($w,phencode) y
    set entries($w,slcselect) z
    set entries($w,plus) +
    set entries($w,pluss) +
    set entries($w,maxG) 0.055
    set entries($w,riseT) 0.00022
    set entries($w,slcprof) "$FSLDIR/data/possum/slcprof"
    set entries($w,numproc) 1
    set entries($w,segs) 10000
    set entries($w,ctt) 0
    set entries($w,comptime) 0
    set entries($w,motion_yn) 0
    set entries($w,pulsechecktest) 1
    set entries($w,b0inh_yn) 0
    set entries($w,b0inhtime_yn) 0
    set entries($w,cover) 100
    set entries($w,flipangle) 90
#tejas - 2.11.12
    set entries($w,seqtype) "epi"
    set entries($w,onlyCheck) 1
    set entries($w,usephantomdims) 0
#    set entries($w,rfavg_yn) 0
    set entries($w,makeb0motion_yn) 0
    set entries($w,custompulse_yn) 0
#tejas-end
    set entries($w,slcsampfactor) 2

    # calculate dependendent quantities from the defaults
    possum:updateTRSLC $w
    possum:updateFOV $w 
    possum:updateechosp $w
    possum:updatecomptime $w 
    
    # set up the GUI widgets
    frame $w.topopts
    frame $w.pulseseq
	array set temp [list a 1 b 2]
	label $w.pulseseq.seqlab -text "Sequence type:" -width 15 -anchor w
	optionMenu2 $w.pulseseq.seqopt entries($w,seqtype) -command "possum:settimings $w; possum:updateTRSLC $w; possum:updateTRSLC2 $w; possum:custompulse $w" "epi" "Echo Planar Imaging (EPI)" "ge" "Gradient Echo (GE)" "custom" "Custom Pulse Sequence"
    pack $w.pulseseq.seqlab $w.pulseseq.seqopt -in $w.pulseseq -side left -anchor w

    frame $w.t
    label $w.t.spc1 -text "" -width 0
    LabelSpinBox $w.t.x -label " TE (s): " -width 8 \
         -textvariable entries($w,te) -range {0.0 10000.0 0.001}
    label $w.t.spc2 -text "" -width 0
    LabelSpinBox $w.t.y -label " TR (s): " -width 8 \
         -textvariable entries($w,tr) -range {0.0 10000.0 0.001} \
	-command "$w.t.y.spin.e validate; possum:updateTRSLC $w" \
	-modifycmd "possum:updateTRSLC $w"
    pack $w.t.spc1 $w.t.x $w.t.spc2 $w.t.y -in $w.t -side left -anchor nw -padx 10 -pady 5
    pack $w.pulseseq $w.t -in $w.topopts -side left -anchor nw

    frame $w.defaultframe
    frame $w.n
    label $w.n.lab -text "Number of Voxels: " -width 15 -anchor w -justify left 
    LabelSpinBox $w.n.x -label " X "  -width 6 \
         -textvariable entries($w,outsize_nx) -range { 1   10000  1 } \
	-command "$w.n.x.spin.e validate; possum:updateFOV $w; possum:updatecomptime $w; possum:updateechosp $w" \
	-modifycmd "possum:updateFOV $w; possum:updatecomptime $w; possum:updateechosp $w"
    LabelSpinBox $w.n.y -label " Y "  -width 6 \
	 -textvariable entries($w,outsize_ny) -range { 1   10000 1 } \
	-command "$w.n.y.spin.e validate; possum:updateFOV $w;possum:updatecomptime $w;possum:updateechosp $w" \
	 -modifycmd " possum:updateFOV $w;possum:updatecomptime $w;possum:updateechosp $w"
    LabelSpinBox $w.n.z -label " Z "  -width 6 \
	 -textvariable entries($w,outsize_nz) -range { 1   10000 1 } \
	-command "$w.n.z.spin.e validate; possum:updateFOV $w;possum:updatecomptime $w;possum:updateechosp $w; possum:updateTRSLC $w" \
	 -modifycmd " possum:updateFOV $w;possum:updatecomptime $w;possum:updateechosp $w; possum:updateTRSLC $w"

    label $w.n.spc1 -text " " -width 3
    label $w.n.nvollab -text "Number of Volumes: " -width 22 -anchor w -justify left 
    LabelSpinBox $w.n.nvolx -label " " -width 6 \
	 -textvariable entries($w,numvol) -range { 1   10000  1 } \
         -command "$w.n.nvolx.spin.e validate; possum:updatecomptime $w" \
	 -modifycmd "possum:updatecomptime $w"
#    pack $w.n.lab $w.n.x $w.n.y $w.n.z -in $w.n -side left -anchor w -padx 3 -pady 3
     pack $w.n.lab $w.n.x $w.n.y $w.n.z $w.n.spc1 $w.n.nvollab $w.n.nvolx -in $w.n -side left -anchor w -padx 3 -pady 3
     
    frame $w.d
    label $w.d.lab -text "Voxel Size (mm): " -width 15 -anchor w -justify left 
    LabelSpinBox $w.d.x -label " X " -width 6 \
	 -textvariable entries($w,outsize_dx) -range { 0.000001  10000.0 0.1 } \
	-command "$w.d.x.spin.e validate; possum:updateFOV $w; possum:updatecomptime $w;possum:updateechosp $w" \
	-modifycmd " possum:updateFOV $w; possum:updatecomptime $w;possum:updateechosp $w"
   
    LabelSpinBox $w.d.y -label " Y "  -width 6 \
	 -textvariable entries($w,outsize_dy) -range { 0.000001   10000.0  0.1 } \
	-command "$w.d.y.spin.e validate; possum:updateFOV $w;possum:updatecomptime $w;possum:updateechosp $w" \
	-modifycmd " possum:updateFOV $w; possum:updatecomptime $w;possum:updateechosp $w"
   
    LabelSpinBox $w.d.z -label " Z "  -width 6 \
	 -textvariable entries($w,outsize_dz) -range { 0.000001   10000.0  0.1 } \
	-command "$w.d.z.spin.e validate; possum:updateFOV $w; possum:updatecomptime $w;possum:updateechosp $w" \
	-modifycmd "possum:updateFOV $w; possum:updatecomptime $w;possum:updateechosp $w"

    label $w.d.spc1 -text " " -width 3
    label $w.d.fliplab -text "Flip angle (deg): " -width 22 -anchor w -justify left 
    LabelSpinBox $w.d.flipx -label " " -width 6 \
	    -textvariable entries($w,flipangle) -range { 0   180  0.1 } 

    pack $w.d.lab $w.d.x $w.d.y $w.d.z $w.d.spc1 $w.d.fliplab $w.d.flipx -in $w.d -side left -anchor w -padx 3 -pady 3 
    
    frame $w.fov
    label $w.fov.lab -text "Field of view (mm): " -width 15 -anchor w -justify left 
    LabelSpinBox $w.fov.x -label " X "  -width 6 \
	 -textvariable entries($w,fov_x) -range { 0.000001   10000.0  0.1 }\
         -command "$w.fov.x.spin.e validate; possum:updateVSIZE $w; possum:updatecomptime $w;possum:updateechosp $w" \
	 -modifycmd "possum:updateVSIZE $w; possum:updatecomptime $w;possum:updateechosp $w"
    LabelSpinBox $w.fov.y -label " Y "  -width 6 \
	 -textvariable entries($w,fov_y) -range { 0.000001   10000.0  0.1 }\
         -command "$w.fov.y.spin.e validate; possum:updateVSIZE $w; possum:updatecomptime $w;possum:updateechosp $w" \
	 -modifycmd "possum:updateVSIZE $w; possum:updatecomptime $w;possum:updateechosp $w"
    LabelSpinBox $w.fov.z -label " Z "  -width 6 \
	 -textvariable entries($w,fov_z) -range { 0.000001   10000.0  0.1 }\
         -command "$w.fov.z.spin.e validate; possum:updateVSIZE $w; possum:updatecomptime $w;possum:updateechosp $w" \
	 -modifycmd "possum:updateVSIZE $w; possum:updatecomptime $w;possum:updateechosp $w"

    label $w.fov.spc1 -text " " -width 3
    label $w.fov.stslclab -text "Starting slice position (mm): " -width 22 -anchor w -justify left 
    LabelSpinBox $w.fov.stslcx -label " " -width 6 \
	    -textvariable entries($w,zstart) -range { 0.0   10000.0  1.0 }

    pack $w.fov.lab $w.fov.x $w.fov.y $w.fov.z $w.fov.spc1 $w.fov.stslclab $w.fov.stslcx -in $w.fov -side left -anchor w -padx 3 -pady 3  

    frame $w.s
	    label $w.s.spc2 -text " " -width 3
	    button $w.s.pulsecheck     -command "Possum:pulsecheck $w 1" -text "Consistency check" -width 15
	    label $w.s.spc3 -text " " -width 5
	    button $w.s.preview -text "Preview slice prescription" -command "possum:previewslices $w" -anchor w
	    label $w.s.spc4 -text " " -width 5
	    button $w.s.createpulse -text "Save Pulse Sequence" -command "possum:savepulsedialog $w" -anchor w
    pack $w.s.spc2 $w.s.preview $w.s.spc3 $w.s.pulsecheck $w.s.spc4 -in $w.s -side left -anchor center -padx 3 -pady 3
    pack $w.s.createpulse -in $w.s

    collapsible frame $w.slp -title "Slice Prescription" -command "$w.nb compute_size; set dummy"
    frame $w.sr
	    label $w.sr.lab -text "Read gradient: " -width 22 -anchor w -justify left 
	    radiobutton $w.sr.x -text "X" -variable entries($w,readgrad) -value x -anchor w
	    radiobutton $w.sr.y -text "Y" -variable entries($w,readgrad) -value y -anchor w
	    radiobutton $w.sr.z -text "Z" -variable entries($w,readgrad) -value z -anchor w
	    $w.sr.x select
	    label $w.sr.spc1 -text " " -width 5
	    label $w.sr.dirlab -text "Slice acquisition order: " -width 20 -anchor w -justify left 
	    radiobutton $w.sr.dirx -text "+" -variable entries($w,plus) -value + -anchor w
	    radiobutton $w.sr.diry -text "-" -variable entries($w,plus) -value - -anchor w
	    $w.sr.dirx select
    pack $w.sr.lab $w.sr.x $w.sr.y $w.sr.z $w.sr.spc1 $w.sr.dirlab $w.sr.dirx $w.sr.diry -in $w.sr -side left -anchor w -padx 3 -pady 3
 
    frame $w.sp    
	    label $w.sp.lab -text "Phase encode gradient: " -width 22 -anchor w -justify left 
	    radiobutton $w.sp.x -text "X" -variable entries($w,phencode) -value x -anchor w
	    radiobutton $w.sp.y -text "Y" -variable entries($w,phencode) -value y -anchor w
	    radiobutton $w.sp.z -text "Z" -variable entries($w,phencode) -value z -anchor w
	    $w.sp.y select
	    label $w.sp.spc2 -text " " -width 5
	    label $w.sp.phlab -text "k-space coverage (%):" -width 20 -anchor w -justify left 
	    LabelSpinBox $w.sp.phx -label "  " -width 6 \
		 -textvariable entries($w,cover) -range { 50   100  0.1 } \
		 -command "$w.ph.x.spin.e validate; possum:updatecomptime $w" \
		 -modifycmd "possum:updatecomptime $w"
    pack $w.sp.lab $w.sp.x $w.sp.y $w.sp.z $w.sp.spc2 $w.sp.phlab $w.sp.phx -in $w.sp -side left -anchor w -padx 3 -pady 3
   
    frame $w.ss
    label $w.ss.lab -text "Slice select gradient: " -width 22 -anchor w -justify left 
    radiobutton $w.ss.x -text "X" -variable entries($w,slcselect) -value x -anchor w
    radiobutton $w.ss.y -text "Y" -variable entries($w,slcselect) -value y -anchor w
    radiobutton $w.ss.z -text "Z" -variable entries($w,slcselect) -value z -anchor w
    $w.ss.z select
    label $w.ss.spc3 -text " " -width 5
    label $w.ss.gaplab -text "Gap between slices (mm):" -width 20 -anchor w -justify left 
    LabelSpinBox $w.ss.gapv -label "  " -width 6 \
	-textvariable entries($w,gap) -range { 0.0   100.0  0.001 } \
	-command "$w.gap.v.spin.e validate; possum:updatecomptime $w" \
	-modifycmd "possum:updatecomptime $w"
    pack $w.ss.lab $w.ss.x $w.ss.y $w.ss.z $w.ss.spc3 $w.ss.gaplab $w.ss.gapv -in $w.ss -side left -anchor w -padx 3 -pady 3

    frame $w.ssf
    LabelSpinBox $w.ssf.fac -label " Object slice resampling factor: " -width 8 \
         -textvariable entries($w,slcsampfactor) -range {1 100} -modifycmd "possum:updatecomptime $w"
    pack $w.ssf.fac -in $w.ssf -side left -anchor w -padx 3 -pady 3

    pack $w.sr $w.sp $w.ss $w.ssf -in $w.slp.b -anchor w -padx 3 -pady 3

    collapsible frame $w.scan -title "Scanner properties" -command "$w.nb compute_size; set dummy"
    
    frame $w.scanfr1
    label $w.scanfr1.maxGlab -text "Maximal gradient strength (T/m):" -width 28 -anchor w
    LabelSpinBox $w.scanfr1.maxGv -label " " -width 10 \
	 -textvariable entries($w,maxG) -range { 0.0   100.0  0.001 }

    label $w.scanfr1.spc1 -text " " -width 10
    label $w.scanfr1.riseTlab -text "Rise time (s): " -width 15 -anchor w
    LabelSpinBox $w.scanfr1.riseTv -width 10 -label " " \
	 -textvariable entries($w,riseT) -range { 0.0   100.0  0.00001 }
    pack $w.scanfr1.maxGlab $w.scanfr1.maxGv $w.scanfr1.spc1 $w.scanfr1.riseTlab $w.scanfr1.riseTv -in $w.scanfr1 -side left -anchor w -pady 3 

    frame $w.bw
    label $w.bw.lab -text "BW (Hz): " -width 28 -anchor w
    LabelSpinBox $w.bw.x -width 10 -label " " \
	    -textvariable entries($w,bw) -range { 0   1000000  10 }\
            -command "$w.bw.x.spin.e validate; possum:updateechosp $w" \
	    -modifycmd "possum:updateechosp $w"

    label $w.bw.spc1 -text " " -width 10
    label $w.bw.echolab -text "Echo spacing (s):   " -width 15 -anchor w
    label $w.bw.echox -textvariable entries($w,echosp) -width 10 
#-readonlybackground white -state readonly
    pack $w.bw.lab $w.bw.x $w.bw.spc1 $w.bw.echolab $w.bw.echox -in $w.bw  -anchor e -side left

    frame $w.slcprof
    label $w.slcprof.lab -text "Slice profile " -width 10 -anchor w
    FileEntry $w.slcprof.dir \
	-textvariable entries($w,slcprof) \
	-title "Select" \
	-width 40 \
	-filedialog directory
    label $w.slcprof.spc1 -text "" -width 2
    button $w.slcprof.prevslcbutton -command "possum:slcprofpreviewdialog $w" \
	    -text "Preview Slice Profile" -width 15
    pack $w.slcprof.lab $w.slcprof.dir $w.slcprof.spc1 $w.slcprof.prevslcbutton -in $w.slcprof  -anchor e -side left -pady 3

    frame $w.trs
    #checkbutton $w.trs.rfavg -text "RF Averaging" -variable entries($w,rfavg_yn) -padx 5 -justify center
    #label $w.trs.spc1 -text "" -width 2
    label $w.trs.lab -text "TRslice (s):" -width 10 -anchor w
    LabelSpinBox $w.trs.z -label " " -width 8 \
         -textvariable entries($w,trslc) -range {0.0 10000.0 0.001} -disabledbackground gray
    checkbutton $w.trs.yn -text "Autoset" -variable entries($w,autotrslc) -command "possum:buttonTRSLC $w" -padx 5
    possum:buttonTRSLC $w
    #pack $w.trs.rfavg $w.trs.spc1 $w.trs.lab $w.trs.z $w.trs.yn -in $w.trs -side left -pady 3
    pack $w.trs.lab $w.trs.z $w.trs.yn -in $w.trs -side left -pady 3

    pack $w.scanfr1 $w.bw $w.slcprof $w.trs -in $w.scan.b -anchor w -padx 3 -pady 3 -expand yes -anchor nw
    pack $w.n $w.d $w.fov $w.s $w.slp $w.scan -in  $w.defaultframe -anchor nw -padx 3 -pady 3
    pack $w.topopts $w.defaultframe -in $pulself -anchor nw -side top -padx 3 -pady 3

    ## Frame for custom pulse sequence
    frame $w.custompulse
	frame $w.inpulse
		label $w.inpulse.lab -text "Pulse-sequence: " -width 15 -justify left -anchor w
		FileEntry $w.inpulse.dir \
			-textvariable entries($w,cuspulse) \
			-title "Custom Pulse-Sequence" \
			-width 40 \
			-filedialog directory \
			-command "possum:checkloadedpulse $w"
	frame $w.customwarn
		label $w.customwarn.lab -text "" -font { Helvetica 10 italic }
	##
	pack $w.inpulse.lab $w.inpulse.dir -in $w.inpulse -side left -anchor w -pady 3
	collapsible frame $w.cusadv -title "Scanner Properties" -command "set dummy"
#		checkbutton $w.rfavg -text "RF Averaging" -variable entries($w,rfavg_yn) -padx 5 -justify center
		frame $w.slcprofcus
			label $w.slcprofcus.lab -text " Slice profile " -width 12 -anchor w
			FileEntry $w.slcprofcus.dir \
				-textvariable entries($w,slcprof) \
				-title "Select" \
				-width 40 \
				-filedialog directory
			label $w.slcprofcus.spc1 -text "" -width 2
			button $w.slcprofcus.prevslcbutton -command "possum:slcprofpreviewdialog $w" \
				-text "Preview Slice Profile" -width 15
		pack $w.slcprofcus.lab $w.slcprofcus.dir $w.slcprofcus.spc1 $w.slcprofcus.prevslcbutton -in $w.slcprofcus  -anchor w -side left -pady 3
#	pack $w.slcprofcus $w.rfavg -in $w.cusadv.b -anchor w -padx 3 -pady 3
	pack $w.slcprofcus -in $w.cusadv.b -anchor w -padx 3 -pady 3
    pack $w.inpulse $w.customwarn $w.cusadv -in $w.custompulse -anchor w -side top -padx 3 -pady 3

    # -------- B0ield -------------
    set guivars($w,lfb0field) [$w.nb getframe b0field]

    #Field strength
    set entries($w,b0strength) 1
    possum:updateb0field $w
    frame $w.b0test
    LabelSpinBox $w.b0test.b0spin -width 8 \
       -textvariable entries($w,b0fieldstrength) -range { 0.0   1000000.0  0.1 } -disabledbackground gray
    LabelFrame $w.b0test.b0fil -text "Field strength                "
    optionMenu2 $w.b0test.b0fil.menu entries($w,b0strength)  -command "possum:updateb0field $w ; possum:updateMRpar $w ; possum:updateb0fieldinh $w ; possum:updateBASEname $w; possum:updateBASEnametime $w; pack forget $w.mrpar.prev; $w.nb compute_size " 0 "1.5 T" 1 "3 T" 2 "Custom field ( T )"
    pack $w.b0test.b0fil.menu
    pack $w.b0test.b0fil -in $w.b0test
    pack $w.b0test -in $guivars($w,lfb0field) -side top -anchor w -padx 3 -pady 3

    #-------------------B0Field------------------------#
    # MR par
    frame $w.mrpar
	    frame $w.mrparopts
	    	    set entries($w,mrpar) "${FSLDIR}/data/possum/MRpar_3T"
		    FileEntry $w.mrparopts.dir \
			-textvariable entries($w,mrpar) \
			-label "MR parameters:            " \
			-title "Select" \
			-width 40 \
			-filedialog directory \
			-command "pack forget $w.mrpar.prev"
		    label $w.mrparopts.spc1 -text "" -width 3
		    button $w.mrparopts.view -command "possum:showMRpar $w"  -text "View MRpar File" -width 22
	    pack $w.mrparopts.dir $w.mrparopts.spc1 $w.mrparopts.view -in $w.mrparopts -anchor w -side left -padx 3 -pady 3
	    frame $w.mrpar.prev -relief raised -borderwidth 1
		label $w.mrpar.prev.contents1 -text "" -justify center
		label $w.mrpar.prev.contents2 -text "" -justify center
		label $w.mrpar.prev.contents3 -text "" -justify center
		label $w.mrpar.prev.contents4 -text "" -justify center
	    pack $w.mrpar.prev.contents1 $w.mrpar.prev.contents2 $w.mrpar.prev.contents3 $w.mrpar.prev.contents4 -in $w.mrpar.prev -anchor w -side left -padx 5 -pady 3
    pack $w.mrparopts -in $w.mrpar -pady 3 -anchor w -side top
    pack $w.mrpar -in $guivars($w,lfb0field) -anchor w -padx 3 -pady 3

    #B0 field inhomogeneities
    frame $w.b0main
    frame $w.b0mainopts
	set entries($w,b0f) "${FSLDIR}/data/possum/b0_ppm.nii.gz"
	possum:updateB0prop $w  $entries($w,b0f)
	FileEntry $w.b0mainopts.dir \
		-textvariable entries($w,b0f) \
		-label "File name:                    " \
		-filetypes IMAGE \
	      	-title "Select" \
		-width 40 \
		-filedialog directory \
		-command  "possum:updateB0prop $w; possum:updatecomptime $w "
	label $w.b0mainopts.spc1 -text "" -width 3
	button $w.b0mainopts.preview -text "Preview Image" -command "possum:previewb0dialog $w" -width 22
    pack $w.b0mainopts.dir $w.b0mainopts.spc1 $w.b0mainopts.preview -in $w.b0mainopts -anchor w -side left -padx 3 -pady 3

    frame $w.b0main.b0u
    radiobutton $w.b0main.b0u.ppm -text "ppm" -variable entries($w,b0units) -value ppm -anchor w
    radiobutton $w.b0main.b0u.tesla -text "Tesla" -variable entries($w,b0units) -value tesla -anchor w
    $w.b0main.b0u.ppm select
    label $w.b0main.b0u.unit -text "Units (for file): "
    pack $w.b0main.b0u.unit $w.b0main.b0u.ppm $w.b0main.b0u.tesla -anchor w -side left

    ## Widgets for standard B0 files
    LabelFrame $w.b0main.b0fi -text "B0 field inhomogeneities (at 1 Tesla)             "
    optionMenu2 $w.b0main.b0fi.menu entries($w,b0inh_yn)  -command "possum:updateb0fieldinh $w; possum:updateBASEname $w ; possum:updatecomptime $w; $w.nb compute_size" 0 "None" 1 "Custom file"
    pack $w.b0main.b0fi.menu
    pack $w.b0main.b0fi -in $w.b0main -side top -anchor w -padx 3 -pady 3

    pack $w.b0main -in $guivars($w,lfb0field) -side top -anchor w -padx 3 -pady 3

    #B0 inhomogeneities changing in time
    frame $w.b0maintime
	frame $w.b0timecourse
	    set entries($w,b0ftimecourse) "${FSLDIR}/data/possum/b0timecourse"
	    FileEntry $w.b0timecourse.dir \
		-textvariable entries($w,b0ftimecourse) \
		-label "B0 time course:            " \
	      	-title "Select" \
		-width 40 \
		-filedialog directory
	    label $w.b0timecourse.spc1 -text "" -width 3
	    button $w.b0timecourse.prev -text "Preview B0 Timecourse" -width 22 -command "possum:previewb0timecoursedialog $w"
	pack $w.b0timecourse.dir $w.b0timecourse.spc1 $w.b0timecourse.prev -in $w.b0timecourse -anchor w -side left -padx 3 -pady 3

	frame $w.b0timeimage
	    set entries($w,b0ftime) "${FSLDIR}/data/possum/b0extra.nii.gz"
	    possum:updateB0timeprop $w  $entries($w,b0ftime)
	    FileEntry $w.b0timeimage.dir \
		-textvariable entries($w,b0ftime) \
		-label "B0 spatial modulation:  " \
		-filetypes IMAGE \
	      	-title "Select" \
		-width 40 \
		-filedialog directory \
	        -command  "possum:updateB0timeprop $w; possum:updatecomptime $w "
	    label $w.b0timeimage.spc1 -text "" -width 3
	    button $w.b0timeimage.prev -text "Preview B0 Spatial Modulation" -width 22 -command "possum:previewb0timedialog $w"
	pack $w.b0timeimage.dir $w.b0timeimage.spc1 $w.b0timeimage.prev -in $w.b0timeimage -anchor w -side left -padx 3 -pady 3

#############	NOT USED CURRENTLY	###############
    FileEntry $w.b0maintime.b0f \
	-textvariable entries($w,b0ftime) \
	-label "B0 spatial modulation:  " \
	-filetypes IMAGE \
      	-title "Select" \
	-width 40 \
	-filedialog directory \
        -command  "possum:updateB0timeprop $w; possum:updatecomptime $w "
    FileEntry $w.b0maintime.b0ftime \
	-textvariable entries($w,b0ftimecourse) \
	-label "B0 time course:  " \
      	-title "Select" \
	-width 40 \
	-filedialog directory \

    frame $w.b0maintime.b0im
    label $w.b0maintime.b0im.label -image "" -text " "
    button $w.b0maintime.b0im.preview -text "Preview Image" -command "possum:previewb0timedialog $w"
    label $w.b0maintime.b0im.b0timewarn -text "" -font { Helvetica 9 italic } -height 2 -justify left -foreground red
########################################################

    frame $w.b0maintime.b0u
    set entries($w,b0unitstime) tesla
    radiobutton $w.b0maintime.b0u.ppm -text "ppm" -variable entries($w,b0unitstime) -value ppm -anchor w
    radiobutton $w.b0maintime.b0u.tesla -text "Tesla" -variable entries($w,b0unitstime) -value tesla -anchor w
    label $w.b0maintime.b0u.unit -text "Units (for file): "
    pack $w.b0maintime.b0u.unit $w.b0maintime.b0u.ppm $w.b0maintime.b0u.tesla -anchor w -side left

    LabelFrame $w.b0maintime.b0fi -text "B0 field inhomogeneities (changing in time)   "
    optionMenu2 $w.b0maintime.b0fi.menu entries($w,b0inhtime_yn)  -command "possum:updateb0fieldinhtime $w; possum:updateBASEnametime $w ; possum:updatecomptime $w ; possum:updatemotion $w; $w.nb compute_size" 0 "None" 1 "Custom file"
    pack $w.b0maintime.b0fi.menu
    pack $w.b0maintime.b0fi -in $w.b0maintime  -side top -anchor w -padx 3 -pady 3
    pack $w.b0maintime -in $guivars($w,lfb0field) -side top -anchor w -padx 3 -pady 3 

    # --------Motion-------------
    set guivars($w,lfmotion) [$w.nb getframe motion]
    set entries($w,mot) "${FSLDIR}/data/possum/motionRzLarge_0.12s"
    FileEntry $w.mot \
	-textvariable entries($w,mot) \
	-label "Motion file  " \
	-title "Select" \
	-width 40 \
	-filedialog directory \
	-command "pack forget $w.motprevimg"
    button $w.motpreview -text "Preview Motion Plot" -command "possum:previewmotionplot $w"
    label $w.motwarn -text "" -font { Helvetica 9 italic } -height 2 -justify left -foreground red
    LabelFrame $w.moti -text ""
    optionMenu2 $w.moti.menu entries($w,motion_yn)  -command "possum:updatemotion $w; possum:updateb0fieldinhtime $w ; possum:updatecomptime $w; $w.nb compute_size " 0 "None" 1 "Custom file"
    label $w.motprevimg -text ""
    pack $w.moti.menu
    pack $w.moti -in $guivars($w,lfmotion) -side top -anchor w -padx 3 -pady 3

    # --------Activation-------------
    set guivars($w,lfactivation) [$w.nb getframe activation]
    set entries($w,activ_yn) 0
    set entries($w,act1) "$FSLDIR/data/possum/activation3D.nii.gz"
    possum:updateACTprop $w $entries($w,act1)
    set entries($w,act2) "$FSLDIR/data/possum/activation3Dtimecourse"

    LabelFrame $w.activ -text ""
    optionMenu2 $w.activ.menu entries($w,activ_yn)  -command "possum:updateactivation $w" 0 "None" 1 "Custom file"
    pack $w.activ.menu
    pack $w.activ -in $guivars($w,lfactivation) -side top -anchor w -padx 3 -pady 3

    frame $w.activation
	frame $w.actt2timecourse
	    FileEntry $w.act1 \
		-textvariable entries($w,act2) \
		-label "T2* time course            " \
		-title "Select" \
		-width 40 \
		-filedialog directory
	    label $w.act1spc -text "" -width 3
	    button $w.act1prev -command "possum:previewt2startimecoursedialog $w" -text "Preview T2* time course" -width 25
	pack $w.act1 $w.act1spc $w.act1prev -in $w.actt2timecourse -anchor w -side left -padx 3 -pady 3
	frame $w.spatialmod
	    FileEntry $w.act2 \
		-textvariable entries($w,act1) \
		-filetypes IMAGE \
		-label "T2* spatial modulation  " \
		-title "Select" \
		-width 40 \
		-filedialog directory \
		-command  "possum:updateACTprop $w"
	    label $w.act2spc -text "" -width 3
	    button $w.act2prev -command "possum:previewactivdialog $w" -text "Preview T2* spatial modulation" -width 25
	pack $w.act2 $w.act2spc $w.act2prev -in $w.spatialmod -anchor w -side left -padx 3 -pady 3
    pack $w.actt2timecourse $w.spatialmod -in $w.activation

    # ------ Noise--------------------
    set guivars($w,lfnoise) [$w.nb getframe noise]
    set entries($w,noise_yn) 0
    set entries($w,noisesnr) 10
    set entries($w,noisesigma) 0
    frame $w.noiseval
    frame $w.noiseval1
    frame $w.noiseval2
    LabelSpinBox $w.noiseval1.snr -label "" -width 8 \
       -textvariable entries($w,noisesnr) -range { 0.0   1000000.0  0.5 } -disabledbackground gray
    LabelSpinBox $w.noiseval2.sigma -label "" -width 8 \
       -textvariable entries($w,noisesigma) -range { 0.0   100000000.0  0.1 } -disabledbackground gray
    radiobutton $w.noiseval1.unitssnr \
       -text "SNR (relative to median object intensity): " \
       -variable entries($w,noiseunits) \
       -value snr -anchor w \
       -command "possum:updatenoiseunits $w " -width 35
    radiobutton $w.noiseval2.unitssigma \
       -text "Absolute intensity (std dev): " \
       -variable entries($w,noiseunits) \
       -value sigma -anchor w \
       -command "possum:updatenoiseunits $w " -width 35
    pack $w.noiseval1.unitssnr $w.noiseval1.snr \
       -in $w.noiseval1 -side left -anchor w -padx 3 -pady 3 
    pack $w.noiseval2.unitssigma $w.noiseval2.sigma \
       -in $w.noiseval2 -side left -anchor w -padx 3 -pady 3 
    pack $w.noiseval1 $w.noiseval2 -in $w.noiseval
    $w.noiseval1.unitssnr select
    $w.noiseval2.sigma configure -state disabled
    LabelFrame $w.noise -text ""
    optionMenu2 $w.noise.menu entries($w,noise_yn)  -command "possum:updatenoise $w" 0 "None" 1 "Thermal (white) noise "
    pack $w.noise.menu
    pack $w.noise -in $guivars($w,lfnoise) -side top -anchor w -padx 3 -pady 3

    #------Output---------------------
    set outputlf [$w.nb getframe output]
    set entries($w,out) "$PWD/simdir"
    FileEntry $w.out \
    -textvariable entries($w,out) \
    -label "Output directory  " \
    -title "Select" \
    -width 50 \
    -filedialog directory

##tejas
    #------Run---------------------
		frame $w.run
		label $w.run.procstext -text "Number of Processors:" -anchor w -justify left 
		LabelSpinBox $w.run.procsspin -label " " -textvariable entries($w,numproc) \
			-range { 1   10000  1 } \
			-command "$w.run.procsspin.spin.e validate; possum:updatecomptime $w" \
			-modifycmd "possum:updatecomptime $w"
		pack $w.run.procstext $w.run.procsspin -in $w.run -side left -anchor nw -padx 3 -pady 5
		possum:updatecomptime $w

		label $w.run.runtimetext -text "Predicted Run time:" -anchor w -justify left 
		entry $w.run.runtimeout -textvariable entries($w,comptime) -width 12 -readonlybackground white -state readonly
		pack $w.run.runtimetext $w.run.runtimeout -in $w.run -side left -anchor nw -padx 3 -pady 5

		# Collapsible frame for advanced options
		collapsible frame $w.runadv -title "Advanced" -command "$w.nb compute_size; set dummy"
			frame $w.runopts
				label $w.runopts.segtext -text "Segment size:" -anchor w -justify left 
				LabelSpinBox $w.runopts.segspin -label " " -width 6 -textvariable entries($w,segs) \
					-range { 1   1000000  1 } \
					-command "$w.runopts.segspin.spin.e validate; possum:updatecomptime $w" \
					-modifycmd "possum:updatecomptime $w"
#				pack $w.runopts.segtext $w.runopts.segspin -in $w.runadv.b -side left -anchor w -padx 3 -pady 5

				label $w.runopts.runtimeusrtext -text "Run time:" -anchor w -justify left 
				LabelSpinBox $w.runopts.runtimeusrspin -label " " \
					-width 4 -textvariable entries($w,ctt) \
					-range { 0   1000000  1 } \
					-command "$w.runopts.runtimeusrspin.spin.e validate"
			pack $w.runopts.segtext $w.runopts.segspin $w.runopts.runtimeusrtext $w.runopts.runtimeusrspin -in $w.runopts -anchor w -side left -anchor w -padx 3 -pady 5

#			frame $w.theorygo
#				button $w.theorygo.generate -command "Possum:generateTheoryVolumeWindow $w" \
#					-text "Spin-Echo Contrast Images" -width 50
#				label $w.theorygo.spc1 -text "" -width 20
#			pack $w.theorygo.spc1 $w.theorygo.generate -in $w.theorygo -pady 5

			frame $w.savescript
				label $w.savescript.status -text "" -width 50 -foreground black
				button $w.savescript.go -text "Save commands to script" -width 50 \
					-command "possum:savescriptdialog $w"
			pack $w.savescript.go -in $w.savescript -pady 5 -anchor w -side top
#		pack $w.runopts $w.theorygo $w.savescript -in $w.runadv.b -anchor w
		pack $w.runopts $w.savescript -in $w.runadv.b -anchor w

    pack $w.out $w.run $w.runadv -in $outputlf -side top -anchor w -pady 3 -padx 5

##tejas-end    

    # Outside the nb part.
    # ---- Pack all of the options ----
    frame $w.f.opts
    pack $w.nb -in $w.f.opts -side top
#    pack $w.np $w.ct $w.advanced -in $w.f.opts  -side left -padx 2
    pack $w.f.opts -in $w.f -side left -padx 8 -pady 6 -expand yes -fill both
   
    # ---- Button Frame ----
    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    button $w.btns.go     -command "Possum:apply $w" \
	    -text "Go" -width 5
    button $w.btns.makedir     -command "Possum:makedir $w" \
	    -text "Make Directory" -width 10
    button $w.btns.cancel    -command "destroy $w" \
	    -text "Exit" -width 5
    button $w.btns.save -command "feat_file:setup_dialog $w a a a [namespace current] *.fsf {Save Possum setup} {possum:write $w} {}" -text "Save"

    button $w.btns.load -command "feat_file:setup_dialog $w a a a [namespace current] *.fsf {Load Possum setup} {possum:load $w} {}" -text "Load"
    button $w.btns.help -command "FmribWebHelp file: ${FSLDIR}/doc/possum/index.html" \
            -text "Help" -width 5
    pack $w.btns.b -side bottom -fill x
    pack $w.btns.go $w.btns.makedir $w.btns.save $w.btns.load $w.btns.cancel $w.btns.help -in $w.btns.b \
	    -side left -expand yes -padx 3 -pady 10 -fill y
    pack $w.f $w.btns -expand yes -fill both
}

proc Possum:pulsecheck { w {onlyCheck 0} } {
    global entries
    set status [ possum:pulsecheck $w $entries($w,obvol) $entries($w,mrpar) $entries($w,seqtype) $entries($w,te) $entries($w,tr) $entries($w,trslc) $entries($w,outsize_nx) $entries($w,outsize_ny) $entries($w,outsize_nz) $entries($w,outsize_dx) $entries($w,outsize_dy) $entries($w,outsize_dz) $entries($w,fov_x)  $entries($w,fov_y)  $entries($w,fov_z)  $entries($w,numvol) $entries($w,zstart) $entries($w,gap) $entries($w,bw) $entries($w,readgrad) $entries($w,phencode) $entries($w,slcselect) $entries($w,plus) $entries($w,maxG)  $entries($w,riseT) $entries($w,b0f) $entries($w,mot)  $entries($w,act1) $entries($w,act2) $entries($w,out) $entries($w,numproc) $entries($w,segs) $entries($w,slcprof) $entries($w,cover) $entries($w,flipangle) $onlyCheck]
    update idletasks
}

proc Possum:apply { w } {
    global entries FSLDIR
    set $entries($w,onlyCheck) 0
    # start by saving the fsf file (with all variables as they are now)
    #if { ! [ file isdirectory $entries($w,out) ] } { 
	#catch { exec sh -c "mkdir $entries($w,out)" } oval
    #}
   
    if { $entries($w,obvol) == "" } {
       puts "The input object not specified."
     return
    } 
    if { $entries($w,mrpar) == "" } {
       puts "The input MR parameters not specified."
       return
    }
    if { $entries($w,slcprof) == "" } {
       puts "The slice profile not specified."
     return
    }
    if { $entries($w,b0inhtime_yn) == 1 && $entries($w,motion_yn) == 1 } {
       puts "Warning: At the moment B0 field changing in time can not be simulated while the object is moving. This will be implemented into POSSUM at a later stage."
       return
    }

    # Custom pulse : Check if pulse sequnce files exist
    if { $entries($w,custompulse_yn) == 1 } {
	foreach i [list "$entries($w,cuspulse)" "$entries($w,cuspulse).info" "$entries($w,cuspulse).posx" "$entries($w,cuspulse).posy" "$entries($w,cuspulse).posz" ] {
		if { [file exists $i] == 0 } {
			puts "Error! Custom pulse file : '$i' not found!"
			return
		}
	}
    }

    # checks if the object is the same size as the b0file 
   if { $entries($w,b0inh_yn) == 1 } {
    if { $entries($w,vcX) != $entries($w,vcXb0) || $entries($w,vcY) != $entries($w,vcYb0) ||  $entries($w,vcZ) != $entries($w,vcZb0) ||  $entries($w,inNx) != $entries($w,inNxb0) ||  $entries($w,inNy) != $entries($w,inNyb0) ||  $entries($w,inNz) != $entries($w,inNzb0) } {
     puts "The object and the B0 file do not match in dimension or voxel size."
     puts "Object dim: $entries($w,inNx), $entries($w,inNy), $entries($w,inNz)"
     puts "B0 dim: $entries($w,inNxb0), $entries($w,inNyb0), $entries($w,inNzb0)"
     puts "Object voxsize: $entries($w,vcX), $entries($w,vcY), $entries($w,vcZ)"
     puts "B0 voxsize: $entries($w,vcXb0), $entries($w,vcYb0), $entries($w,vcZb0)"
     return
    }
   }
    # checks if the object is the same size as the b0time file
    if { $entries($w,b0inhtime_yn) == 1 } {
    if { $entries($w,vcX) != $entries($w,vcXb0time) || $entries($w,vcY) != $entries($w,vcYb0time) ||  $entries($w,vcZ) != $entries($w,vcZb0time) ||  $entries($w,inNx) != $entries($w,inNxb0time) ||  $entries($w,inNy) != $entries($w,inNyb0time) ||  $entries($w,inNz) != $entries($w,inNzb0time) } {
     puts "The object and the B0 file (time changing) do not match in dimension or voxel size."
     puts "Object dim: $entries($w,inNx), $entries($w,inNy), $entries($w,inNz)"
     puts "B0 dim: $entries($w,inNxb0time), $entries($w,inNyb0time), $entries($w,inNzb0time)"
     puts "Object voxsize: $entries($w,vcX), $entries($w,vcY), $entries($w,vcZ)"
     puts "B0 voxsize: $entries($w,vcXb0time), $entries($w,vcYb0time), $entries($w,vcZb0time)"
     return
    }
   }
   if { $entries($w,activ_yn) == 1 } {
    if { $entries($w,vcX) != $entries($w,vcXact) || $entries($w,vcY) != $entries($w,vcYact) ||  $entries($w,vcZ) != $entries($w,vcZact) ||  $entries($w,inNx) != $entries($w,inNxact) ||  $entries($w,inNy) != $entries($w,inNyact) ||  $entries($w,inNz) != $entries($w,inNzact) } {
     puts "The object and the activation file do not match in dimension or voxel size."
     puts "Object dim: $entries($w,inNx), $entries($w,inNy), $entries($w,inNz)"
     puts "Activation dim: $entries($w,inNxact), $entries($w,inNyact), $entries($w,inNzact)"
     puts "Object voxsize: $entries($w,vcX), $entries($w,vcY), $entries($w,vcZ)"
     puts "Activation voxsize: $entries($w,vcXact), $entries($w,vcYact), $entries($w,vcZact)"
     return
    }
   }
   #checks if the output voxel size is smaller than the input voxel size
   if { $entries($w,vcX) > $entries($w,outsize_dx) || $entries($w,vcY) > $entries($w,outsize_dy) ||  $entries($w,vcZ) > $entries($w,outsize_dz)} {
     puts "The input object voxel size (every direction) should not be bigger than the output image voxel size."
     puts "The input object voxel size: $entries($w,vcX), $entries($w,vcY), $entries($w,vcZ) "
     puts "The output image voxel size: $entries($w,outsize_dx), $entries($w,outsize_dy), $entries($w,outsize_dz)"
     return
    }
   #check if pulse is properly set up
   Possum:pulsecheck $w
   if { $entries($w,pulsechecktest) == 0 } {
	puts "Error in pulse sequence!"
	return
   }

    puts "Creating the POSSUM directory..."
    puts ""
    if { $entries($w,out)  == "" } {
       puts "The output directory not specified."
     exit
    } else { 
	new_file $entries($w,out)
	catch { exec sh -c "mkdir $entries($w,out)" } oval
    }
    possum:write $w $entries($w,out)/possum.fsf

    # now do some logic to figure out the parameters to pass on
    if { $entries($w,b0inh_yn) == 0 } { 
	set b0file "" 
    } else {
	set b0file $entries($w,b0f)
    }
    if { $entries($w,b0inhtime_yn) == 0 } { 
	set b0filetime "" 
	set b0filetimecourse "" 
    } else {
	set b0filetime $entries($w,b0ftime)
	set b0filetimecourse $entries($w,b0ftimecourse)
    }
    if { $entries($w,motion_yn) == 0 } { 
	set motfile "${FSLDIR}/data/possum/zeromotion" 
    } else {
	set motfile $entries($w,mot)
    }
    if { $entries($w,activ_yn) == 0 } { 
	set act1file "" 
	set act2file "" 
    } else {
	set act1file $entries($w,act1)
	set act2file $entries($w,act2)
    }
    set filename "$entries($w,out)/noise"
    set log [open "$filename" w]
    if { $entries($w,noiseunits) == "snr" && $entries($w,noise_yn) == 1 } { 
	puts $log "snr $entries($w,noisesnr) "
    } else {
	puts $log "sigma $entries($w,noisesigma) "
    }
    close $log
    set status [ possum:proc $w $entries($w,proctime) $entries($w,obvol) $entries($w,mrpar) $entries($w,seqtype) $entries($w,te) $entries($w,tr) $entries($w,trslc) $entries($w,outsize_nx) $entries($w,outsize_ny) $entries($w,outsize_nz) $entries($w,outsize_dx) $entries($w,outsize_dy) $entries($w,outsize_dz) $entries($w,fov_x)  $entries($w,fov_y)  $entries($w,fov_z)  $entries($w,numvol) $entries($w,zstart) $entries($w,gap) $entries($w,bw) $entries($w,readgrad) $entries($w,phencode) $entries($w,slcselect) $entries($w,plus) $entries($w,maxG)  $entries($w,riseT) $b0file $entries($w,b0fieldstrength) $entries($w,b0units)  $b0filetime $b0filetimecourse $entries($w,b0unitstime) $motfile $act1file $act2file $entries($w,out) $entries($w,numproc) $entries($w,segs) $entries($w,slcprof) $entries($w,cover) $entries($w,flipangle) $entries($w,slcsampfactor) ]
    update idletasks
    puts "Job submitted."
    puts ""
    puts "You can follow the POSSUM process by looking at the possum.log file."
    puts ""
    puts "If you want to see the individual processes see the logs directory."
}

proc Possum:makedir { w } {
    global entries FSLDIR

    # start by saving the fsf file (with all variables as they are now)
    #if { ! [ file isdirectory $entries($w,out) ] } { 
	#catch { exec sh -c "mkdir $entries($w,out)" } oval
    #}
    if { $entries($w,obvol) == "" } {
       puts "The input object not specified."
     return
    } 
    if { $entries($w,mrpar) == "" } {
       puts "The input MR parameters not specified."
       return
    }
    if { $entries($w,slcprof) == "" } {
       puts "The slice profile not specified."
     return
    }
    # Custom pulse : Check if pulse sequnce files exist
    if { $entries($w,custompulse_yn) == 1 } {
	foreach i [list "$entries($w,cuspulse)" "$entries($w,cuspulse).info" "$entries($w,cuspulse).posx" "$entries($w,cuspulse).posy" "$entries($w,cuspulse).posz" ] {
		if { [file exists $i] == 0 } {
			puts "Error! Custom pulse file : '$i' not found!"
			return
		}
	}
    }
    if { $entries($w,b0inhtime_yn) == 1 && $entries($w,motion_yn) == 1 } {
       puts "Warning: At the moment B0 field changing in time can not be simulated while the object is moving. This will be implemented into POSSUM at a later stage."
       return
    }
    # checks if the object is the same size as the b0file 
   if { $entries($w,b0inh_yn) == 1 } {
    if { $entries($w,vcX) != $entries($w,vcXb0) || $entries($w,vcY) != $entries($w,vcYb0) ||  $entries($w,vcZ) != $entries($w,vcZb0) ||  $entries($w,inNx) != $entries($w,inNxb0) ||  $entries($w,inNy) != $entries($w,inNyb0) ||  $entries($w,inNz) != $entries($w,inNzb0) } {
     puts "The object and the B0 file do not match in dimension or voxel size."
     puts "Object dim: $entries($w,inNx), $entries($w,inNy), $entries($w,inNz)"
     puts "B0 dim: $entries($w,inNxb0), $entries($w,inNyb0), $entries($w,inNzb0)"
     puts "Object voxsize: $entries($w,vcX), $entries($w,vcY), $entries($w,vcZ)"
     puts "B0 voxsize: $entries($w,vcXb0), $entries($w,vcYb0), $entries($w,vcZb0)"
     return
    }
   }
    # checks if the object is the same size as the b0time file
    if { $entries($w,b0inhtime_yn) == 1 } {
    if { $entries($w,vcX) != $entries($w,vcXb0time) || $entries($w,vcY) != $entries($w,vcYb0time) ||  $entries($w,vcZ) != $entries($w,vcZb0time) ||  $entries($w,inNx) != $entries($w,inNxb0time) ||  $entries($w,inNy) != $entries($w,inNyb0time) ||  $entries($w,inNz) != $entries($w,inNzb0time) } {
     puts "The object and the B0 file (time changing) do not match in dimension or voxel size."
     puts "Object dim: $entries($w,inNx), $entries($w,inNy), $entries($w,inNz)"
     puts "B0 dim: $entries($w,inNxb0time), $entries($w,inNyb0time), $entries($w,inNzb0time)"
     puts "Object voxsize: $entries($w,vcX), $entries($w,vcY), $entries($w,vcZ)"
     puts "B0 voxsize: $entries($w,vcXb0time), $entries($w,vcYb0time), $entries($w,vcZb0time)"
     return
    }
   }
   if { $entries($w,activ_yn) == 1 } {
    if { $entries($w,vcX) != $entries($w,vcXact) || $entries($w,vcY) != $entries($w,vcYact) ||  $entries($w,vcZ) != $entries($w,vcZact) ||  $entries($w,inNx) != $entries($w,inNxact) ||  $entries($w,inNy) != $entries($w,inNyact) ||  $entries($w,inNz) != $entries($w,inNzact) } {
     puts "The object and the activation file do not match in dimension or voxel size."
     puts "Object dim: $entries($w,inNx), $entries($w,inNy), $entries($w,inNz)"
     puts "Activation dim: $entries($w,inNxact), $entries($w,inNyact), $entries($w,inNzact)"
     puts "Object voxsize: $entries($w,vcX), $entries($w,vcY), $entries($w,vcZ)"
     puts "Activation voxsize: $entries($w,vcXact), $entries($w,vcYact), $entries($w,vcZact)"
     return
    }
   }
   #checks if the output voxel size is smaller than the input voxel size
   if { $entries($w,vcX) > $entries($w,outsize_dx) || $entries($w,vcY) > $entries($w,outsize_dy) ||  $entries($w,vcZ) > $entries($w,outsize_dz)} {
     puts "The input object voxel size (every direction) should not be bigger than the output image voxel size."
     puts "The input object voxel size: $entries($w,vcX), $entries($w,vcY), $entries($w,vcZ) "
puts "The output image voxel size: $entries($w,outsize_dx), $entries($w,outsize_dy), $entries($w,outsize_dz)"
     return
    }
   #check if pulse is properly set up
   Possum:pulsecheck $w
   if { $entries($w,pulsechecktest) == 0 } {
	puts "Error in pulse sequence!"
	return
   }
    if { $entries($w,out)  == "" } {
       puts "The output directory not specified."
     exit
    } else { 
	new_file $entries($w,out)
	catch { exec sh -c "mkdir $entries($w,out)" } oval
    }
    possum:write $w $entries($w,out)/possum.fsf

    # now do some logic to figure out the parameters to pass on
    if { $entries($w,b0inh_yn) == 0 } { 
	set b0file "" 
    } else {
	set b0file $entries($w,b0f)
    }
    if { $entries($w,b0inhtime_yn) == 0 } { 
	set b0filetime "" 
	set b0filetimecourse "" 
    } else {
	set b0filetime $entries($w,b0ftime)
	set b0filetimecourse $entries($w,b0ftimecourse)
    }
    if { $entries($w,motion_yn) == 0 } { 
	set motfile "${FSLDIR}/data/possum/zeromotion" 
    } else {
	set motfile $entries($w,mot)
    }
    if { $entries($w,activ_yn) == 0 } { 
	set act1file "" 
	set act2file "" 
    } else {
	set act1file $entries($w,act1)
	set act2file $entries($w,act2)
    }
    set filename "$entries($w,out)/noise"
    set log [open "$filename" w]
    if { $entries($w,noiseunits) == "snr" && $entries($w,noise_yn) == 1 } { 
	puts $log "snr $entries($w,noisesnr) "
    } else {
	puts $log "sigma $entries($w,noisesigma) "
    }
    close $log
    set status [ possum:procmakedir $w $entries($w,proctime) $entries($w,obvol) $entries($w,mrpar) $entries($w,seqtype) $entries($w,te) $entries($w,tr) $entries($w,trslc) $entries($w,outsize_nx) $entries($w,outsize_ny) $entries($w,outsize_nz) $entries($w,outsize_dx) $entries($w,outsize_dy) $entries($w,outsize_dz) $entries($w,fov_x)  $entries($w,fov_y)  $entries($w,fov_z)  $entries($w,numvol) $entries($w,zstart) $entries($w,gap) $entries($w,bw) $entries($w,readgrad) $entries($w,phencode) $entries($w,slcselect) $entries($w,plus) $entries($w,maxG)  $entries($w,riseT) $b0file $entries($w,b0fieldstrength) $entries($w,b0units)  $b0filetime $b0filetimecourse $entries($w,b0unitstime) $motfile $act1file $act2file $entries($w,out) $entries($w,numproc) $entries($w,segs) $entries($w,slcprof) $entries($w,cover) $entries($w,flipangle) $entries($w,slcsampfactor)]
    update idletasks
   
}

proc possum:previewimage { w } {
    global entries FSLDIR
    set convertcom "${FSLDIR}/bin/pngappend"
    set filenm $entries($w,obvol)
    set validim 0
    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "${FSLDIR}/bin/slicer $filenm -a ${tmpnam}.png" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png ${tmpnam}" } oval
    }
    if { $validim == 1 } {
	$w.objim.label configure -image $graphpic
    } else { 
	$w.objim.label configure -image "" -text "Could not generate preview"
    }
    $w.nb compute_size
}

proc possum:previewimagedialog { w } {
    global entries FSLDIR
    set count 0
    set w1 ".dialog[incr count]"
    while { [ winfo exists $w1 ] } {
        set w1 ".dialog[incr count]"
    }
    toplevel $w1
    wm title $w1 "Preview of '[exec sh -c "basename $entries($w,obvol)"]'"
    wm iconname $w1 "SlicePreview"
    wm iconbitmap $w1 @${FSLDIR}/tcl/fmrib.xbm
    frame $w1.sprev
    label $w1.sprev.label -image "" -text "\n    Generating preview ... please wait    \n\n"
    pack $w1.sprev.label -in $w1.sprev
    pack $w1.sprev -in $w1
    update idletasks

    set convertcom "${FSLDIR}/bin/pngappend"
    set filenm $entries($w,obvol)
    set validim 0
    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "${FSLDIR}/bin/slicer $filenm -a ${tmpnam}.png" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png ${tmpnam}" } oval
    }
    if { $validim == 1 } {
	$w1.sprev.label configure -image $graphpic
    } else { 
	$w1.sprev.label configure -image "" -text "Could not generate preview" -width 50 -height 5 -foreground red
    }
}

proc possum:previewimage_steve { w } {
    global entries FSLDIR
    set convertcom "${FSLDIR}/bin/pngappend"
    set filenm $entries($w,obvol)
    set validim 0
    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam /tmp/possum" } tmpnam
	if { [ exec sh -c "${FSLDIR}/bin/fslnvols $filenm 2> /dev/null" ] == 3 } {
	    catch { exec sh -c "${FSLDIR}/bin/fslsplit $filenm $tmpnam" } oval
	    catch { exec sh -c "${FSLDIR}/bin/overlay 0 0 ${tmpnam}0001 0 1 ${tmpnam}0000 0.5 1 ${tmpnam}0002 0.5 1 ${tmpnam}out" } oval
	    set filenm ${tmpnam}out
	}
	catch { exec sh -c "${FSLDIR}/bin/slicer $filenm -a ${tmpnam}.png" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png ${tmpnam}" } oval
    }
    if { $validim == 1 } {
	$w.objim.label configure -image $graphpic
    } else { 
	$w.objim.label configure -image "" -text "Could not generate preview"
    }
    $w.nb compute_size
}

proc possum:updateFOV {w} {
    global entries
    set entries($w,fov_x) [ expr $entries($w,outsize_nx) *$entries($w,outsize_dx) ]
    set entries($w,fov_y) [ expr $entries($w,outsize_ny) *$entries($w,outsize_dy) ]
    set entries($w,fov_z) [ expr $entries($w,outsize_nz) *$entries($w,outsize_dz) ]
}

proc possum:updateVSIZE {w} {
    global entries
    set entries($w,outsize_dx) [ expr $entries($w,fov_x) * 1.0 / $entries($w,outsize_nx) ]
    set entries($w,outsize_dy) [ expr $entries($w,fov_y) * 1.0 / $entries($w,outsize_ny) ]
    set entries($w,outsize_dz) [ expr $entries($w,fov_z) * 1.0 / $entries($w,outsize_nz) ]
}

proc possum:updateTRSLC {w} {
    global entries guivars FSLDIR
#tejas-edit
    if { $entries($w,autotrslc) == 1 && $entries($w,seqtype) == "epi"} {
	set tmp [ expr $entries($w,tr)*1.0/$entries($w,outsize_nz) ]
        set entries($w,trslc) [ possum:modTRslice $tmp ]  
    }
    if { $entries($w,seqtype) != "epi" } {
	set entries($w,trslc) 0
#	pack forget $w.t.trs
	$w.trs.z configure -state disabled
    }
}

proc possum:buttonTRSLC {w} {
    global entries guivars FSLDIR
    if { $entries($w,autotrslc) == 1 } {
	$w.trs.z configure -state disabled
    } else {
	$w.trs.z configure -state normal
    }
    possum:updateTRSLC $w
}

proc possum:updateTRSLC2 {w} {
    global entries guivars FSLDIR
    if {$entries($w,seqtype) == "epi"} {
#	pack $w.t.trs -in $w.t -side left -anchor w -padx 3 -pady 3
	$w.trs.z configure -state normal
    }
}
#tejas-end

proc possum:updateb0field { w } {
    global entries guivars FSLDIR
    pack forget $w.b0test.b0spin
    if { $entries($w,b0strength) == 2 } {
        pack $w.b0test.b0spin -in $w.b0test  -side left -anchor w -padx 3 -pady 3
    } 
    if { $entries($w,b0strength) == 1 } {
	set entries($w,b0fieldstrength) 3.0
    } 
    if { $entries($w,b0strength) == 0 } {
	set entries($w,b0fieldstrength) 1.5
    }
}

proc possum:updateMRpar { w } {
    global entries guivars FSLDIR
    if { $entries($w,b0strength) == 2 } {
	set entries($w,mrpar) ""
    } 
    if { $entries($w,b0strength) == 1 } {
	set entries($w,mrpar) "${FSLDIR}/data/possum/MRpar_3T"
    } 
    if { $entries($w,b0strength) == 0 } {
	set entries($w,mrpar) "${FSLDIR}/data/possum/MRpar_1.5T"
    }
}

proc possum:updateb0fieldinh { w } {
    global entries guivars FSLDIR
    if { $entries($w,b0inh_yn)  == 1} {
	pack $w.b0mainopts -in $w.b0main -side top -anchor w -padx 3 -pady 3
        pack $w.b0main.b0u -in $w.b0main -anchor w -padx 3 -pady 3
    } else {
        pack forget $w.b0main.b0u
	pack forget $w.b0mainopts
    }
}

proc possum:updateb0fieldinhtime { w } {
    global entries guivars FSLDIR
    if { $entries($w,b0inhtime_yn)} {
	pack $w.b0timecourse -in $w.b0maintime -anchor w -padx 3 -pady 3
	pack $w.b0timeimage -in $w.b0maintime -anchor w -padx 3 -pady 3
	pack $w.b0maintime.b0u $w.b0maintime.b0im -in $w.b0maintime -anchor w -padx 3 -pady 3
	if { $entries($w,motion_yn) == 1 } {
		$w.b0maintime.b0im.b0timewarn configure \
		-text "Note:\tAt the moment B0 field changing in time can not be simulated while the object is moving.\n\tThis will be implemented into POSSUM at a later stage."
		pack $w.b0maintime.b0im.b0timewarn -in $w.b0maintime.b0im

		$w.motwarn configure \
		-text "Note:\tAt the moment B0 field changing in time can not be simulated while the object is moving.\n\tThis will be implemented into POSSUM at a later stage."
		pack $w.motwarn -in $guivars($w,lfmotion) -anchor w -side top -padx 3 -pady 3		
		
	}
    } else {
        pack forget $w.b0maintime.b0u
        pack forget $w.b0maintime.b0im
        pack forget $w.b0timecourse
	pack forget $w.b0timeimage
	$w.motwarn configure -text ""
    }
}

proc possum:updateBASEname { w } {
    global entries guivars FSLDIR
    if { $entries($w,b0inh_yn)} {
	set entries($w,b0f) "${FSLDIR}/data/possum/b0_ppm.nii.gz"
	set entries($w,b0units) "ppm"
        possum:updateB0prop $w $entries($w,b0f)
    } else {
	set entries($w,b0f) ""
    }
}

proc possum:updateBASEnametime { w } {
    global entries guivars FSLDIR
    if { $entries($w,b0inhtime_yn)} {
	set entries($w,b0ftime) "${FSLDIR}/data/possum/b0extra.nii.gz"
	set entries($w,b0ftimecourse) "${FSLDIR}/data/possum/b0timecourse"
        possum:updateB0timeprop $w $entries($w,b0ftime)
    } else {
	set entries($w,b0ftime) ""
	set entries($w,b0ftimecourse) ""
    }
}

proc possum:previewb0 { w } {
    global entries FSLDIR
    set convertcom "${FSLDIR}/bin/pngappend"
    set filenm "$entries($w,b0f)"
    set validim 0
    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "${FSLDIR}/bin/slicer $filenm -a ${tmpnam}.png" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png ${tmpnam}" } oval
    }
    if { $validim == 1 } {
	$w.b0main.b0im.label configure -image $graphpic
    } else { 
	$w.b0main.b0im.label configure -image "" -text "Could not generate preview"
    }
    $w.nb compute_size
}

proc possum:previewb0dialog { w } {
    global entries FSLDIR
    set count 0
    set w1 ".dialog[incr count]"
    while { [ winfo exists $w1 ] } {
        set w1 ".dialog[incr count]"
    }
    toplevel $w1
    wm title $w1 "Preview of '[exec sh -c "basename $entries($w,b0f)"]'"
    wm iconname $w1 "SlicePreview"
    wm iconbitmap $w1 @${FSLDIR}/tcl/fmrib.xbm
    frame $w1.sprev
    label $w1.sprev.label -image "" -text "\n    Generating preview ... please wait    \n\n"
    pack $w1.sprev.label -in $w1.sprev
    pack $w1.sprev -in $w1
    update idletasks

    set convertcom "${FSLDIR}/bin/pngappend"
    set filenm "$entries($w,b0f)"
    set validim 0
    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "${FSLDIR}/bin/slicer $filenm -a ${tmpnam}.png" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png ${tmpnam}" } oval
    }
    if { $validim == 1 } then {
	$w1.sprev.label configure -image $graphpic
    } else { 
	$w1.sprev.label configure -image "" -text "Could not generate preview" -height 5 -width 50 -foreground red
    }

}

proc possum:previewb0time { w } {
    global entries FSLDIR
    set convertcom "${FSLDIR}/bin/pngappend"
    set filenm "$entries($w,b0ftime)"
    set validim 0
    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "${FSLDIR}/bin/slicer $filenm -a ${tmpnam}.png" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png ${tmpnam}" } oval
    }
    if { $validim == 1 } {
	$w.b0maintime.b0im.label configure -image $graphpic
    } else { 
	$w.b0maintime.b0im.label configure -image "" -text "Could not generate preview"
    }
    $w.nb compute_size
}

#tejas-add
proc possum:previewb0timedialog { w } {
    global entries FSLDIR
    set count 0
    set w1 ".dialog[incr count]"
    while { [ winfo exists $w1 ] } {
        set w1 ".dialog[incr count]"
    }
    toplevel $w1
    wm title $w1 "Preview of '[exec sh -c "basename $entries($w,b0ftime)"]'"
    wm iconname $w1 "SlicePreview"
    wm iconbitmap $w1 @${FSLDIR}/tcl/fmrib.xbm
    frame $w1.sprev
    label $w1.sprev.label -image "" -text "\n    Generating preview ... please wait    \n\n"
    pack $w1.sprev.label -in $w1.sprev
    pack $w1.sprev -in $w1
    update idletasks

    set convertcom "${FSLDIR}/bin/pngappend"
    set filenm "$entries($w,b0ftime)"
    set validim 0
    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "${FSLDIR}/bin/slicer $filenm -a ${tmpnam}.png" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png ${tmpnam}" } oval
    }
    if { $validim == 1 } {
	$w1.sprev.label configure -image $graphpic
    } else { 
	$w1.sprev.label configure -image "" -text "Could not generate preview" -height 5 -width 50 -foreground red
    }
}

proc possum:previewactivdialog { w } {
    global entries FSLDIR
    global entries FSLDIR
    set count 0
    set w1 ".dialog[incr count]"
    while { [ winfo exists $w1 ] } {
        set w1 ".dialog[incr count]"
    }
    toplevel $w1
    wm title $w1 "Preview of Activation Timecourse"
    wm iconname $w1 "SlicePreview"
    wm iconbitmap $w1 @${FSLDIR}/tcl/fmrib.xbm
    frame $w1.sprev
    label $w1.sprev.label -image "" -text "\n    Generating preview ... please wait    \n\n"
    pack $w1.sprev.label -in $w1.sprev
    pack $w1.sprev -in $w1
    update idletasks

    set convertcom "${FSLDIR}/bin/pngappend"
    set filenm $entries($w,act1)
    set validim 0
    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "${FSLDIR}/bin/imtest $entries($w,obvol)" } oval
	if { $oval == 1 } {
	    catch { exec sh -c "${FSLDIR}/bin/slicer $filenm $entries($w,obvol) -a ${tmpnam}.png" } oval 
	} else {
	    catch { exec sh -c "${FSLDIR}/bin/slicer $filenm -a ${tmpnam}.png" } oval 
	}
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png ${tmpnam}" } oval
    }
    if { $validim == 1 } {
	$w1.sprev.label configure -image $graphpic
    } else { 
	$w1.sprev.label configure -image "" -text "Could not generate preview" -height 5 -width 50 -foreground red
    }
}

proc possum:previewmotionplotdialog { w } {
    global entries FSLDIR
    set count 0
    set w1 ".dialog[incr count]"
    while { [ winfo exists $w1 ] } {
        set w1 ".dialog[incr count]"
    }
    toplevel $w1
    set filenm $entries($w,mot)
    wm title $w1 "Preview of '[exec sh -c "basename $filenm"]'"
    wm iconname $w1 "SlicePreview"
    wm iconbitmap $w1 @${FSLDIR}/tcl/fmrib.xbm
    frame $w1.sprev
    label $w1.sprev.label -image "" -text "\n    Generating preview ... please wait    \n\n"
    pack $w1.sprev.label -in $w1.sprev
    pack $w1.sprev -in $w1
    update idletasks

    set convertcom "${FSLDIR}/bin/pngappend"
    set validim 0
#    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    set oval [file exists $filenm]
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "${FSLDIR}/bin/fsl_tsplot -i $filenm -o ${tmpnam}.png --start=2 --finish=7 -x 'Motion matrix row' -y Magnitude -w 500" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png" } oval
    }
    if { $validim == 1 } {
	$w1.sprev.label configure -image $graphpic
    } else { 
	$w.objim.label configure -image "" -text "Could not generate preview"
    }
}

proc possum:previewmotionplot { w } {
    global entries guivars FSLDIR
    set convertcom "${FSLDIR}/bin/pngappend"
    set validim 0
    set width 100
    set filenm $entries($w,mot)
    if { [file exists $filenm] } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "cat $filenm | wc -l" } sizeoffile
#	if { $sizeoffile > 10 } { set width 700 }
	set legend [list Tx(m) Ty(m) Tz(m) Rx(rad) Ry(rad) Rz(rad) ]
	foreach i $legend {
		catch { exec sh -c "echo '$i' >> $tmpnam.legend" } forEachCatch
	}
	catch { exec sh -c "${FSLDIR}/bin/fsl_tsplot -i $filenm -o ${tmpnam}.png -l $tmpnam.legend --start=2 --finish=7 -x  'Motion matrix row' -y Magnitude -w 500" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png $tmpnam.legend" } oval
    }
    if { $validim == 1 } {
	$w.motprevimg configure -image $graphpic
	pack $w.motprevimg -in  $guivars($w,lfmotion) -side top -anchor w -padx 3 -pady 3
    } else { 
	$w.motprevimg configure -image "" -text "Could not generate preview" -height 5 -width 50 -foreground red
	pack $w.motprevimg -in  $guivars($w,lfmotion) -side top -anchor w -padx 3 -pady 3
    }
    $w.nb compute_size
}

proc possum:slcprofpreviewdialog { w } {
    global entries FSLDIR
    set count 0
    set w1 ".dialog[incr count]"
    while { [ winfo exists $w1 ] } {
        set w1 ".dialog[incr count]"
    }
    toplevel $w1
    set filenm $entries($w,slcprof)
    wm title $w1 "Preview of '[exec sh -c "basename $filenm"]'"
    wm iconname $w1 "SlicePreview"
    wm iconbitmap $w1 @${FSLDIR}/tcl/fmrib.xbm
    frame $w1.sprev
    label $w1.sprev.label -image "" -text "\n    Generating preview ... please wait    \n\n"
    pack $w1.sprev.label -in $w1.sprev
    pack $w1.sprev -in $w1
    update idletasks

    set convertcom "${FSLDIR}/bin/pngappend"
    set validim 0
#    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    set imtest [file exists $filenm]
    if { $imtest == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "${FSLDIR}/bin/fsl_tsplot -i $filenm -o ${tmpnam}.png --start=2 --finish=2 -x Steps -y Magnitude -w 500" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png" } oval
    }
    if { $validim == 1 } {
	$w1.sprev.label configure -image $graphpic
    } else { 
	$w1.sprev.label configure -image "" -text "Could not generate preview" -width 50 -height 5 -foreground red
    }
}

proc possum:slcprofpreview { w } {
    global entries FSLDIR
    set convertcom "${FSLDIR}/bin/pngappend"
    set validim 0
    set width 100
    set filenm $entries($w,slcprof)
    if { [file exists $filenm] } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "cat $filenm | wc -l" } sizeoffile
#	if { $sizeoffile > 10 } { set width 700 }
	catch { exec sh -c "${FSLDIR}/bin/fsl_tsplot -i $filenm -o ${tmpnam}.png -l $FSLDIR/tcl/motlegend --start=2 --finish=7 -x Steps -y Magnitude -w 500" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png" } oval
    }
    if { $validim == 1 } {
	pack $w.slcprof.prev -in $w.slcprof -anchor e -side left -pady 3
	$w.slcprof.prev configure -image $graphpic
    } else {
	pack $w.slcprof.prev -in $w.slcprof -anchor e -side left -pady 3
	$w.slcprof.prev configure -image "" -text "Could not generate preview"
    }
    $w.nb compute_size
}

proc possum:previewb0timecoursedialog { w } {
	global entries FSLDIR POSSUMDIR
	set count 0
	set w1 ".dialog[incr count]"
	while { [ winfo exists $w1 ] } {
		set w1 ".dialog[incr count]"
	}
	toplevel $w1
	set filenm $entries($w,b0ftimecourse)
	wm title $w1 "Preview of '[exec sh -c "basename $filenm"]'"
	wm iconname $w1 "SlicePreview"
	wm iconbitmap $w1 @${FSLDIR}/tcl/fmrib.xbm
	frame $w1.sprev
	label $w1.sprev.label -image "" -text "\n    Generating preview ... please wait    \n\n"
	pack $w1.sprev.label -in $w1.sprev
	pack $w1.sprev -in $w1
	update idletasks

	set convertcom "${FSLDIR}/bin/pngappend"
	set validim 0
	if { [file exists $filenm] == 0 } { set oval 0 } else { set oval 1 }
	if { $oval == 1 } {
		catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
		catch { exec sh -c "${FSLDIR}/bin/fsl_tsplot -i $filenm -o ${tmpnam}.png --start=2 --finish=2 -x Steps -y Magnitude -w 1000" } oval
		catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
		if { [ file exists ${tmpnam}.gif ] } {
			set graphpic [image create photo -file ${tmpnam}.gif ]
			set validim 1
		}
		catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png" } oval
	}
	if { $validim == 1 } {
		$w1.sprev.label configure -image $graphpic
	} else { 
		$w1.sprev.label configure -image "" -text "Could not generate preview" -height 5 -width 50 -foreground red
	}
}

proc possum:previewt2startimecoursedialog { w } {
	global entries FSLDIR POSSUMDIR
	set count 0
	set w1 ".dialog[incr count]"
	while { [ winfo exists $w1 ] } {
		set w1 ".dialog[incr count]"
	}
	toplevel $w1
	set filenm $entries($w,act2)
	wm title $w1 "Preview of '[exec sh -c "basename $filenm"]'"
	wm iconname $w1 "SlicePreview"
	wm iconbitmap $w1 @${FSLDIR}/tcl/fmrib.xbm
	frame $w1.sprev
	label $w1.sprev.label -image "" -text "\n    Generating preview ... please wait    \n\n"
	pack $w1.sprev.label -in $w1.sprev
	pack $w1.sprev -in $w1
	update idletasks

	set convertcom "${FSLDIR}/bin/pngappend"
	set validim 0
	if { [file exists $filenm] == 0 } { set oval 0 } else { set oval 1 }
	if { $oval == 1 } {
		catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
		catch { exec sh -c "${FSLDIR}/bin/fsl_tsplot -i $filenm -o ${tmpnam}.png --start=2 --finish=2 -x Steps -y Magnitude -w 500" } oval
		catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
		if { [ file exists ${tmpnam}.gif ] } {
			set graphpic [image create photo -file ${tmpnam}.gif ]
			set validim 1
		}
		catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png" } oval
	}
	if { $validim == 1 } {
		$w1.sprev.label configure -image $graphpic
	} else { 
		$w1.sprev.label configure -image "" -text "Could not generate preview" -height 5 -width 50 -foreground red
	}
}

###	Add a new window for T1 weighted images
proc Possum:generateTheoryVolumeWindow { w } {
    global entries FSLDIR PWD
    package require BWLabelSpinBox
    set count 0
    set w1 ".dialog[incr count]"
    while { [ winfo exists $w1 ] } {
        set w1 ".dialog[incr count]"
    }
    toplevel $w1
    wm title $w1 "Theory Volume"
    wm iconname $w1 "Possum"
    wm iconbitmap $w1 @${FSLDIR}/tcl/fmrib.xbm
##    frame $w1.f
    NoteBook $w1.nb -side top -bd 2 -tabpady {5 5} -arcradius 3
    $w1.nb insert 0 object -text "Object"
    $w1.nb insert 1 imdims -text "Image Dimensions"
    $w1.nb insert 2 mrpar -text "MR Parameters"
    $w1.nb insert 3 noise -text "Noise"
    $w1.nb insert 4 run -text "Output"
    $w1.nb raise object

    set entries($w1,obvol) $entries($w,obvol)
    set entries($w1,mrpar) "${FSLDIR}/data/possum/MRpar_1.5T"
    set entries($w1,theoryout) "$PWD/SEvolume.nii.gz"
    set entries($w1,mrparcontents) ""

    frame $w1.opts
	    #---input object-----#
	    set objlf [$w1.nb getframe object]
	    frame $w1.object
		FileEntry $w1.object.dir \
		-textvariable entries($w1,obvol) \
		-filetypes IMAGE \
		-label "Input Object      " \
		-title "Select" \
		-width 40 \
		-filedialog directory \
		-command "$w1.status.lab configure -text '' -foreground black;possum:updatetheoryobj $w1"
	    pack $w1.object.dir -in $w1.object -side top -anchor w -pady 3 -padx 5
	    frame $w1.objim
		button $w1.objim.preview -text "Preview Image" -command "possum:previewimagedialog $w1"
	    pack $w1.objim.preview -in $w1.objim -pady 10
	    pack $w1.object $w1.objim -in $objlf -side top -anchor w -pady 3 -padx 5

	    #----image dims----#
	    possum:updatetheoryobj $w1
	    set entries($w1,outsize_nx) $entries($w,inNx)
	    set entries($w1,outsize_ny) $entries($w,inNy)
	    set entries($w1,outsize_nz) $entries($w,inNz)
	    set entries($w1,outsize_dx) [string range $entries($w1,vcX) 0 2]
	    set entries($w1,outsize_dy) [string range $entries($w1,vcY) 0 2]
	    set entries($w1,outsize_dz) [string range $entries($w1,vcZ) 0 2]
	    set entries($w1,zstart) 0

	    set entries($w1,te) 0.01
	    set entries($w1,tr) 0.7
	    set entries($w1,noisesigma) 0

	    possum:updateFOV $w1

	    set imdimlf [$w1.nb getframe imdims]
	    frame $w1.n
		    label $w1.n.lab -text "Number of Voxels: " -width 15 -anchor w -justify left 
		    LabelSpinBox $w1.n.x -label " X "  -width 6 \
			 -textvariable entries($w1,outsize_nx) \
			-command "$w1.n.x.spin.e validate; possum:updateFOV $w1" \
			-modifycmd "possum:updateFOV $w1" \
			-range { 1   10000  1 }
		    LabelSpinBox $w1.n.y -label " Y "  -width 6 \
			 -textvariable entries($w1,outsize_ny) \
			-command "$w1.n.y.spin.e validate; possum:updateFOV $w1" \
			-modifycmd " possum:updateFOV $w1" \
			-range { 1   10000  1 }
		    LabelSpinBox $w1.n.z -label " Z "  -width 6 \
			 -textvariable entries($w1,outsize_nz) \
			 -command "$w1.n.z.spin.e validate; possum:updateFOV $w1" \
			 -modifycmd " possum:updateFOV $w1" \
			-range { 1   10000  1 }
	     pack $w1.n.lab $w1.n.x $w1.n.y $w1.n.z -in $w1.n -side left -anchor w -padx 3 -pady 3

	    frame $w1.d
		    label $w1.d.lab -text "Voxel Size (mm): " -width 15 -anchor w -justify left 
		    LabelSpinBox $w1.d.x -label " X " -width 6 \
			 -textvariable entries($w1,outsize_dx) \
			-command "$w1.d.x.spin.e validate; possum:updateFOV $w1" \
			-modifycmd "possum:updateFOV $w1" \
			-range { 0.000001   10000.0  0.1 }
		   
		    LabelSpinBox $w1.d.y -label " Y "  -width 6 \
			 -textvariable entries($w1,outsize_dy) \
			-command "$w1.d.y.spin.e validate; possum:updateFOV $w1" \
			-modifycmd "possum:updateFOV $w1" \
			-range { 0.000001   10000.0  0.1 }
		   
		    LabelSpinBox $w1.d.z -label " Z "  -width 6 \
			 -textvariable entries($w1,outsize_dz) \
			-command "$w1.d.z.spin.e validate; possum:updateFOV $w1" \
			-modifycmd "possum:updateFOV $w1" \
			-range { 0.000001   10000.0  0.1 }
	     pack $w1.d.lab $w1.d.x $w1.d.y $w1.d.z -in $w1.d -side left -anchor w -padx 3 -pady 3

	    frame $w1.fov
		    label $w1.fov.lab -text "Field of view (mm): " -width 15 -anchor w -justify left 
		    LabelSpinBox $w1.fov.x -label " X " \
			 -textvariable entries($w1,fov_x) -width 6 \
			 -command "$w1.fov.x.spin.e validate; possum:updateVSIZE $w1" \
			 -modifycmd "possum:updateVSIZE $w1" \
			-range { 0.000001   10000.0  0.1 }
		    LabelSpinBox $w1.fov.y -label " Y "  -width 6 \
			 -textvariable entries($w1,fov_y) \
			 -command "$w1.fov.y.spin.e validate; possum:updateVSIZE $w1" \
			 -modifycmd "possum:updateVSIZE $w1" \
			-range { 0.000001   10000.0  0.1 }
		    LabelSpinBox $w1.fov.z -label " Z "  -width 6 \
			 -textvariable entries($w1,fov_z) \
			 -command "$w1.fov.z.spin.e validate; possum:updateVSIZE $w1" \
			 -modifycmd "possum:updateVSIZE $w1" \
			-range { 0.000001   10000.0  0.1 }
	     pack $w1.fov.lab $w1.fov.x $w1.fov.y $w1.fov.z -in $w1.fov -side left -anchor w -padx 3 -pady 3

	     pack $w1.n $w1.d $w1.fov -in $imdimlf -side top -anchor w -pady 3 -padx 5

	     #----mrpar----#
	     set mrparlf [$w1.nb getframe mrpar]
	     frame $w1.t
		     label $w1.t.spc1 -text "" -width 4
		     label $w1.t.lab -text "" -width 0
		     LabelSpinBox $w1.t.x -label "TE (s)   " -width 8 \
		     	-textvariable entries($w1,te)
		     label $w1.t.spc2 -text "" -width 3
		     LabelSpinBox $w1.t.y -label "TR (s)  " -width 8 \
		     	-textvariable entries($w1,tr) \
		     	-command "$w1.t.y.spin.e validate"
	     pack $w1.t.spc1 $w1.t.lab $w1.t.x $w1.t.spc2 $w1.t.y -in $w1.t -side left -anchor center -padx 3 -pady 10

	     frame $w1.mrpar
		FileEntry $w1.mrpar.dir \
			-textvariable entries($w1,mrpar) \
			-label "MR parameters  " \
			-title "Select" \
			-width 40 \
			-filedialog directory \
			-command "pack forget $w1.mrpar.prev"
		button $w1.mrpar.view -width 20 -command "possum:showMRpar $w1 1"  -text "View MRpar File"
		frame $w1.mrpar.prev -relief raised -borderwidth 1
			label $w1.mrpar.prev.contents1 -text "" -justify center
			label $w1.mrpar.prev.contents2 -text "" -justify center
			label $w1.mrpar.prev.contents3 -text "" -justify center
			label $w1.mrpar.prev.contents4 -text "" -justify center
	     pack $w1.mrpar.prev.contents1 $w1.mrpar.prev.contents2 $w1.mrpar.prev.contents3 $w1.mrpar.prev.contents4 -in $w1.mrpar.prev -anchor w -side left -padx 5 -pady 3
	     pack $w1.mrpar.dir -in $w1.mrpar -anchor w -padx 3 -pady 3
	     pack $w1.mrpar.view -in $w1.mrpar -anchor center -padx 3 -pady 3
	     pack $w1.t $w1.mrpar -in $mrparlf -side top -anchor w -pady 3 -padx 5

	    #----noise----#
	     set noiself [$w1.nb getframe noise]
	     frame $w1.noise
	     LabelSpinBox $w1.noise.snr -label "Sigma: " -width 8 \
		-textvariable entries($w1,noisesigma) -range { 0.0   100000000.0  0.1 } -disabledbackground gray
	     pack $w1.noise.snr -in $w1.noise -anchor w -padx 3 -pady 3
	     pack $w1.noise -in $noiself -side top -anchor w -pady 3 -padx 5

	    #----run----#
	    set runlf [$w1.nb getframe run]
	    frame $w1.out
	    	FileEntry $w1.out.dir \
			-textvariable entries($w1,theoryout) \
			-label "Output file  " \
			-title "Select" \
			-width 50 \
			-filedialog directory
	    	label $w1.out.donelab -text "" -width 50
	    	pack $w1.out.dir $w1.out.donelab -in $w1.out -padx 5 -pady 5
	    pack $w1.out -in $runlf -side top -anchor w -pady 3 -padx 5

    #----outside nb-----#
    frame $w1.btns -relief raised -borderwidth 1
    button $w1.btns.go     -command "possum:theoryGo $w1" \
	-text "Go" -width 10
    button $w1.btns.cancel    -command "destroy $w1" \
	    -text "Exit" -width 10
    pack $w1.btns.go $w1.btns.cancel -in $w1.btns -side left -expand yes -padx 3 -pady 10 -fill y

    frame $w1.status -relief raised -borderwidth 1
	label $w1.status.lab -text ""
    pack $w1.status.lab -in $w1.status -padx 3

    pack $w1.nb $w1.btns $w1.status -expand yes -fill both
##    pack $w1.f -in $w1 -side left -padx 8 -pady 6 -expand yes -fill both
}

proc possum:updatetheoryobj { w args } {
    global entries FSLDIR
    if { $entries($w,obvol) == "" || [file exists $entries($w,obvol)] == 0 } {
	$w.status.lab configure -text "File $entries($w,obvol) not found!" -foreground red
	return
    }
    set entries($w,vcX) [ exec sh -c "$FSLDIR/bin/fslval $entries($w,obvol) pixdim1" ]  
    set entries($w,vcY) [ exec sh -c "$FSLDIR/bin/fslval $entries($w,obvol) pixdim2" ]  
    set entries($w,vcZ) [ exec sh -c "$FSLDIR/bin/fslval $entries($w,obvol) pixdim3" ]  
    set entries($w,inNx) [ exec sh -c "$FSLDIR/bin/fslval $entries($w,obvol) dim1" ]  
    set entries($w,inNy) [ exec sh -c "$FSLDIR/bin/fslval $entries($w,obvol) dim2" ]  
    set entries($w,inNz) [ exec sh -c "$FSLDIR/bin/fslval $entries($w,obvol) dim3" ] 
    set entries($w,inNt) [ exec sh -c "$FSLDIR/bin/fslval $entries($w,obvol) dim4" ]
}
   
proc possum:showMRpar { w {theory 0} } {
    global entries FSLDIR
    $w.mrpar.prev.contents1 configure -text ""
    $w.mrpar.prev.contents2 configure -text ""
    $w.mrpar.prev.contents3 configure -text ""
    $w.mrpar.prev.contents4 configure -text ""

    if { $entries($w,mrpar) == "" || [file exists $entries($w,mrpar)] == 0 } {
	set mrparcontents "MRpar file $entries($w,mrpar) not found!"
	$w.mrpar.prev.contents1 configure -text "$mrparcontents" -foreground red
	pack $w.mrpar.prev -in $w.mrpar -anchor center -padx 3 -pady 3
	$w.nb compute_size
	return
    }
    catch { exec sh -c "cat $entries($w,mrpar)" } mrparcontents
    if { $mrparcontents == "" } {
	set mrparcontents "File Empty!"
	$w.mrpar.prev.contents1 configure -text "$mrparcontents" -foreground red
	pack $w.mrpar.prev -in $w.mrpar -anchor center -padx 3 -pady 3
	$w.nb compute_size
	return
    } else {
	if { $theory == 0 } {
		set t2title "T2*(s)"
		for { set i 1 } { $i <= 4 } { incr i } {
			catch { exec sh -c "awk '{print $$i}' $entries($w,mrpar)" } mrparcontents
			switch $i {
					1 {set mrparcontents "T1(s)\n$mrparcontents"}
					2 {set mrparcontents "$t2title\n$mrparcontents"}
					3 {set mrparcontents "Proton Density\n$mrparcontents"}
					4 {set mrparcontents "Chemical Shift\n$mrparcontents"}
			}
	    		$w.mrpar.prev.contents$i configure -text "$mrparcontents" -foreground black
		}
	} else {
		set t2title "T2(s)"
		for { set i 1 } { $i <= 3 } { incr i } {
			catch { exec sh -c "awk '{print $$i}' $entries($w,mrpar)" } mrparcontents
			switch $i {
					1 {set mrparcontents "T1(s)\n$mrparcontents"}
					2 {set mrparcontents "$t2title\n$mrparcontents"}
					3 {set mrparcontents "Proton Density\n$mrparcontents"}
			}
	    		$w.mrpar.prev.contents$i configure -text "$mrparcontents" -foreground black
		}
	}
	pack $w.mrpar.prev -in $w.mrpar -anchor center -padx 3 -pady 3
	$w.nb compute_size
    }
}

proc possum:theoryGo { w } {
	$w.status.lab configure -text "Generating image..." -foreground black
	update idletasks
	Possum:generateTheoryVolume $w
}

proc Possum:generateTheoryVolume { w } {
    global entries FSLDIR POSSUMDIR
    set mrpar $entries($w,mrpar)
    set phantom $entries($w,obvol)
    set output $entries($w,theoryout)
    set sigma $entries($w,noisesigma)

    if { $mrpar == "" || [file exists $mrpar] == 0} then {
	set msg "Error! MRpar file not found!"
	$w.status.lab configure -text $msg -foreground red
	return 0
    }
    if { $phantom == "" || [file exists $phantom] == 0} then {
	set msg "Error! Input Object not found!"
	$w.status.lab configure -text $msg -foreground red
	return 0
    }
    if { $entries($w,te) == 0 || $entries($w,tr) == 0 } then {
	set msg "Error! TE/TR can't be zero!"
	$w.status.lab configure -text $msg -foreground red
	return 0;
    }
    if { $output == "" || [ file exists [exec sh -c "dirname $output"] ] == 0 } then {
	set msg "Error! output directory doesn't exist!"
	$w.status.lab configure -text $msg -foreground red
	return 0
    }
    set outwd [file dirname $output]
    catch { exec sh -c "${FSLDIR}/bin/fslval $phantom dim4" } phsize
    set 1 {$1}
    catch { exec sh -c "wc -l $mrpar | awk '{print $1}'" } size
    if { $phsize != $size } {
	set msg "Error! Input object doesn't have enough MRpar values."
	$w.status.lab configure -text $msg -foreground red
	return 0
    }

    set imdimx $entries($w,outsize_nx)
    set imdimy $entries($w,outsize_ny)
    set imdimz $entries($w,outsize_nz)
    set voxdimx $entries($w,outsize_dx)
    set voxdimy $entries($w,outsize_dy)
    set voxdimz $entries($w,outsize_dz)
    set zstart [ expr round($entries($w,zstart)/$voxdimz) ]

    regsub -all ".nii.gz" $output "" outfn
    set logfile "${outfn}.log"
    set command "${POSSUMDIR}/bin/tcalc -i $phantom -o $outfn --mrpar=$mrpar --te=$entries($w,te) --tr=$entries($w,tr) --nx=$imdimx --ny=$imdimy --nz=$imdimz --dx=$voxdimx --dy=$voxdimy --dz=$voxdimz --zstart=$zstart --sigma=$sigma -v"
    catch { exec sh -c "$command 3>&1 1> $logfile 2>&3-" } oval
    if { $oval == "" } {
	set msg "Volume generated : '${output}'"
    } else {
	puts $oval
	set msg "Could not generate image! Check $logfile for details"
	$w.status.lab configure -foreground red
    }
    $w.status.lab configure -text $msg
}

proc possum:savescriptdialog { w } {
	global entries FSLDIR PWD
	set count 0
	set w1 ".dialog[incr count]"
	while { [ winfo exists $w1 ] } {
		set w1 ".dialog[incr count]"
	}
	toplevel $w1
	wm title $w1 "Save Commands to script"
	wm iconname $w1 "Possum"
	wm iconbitmap $w1 @${FSLDIR}/tcl/fmrib.xbm

	set entries($w1,scriptpath) "$PWD/possum-script"
	frame $w1.opts -relief raised -borderwidth 1
		FileEntry $w1.opts.scriptdir \
			-textvariable entries($w1,scriptpath) \
			-filetypes IMAGE \
			-label "Output File: " \
			-title "Select" \
			-width 40 \
			-filedialog directory
	pack $w1.opts.scriptdir -in $w1.opts -padx 3 -pady 5
	frame $w1.bttns -relief raised -borderwidth 1
		button $w1.bttns.go -text "Go" -command "$w1.warn.lab configure -text '' -foreground black; update idletasks; Possum:savescripts $w $w1" -width 10
		button $w1.bttns.cancel -text "Exit" -command "destroy $w1" -width 10
	pack $w1.bttns.go $w1.bttns.cancel -in $w1.bttns -side left -expand yes -padx 3 -pady 10 -fill y
	frame $w1.warn -relief raised -borderwidth 1
		label  $w1.warn.lab -text ""
	pack $w1.warn.lab -in $w1.warn -padx 3 -pady 3

	pack $w1.opts $w1.bttns $w1.warn -expand yes -fill both				
}

proc Possum:savescripts { w w1 } {
	global entries
	if { $entries($w1,scriptpath) == "" } {
		$w1.warn.lab configure -text "Script path cannot be empty!" -foreground red
		return
	}
	set status [possum:savescripts $w $entries($w,seqtype) $entries($w,te) $entries($w,tr) $entries($w,trslc) $entries($w,outsize_nx) $entries($w,outsize_ny) $entries($w,outsize_nz) $entries($w,outsize_dx) $entries($w,outsize_dy) $entries($w,outsize_dz) $entries($w,zstart) $entries($w,flipangle) $entries($w,numvol) $entries($w,maxG) $entries($w,riseT) $entries($w,bw) $entries($w,gap) $entries($w,slcselect) $entries($w,phencode) $entries($w,readgrad) $entries($w,cover) $entries($w,pluss) $entries($w,pluss) $entries($w1,scriptpath) $entries($w1,scriptpath)]
	if { $status != "" } {
		$w1.warn.lab configure -text "$status" -foreground red
	} else {
		$w1.warn.lab configure -text "Script saved! '$entries($w1,scriptpath)'" -foreground black
	}
}

proc possum:savescripts { w seqtype te tr trslc nx ny nz dx dy dz zstart angle nvol maxg riset bw gap slcdir phasedir readdir coverage plus pluss scriptpath outpath} {
	global entries FSLDIR POSSUMDIR
	## Read possum script template	##CHECK PATH!	## CHANGE SCRIPT NAME
	if { [ file exists [ file dirname $scriptpath] ] == 0 } { return "Script path cannot be accessed" }

	set templatepath "${FSLDIR}/data/possum/possum_script_template"
	if { [ file exists $templatepath ] == 0 } { return "Template file not found!" }

	## Check all the compulsary files
	set fileids [ list "Input Object" "Output Directory" "Motion file" "MRpar file" "Slice Profile file" ]
	set filevals [ list $entries($w,obvol) [file dirname $outpath] $entries($w,mot) $entries($w,mrpar) $entries($w,slcprof) ]
	foreach i $fileids j $filevals {
		if { $j == "" || [file exists $j] == 0 } {
			puts "(Not found!) $i : $j"
			return "$i not found!"
		}
	}
	if { $entries($w,custompulse_yn) == 1 } {
		return "Error! Script cannot be generated for custom pulse sequence!"
	}
	if { $entries($w,numproc) == 0 } {
		return "Error! Number of processors cannot be 0!"
	}

	## Check pulse sequence
	Possum:pulsecheck $w
	if { $entries($w,pulsechecktest) == 0 } { 
		return "Check pulse sequence!"
	}

	## Open script for writing
	set filenm [open $scriptpath w]
	catch { exec sh -c "cat $templatepath" } script

	## Replace date at the top of the script
	regsub -all "<date>" $script "[clock format [clock seconds] -format %d-%m-%Y]" script

	if { $entries($w,motion_yn) == 1 } { set motion $entries($w,mot) } else { set motion "${FSLDIR}/data/possum/zeromotion" }

	set pulsedir "$entries($w,out)/pulse"

	set plcholders [list brain outdir motion mrpar slcprof pulse nprocs segs ntime ]
	set filevals   [list $entries($w,obvol) $entries($w,out) $motion $entries($w,mrpar) $entries($w,slcprof) $pulsedir $entries($w,numproc) $entries($w,segs) 2000 ]

	foreach key $plcholders val $filevals {
		regsub -all "<$key>" $script "$val" script
	}

	set plcholders  [list seqtype te tr trslc nx ny nz dx dy dz zstart angle nvol maxg riset bw gap slcdir phasedir readdir coverage]
	set pulseparams [list $seqtype $te $tr $trslc $nx $ny $nz [expr $dx * 0.001] [expr $dy * 0.001] [expr $dz * 0.001] [expr $zstart * 0.001] $angle $nvol $maxg $riset $bw $gap "${slcdir}${plus}" "${phasedir}${pluss}" "${readdir}${pluss}" $coverage]
	foreach key $plcholders val $pulseparams {
		regsub -all "<$key>" $script "$val" script
	}

	## Possum command optional arguments
	set possumopts ""

	# B0 Field inh
	if { $entries($w,b0inh_yn)  == 1} {
		set b0opts "B0=$entries($w,b0f)"
		if { $entries($w,b0units) == "ppm" } {
			set b0opts "$b0opts\n${FSLDIR}/bin/fslmaths \$B0 -mul $entries($w,b0fieldstrength) -div 1000000 \$OUTDIR/b0"
		} else {
			set b0opts "$b0opts\ncp \$B0 \$OUTDIR/b0.nii.gz"
		}
		if { $entries($w,motion_yn) == 1 } {
			set b0opts "$b0opts\n${FSLDIR}/bin/fslroi \$OUTDIR/b0 \$OUTDIR/b0x_dx.nii.gz 8 1"
			set b0opts "$b0opts\n${FSLDIR}/bin/fslroi \$OUTDIR/b0 \$OUTDIR/b0x_dy.nii.gz 7 1"
			set b0opts "$b0opts\n${FSLDIR}/bin/fslroi \$OUTDIR/b0 \$OUTDIR/b0x_dz.nii.gz 6 1"
			set b0opts "$b0opts\n${FSLDIR}/bin/fslroi \$OUTDIR/b0 \$OUTDIR/b0y_dx.nii.gz 5 1"
			set b0opts "$b0opts\n${FSLDIR}/bin/fslroi \$OUTDIR/b0 \$OUTDIR/b0y_dy.nii.gz 4 1"
			set b0opts "$b0opts\n${FSLDIR}/bin/fslroi \$OUTDIR/b0 \$OUTDIR/b0y_dz.nii.gz 3 1"
			set b0opts "$b0opts\n${FSLDIR}/bin/fslroi \$OUTDIR/b0 \$OUTDIR/b0z_dx.nii.gz 2 1"
			set b0opts "$b0opts\n${FSLDIR}/bin/fslroi \$OUTDIR/b0 \$OUTDIR/b0z_dy.nii.gz 1 1"
			set b0opts "$b0opts\n${FSLDIR}/bin/fslroi \$OUTDIR/b0 \$OUTDIR/b0z_dz.nii.gz 0 1"
		} else {
			set b0opts "$b0opts\n${FSLDIR}/bin/fslroi \$OUTDIR/b0 \$OUTDIR/b0z_dz.nii.gz 0 1"
		}
		regsub -all "<b0file>" $script "$b0opts\n" script
		set possumopts "$possumopts -b \$OUTDIR/b0"
	} else {
		regsub -all "<b0file>" $script "" script
	}

	## B0 Field inh (changing in time)
	if {  $entries($w,b0inhtime_yn) == 1 } {
		regsub -all "<b0extra>" $script "B0EXTRA=$entries($w,b0ftime)\nB0EXTRATIMECOURSE=$entries($w,b0ftimecourse)" script
		set possumopts "$possumopts --b0extra=\$B0EXTRA --b0time=\$B0EXTRATIMECOURSE"
	} else {
		regsub -all "<b0extra>" $script "" script
	}

	## Activation
	if { $entries($w,activ_yn) == 1 } {
		regsub -all "<activ-timecourse>" $script "ACTTIMECOURSE=$entries($w,act2)" script
		regsub -all "<activ-spatmod>" $script "ACTSPATMOD=$entries($w,act1)" script
		set possumopts "$possumopts -a \$ACTSPATMOD -t \$ACTTIMECOURSE"
	} else {
		regsub -all "<activ-timecourse>" $script "" script
		regsub -all "<activ-spatmod>" $script "" script
	}

	## Noise
	regsub -all "<sigma>" $script "$entries($w,noisesigma)" script
	regsub -all "<snr>" $script "$entries($w,noisesnr)" script
	if { $entries($w,noise_yn) == 1 } {
		if { $entries($w,noiseunits) == "snr" } {
			regsub -all "<noise>" $script "echo snr \$SNR > \$OUTDIR/noise" script
			regsub -all "<noise-comment>" $script "\#\#Noise added with SIGMA value" script
		} else {
			regsub -all "<noise>" $script "echo sigma \$SIGMA > \$OUTDIR/noise" script
			regsub -all "<noise-comment>" $script "\#\#Noise added with SNR value" script
		}
	} else {
		regsub -all "<noise>" $script "echo sigma 0 > \$OUTDIR/noise" script
		regsub -all "<noise-comment>" $script "\#\#No Noise added" script
	}

	## RF averaging
#	if { $entries($w,rfavg_yn) == 1} {
#		set possumopts "$possumopts --rfavg"
#	}

	## Add POSSUM options to possum-command
	regsub -all "<possumopts>" $script "$possumopts" script

	## Remove unnecessary blank spaces
	regsub -all "\n\n\n" $script "" script

	## Write script to file
	puts $filenm $script
	close $filenm

	## temp
	return
}

proc possum:settimings { w } {
	global entries
	if { $entries($w,seqtype) == "epi" } {
		set entries($w,te) 0.030
		set entries($w,tr) 3
		set entries($w,trslc) 0.12
	} elseif {$entries($w,seqtype) == "ge" } {
		set entries($w,te) 0.01
		set entries($w,tr) 0.7
		set entries($w,trslc) 0
	} else {
		set entries($w,te) 0
		set entries($w,tr) 0
		set entries($w,trslc) 0
	}
	update idletasks
}

proc possum:previewactiv { w } {
    global entries FSLDIR
    set convertcom "${FSLDIR}/bin/pngappend"
    set filenm $entries($w,act1)
    set validim 0
    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	catch { exec sh -c "${FSLDIR}/bin/imtest $entries($w,obvol)" } oval
	if { $oval == 1 } {
	    catch { exec sh -c "${FSLDIR}/bin/slicer $filenm $entries($w,obvol) -a ${tmpnam}.png" } oval 
	} else {
	    catch { exec sh -c "${FSLDIR}/bin/slicer $filenm -a ${tmpnam}.png" } oval 
	}
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png ${tmpnam}" } oval
    }
    if { $validim == 1 } {
	$w.activim.label configure -image $graphpic
    } else { 
	$w.activim.label configure -image "" -text "Could not generate preview"
    }
    $w.nb compute_size
}

proc possum:savepulsedialog { w } {
    global entries FSLDIR PWD
    set count 0
    set w1 ".dialog[incr count]"
    while { [ winfo exists $w1 ] } {
        set w1 ".dialog[incr count]"
    }
    toplevel $w1
    wm title $w1 "Save Pulse-Sequence"
    wm iconname $w1 "Possum"
    wm iconbitmap $w1 @${FSLDIR}/tcl/fmrib.xbm

    set entries($w1,pulse) "$PWD/pulse"
    frame $w1.opts -relief raised -borderwidth 1
	    FileEntry $w1.opts.dir \
		-textvariable entries($w1,pulse) \
		-filetypes IMAGE \
		-label "Output File: " \
		-title "Select" \
		-width 40 \
		-filedialog directory
    pack $w1.opts.dir -in $w1.opts -padx 3 -pady 5
    frame $w1.bttns -relief raised -borderwidth 1
	    button $w1.bttns.go -text "Go" -command "possum:savepulse $w $w1" -width 10
	    button $w1.bttns.cancel -text "Exit" -command "destroy $w1" -width 10
    pack $w1.bttns.go $w1.bttns.cancel -in $w1.bttns -side left -expand yes -padx 3 -pady 10 -fill y
    frame $w1.warn -relief raised -borderwidth 1
	label  $w1.warn.lab -text ""
    pack $w1.warn.lab -in $w1.warn -padx 3 -pady 3

    pack $w1.opts $w1.bttns $w1.warn -expand yes -fill both
}

proc possum:savepulse { w  w1 } {
	global entries FSLDIR POSSUMDIR PWD
        $w1.warn.lab configure -text "Generating pulse-sequence..." -foreground black
	update idletasks
	set dx [ expr $entries($w,outsize_dx) * 0.001 ]
	set dy [ expr $entries($w,outsize_dy) * 0.001 ]
	set dz [ expr $entries($w,outsize_dz) * 0.001 ]
	set zs [ expr $entries($w,zstart) * 0.001 ]
	set gap [ expr $entries($w,gap) * 0.001 ]

	set outwd [file dirname $entries($w1,pulse) ]

	if { [ file exists $outwd ] == 0 } {
		$w1.warn.lab configure -text "Error: '$outwd' - No such directory!"
		return
	}

	if { [ file exists $entries($w,obvol) ] == 0 } {
		$w1.warn.lab configure -text "Error: '$entries($w,obvol)' - No such file!"
		return
	}

	set pulsecom "${POSSUMDIR}/bin/pulse -i $entries($w,obvol) -o $entries($w1,pulse) --seq=$entries($w,seqtype) --te=$entries($w,te) --tr=$entries($w,tr) --trslc=$entries($w,trslc) --nx=$entries($w,outsize_nx) --ny=$entries($w,outsize_ny) --numslc=$entries($w,outsize_nz) --dx=$dx --dy=$dy --slcthk=$dz --numvol=$entries($w,numvol) --zstart=$zs --gap=$gap --bw=$entries($w,bw) --readdir=$entries($w,readgrad)$entries($w,pluss) --phasedir=$entries($w,phencode)$entries($w,pluss) --slcdir=$entries($w,slcselect)$entries($w,plus) --maxG=$entries($w,maxG) --riset=$entries($w,riseT) --angle=$entries($w,flipangle) --cover=$entries($w,cover) -v"
	Possum:pulsecheck $w
	if { $entries($w,pulsechecktest) == 0 } {
		$w1.warn.lab configure -text "Check parameters!" -foreground red
		return
	}
	# Temporarily use pipe 3 to store stderr (2). This is because stderr cannot be piped, only stdout can.
	catch { exec sh -c "$pulsecom 3>&1 1 > $outwd/pulse.log 2>&3- | tee -a $outwd/pulse.log" } pulsestatus
	if { $pulsestatus != "" } {
		$w1.warn.lab configure -text "Error: Pulse-Sequence not generated. Check pulse.log" -foreground red
	} else {
		$w1.warn.lab configure -text "Pulse-Sequence Saved!"
	}
}

proc possum:custompulse { w } {
	global entries FSLDIR POSSUMDIR PWD
	set pulself [$w.nb getframe pulse]
	if { $entries($w,seqtype) == "custom" } {
		set entries($w,custompulse_yn) 1
		set entries($w,cuspulse) "$PWD/pulse"
		pack forget $w.t $w.defaultframe $w.customwarn.lab
		pack $w.custompulse -in $pulself -anchor w -padx 3 -pady 3
	} else {
		set entries($w,custompulse_yn) 0
		pack forget pack $w.custompulse
		pack $w.t -in $w.topopts
		pack $w.topopts $w.defaultframe -in $pulself -anchor nw -side top -padx 3 -pady 3
	}
}

proc possum:checkloadedpulse { w dummy } {
	global entries FSLDIR
	set loaded $entries($w,cuspulse)
	set filelist [list "$loaded" "$loaded.info" "$loaded.posx" "$loaded.posy" "$loaded.posz" ]
	set msg ""
	foreach i $filelist {
		if { [file exists $i] == 0 } {
			set msg "$msg, [file tail $i]"
		}
	}
	if { $msg != "" } {
		$w.customwarn.lab configure -text "Required File(s) $msg not found" -foreground red
		pack $w.customwarn.lab -in $w.customwarn -side top -anchor center -padx 3 -pady 3
	} else {
		pack forget $w.customwarn.lab
	}
}

proc possum:previewslices { w } {
    global entries FSLDIR
    set count 0
    set w1 ".dialog[incr count]"
    while { [ winfo exists $w1 ] } {
        set w1 ".dialog[incr count]"
    }
    toplevel $w1
    wm title $w1 "Preview of Slice Prescription"
    wm iconname $w1 "SlicePreview"
    wm iconbitmap $w1 @${FSLDIR}/tcl/fmrib.xbm
    frame $w1.sprev
    label $w1.sprev.label -image "" -text "\n    Generating preview ... please wait    \n\n"
    pack $w1.sprev.label -in $w1.sprev
    pack $w1.sprev -in $w1
    # force this message to popup now
    update
    set convertcom "${FSLDIR}/bin/pngappend"
    set filenm $entries($w,obvol)
    set validim 0
    catch { exec sh -c "${FSLDIR}/bin/imtest $filenm" } oval
    if { $oval == 1 } {
	catch { exec sh -c "${FSLDIR}/bin/tmpnam" } tmpnam
	set slcselect $entries($w,slcselect)
	catch { exec sh -c "${FSLDIR}/bin/fslroi $filenm ${tmpnam}_sw 0 1" } oval
	catch { exec sh -c "${FSLDIR}/bin/fslswapdim ${tmpnam}_sw $entries($w,readgrad) $entries($w,phencode) $slcselect ${tmpnam}_sw" } oval
	catch { exec sh -c "${FSLDIR}/bin/fslval ${tmpnam}_sw dim1" } in_nx
	catch { exec sh -c "${FSLDIR}/bin/fslval ${tmpnam}_sw dim2" } in_ny
	catch { exec sh -c "${FSLDIR}/bin/fslval ${tmpnam}_sw pixdim1" } dx	
	catch { exec sh -c "${FSLDIR}/bin/fslval ${tmpnam}_sw pixdim2" } dy	
	catch { exec sh -c "${FSLDIR}/bin/fslval ${tmpnam}_sw pixdim3" } dz	
	catch { exec sh -c "${FSLDIR}/bin/fslstats ${tmpnam}_sw -r | awk '{ print \$2 }'" } imax
	set imax [ expr $imax*6 ]
	set zstart [ expr round($entries($w,zstart)/$dz) ]
	set xsize [ expr  round($entries($w,outsize_nx)*$entries($w,outsize_dx)/$dx) ]
	set ysize [ expr  round($entries($w,outsize_ny)*$entries($w,outsize_dy)/$dy) ]
	set zsize [ expr  round($entries($w,outsize_nz)*($entries($w,outsize_dz)+$entries($w,gap))/$dz) ]
	set xstart [ expr round(($in_nx-$xsize)/2.0) ]
	if { $xstart < 0 } { set xstart 0 }
	set ystart [ expr round(($in_ny-$ysize)/2.0) ]
	if { $ystart < 0 } { set ystart 0 }
	# puts "roi $xstart $xsize $ystart $ysize $zstart $zsize"
	catch { exec sh -c "${FSLDIR}/bin/fslmaths ${tmpnam}_sw -roi $xstart $xsize $ystart $ysize $zstart $zsize 0 1 -mul 5 -add ${tmpnam}_sw ${tmpnam}_sw" } oval	
	catch { exec sh -c "${FSLDIR}/bin/slicer ${tmpnam}_sw -i 0 $imax -a ${tmpnam}.png" } oval
	catch { exec sh -c "$convertcom ${tmpnam}.png ${tmpnam}.gif" } oval
	if { [ file exists ${tmpnam}.gif ] } {
	    set graphpic [image create photo -file ${tmpnam}.gif ]
	    set validim 1
	}
	catch { exec sh -c "rm -f ${tmpnam}.gif ${tmpnam}.png ${tmpnam}_sw* ${tmpnam}" } oval
    }
    if { $validim == 1 } {
	$w1.sprev.label configure -image $graphpic
    } else { 
	$w1.sprev.label configure -image "" -text "Could not generate preview"
    }
    pack forget $w1.sprev.label
    pack forget $w1.sprev
    pack $w1.sprev.label -in $w1.sprev
    button $w1.cancel -command "destroy $w1" -text "Dismiss"
    pack $w1.sprev $w1.cancel -in $w1
    update
}

proc possum:updatemotion { w } {
    global entries guivars FSLDIR
    if { $entries($w,motion_yn) == 1 } {
	pack $w.mot $w.motpreview $w.motprevimg -in  $guivars($w,lfmotion) -side top -anchor w -padx 3 -pady 3
	set entries($w,slcsampfactor) 6
    } else {
	pack forget $w.mot $w.motpreview $w.motprevimg $w.b0maintime.b0im.b0timewarn $w.motwarn
	set entries($w,slcsampfactor) 2
    }
    $w.nb compute_size
}

#proc possum:updaterfavg { w } {
#	global entries guivars FSLDIR
#	if { $entries($w,rfavg_yn) == 1 } {
# 		$w.trs.rfavg select
# 		$w.rfavg select
# 	} else {
# 		$w.trs.rfavg deselect
# 		$w.rfavg deselect
# 	}
# }

proc possum:updateactivation { w } {
    global entries guivars FSLDIR
    if { $entries($w,activ_yn) == 1 } {
	pack $w.activation -in $guivars($w,lfactivation) -anchor w
    } else {
	pack forget $w.activation
    }
}

proc possum:updatenoise { w } {
    global entries guivars FSLDIR
    if { $entries($w,noise_yn) == 1 } {
	pack $w.noiseval -in $guivars($w,lfnoise) -side top -anchor w -padx 3 -pady 3
    } else {
	pack forget $w.noiseval
    }
}
proc possum:updatenoiseunits { w } {
    global entries guivars FSLDIR
    if { $entries($w,noiseunits) == "snr" } {
	$w.noiseval2.sigma configure -state disabled ; $w.noiseval1.snr configure -state normal
    } else {
        $w.noiseval2.sigma configure -state normal ; $w.noiseval1.snr configure -state disabled
    }
}

proc possum:twosigfigs { num } {
    set pten [ expr log10($num) ]
    set pten [ expr floor($pten) - 1 ]
    set pten [ expr exp($pten*log(10)) ]
    set tsf [ expr round($num / $pten) ]
    set tsf [ expr $tsf*$pten ]
    return $tsf
}

proc possum:modTRslice { num } {
    set tsf [ expr $num*0.999 ]
    return $tsf
}

proc possum:updateOBprop { w { filename foo } } {
 global entries FSLDIR
    set filename $entries($w,obvol)
    if { $filename != "" } {
	    set entries($w,vcX) [ exec sh -c "$FSLDIR/bin/fslval $filename pixdim1" ]  
	    set entries($w,vcY) [ exec sh -c "$FSLDIR/bin/fslval $filename pixdim2" ]  
	    set entries($w,vcZ) [ exec sh -c "$FSLDIR/bin/fslval $filename pixdim3" ]  
	    set entries($w,inNx) [ exec sh -c "$FSLDIR/bin/fslval $filename dim1" ]  
	    set entries($w,inNy) [ exec sh -c "$FSLDIR/bin/fslval $filename dim2" ]  
	    set entries($w,inNz) [ exec sh -c "$FSLDIR/bin/fslval $filename dim3" ] 
	    set entries($w,inNt) [ exec sh -c "$FSLDIR/bin/fslval $filename dim4" ] 
    }
    return 0
}

proc possum:updateB0prop { w { filename foo } } {
 global entries FSLDIR
    set filename "$entries($w,b0f)"
    set entries($w,vcXb0) [ exec sh -c "$FSLDIR/bin/fslval $filename pixdim1" ]  
    set entries($w,vcYb0) [ exec sh -c "$FSLDIR/bin/fslval $filename pixdim2" ]  
    set entries($w,vcZb0) [ exec sh -c "$FSLDIR/bin/fslval $filename pixdim3" ]  
    set entries($w,inNxb0) [ exec sh -c "$FSLDIR/bin/fslval $filename dim1" ]  
    set entries($w,inNyb0) [ exec sh -c "$FSLDIR/bin/fslval $filename dim2" ]  
    set entries($w,inNzb0) [ exec sh -c "$FSLDIR/bin/fslval $filename dim3" ] 
    set entries($w,inNtb0) [ exec sh -c "$FSLDIR/bin/fslval $filename dim4" ]
    return 0
}

proc possum:updateB0timeprop { w { filename foo } } {
 global entries FSLDIR
    set filename "$entries($w,b0ftime)"
    set entries($w,vcXb0time) [ exec sh -c "$FSLDIR/bin/fslval $filename pixdim1" ]  
    set entries($w,vcYb0time) [ exec sh -c "$FSLDIR/bin/fslval $filename pixdim2" ]  
    set entries($w,vcZb0time) [ exec sh -c "$FSLDIR/bin/fslval $filename pixdim3" ]  
    set entries($w,inNxb0time) [ exec sh -c "$FSLDIR/bin/fslval $filename dim1" ]  
    set entries($w,inNyb0time) [ exec sh -c "$FSLDIR/bin/fslval $filename dim2" ]  
    set entries($w,inNzb0time) [ exec sh -c "$FSLDIR/bin/fslval $filename dim3" ] 
    return 0
}

proc possum:updateACTprop { w { filename foo } } {
 global entries FSLDIR
    set filename $entries($w,act1)
    set entries($w,vcXact) [ exec sh -c "$FSLDIR/bin/fslval $filename pixdim1" ]  
    set entries($w,vcYact) [ exec sh -c "$FSLDIR/bin/fslval $filename pixdim2" ]  
    set entries($w,vcZact) [ exec sh -c "$FSLDIR/bin/fslval $filename pixdim3" ]  
    set entries($w,inNxact) [ exec sh -c "$FSLDIR/bin/fslval $filename dim1" ]  
    set entries($w,inNyact) [ exec sh -c "$FSLDIR/bin/fslval $filename dim2" ]  
    set entries($w,inNzact) [ exec sh -c "$FSLDIR/bin/fslval $filename dim3" ] 
    set entries($w,inNtact) [ exec sh -c "$FSLDIR/bin/fslval $filename dim4" ] 
    return 0
}

proc possum:updatecomptime { w { filename foo } } {
    global entries FSLDIR
    if { $entries($w,seqtype) == "epi" } {

	    # Quick fix : Multiply computed time by factor
	    set multiplier 2
	    # number of voxels per 1mm3
	    set Nvpv [ expr $entries($w,slcsampfactor) / ( $entries($w,vcX) * $entries($w,vcY) * $entries($w,vcZ) )  ]
	    set dimX [ expr $entries($w,vcX) * $entries($w,inNx) ] 
	    set dimY [ expr $entries($w,vcY) * $entries($w,inNy) ]
	    set dimZ [ expr $entries($w,vcZ) * $entries($w,inNz) ]
	    if { $entries($w,motion_yn) == 0 } {
	       set Zmaxm 0
	       set Nfev 12
		if { $entries($w,b0inh_yn) == 1 } {
		    set Nfev 94
		}
	    } else {
	       set Zmaxm 5
	       set Nfev 91
	    }
	    set Zmaxm [ expr $Zmaxm * $multiplier]
	    set Nfev [ expr $Nfev * $multiplier]
	    # max mm extra due to slc profile
	    set SlcP [expr $entries($w,outsize_dz) * 1.5 ]
	    # number of events in the pulse sequence
	    set Nevent [expr $entries($w,outsize_nx) * $entries($w,outsize_ny) * ( ($Zmaxm + $SlcP ) / $entries($w,outsize_dz) + 1 ) * $entries($w,numvol) * $entries($w,cover)/100.0 ]
	    # number of voxels
	    set gap [expr $entries($w,gap) * ($entries($w,outsize_nz) - 1 ) ]
	#set proportion of the non-zero voxels
	    set P 0.4
	    set Nvoxel [ expr $Nvpv * $dimX * $dimY * ( $entries($w,outsize_nz) * $entries($w,outsize_dz) + $Zmaxm + $SlcP + $gap ) * $entries($w,inNt) ]
	    # number of flops per event and per voxel
	    # total number of flops 
	    set Nflops [ expr $Nevent * $Nvoxel * $Nfev ]
	    # time for N computer proc capable of 1 Giga Flop (in seconds)
	    set tottime [ expr $Nflops * 0.000000001 / $entries($w,numproc) ]
	    if { $tottime < 1 } { set tottime 1 }
	    set entries($w,proctime) [ expr int ( $tottime / 60 ) ]
	    if { $tottime > 86400 } {
		set entries($w,comptime) "[ possum:twosigfigs [ expr $tottime / 86400 ] ] days"
		return
	    }
	    if { $tottime > 3600 } {
		set entries($w,comptime) "[ possum:twosigfigs [ expr $tottime / 3600 ] ] hours"
		return
	    }
	    if { $tottime > 60 } {
		set entries($w,comptime) "[ possum:twosigfigs [ expr $tottime / 60 ] ] minutes"
		return
	    } else {
		set entries($w,comptime) "[ possum:twosigfigs $tottime ] seconds"
	    }
    } elseif { $entries($w,seqtype) == "ge" } {
	    # number of voxels per 1mm3
	    set Nvpv [ expr 1 / ( $entries($w,vcX) * $entries($w,vcY) * $entries($w,vcZ) )  ]
	    set dimX [ expr $entries($w,vcX) * $entries($w,inNx) ] 
	    set dimY [ expr $entries($w,vcY) * $entries($w,inNy) ]
	    set dimZ [ expr $entries($w,vcZ) * $entries($w,inNz) ]
	    if { $entries($w,motion_yn) == 0 } {
	       set Zmaxm 0
	       set Nfev 50
		if { $entries($w,b0inh_yn) == 1 } {
		    set Nfev 150
		}
	    } else {
	       set Zmaxm 5
	       set Nfev 380
	    }
	    # max mm extra due to slc profile
	    set SlcP [expr $entries($w,outsize_dz) * 1.5 ]
	    # number of events in the pulse sequence
	    set Nevent [expr $entries($w,outsize_nx) * $entries($w,outsize_ny) * ( ($Zmaxm + $SlcP ) / $entries($w,outsize_dz) + 1 ) * $entries($w,numvol) * $entries($w,cover)/100.0 ]
	    # number of voxels
	    set gap [expr $entries($w,gap) * ($entries($w,outsize_nz) - 1 ) ]
	#set proportion of the non-zero voxels
	    set P 0.4
	    set Nvoxel [ expr $Nvpv * $dimX * $dimY * ( $entries($w,outsize_nz) * $entries($w,outsize_dz) + $Zmaxm + $SlcP + $gap ) * $entries($w,inNt) ]
	    # number of flops per event and per voxel
	    # total number of flops 
	    set Nflops [ expr $Nevent * $Nvoxel * $Nfev ]
	    # time for N computer proc capable of 1 Giga Flop (in seconds)
	    set tottime [ expr $Nflops * 0.000000001 / $entries($w,numproc) ]
	    if { $tottime < 1 } { set tottime 1 }
	    set entries($w,proctime) [ expr int ( $tottime / 60 ) ]
	    if { $tottime > 86400 } {
		set entries($w,comptime) "[ possum:twosigfigs [ expr $tottime / 86400 ] ] days"
		return
	    }
	    if { $tottime > 3600 } {
		set entries($w,comptime) "[ possum:twosigfigs [ expr $tottime / 3600 ] ] hours"
		return
	    }
	    if { $tottime > 60 } {
		set entries($w,comptime) "[ possum:twosigfigs [ expr $tottime / 60 ] ] minutes"
		return
	    } else {
		set entries($w,comptime) "[ possum:twosigfigs $tottime ] seconds"
	    }
    }
}

proc possum:updateechosp { w } {
global entries 
    set dx [ expr $entries($w,outsize_dx) * 0.001 ]
    set dy [ expr $entries($w,outsize_dy) * 0.001 ]
    set dz [ expr $entries($w,outsize_dz)* 0.001 ]
    set zs [ expr $entries($w,zstart) * 0.001 ]
    # checks that the Pulse Sequence parameters are appropriate
    set gammabar 42580000
    set tana [expr $entries($w,maxG)/$entries($w,riseT) ]
    set dtx [ expr 1.0/$entries($w,bw)]
    set dkx [ expr 1.0/($entries($w,outsize_nx)*${dx})]
    set dky [ expr 1.0/($entries($w,outsize_ny)*${dy})]
    set Gx  [ expr ${dkx}/(${gammabar}*${dtx})]
    set dt  [ expr ${Gx}/$tana]
    set dty [ expr sqrt(4*${dky}/(${gammabar}*${tana}))]
    set entries($w,echosp) [expr round(10000000*(($entries($w,outsize_nx)-1)*$dtx+2*$dt))/10000000.0]
  return 0
}

proc possum:pulsecheck { w obvol mrpar seqtype te tr trslc outsize_nx outsize_ny outsize_nz outsize_dx outsize_dy outsize_dz fov_x fov_y fov_z numvol zstart gap bw readdir phasedir slcdir plus maxG riseT b0f mot act1 act2 out numproc segs slcprof cover flipangle {onlyCheck 0}} {
    global entries FSLDIR
    set chkwarning ""
    set chkmessage ""
    if { $seqtype == "epi" } {
	    set dx [ expr $outsize_dx * 0.001 ]
	    set dy [ expr $outsize_dy * 0.001 ]
	    set dz [ expr $outsize_dz * 0.001 ]
	    set zs [ expr $zstart * 0.001 ]
	    # checks that the Pulse Sequence parameters are appropriate
	    set gammabar 42580000
	    set Gz 7.128*1e-03
	    set tana [expr ${maxG}/${riseT} ]
	    set dtz [expr $Gz/$tana ]
	    set rft [expr 4*0.001 ]
	    set dtz1 [expr sqrt($Gz*($dtz+$rft)*2/$tana)]
	    set Gz1 [expr $dtz1*$tana/2]
	    set TA [expr $rft/2+$dtz+$dtz1 ]
	    set dtx [ expr 1.0/$bw]
	    set dkx [ expr 1.0/(${outsize_nx}*${dx})]
	    set dky [ expr 1.0/(${outsize_ny}*${dy})]
	    set Gx  [ expr ${dkx}/(${gammabar}*${dtx})]
	    set dt  [ expr ${Gx}/$tana]
	    set dty [ expr sqrt(4*${dky}/(${gammabar}*${tana}))]
	    set Gy  [expr ${dty}*${tana}/2]
	    set dtx1 [expr sqrt(${Gx}*($dt+${outsize_nx}*$dtx)*2/$tana)] 
	    set Gx1 [ expr $dtx1*$tana/2]
	    set dty1 [expr sqrt(${outsize_ny}/2)*$dty]
	    set Gy1 [expr $dty1*$tana/2]
	    # Takes into account that the kspace can be partial 
	    if { $cover == 100 } {
		set bottomkspace [expr $outsize_ny/2]
	    } else {
		set bottomkspace [expr ($cover-50)*$outsize_ny/100.0 ]
	    }
	    set TEl [expr $bottomkspace*(2*$dt+($outsize_nx-1)*$dtx)+($dt+${outsize_nx}/2*$dtx)]
	    set TEr [expr (${outsize_ny}/2-1)*(2*$dt+(${outsize_nx}-1)*$dtx)+($dt+(${outsize_nx}/2-1)*$dtx)]
	    set TD [expr $te - $TEl ]
	    set TC [expr $TD - $dtx1 ]
	    set TB [expr $TC - $dty1 ]
	    set TF [expr $te + $TEr ]
	    set tcrush [expr 100 * $riseT]
	    set TG [expr $TF + 2 * $riseT + $tcrush ]
	    set tmpSLC [expr $trslc*$outsize_nz ]
	    set count 0
	    set w0 ".dialog[incr count]"
	    while { [ winfo exists $w0 ] } {
		set w0 ".dialog[incr count]"
	    }
	    #making sure that the scan FOV does not exceed the given FOV in slcseldir
	    set maxFOV [expr $zstart + $fov_z]
	    set dimX [ expr $entries($w,vcX) * $entries($w,inNx) ] 
	    set dimY [ expr $entries($w,vcY) * $entries($w,inNy) ]
	    set dimZ [ expr $entries($w,vcZ) * $entries($w,inNz) ]
	    set dim $dimZ
	    if { $slcdir == "x" } {
		set dim $dimX
	    } elseif { $slcdir == "y" } { 
		set dim $dimY
	    } else {
		set dim $dimZ
	    }
	    if { $maxFOV > $dim } {
		set newZstart [ expr $zstart - $maxFOV + $dim ]
		if { $newZstart < 0 } {
		    set newZstart 0
		    set newFOV $dim
		    set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
		    set entries($w,pulsechecktest) 0 
		    set chkmessage "$chkmessage The selected slices do not fit within the input object. Try changing the starting slice position to $newZstart and the field of view (Z) to less or equal to $newFOV.\n"
		} else {
		    set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
		    set entries($w,pulsechecktest) 0 
		    set chkmessage "$chkmessage The selected slices do not fit within the input object. Try changing the starting slice position to $newZstart.\n"
		}  
	    } elseif {$tr < $tmpSLC} {
		set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
		set entries($w,pulsechecktest) 0 
		set chkmessage "$chkmessage Try changing the TR to $tmpSLC.\n"
	    } elseif { $TB < 0 } {
		set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
		set entries($w,pulsechecktest) 0 
		set tmpB [expr $te -$TA]
		set tmpC [expr ($outsize_ny*($outsize_nx-1)+$outsize_nx)/2]
		set tmpA [expr $dkx*($outsize_ny+1)/($tana*$gammabar)]
		set tmpABC [expr $tmpB*$tmpB-4*$tmpA*$tmpC]
		set newTE [expr $TEl+$dtx1+$dty1+$TA ]
		if { $tmpABC < 0 } {
		    set chkmessage "$chkmessage Try changing the TE to greater than $newTE s.\n"
		} else {
		    set BW1 [expr ($tmpB-sqrt($tmpABC))/(2*$tmpA)]
		    set BW2 [expr ($tmpB+sqrt($tmpABC))/(2*$tmpA)]
		    if { $BW1 < 0 || $BW2 < 0 } {
			set chkmessage "$chkmessage Try changing the TE to greater than $newTE s. \n"
		    } else {
			set chkmessage "$chkmessage Try changing the BW to a value between $BW1 and $BW2\nOR\nTry changing the TE to greater than $newTE s.\n"
		    }
		}
	    } elseif { $TG > $trslc } {
		set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
		set entries($w,pulsechecktest) 0 
		set tmpB [expr $trslc-$te-2*$riseT-$tcrush]
		set tmpC [expr ($outsize_ny/2-1)*($outsize_nx-1)+$outsize_nx/2-1]
		set tmpA [expr $dkx*($outsize_ny-1)/($tana*$gammabar)]
		set tmpABC [expr $tmpB*$tmpB-4*$tmpA*$tmpC]
		set newTRslc [expr $te+$TEr+2*$riseT+$tcrush]
		set newTRval [expr $newTRslc/0.999*$outsize_nz]
		if { $tmpABC < 0 } {
		    set newTRslc [expr $te+$TEr+2*$riseT+$tcrush]
		    set newTRval [expr $newTRslc/0.999*$outsize_nz]
		    if { $entries($w,autotrslc) == 1 } {
			set chkmessage "$chkmessage Try changing TR to be greater than $newTRval s.\n"
		    } else {
			set chkmessage "$chkmessage Try changing the TRslc value to greater than $newTRslc s.\n"
		    }
		} else {
		    set BW1 [expr ($tmpB-sqrt($tmpABC))/(2*$tmpA)]
		    set BW2 [expr ($tmpB+sqrt($tmpABC))/(2*$tmpA)]
		    if { $BW1 < 0 || $BW2 < 0 } {
			if { $entries($w,autotrslc) == 1 } {
			    set chkmessage "$chkmessage Try changing TR to be greater than $newTRval s.\n"
			} else {
			    set chkmessage "$chkmessage Try changing the TRslc value to greater than $newTRslc s.\n"
			}
		    } else {
			if { $entries($w,autotrslc) == 1 } {
			    set chkmessage "$chkmessage Try changing the BW to a value between $BW1 and $BW2\nOR\nTry changing TR to be greater than $newTRval s.\n"
			} else {
			    set chkmessage "$chkmessage Try changing the BW to a value between $BW1 and $BW2\nOR\nTry changing the TRslc value to greater than $newTRslc s.\n"
			}
		    }
		}
	    } else {
		if { $readdir == $phasedir || $readdir == $slcdir || $phasedir == $slcdir } {
		    set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
		    set entries($w,pulsechecktest) 0 
		    set chkmessage "$chkmessage Read-, phase-, and slice- directions must be different.\n"
		} else {
		    set chkmessage "$chkmessage All is well with the pulse sequence set up.\n"
		    set entries($w,pulsechecktest) 1
		}
	    }
	if { $entries($w,pulsechecktest) != 1 || $onlyCheck == 1 } {
	    	toplevel $w0
	    	wm title $w0 "Pulse Sequence Check"
	    	wm iconname $w0 "PulseCheck"
	    	wm iconbitmap $w0 @${FSLDIR}/tcl/fmrib.xbm

		label $w0.msg1 -text "$chkwarning" -font {Helvetica 12 bold} -foreground red
		label $w0.msg2 -text "$chkmessage" -font {Helvetica 12 bold} -foreground black
		button $w0.cancel -command "destroy $w0" -text "Dismiss"
		pack $w0.msg1 $w0.msg2 $w0.cancel -in $w0
	}

    } elseif { $seqtype == "ge" } {
	    set dx [ expr $outsize_dx * 0.001 ]
	    set dy [ expr $outsize_dy * 0.001 ]
	    set dz [ expr $outsize_dz * 0.001 ]
	    set zs [ expr $zstart * 0.001 ]
	    # checks that the Pulse Sequence parameters are appropriate
	    set gammabar 42580000
	    set Gz 7.128*1e-03
	    set tana [expr ${maxG}/${riseT} ]
	    set dtz [expr $Gz/$tana ]
	    set rft [expr 4*0.001 ]
####
	    set dtz1 [expr sqrt($Gz*($dtz+$rft)*2/$tana)]
	    set Gz1 [expr $dtz1*$tana/2]
	    set TA [expr $rft/2+$dtz+$dtz1 ]
	    set dtx [ expr 1.0/$bw]
	    set dkx [ expr 1.0/(${outsize_nx}*${dx})]
	    set dky [ expr 1.0/(${outsize_ny}*${dy})]
	    set Gx  [ expr ${dkx}/(${gammabar}*${dtx})]
	    set dt  [ expr ${Gx}/$tana]
	    set dty [ expr sqrt(4*${dky}/(${gammabar}*${tana}))]
	    set Gy  [expr ${dty}*${tana}/2]
	    set dtx1 [expr sqrt(${Gx}*($dt+${outsize_nx}*$dtx)*2/$tana)] 
	    set Gx1 [ expr $dtx1*$tana/2]
	    set dty1 [expr sqrt(${outsize_ny}/2)*$dty]
	    set Gy1 [expr $dty1*$tana/2]
	    # Takes into account that the kspace can be partial 
	    if { $cover == 100 } {
		set bottomkspace [expr $outsize_ny/2]
	    } else {
		set bottomkspace [expr ($cover-50)*$outsize_ny/100.0 ]
	    }
	    set TEl [expr sqrt($outsize_ny/2)*$dty]
	    set TEr [expr $dt+($outsize_nx/2-1)*$dtx]
	    set TD [expr $te - $TEl ]
	    set TC [expr $TD - $dtx1 ]
	    set TB [expr $TC - $dty1 ]
	    set TF [expr $te + $TEr ]
	    set tcrush [expr 100 * $riseT]
	    set TG [expr $TF + 2 * $riseT + $tcrush ]
	    set tmpSLC [expr $trslc*$outsize_nz ]
#####

	    set count 0
	    set w0 ".dialog[incr count]"
	    while { [ winfo exists $w0 ] } {
		set w0 ".dialog[incr count]"
	    }
	    #making sure that the scan FOV does not exceed the given FOV in slcseldir
	    set maxFOV [expr $zstart + $fov_z]
	    set dimX [ expr $entries($w,vcX) * $entries($w,inNx) ] 
	    set dimY [ expr $entries($w,vcY) * $entries($w,inNy) ]
	    set dimZ [ expr $entries($w,vcZ) * $entries($w,inNz) ]
	    set dim $dimZ
	    if { $slcdir == "x" } {
		set dim $dimX
	    } elseif { $slcdir == "y" } { 
		set dim $dimY
	    } else {
		set dim $dimZ
	    }

	    if { $maxFOV > $dim } {
		set newZstart [ expr $zstart - $maxFOV + $dim ]
		if { $newZstart < 0 } {
		    set newZstart 0
		    set newFOV $dim
		    set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
		    set entries($w,pulsechecktest) 0 
		    set chkmessage "$chkmessage The selected slices do not fit within the input object. Try changing the starting slice position to $newZstart and the field of view (Z) to less or equal to $newFOV.\n"
		} else {
		    set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
		    set entries($w,pulsechecktest) 0 
		    set chkmessage "$chkmessage The selected slices do not fit within the input object. Try changing the starting slice position to $newZstart.\n"
		}  
	    } elseif { $TD<0 || $TB<0 || $TA<0 || $TC<0 } {
		    set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
		    set entries($w,pulsechecktest) 0 
		    set chkmessage "$chkmessage TE is not long enough to accomodate for the resX, resY, and the BW\n"
	    } elseif { $TF>$tr || $TG> $tr } {
		    set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
		    set entries($w,pulsechecktest) 0 
		    set chkmessage "$chkmessage TR is not long enough.\n"
	    } else {
		    if { $readdir == $phasedir || $readdir == $slcdir || $phasedir == $slcdir } {
		    set chkwarning "$chkwarning WARNING: Pulse sequence parameters did not pass the consistency check.\n\n"
		    set entries($w,pulsechecktest) 0 
		    set chkmessage "$chkmessage Read-, phase-, and slice- directions must be different.\n"
		} else {
		    set chkmessage "$chkmessage All is well with the pulse sequence set up.\n"
		    set entries($w,pulsechecktest) 1
	    }
	}
	if { $entries($w,pulsechecktest) != 1 || $onlyCheck == 1 } {
		toplevel $w0
		wm title $w0 "Pulse Sequence Check"
		wm iconname $w0 "PulseCheck"
		wm iconbitmap $w0 @${FSLDIR}/tcl/fmrib.xbm

		label $w0.msg1 -text "$chkwarning" -font {Helvetica 12 bold} -foreground red
		label $w0.msg2 -text "$chkmessage" -font {Helvetica 12 bold} -foreground black
		button $w0.cancel -command "destroy $w0" -text "Dismiss"
		pack $w0.msg1 $w0.msg2 $w0.cancel -in $w0
	}
    }
    return 0
}


proc possum:write { w filename } {
    global entries FSLDIR
set channel [ open ${filename} "w" ]
puts $channel "set entries(\$w,act1) \"$entries($w,act1)\""
puts $channel "set entries(\$w,act2) \"$entries($w,act2)\""
puts $channel "set entries(\$w,activ_yn) \"$entries($w,activ_yn)\""
puts $channel "set entries(\$w,autotrslc) \"$entries($w,autotrslc)\""
puts $channel "set entries(\$w,b0f) \"$entries($w,b0f)\""
puts $channel "set entries(\$w,b0fieldstrength) \"$entries($w,b0fieldstrength)\""
puts $channel "set entries(\$w,b0ftime) \"$entries($w,b0ftime)\""
puts $channel "set entries(\$w,b0ftimecourse) \"$entries($w,b0ftimecourse)\""
puts $channel "set entries(\$w,b0inh_yn) \"$entries($w,b0inh_yn)\""
puts $channel "set entries(\$w,b0inhtime_yn) \"$entries($w,b0inhtime_yn)\""
puts $channel "set entries(\$w,b0strength) \"$entries($w,b0strength)\""
puts $channel "set entries(\$w,b0units) \"$entries($w,b0units)\""
puts $channel "set entries(\$w,b0unitstime) \"$entries($w,b0unitstime)\""
puts $channel "set entries(\$w,bw) \"$entries($w,bw)\""
puts $channel "set entries(\$w,comptime) \"$entries($w,comptime)\""
puts $channel "set entries(\$w,cover) \"$entries($w,cover)\""
puts $channel "set entries(\$w,flipangle) \"$entries($w,flipangle)\""
puts $channel "set entries(\$w,echosp) \"$entries($w,echosp)\""
puts $channel "set entries(\$w,fov_x) \"$entries($w,fov_x)\""
puts $channel "set entries(\$w,fov_y) \"$entries($w,fov_y)\""
puts $channel "set entries(\$w,fov_z) \"$entries($w,fov_z)\""
puts $channel "set entries(\$w,gap) \"$entries($w,gap)\""
puts $channel "set entries(\$w,inNt) \"$entries($w,inNt)\""
puts $channel "set entries(\$w,inNtact) \"$entries($w,inNtact)\""
puts $channel "set entries(\$w,inNtb0) \"$entries($w,inNtb0)\""
puts $channel "set entries(\$w,inNx) \"$entries($w,inNx)\""
puts $channel "set entries(\$w,inNxact) \"$entries($w,inNxact)\""
puts $channel "set entries(\$w,inNxb0) \"$entries($w,inNxb0)\""
puts $channel "set entries(\$w,inNxb0time) \"$entries($w,inNxb0time)\""
puts $channel "set entries(\$w,inNy) \"$entries($w,inNy)\""
puts $channel "set entries(\$w,inNyact) \"$entries($w,inNyact)\""
puts $channel "set entries(\$w,inNyb0) \"$entries($w,inNyb0)\""
puts $channel "set entries(\$w,inNyb0time) \"$entries($w,inNyb0time)\""
puts $channel "set entries(\$w,inNz) \"$entries($w,inNz)\""
puts $channel "set entries(\$w,inNzact) \"$entries($w,inNzact)\""
puts $channel "set entries(\$w,inNzb0) \"$entries($w,inNzb0)\""
puts $channel "set entries(\$w,inNzb0time) \"$entries($w,inNzb0time)\""
puts $channel "set entries(\$w,maxG) \"$entries($w,maxG)\""
puts $channel "set entries(\$w,mot) \"$entries($w,mot)\""
puts $channel "set entries(\$w,motion_yn) \"$entries($w,motion_yn)\""
puts $channel "set entries(\$w,mrpar) \"$entries($w,mrpar)\""
puts $channel "set entries(\$w,noise_yn) \"$entries($w,noise_yn)\""
puts $channel "set entries(\$w,noisesigma) \"$entries($w,noisesigma)\""
puts $channel "set entries(\$w,noisesnr) \"$entries($w,noisesnr)\""
puts $channel "set entries(\$w,noiseunits) \"$entries($w,noiseunits)\""
puts $channel "set entries(\$w,numproc) \"$entries($w,numproc)\""
puts $channel "set entries(\$w,segs) \"$entries($w,segs)\""
puts $channel "set entries(\$w,numvol) \"$entries($w,numvol)\""
puts $channel "set entries(\$w,obvol) \"$entries($w,obvol)\""
puts $channel "set entries(\$w,out) \"$entries($w,out)\""
puts $channel "set entries(\$w,outsize_dx) \"$entries($w,outsize_dx)\""
puts $channel "set entries(\$w,outsize_dy) \"$entries($w,outsize_dy)\""
puts $channel "set entries(\$w,outsize_dz) \"$entries($w,outsize_dz)\""
puts $channel "set entries(\$w,outsize_nx) \"$entries($w,outsize_nx)\""
puts $channel "set entries(\$w,outsize_ny) \"$entries($w,outsize_ny)\""
puts $channel "set entries(\$w,outsize_nz) \"$entries($w,outsize_nz)\""
puts $channel "set entries(\$w,phencode) \"$entries($w,phencode)\""
puts $channel "set entries(\$w,plus) \"$entries($w,plus)\""
puts $channel "set entries(\$w,proctime) \"$entries($w,proctime)\""
puts $channel "set entries(\$w,pulsechecktest) \"$entries($w,pulsechecktest)\""
puts $channel "set entries(\$w,readgrad) \"$entries($w,readgrad)\""
puts $channel "set entries(\$w,riseT) \"$entries($w,riseT)\""
puts $channel "set entries(\$w,slcprof) \"$entries($w,slcprof)\""
puts $channel "set entries(\$w,slcselect) \"$entries($w,slcselect)\""
puts $channel "set entries(\$w,te) \"$entries($w,te)\""
puts $channel "set entries(\$w,tr) \"$entries($w,tr)\""
puts $channel "set entries(\$w,trslc) \"$entries($w,trslc)\""
puts $channel "set entries(\$w,vcX) \"$entries($w,vcX)\""
puts $channel "set entries(\$w,vcXact) \"$entries($w,vcXact)\""
puts $channel "set entries(\$w,vcXb0) \"$entries($w,vcXb0)\""
puts $channel "set entries(\$w,vcXb0time) \"$entries($w,vcXb0time)\""
puts $channel "set entries(\$w,vcY) \"$entries($w,vcY)\""
puts $channel "set entries(\$w,vcYact) \"$entries($w,vcYact)\""
puts $channel "set entries(\$w,vcYb0) \"$entries($w,vcYb0)\""
puts $channel "set entries(\$w,vcYb0time) \"$entries($w,vcYb0time)\""
puts $channel "set entries(\$w,vcZ) \"$entries($w,vcZ)\""
puts $channel "set entries(\$w,vcZact) \"$entries($w,vcZact)\""
puts $channel "set entries(\$w,vcZb0) \"$entries($w,vcZb0)\""
puts $channel "set entries(\$w,vcZb0time) \"$entries($w,vcZb0time)\""
puts $channel "set entries(\$w,zstart) \"$entries($w,zstart)\""
puts $channel "set entries(\$w,seqtype) \"$entries($w,seqtype)\""
#puts $channel "set entries(\$w,rfavg_yn) \"$entries($w,rfavg_yn)\""
puts $channel "set entries(\$w,custompulse_yn) \"$entries($w,custompulse_yn)\""
puts $channel "set entries(\$w,cuspulse) \"$entries($w,cuspulse)\""
puts $channel "set entries(\$w,slcsampfactor) \"$entries($w,slcsampfactor)\""
close $channel
}

proc possum:load { w filename } {
    global entries FSLDIR
    source ${filename}
    possum:updateFOV $w
    possum:updateVSIZE $w
    possum:updateTRSLC $w
    possum:updateb0field $w
    possum:updateb0fieldinh $w
    possum:updateb0fieldinhtime $w
    possum:updatemotion $w
    possum:updateactivation $w
    possum:updatenoise $w
    possum:updatenoiseunits $w
    possum:updatecomptime $w 
    possum:updateechosp $w
    possum:custompulse $w
#    possum:updaterfavg $w
}
proc possum:procmakedir { w comptime obvol mrpar seqtype te tr trslc outsize_nx outsize_ny outsize_nz outsize_dx outsize_dy outsize_dz fov_x fov_y fov_z numvol zstart gap bw readdir phasedir slcdir plus maxG riseT b0f b0fieldstrength b0units b0extra b0timecourse b0extraunits mot act1 act2 out numproc segs slcprof cover flipangle slcsampfactor } {
    global entries FSLDIR POSSUMDIR
    set dx [ expr $outsize_dx * 0.001 ]
    set dy [ expr $outsize_dy * 0.001 ]
    set dz [ expr $outsize_dz * 0.001 ]
    set zs [ expr $zstart * 0.001 ]
    set gap [ expr $gap * 0.001 ]
  
    catch { exec sh -c "${FSLDIR}/bin/fslval $obvol dim1" } dim1
    catch { exec sh -c "${FSLDIR}/bin/fslval $obvol dim2" } dim2
    catch { exec sh -c "${FSLDIR}/bin/fslval $obvol dim3" } dim3
    catch { exec sh -c "${FSLDIR}/bin/fslval $obvol dim4" } dim4
    catch { exec sh -c "${FSLDIR}/bin/fslval $obvol pixdim1" } pixdim1
    catch { exec sh -c "${FSLDIR}/bin/fslval $obvol pixdim2" } pixdim2
    catch { exec sh -c "${FSLDIR}/bin/fslval $obvol pixdim3" } pixdim3
    catch { exec sh -c "${FSLDIR}/bin/fslval $obvol pixdim4" } pixdim4
    set newdim1 $dim1
    set newdim2 $dim2
    set newdim3 $dim3
    set newpixdim1 $pixdim1
    set newpixdim2 $pixdim2
    set newpixdim3 $pixdim3
    if { $slcdir == "x" } {
	set newpixdim1 [ expr $outsize_dx / $slcsampfactor ]
	set newdim1 [ expr round( $dim1 * $pixdim1 / $newpixdim1 ) ]
    }
    if { $slcdir == "y" } {
	set newpixdim2 [ expr $outsize_dy / $slcsampfactor ]
	set newdim2 [ expr round( $dim2 * $pixdim2 / $newpixdim2 ) ]
    }
    if { $slcdir == "z" } {
	set newpixdim3 [ expr $outsize_dz / $slcsampfactor ]
	set newdim3 [ expr round( $dim3 * $pixdim3 / $newpixdim3 ) ]
    }
    puts "${FSLDIR}/bin/fslcreatehd $newdim1 $newdim2 $newdim3 $dim4 $newpixdim1 $newpixdim2 $newpixdim3 $pixdim4 0 0 0 16 $out/brainref"
    puts "${FSLDIR}/bin/flirt -in $obvol -ref $out/brainref -applyxfm -out $out/brain"
    catch { exec sh -c "${FSLDIR}/bin/fslcreatehd $dim1 $dim2 $newdim3 $dim4 $pixdim1 $pixdim2 $newpixdim3 $pixdim4 0 0 0 16 $out/brainref" } oval
    catch { exec sh -c "${FSLDIR}/bin/flirt -in $obvol -ref $out/brainref -applyxfm -out $out/brain" } oval
    catch { exec sh -c "${FSLDIR}/bin/imrm $out/brainref" } oval
#    catch { exec sh -c "${FSLDIR}/bin/imcp $obvol $out/brain" } oval
    catch { exec sh -c "${FSLDIR}/bin/imrm $out/braintmp" } oval
    catch { exec sh -c "cp $mrpar $out/MRpar" } oval
    catch { exec sh -c "cp $slcprof $out/slcprof" } oval
    catch { exec sh -c "cp $mot $out/motion" } oval
#    catch { exec sh -c "${FSLDIR}/bin/imcp $act1 $out/T2" } oval
    catch { exec sh -c "${FSLDIR}/bin/flirt -in $act1 -ref $out/brainref -applyxfm -out $out/T2" } oval
    catch { exec sh -c "cp $act2 $out/T2timecourse" } oval

    if { $b0f != "" && $mot == "${FSLDIR}/data/possum/zeromotion" } {
	catch { exec sh -c "${FSLDIR}/bin/flirt -in ${b0f} -ref $out/brainref -applyxfm -out $out/b0newref" } oval
	catch { exec sh -c "${FSLDIR}/bin/fslroi $out/b0newref ${out}/b0z_dz.nii.gz 0 1" } oval
	#	catch { exec sh -c "${FSLDIR}/bin/fslroi ${b0f} ${out}/b0z_dz.nii.gz 0 1" } oval
       if { $b0units == "ppm" } {
	  catch { exec sh -c "${FSLDIR}/bin/fslmaths $out/b0z_dz -mul $b0fieldstrength -div 1000000 $out/b0z_dz" } oval
       }
    }

    if { $b0f != "" && $mot != "${FSLDIR}/data/possum/zeromotion" } {
	catch { exec sh -c "${FSLDIR}/bin/flirt -in ${b0f} -ref $out/brainref -applyxfm -out $out/b0inh" } oval
        if { $b0units == "ppm" } {
	    catch { exec sh -c "${FSLDIR}/bin/fslmaths $out/b0inh -mul $b0fieldstrength -div 1000000 $out/b0inh" } oval
	}
##check dim
	if { [exec sh -c "${FSLDIR}/bin/fslval $out/b0inh dim4"] == 9 } {
		catch { exec sh -c "${FSLDIR}/bin/fslroi $out/b0inh ${out}/b0x_dx.nii.gz 8 1" } oval
		catch { exec sh -c "${FSLDIR}/bin/fslroi $out/b0inh ${out}/b0x_dy.nii.gz 7 1" } oval
		catch { exec sh -c "${FSLDIR}/bin/fslroi $out/b0inh ${out}/b0x_dz.nii.gz 6 1" } oval
		catch { exec sh -c "${FSLDIR}/bin/fslroi $out/b0inh ${out}/b0y_dx.nii.gz 5 1" } oval
		catch { exec sh -c "${FSLDIR}/bin/fslroi $out/b0inh ${out}/b0y_dy.nii.gz 4 1" } oval
		catch { exec sh -c "${FSLDIR}/bin/fslroi $out/b0inh ${out}/b0y_dz.nii.gz 3 1" } oval
		catch { exec sh -c "${FSLDIR}/bin/fslroi $out/b0inh ${out}/b0z_dx.nii.gz 2 1" } oval
		catch { exec sh -c "${FSLDIR}/bin/fslroi $out/b0inh ${out}/b0z_dy.nii.gz 1 1" } oval
		catch { exec sh -c "${FSLDIR}/bin/fslroi $out/b0inh ${out}/b0z_dz.nii.gz 0 1" } oval
		catch { exec sh -c "${FSLDIR}/bin/imrm $out/b0inh" } oval
	} else {
		puts "Error: B0 field size is incorrect!"
		return
	}
    }
    if { $b0extra != "" && $mot == "${FSLDIR}/data/possum/zeromotion" } {
	catch { exec sh -c "${FSLDIR}/bin/flirt -in $b0extra -ref $out/brainref -applyxfm -out $out/b0extra" } oval
	catch { exec sh -c "cp $b0timecourse $out/b0timecourse" } oval
	if { $b0extraunits == "ppm" } {
	    catch { exec sh -c "${FSLDIR}/bin/fslmaths $out/b0extra -mul $b0fieldstrength -div 1000000 $out/b0extra" } oval
	    catch { exec sh -c "cp $b0timecourse $out/b0timecourse" } oval
	}
    }
    if { $b0extra != "" && $mot != "${FSLDIR}/data/possum/zeromotion" } {
       puts "Warning: At the moment B0 field changing in time can not be simulated while the object is moving. This will be implemented into POSSUM at a later stage."
       return
    }
    if { $entries($w,custompulse_yn) == 0 } {
	    set seq $seqtype
	    if { $seq == "epi" } {
		    set pulsecom  "${POSSUMDIR}/bin/pulse -i $out/brain -o ${out}/pulse --te=${te} --tr=${tr} --trslc=${trslc} --nx=${outsize_nx} --ny=${outsize_ny} --dx=${dx} --dy=${dy} --maxG=${maxG} --riset=${riseT} --bw=${bw} --numvol=${numvol} --numslc=${outsize_nz} --slcthk=${dz} --zstart=${zs} --seq=${seq} --slcdir=${slcdir}${plus} --readdir=${readdir}$entries($w,pluss) --phasedir=${phasedir}$entries($w,pluss) --gap=$gap -v --cover=$cover --angle=$flipangle"
	    } elseif { $seq == "ge" } {
		    set pulsecom  "${POSSUMDIR}/bin/pulse -i $out/brain -o ${out}/pulse --te=${te} --tr=${tr} --nx=${outsize_nx} --ny=${outsize_ny} --dx=${dx} --dy=${dy} --maxG=${maxG} --riset=${riseT} --bw=${bw} --numvol=${numvol} --numslc=${outsize_nz} --slcthk=${dz} --zstart=${zs} --seq=${seq} --slcdir=${slcdir}${plus} --readdir=${readdir}$entries($w,pluss) --phasedir=${phasedir}$entries($w,pluss) --gap=$gap -v --cover=$cover --angle=$flipangle"
	    }
	    catch { exec sh -c "echo $pulsecom >> $out/pulse.com" } oval
	    catch { exec sh -c "echo $pulsecom >> $out/possum.log" } oval
	    Possum:pulsecheck $w
	    if { $entries($w,pulsechecktest) == 0 } {
		return
	    }
	    puts $pulsecom
	    fsl:exec "$pulsecom >> $out/possum.log 2>&1" 
    } else {
	    set pulsebasename $entries($w,cuspulse)
	    set filelist [list "" ".info" ".posx" ".posy" ".posz" ".readme"]
	    foreach ext $filelist {
		if { [file exists $pulsebasename$ext ] != 0 } {
			catch { exec sh -c "cp -rv ${pulsebasename}${ext} ${out}/pulse${ext} >> possum.log" } copystatus
		} else {
			puts "File '$pulsebasename$ext' not found!"
			return 0
		}
	    }
   }
   return 0
}

proc possum:proc { w comptime obvol mrpar seqtype te tr trslc outsize_nx outsize_ny outsize_nz outsize_dx outsize_dy outsize_dz fov_x fov_y fov_z numvol zstart gap bw readdir phasedir slcdir plus maxG riseT b0f b0fieldstrength b0units b0extra b0timecourse b0extraunits mot act1 act2 out numproc segs slcprof cover flipangle slcsampfactor } {
    global entries FSLDIR POSSUMDIR

    possum:procmakedir $w $comptime $obvol $mrpar $seqtype $te $tr $trslc $outsize_nx $outsize_ny $outsize_nz $outsize_dx $outsize_dy $outsize_dz $fov_x $fov_y $fov_z $numvol $zstart $gap $bw $readdir $phasedir $slcdir $plus $maxG $riseT $b0f $b0fieldstrength $b0units $b0extra $b0timecourse $b0extraunits $mot $act1 $act2 $out $numproc $segs $slcprof $cover $flipangle $slcsampfactor

    if { $entries($w,ctt) != 0 } {
       set comptime $entries($w,ctt)
    }
    # Execute possum
#    if { $entries($w,rfavg_yn) == 0 } {
    	set possumcom  "${POSSUMDIR}/bin/possumX $out -n $numproc -t $comptime -s $segs"
#    } else {
#	set possumcom  "${POSSUMDIR}/bin/possumX $out -n $numproc -t $comptime -s $segs -r"
#    }
    catch { exec sh -c "echo $possumcom >> $out/possum.log" } oval
    fsl:exec "$possumcom  >> $out/possum.log 2>&1"
    return 0
}

wm withdraw .
possum .rename
tkwait window .rename

