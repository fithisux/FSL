#

# FMRIB_Prepare_Fieldmap - the GUI for fmrib_prepare_fieldmap
#
# Mark Jenkinson, FMRIB Image Analysis Group
#
# Copyright (C) 2008 University of Oxford
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


# setup
source [ file dirname [ info script ] ]/fslstart.tcl

proc prepfmap { w } {

    global entries FSLDIR PWD

    # ---- Set up Frames ----
    toplevel $w
    wm title $w "FSL Prepare Fieldmap"
    wm iconname $w "FSL Prepare Fieldmap"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm
    frame $w.f

    set inFmrib [ string last fmrib.ox.ac.uk [ info hostname ] ]
    # ----- Scanner ------

    TitleFrame $w.f.scanner -text "Scanner" -relief groove 
    set lfscanner [ $w.f.scanner getframe ]

    set entries($w,scanner) "SIEMENS"
    puts $inFmrib
    if { $inFmrib == -1 } {
	radiobutton $w.f.ocmr -text "  Siemens" \
	    -variable entries($w,scanner) -value SIEMENS -anchor w -width 18
	pack $w.f.ocmr -in $lfscanner -side left -anchor w -pady 3 -padx 5 -expand yes -fill x
    } else {
	radiobutton $w.f.ocmr -text "  Siemens (OCMR & New FMRIB)" \
	    -variable entries($w,scanner) -value SIEMENS -anchor w -width 18
	radiobutton $w.f.fmrib -text "  Varian (Old FMRIB)" \
	    -variable entries($w,scanner) -value VARIAN -anchor w
	pack $w.f.ocmr $w.f.fmrib -in $lfscanner -side left -anchor w -pady 3 -padx 5 -expand yes -fill x
    }
    # ----- Inputs ------

    TitleFrame $w.f.input -text "Input Images"  -relief groove 
    set lfinput [ $w.f.input getframe ]


    set entries($w,phim) ""
    set entries($w,magim) ""
    

    FileEntry  $w.f.phase -textvariable entries($w,phim) -label "  Phase Image                             " -title "Select" -width 40 -filedialog directory  -filetypes IMAGE

    FileEntry  $w.f.mag -textvariable entries($w,magim) -label  "  Magnitude Image (Brain Extracted) " -title "Select" -width 40 -filedialog directory  -filetypes IMAGE

    pack $w.f.phase $w.f.mag -in $lfinput -side top -anchor w -pady 3 -padx 5

    # ---- Delta TE ----

    TitleFrame $w.f.dTE -text "Fieldmap Sequence Parameters" -relief groove 
    set lfdte [ $w.f.dTE getframe ]

    set entries($w,deltate) 2.5

    LabelSpinBox $w.f.deltaTE -label "  Difference of Echo Times (in milliseconds)  " -textvariable entries($w,deltate) -range {0.0 100.0 0.01 }

    pack $w.f.deltaTE -in $lfdte -side top -anchor w -padx 5 -pady 5

    # ----- Output ------

    TitleFrame $w.f.output -text "Output" -relief groove 
    set lfoutput [ $w.f.output getframe ]

    set entries($w,fmap) ""
    
    FileEntry  $w.f.fmap -textvariable entries($w,fmap) -label  "  Fieldmap image (rad/s)    " -title "Select" -width 40 -filedialog directory  -filetypes IMAGE

    pack $w.f.fmap -in $lfoutput -side top -anchor w -pady 3 -padx 5

    # ---- Global Pack -----

    pack $w.f.scanner $w.f.input $w.f.dTE $w.f.output -in $w.f -side top -anchor w -pady 0 -padx 5 -expand yes -fill x

    # ---- Button Frame ----

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.apply     -command "Prepfmap:apply $w" \
	    -text "Go" -width 5

    button $w.cancel    -command "destroy $w" \
	    -text "Exit" -width 5

    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/redirects/flirt.html" \
	    -text "Help" -width 5

    pack $w.btns.b -side bottom -fill x
#    pack $w.apply $w.cancel $w.help -in $w.btns.b
    pack $w.apply $w.cancel -in $w.btns.b \
	    -side left -expand yes -padx 3 -pady 10 -fill y
    
    pack $w.f $w.btns -expand yes -fill both

}


proc Prepfmap:apply { w } {
    global entries

    catch { prepfmap:proc $entries($w,scanner) $entries($w,phim) $entries($w,magim) $entries($w,deltate) $entries($w,fmap) }

    update idletasks
    puts "Done"
}


proc prepfmap:proc { scanner phim magim deltate fmap } {

    global FSLDIR

    # Do pop-up

    set count 0
    set w1 ".dialog[incr count]"
    while { [ winfo exists $w1 ] } {
        set w1 ".dialog[incr count]"
    }
    toplevel $w1
    wm title $w1 "Preparing Fieldmap"
    wm iconname $w1 "PrepareFieldmapOutput"
    wm iconbitmap $w1 @${FSLDIR}/tcl/fmrib.xbm
    frame $w1.sprev
    label $w1.sprev.label -text "\n    Running script ... please wait    \n\n"
    pack $w1.sprev.label -in $w1.sprev
    pack $w1.sprev -in $w1
    # force this message to popup now
    update

    # run command
    set thecommand "${FSLDIR}/bin/fsl_prepare_fieldmap $scanner $phim $magim $fmap $deltate"
    puts $thecommand
    set ret [ catch { exec sh -c "$thecommand" } otxt ]

    # show output
    pack forget $w1.sprev.label
    pack forget $w1.sprev
    pack $w1.sprev.label -in $w1.sprev
    set tottxt "Output from script:

$otxt"
    $w1.sprev.label configure -text "$tottxt"
    button $w1.cancel -command "destroy $w1" -text "Dismiss"
    pack $w1.sprev $w1.cancel -in $w1
    update
    return 0
}

# Call GUI

wm withdraw .
prepfmap .rename
tkwait window .rename

