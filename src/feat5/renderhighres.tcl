#{{{ copyright and setup

#   renderhighres - render stats onto high res structurals
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 1999-2005 University of Oxford
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
#bind Button <1> "focus %W ; [bind Button <1>] "
#tk_focusFollowsMouse
set VARS(history) {}

#}}}
#{{{ feat5:renderhighres

proc feat5:renderhighres { w } {
    global PXHOME FSLDIR VARS argc argv PWD entries rendervars

    #{{{ setup main window

toplevel $w

wm title      $w "Upsample+overlay FEAT stats"
wm iconname   $w "FEAT upsample"
wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

frame $w.f
pack $w.f -in $w -side top



FileEntry $w.f.featdir -textvariable entries($w,1) -label "Select a FEAT directory" -title "Select a FEAT directory"  -width 20 -filedialog directory -filetypes "*.feat" -dirasfile "design.fsf" 
#FSLFileEntry $w.f.featdir -varia entries($w,1)-pattern "*.feat" -directory "~" -label "Select a FEAT directory" -title "Select a FEAT directory" -width 20 -dirasfile "design.fsf" -filterhist VARS(history)

frame $w.f.space
label $w.f.space.label -text "Space to upsample to: "
set rendervars(space) standard
radiobutton $w.f.space.highres  -text "main structural"  -value highres  -variable rendervars(space)
radiobutton $w.f.space.standard -text "standard"         -value standard -variable rendervars(space)
pack $w.f.space.label $w.f.space.highres $w.f.space.standard -in $w.f.space -side left -anchor w

frame $w.f.background
label $w.f.background.label -text "Background image: "
set rendervars(background) highres
radiobutton $w.f.background.highres  -text "main structural"  -value highres  -variable rendervars(background)
radiobutton $w.f.background.standard -text "standard"         -value standard -variable rendervars(background)
pack $w.f.background.label $w.f.background.highres $w.f.background.standard -in $w.f.background -side left -anchor w

frame $w.f.autothresh
label $w.f.autothresh.label -text "Colour-code using existing post-threshold range"
set rendervars(autothresh) 1
checkbutton $w.f.autothresh.yn -variable rendervars(autothresh) -command "feat5:autothresh_update $w"
pack $w.f.autothresh.label $w.f.autothresh.yn -in $w.f.autothresh -side left -anchor w

frame $w.f.minmax
set rendervars(zmin) 1


LabelSpinBox $w.f.zmin -label "  Colour-code range:  Min"  -textvariable rendervars(zmin) -range {0 10000 1 } -width 4
#tixControl $w.f.zmin -label "  Colour-code range:  Min" -variable rendervars(zmin) -step 1 -min 0 -selectmode immediate -options { entry.width 4 }
set rendervars(zmax) 15

LabelSpinBox $w.f.zmax -label "Max"  -textvariable rendervars(zmax) -range {0 10000 1 } -width 4
#tixControl $w.f.zmax -label "Max" -variable rendervars(zmax) -step 1 -min 0 -selectmode immediate -options { entry.width 4 }
pack $w.f.zmin $w.f.zmax -in $w.f.minmax -side left -anchor w -padx 5

pack $w.f.featdir $w.f.space $w.f.background $w.f.autothresh -in $w.f -padx 5 -pady 5 -anchor w

#}}}
    #{{{ Button Frame

frame $w.btns
frame $w.btns.b -relief raised -borderwidth 1
    
button $w.apply -command "feat5:renderhighres_proc $w" -text "Go"

button $w.cancel -command "destroy $w" -text "Exit"

pack $w.btns.b -side bottom -fill x -padx 3 -pady 3
pack $w.apply $w.cancel -in $w.btns.b -side left -expand yes -padx 3 -pady 3 -fill y

pack $w.btns -in $w -side bottom -fill x

#}}}
}

#}}}
#{{{ feat5:autothresh_update

proc feat5:autothresh_update { w } {
    global rendervars

    pack forget $w.f.minmax

    if { $rendervars(autothresh) == 0 } {
	pack $w.f.minmax -in $w.f -padx 5 -pady 5 -anchor w -after $w.f.autothresh
    }

}

#}}}
#{{{ feat5:renderhighres_proc

proc feat5:renderhighres_proc { w } {

    global FSLDIR entries rendervars

    set FD $entries($w,1)

    puts "Processing $FD"

    fsl:exec "${FSLDIR}/bin/renderhighres $FD $rendervars(space) $rendervars(background) $rendervars(autothresh) $rendervars(zmin) $rendervars(zmax)"

    puts "Done"
}

#}}}
#{{{ call GUI and wait

wm withdraw .
feat5:renderhighres .r
tkwait window .r

#}}}

