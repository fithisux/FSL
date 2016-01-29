#

#{{{ setup

#   susan.tcl - GUI for susan smoothing
#
#   Stephen Smith, Matthew Webster FMRIB Image Analysis Group
#
#   Copyright (C) 1999-2006 University of Oxford
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

#{{{ proc susan

proc susan { w } {

    global susanvars usanentries FSLDIR PWD argc argv TN HOME tempSpin
 
    #{{{ setup main window

toplevel $w

wm title $w "SUSAN   Structure-Preserving Noise Reduction"
wm iconname $w "SUSAN"
wm iconbitmap $w @$FSLDIR/tcl/fmrib.xbm

frame $w.f

#}}}
    #{{{ setup variable defaults

set susanvars(dt)  3
set susanvars(um)  1
set susanvars(dim) 3
set susanvars(bt) 0

#}}}
    #{{{ input and output images

if { $argc > 0 && [ string length [ lindex $argv 0 ] ] > 0 } {
    set inputname [ imglob [ lindex $argv 0 ] ] 
    if { [ imtest $inputname ] } {
	if { [ string first / $inputname ] == 0 || [ string first ~ $inputname ] == 0 } {
	    set susanvars($w,input) $inputname
	} else {
	    set susanvars(input) ${PWD}/$inputname
	}
	set susanvars(output) $susanvars(input)_susan
	susan:setrange $w dummy
    }
}

FileEntry $w.f.input -textvariable susanvars(input) -label "Input image  " -title "Select the input image" -width 55 -filedialog directory -filetypes IMAGE -command "susan:setrange $w"

FileEntry $w.f.output -textvariable susanvars(output) -label "Output image" -title "Select the output image" -width 55 -filedialog directory -filetypes IMAGE

pack $w.f.input $w.f.output -in $w.f -side top -padx 5 -pady 5 -anchor w

#}}}
    #{{{ top row

frame $w.f.top

frame $w.f.op
label $w.f.oplabel -text "Dimensionality" 
optionMenu2 $w.f.opmenu susanvars(dim) 2 "2D" 3 "3D"
pack $w.f.oplabel $w.f.opmenu -in $w.f.op -side top -side left

LabelSpinBox   $w.f.bt -textvariable susanvars(bt) -label "Brightness threshold" -range {0.0 1e20 1} -width 8

LabelSpinBox   $w.f.dt -textvariable susanvars(dt) -label "Mask SD" -range {0.0 10000 1} -width 5 

pack $w.f.top -in $w.f -side top -padx 3 -pady 3 -expand yes -anchor w

pack $w.f.op $w.f.bt $w.f.dt -in $w.f.top -side left -padx 5 -pady 3 -expand yes -anchor w

#}}}
    #{{{ help text

label $w.f.label -wraplength 550 -text "Mask SD (approx HWHM) is in current length units (e.g. mm), not voxels.\nSet Mask SD to 0 for a (much faster) flat 3x3 voxels or 3x3x3 voxels mask." -fg "#303030"

pack $w.f.label -in $w.f -side top -anchor s

#}}}
    #{{{ advanced options

collapsible frame $w.f.opts -title "Advanced options"

#{{{ use median?

frame $w.f.um

label $w.f.um.label -text "Use median when no neighbourhood is found "

checkbutton $w.f.um.but -variable susanvars(um)

pack $w.f.um.label $w.f.um.but -in $w.f.um -side left

pack $w.f.um  -in $w.f.opts.b -anchor w -pady 5

#}}}
#{{{ more USANs

set susanvars(maxusans) 2
set susanvars(usans) 0
LabelSpinBox  $w.f.usans -textvariable susanvars(usans) -label  "Separate images to find USAN from " -range " 0 $susanvars(maxusans) 1 " -command " $w.f.usans.spin.e validate; susan:updateusan $w" -modifycmd "susan:updateusan $w"
#tixControl $w.f.usans -label "Separate images to find USAN from " -variable susanvars(usans) -step 1 -min 0 -max $susanvars(maxusans) -selectmode immediate -command "susan:updateusan $w"

pack $w.f.usans -in $w.f.opts.b -anchor w -pady 5

set i 1
while { $i <= $susanvars(maxusans) } {
    frame $w.f.usanentries($i)


FileEntry $w.f.ue$i -textvariable usanentries($w,$i) -label "USAN image $i" -title "Select the USAN image" -width 40 -filedialog directory -filetypes IMAGE 
#FSLFileEntry $w.f.ue$i -variable usanentries($w,$i) -pattern "IMAGE" -directory $PWD -label "USAN image $i" -title "Select the USAN image" -width 40 -filterhist VARS(history)
	
    set susanvars(ubt,$i) 1

LabelSpinBox $w.f.ubt$i -textvariable susanvars(ubt,$i) -label "Brightness threshold" -range {0.0 1e20 1} -width 5
#tixControl $w.f.ubt$i -label "Brightness threshold" -variable susanvars(ubt,$i) -step 1 -min 0 -selectmode immediate -options {entry.width 5}

    pack $w.f.ue$i $w.f.ubt$i -in $w.f.usanentries($i) -padx 3 -pady 3 -side left
    incr i 1
}

#}}}

pack $w.f.opts -in $w.f -side bottom -anchor w -pady 5

#}}}
    #{{{ button frame

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.apply -command "susan:apply $w keep" -text "Go" -width 5
    bind $w.apply <Return> {
	[winfo toplevel %W].apply invoke
    }
	    
    button $w.cancel    -command "susan:destroy $w" \
	    -text "Exit" -width 5
    bind $w.cancel <Return> {
	[winfo toplevel %W].cancel invoke
    }

    button $w.help -command "FmribWebHelp file: $FSLDIR/doc/redirects/susan.html" \
            -text "Help" -width 5
    bind $w.help <Return> {
        [winfo toplevel %W].help invoke
    }

    pack $w.btns.b -side bottom -fill x -padx 3 -pady 3
    pack $w.apply $w.cancel $w.help -in $w.btns.b \
	    -side left -expand yes -padx 3 -pady 3 -fill y

#}}}
		
    pack $w.f $w.btns -expand yes -fill both
}

#}}}
#{{{ proc susan:setrange

proc susan:setrange { w  {dummy "" }} {

    global susanvars FSLDIR

    set susanvars(input)  [ remove_ext $susanvars(input) ]
    set susanvars(output) [ remove_ext $susanvars(input) ]_susan
    
    if { ! [ catch { exec sh -c "${FSLDIR}/bin/fslstats $susanvars(input) -r" } minmax ] } {
	set min [ lindex $minmax 0 ]
	set max [ lindex $minmax 1 ]
	set susanvars(bt) [ expr ( $max - $min ) / 10.0 ]
    }
}

#}}}
#{{{ proc susan:apply

proc susan:apply { w dialog } {

    global susanvars usanentries TN

    foreach v { bt dt usans} {
	$w.f.$v.spin.e validate
    }
    susan:updateusan $w
    #{{{ process USAN entries

set usanlist ""

set i 1
while { $i <= $susanvars(usans) } {
    lappend usanlist [ remove_ext $usanentries($w,$i) ]
    lappend usanlist $susanvars(ubt,$i)
    incr i 1
}

#}}}
    susan_proc $susanvars(input) $susanvars(bt) $susanvars(output) $susanvars(dt) $susanvars(dim) $susanvars(um) $susanvars(usans) $usanlist
    
    update idletasks
    
    if {$dialog == "destroy"} {
	susan:destroy $w
    }
}

#}}}
#{{{ proc susan:destroy

proc susan:destroy { w } {
    destroy $w
}

#}}}
#{{{ updates

proc susan:updateusan { w } {
    global susanvars

    if { $susanvars(usans) == 0 } {
	pack $w.f.bt -in $w.f.top -before $w.f.dt -side left -padx 5 -pady 3 -expand yes -anchor w
    } else {
	pack forget $w.f.bt
    }

     set i 1
     while { $i <= $susanvars(maxusans) } { 
     pack forget $w.f.usanentries($i)
	incr i 1
    }

    set i 1
    while { $i <= $susanvars(usans) } {
	pack $w.f.usanentries($i) -in $w.f.opts.b -anchor w -pady 5
	incr i 1
    }
}

#}}}

#{{{ tail end

wm withdraw .
susan .rename
tkwait window .rename

#}}}
