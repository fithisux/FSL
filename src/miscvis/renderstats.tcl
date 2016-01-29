#
#   renderstats - colour rendering of two or more images
#
#   Stephen Smith and Matthew Webster, FMRIB Image Analysis Group
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

#{{{ setup

source [ file dirname [ info script ] ]/fslstart.tcl

#}}}
#{{{ render

proc render { w } {
    global entries fentries vars FSLDIR argc argv PWD USER
 
    #{{{ setup main window

toplevel $w

wm title $w "Statistics Colour Rendering"
wm iconname $w "Render"
wm iconbitmap $w @$FSLDIR/tcl/fmrib.xbm

frame $w.f

#}}}
    #{{{ Stats

set vars($w,maxstatss) 2
set vars($w,statss) 1

LabelSpinBox $w.f.statss -label " Number of stats images " -textvariable  vars($w,statss) -range "1 $vars($w,maxstatss) 1" -modifycmd "render:updatestats $w" -command "$w.f.stass.spin.e validate;render:updatestats $w"
#tixControl $w.f.statss -label " Number of stats images " -variable vars($w,statss) -step 1 -min 1 -max $vars($w,maxstatss) -selectmode immediate -command "render:updatestats $w"

pack $w.f.statss -in $w.f -anchor w -pady 5

set i 1
while { $i <= $vars($w,maxstatss) } {
    set vars($w,statsmax,$i) "max"
    frame $w.f.statsentries($i)

FileEntry $w.f.sel$i -textvariable fentries($w,$i) -label " Stats image ${i} " -title "Select" -width 35 -filetypes IMAGE
$w.f.sel$i.labf configure -width 15
#FSLFileEntry $w.f.sel$i -variable fentries($w,$i)-pattern "IMAGE"-directory $PWD-label " Stats image ${i} "-labelwidth 17 -title "Select"-width 35 -filterhist VARS(history)
    

    label $w.f.uelmin($i) -text " Min"
    entry $w.f.uemin($i) -textvariable vars($w,statsmin,$i) -width 7
    label $w.f.uelmax($i) -text " Max"
    entry $w.f.uemax($i) -textvariable vars($w,statsmax,$i) -width 7
    pack $w.f.sel$i $w.f.uelmin($i) $w.f.uemin($i) $w.f.uelmax($i) $w.f.uemax($i) -in $w.f.statsentries($i) -padx 3 -pady 3 -side left
    incr i 1
}

render:updatestats $w 

$w.f.uelmin(1) configure -fg "red"
$w.f.uelmax(1) configure -fg "yellow"
$w.f.uelmin(2) configure -fg "blue"
$w.f.uelmax(2) configure -fg "cyan"

#}}}
    #{{{ Background image

set vars($w,statsmin,0) "min"
set vars($w,statsmax,0) "max"

frame $w.f.statsentries(0)

FileEntry $w.f.sel0 -textvariable fentries($w,0) -label " Background image " -title "Select" -width 35 -filetypes IMAGE
$w.f.sel0.labf configure -width 15
#FSLFileEntry $w.f.sel0 -variable fentries($w,0) -pattern "IMAGE" -directory $PWD -label " Background image " -labelwidth 17 -title "Select" -width 35 -filterhist VARS(history)
    

label $w.f.uelmin(0) -text " Min"
entry $w.f.uemin(0) -textvariable vars($w,statsmin,0) -width 7
label $w.f.uelmax(0) -text " Max"
entry $w.f.uemax(0) -textvariable vars($w,statsmax,0) -width 7
pack $w.f.sel0 $w.f.uelmin(0) $w.f.uemin(0) $w.f.uelmax(0) $w.f.uemax(0) -in $w.f.statsentries(0) -padx 3 -pady 3 -side left

$w.f.uelmax(0) configure -fg "white"

#}}}
    #{{{ Output image

    frame $w.f.statsentries(3)

FileEntry $w.f.sel3 -textvariable fentries($w,3) -label " Output image " -title "Select" -width 35 -filetypes IMAGE
$w.f.sel3.labf configure -width 15
#FSLFileEntry $w.f.sel3 -variable fentries($w,3) -pattern "IMAGE" -directory $PWD -label " Output image " -labelwidth 17 -title "Select" -width 35 -filterhist VARS(history)

    pack $w.f.sel3 -in $w.f.statsentries(3) -padx 3 -pady 3 -side left


#}}}
    #{{{ Colourmap type

frame $w.f.types

optionMenu2 $w.f.type vars($w,type) 0 "Solid colours" 1 "Transparent colours"
#tixOptionMenu $w.f.type -label " Colourmap type:  " -variable vars($w,type)
#$w.f.type add command 0 -label "Solid colours"
#$w.f.type add command 1 -label "Transparent colours"
set vars($w,type) 1

set vars($w,checker) 0
label $w.f.checkerlabel -text " Apply checkerboard mask to colour overlay"
checkbutton $w.f.checker -variable vars($w,checker)

pack $w.f.type $w.f.checkerlabel $w.f.checker -in $w.f.types -padx 3 -pady 3 -side left

#}}}
    #{{{ Image type:

set vars($w,output) 0

optionMenu2 $w.f.output vars($w,output) 0 "Floating point" 1 "Integer (required for 3D rendering in MEDx)"
#tixOptionMenu $w.f.output -label "  Output data type:  " -variable vars($w,output)
#$w.f.output add command 0 -label "Floating point"
#$w.f.output add command 1 -label "Integer (required for 3D rendering in MEDx)"

pack $w.f.output $w.f.types -anchor w -side bottom
    pack $w.f.statsentries(3) -anchor w -side bottom

pack $w.f.statsentries(0) -anchor w -side bottom

#}}}
    #{{{ button frame

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    pack $w.btns.b -side bottom -fill x -padx 3 -pady 3
    
    button $w.apply -command "render:apply $w" -text "Go" -width 5
    bind $w.apply <Return> {
	[winfo toplevel %W].apply invoke
    }
    pack $w.apply -in $w.btns.b -side left -expand yes -padx 3 -pady 3 -fill y
	    
    set vars($w,cmap) 0


    button $w.cancel    -command "render:destroy $w" \
	    -text "Exit" -width 5
    bind $w.cancel <Return> {
	[winfo toplevel %W].cancel invoke
    }

    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/redirects/miscvis.html" \
            -text "Help" -width 5
    bind $w.help <Return> {
        [winfo toplevel %W].help invoke
    }
    pack $w.cancel $w.help -in $w.btns.b -side left -expand yes -padx 3 -pady 3 -fill y

#}}}
		
    pack $w.f $w.btns -expand yes -fill both
}

#}}}
#{{{ render:apply

proc render:apply { w } {
    global vars entries fentries FSLDIR


    if {  $vars($w,cmap) } {
	#{{{ setup display range for current image, if CALmin/max are set

MxGetCurrentImage OUT

MxGetHeaderInfo $OUT 0 OP

if { [ keylget OP AVWDsr.Value.IDCalMin.Value dmin ] && [ keylget OP AVWDsr.Value.IDCalMax.Value dmax ] } {
    puts "Display range $dmin $dmax"
    MxSetDisplayRange $OUT V $dmin $dmax
}

#}}}
    } else {
	#{{{ start setting up thecommand

set thecommand "$FSLDIR/bin/overlay $vars($w,type) $vars($w,output)"

if { $vars($w,checker) } {
    set thecommand "$thecommand -c"
}

#}}}
	#{{{ process entries

set i 0
while { $i <= $vars($w,statss) } {

    #{{{ if INMEDX save pages to file



#}}}
    #{{{ check dimensions are consistent

    set xdim($i) [ exec sh -c "$FSLDIR/bin/fslval $fentries($w,$i) dim1" ]
    set ydim($i) [ exec sh -c "$FSLDIR/bin/fslval $fentries($w,$i) dim2" ]
    set zdim($i) [ exec sh -c "$FSLDIR/bin/fslval $fentries($w,$i) dim3" ]

    if { $i > 0 } {
	if { $zdim($i) != $zdim(0) || $ydim($i) != $ydim(0) || $xdim($i) != $xdim(0) } {
	    MxPause "Not all of the stats images are of the same size as the input image"
	    return
	}
    }

#}}}
    #{{{ setup thresholds

    if { $vars($w,statsmin,$i) == "" || $vars($w,statsmax,$i) == "" } {
	MxPause "Not all of the thresholds have been set"
	return
    }

    set RANGE -R
    if { $i == 0 } {
	set RANGE -r
    }

    catch { exec sh -c "${FSLDIR}/bin/fslstats $fentries($w,$i) $RANGE" } minmax

    if { $vars($w,statsmin,$i) == "min" } {
	set vars($w,statsmin,$i) [ lindex $minmax 0 ]
    }
    if { $vars($w,statsmax,$i) == "max" } {
	set vars($w,statsmax,$i) [ lindex $minmax 1 ]
    }

#}}}

    set thecommand "$thecommand [ remove_ext $fentries($w,$i) ] $vars($w,statsmin,$i) $vars($w,statsmax,$i)"

    incr i 1
}

#}}}
	#{{{ run the program

    if { $fentries($w,3) == "" || ! [ file writable [ file dirname $fentries($w,3) ] ] } {
	MxPause "Please select writable output image"
	return
    }

set thecommand "$thecommand [ remove_ext $fentries($w,3) ]"

puts $thecommand
set result [ catch { exec sh -c $thecommand } ErrMsg ]
if {$result != 0} {
    puts "$ErrMsg"
}

#}}}
}

    #{{{ load LUTS



#}}}
    
    set vars($w,cmap) 0
    puts Done
    update idletasks
}

#}}}
#{{{ render:destroy

proc render:destroy { w } {
    destroy $w
}

#}}}
#{{{ render:updatestats

proc render:updatestats { w } {
    global vars

    set i 1
    while { $i <= $vars($w,maxstatss) } {
	pack forget $w.f.statsentries($i)
	incr i 1
    }

    set i 1
    while { $i <= $vars($w,statss) } {
	pack $w.f.statsentries($i) -in $w.f -anchor w
	incr i 1
    }
}

#}}}
#{{{ tail end

wm withdraw .
render .rename
tkwait window .rename

#}}}
