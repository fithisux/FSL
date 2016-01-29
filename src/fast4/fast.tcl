#

#   fast.tcl - GUI for FAST - FMRIB's Automated Segmentation Tool
#
#   Stephen Smith, Matthew Webster and Yongyue Zhang, FMRIB Image Analysis Group
#
#   Copyright (C) 2001-2006 University of Oxford
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

#{{{ setups

source [ file dirname [ info script ] ]/fslstart.tcl

set VARS(history) {}

set MC 5

catch { exec sh -c "${FSLDIR}/bin/fast 2>&1 | grep Version | awk '{ print \$4 }' -" } VERSION

#}}}
#{{{ fast

proc fast { w } {

    #{{{ vars and setup

global vars FSLDIR argc argv PWD HOME VERSION entries MC

toplevel $w
wm title $w "FAST - FMRIB's Automated Segmentation Tool - v$VERSION"
wm iconname $w "FAST"
wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

frame $w.f

set entries($w,$MC) ""

#}}}
    TitleFrame $w.f.input -text "Input" -relief groove 
    set lfinput [ $w.f.input getframe ]
    #{{{ input channels

set vars(channels) 1

LabelSpinBox $lfinput.channels -label "Number of input channels " -textvariable vars(channels) -range " 1 $MC 1 " -command " $lfinput.channels.spin.e validate; fast:updateinputs $w" -modifycmd "fast:updateinputs $w"

pack $lfinput.channels -in $lfinput -side top -anchor w -padx 5 -pady 5

#}}}
    #{{{ input image (if not in medx)

set i 0
while { $i < $MC } {

    if { $i == 0 && $argc > 0 && [ string length [ lindex $argv 0 ] ] > 0 } {
	set inputname [ imglob -extension [ lindex $argv 0 ] ]
	if { [ string first / $inputname ] == 0 || [ string first ~ $inputname ] == 0 } {
	    set entries($w,$i) $inputname
	} else {
	    set entries($w,$i) ${PWD}/$inputname
	}
	set entries($w,$MC) $entries($w,$i)
    }

    FileEntry $lfinput.input$i -textvariable entries($w,$i) -label "Input image   " -title "Select the input image"  -width 40 -filedialog directory  -filetypes IMAGE -command "fast:select $w" 

    incr i 1
}

#}}}
    #{{{ input options

frame $lfinput.inopts

label $lfinput.inopts.typelabel -text "Image type"

set vars(type) 1
optionMenu2  $lfinput.inopts.type vars(type) 1 "T1-weighted" 2 "T2-weighted" 3 "Proton density"

pack $lfinput.inopts.typelabel $lfinput.inopts.type -in $lfinput.inopts -side left

#}}}
    TitleFrame $w.f.output -text "Output" -relief groove 
    set lfoutput [ $w.f.output getframe ]
    #{{{ output image (if not in medx)

FileEntry $lfoutput.output -textvariable entries($w,$MC) -label "Output image(s) basename" -title "Select a basename for the output image(s)" -width 40 -filedialog directory  -filetypes * 

    pack $lfoutput.output -in $lfoutput -side top -padx 5 -pady 5 -anchor w

#}}}
    #{{{ output classes

set vars(classes) 3


LabelSpinBox $lfoutput.classes -label "Number of classes" -textvariable vars(classes) -range " 2 6 1 " -command "$lfoutput.classes.spin.e validate;fast:updateinputs $w" -modifycmd  "fast:updateinputs $w"

pack $lfoutput.classes -in $lfoutput -side top -anchor w -padx 5 -pady 5

#}}}
    #{{{ output options

frame $lfoutput.outoptsA

label $lfoutput.outoptsA.label -text "Output images:"

frame $lfoutput.outoptsB

label $lfoutput.outoptsB.seglabel -text "Binary segmentation: Also output one image per class"
set vars(segall) 0
checkbutton $lfoutput.outoptsB.segall -variable vars(segall)

frame $lfoutput.outoptsC

set vars(prob) 0
# label $lfoutput.outoptsC.problabel -text "Probability maps"
# checkbutton $lfoutput.outoptsC.prob -variable vars(prob)

label $lfoutput.outoptsC.pvlabel -text "Partial volume maps"
set vars(pv) 1
checkbutton $lfoutput.outoptsC.pv -variable vars(pv)

label $lfoutput.outoptsC.restoredlabel -text "Restored input"
set vars(restored) 0
checkbutton $lfoutput.outoptsC.restored -variable vars(restored)

label $lfoutput.outoptsC.biaslabel -text "Estimated Bias field"
set vars(bias) 0
checkbutton $lfoutput.outoptsC.bias -variable vars(bias)

pack $lfoutput.outoptsA.label -side left
pack $lfoutput.outoptsB.seglabel $lfoutput.outoptsB.segall -side left
pack 	$lfoutput.outoptsC.pvlabel $lfoutput.outoptsC.pv \
	$lfoutput.outoptsC.restoredlabel $lfoutput.outoptsC.restored \
	$lfoutput.outoptsC.biaslabel $lfoutput.outoptsC.bias \
	-in $lfoutput.outoptsC -side left

pack $lfoutput.outoptsA $lfoutput.outoptsB $lfoutput.outoptsC -in $lfoutput -side top -anchor w -padx 5

#}}}

    #{{{ Advanced Options

collapsible frame $w.f.opts -title "Advanced options"

TitleFrame $w.f.opts.b.f -text "Advanced" -relief groove 
set advf [ $w.f.opts.b.f getframe ]

#{{{ use a-priori probability map

frame $advf.apriori

set vars(apriori) 0
checkbutton $advf.apriori.yn -variable vars(apriori) -command "fast:updateapriori $w"
label $advf.apriori.label -text "Use a-priori probability maps for initialisation"

set vars(aprioriTransform) "${FSLDIR}/etc/flirtsch/ident.mat"
FileEntry $advf.apriori.transform -textvariable vars(aprioriTransform) -label "Standard to Input FLIRT transform" -title "Select a flirt transformation matrix" -width 30 -filedialog directory  -filetypes *.mat

set vars(apriori_final) 0
checkbutton $advf.apriori.yn_final -variable vars(apriori_final)
label $advf.apriori.label_final -text "and for final segmentation"

pack $advf.apriori.transform -side bottom
pack $advf.apriori.label $advf.apriori.yn -side left

#}}}
#{{{ initial manual segment

set vars(manualSegment) ""

FileEntry $advf.manualSegment -textvariable vars(manualSegment) -label "Use file of initial tissue-type means" -title "Select a file of initial tissue-type means" -width 30 -filedialog directory  -filetypes * 


#

#}}}

#{{{ MRF parameter

set vars(mrfbeta) 0.1

LabelSpinBox $advf.mrfbeta -label "Main MRF parameter " -textvariable vars(mrfbeta) -range " 0.0 1000.0 0.1 " -command "$advf.mrfbeta.spin.e validate;fast:updateinputs $w" -modifycmd  "fast:updateinputs $w"

pack $advf.mrfbeta -in $advf -side top -anchor w -padx 5 -pady 5

#}}}

#{{{ bias field iterations

set vars(biasiter) 4

LabelSpinBox $advf.biasiter -label "Number of iterations for bias field removal " -textvariable vars(biasiter) -range " 1 100 1 " -command "$advf.biasiter.spin.e validate;fast:updateinputs $w" -modifycmd  "fast:updateinputs $w"

pack $advf.biasiter -in $advf -side top -anchor w -padx 5 -pady 5

#}}}

#{{{ bias field smoothing

set vars(biassmooth) 20.0

LabelSpinBox $advf.biassmooth -label "Bias field smoothing (FWHM in mm) " -textvariable vars(biassmooth) -range " 0.1 1000.0 0.1 " -command "$advf.biassmooth.spin.e validate;fast:updateinputs $w" -modifycmd  "fast:updateinputs $w"

pack $advf.biassmooth -in $advf -side top -anchor w -padx 5 -pady 5

#}}}


pack $advf.apriori -in $advf -side top -anchor w -padx 5 -pady 5

pack $w.f.opts.b.f -in $w.f.opts.b -anchor w

#}}}

    fast:updateinputs $w
    pack $w.f.input $w.f.output $w.f.opts -in $w.f -side top -padx 5 -pady 0 -anchor w

    #{{{ Button Frame

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
 
    button $w.apply     -command "fast:apply $w" \
        -text "Go" -width 5
 
    button $w.cancel    -command "destroy $w" \
        -text "Exit" -width 5
 
    button $w.help -command "FmribWebHelp file: ${FSLDIR}/doc/redirects/fast.html" \
            -text "Help" -width 5

    pack $w.btns.b -side bottom -fill x
    pack $w.apply $w.cancel $w.help -in $w.btns.b -side left -expand yes -padx 3 -pady 10 -fill y
 
    pack $w.f $w.btns -expand yes -fill both -padx 5 -pady 5

#}}}
}

#}}}
#{{{ fast:updateinputs

proc fast:updateinputs { w } {

    global vars MC

    set lfinput [ $w.f.input getframe ]
    set advf [ $w.f.opts.b.f getframe ]

    set i 0
    while { $i < $MC } {
	pack forget $lfinput.input$i
	incr i 1
    }

    if { $vars(channels) == 1 } {

	pack $lfinput.inopts -in $lfinput -side top -padx 5 -pady 5 -anchor w

	$lfinput.input0 configure -label "Input image"
	pack $lfinput.input0 -in $lfinput -side top -padx 5 -pady 5 -anchor w -after $lfinput.channels

	pack $advf.manualSegment -in $advf -side top -anchor w -padx 5 -pady 5 -after $advf.apriori 

    } else {

	pack forget $lfinput.inopts

	set i [ expr $vars(channels) - 1 ]
	while { $i >= 0 } {
	    $lfinput.input$i configure -label "Input image [ expr $i + 1 ]"
	    pack $lfinput.input$i -in $lfinput -side top -padx 5 -pady 5 -anchor w -after $lfinput.channels
	    incr i -1
	}

	pack forget $advf.manualSegment
    }
}

#}}}
#{{{ fast:updateapriori

proc fast:updateapriori { w } {

    global vars

    set vars(apriori_final) 0

    set advf [ $w.f.opts.b.f getframe ]

    pack forget $advf.apriori.yn_final $advf.apriori.label_final

    if { $vars(apriori) == 1 } {
	pack $advf.apriori.label_final $advf.apriori.yn_final -in $advf.apriori -side left
    }
}

#}}}
#{{{ fast:select

proc fast:select { w { dummy "" } } {

    global vars entries MC

    set entries($w,0) [ remove_ext $entries($w,0) ]

    set entries($w,$MC) $entries($w,0)

#    if { [ string length $entries($w,$MC) ] == 0 } {
#	set entries($w,$MC) [ file rootname $entries($w,0) ]
#    }

}

#}}}
#{{{ fast:apply

proc fast:apply { w } {

    global vars entries MC

    set inlist ""
    for { set i 0 } { $i < $vars(channels) } { incr i 1 } {
	lappend inlist $entries($w,$i)
    }

    fast:proc $vars(channels) $inlist $vars(type) $entries($w,$MC) $vars(classes) $vars(segall) $vars(prob) $vars(pv) $vars(restored) $vars(bias) $vars(apriori) $vars(apriori_final) $vars(manualSegment) $vars(biasiter) $vars(biassmooth) $vars(mrfbeta)

    update idletasks
}

#}}}
#{{{ fast:proc

proc fast:proc { channels inlist type output classes segall prob pv restored bias apriori apriori_final manualSegment biasiter biassmooth mrfbeta } {

    #{{{ setup for running fast 

global FSLDIR HOME vars

set output [ file rootname $output ]

#}}}
    #{{{ run command

if { $channels == 1 } {
    set thecommand "$FSLDIR/bin/fast -t $type"
} else {
    set thecommand "$FSLDIR/bin/fast -S $channels"
}

set thecommand "$thecommand -n $classes"
set thecommand "${thecommand} -H $mrfbeta"
set thecommand "${thecommand} -I $biasiter"
set thecommand "${thecommand} -l $biassmooth"


if { $segall } {
    set thecommand "${thecommand} -g"
}

if { $prob } {
    set thecommand "${thecommand} -p" 
}

if { ! $pv } {
    set thecommand "${thecommand} --nopve"
}

if { $restored } {
    set thecommand "${thecommand} -B" 
}

if { $bias } {
    set thecommand "${thecommand} -b"
}

if { $apriori_final } {
    set thecommand "${thecommand} -a $vars(aprioriTransform) -P"
} elseif { $apriori } {
    set thecommand "${thecommand} -a $vars(aprioriTransform)"
}

if { $manualSegment != "" } {
    set thecommand "${thecommand} -s $manualSegment" 
}

set thecommand "$thecommand -o $output $inlist"

fsl:exec "$thecommand" -t 20

puts "Finished"

#}}}
}

#}}}
#{{{ call GUI

wm withdraw .
fast .rename
tkwait window .rename

#}}}

