#{{{ copyright and setup

#   make_flobs - GUI for setting up FLOBS basis set
#
#   Stephen Smith, Mark Woolrich and Matthew Webster, FMRIB Image Analysis Group
#
#   Copyright (C) 2006 University of Oxford
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
#{{{ make_flobs

proc make_flobs { w } {
    global PXHOME FSLDIR VARS argc argv PWD entries flobsvars

    #{{{ setup main window

toplevel $w -visual truecolor

wm title      $w "Make FLOBS"
wm iconname   $w "Make FLOBS"
wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

frame $w.f
pack $w.f -in $w -side top

set graphpic [ image create photo -file ${FSLDIR}/tcl/halfcoshrfparams.ppm ]
button $w.f.pic -image $graphpic -borderwidth 1

set flobsvars(m1min) 0
set flobsvars(m1max) 2
set flobsvars(m2min) 3
set flobsvars(m2max) 8
set flobsvars(m3min) 3
set flobsvars(m3max) 8
set flobsvars(m4min) 3
set flobsvars(m4max) 8
set flobsvars(cmin) 0
set flobsvars(cmax) 0.3

frame $w.f.g
LabelSpinBox  $w.f.g.m1min -textvariable flobsvars(m1min) -label "m1   Min" -range " 0.0 100000 0.5 "
grid $w.f.g.m1min -in $w.f.g -column 0 -row 0
LabelSpinBox  $w.f.g.m1max -textvariable flobsvars(m1max) -label "   Max"   -range " 0.0 100000 0.5 "
grid $w.f.g.m1max -in $w.f.g -column 1 -row 0
LabelSpinBox  $w.f.g.m2min -textvariable flobsvars(m2min) -label "m2   Min" -range " 0.0 100000 0.5 "
grid $w.f.g.m2min -in $w.f.g -column 0 -row 1
LabelSpinBox  $w.f.g.m2max -textvariable flobsvars(m2max) -label "   Max"   -range " 0.0 100000 0.5 "
grid $w.f.g.m2max -in $w.f.g -column 1 -row 1
LabelSpinBox  $w.f.g.m3min -textvariable flobsvars(m3min) -label "m3   Min" -range " 0.0 100000 0.5 "
grid $w.f.g.m3min -in $w.f.g -column 0 -row 2
LabelSpinBox  $w.f.g.m3max -textvariable flobsvars(m3max) -label "   Max"   -range " 0.0 100000 0.5 "
grid $w.f.g.m3max -in $w.f.g -column 1 -row 2
LabelSpinBox  $w.f.g.m4min -textvariable flobsvars(m4min) -label "m4   Min" -range " 0.0 100000 0.5 "
grid $w.f.g.m4min -in $w.f.g -column 0 -row 3
LabelSpinBox  $w.f.g.m4max -textvariable flobsvars(m4max) -label "   Max"   -range " 0.0 100000 0.5 "
grid $w.f.g.m4max -in $w.f.g -column 1 -row 3
LabelSpinBox  $w.f.g.cmin -textvariable flobsvars(cmin) -label "c    Min" -range " 0.0 100000 0.1 "
grid $w.f.g.cmin -in $w.f.g -column 0 -row 4
LabelSpinBox  $w.f.g.cmax -textvariable flobsvars(cmax) -label "   Max"   -range " 0.0 100000 0.1 "
grid $w.f.g.cmax -in $w.f.g -column 1 -row 4

set flobsvars(nflobs) 3
LabelSpinBox  $w.f.nflobs -textvariable flobsvars(nflobs) -label "Number of basis functions" -range " 1 100000 1 "

FileEntry $w.f.featdir -textvariable entries($w,1) -label "Output directory" -title "Output directory" -width 20 -filedialog directory  -filetypes *.flobs 

pack $w.f.pic $w.f.g $w.f.nflobs $w.f.featdir -in $w.f -padx 5 -pady 5 -anchor w

#}}}
    #{{{ Button Frame

frame $w.f.btns
    
button $w.f.btns.preview -command "make_flobs:apply $w 0" -text "Preview"

button $w.f.btns.go -command "make_flobs:apply $w 1" -text "Go"

button $w.f.btns.cancel -command "destroy $w" -text "Exit"

button $w.f.btns.help -command "FmribWebHelp file: ${FSLDIR}/doc/redirects/flobs.html" -text "Help" -width 5
 
pack $w.f.btns.preview $w.f.btns.go $w.f.btns.cancel $w.f.btns.help -in $w.f.btns -side left -expand yes -padx 3 -pady 3 -fill y

pack $w.f.btns -in $w.f -side bottom -fill x -padx 3 -pady 3

#}}}
}

#}}}
#{{{ make_flobs:apply

proc make_flobs:apply { w dooutput } {
    global FSLDIR FSLSLASH entries flobsvars

    if { $dooutput && $entries($w,1) == "" } {
	MxPause "Please select output directory"
	return 1
    }

    set tmpdir [ exec sh -c "${FSLDIR}/bin/tmpnam /tmp/flobs" ]

    #{{{ write params file

    set channel [ open ${tmpdir}.txt "w" ]

    puts $channel "$flobsvars(m1min) $flobsvars(m1max)"
    puts $channel "$flobsvars(m2min) $flobsvars(m2max)"
    puts $channel "$flobsvars(m3min) $flobsvars(m3max)"
    puts $channel "$flobsvars(m4min) $flobsvars(m4max)"
    puts $channel "0 0"

    puts $channel "0 0"
    puts $channel "$flobsvars(cmin)  $flobsvars(cmax)"
    puts $channel "0 0"

    close $channel

#}}}

    set window_length [ expr int( $flobsvars(m1max) + $flobsvars(m2max) + $flobsvars(m3max) + $flobsvars(m4max) + 2 )  ]

    # set number of samples
    set nhs 1000
    if { ! $dooutput } {
	set nhs 200
    }

    catch { exec sh -c "${FSLDIR}/bin/halfcosbasis --hf=${tmpdir}.txt --nbfs=$flobsvars(nflobs) --ns=$window_length --logdir=${tmpdir}.flobs --nhs=$nhs --res=0.05" } errmsg
    set outdir ${tmpdir}.flobs

    catch { exec sh -c "/bin/mv ${tmpdir}.txt ${tmpdir}.flobs/params.txt" } errmsg

    if { $dooutput } {
	set outdir [ new_filename [ file rootname $entries($w,1) ].flobs ]
	catch { exec sh -c "/bin/mv ${tmpdir}.flobs $outdir" } errmsg
    }

    #{{{ write web page report

set channel [ open $outdir/report.html "w" ]

puts $channel "<HTML>

<TITLE>FLOBS Report</TITLE>

<BODY BACKGROUND=\"file:${FSLSLASH}${FSLDIR}/doc/images/fsl-bg.jpg\">

<hr><CENTER>
<H1>Basis Set Report</H1>
[ exec date ]
"

if { $dooutput } {
    puts $channel "<br>FLOBS output is
<br>$outdir
<br>To use this in FEAT, select the <b>Optimal/custom basis functions</b> convolution option in the <b>Stats</b> section of the FEAT GUI. Using the file selector which then appears, choose this FLOBS output directory. Make this HRF setting and file selection for each relevant <b>original EV</b>.
"
}

puts $channel "<hr>

m1=\[$flobsvars(m1min)-$flobsvars(m1max)s\] &nbsp;&nbsp; m2=\[$flobsvars(m2min)-$flobsvars(m2max)s\] &nbsp;&nbsp; m3=\[$flobsvars(m3min)-$flobsvars(m3max)s\] &nbsp;&nbsp; m4=\[$flobsvars(m4min)-$flobsvars(m4max)s\] &nbsp;&nbsp; c=\[$flobsvars(cmin)-$flobsvars(cmax)\]

<p><IMG BORDER=0 SRC=\"hrfsamps.png\">
<p><IMG BORDER=0 SRC=\"eigenvalues.png\"> &nbsp;&nbsp; <IMG BORDER=0 SRC=\"hrfbasisfns.png\">

</CENTER>

<p>
Basis set generated by FLOBS (FMRIB's Linear Optimal Basis Set)<br>
M.W. Woolrich, T.E.J. Behrens, and S.M. Smith. Constrained linear basis sets for HRF modelling using Variational Bayes. NeuroImage, 21:4(1748-1761) 2004.

<HR><FONT SIZE=1>This page produced automatically by Make_flobs - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL</A>.</FONT>

</BODY></HTML>
"

close $channel

#}}}

    FmribWebHelp file: $outdir/report.html
}

#}}}
#{{{ call GUI and wait

wm withdraw .
make_flobs .r
tkwait window .r

#}}}
