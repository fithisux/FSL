#{{{ copyright and setup

#   featquery - apply masking etc to get out stats from FEAT runs
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 2002-2008 University of Oxford
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
#{{{ featquery_select_label

proc featquery_select_label { } {
    global FSLDIR fmri atlasname atlaslabelcount atlaslabelid atlaslabelname atlasimage

    set count 0
    set w0 ".dialog[incr count]"
    while { [ winfo exists $w0 ] } {
        set w0 ".dialog[incr count]"
    }

    toplevel $w0

    wm iconname $w0 $atlasname($fmri(atlasmask))
    wm title $w0 $atlasname($fmri(atlasmask))

    frame $w0.f
    pack $w0.f -expand yes -fill both -in $w0 -side top

    canvas $w0.f.viewport -yscrollcommand "$w0.f.ysbar set"
    scrollbar $w0.f.ysbar -command "$w0.f.viewport yview" -orient vertical
    frame $w0.f.viewport.f
    $w0.f.viewport create window 0 0 -anchor nw -window $w0.f.viewport.f
    bind $w0.f.viewport.f <Configure> "feat5:scrollform_resize $w0 $w0.f.viewport"
    pack $w0.f.viewport -side left -fill both -expand true -in $w0.f

    # test for whether this is a label image, in which case omit the first ("0") label from the button list
    set theatlasimage "[ string range $atlasimage($fmri(atlasmask)) 0 [ expr [ string last "mm" $atlasimage($fmri(atlasmask)) ] - 2 ] ]${fmri(regres)}mm"
    set nvols [ fsl:exec "${FSLDIR}/bin/fslnvols ${FSLDIR}/data/atlases/$theatlasimage" -n ]
    set starti 0
    if { $nvols == 1 } {
	set starti 1
    }

    for { set i $starti } { $i < $atlaslabelcount($fmri(atlasmask)) } { incr i 1 } {
	button $w0.button$i -text "$atlaslabelid($fmri(atlasmask),$i) $atlaslabelname($fmri(atlasmask),$i)" -command "set fmri(atlaslabel) $i ; featquery_update ; destroy $w0" -anchor w
	pack $w0.button$i -in $w0.f.viewport.f -side top -expand yes -fill both -padx 0 -pady 0
    }
}

#}}}
#{{{ featquery GUI

proc featquery { w } {
    global fmri PXHOME FSLDIR VARS argc argv PWD feat_files query vars n_atlases atlasname atlasimage atlaslabelcount atlaslabelid atlaslabelname atlasmenu

    toplevel $w
    wm title      $w "FEATQuery"
    wm iconname   $w "FEATQuery"
    wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

    set fmri(multiple) 1
    set fmri(level) 1
    set fmri(analysis) 4
    set fmri(inputtype) 1
    set fmri(regres) 2

    #{{{ select FEAT directories

TitleFrame $w.select -text "Input FEAT directories" -relief groove 
set fmri(selectf) [ $w.select getframe ]

LabelSpinBox $fmri(selectf).number -label "Number of FEAT directories " -textvariable fmri(multiple) -range {1 10000 1 } -width 3

button $fmri(selectf).button -text "Select" -command "feat5:multiple_select $w -1 \"Select FEAT directories\" "

pack $fmri(selectf).number $fmri(selectf).button -in $fmri(selectf) -padx 5 -pady 5 -side left

#}}}
    #{{{ setup mask or coordinates

TitleFrame $w.roi -text "Input ROI selection" -relief groove 
set fmri(roif) [ $w.roi getframe ]

#{{{ atlas

frame $fmri(roif).atlas

set fmri(atlasmask) 0
set fmri(oldatlasmask) 0
set fmri(atlaslabel) -1

label $fmri(roif).atlas.label -text "Use atlas:"

set atlasmenu ""
parseatlases

eval "optionMenu2 $fmri(roif).atlas.menu fmri(atlasmask) -command \"featquery_update\" 0 None $atlasmenu"

pack $fmri(roif).atlas.label $fmri(roif).atlas.menu -in $fmri(roif).atlas -padx 5 -side left

frame $fmri(roif).atlas2
label $fmri(roif).atlas2.labelselected -text ""
button $fmri(roif).atlas2.button -command "featquery_select_label" -text "Select label"
pack $fmri(roif).atlas2.button $fmri(roif).atlas2.labelselected -in $fmri(roif).atlas2 -padx 5 -side left

#}}}
#{{{ mask

frame $fmri(roif).mask
set fmri(maskf) $fmri(roif).mask

FileEntry $fmri(maskf).mask -textvariable fmri(mask) -label "Mask image " -title "Select the mask image" -width 30 -filedialog directory -filetypes IMAGE

frame $fmri(maskf).type
set fmri(masktype) 1
radiobutton $fmri(maskf).type.image -text "Use mask image"   -value 1 -variable fmri(masktype) -command "featquery_update"
radiobutton $fmri(maskf).type.coord -text "Use co-ordinates" -value 2 -variable fmri(masktype) -command "featquery_update"
pack $fmri(maskf).type.image $fmri(maskf).type.coord -in $fmri(maskf).type -padx 5 -side left -expand yes

pack $fmri(maskf).mask $fmri(maskf).type -in $fmri(maskf) -padx 5 -pady 5 -side top

frame $fmri(maskf).coords

set fmri(coordtype) -vox
radiobutton $fmri(maskf).coords.vox -text "vox" -value -vox -variable fmri(coordtype)
radiobutton $fmri(maskf).coords.mm  -text "mm"  -value -mm  -variable fmri(coordtype)

set fmri(cX) 0
set fmri(cY) 0
set fmri(cZ) 0

LabelSpinBox $fmri(maskf).coords.cX -label "  X " -textvariable fmri(cX) -range {-10000 10000 1 } -width 5
LabelSpinBox $fmri(maskf).coords.cY -label "Y " -textvariable fmri(cY) -range {-10000 10000 1 } -width 5
LabelSpinBox $fmri(maskf).coords.cZ -label "Z " -textvariable fmri(cZ) -range {-10000 10000 1 } -width 5

pack $fmri(maskf).coords.vox $fmri(maskf).coords.mm $fmri(maskf).coords.cX $fmri(maskf).coords.cY $fmri(maskf).coords.cZ -in $fmri(maskf).coords -padx 5 -side left -expand yes

#}}}

pack $fmri(roif).atlas $fmri(maskf) -in $fmri(roif)

#}}}
    #{{{ options

TitleFrame $w.opts -text "Output options" -relief groove 
set fmri(optsf) [ $w.opts getframe ]

#{{{ atlas

frame $fmri(optsf).atlas

set fmri(atlas) 0

label $fmri(optsf).atlas.label -text "Use atlas:"

eval "optionMenu2 $fmri(optsf).atlas.menu fmri(atlas) 0 None $atlasmenu"

pack $fmri(optsf).atlas.label $fmri(optsf).atlas.menu -in $fmri(optsf).atlas -padx 5 -side left

#}}}
#{{{ percent

frame $fmri(optsf).percent

set fmri(showpercent) 0
checkbutton $fmri(optsf).percent.yn -variable fmri(showpercent) -text "Convert PE/COPE values to % (VARCOPE to %^2)"

pack $fmri(optsf).percent.yn -in $fmri(optsf).percent -padx 5 -side left

#}}}
#{{{ mask weighting

frame $fmri(optsf).maskweight

set fmri(maskweight_yn) 0
checkbutton $fmri(optsf).maskweight.yn -variable fmri(maskweight_yn) -text "Do not binarise mask (allow weighting)" -command featquery_update

pack $fmri(optsf).maskweight.yn -in $fmri(optsf).maskweight -padx 5 -side left

#}}}
#{{{ interp thresholding

frame $fmri(optsf).ithresh

set fmri(ithresh_yn) 0
checkbutton $fmri(optsf).ithresh.yn -variable fmri(ithresh_yn) -text "Change post-interpolation thresholding of mask" -command "featquery_update"

set fmri(ithresh) 0.5

LabelSpinBox  $fmri(optsf).ithresh.thresh -textvariable fmri(ithresh)    -range {0.0 10000 0.1 } -width 5

pack $fmri(optsf).ithresh.yn $fmri(optsf).ithresh.thresh  -in $fmri(optsf).ithresh -padx 5 -side left

#}}}
#{{{ thresholding

frame $fmri(optsf).thresh

set fmri(statsthresh_yn) 0
checkbutton $fmri(optsf).thresh.yn -variable fmri(statsthresh_yn) -text "Threshold stats images as well as masking" -command "featquery_update"

set fmri(statsthresh) 0
LabelSpinBox  $fmri(optsf).thresh.thresh -textvariable fmri(statsthresh) -range {-10000.0 10000 1 } -width 5

pack $fmri(optsf).thresh.yn $fmri(optsf).thresh.thresh  -in $fmri(optsf).thresh -padx 5 -side left

#}}}
#{{{ tsplot

frame $fmri(optsf).tsplot

set fmri(tsplot_yn) 1
checkbutton $fmri(optsf).tsplot.yn -variable fmri(tsplot_yn) -text "Create timeseries plots"

pack $fmri(optsf).tsplot.yn -in $fmri(optsf).tsplot -padx 5 -side left

#}}}
#{{{ web browser popup

frame $fmri(optsf).popup

set fmri(fqpopup) 1
checkbutton $fmri(optsf).popup.yn -variable fmri(fqpopup) -text "Popup results in web browser"

pack $fmri(optsf).popup.yn -in $fmri(optsf).popup -padx 5 -side left

#}}}
#{{{ output directory name

frame $fmri(optsf).name
label $fmri(optsf).name.label -text "Featquery output directory name"
set fmri(output) "featquery"
entry $fmri(optsf).name.entry -textvariable fmri(output) -width 15
pack $fmri(optsf).name.label $fmri(optsf).name.entry -in $fmri(optsf).name -padx 5 -side left

#}}}

pack $fmri(optsf).atlas $fmri(optsf).percent $fmri(optsf).maskweight $fmri(optsf).ithresh $fmri(optsf).thresh $fmri(optsf).tsplot $fmri(optsf).popup $fmri(optsf).name -in $fmri(optsf) -side top -anchor w

#}}}
    #{{{ bottom buttons

frame $w.btns -borderwidth 1

button $w.btns.apply  -text "Go"   -command "featquery_proc $w"
button $w.btns.cancel -text "Exit" -command "destroy $w"
button $w.btns.help   -text "Help" -command "FmribWebHelp file: ${FSLDIR}/doc/redirects/featquery.html"

pack $w.btns.apply $w.btns.cancel $w.btns.help -in $w.btns -padx 5 -pady 5 -side left -expand yes

#}}}

    pack $w.select $w.roi $w.opts $w.btns -in $w -padx 5 -pady 5 -fill x 
    featquery_update
}

#}}}
#{{{ featquery_update

proc featquery_update { } {
    global fmri atlaslabelname

    if { $fmri(atlasmask) != $fmri(oldatlasmask) } {
	set fmri(atlaslabel) -1
    }

    pack forget $fmri(maskf) $fmri(maskf).coords $fmri(roif).atlas2
    if { $fmri(atlasmask) == 0 } {
	pack $fmri(maskf) -in $fmri(roif) -after $fmri(roif).atlas 
	if { $fmri(masktype) == 2 } {
	    pack $fmri(maskf).coords -in $fmri(maskf) -padx 5 -pady 5 -side top
	}
    } else {
	set fmri(masktype) 1
	set fmri(maskweight_yn) 1
	if { $fmri(atlaslabel) >= 0 } {
	    $fmri(roif).atlas2.labelselected configure -text "$atlaslabelname($fmri(atlasmask),$fmri(atlaslabel))"
	} else {
	    $fmri(roif).atlas2.labelselected configure -text ""
	}
	pack $fmri(roif).atlas2 -in $fmri(roif) -side top
    }

    pack forget $fmri(optsf).ithresh
    if { $fmri(maskweight_yn) == 0 } {
	pack $fmri(optsf).ithresh -in $fmri(optsf) -after $fmri(optsf).maskweight -side top -anchor w
    }

    pack forget $fmri(optsf).thresh.thresh
    if { $fmri(statsthresh_yn) } {
	pack $fmri(optsf).thresh.thresh -in $fmri(optsf).thresh -after $fmri(optsf).thresh.yn -padx 5 -side left
    }

    pack forget $fmri(optsf).ithresh.thresh
    if { $fmri(ithresh_yn) } {
	pack $fmri(optsf).ithresh.thresh -in $fmri(optsf).ithresh -after $fmri(optsf).ithresh.yn -padx 5 -side left
    }

    set fmri(oldatlasmask) $fmri(atlasmask)
}

#}}}
#{{{ featquery_whichstats

proc featquery_whichstats { w } {
    global fmri feat_files PWD FSLDIR atlasimage

    if { [ winfo exists $w.stats ] } {
	destroy $w.stats
    }

    TitleFrame $w.stats -text "Stats images of interest" -relief groove 
    set fmri(statsf) [ $w.stats getframe ]
    set w0 $fmri(statsf)
    pack $w.stats -in $w -after $w.select -padx 5 -pady 5 -fill x

    #{{{ setup scrollbar viewport as $v

    frame $w0.f
    pack $w0.f -in $w0 -side top -anchor w
    canvas $w0.f.viewport -yscrollcommand "$w0.f.ysbar set" -xscrollcommand "$w0.xsbar set" -borderwidth 0
    scrollbar $w0.xsbar -command "$w0.f.viewport xview" -orient horizontal
    scrollbar $w0.f.ysbar -command "$w0.f.viewport yview" -orient vertical
    frame $w0.f.viewport.f
    $w0.f.viewport create window 0 0 -anchor nw -window $w0.f.viewport.f
    bind $w0.f.viewport.f <Configure> "feat5:scrollform_resize $w0 $w0.f.viewport"
    pack $w0.f.viewport -side left -fill both -expand true -in $w0.f

    set v $w0.f.viewport.f

#}}}

    cd $feat_files(1)

    set row 0

    set fmri(nstats) 0
    set fmri(statslist) ""

    foreach statstype { stats/pe stats/cope stats/varcope stats/tstat stats/fstat stats/zstat stats/zfstat thresh_zstat thresh_zfstat } {
    
	set statslist [ lsort -dictionary [ imglob ${statstype}*.* ] ] 
    
	if { [ llength $statslist ] > 0 } {

	    label $v.${row}_0 -text "$statstype  "
	    grid $v.${row}_0 -in $v -column 0 -row $row
	    
	    set column 1

	    foreach stats $statslist {
		set fmri(dostats.$fmri(nstats)) 0
		set i [ string trimleft $stats "abcdefghijklmnopqrstuvwxyz_/" ]
		checkbutton $v.${row}_$column -variable fmri(dostats.$fmri(nstats)) -text "$i "
		grid $v.${row}_$column -in $v -column $column -row $row
		incr column 1

		set fmri(statslist) "$fmri(statslist) $stats"
		incr fmri(nstats) 1
	    }
	    
	    incr row 1
	}
	
    }

    if { [ imtest reg/standard ] } {
	set fmri(regres) [ expr round ( abs ( [ exec sh -c "$FSLDIR/bin/fslval reg/standard pixdim1" ] ) ) ]
    }

    cd $PWD
}

#}}}
#{{{ featquery_proc

proc featquery_proc { w } {
    global fmri feat_files FSLDIR n_atlases atlasname atlasimage atlaslabelcount atlaslabelid atlaslabelname atlasmenu

    set mn ""
    if { $fmri(atlasmask) > 0 } {
	set atlasimage($fmri(atlasmask)) "[ string range $atlasimage($fmri(atlasmask)) 0 [ expr [ string last "mm" $atlasimage($fmri(atlasmask)) ] - 2 ] ]${fmri(regres)}mm"
	set fmri(mask) [ fsl:exec "${FSLDIR}/bin/tmpnam /tmp/featquery" -n ]
	set nvols [ fsl:exec "${FSLDIR}/bin/fslnvols ${FSLDIR}/data/atlases/$atlasimage($fmri(atlasmask))" -n ]
	if { $nvols == 1 } {
	    fsl:exec "${FSLDIR}/bin/fslmaths ${FSLDIR}/data/atlases/$atlasimage($fmri(atlasmask)) -thr $fmri(atlaslabel) -uthr $fmri(atlaslabel) -bin $fmri(mask)"
	} else {
	    fsl:exec "${FSLDIR}/bin/fslroi ${FSLDIR}/data/atlases/$atlasimage($fmri(atlasmask)) $fmri(mask) $fmri(atlaslabel) 1"
	}
	fsl:echo ${fmri(mask)}.name "$atlaslabelname($fmri(atlasmask),$fmri(atlaslabel))"
    }

    if { ( $fmri(mask) == "" ) ||
	 ( ! [ imtest $fmri(mask) ] && ! [ imtest $feat_files(1)/$fmri(mask) ] ) } {
	puts "\n\nYou need to set the mask image!\n\n"
	return
    }

    set thecommand "${FSLDIR}/bin/featquery $fmri(multiple)"

    for { set n 1 } { $n <= $fmri(multiple) } { incr n 1 } {
	set thecommand "$thecommand $feat_files($n)"
    }

    set nstats 0
    set statslist ""
    for { set n 0 } { $n < $fmri(nstats) } { incr n 1 } {
	if { $fmri(dostats.$n) == 1 } {
	    set statslist "$statslist [ lindex $fmri(statslist) $n ]"
	    incr nstats 1
	}
    }
    set thecommand "$thecommand $nstats $statslist $fmri(output)"

    if { $fmri(atlas) > 0 } {
    	set thecommand "$thecommand -a $fmri(atlas)"
    }

    if { $fmri(showpercent) == 1 } {
    	set thecommand "$thecommand -p"
    }

    if { $fmri(statsthresh_yn) } {
	set thecommand "$thecommand -t $fmri(statsthresh)"
    }

    if { $fmri(ithresh_yn) } {
	set thecommand "$thecommand -i $fmri(ithresh)"
    }

    if { $fmri(tsplot_yn) } {
	set thecommand "$thecommand -s"
    }

    if { $fmri(maskweight_yn) } {
	set thecommand "$thecommand -w"
    }

    if { $fmri(fqpopup) } {
	set thecommand "$thecommand -b"
    }

    set thecommand "$thecommand $fmri(mask)"

    if { $fmri(masktype) == 2 } {
	set thecommand "$thecommand $fmri(coordtype) $fmri(cX) $fmri(cY) $fmri(cZ)"
    }

    puts $thecommand
    catch { exec sh -c "$thecommand" & } junk
}

#}}}
#{{{ call GUI and wait

wm withdraw .
featquery .r
tkwait window .r

#}}}

