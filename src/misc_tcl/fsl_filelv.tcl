proc feat_file:setup_dialogFid { w output_file output_filhis output_selhis name_space pattern title command dirasfile} {
    global data env 
    upvar 1 $output_file outputfile
    #old version of FSLFile_new for standalone
    upvar 1 $output_filhis filterhist
    upvar 1 $output_selhis selhis  

     if { [ info exists env(w_filesel) ] } {
	if { [ winfo exists $env(w_filesel) ] } {
            raise $env(w_filesel)
	    return 0
	}
    }
    set count 0
    set w0 ".wdialog[incr count]"
    while { [ winfo exists $w0 ] } {
        set w0 ".wdialog[incr count]"
    }
    set env(w_filesel) $w0
    toplevel $w0 
    wm title $w0 $title
    wm geometry $w0 600x300
    set data(-dirasfile) {}
    if { [ info exists dirasfile ] }  { set data(-dirasfile) $dirasfile }
    #FilterFrame
    frame $w0.f1 -border 10
    LabelComboBox   $w0.f1.filter -command "feat_file:filter $w0" -modifycmd "feat_file:filter $w0" -label "Filter:" -width 50
    set data(w:filter)  $w0.f1.filter
    pack $data(w:filter) -side top -expand yes -fill both
    #CenterFrame
    PanedWindow $w0.f2 -side top
    #Left
    set dir [$w0.f2 add -weight 10 ]
    $dir config -relief flat 
    label $dir.lab -text "Directories:" -underline 0
    set data(w:dirlist) [ listbox $dir.list -selectmode single -yscrollcommand "$dir.sb set" -width 4 -height 10 ]
    scrollbar $dir.sb -command "$data(w:dirlist) yview "
    pack $dir.lab -side top -fill x -padx 10
    pack $dir.sb -side right -fill y -padx {0 10}
    pack  $data(w:dirlist) -side right -expand yes -fill both -padx {10 0}
    # Right
    set file [$w0.f2 add -weight 7]
    $file config -relief flat
    if { $data(-dirasfile) != {} } { label $file.lab -text "Filtered Directories:" -underline 2 } else {
    label $file.lab -text "Files:" -underline 2 }
    set data(w:filelist) [ listbox $file.list -selectmode single -yscrollcommand "$file.sb set" -width 4 -height 10 ]
    scrollbar $file.sb -command "$data(w:filelist) yview "
    pack $file.lab -side top 
    pack $file.sb -side right -fill y -padx {0 10}
    pack  $data(w:filelist) -side right -expand yes -fill both -padx {10 0}
    # Far Right
    set info [$w0.f2 add -weight 13]
    $info config -relief flat
    label $info.lab -text "Info:" -underline 2 
    set data(w:infolist)  [ listbox $info.list -selectmode single -yscrollcommand "$info.sb set" -width 4 -height 10 ]
    scrollbar $info.sb -command "$data(w:infolist) yview "
    pack $info.lab -side top -fill x -padx 10
    pack $info.sb -side right -fill y -padx {0 10}
    pack $data(w:infolist) -side right -expand yes -fill both -padx {10 0}
    #SelectionFrame
    frame $w0.f3 -border 10
    LabelComboBox   $w0.f3.selection -command "feat_file:invoke $w0 $output_file $output_filhis $output_selhis $name_space {$command}"  -modifycmd "feat_file:invoke $w0 $output_file $output_filhis  $output_selhis $name_space {$command}; return -code break" -label "Selection:" 
    set data(w:selection) $w0.f3.selection
    pack $data(w:selection) -side top -fill both
    #Buttons
    frame $w0.f4
    button $w0.f4.but_ok -text "Ok" -command "feat_file:invoke $w0 $output_file $output_filhis $output_selhis $name_space {$command}"
    button $w0.f4.but_filt -text "Filter" -command "feat_file:filter $w0" 
    button $w0.f4.but_cancel -text "Cancel" -command "destroy $w0"
    pack $w0.f4.but_ok  $w0.f4.but_filt  $w0.f4.but_cancel -side left -padx 10
    #ConstructFrame
    pack $w0.f1 -in $w0 -side top -fill x
    pack  $w0.f4 -side bottom -pady 10
    pack $w0.f3 -in $w0 -side bottom -fill x
    pack $w0.f2 -in $w0 -side top -fill both -expand yes
    #Initialise Variables
    set data(flag)      0
    set data(fakeDir)   0
    set data(-browsecmd) {}
    set data(-selhist) {}
    set data(-filterhist) {}
   
    #if { [ info exists $output_filhis ] } { set data(-filterhist) $filterhist }
    if { [ info exists filterhist ] } { set data(-filterhist) $filterhist }
    if { [ info exists selhis ] } { set data(-selhist) $selhis }
    set data(-pattern) $pattern

    if { [info exists outputfile] } { 
	if { $outputfile != [ $data(w:selection).combo.e get ] } {
           $data(w:selection) configure -text [rmSlash $outputfile]
	    if { [file isdirectory $outputfile] && [file exists $outputfile] } { 
		if { [ lindex $data(-filterhist) end ] != "$outputfile/$data(-pattern)" } { lappend data(-filterhist) $outputfile/$data(-pattern) }            
            } else {
		if {[file exists [file dir $outputfile]] && [file isdirectory [file dir $outputfile]] && [ lindex $data(-filterhist) end ] != "[file dir $outputfile]/$data(-pattern)"  } {  
               lappend data(-filterhist) [file dir $outputfile]/$data(-pattern) 
            }}
            if { [ lindex $data(-selhist) end ] != "$outputfile" } { lappend data(-selhist) $outputfile }
	}
    }

    if { $data(-selhist) != {} } { 
       $data(w:selection) configure -values  $data(-selhist) 
	$data(w:selection).combo setvalue last
    }

    #If there is a previous filter history obtain values from pattern etc from it 
    if { $data(-filterhist) != {} } { 
       $data(w:filter) configure -values  $data(-filterhist) 
       $data(w:filter).combo setvalue last
       feat_file:InterpFilter $w0 [ $data(w:filter).combo.e get ]
    } else {
       if {[info exists env(PWD)]} { set data(-directory) $env(PWD) } else { set data(-directory) [pwd] }
       set filter [feat_file:GetFilter $w0 $data(-directory) $data(-pattern)]
       $data(w:filter) configure -text "$filter"
       $data(w:filter) configure -values [feat_file:InterpFilter $w $filter]
    }
    feat_file:LoadDirIntoLists w0

    if { [info exists outputfile] && ![file isdirectory $outputfile] && [file exists $outputfile]} { 
                   set i 0
		while {[file tail $outputfile ] != [ $file.list get $i ] && $i < 1000 } {
                    $file.list selection clear $i
		    set i [ expr $i + 1 ]
                    $file.list selection set $i
               }     
            }
    bind  $data(w:dirlist) <Double-1>      "feat_file:InvokeDirFid $w0"
    bind  $data(w:dirlist) <<ListboxSelect>> "feat_file:SelectDirFid $w0"
    bind  $data(w:filelist) <<ListboxSelect>> "feat_file:SelectFileFid $w0" 
    bind  $data(w:filelist) <Double-1> "feat_file:invoke $w0 $output_file $output_filhis $output_selhis $name_space {$command}"
    bind  $w0 <Destroy> "set env(w_filesel) -1"
}

proc feat_file:SelectDirFid {w} {
    global data
    if {$data(fakeDir) > 0} {
	incr data(fakeDir) -1
	$data(w:dirlist) select clear 0 end
	$data(w:dirlist) activate -1
	return
    }

    if {$data(flag)} {
	return
    }
    set data(flag) 1
    set subdir [$data(w:dirlist) curselection]
    if {$subdir != ""} {
    set subdir [ $data(w:dirlist) get  [ $data(w:dirlist) curselection ] ]
    }
    if {$subdir == {}} {
	set subdir "."
    }

    set filter [feat_file:GetFilter $w $data(-directory) $subdir/$data(-pattern)]
    $data(w:filter) configure -text $filter    

    if {[file exists $data(-directory)/$subdir/info.patient]==1} {
	$data(w:infolist) select clear 0 end
  
        set infoout [open "$data(-directory)/$subdir/info.patient"]
	gets $infoout help
	gets $infoout help
	gets $infoout help
	while {[gets $infoout help] >= 0 } {
	    $data(w:infolist) insert end $help
	}
	close $infoout
    } else {
	    $data(w:infolist) select clear 0 end
    }
    set data(flag) 0
}


proc feat_file:InvokeDirFid {w} {
    global data

    feat_file:SelectDirFid $w

    set theDir [$data(w:dirlist) get active]

    set OWD [ pwd ]
    cd $data(-directory)/$theDir
    set data(-directory) [ pwd ]
    cd $OWD

    $data(w:dirlist) select clear 0 end

    feat_file:InterpFilter $w [feat_file:GetFilter $w $data(-directory) $data(-pattern)]

    feat_file:LoadDir $w
    #update selection if go back via ..
    $data(w:selection) configure -text "$data(-directory)"
    
}

proc feat_file:SelectFileFid {w} {
    global data

    if {$data(flag)} {
	return
    }
    set data(flag) 1

    # Reset the "Filter:" box to the current directory:
    $data(w:dirlist) select clear 0 end
    set filter [feat_file:GetFilter $w $data(-directory) $data(-pattern)]
    $data(w:filter) configure -text $filter
    # Now select the file

    set selected [$data(w:filelist) curselection]
    if {$selected != ""} {
    set selected [ $data(w:filelist) get  [ $data(w:filelist) curselection ] ]
    }

    if {$selected  != {}} {
	# Make sure that the selection is not empty!
	if {$data(-directory) == "/"} {
            $data(w:selection) configure -text "/$selected"
	    set data(-value) /$selected
	} else {
            $data(w:selection) configure -text "$data(-directory)/$selected"
	    set data(-value) $data(-directory)/$selected
	}
	if {$data(-browsecmd) != {}} {
	    #tixEvalCmdBinding $w $data(-browsecmd) {} [$data(w:selection) cget -value]
	}
   }
    set data(flag) 0
    $data(w:infolist) select clear 0 end
    foreach {infofile} {comment sequence options seq_params scan_range matrix misc fmri_protocol} {
	
	if {[file exists $data(-directory)/$selected/info.$infofile]==1} {
	    set infoout [open "$data(-directory)/$selected/info.$infofile"]
	    gets $infoout help
	    gets $infoout help
	    gets $infoout help
	    while {[gets $infoout help] >= 0 } {
		$data(w:infolist) insert end $help
	    }
	    close $infoout
	}
    }
}
