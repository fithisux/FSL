set FSLDIR        $env(FSLDIR)
set FSLOUTPUTTYPE $env(FSLOUTPUTTYPE)
set FSLTCLSH      $env(FSLTCLSH)
set FSLWISH       $env(FSLWISH)

set USER       $env(USER)
set HOME       [ exec sh -c " cd ; pwd " ]
set PWD        [ exec sh -c " pwd " ]
set PROCID     [ pid ]
set HOSTNAME   [ exec hostname ]
set OS         [ exec uname ]

set auto_path [ linsert $auto_path 0 $FSLDIR/tcl/ ]
set MYEXEC    [ info nameofexecutable ]
set MYSHELL   [ file tail $MYEXEC ]

# test for whether we are in tclsh or wish
if { [  string match -nocase *wish* $MYSHELL ] } {
    package require Tk
#    tk_focusFollowsMouse
#bind Button <Enter> { focus %W ; tk::ButtonEnter %W }
    bind Button <1> { focus %W ; tk::ButtonDown %W }
    option add *Font {Helvetica 10} widgetDefault
    option add *NoteBook.font {Helvetica 10 bold}
    option add *Foreground black
    option add *background lightgray 
    option add *highlightBackground lightgray
    option add *FileEntry.Background grey95
    option add *FileEntry.e.background grey95
    option add *LabelComboBox.combo.e.background grey95
    option add *Scrollbar.Background grey60
    option add *highlightColor grey86
    option add *Checkbutton*selectColor #ffff00  
    option add *Radiobutton*selectColor #ffff00  
    option add *Listbox.background grey95
    option add *activeforeground #b0b0ff
    option add *activeBackground grey86
    option add *insertBackground #ffb0b0 
    option add *troughColor grey95
    option add *selectBackground #b0b0ff
    option add *disabledBackground #ffb0b0
    option add *disabledForeground grey77
    option add *Entry.Background grey95 
    option add *Entry.width 4
    #New options for Spinbox...
    option add *Entry.validate focusout 
    option add *SpinBox.vcmd2  {validNum %W %V %P %s [lindex [[string range %W 0 end-2] cget -range] 0] [lindex [[string range %W 0 end-2] cget -range] 1]} 
    option add *SpinBox.invcmd2 { set [%W cget -textvariable] $tempSpin; %W config -validate %v }
    lappend auto_path $FSLDIR/tcl/BWidget
    package require BWidget
    package require BWFileEntry
    package require BWLabelComboBox 
    package require BWLabelSpinBox 
}

#Choose load_varian specific OR FSLPARALLEL
set FSLPARALLEL 0
if { [ info exists env(SGE_ROOT) ] && $env(SGE_ROOT) != "" } { set FSLPARALLEL 1 }

# what type of OS?
set OSFLAVOUR unix
set gui_ext ""
set FSLSLASH ""
set UNAME [ exec uname ]
if { $UNAME == "Darwin" } {
    set OSFLAVOUR macos
    set gui_ext "_gui"
} elseif { [ string compare CYGWIN [ string range $UNAME 0 5 ] ] == 0 } {
    set OSFLAVOUR cygwin
    set gui_ext "_gui"
    set FSLSLASH "C:/cygwin"
}

