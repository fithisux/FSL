#
# InputFromWidget - Generic Input Widget
#	Will get any number of inputs from the user.
#
#	Input:
#		Title - title for window
#		NumInputs - number of inputs
#		Text1 - text for the first input field
#		Value1 - value entered for the first input field
#		Width1 - width of the first input fields
#		...
#		Textn - text for the nth input field
#		Valuen - value entered for the nth input field
#		Widthn - width of the nth input fields
#
#	Side Effect:	Each of the Value arguments may be changed
#			by the user's inputs
#
proc InputFromWidget {Title NumInputs args} {

	#
	# Make an independent window widget and title it
	#
	toplevel .main -borderwidth 5
	wm title .main $Title
	wm iconname .main $Title

	#
	# Make the Each input widget
	#
	for {set i 0} {$i < [expr $NumInputs*3]} {incr i 3} {
		frame .main.frame$i
		label .main.frame$i.text -text [lindex $args $i]
		entry .main.frame$i.input -width [lindex $args \
		    [expr $i + 2]] -relief sunken -textvariable \
		    [lindex $args [expr $i + 1]]
	}

	#
	# Make an OK button.  Use the command to destroy this
	#	widget and effectively quit this procedure
	#
	button .main.ok -text OK -command {destroy .main} -width 5

	#
	# Pack them together.  This also displays them on the screen.
	#
	for {set i 0} {$i < [expr $NumInputs*3]} {incr i 3} {
		pack .main.frame$i -side top
		pack .main.frame$i.text .main.frame$i.input \
			-side left -padx 1m -pady 2m
	}
	pack .main.ok -side top

	#
	# Wait until the main window is destroyed.  This
	#	happens when the user clicks on the OK button.
	#
	tkwait window .main
	update idletasks
}

#
# Update the user about what's happening...
#
#	Side Effect:	Creates/Overwrites global variable UpdateMessage
#
proc ScriptUpdate {Message} {
	global UpdateMessage

	set UpdateMessage $Message

	#
	# Make an independent window if necessary
	#
	if {[set exists [winfo exists .update]] == 0} {
		toplevel .update -borderwidth 5
		wm title .update "Script Update"
		wm iconname .update "Script Update"

		#
		# Add the message
		#
		set font [FindFont "-adobe-helvetica-bold-r" 18]
		message .update.msg -width 8c -justify center \
			-font $font -textvariable UpdateMessage

		pack .update.msg
	}

	update idletasks
}

#
# Get rid of the Update window...
#
proc CancelScriptUpdate {} {

	if [winfo exists .update] {
		destroy .update
	}

	update idletasks
}


#
# Yes/No Type of Widget
#       Will get a yes no type of response from the user
#
#       Input:
#               Question - Question to pose to the user
#               yes - text for the "yes" button
#               no - text for the "no" button
#
#	Side Effect:	Creates/Overwrites global array variable YesNoWidget
#
proc YesNoWidget {Question yes no} {
	global YesNoWidget

	#
	# Make an independent window widget and title it
	#
	toplevel .main -borderwidth 5
	wm title .main "Script Q/A"
	wm iconname .main "Script Q/A"

	#
	# If the user has positioned us previously - let's
	#	put ourselves there
	#
	if {[lsearch [array names YesNoWidget] position] != -1} {
	    # get out just the position part of the geometry...
	    set plus [string first "+" $YesNoWidget(position)]
	    set minus [string first "-" $YesNoWidget(position)]
	    if {$plus == -1 || ($minus < $plus && $minus != -1)} {
		set start $minus
	    } else { set start $plus }
	    set position [string range $YesNoWidget(position) $start \
		[string length $YesNoWidget(position)]]
	    wm geometry .main $position
	}

	#
	# Display the Question
	#
	frame .main.f
	set font [FindFont "-adobe-helvetica-bold-r" 18]
	message .main.msg -text $Question -width 300 \
		-font $font -justify center
	pack .main.msg -in .main.f -side left
	pack .main.f -side top -fill x

	#
	# Make the buttons
	#
	frame .main.buttons
	button .main.buttons.yes -text $yes -command \
		{set YesNoWidget(reply) 1}
	button .main.buttons.no -text $no -command \
		{set YesNoWidget(reply) 0}
	pack .main.buttons -side top
	pack .main.buttons.yes .main.buttons.no -side left \
		-padx 1m -pady 2m
 
	#
	# Wait until something happens to the YesNoWidget(reply) variable
	#
	tkwait variable YesNoWidget(reply)

	#
	# remember our current position
	#
	set YesNoWidget(position) [wm geometry .main]

	#
	# destroy the window and return the value of the reply
	#
	destroy .main
	update idletasks

	return $YesNoWidget(reply)
}


#
# Define MxPause
#
proc MxPause { string } {
    global prompt MxPause

    set w [toplevel .prompt -borderwidth 5]
    wm title $w " "
    wm iconname $w " "

    #
    # If the user has positioned us previously - let's
    #       put ourselves there
    #
    if {[lsearch [array names MxPause] position] != -1} {
        # get out just the position part of the geometry...
        set position [string range $MxPause(position) \
            [string first "+" $MxPause(position)] \
            [string length $MxPause(position)]]
        wm geometry $w $position
    }

    frame $w.f
    #set font [FindFont "-adobe-helvetica-bold-r" 18]
    #message $w.msg -text $string -width 300 -font $font -justify center
    message $w.msg -text $string -width 500
    set b [frame $w.buttons -bd 1]

    pack $w.msg -in $w.f -side left
    pack $w.f $w.buttons -side top -fill x

    button $b.ok -text OK -command "MxPause:done $w" -width 5
    pack $b.ok -side top
    focus $b.ok
 
    tkwait visibility $w
    tkwait window $w

    update idletasks

    return 0
}

proc MxPause:done {w} {
    global MxPause

    set MxPause(position) [wm geometry $w]
    destroy $w
}

proc FindFont { pattern size } {

    set any [catch {exec xlsfonts -fn "$pattern*" | grep -v " "} fonts]
#    set any [catch {exec xlsfonts -fn "$pattern*"} fonts]

    if {$any == 1} {
        return fixed
    }

    # Find the font with the size that is closest
    # to the desired size
    set diff 999
    set rfont null
    foreach f $fonts {
	set s [ lindex [ split $f - ] 7 ]
	set d [ expr abs ( $size - $s ) ]
	if {$d < $diff} { set diff $d; set rfont $f }
    }

    if {$rfont != "null"} { return $rfont } else { return fixed }
}
