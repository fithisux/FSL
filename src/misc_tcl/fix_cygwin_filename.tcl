
proc fix_cygwin_filename { f } {

    global OSFLAVOUR FSLSLASH

    if { $OSFLAVOUR == "cygwin" } {

	regsub -all {\\} $f / f

	set s [ string last $FSLSLASH $f ]
	if { $s > -1 } {
	    set s [ expr $s + [ string length $FSLSLASH ] - 1 ]
	    set f [ string replace $f 0 $s ]
	}

	# the following should work instead of the above, but it doesn't....
	# set f [ exec sh -c "cygpath -u $f" ]

    }

    return $f
}

