
proc FmribWebHelp { prefix file } {

    global OSFLAVOUR

    regsub -all "//" $file "/" file

    if { $OSFLAVOUR == "macos" } {

        catch { exec sh -c "open $file" & }

    } elseif { $OSFLAVOUR == "cygwin" } {

	set url [ exec sh -c "cygpath -w $file" ]
	eval exec [auto_execok start] {"${prefix}//$url"}

    } else {

	foreach executable {firefox mozilla netscape iexplorer opera konquerer galeon amaya mosaic lynx w3m links browsex elinks} {
	    set executable [auto_execok $executable]
	    if [string length $executable] {
		catch { exec sh -c "$executable -remote \"openURL(${prefix}//${file},new-window)\" || $executable ${prefix}//${file} " 2> /dev/null & }
		break
	    }
	}
    }
}

