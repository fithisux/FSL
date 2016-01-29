#usage optionMenu2 .m someVar [-command "command" ] h1 "Hello1" h2 "Hello2" h3 "Hello3"
#NB This includes an optional [ -command variable ] which MUST be placed after path/variable
#and before arguement list. This is NOT a robust implementation - don't even consider
#trying to use configure -command (although you can use set $w.menu.command at at pinch)
#This uses a fairly "clunky" method to remove the trace (by unsetting and resetting variable)
#due to a more orthodox vdelete not working correctly. This requires the pertinent variables
 proc optionMenu2 {w varName firstValue firstText args} {
    # Initialize the variable
    upvar #0 $varName var
    if { ![info exists var] } {
            set var $firstValue
    }

    # Create an array of the displayed text for each value
    set map($firstValue) $firstText
    array set map $args

    # Create the menubutton and associated menu
    menubutton $w \
        -text $map($var) \
        -indicatoron 1 \
        -menu $w.menu \
        -relief raised \
        -bd 2 \
        -highlightthickness 2 \
        -anchor c  \
        -direction below
    menu $w.menu -tearoff 0

    #If the next arguement is -command, set up command call else setup first menu item
    if { ![string compare $firstValue "-command"] } { set command $firstText } else {
	
    set command "" 
    $w.menu add radiobutton \
        -label      $firstText \
        -value      $firstValue \
        -variable   $varName \
        -indicatoron false
    }
    foreach {value text} $args {
        $w.menu add radiobutton \
              -label      $text \
            -value      $value \
            -indicatoron false \
            -variable   $varName \
            -command $command
    }
    # This procedure will keep the menubutton's text sync'ed with the variable's content.
    #$var($varIndex) added to allow for array element passing...
    proc optionMenuTrace {w valueText varName varIndex varOp} {
        upvar #0 $varName var
        array set map $valueText
        $w configure -text $map($var($varIndex))   
    }

    proc optionMenuRemoveTrace {varName} {
    global fmri
        set temp [set [set varName]]
        unset $varName
        set $varName $temp
    }
    # Define a trace so the menubutton's text is sync'ed when the variable changes.
    set valueText [array get map]
    uplevel #0 [list trace add variable $varName write [list optionMenuTrace $w $valueText]]
    # Remove the trace when the menubutton is destroyed.
    #bind $w <Destroy> [list trace vdelete $varName w [list optionMenuTrace $w $valueText]]
    bind $w <Destroy> "optionMenuRemoveTrace $varName"

    # Return the menu's path
    return $w.menu
} 


# usage -vcmd {validNum %W %V %P %s +10 +92}
#NB This uses the global variable tempSpin
    proc validNum {win event X oldX min max } {
        global tempSpin
        # Make sure min<=max
        if {$min > $max} {
            set tmp $min; set min $max; set max $tmp
        }
        set default $min
	if { [ expr $min + $max ] == 0 } { set default 0 }
        if {[string is int -strict $min] && [string is int -strict $max] } {set strongCheck [expr {[string is int -strict $X]}] } else {set strongCheck [expr {[string is double -strict $X]}] }
        if { !$strongCheck } { set tempSpin $default }
        if { $strongCheck && $X < $min } { set tempSpin $min }
        if { $strongCheck && $X > $max } { set tempSpin $max }
	if {[string is int -strict $min] && [string is int -strict $max] } { set strongCheck [expr {[string is int -strict $X] && ($X >= $min) && ($X <= $max)}]} else { set strongCheck [expr {[string is double -strict $X] && ($X >= $min) && ($X <= $max)}]}  
        return $strongCheck
    }
