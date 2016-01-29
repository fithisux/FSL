# ----------------------------------------------------------------------
#  EXAMPLE: use "wm" commands to manage a balloon help window
# ----------------------------------------------------------------------
#  Effective Tcl/Tk Programming
#    Mark Harrison, AsiaInfo Software Research and Development
#    Michael McLennan, Bell Labs Innovations for Lucent Technologies
#    Addison-Wesley Professional Computing Series
# ======================================================================
#  Copyright (c) 1996-1997  Lucent Technologies Inc. and Mark Harrison
# ======================================================================

option add *Balloonhelp*background white widgetDefault
option add *Balloonhelp*foreground black widgetDefault
option add *Balloonhelp.info.wrapLength 500 widgetDefault
option add *Balloonhelp.info.justify left widgetDefault
#option add *Balloonhelp.info.font -*-lucida-medium-r-normal-sans-*-120-* widgetDefault

toplevel .balloonhelp -class Balloonhelp -background black -borderwidth 1 -relief flat

label .balloonhelp.arrow -anchor nw -bitmap @[file join ${FSLDIR}/tcl/arrow.xbm] -bg #b0ffb0 -fg black
pack .balloonhelp.arrow -side left -fill y

label .balloonhelp.info -bg #b0ffb0 -fg black
pack .balloonhelp.info -side left -fill y

wm overrideredirect .balloonhelp 1
wm withdraw .balloonhelp
# ----------------------------------------------------------------------
#  USAGE:  balloonhelp_for <win> <mesg>
#
#  Adds balloon help to the widget named <win>.  Whenever the mouse
#  pointer enters this widget and rests within it for a short delay,
#  a balloon help window will automatically appear showing the
#  help <mesg>.
# ----------------------------------------------------------------------
proc balloonhelp_for {win mesg} {
    global bhInfo
    set bhInfo($win) $mesg
    set bhInfo(click) 0

       bind $win <B1-Leave> {set bhInfo(click) 1}
    bind $win <Enter> {if {$bhInfo(click)!=1} { balloonhelp_pending %W } else {balloonhelp_cancel;set bhInfo(click) 0 } }
    bind $win <Leave> {balloonhelp_cancel; set bhInfo(click) 0}
 
}

# ----------------------------------------------------------------------
#  USAGE:  balloonhelp_control <state>
#
#  Turns balloon help on/off for the entire application.
# ----------------------------------------------------------------------
# modified for auto-toggle
set bhInfo(active) 1

proc balloonhelp_control {state} {
    global bhInfo
    if {$state==2 && $bhInfo(active)} {set state 0}
    if {$state==2 && !$bhInfo(active)} {set state 1}
    if {$state==1} {set bhInfo(active) 1} 
    if {$state==0} { 
        balloonhelp_cancel
        set bhInfo(active) 0
    }

}

# ----------------------------------------------------------------------
#  USAGE:  balloonhelp_pending <win>
#
#  Used internally to mark the point in time when a widget is first
#  touched.  Sets up an "after" event so that balloon help will appear
#  if the mouse remains within the current window.
# ----------------------------------------------------------------------
proc balloonhelp_pending {win} {
    global bhInfo

    balloonhelp_cancel
    set bhInfo(pending) [after 3000 [list balloonhelp_show $win]]
}

# ----------------------------------------------------------------------
#  USAGE:  balloonhelp_cancel
#
#  Used internally to mark the point in time when the mouse pointer
#  leaves a widget.  Cancels any pending balloon help requested earlier
#  and hides the balloon help window, in case it is visible.
# ----------------------------------------------------------------------
proc balloonhelp_cancel {} {
    global bhInfo

    if {[info exists bhInfo(pending)]} {
        after cancel $bhInfo(pending)
        unset bhInfo(pending)
    }
    wm withdraw .balloonhelp
}

# ----------------------------------------------------------------------
#  USAGE:  balloonhelp_show <win>
#
#  Used internally to display the balloon help window for the
#  specified <win>.
# ----------------------------------------------------------------------
proc balloonhelp_show {win} {
    global bhInfo
  
    if { $bhInfo(active) && [ winfo exists $win ] } {
    .balloonhelp.info configure -text $bhInfo($win)
#Define initial position + screensize
    set x [expr [winfo rootx $win]+[winfo width $win]]
    set y [expr [winfo rooty $win]]
    set scrwidth  [winfo vrootwidth  $win]
    set scrheight [winfo vrootheight $win]
#Create Initial window on bottom pixel
    set x1 [expr {$scrwidth - 1}]
    set y1 [expr {$scrheight - 1}]
    wm geometry .balloonhelp +$x1+$y1
    wm deiconify .balloonhelp
    raise .balloonhelp
    update
    set width  [winfo reqwidth .balloonhelp]
    set height [winfo reqheight .balloonhelp]
    set scrwidth  [winfo vrootwidth  $win]
    set scrheight [winfo vrootheight $win]
    # Reset the balloon if it is off the screen...

    if {[expr {$x + $width}] > $scrwidth} {set x [expr {$x-$width-[winfo width $win]}] }
    if {$x < 0} {set x 0}
    if {[expr {$y + $height}] > $scrheight} {set y [expr {$scrheight - $height}]}
    if {$y < 0} {set y 0}
    # Redisplay the balloon
        wm geometry .balloonhelp +$x+$y
        wm deiconify .balloonhelp
        after idle raise .balloonhelp
        after idle update
    }
    if { [ info exists bhInfo(pending) ] } { unset bhInfo(pending) }
}

