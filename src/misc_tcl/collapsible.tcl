 #########################################################################
 #                                                                       #
 # Copyright (C) 1993 by General Electric company.  All rights reserved. #
 #                                                                       #
 # Permission to use, copy, modify, and distribute this                  #
 # software and its documentation for any purpose and without            #
 # fee is hereby granted, provided that the above copyright              #
 # notice appear in all copies and that both that copyright              #
 # notice and this permission notice appear in supporting                #
 # documentation, and that the name of General Electric not be used in   #
 # advertising or publicity pertaining to distribution of the            #
 # software without specific, written prior permission.                  #
 #                                                                       #
 # General Electric makes no representations about the suitability of    #
 # this software for any purpose.  It is provided ``as is''              #
 # without express or implied warranty.                                  #
 #                                                                       #
 # This work was supported in part by the DARPA Initiative in Concurrent #
 # Engineering (DICE) through DARPA Contracts MDA972-88-C-0047 and       #
 # MDA972-92-C-0027.                                                     #
 #                                                                       #
 # This work was supported in part by the Tri-Services Microwave and     #
 # Millimeter-Wave Advanced Computational Environment (MMACE) program    #
 # under Naval Research Laboratory contract N00014-92-C-2044.            #
 #                                                                       #
 #########################################################################


# File: collapsible.tcl
#
# Description:
#	Procedure to define a `collapsible' object.
#
# Global variables:
#c	collapsible_priv(visible,$w)
#		1 if the collapsible widget $w is visible, and 0 otherwise.

 # $Id: collapsible.tcl,v 1.1 2003/02/21 16:33:00 mark Exp $
 # $Source: /usr/local/share/sources/misc_tcl/collapsible.tcl,v $
 # $Log: collapsible.tcl,v $
 # Revision 1.1  2003/02/21 16:33:00  mark
 # Initial revision
 #
# Revision 1.2  1995/10/04  14:06:22  mikeq
# *** empty log message ***
#
 # Revision 1.1  1995/09/07 22:13:05  mikeq
 # Initial revision
 #
# Revision 1.1  95/05/25  11:39:36  11:39:36  mikeq (Mike Quesnel)
# Initial revision

 # Revision 1.12  1994/10/27  18:29:42  kennykb
 # Release 2.0 -- 10-27-94.  To be uploaded to archive sites.
 #
 # Revision 1.11  1993/11/01  18:20:46  kennykb
 # Beta release to be announced on comp.lang.tcl
 #
 # Revision 1.10  1993/10/27  15:52:49  kennykb
 # Package for alpha release to the Net, and for MMACE 931101 release.
 #
 # Revision 1.9  1993/10/20  19:10:47  kennykb
 # Alpha release #1 was thawed for bug fixes in tk 3.3.  Now frozen again at this
 # point.
 #
 # Revision 1.8  1993/10/20  18:40:24  kennykb
 # Repaired copyright notice so that it doesn't look like structured commentary.
 #
 # Revision 1.7  1993/10/14  18:15:42  kennykb
 # Cleaned up alignment of log messages, to avoid problems extracting
 # structured commentary.
 #
 # Revision 1.6  1993/10/14  18:06:59  kennykb
 # Added GE legal notice to head of file in preparation for release.
 #
 # Revision 1.5  1993/10/14  14:02:02  kennykb
 # Alpha release #1 frozen at this point.
 #
 # Revision 1.4  1993/07/21  19:44:36  kennykb
 # Finished cleaning up structured commentary.
 #
 # Revision 1.3  1993/07/20  13:12:03  kennykb
 # Made `choicebox', `collapsible', and `debug' conform with naming and
 # commentary conventions
 #
 # Revision 1.2  1993/07/16  15:58:00  kennykb
 # Renamed all commands that start with `wiget.' to either `widget_' or
 # `widget:'.  Added the code for creating composite widgets.
 #
 # Revision 1.1  1993/06/03  15:26:18  kennykb
 # Initial revision
 #

# Procedure: collapsible
#
# Synopsis:
#	Create a `collapsible' widget which may be added to or removed from
#	a larger display under user control.
#
# Usage:
#c	collapsible widgetType pathName ?option value?...
#
# Parameters:
#c	widgetType
#		A widget creation command such as `button' or `listbox'.
#
#c	pathName
#		Name of the frame in which the collapsible widget will be
#		created.
#
# Options:
#	Name:			title
#	Class:			Title
#	Command-Line Switch:	-title
#		Title of the collapsible widget.  Default is the null string.
#
#	Name:			visible
#	Class:			Visible
#	Command-Line Switch:	-visible
#		A Boolean value (0 or 1) giving the default state of the
#		visibility of the widget.  Default is 0.
#
#	Name:			command
#	Class:			Command
#	Command-Line Switch:	-command
#		Proc that will be executed when the frame is collapsed or
#		opened.
#
#	Other options are those accepted by the `widgetType' command.
#
# Description:
#	A `collapsible' widget is one that is packed into its parent's frame
#	only on demand.  It appears on the screen as an arrow (pointing right
#	or downward) at upper left, a title at upper center, and the widget 
#	itself at lower right.  The arrow is a button. Invoking the button
#	toggles the state of the widget between being visible and invisible
#	by packing it or unpacking it.
#
# Bugs:
#	- If a widget becomes invisible while it or one of its components has
#	the keyboard focus, it doesn't relinquish the focus.  It ought to.
#
#	- The collabsible processor does not honor the `configure' widget
#	command.
#
#	- The collapsible processor does not use the widget creation
#	definitions to make itself a first-class widget.

option add *Collapsible.visible 0 widgetDefault
option add *Collapsible.title {} widgetDefault
option add *Collapsible.command {} widgetDefault
option add *Collapsible.a.relief flat widgetDefault

set bitmapdir ${FSLDIR}/tcl

proc collapsible {type w args} {
	global collapsible_priv
	global bitmapdir

	frame $w -class Collapsible

	set title [option get $w title Title]
	set visible [option get $w visible Visible]
	set command [option get $w command Command]
	set fargs {}

	while {[llength $args] >= 2} {
		set option [lindex $args 0]
		set value [lindex $args 1]
		set args [lrange $args 2 end]
		case $option in {
			-title { set title $value }
			-visible { set visible $value }
			-command { set command $value }
			default {
				lappend fargs $option $value
			}
		}
	}

	set collapsible_priv(visible,$w) $visible
	set collapsible_priv(command,$w) $command

	pack append $w \
		[button $w.a \
			-command "collapsible_toggle $w" \
			-borderwidth 0 \
			-bitmap @$bitmapdir/triangler.xbm] \
				{left frame n} \
		[label $w.t -text $title] \
				{top frame w}

	eval [list $type $w.b] $fargs	

	if {$visible} {
		collapsible_show $w
	}
	return $w
}

# Procedure:	collapsible_toggle
#
# Synopsis:
#	Toggle the state of a collapsible widget when a user requests the
#	change.
#
# Usage:
#c	collapsible_toggle pathName
#
# Parameters:
#c	pathName
#		Path name of a collapsible widget
#
# Return value:
#	None.
#
# Description:
#	collapsible_toggle toggles the state of a collapsible widget.
#	If it is invisible, it makes it visible, and vice versa.

proc collapsible_toggle w {
	global collapsible_priv

	if {$collapsible_priv(visible,$w)} {
		collapsible_hide $w
	} else {
		collapsible_show $w
	}
	if {[string length $collapsible_priv(command,$w)]} {
		eval $collapsible_priv(command,$w) $collapsible_priv(visible,$w)
	}
}

# Procedure:	collapsible_hide
#
# Synopsis:
#	Make a collapsible widget invisible.
#
# Usage:
#c	collapsible_hide pathName
#
# Parameters:
#c	pathName
#		Path name of a collapsible widget
#
# Return value:
#	None.
#
# Description:
#	collapsible_hide makes a collapsible widget invisible.

proc collapsible_hide w {
	global collapsible_priv
	global bitmapdir
	set collapsible_priv(visible,$w) 0
	if {[string compare $w.a [focus]] == 0} {
		$w.a config -bitmap @$bitmapdir/triangler.xbm
	} else {
		$w.a config -bitmap @$bitmapdir/triangler.xbm
	}
	pack unpack $w.b
	# NOTE: Need to defocus any descendant of $w.b!
}

# Procedure: collapsible_show
#
# Synopsis:
#	Make a collapsible widget visible
#
# Usage:
#c	collapsible_show pathName
#
# Parameters:
#c	pathName
#		Path name of a collapsible widget.
#
# Return Value:
#	None.
#
# Description:
#	collapsible_show makes a collapsible widget visible.

proc collapsible_show w {
	global collapsible_priv
	global bitmapdir
	set collapsible_priv(visible,$w) 1
	if {[string compare $w.a [focus]] == 0} {
		$w.a config -bitmap @$bitmapdir/filltriangled.xbm
	} else {
		$w.a config -bitmap @$bitmapdir/triangled.xbm
	}
	pack append $w $w.b {top expand fill}
}

