#* 
#* ------------------------------------------------------------------
#* labelcombobox.tcl - Labeled ComboBox
#* Created by Robert Heller on Thu Feb 16 10:03:17 2006
#* ------------------------------------------------------------------
#* Modification History: $Log: labelcombobox.tcl,v $
#* Modification History: Revision 1.1  2006/11/28 15:34:01  mwebster
#* Modification History: updates
#* Modification History:
#* Modification History: Revision 1.1  2006/02/16 15:20:02  heller
#* Modification History: Added LabelComboBox
#* Modification History:
#* Modification History: Revision 1.1  2002/07/28 14:03:50  heller
#* Modification History: Add it copyright notice headers
#* Modification History:
#* ------------------------------------------------------------------
#* Contents:
#* ------------------------------------------------------------------
#*  
#*     Generic Project
#*     Copyright (C) 2005  Robert Heller D/B/A Deepwoods Software
#* 			51 Locke Hill Road
#* 			Wendell, MA 01379-9728
#* 
#*     This program is free software; you can redistribute it and/or modify
#*     it under the terms of the GNU General Public License as published by
#*     the Free Software Foundation; either version 2 of the License, or
#*     (at your option) any later version.
#* 
#*     This program is distributed in the hope that it will be useful,
#*     but WITHOUT ANY WARRANTY; without even the implied warranty of
#*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#*     GNU General Public License for more details.
#* 
#*     You should have received a copy of the GNU General Public License
#*     along with this program; if not, write to the Free Software
#*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#* 
#*  
#* 
# ------------------------------------------------------------------------------
#  $Id: labelcombobox.tcl,v 1.1 2006/11/28 15:34:01 mwebster Exp $
# ------------------------------------------------------------------------------
#  Index of commands:
#     - LabelComboBox::create
#     - LabelComboBox::configure
#     - LabelComboBox::cget
#     - LabelComboBox::bind
#     - LabelComboBox::get
#     - LabelComboBox::getlistbox
#     - LabelComboBox::getvalue
#     - LabelComboBox::icursor
#     - LabelComboBox::post
#     - LabelComboBox::setvalue
#     - LabelComboBox::unpost
# ------------------------------------------------------------------------------

namespace eval LabelComboBox {
    Widget::define LabelComboBox labelcombobox ComboBox LabelFrame

    Widget::bwinclude LabelComboBox LabelFrame .labf \
        remove {-relief -borderwidth -focus} \
        rename {-text -label} \
        prefix {label -justify -width -anchor -height -font -textvariable}

    Widget::bwinclude LabelComboBox ComboBox .combo \
        remove {-fg -bg} \
        rename {-foreground -comboboxfg -background -comboboxbg 
		-height -comboboxheight -listboxwidth comboboxlistboxwidth}

    Widget::addmap LabelComboBox "" :cmd {-background {}}

    Widget::syncoptions LabelComboBox ComboBox .combo {-autocomplete -bwlistbox
						       -expand -hottrack 
						       -images -modifycmd
						       -postcommand -values {}}
    Widget::syncoptions LabelComboBox LabelFrame .labf {-label -text -underline {}}

    ::bind BwLabelComboBox <FocusIn> [list focus %W.labf]
    ::bind BwLabelComboBox <Destroy> [list LabelComboBox::_destroy %W]
}


# ------------------------------------------------------------------------------
#  Command LabelComboBox::create
# ------------------------------------------------------------------------------
proc LabelComboBox::create { path args } {
    array set maps [list LabelComboBox {} :cmd {} .labf {} .combo {}]
    array set maps [Widget::parseArgs LabelComboBox $args]

    eval [list frame $path] $maps(:cmd) -class LabelComboBox \
	    -relief flat -bd 0 -highlightthickness 0 -takefocus 0
    Widget::initFromODB LabelComboBox $path $maps(LabelComboBox)
	
    set labf  [eval [list LabelFrame::create $path.labf] $maps(.labf) \
                   [list -relief flat -borderwidth 0 -focus $path.combo]]
    set subf  [LabelFrame::getframe $labf]
    set combo [eval [list ComboBox::create $path.combo] $maps(.combo)]

    pack $combo -in $subf -fill both -expand yes
    pack $labf  -fill both -expand yes

    bindtags $path [list $path BwLabelComboBox [winfo toplevel $path] all]

    return [Widget::create LabelComboBox $path]
}


# ------------------------------------------------------------------------------
#  Command LabelComboBox::configure
# ------------------------------------------------------------------------------
proc LabelComboBox::configure { path args } {
    return [Widget::configure $path $args]
}


# ------------------------------------------------------------------------------
#  Command LabelComboBox::cget
# ------------------------------------------------------------------------------
proc LabelComboBox::cget { path option } {
    return [Widget::cget $path $option]
}


# ------------------------------------------------------------------------------
#  Command LabelComboBox::bind
# ------------------------------------------------------------------------------
proc LabelComboBox::bind { path args } {
    return [eval [list ::bind $path.combo] $args]
}

# ------------------------------------------------------------------------------
#  Command LabelComboBox::get
# ------------------------------------------------------------------------------
proc LabelComboBox::get { path args } {
    return [eval [list ::get $path.combo] $args]
}

# ------------------------------------------------------------------------------
#  Command LabelComboBox::getlistbox
# ------------------------------------------------------------------------------
proc LabelComboBox::getlistbox { path args } {
    return [eval [list ::getlistbox $path.combo] $args]
}

# ------------------------------------------------------------------------------
#  Command LabelComboBox::getvalue
# ------------------------------------------------------------------------------
proc LabelComboBox::getvalue { path args } {
    return [eval [list ::getvalue $path.combo] $args]
}

# ------------------------------------------------------------------------------
#  Command LabelComboBox::icursor
# ------------------------------------------------------------------------------
proc LabelComboBox::icursor { path args } {
    return [eval [list ::icursor $path.combo] $args]
}

# ------------------------------------------------------------------------------
#  Command LabelComboBox::post
# ------------------------------------------------------------------------------
proc LabelComboBox::post { path args } {
    return [eval [list ::post $path.combo] $args]
}

# ------------------------------------------------------------------------------
#  Command LabelComboBox::setvalue
# ------------------------------------------------------------------------------
proc LabelComboBox::setvalue { path args } {
    return [eval [list ::setvalue $path.combo] $args]
}

# ------------------------------------------------------------------------------
#  Command LabelComboBox::unpost
# ------------------------------------------------------------------------------
proc LabelComboBox::unpost { path args } {
    return [eval [list ::unpost $path.combo] $args]
}


#------------------------------------------------------------------------------
#  Command LabelComboBox::_path_command
#------------------------------------------------------------------------------
proc LabelComboBox::_path_command { path cmd larg } {
    if { [string equal $cmd "configure"] ||
         [string equal $cmd "cget"] ||
         [string equal $cmd "bind"] ||
	 [string equal $cmd "get"] ||
	 [string equal $cmd "getlistbox"] ||
	 [string equal $cmd "getvalue"] ||
	 [string equal $cmd "icursor"] ||
	 [string equal $cmd "post"] ||
	 [string equal $cmd "setvalue"] ||
	 [string equal $cmd "unpost"]} {
        return [eval [list LabelComboBox::$cmd $path] $larg]
    } else {
        return [eval [list $path.e:cmd $cmd] $larg]
    }
}


proc LabelComboBox::_destroy { path } {
    Widget::destroy $path
}

package provide BWLabelComboBox 1.0.0
