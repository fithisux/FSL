#* 
#* ------------------------------------------------------------------
#* BWExtras.tcl -- Assorted extra composite widgets
#* ------------------------------------------------------------------
#* labelspinbox.tcl - Labeled SpinBox
#* Created by Robert Heller on Wed Feb 15 21:21:45 2006
#* ------------------------------------------------------------------
#* Modification History: $Log: labelspinbox.tcl,v $
#* Modification History: Revision 1.1  2006/11/28 15:34:01  mwebster
#* Modification History: updates
#* Modification History:
#* Modification History: Revision 1.1.1.1  2006/02/16 14:58:07  heller
#* Modification History: Imported sources
#* Modification History:
#* ------------------------------------------------------------------
#* Contents:
#* ------------------------------------------------------------------
#*  
#* Copyright (c) 2006, Robert Heller
#* All rights reserved.
#* 
#* Redistribution and use in source and binary forms, with or without
#* modification, are permitted provided that the following conditions are
#* met:
#* 
#*     * Redistributions of source code must retain the above copyright
#*       notice, this list of conditions and the following disclaimer.
#*     * Redistributions in binary form must reproduce the above copyright
#*       notice, this list of conditions and the following disclaimer in the
#*       documentation and/or other materials provided with the distribution.
#*     * Neither the name of the Deepwoods Software nor the names of its
#*       contributors may be used to endorse or promote products derived from
#*       this software without specific prior written permission.
#* 
#* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
#* IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
#* TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#* PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#* 
#* 
# ------------------------------------------------------------------------------
#  $Id: labelspinbox.tcl,v 1.1 2006/11/28 15:34:01 mwebster Exp $
# ------------------------------------------------------------------------------
#  Index of commands:
#     - LabelSpinBox::create
#     - LabelSpinBox::configure
#     - LabelSpinBox::cget
#     - LabelSpinBox::bind
#     - LabelSpinBox::setvalue
#     - LabelSpinBox::getvalue
# ------------------------------------------------------------------------------

namespace eval LabelSpinBox {
    Widget::define LabelSpinBox labelspinbox SpinBox LabelFrame

    Widget::bwinclude LabelSpinBox LabelFrame .labf \
        remove {-relief -borderwidth -focus} \
        rename {-text -label} \
        prefix {label -justify -width -anchor -height -font -textvariable}

    Widget::bwinclude LabelSpinBox SpinBox .spin \
        remove {-fg -bg} \
        rename {-foreground -spinboxfg -background -spinboxbg}

    Widget::addmap LabelSpinBox "" :cmd {-background {}}

    Widget::syncoptions LabelSpinBox SpinBox .spin {-range -values {}}
    Widget::syncoptions LabelSpinBox LabelFrame .labf {-label -text -underline {}}

    ::bind BwLabelSpinBox <FocusIn> [list focus %W.labf]
    ::bind BwLabelSpinBox <Destroy> [list LabelSpinBox::_destroy %W]
}


# ------------------------------------------------------------------------------
#  Command LabelSpinBox::create
# ------------------------------------------------------------------------------
proc LabelSpinBox::create { path args } {
    array set maps [list LabelSpinBox {} :cmd {} .labf {} .spin {}]
    array set maps [Widget::parseArgs LabelSpinBox $args]

    eval [list frame $path] $maps(:cmd) -class LabelSpinBox \
	    -relief flat -bd 0 -highlightthickness 0 -takefocus 0
    Widget::initFromODB LabelSpinBox $path $maps(LabelSpinBox)
	
    set labf  [eval [list LabelFrame::create $path.labf] $maps(.labf) \
                   [list -relief flat -borderwidth 0 -focus $path.spin]]
    set subf  [LabelFrame::getframe $labf]
    set spin [eval [list SpinBox::create $path.spin] $maps(.spin)]

    pack $spin -in $subf -fill both -expand yes
    pack $labf  -fill both -expand yes

    bindtags $path [list $path BwLabelSpinBox [winfo toplevel $path] all]

    return [Widget::create LabelSpinBox $path]
}


# ------------------------------------------------------------------------------
#  Command LabelSpinBox::configure
# ------------------------------------------------------------------------------
proc LabelSpinBox::configure { path args } {
    return [Widget::configure $path $args]
}


# ------------------------------------------------------------------------------
#  Command LabelSpinBox::cget
# ------------------------------------------------------------------------------
proc LabelSpinBox::cget { path option } {
    return [Widget::cget $path $option]
}


# ------------------------------------------------------------------------------
#  Command LabelSpinBox::setvalue
# ------------------------------------------------------------------------------
proc LabelSpinBox::setvalue { path args } {
    return [eval [list ::setvalue $path.spin] $args]
}

# ------------------------------------------------------------------------------
#  Command LabelSpinBox::getvalue
# ------------------------------------------------------------------------------
proc LabelSpinBox::getvalue { path args } {
    return [eval [list ::getvalue $path.spin] $args]
}

# ------------------------------------------------------------------------------
#  Command LabelSpinBox::bind
# ------------------------------------------------------------------------------
proc LabelSpinBox::bind { path args } {
    return [eval [list ::bind $path.spin] $args]
}


#------------------------------------------------------------------------------
#  Command LabelSpinBox::_path_command
#------------------------------------------------------------------------------
proc LabelSpinBox::_path_command { path cmd larg } {
    if { [string equal $cmd "configure"] ||
         [string equal $cmd "cget"] ||
         [string equal $cmd "bind"] ||
	 [string equal $cmd "setvalue"] ||
	 [string equal $cmd "getvalue"]} {
        return [eval [list LabelSpinBox::$cmd $path] $larg]
    } else {
        return [eval [list $path.e:cmd $cmd] $larg]
    }
}


proc LabelSpinBox::_destroy { path } {
    Widget::destroy $path
}

package provide BWLabelSpinBox 1.0.0
