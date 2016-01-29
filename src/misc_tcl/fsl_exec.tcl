#

#   fsl:exec execute a shell command with logging all setup
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 2000-2007 University of Oxford
#
#   Part of FSL - FMRIB's Software Library
#   http://www.fmrib.ox.ac.uk/fsl
#   fsl@fmrib.ox.ac.uk
#   
#   Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
#   Imaging of the Brain), Department of Clinical Neurology, Oxford
#   University, Oxford, UK
#   
#   
#   LICENCE
#   
#   FMRIB Software Library, Release 5.0 (c) 2012, The University of
#   Oxford (the "Software")
#   
#   The Software remains the property of the University of Oxford ("the
#   University").
#   
#   The Software is distributed "AS IS" under this Licence solely for
#   non-commercial use in the hope that it will be useful, but in order
#   that the University as a charitable foundation protects its assets for
#   the benefit of its educational and research purposes, the University
#   makes clear that no condition is made or to be implied, nor is any
#   warranty given or to be implied, as to the accuracy of the Software,
#   or that it will be suitable for any particular purpose or for use
#   under any specific conditions. Furthermore, the University disclaims
#   all responsibility for the use which is made of the Software. It
#   further disclaims any liability for the outcomes arising from using
#   the Software.
#   
#   The Licensee agrees to indemnify the University and hold the
#   University harmless from and against any and all claims, damages and
#   liabilities asserted by third parties (including claims for
#   negligence) which arise directly or indirectly from the use of the
#   Software or the sale of any products based on the Software.
#   
#   No part of the Software may be reproduced, modified, transmitted or
#   transferred in any form or by any means, electronic or mechanical,
#   without the express permission of the University. The permission of
#   the University is not required if the said reproduction, modification,
#   transmission or transference is done without financial return, the
#   conditions of this Licence are imposed upon the receiver of the
#   product, and all original and amended source code is included in any
#   transmitted product. You may be held legally responsible for any
#   copyright infringement that is caused or encouraged by your failure to
#   abide by these terms and conditions.
#   
#   You are not permitted under this Licence to use this Software
#   commercially. Use for which any financial return is received shall be
#   defined as commercial use, and includes (1) integration of all or part
#   of the source code or the Software into a product for sale or license
#   by or on behalf of Licensee to third parties or (2) use of the
#   Software or any derivative of it for research with the final aim of
#   developing software products for sale or license to a third party or
#   (3) use of the Software or any derivative of it for research with the
#   final aim of developing non-software products for sale or license to a
#   third party, or (4) use of the Software to provide any service to an
#   external organisation for which payment is received. If you are
#   interested in using the Software commercially, please contact Isis
#   Innovation Limited ("Isis"), the technology transfer company of the
#   University, to negotiate a licence. Contact details are:
#   innovation@isis.ox.ac.uk quoting reference DE/9564.

#
# usage: fsl:echo <filename> "string" [-o]
# -o : overwrite instead of appending
#
proc fsl:echo { thefilename thestring args } {
    if { $thefilename != "" } {
	if { [ llength $args ] != 0 } {
	    catch { exec sh -c "echo \"$thestring\" >  $thefilename" >& /dev/null } putserr
	} else {
	    catch { exec sh -c "echo \"$thestring\" >> $thefilename" >& /dev/null } putserr
	}
    } else {
	catch { puts $thestring } putserr
    }
    if { [ glob -nocomplain logs/* ] != "" } {
	catch { exec sh -c "cat logs/* > report_log.html" } putserr
    }
}


# usage: fsl:exec "the command" [options]
# note the double-quotes round "the command"
# -n                : don't output logging information
# -b <job_duration> : submit to SGE (if available) and return immediately; otherwise wait to return until <command> has finished
#                     job_duration is in minutes, set 0 for unknown
# -f                : "the command" is a filename to be treated as a parallel batch job (one command per line)
# -h <job-hold-id>  : ID of SGE job that must finish before this job gets run
# -N <job_name>     : the job will have the name <job_name> on the SGE submission
# -l <logdir>       : output logs in <logdir>
# -i                : ignore errors ( return exit status 0 to calling script )

proc fsl:exec { thecommand args } {

    global logout comout FSLDIR FD SGEID

    # process args
    set do_logout 1
    set do_runfsl 0
    set do_runfsl_command ""
    set do_parallel ""
    set job_holds ""
    set job_name ""
    set runtime 300
    set ignoreErrors 0
    set args [ join $args ]
    set logdir ""
    for { set argindex 0 } { $argindex < [ llength $args ] } { incr argindex 1 } {
	set thearg [ lindex $args $argindex ]

	if { $thearg == "-n" } {
	    set do_logout 0
	}

	if { $thearg == "-f" } {
	    set do_parallel "-t"
	}

	if { $thearg == "-h" } {
	    incr argindex 1
	    set theid [ lindex $args $argindex ]
	    if { $theid > 0 } {
		if { $job_holds == "" } {
		    set job_holds "-j $theid"
		}  else {
		    set job_holds "${job_holds},$theid"
		}
	    }
	}

	if { $thearg == "-N" } {
	    incr argindex 1
	    set job_name "-N [ lindex $args $argindex ]"
	}

	if { $thearg == "-l" } {
	    incr argindex 1
	    set logdir "-l [ lindex $args $argindex ]"
	}

	if { $thearg == "-b" } {
	    set do_runfsl 1
	    incr argindex 1
	    set runtime [ lindex $args $argindex ]
	}

	if { $thearg == "-i" } {
	    set ignoreErrors 1
	}


    }

    # add runfsl call if required
    if { $do_runfsl } {
	set thecommand "${FSLDIR}/bin/fsl_sub -T $runtime $logdir $job_name $job_holds $do_parallel $thecommand"
    }

    # save just the command to report.com if it has been setup
    if { [ info exists comout ] && $comout != "" } {
	fsl:echo $comout "$thecommand\n"
    }

    # if logout doesn't exist, set empty so that logging will go to screen
    if { ! [ info exists logout ] } {
	set logout ""
    }

    # run and log the actual command
    if { $do_logout } {
	fsl:echo $logout "\n$thecommand"
	if { $logout != "" } {
	    set errorCode [ catch { exec sh -c $thecommand >>& $logout } errmsg ]
	    set actualCode [lindex $::errorCode 0]
	    if { $errorCode != 0 && $actualCode ne "NONE" &&  $ignoreErrors != 1 } { 
	    } else {
	    # Trim errmsg to final line
		set errmsg [ exec sh -c "tail -n 1 $logout" ]
	    }
	} else {
	    set errorCode [ catch { exec sh -c $thecommand } errmsg ]
	    set actualCode [lindex $::errorCode 0]
	    catch { puts $errmsg } putserr
	}
    } else {
	set errorCode [ catch { exec sh -c $thecommand } errmsg ]
	set actualCode [lindex $::errorCode 0]
    }

    if { $errorCode != 0 && $actualCode ne "NONE" &&  $ignoreErrors != 1 } { 
	if { $do_logout } {
	    fsl:echo $logout "\nFATAL ERROR ENCOUNTERED:\nCOMMAND:\n$thecommand\nERROR MESSAGE:\n$errmsg\nEND OF ERROR MESSAGE" 	
	}
    } else {
       	set errorCode 0
    }


    # now return errmsg in case the exec call needs to know the output
    return -code $errorCode $errmsg 
}

