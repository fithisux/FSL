#{{{ copyright and setup 

#   FEAT TCL FSLIO wrappers
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 2004 University of Oxford
#
#   TCLCOPYRIGHT

#}}}

proc imcp { args } {
    global FSLDIR

    regsub -all "\{" $args "" cleanedargs
    regsub -all "\}" $cleanedargs "" cleanedargs

    return [ exec sh -c "${FSLDIR}/bin/imcp $args" ]
}

proc imglob { args } {
    global FSLDIR

    regsub -all "\{" $args "" cleanedargs
    regsub -all "\}" $cleanedargs "" cleanedargs

    return [ exec sh -c "${FSLDIR}/bin/imglob $args" ]
}

proc imln { args } {
    global FSLDIR

    regsub -all "\{" $args "" cleanedargs
    regsub -all "\}" $cleanedargs "" cleanedargs

    return [ exec sh -c "${FSLDIR}/bin/imln $args" ]
}

proc immv { args } {
    global FSLDIR

    regsub -all "\{" $args "" cleanedargs
    regsub -all "\}" $cleanedargs "" cleanedargs

    return [ exec sh -c "${FSLDIR}/bin/immv $args" ]
}

proc imrm { args } {
    global FSLDIR

    regsub -all "\{" $args "" cleanedargs
    regsub -all "\}" $cleanedargs "" cleanedargs

    return [ exec sh -c "${FSLDIR}/bin/imrm $args" ]
}

proc imtest { args } {
    global FSLDIR

    regsub -all "\{" $args "" cleanedargs
    regsub -all "\}" $cleanedargs "" cleanedargs

    return [ exec sh -c "${FSLDIR}/bin/imtest $args" ]
}

proc remove_ext { args } {
    global FSLDIR

    regsub -all "\{" $args "" cleanedargs
    regsub -all "\}" $cleanedargs "" cleanedargs

    if { $cleanedargs == "" } { return 1 }
    return [ exec sh -c "${FSLDIR}/bin/remove_ext $cleanedargs" ]
}

