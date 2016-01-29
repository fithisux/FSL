# FSL configuration file 
#  - to be sourced by the user, typically in .bashrc or equivalent
#  - note that the user should set 

# Written by Mark Jenkinson
#  FMRIB Analysis Group, University of Oxford

# SHBASECOPYRIGHT


#### Set up standard FSL user environment variables ####

# The following variable selects the default output image type
# Legal values are:  ANALYZE  NIFTI  NIFTI_PAIR  ANALYZE_GZ  NIFTI_GZ  NIFTI_PAIR_GZ
# This would typically be overwritten in ${HOME}/.fslconf/fsl.sh if the user wished
#  to write files with a different format
FSLOUTPUTTYPE=NIFTI_GZ
export FSLOUTPUTTYPE

# Comment out the definition of FSLMULTIFILEQUIT to enable 
#  FSL programs to soldier on after detecting multiple image
#  files with the same basename ( e.g. epi.hdr and epi.nii )
FSLMULTIFILEQUIT=TRUE ; export FSLMULTIFILEQUIT


# The following variables specify paths for programs and can be changed
#  or replaced by different programs ( e.g. FSLDISPLAY=open   for MacOSX)

FSLTCLSH=$FSLDIR/bin/fsltclsh
FSLWISH=$FSLDIR/bin/fslwish

export FSLTCLSH FSLWISH 

# The following variables are used for running code in parallel across
#  several machines ( i.e. for FDT )

FSLLOCKDIR=
FSLMACHINELIST=
FSLREMOTECALL=

export FSLLOCKDIR FSLMACHINELIST FSLREMOTECALL

# The following variables are used to configure CUDA capable queues - if you
# are using Grid Engine then this queue will be used to enqueue tasks that
# support execution on NVIDIA CUDA hardware.
FSLGECUDAQ=cuda.q

export FSLGECUDAQ

# Set up development variables (not for the faint-hearted)

FSLCONFDIR=$FSLDIR/config
FSLMACHTYPE=`$FSLDIR/etc/fslconf/fslmachtype.sh`

export FSLCONFDIR FSLMACHTYPE


###################################################
### Add other global environment variables here ###
###      or change the definitions above        ###
###################################################


# USER VARIABLES HERE


###################################################
####    DO NOT ADD ANYTHING BELOW THIS LINE    ####
###################################################

if [ -f /usr/local/etc/fslconf/fsl.sh ] ; then
  . /usr/local/etc/fslconf/fsl.sh ;
fi


if [ -f /etc/fslconf/fsl.sh ] ; then
  . /etc/fslconf/fsl.sh ;
fi


if [ -f "${HOME}/.fslconf/fsl.sh" ] ; then
  . "${HOME}/.fslconf/fsl.sh" ;
fi
