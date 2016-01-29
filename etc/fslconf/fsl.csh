# FSL configuration file 
#  - to be sourced by the user, typically in .bashrc or equivalent
#  - note that the user should set 

# Written by Mark Jenkinson
#  FMRIB Analysis Group, University of Oxford

# SHBASECOPYRIGHT


#### Set up standard FSL user environment variables ####

# The following variable selects the default output image type
# Legal values are:  ANALYZE  NIFTI  NIFTI_PAIR  ANALYZE_GZ  NIFTI_GZ  NIFTI_PAIR_GZ
# This would typically be overwritten in ${HOME}/.fslconf/fsl.csh if the user wished
#  to write files with a different format
setenv FSLOUTPUTTYPE NIFTI_GZ

# Comment out the definition of FSLMULTIFILEQUIT to enable 
#  FSL programs to soldier on after detecting multiple image
#  files with the same basename ( e.g. epi.hdr and epi.nii )
setenv FSLMULTIFILEQUIT TRUE


# The following variables specify paths for programs and can be changed
#  or replaced by different programs ( e.g. FSLDISPLAY=open   for MacOSX)

setenv FSLTCLSH $FSLDIR/bin/fsltclsh
setenv FSLWISH $FSLDIR/bin/fslwish

# The following variables are used for running code in parallel across
#  several machines ( i.e. for FDT )

setenv FSLLOCKDIR 
setenv FSLMACHINELIST 
setenv FSLREMOTECALL 

# The following variables are used to configure CUDA capable queues - if you
# are using Grid Engine then this queue will be used to enqueue tasks that
# support execution on NVIDIA CUDA hardware.
setenv FSLGECUDAQ cuda.q

# Set up development variables (not for the faint-hearted)
# Uncomment the following if you wish to compile FSL source code
#setenv FSLCONFDIR $FSLDIR/config
#setenv FSLMACHTYPE `$FSLDIR/etc/fslconf/fslmachtype.sh`


###################################################
### Add other global environment variables here ###
###      or change the definitions above        ###
###################################################


# USER VARIABLES HERE


###################################################
####    DO NOT ADD ANYTHING BELOW THIS LINE    ####
###################################################

if ( -f /usr/local/etc/fslconf/fsl.csh ) then
  source /usr/local/etc/fslconf/fsl.csh
endif


if ( -f /etc/fslconf/fsl.csh ) then
  source /etc/fslconf/fsl.csh
endif


if ( -f "${HOME}/.fslconf/fsl.csh" ) then
  source "${HOME}/.fslconf/fsl.csh"
endif
