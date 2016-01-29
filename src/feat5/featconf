
#########################################################################
# Misc
#########################################################################

# Balloon help
set fmri(help_yn) 1

# Run Featwatcher
set fmri(featwatcher_yn) 1

# Brain/background threshold, %
set fmri(brain_thresh) 10

# Critical z for design efficiency calculation
set fmri(critical_z) 5.3

# Noise level
set fmri(noise) 0.66

# Noise AR(1)
set fmri(noisear) 0.34

# 0 : Overwrite original post-stats results
# 1 : Copy original FEAT directory for new Contrasts, Thresholding, Rendering
set fmri(newdir_yn) 0

# Cleanup first-level standard-space images
set fmri(sscleanup_yn) 0


#########################################################################
# Data
#########################################################################

# Delete volumes
set fmri(ndelete) 0

# TR(s)
set fmri(tr) 3.0

# High pass filter cutoff
set fmri(paradigm_hp) 100

#########################################################################
# Pre-stats
#########################################################################

# Motion correction
# 0 : None
# 1 : MCFLIRT
set fmri(mc) 1

# B1 unwarping
set fmri(regunwarp_yn) 0

# Alternate reference
set fmri(alternateReference_yn) 0

# EPI dwell time (ms)
set fmri(dwell) 0.7

# EPI TE (ms)
set fmri(te) 35

# % Signal loss threshold
set fmri(signallossthresh) 10

# Unwarp direction
set fmri(unwarp_dir) y-

# Slice timing correction
# 0 : None
# 1 : Regular up (0, 1, 2, 3, ...)
# 2 : Regular down
# 3 : Use slice order file
# 4 : Use slice timings file
# 5 : Interleaved (0, 2, 4 ... 1, 3, 5 ... )
set fmri(st) 0

# BET brain extraction
set fmri(bet_yn) 1

# Spatial smoothing FWHM (mm)
set fmri(smooth) 5

# Intensity normalization
set fmri(norm_yn) 0

# Perfusion subtraction
set fmri(perfsub_yn) 0

# Perfusion tag/control order
set fmri(tagfirst) 1

# Highpass temporal filtering
set fmri(temphp_yn) 1

# Lowpass temporal filtering
set fmri(templp_yn) 0

# MELODIC ICA data exploration
set fmri(melodic_yn) 0


#########################################################################
# Stats
#########################################################################

# Carry out prewhitening?
set fmri(prewhiten_yn) 1

# Add motion parameters to model
# 0 : No
# 1 : Yes
set fmri(motionevs) 0
set fmri(motionevsbeta) ""

# General confound matrix (any number of EVs)
set fmri(confoundevs) 0

# Robust outlier detection in FLAME
set fmri(robust_yn) 0

# Higher-level modelling
# 3 : Fixed effects
# 0 : Mixed Effects: Simple OLS
# 2 : Mixed Effects: FLAME 1
# 1 : Mixed Effects: FLAME 1+2
set fmri(mixed_yn) 2

# HRF Convolution
# 0 : None
# 1 : Gaussian
# 2 : Gamma
# 3 : Double-Gamma HRF
# 4 : Gamma basis functions
# 5 : Sine basis functions
# 6 : FIR basis functions
# 7 : Optimal/custom basis functions
set fmri(default_convolve) 2

# Convolve phase
set fmri(default_convolve_phase) 0

# Gauss sigma
set fmri(default_gausssigma) 2.8

# Gauss delay
set fmri(default_gaussdelay) 5

# Gamma sigma
set fmri(default_gammasigma) 3

# Gamma delay
set fmri(default_gammadelay) 6

# Optimal/custom HRF convolution file
set fmri(default_bfcustom) "${FSLDIR}/etc/default_flobs.flobs/hrfbasisfns.txt"

# Add temporal derivative
set fmri(default_deriv_yn) 1

# Number of additional voxel-dependent EVs
set fmri(evs_vox) 0

#########################################################################
# Post-stats
#########################################################################

# Thresholding
# 0 : None
# 1 : Uncorrected
# 2 : Voxel
# 3 : Cluster
set fmri(thresh) 3

# Contrast masking - use >0 instead of thresholding?
set fmri(conmask_zerothresh_yn) 0

# P threshold
set fmri(prob_thresh) 0.05

# Z threshold
set fmri(z_thresh) 2.3

# Z min/max for colour rendering
# 0 : Use actual Z min/max
# 1 : Use preset Z min/max
set fmri(zdisplay) 0

# Z min in colour rendering
set fmri(zmin) 2

# Z max in colour rendering
set fmri(zmax) 8

# Colour rendering type
# 0 : Solid blobs
# 1 : Transparent blobs
set fmri(rendertype) 1

# Background image for higher-level stats overlays
# 1 : Mean highres
# 2 : First highres
# 3 : Mean functional
# 4 : First functional
# 5 : Standard space template
set fmri(bgimage) 1

# Create time series plots
set fmri(tsplot_yn) 1

#########################################################################
# Registration
#########################################################################

set fmri(reginitial_highres_yn) 0
set fmri(reghighres_yn) 0
set fmri(regstandard_yn) 1

# Search space for registration to initial structural
# 0   : No search
# 90  : Normal search
# 180 : Full search
set fmri(reginitial_highres_search) 90

# Degrees of Freedom for registration to initial structural
# options: 3/6/7/9/12
set fmri(reginitial_highres_dof) 3

# Search space for registration to main structural
set fmri(reghighres_search) 90

# Degrees of Freedom for registration to main structural
set fmri(reghighres_dof) BBR

# Standard image
set fmri(regstandard) "${FSLDIR}/data/standard/MNI152_T1_2mm_brain"

# Search space for registration to standard space
set fmri(regstandard_search) 90

# Degrees of Freedom for registration to standard space
set fmri(regstandard_dof) 12

# Do nonlinear registration from structural to standard space?
set fmri(regstandard_nonlinear_yn) 0

# Control nonlinear warp field resolution
set fmri(regstandard_nonlinear_warpres) 10


#########################################################################
# Non-GUI options
#########################################################################

# Alternative example_func image (not derived from input 4D dataset)
set fmri(alternative_example_func) ""

# Alternative (to BETting) mask image
set fmri(alternative_mask) ""

# Initial structural space registration initialisation transform
set fmri(init_initial_highres) ""

# Structural space registration initialisation transform
set fmri(init_highres) ""

# Standard space registration initialisation transform
set fmri(init_standard) ""

# For full FEAT analysis: overwrite existing .feat output dir?
set fmri(overwrite_yn) 0

# For custom FNIRT configuration change the following
set fmri(fnirt_config) "T1_2_MNI152_2mm"

#########################################################################

