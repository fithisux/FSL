/*  costfns.h

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2001 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 5.0 (c) 2012, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Isis
    Innovation Limited ("Isis"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    innovation@isis.ox.ac.uk quoting reference DE/9564. */


#if !defined(__costfns_h)
#define __costfns_h

#include "newimageall.h"

using namespace NEWIMAGE;

namespace NEWIMAGE {

  enum costfns { Woods, CorrRatio, MutualInfo, NormCorr, NormMI, LeastSq, LabelDiff,
		 NormCorrSinc, BBR, Unknown };

  costfns costfn_type(const string& cname);

  class Costfn {
  public:
    const volume<float> &refvol;
    const volume<float> &testvol;
    const volume<float> &rweight;
    const volume<float> &tweight;
    volume<float> wmseg;  // WM segmentation (or any useful boundary) for BBR
    volume<float> fmap;  // fieldmap for BBR
    volume<float> fmap_mask;  // fieldmap mask for BBR
    //volume4D<float> nonlin_basis;
    mutable volume4D<float> debugvol;
  private:
    int *bindex;
    int no_bins;
    ColumnVector plnp;
    int *jointhist;
    int *marghist1;
    int *marghist2;
    float *fjointhist;
    float *fmarghist1;
    float *fmarghist2;
    mutable int p_count;
    costfns p_costtype;
    bool validweights;
    float bin_a0;
    float bin_a1;
    float bbr_dist;  // in mm
    float bbr_offset;
    float bbr_slope;
    Matrix bbr_pts;  // in mm coords
    Matrix bbr_norms;  // in mm coords pointing from wm to other
    float *gm_coord_x;  // in mm coords
    float *gm_coord_y;
    float *gm_coord_z;
    float *wm_coord_x;
    float *wm_coord_y;
    float *wm_coord_z;
    int no_coords;
    int vertex_step;
    int pe_dir;  // for fieldmap application
    string bbr_type;
    bool debug_mode;
  public: 
    float smoothsize;
    float fuzzyfrac;

  public:
    // Publicly available calls
    Costfn(const volume<float>& refv, const volume<float>& inputv);
    Costfn(const volume<float>& refv, const volume<float>& inputv,
	   const volume<float>& refweight, const volume<float>& inweight);
    ~Costfn();

    void set_debug_mode(bool debug_flag=true);
    void set_costfn(const costfns& costtype) { p_costtype = costtype; }
    costfns get_costfn(void) { return p_costtype; }
    void set_no_bins(int n_bins);
    int set_bbr_seg(const volume<float>& bbr_seg);
    int set_bbr_coords(const Matrix& coords, const Matrix& norms);
    int set_bbr_type(const string& typenm);
    int set_bbr_step(int step);
    int set_bbr_slope(float slope);
    int set_bbr_fmap(const volume<float>& fieldmap, int phase_encode_direction);
    int set_bbr_fmap(const volume<float>& fieldmap, const volume<float>& fieldmap_mask, int phase_encode_direction);
    int count() const { return p_count; }

    // General cost function call
    float cost(const Matrix& affmat) const;    // affmat is voxel to voxel
    // affmat is voxel to voxel and non-linear parameters are arbitrary
    float cost(const Matrix& affmat, const ColumnVector& nonlin_params) const;
    // in the following, all warps are mm to mm
    float cost(const volume4D<float>& warp) const;
    float cost_gradient(volume4D<float>& gradvec,
			const volume4D<float>& warp, bool nullbc) const; 

    // some basic entropy calls
    float ref_entropy(const Matrix& aff) const;
    float test_entropy(const Matrix& aff) const;
    float joint_entropy(const Matrix& aff) const;

    volume<float> image_mapper(const Matrix& affmat) const;    // affmat is voxel to voxel
    Matrix mappingfn(const Matrix& affmat) const;    // affmat is voxel to voxel
    float get_bin_intensity(int bin_number) const;
    float get_bin_number(float intensity) const;
    bool is_bbr_set(void) const;

    // a resampling function (since it is logical to keep it with the general cost processing for bbr)
    float bbr_resamp(const Matrix& aff, const ColumnVector& nonlin_params, volume<float>& resampvol) const;

  private:
    // Prevent default behaviours
    Costfn();
    Costfn operator=(const Costfn&);
    Costfn(const Costfn&);

    // Internal functions available
    float normcorr(const Matrix& aff) const; 
    float normcorr_smoothed(const Matrix& aff) const;
    float normcorr_smoothed_sinc(const Matrix& aff) const;
    float normcorr_fully_weighted(const Matrix& aff,
				  const volume<float>& refweight, 
				  const volume<float>& testweight) const;
    
    float leastsquares(const Matrix& aff) const;
    float leastsquares_smoothed(const Matrix& aff) const;
    float leastsquares_fully_weighted(const Matrix& aff, 
				      const volume<float>& refweight, 
				      const volume<float>& testweight) const;
    
    float labeldiff(const Matrix& aff) const;
    float labeldiff_smoothed(const Matrix& aff) const;
    float labeldiff_fully_weighted(const Matrix& aff, 
				      const volume<float>& refweight, 
				      const volume<float>& testweight) const;
    
    float bbr(const Matrix& aff) const;
    float bbr(const Matrix& aff, const ColumnVector& nonlin_params) const;
    float bbr(const Matrix& aff, const ColumnVector& nonlin_params, 
	      volume<float>& resampvol, bool resampling_required) const;
    float fmap_extrap(const double& x_vox, const double& y_vox, const double& z_vox, const ColumnVector& v_pe) const;
    int vox_coord_calc(ColumnVector& tvc, ColumnVector& rvc, const Matrix& aff, const ColumnVector& nonlin_params, 
		       const Matrix& iaffbig, const Matrix& mm2vox, const ColumnVector& pe_dir_vec) const;

    float woods_fn(const Matrix& aff) const; 
    float woods_fn_smoothed(const Matrix& aff) const; 
    
    float corr_ratio(const Matrix& aff) const; 
    float corr_ratio_smoothed(const Matrix& aff) const; 
    float corr_ratio_fully_weighted(const Matrix& aff,
				    const volume<float>& refweight, 
				    const volume<float>& testweight) const;
    float corr_ratio_fully_weighted(const volume4D<float>& warpvol,
				    const volume<float>& refweight, 
				    const volume<float>& testweight) const;
    float corr_ratio_gradient_fully_weighted(volume4D<float>& gradvec,
					     const volume4D<float>& warpvol,
					     const volume<float>& refweight, 
					     const volume<float>& testweight,
					     bool nullbc) const;
     
    float mutual_info(const Matrix& aff) const;
    float mutual_info_smoothed(const Matrix& aff) const;
    float mutual_info_fully_weighted(const Matrix& aff,
				     const volume<float>& refweight, 
				     const volume<float>& testweight) const;
    
    float normalised_mutual_info(const Matrix& aff) const;
    float normalised_mutual_info_smoothed(const Matrix& aff) const;
    float normalised_mutual_info_fully_weighted(const Matrix& aff,
						const volume<float>& refweight, 
						const volume<float>& testweight) const;
    
    float cost(const Matrix& affmat,
	       const volume<float>& refweight, 
	       const volume<float>& testweight) const;

    float cost(const Matrix& affmat,
	       const ColumnVector& nonlin_params,
	       const volume<float>& refweight, 
	       const volume<float>& testweight) const;
    
    float cost(const volume4D<float>& warp,
	       const volume<float>& refweight, 
	       const volume<float>& testweight) const;
    
    float cost_gradient(volume4D<float>& gradvec,
			const volume4D<float>& warp,
			const volume<float>& refweight, 
			const volume<float>& testweight,
			bool nullbc=false) const;

    int p_corr_ratio_image_mapper(volume<float>& vout,
				  Matrix& mappingfn,
				  const volume<float>& vref, 
				  const volume<float>& vtest,
				  const volume<float>& refweight, 
				  const volume<float>& testweight,
				  int *bindex, const Matrix& aff,
				  const int no_bins, const float smoothsize) const;
    };


   //////////////////////////////////////////////////////////////////////////


}

#endif







