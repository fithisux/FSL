/*  fibre_gpu.h

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

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

#ifndef __FIBRE_GPU_
#define __FIBRE_GPU_

class FibreGPU{
public:  
    	float m_th;
	float m_th_prop;
	float m_th_prior;
	int m_th_acc; 
    	int m_th_rej;

    	float m_ph;
	float m_ph_prop;
	float m_ph_prior;
	int m_ph_acc;
    	int m_ph_rej; 

    	float m_f;
	float m_f_prop;
	float m_f_prior;
	int m_f_acc;
    	int m_f_rej;

    	//float m_lam;
    	//float m_lam_prop;
    	//float m_lam_prior;
    
    	float m_prior_en;
    	bool m_lam_jump;
    	//float m_d;
};

class MultifibreGPU{
public:    
    
    	float m_f0;
    	float m_f0_prop;
    	float m_f0_prior;
	int m_f0_acc;
    	int m_f0_rej;

    	float m_tau;
    	float m_tau_prop;
    	float m_tau_prior;
	int m_tau_acc;
    	int m_tau_rej;

	float m_S0;
    	float m_S0_prop;
    	float m_S0_prior;
    	int m_S0_acc;
    	int m_S0_rej;

    	float m_d; 			
    	float m_d_prop;
    	float m_d_prior; 
    	int m_d_acc;
    	int m_d_rej;

    	float m_dstd;
    	float m_dstd_prop;
    	float m_dstd_prior;	 
    	int m_dstd_acc;
    	int m_dstd_rej;

	float m_R;
    	float m_R_prop;
    	float m_R_prior;	 
    	int m_R_acc;
    	int m_R_rej;

    	float m_prior_en;		
    	float m_likelihood_en;
    	float m_energy;
};

#endif
