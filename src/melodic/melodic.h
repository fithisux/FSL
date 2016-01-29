/*  MELODIC - Multivariate exploratory linear optimized decomposition into 
                 independent components
    
    melodic.h - main program header

    Christian F. Beckmann and Matthew Webster, FMRIB Analysis Group
    
    Copyright (C) 1999-2013 University of Oxford */

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

#ifndef __MELODIC_h
#define __MELODIC_h

#include <strstream>

#ifdef __APPLE__
#include <mach/mach.h>
#define memmsg(msg) { \
  MelodicOptions&opt = MelodicOptions::getInstance(); \
  struct task_basic_info t_info; \
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT; \
  if (KERN_SUCCESS == task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t) &t_info, &t_info_count)) \
	{ \
		if(opt.debug.value()) {\
		    cout << " MEM: " << msg << " res: " << t_info.resident_size/1000000 << " virt: " << t_info.virtual_size/1000000 << "\n"; \
		    cout.flush(); \
		 } \
	} \
}
#else
#define memmsg(msg) { \
  MelodicOptions&opt = MelodicOptions::getInstance(); \
  if(opt.debug.value()) {\
		cout << msg; \
		cout.flush(); \
  } \
}
#endif



// a simple message macro that takes care of cout and log
#define message(msg) { \
  MelodicOptions& opt = MelodicOptions::getInstance(); \
  if(opt.verbose.value()) \
    { \
      cout << msg; \
    } \
  Log& logger  =   LogSingleton::getInstance(); \
  logger.str() << msg; \
  cout.flush(); \
}

#define dbgmsg(msg) { \
  MelodicOptions&opt = MelodicOptions::getInstance(); \
  if(opt.debug.value()) {\
		cout << msg; } \
}

#define outMsize(msg,Mat) { \
  MelodicOptions& opt = MelodicOptions::getInstance();		\
  if(opt.debug.value())						\
    cerr << "     " << msg << "  " <<Mat.Nrows() << " x " << Mat.Ncols() << endl;	\
}

namespace Melodic{

const string version = "3.14";  

// The two strings below specify the title and example usage that is	
// printed out as the help or usage message
const string title=string("MELODIC (Version ")+version+")"+
		string("\n Multivariate Exploratory Linear Optimised Decomposition into Independent Components\n")+
		string("\nAuthor: Christian F. Beckmann \n Copyright(c) 2001-2013 University of Oxford");

const string usageexmpl=string(" melodic -i <filename> <options>")+
		   string("\n \t \t to run melodic")+
	//	   string("\n melodic -i <filename> --mix=melodic_mix")+
	//	   string(" --filter=\"string of component numbers\"")+
	//	   string("\n \t \t to remove estimated ICs from input")+
		   string("\n melodic -i <filename> --ICs=melodic_IC")+
		   string(" --mix=melodic_mix <options>")+
		   string("\n \t \t to run Mixture Model based inference on estimated ICs")+
		   string("\n melodic --help ");
}

#endif
