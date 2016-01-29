/*  prelude.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000/2001 University of Oxford  */

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

#include "unwarpfns.h"
#include "utils/options.h"

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

using namespace Utilities;

string title="prelude (Version 2.0 in c# minor)\nPhase Region Expanding Labeller for Unwrapping Discrete Estimates\nCopyright(c) 2001, University of Oxford (Mark Jenkinson)";
string examples="prelude -c <rawphase> -o <unwrappedphase> [options]\nprelude -p <phasevol> -a <absvol> -o <unwrappedphase> [options]";

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> labelslices(string("--labelslices"), false,
		    string("does label processing in 2D (slice at a time)"),
		    false, no_argument);
Option<bool> allslices(string("-s,--slices"), false,
		       string("does all processing in 2D (slice at a time)"),
		       false, no_argument);
Option<bool> force3D(string("-f,--force3D"), false,
		     string("forces all processing to be full 3D"),
		     false, no_argument);
Option<bool> removeramps(string("--removeramps"), false,
		     string("remove phase ramps during unwrapping"),
		     false, no_argument);
Option<int> num_phase_split(string("-n,--numphasesplit"), 8,
			      string("number of phase partitions to use"),
			      false, requires_argument);
Option<int> startimno(string("--start"), 0,
		      string("first image number to process (default 0)"),
		      false, requires_argument);
Option<int> endimno(string("--end"), 0,
		    string("final image number to process (default Inf)"),
		    false, requires_argument);
Option<float> thresh(string("-t,--thresh"), 0,
		    string("intensity threshold for masking"),
		    false, requires_argument);
Option<string> complexvol(string("-c,--complex"), string(""),
			  string("filename of complex phase input volume"),
			  false, requires_argument);
Option<string> absvol(string("-a,--abs"), string(""),
		      string("filename of absolute input volume"),
		      false, requires_argument);
Option<string> phasevol(string("-p,--phase"), string(""),
			string("filename of raw phase input volume"),
			false, requires_argument);
Option<string> uphasevol(string("-o,--out,-u,--unwrap"), string(""),
			 string("filename for saving the unwrapped phase output"),
			 true, requires_argument);
Option<string> rawphasevol(string("-r,--rawphase"), string(""),
			   string("filename for saving the raw phase output"),
			   false, requires_argument);
Option<string> labelvol(string("-l,--labels"), string(""),
			   string("filename for saving the area labels output"),
			   false, requires_argument);
Option<string> maskvol(string("--savemask"), string(""),
			   string("filename for saving the mask volume"),
			   false, requires_argument);
Option<string> inmask(string("-m,--mask"), string(""),
		      string("filename of mask input volume"),
		      false, requires_argument);

////////////////////////////////////////////////////////////////////////////


void wrap(volume<float>& phasevol)
{
  for (int z=phasevol.minz(); z<=phasevol.maxz(); z++) {
    for (int y=phasevol.miny(); y<=phasevol.maxy(); y++) {
      for (int x=phasevol.minx(); x<=phasevol.maxx(); x++) {
	phasevol(x,y,z) = wrap(phasevol(x,y,z));
      }
    }
  }
}


int do_unwrapping() 
{

  volume4D<float> phasemaps, absmaps, masks;

  if (complexvol.set()) {
    volume4D<float> rvol, ivol;
    read_complexvolume4D(rvol,ivol,complexvol.value());
    for (int n=0; n<rvol.tsize(); n++) {
      phasemaps.addvolume(phase(rvol[n],ivol[n]));
      absmaps.addvolume(abs(rvol[n],ivol[n]));
    }
    phasemaps.setTR(rvol.TR());
  } else {
    if (verbose.value()) cout << "Loading volumes" << endl;
    read_volume4D(phasemaps,phasevol.value());
    if (verbose.value()) cout << "Phase loaded" << endl;
    read_volume4D(absmaps,absvol.value());
    if (verbose.value()) cout << "Magnitude loaded" << endl;
  }

  bool threshset = thresh.set();
  float threshval = thresh.value();
  int imgstart=0, imgend, maskend=0;
  imgend = phasemaps.maxt(); 
  if (startimno.set())  imgstart = startimno.value();
  if (endimno.set())    imgend = endimno.value();
  if (inmask.set()) {
    read_volume4D(masks,inmask.value());
    if (verbose.value()) cout << "Mask loaded" << endl;
    maskend = masks.maxt();
    if (!threshset) {
      threshset = true;
      threshval = 0.5;
    }
  }

  // sanity check on phase values
  for (int n=imgstart; n<=imgend; n++) {
    if ((phasemaps[n].max() - phasemaps[n].min())>3.0*M_PI) {
      cerr << "ERROR: input phase image exceeds allowable phase range."
	   << endl << "Allowable range is 6.283 radians.  Image range is: " 
	   << phasemaps[n].max() - phasemaps[n].min() << " radians."
	   << endl << "Aborting." << endl;
      return -1;
    }
    if ( (phasemaps[n].max()>1.001*M_PI) || (phasemaps[n].min()<-1.001*M_PI) )
      {
	if (verbose.value()) 
	  { cout << "Rewrapping phase range to [-pi,pi]" << endl; }
	wrap(phasemaps[n]);
      }
  }
  

  bool label_2D = labelslices.value();
  bool unwrap_2D = allslices.value();
  // The automatic behaviour is to use 2D labels for large images
  if ((phasemaps[imgstart].xsize() > 64) || (phasemaps[imgstart].ysize() > 64))
    {
      label_2D = true;
    }
  if (force3D.value()) { label_2D = false;  unwrap_2D = false; }

  volume4D<float> uphase;

  for (int n=imgstart; n<=imgend; n++) {

    // Make mask
    volume<float> mask;
    if (inmask.set()) {
      mask = masks[Min(n,maskend)];
    } else {
      mask = absmaps[n];
    }

    if (!threshset) { threshval = basic_mask_threshold(mask); }

    if (unwrap_2D) {
      mask = make_head_mask2D(mask,threshval);
    } else {
      mask = make_head_mask(mask,threshval);
    }
    if (maskvol.set()) { save_volume(mask,maskvol.value()); }

    // Remove phase ramps
    ColumnVector ramps;
    if (removeramps.value()) {
      if (verbose.value()) { cout << "Removing linear ramps" << endl; }
      ramps = estimate_linear_ramps(phasemaps[n],mask);
      if (verbose.value()) { cout << "Estimated linear ramps" << endl; }
      remove_linear_ramps(ramps,phasemaps[n],mask);
      if (verbose.value()) { cout << "Removed linear ramps" << endl; }
      if (verbose.value()) { cout << "Ramp values: " << ramps << endl; }
      if (verbose.value()) { save_volume(phasemaps[n],"grot"); }
    }

    // Make labels
    volume<int> label;
    if (verbose.value()) 
      {cout << "Number of phase splits = " << num_phase_split.value() << endl;}
    if (unwrap_2D) {
      label = find_phase_labels2D(phasemaps[n],mask,
				  num_phase_split.value(),false);
    } else if (label_2D) {
      label = find_phase_labels2D(phasemaps[n],mask,
				  num_phase_split.value(),true);
    } else {
      label = find_phase_labels(phasemaps[n],mask,num_phase_split.value());
    }
    if (labelvol.set()) { save_volume(label,labelvol.value()); }

    // Unwrap phase
    volume<float> uph;
    if (unwrap_2D) {
      uph = unwrap2D(phasemaps[n],label,verbose.value());
    } else {
      uph = unwrap(phasemaps[n],label,verbose.value());
    }
 
    // Restore phase ramps
    if (removeramps.value()) {
      restore_linear_ramps(ramps,uph,mask);
      restore_linear_ramps(ramps,phasemaps[n],mask);
      if (verbose.value()) { cout << "Removed linear ramps" << endl; }
    }

    uphase.addvolume(uph);
  }

  // Save outputs
  uphase.setTR(phasemaps.TR());
  save_volume4D(uphase,uphasevol.value());

  if (rawphasevol.set()) {
    save_volume4D(phasemaps,rawphasevol.value());
  }

  return 0;
}




int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    options.add(complexvol);
    options.add(absvol);
    options.add(phasevol);
    options.add(uphasevol);
    options.add(num_phase_split);
    options.add(labelslices);
    options.add(allslices);
    options.add(force3D);
    options.add(thresh);
    options.add(inmask);
    options.add(startimno);
    options.add(endimno);
    options.add(maskvol);
    options.add(rawphasevol);
    options.add(labelvol);
    options.add(removeramps);
    options.add(verbose);
    options.add(help);
    
    options.parse_command_line(argc, argv);

    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
    
    if ( ( (complexvol.set()) && (absvol.set()) ) || 
	 ( (complexvol.set()) && (phasevol.set()) ) ) 
      {
	options.usage();
	cerr << endl 
	     << "Cannot specify both --complex AND --phase or --abs."
	     << endl;
	exit(EXIT_FAILURE);
      }
 
    if ( ( (phasevol.set()) && (absvol.unset()) ) || 
	 ( (phasevol.unset()) && (absvol.set()) ) ) 
      {
	options.usage();
	cerr << endl 
	     << "Both --phase AND --abs must be used together." 
	     << endl;
	exit(EXIT_FAILURE);
      }
    
    if (num_phase_split.value() < 2)
      {
	options.usage();
	cerr << endl 
	     << "Always set --numphasesplit greater than 1." 
	     << endl << "NOT " << num_phase_split.value() << endl;
	exit(EXIT_FAILURE);
      }
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  return do_unwrapping();
}






