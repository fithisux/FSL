/*  slicetimer.cc - FMRIB's Slice Timing Utility
    
    Peter Bannister, Stephen Smith and Matthew Webster, FMRIB Image Analysis Group
    
    Copyright (C) 2001-2009 University of Oxford  */

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


#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>

#include "miscmaths/optimise.h"
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "utils/options.h"
#include "miscmaths/kernel.h"

using namespace MISCMATHS;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace Utilities;

string title="slicetimer \nFMRIB's Interpolation for Slice Timing\nCopyright(c) 2001-2009, University of Oxford";
string examples="slicetimer -i <timeseries> [-o <corrected_timeseries>] [options]\n";

Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> odd(string("--odd"), false, 
		     string("use interleaved acquisition"), 
		     false, no_argument);
Option<bool> down(string("--down"), false, 
		     string("reverse slice indexing (default is: slices were acquired bottom-up)"), 
		     false, no_argument);
Option<string> inputname(string("-i,--in"), string(""),
			 string("filename of input timeseries"),
			 true, requires_argument);
Option<string> outputname(string("-o,--out"), string(""),
			  string("filename of output timeseries"),
			  false, requires_argument);
Option<string> tcustom(string("--tcustom"), string(""),
			  string("filename of single-column slice timings, in fractions of TR, +ve values shift slices forwards in time."),
			  false, requires_argument);
Option<float>  tglobal(string("--tglobal"), 0,
			  string("global shift in fraction of TR, (default is 0)"),
			  false, requires_argument);
Option<string> ocustom(string("--ocustom"), string(""),
			  string("filename of single-column custom interleave order file (first slice is referred to as 1 not 0)"),
			  false, requires_argument);
Option<int> direction(string("-d,--direction"), 3,
		      string("direction of slice acquisition (x=1,y=2,z=3) - default is z"),
		      false, requires_argument);
Option<float> repeat(string("-r,--repeat"), 3.0,
		      string("Specify TR of data - default is 3s"),
		      false, requires_argument);

int do_slice_correction()
{
  volume4D<float> timeseries;
  Matrix timings;
  
  int no_volumes, no_slices;
  float repeat_time, offset=0, slice_spacing;

  if (inputname.set()) {
    if (verbose.value()) { cout << "Reading input volume" << endl; }
    read_volume4D(timeseries,inputname.value());
    if (!outputname.set())
      outputname.set_value(inputname.value() + "_st");
  } else if (outputname.set()) {
    cerr << "Must specify an input volume (-i or --in) to generate corrected data." 
	 << endl;
    return -1;
  }

  no_volumes = timeseries.tsize();
  repeat_time = timeseries.tdim();
  if (repeat_time ==0){
    // cerr << "Zero TR in file header - fixing ... " ;
    repeat_time = repeat.value();
  }

  if (direction.value() == 1) timeseries.swapdimensions(3,2,1); // Flip z and x
  if (direction.value() == 2) timeseries.swapdimensions(1,3,2); // Flip z and y
  no_slices = timeseries.zsize();
  slice_spacing = repeat_time / no_slices;

  if (tcustom.set()) {
    timings.ReSize(1, timeseries.zsize());
    timings = read_ascii_matrix(tcustom.value(), timeseries.zsize(), 1);
  }  else if (ocustom.set()) {
    timings.ReSize(1, timeseries.zsize());
    timings = read_ascii_matrix(ocustom.value(), timeseries.zsize(), 1);
  }

  ColumnVector userkernel = sinckernel1D("hanning", 7, 1201);
  // for(int i=1; i<=1201; i++) cout << i << " " << userkernel(i) << endl;

  float recenter = (((float) no_slices)/2 - 0.5)/ no_slices; // only valid for slice-count-based corrections
  
  for (int slice=1; slice<=no_slices; slice++) {

    if (tglobal.set()) {
      offset = tglobal.value();
    } else if (tcustom.set()) {
      offset = timings(slice, 1);
    } else if (ocustom.set()) { 
      int local_count=1;
      while (local_count <= no_slices) {
	if (timings(local_count, 1) == slice) {
	  offset = recenter -(local_count -1)* (slice_spacing / repeat_time);
	  local_count = no_slices + 1;
	} else 
	  local_count++;
      }
    } else if (odd.value()) { // acquisition order: 1,3,5, ..., 2,4,6 ...
      if ((slice % 2) == 0) // even
	offset = recenter - ( ceil((float)no_slices / 2) + ((slice -1)/ 2)) * (slice_spacing / repeat_time);
      else
	offset = recenter - ((slice -1) / 2) * (slice_spacing / repeat_time);
    } else if (down.value()) {
      offset = recenter - (no_slices - slice) * (slice_spacing / repeat_time);
    } else { 
      offset = recenter - (slice -1) * (slice_spacing / repeat_time);
    }

    for (int x_pos = 0; x_pos < timeseries. xsize(); x_pos++)
      for (int y_pos = 0; y_pos < timeseries. ysize(); y_pos++){
	ColumnVector voxeltimeseries = timeseries.voxelts(x_pos,y_pos,slice-1);
	ColumnVector interpseries = voxeltimeseries;
	for (int time_step=1; time_step <= no_volumes; time_step++){
	  // interpseries(time_step) = interpolate_1d(voxeltimeseries, time_step - offset);
	  interpseries(time_step) = kernelinterpolation_1d(voxeltimeseries, time_step + offset, userkernel, 7);
	}
	timeseries.setvoxelts(interpseries,x_pos,y_pos,slice-1);
      }
    
    if (verbose.value())
      cerr << "Slice " << slice << " offset " << offset << endl;
  }
  
  if (direction.value() == 1) timeseries.swapdimensions(3,2,1); // reverse Flip z and x
  if (direction.value() == 2) timeseries.swapdimensions(1,3,2); // reverse Flip z and y

  save_volume4D(timeseries, outputname.value());
  return 0;
}

int main (int argc,char** argv)
{
  Tracer tr("main");
  
  OptionParser options(title, examples);

  try {
    options.add(inputname);
    options.add(outputname);
    options.add(help);
    options.add(verbose);
    options.add(down);
    options.add(repeat);
    options.add(direction);
    options.add(odd);
    options.add(tcustom);
    options.add(tglobal);
    options.add(ocustom);

    options.parse_command_line(argc, argv);

    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
        
    if ( inputname.unset()) 
      {
	options.usage();
	cerr << endl 
	     << "--in or -i MUST be used." 
	     << endl;
	exit(EXIT_FAILURE);
      }
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  int retval = do_slice_correction();
  
  if (retval!=0) {
    cerr << endl << endl << "Error detected: try -h for help" << endl;
  }

  return retval;
}

