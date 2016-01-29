/*  readoptions.h

    Michael Chappell - FMRIB Image Analysis Group

    Copyright (C) 2009 University of Oxford  */

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

#if !defined(ReadOptions_h)
#define ReadOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;

namespace OXASL {

class ReadOptions {
 public:
  static ReadOptions& getInstance();
  ~ReadOptions() { delete ropt; }
  
  Option<bool> help;

  Option<string> datafile;
  Option<string> maskfile;
  Option<int> ntis;
  
  Option<string> inblockform;
  Option<string> inaslform;
  Option<bool> ispairs;

  Option<bool> splitpairs;
  Option<bool> tcdiff;
  Option<bool> surrtcdiff;
  
  Option<string> outblockform;
  Option<string> out;

  Option<string> meanout;
  Option<string> splitout;
  
  Option<string> epochout;
  Option<int> epochlen;
  Option<int> epochover;
  Option<string> epochunit;

  Option<string> deconvout;
  Option<string> aif;

  void parse_command_line(int argc, char** argv);
  
 private:
  ReadOptions();  
  const ReadOptions& operator=(ReadOptions&);
  ReadOptions(ReadOptions&);

  OptionParser options; 
      
  static ReadOptions* ropt;
  
};

 inline ReadOptions& ReadOptions::getInstance(){
   if(ropt == NULL)
     ropt = new ReadOptions();
   
   return *ropt;
 }

 inline ReadOptions::ReadOptions() :

help(string("-h,--help"), false,
		    string("display this message"),
		    false, no_argument),
   //input files
   datafile(string("--data,--datafile"), string("ASL datafile"),
			  string("data file"),
		     true, requires_argument),  
   maskfile(string("--mask"), string("maskfile\n"),
	    string("mask"),
	    false, requires_argument),
   //input file information
   ntis(string("--ntis"),0,
	string("Number of TIs in file"),
	true, requires_argument),
   inblockform(string("--ibf,--inblockform"),string("rpt"),
	       string("Input block format:\n          rpt - blocks of measurements that include all TIs\n          tis - blocks of repeated measurements at a single TI"),
	       false,requires_argument),
   inaslform(string("--iaf,--inaslform"),string("diff"),
	     string("ASL data form:\n          diff - differenced data {default}\n          tc - Tag-Control pairs\n          ct - Control-Tag pairs\n"),
	     false,requires_argument),
   ispairs(string("--pairs,--inpairs"),false,
	   string("Data contains adjacent pairs of measuremnts (e.g. Tag, Control)"),
	   false,no_argument),
   

   //asaq(string("--asaq"),false,
   //	string("Data is as aquired: same as --blocked --pairs"),
   //	false,no_argument),

   // manipulation options
   splitpairs(string("--spairs"),false,
	      string("Split the pairs within the data, e.g. to separate tag and control images in output"),
	      false,no_argument),
   tcdiff(string("--diff"), false,
	   string("Take the difference between the pairs, i.e. Tag control difference\n"),
	   false,no_argument),
   surrtcdiff(string("--surrdiff"), false,
	      string("Do surround subtraction on the pairs"),
	      false,no_argument),

   //basic output
   outblockform(string("--obf,--outblockform"),string("notset"),
	       string("Output block format (for --out=):\n          rpt - blocks of measurements that include all TIs\n          tis - blocks of repeated measurements at a single TI\n          Default is same as input block format (--ibf)"),
		  false,requires_argument),
   out(string("--out"),string("Out filename"),
       string("Output data file"),
       false,requires_argument),

   // other output options
   meanout(string("--mean"),string(""),
	  string("Output ASL data having taken mean at each TI to file"),
	   false, requires_argument),
   splitout(string("--split"),string(""),
	    string("Split data into separate files each each TI, specify filename root\n"),
	    false, requires_argument),

   epochout(string("--epoch"),string(""),
	    string("Output epochs of ASL data (takes mean at each TI within the epoch)"),
	    false, requires_argument),
   epochlen(string("--elen,--epochlen"),1,
	    string("Length of epochs in number of repeats"),
	    false, requires_argument),
   epochover(string("--eol,--epochol"),0,
	     string("Ammount of overlap between epochs in number of repeats"),
	     false, requires_argument),
   epochunit(string("--eunit,--epochunit"),string("rpt"),
	     string("Epochs to be determined over:\n          rpt - repeats in the data {default}\n          tis - TIs in the data\n"),
	     false,requires_argument),

   deconvout(string("--deconv"),string(""),
	     string("Deconvolution of data with arterial input functions"),
	     false,requires_argument),
   aif(string("--aif"),string(""),
       string("Arterial input functions for deconvolution (4D volume, one aif for each voxel within mask)"),
       false,requires_argument),

  
   options("asl_file","asl_file --data=<asldata> --ibf=rpt --iaf=tc --diff --out=<diffdata>\n")
   {
     try {
       options.add(help);

       options.add(datafile);
       options.add(maskfile);

       options.add(ntis);
       options.add(inblockform);
       options.add(inaslform);
       options.add(ispairs);
 
       options.add(splitpairs);
       options.add(tcdiff);
       options.add(surrtcdiff);

       options.add(out);
       options.add(outblockform);

       options.add(meanout);
       options.add(splitout);

       options.add(epochout);
       options.add(epochlen);
       options.add(epochover);
       options.add(epochunit);

       options.add(deconvout);
       options.add(aif);

        }
     catch(X_OptionError& e) {
       options.usage();
       cerr << endl << e.what() << endl;
     } 
     catch(std::exception &e) {
       cerr << e.what() << endl;
     }    
     
   }

}
#endif



