/*  ContrastMgrOptions.cc

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

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

#define WANT_STREAM
#define WANT_MATH

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "newmatap.h"
#include "newmatio.h"
#include "ContrastMgrOptions.h"

using namespace NEWMAT;
namespace FILM {

ContrastMgrOptions* ContrastMgrOptions::gopt = NULL;

void ContrastMgrOptions::parse_command_line(int argc, char** argv, Log& logger)
{
  if(argc<2){
    print_usage(argc,argv);
    exit(1);
  }
  
  int inp = 1;
  int n=1;
  string arg;
  char first;
  
  while (n<argc) {
    arg=argv[n];

    if (arg.size()<1) { n++; continue; }
    first = arg[0];
    if (first!='-') {
      if(inp == 1)
	gopt->datadir = arg;
      else if(inp == 2)
	gopt->contrastfname = arg;      
      else
	{
	  cerr << "Mismatched argument " << arg << endl;
	  break;
	}
      n++;
      inp++;
      continue;
    }
    
    // put options without arguments here
    if ( arg == "-help" ) {
      print_usage(argc,argv);
      exit(0);
    } else if ( arg == "-v" ) {
      gopt->verbose = true;
      n++;
      continue;
    } else if ( arg == "-d" ) {
      gopt->debug = true;
      n++;
      continue;
    }

    if ( arg == "-cope") {
      gopt->copenumber = atoi(argv[n+1]);
      n+=2;
      continue;
    } 
    else if ( arg == "-suffix") {
      gopt->suffix = argv[n+1];
      n+=2;
      continue;
    } 
    else if (arg == "-f") {
      gopt->fcontrastfname = argv[n+1];
      n+=2;
      continue;
    } 

    if (n+1>=argc) 
      { 
	cerr << "Lacking argument to option " << arg << endl;
	break; 
      }
  }  // while (n<argc)

  logger.setDir(gopt->datadir,"cmlogfile");
  for(int a = 0; a < argc; a++)
    logger.str() << argv[a] << " ";
  logger.str() << endl << "---------------------------------------------" << endl << endl;

  if (gopt->datadir.size()<1) {
    print_usage(argc,argv);
    throw Exception("Data directory not found");
  }
  if (gopt->contrastfname.size()<1) {
    print_usage(argc,argv);
    throw Exception("Contrast filename not found");
  }
 
}

void ContrastMgrOptions::print_usage(int argc, char *argv[])
{

  cout << "Usage: " << argv[0] << " [options] <datadir> <tcontrastfile>\n\n"
       << "  Available options are:\n"
       << "        -v                                 (outputs full data)\n"
       << "        -cope <num>                        (contrast number to start labelling copes from, default is 1)\n"
       << "        -suffix <string>                   (suffix to put on the end of the cope filename before the contrast number, default is nothing)\n"
       << "        -f <fcontrastfile>                 (fcontrast file)\n"
       << "        -d                                 (debug)\n"
       << "        -help\n\n";
}

#ifndef NO_NAMESPACE
}
#endif
