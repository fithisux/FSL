/*  FilmGlsOptionsRes.cc

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
#include "FilmGlsOptionsRes.h"

#ifndef NO_NAMESPACE
using namespace NEWMAT;
namespace FILM {
#endif

FilmGlsOptionsRes* FilmGlsOptionsRes::gopt = NULL;

void FilmGlsOptionsRes::parse_command_line(int argc, char** argv, ofstream& logfile)
{

  if(argc<2){
    print_usage(argc,argv);
    exit(1);
  }
  
  int inp = 1;
  int n=1;
  string arg;
  char first;
  
  logfile << "Arguments: " << endl;
  while (n<argc) {
    arg=argv[n];

    // output argument to logfile:
    logfile << arg << " ";

    if (arg.size()<1) { n++; continue; }
    first = arg[0];
    if (first!='-') {
      if(inp == 1)
	gopt->inputfname = arg;
      else if(inp == 2)
	gopt->paradigmfname = arg;
      else if(inp == 3)
	gopt->thresh = atoi(argv[n]);
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
    } else if ( arg == "-vdms" ) {
      gopt->verbose_dms = true;
      n++;
      continue;
    } else if ( arg == "-ar" ) {
      gopt->fitAutoRegressiveModel = true;
      gopt->tukey = false;
      n++;
      continue;
    } else if ( arg == "-pava" ) {
      gopt->pava = true;
      gopt->tukey = false;
      n++;
      continue;
    } else if ( arg == "-hfr" ) {
      gopt->highfreqremoval=true;
      n++;
      continue;
    } else if ( arg == "-sa" ) {
      gopt->smoothACEst = true;
      n++;
      continue;
    } else if ( arg == "-det" ) {
      gopt->detrend = true;
      n++;
      continue;
    } else if ( arg == "-noest" ) {
      gopt->noest = true;
      n++;
      continue;
    } else if ( arg == "-justpw" ) {
      gopt->justprewhiten = true;
      n++;
      continue;
    } 



    if (n+1>=argc) 
      { 
	cerr << "Lacking argument to option " << arg << endl;
	break; 
      }

    // put options with 1 argument here
    if ( arg == "-ms") {
      gopt->ms = atoi(argv[n+1]);
      logfile << argv[n+1] << " ";
      n+=2;
    } else if ( arg == "-iters") {
      gopt->numiters = atoi(argv[n+1]);
      logfile << argv[n+1] << " ";
      n+=2;
    } else if ( arg == "-rn" ) {
      gopt->datadir = argv[n+1];      
      logfile << argv[n+1] << " ";
      n+=2;
    } else if ( arg == "-tukey" ) {
      gopt->tukeysize = atoi(argv[n+1]);
      gopt->tukey = true;
      logfile << argv[n+1] << " ";
      n+=2;
    } else if ( arg == "-mt" ) {
      gopt->multitapersize = atoi(argv[n+1]);
      gopt->tukey = false;
      gopt->multitaper = true;
      logfile << argv[n+1] << " ";
      n+=2;
    }
    else if ( arg == "-epith" ) {
      gopt->epith = atoi(argv[n+1]);      
      logfile << argv[n+1] << " ";
      n+=2;
    } else { 
      cerr << "Unrecognised option " << arg << endl;
      n++;
    } 
  }  // while (n<argc)

  logfile << endl << "---------------------------------------------" << endl << endl;

  if (gopt->inputfname.size()<1) {
    print_usage(argc,argv);
    throw Exception("Input filename not found");
  }
  if (gopt->paradigmfname.size()<1) {
    print_usage(argc,argv);
    throw Exception("Paradigm filename not found");
  }
  if (gopt->thresh < 0) {
    print_usage(argc,argv);
    throw Exception("Invalid thresh");
  }
}

void FilmGlsOptionsRes::print_usage(int argc, char *argv[])
{
  cout << "Usage: " << argv[0] << " [options] <groupfile> <paradigmfile> <thresh>\n\n"
       << "  Available options are:\n"
       << "        -zout <zscoresfile>                (default is "
       << gopt->zscoresfname << ")\n"
       << "        -sa                                (smooths auto corr estimates)\n"
       << "        -ms <num>                          (susan mask size)\n"
       << "        -iters <num>                       (number of prewhitening iterations to use - default is "
       << gopt->numiters << ")\n"
       << "        -hfr                               (perform high frequency artefact removal)\n"
       << "        -epith <num>                       (set susan brightness threshold - otherwise it is estimated)\n"
       << "        -v                                 (outputs full data)\n"
       << "        -vdms                              (outputs prewhiten design matrices)\n"
       << "        -ar                                (fits autoregressive model - default " 
       << "is to use tukey with M=2*sqrt(numvols))\n"
       << "        -tukey <num>                       (uses tukey window to estimate autocorr with window size num - default "
       << "is to use tukey with M=2*sqrt(numvols))\n"
       << "        -mt <num>                          (uses multitapering with slepian tapers and num is the time-bandwidth product - default " 
       << "is to use tukey with M=2*sqrt(numvols))\n"
       << "        -pava                              (estimates autocorr using PAVA - default " 
       << "is to use tukey with M=2*sqrt(numvols))\n"
       << "        -noest                             (do not estimate auto corrs)\n"
       << "        -det                               (linear detrend data)\n"
       << "        -rn <dir>                          (directory name to store results in, default is "
       << gopt->datadir << ")\n"
       << "        -justpw                            (just prewhiten and output prewhitened data)\n"
       << "        -help\n\n";
}

#ifndef NO_NAMESPACE
}
#endif
