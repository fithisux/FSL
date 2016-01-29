// Utility for resampling of DTI data using fields/data
// estimated by topup.
//
// applytopup.cpp
//
// Jesper Andersson, FMRIB Image Analysis Group
//
// Copyright (C) 2012 University of Oxford

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

#include <cstring>
#include "utils/options.h"
#include "newmat.h"
#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS           // To allow us to use .set_sform etc
#endif
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
#include "warpfns/warpfns.h"
#include "basisfield/basisfield.h"
#include "basisfield/splinefield.h"
#include "warpfns/fnirt_file_reader.h"
#include "topup_file_io.h"
#include "displacement_vector.h"
#include "applytopup.h"

using namespace Utilities;
using namespace TOPUP;

////////////////////////////////////////////////////////////////////////////

// COMMAND LINE OPTIONS

string title="applytopup (Version 1.0)\nCopyright(c) 2009, University of Oxford (Jesper Andersson)";
string examples=string("applytopup -in=topdn,botup --topup=mytu --inindex=1,2 --out=hifi\n") +
                string("applytopup -in=topdn --topup=mytu --inindex=1 --method=jac --interp=spline --out=hifi\n");

Utilities::Option<bool> verbose(string("-v,--verbose"), false, 
		                string("switch on diagnostic messages"), 
		                false, Utilities::no_argument);
Utilities::Option<bool> help(string("-h,--help"), false,
		             string("display this message"),
		             false, Utilities::no_argument);
Utilities::Option<string> interpstr(string("-n,--interp"), string("spline"),
		                    string("interpolation method {trilinear,spline}, default=spline"),
		                    false, Utilities::requires_argument);
Utilities::Option<string> infname(string("-i,--imain"), string(""),
		                  string("comma separated list of names of input image (to be corrected)"),
		                  true, Utilities::requires_argument);
Utilities::Option<string> datain(string("-a,--datain"), string(""),
		                 string("name of text file with PE directions/times"),
		                 true, Utilities::requires_argument);
Utilities::Option<string> inindex(string("-x,--inindex"), string(""),
		                  string("comma separated list of indicies into --datain of the input image (to be corrected)"),
		                  true, Utilities::requires_argument);
Utilities::Option<string> outfname(string("-o,--out"), string(""),
		                   string("basename for output (warped) image"),
		                   true, Utilities::requires_argument);
Utilities::Option<string> datatype(string("-d,--datatype"), string(""),
                                   string("Force output data type [char short int float double]."),
                                   false, Utilities::requires_argument);
Utilities::Option<string> method(string("-m,--method"), string("lsr"),
                                 string("Use jacobian modulation (jac) or least-squares resampling (lsr), default=lsr."),
                                 false, Utilities::requires_argument);
Utilities::Option<string> topupfname(string("-t,--topup"), string(""),
			             string("name of field/movements (from topup)"),
			             true, Utilities::requires_argument);

int main(int argc, char *argv[])
{

  Tracer tr("main");

  OptionParser options(title, examples);

  try {
    options.add(infname);
    options.add(datain);
    options.add(inindex);
    options.add(topupfname);
    options.add(outfname);
    options.add(method);
    options.add(interpstr);
    options.add(datatype);
    options.add(verbose);
    options.add(help);

    int i=options.parse_command_line(argc, argv);
    if (i < argc) {
      for (; i<argc; i++) {
        cerr << "Unknown input: " << argv[i] << endl;
      }
      exit(EXIT_FAILURE);
    }

    if (help.value() || !options.check_compulsory_arguments(true)) {
      options.usage();
      exit(EXIT_FAILURE);
    }
  }  
  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } 
  catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  return(applytopup());
}

namespace TOPUP {

int applytopup()
{
  // Check that input is kosher

  // Assert value for data-type
  short  dtypecode = DT_FLOAT;
  if (datatype.set()) {
    if (datatype.value()==string("char")) dtypecode = DT_UNSIGNED_CHAR;
    else if (datatype.value()==string("short")) dtypecode = DT_SIGNED_SHORT;
    else if (datatype.value()==string("int")) dtypecode = DT_SIGNED_INT;
    else if (datatype.value()==string("float")) dtypecode = DT_FLOAT;
    else if (datatype.value()==string("double")) dtypecode = DT_DOUBLE;
    else throw ApplyTopupException(string("Unknown data type ")+datatype.value());
  }
  NEWIMAGE::interpolation  interp = NEWIMAGE::spline;
  if (interpstr.set()) {
    if (interpstr.value()==string("trilinear")) interp = NEWIMAGE::trilinear;
    else if (interpstr.value()==string("spline")) interp = NEWIMAGE::spline;  
    else throw ApplyTopupException(string("Unknown interpolation type ")+interpstr.value());
  }

  // Parse list and read input files

  // Parse and check comma separated lists
  std::vector<std::string> infnames = parse_commaseparated_list(infname.value());
  std::vector<unsigned int> inindices = parse_commaseparated_numbers(inindex.value());
  if (infnames.size() != inindices.size()) throw ApplyTopupException("Mismatched --in and --inindex lists");

  // Read input image files
  std::vector<NEWIMAGE::volume4D<float> >  scans(infnames.size());
  for (unsigned int i=0; i<infnames.size(); i++) {
    NEWIMAGE::read_volume4D(scans[i],infnames[i]);
    if (i) {
      if (!NEWIMAGE::samesize(scans[i],scans[0])) throw ApplyTopupException(string("Mismatched input files ")+infnames[i]+string(" and ")+infnames[0]);
    }
  }

  // Read and parse file with phase encode vectors and readout times
  TopupDatafileReader  datafile(datain.value());
  for (unsigned int i=0; i<inindices.size(); i++) if (!inindices[i] || inindices[i]>datafile.N()) throw ApplyTopupException("Invalid entry in --inindex list");

  // Read and parse field and movements estimated by topup
  TopupFileReader  topupfile(topupfname.value());

  NEWIMAGE::volume4D<float>  ovol;
  NEWIMAGE::volume<char>   mask(scans[0][0].xsize(),scans[0][0].ysize(),scans[0][0].zsize());
  copybasicproperties(scans[0][0],mask);
  mask = 1;

  if (method.value() == "jac") { // Use simple resampling and jacobian modulation
    ovol = jac_resample(datafile,topupfile,inindices,mask,scans,interp);
  }
  else if (method.value() == "lsr") { // Use least squares restoration
    // First resample everything using the movement parameters
    resample_using_movement_parameters(datafile,topupfile,inindices,mask,scans,interp);

    // Then do LSR resampling
    ovol = lsr_resample(datafile,topupfile,inindices,mask,scans);
  }
  else if (method.value() == "vb2D") { // Use VB based least squares restoration
    // First resample everything using the movement parameters
    resample_using_movement_parameters(datafile,topupfile,inindices,mask,scans,interp);

    // Then do LSR resamling using VB to infer on the regularisation
    ovol = vb_resample_2D(datafile,topupfile,inindices,mask,scans);
  }
  else if (method.value() == "vb3D") { // Use VB based least squares restoration
    // First resample everything using the movement parameters
    resample_using_movement_parameters(datafile,topupfile,inindices,mask,scans,interp);

    // Then do LSR resamling using VB to infer on the regularisation
    ovol = vb_resample_3D(datafile,topupfile,inindices,mask,scans);
  }
  else if (method.value() == "vb4D") { // Use VB based least squares restoration
    // First resample everything using the movement parameters
    resample_using_movement_parameters(datafile,topupfile,inindices,mask,scans,interp);

    // Then do LSR resamling using VB to infer on the regularisation
    ovol = vb_resample_4D(datafile,topupfile,inindices,mask,scans);
  }
  else {
    // Dodgy!!
  }
  // Mask out other dodgy bits
  for (int t=0; t<ovol.tsize(); t++) {
    for (int k=0; k<ovol.zsize(); k++) {
      for (int j=0; j<ovol.ysize(); j++) {
        for (int i=0; i<ovol.xsize(); i++) {
          ovol(i,j,k,t) = (mask(i,j,k)) ? ovol(i,j,k,t) : 0.0;
        }
      }
    }
  }

  // Save the results
  ovol.setDisplayMaximumMinimum(0,0);
  if (datatype.set()) {
    save_volume4D_dtype(ovol,outfname.value(),dtypecode);
  }
  else {
    save_volume4D_dtype(ovol,outfname.value(),NEWIMAGE::dtype(infnames[0]));
  }
    
  return(EXIT_SUCCESS);
}

NEWIMAGE::volume4D<float> vb_resample_4D(const TopupDatafileReader&                       datafile, 
                                         const TopupFileReader&                           topupfile, 
                                         const std::vector<unsigned int>&                 inindices, 
                                         NEWIMAGE::volume<char>&                          mask,
                                         std::vector<NEWIMAGE::volume4D<float> >&         scans)
{
  NEWIMAGE::volume<float>   field = topupfile.FieldAsVolume();
  TopupCollections          collections(inindices,datafile);
  NEWIMAGE::volume4D<float> ovol = scans[0];
  ovol = 0.0;

  cout << "NCollections = " << collections.NCollections() << endl;
  cout << "NScans = " << collections.NScans(0) << endl;

  for (unsigned int c=0; c<collections.NCollections(); c++) {
    // Find PE direction for this collection
    bool row = false;
    unsigned int pesz = scans[0].ysize();
    unsigned int fesz = scans[0].xsize();
    if (datafile.PhaseEncodeVector(collections.IndexAt(c,0))(1)) { // If PE in x-dir for this collection
      row = true;
      pesz = scans[0].xsize();
      fesz = scans[0].ysize();
    }
    int m=pesz;
    int n=2*m;
    int N=fesz*scans[0].zsize()*scans[0].tsize();
    std::vector<std::vector<MISCMATHS::SpMat<double> > >           K(scans[0].zsize());    // ALL K-matrices
    std::vector<std::vector<MISCMATHS::SpMat<double> > >           KtK(scans[0].zsize());  // ALL KtK-matrices
    std::vector<std::vector<std::vector<NEWMAT::ColumnVector> > >  y(scans[0].tsize());    // All data for all diffusion gradients
    std::vector<std::vector<std::vector<double> > >                yty(scans[0].tsize());  // All yty for all diffusion gradients
    std::vector<double> sf(2,0.0);
    if (row) { sf[0] = scans[0].xdim()/scans[0].ydim(); sf[1] = scans[0].xdim()/scans[0].zdim(); }
    else { sf[0] = scans[0].ydim()/scans[0].xdim(); sf[1] = scans[0].ydim()/scans[0].zdim(); }   
    std::vector<std::vector<SpMat<double> > >  Lijk = GetLijk(pesz,sf);               // Regularisation sub-matrices
    MISCMATHS::SpMat<double>                   Lii = GetLii(pesz);                    // For initialization of mu
    // Pre-compute K- and KtK-matrices
    for (int sl=0; sl<scans[0].zsize(); sl++) {
      K[sl].resize(fesz);
      KtK[sl].resize(fesz);
      for (unsigned int i=0; i<fesz; i++) {
	DispVec         dv(pesz);
	if (row) dv.SetFromRow(field,sl,i);
	else dv.SetFromColumn(field,sl,i);
        for (unsigned int a=0; a<collections.NScans(c); a++) { // Loop over acquisitions in collection
  	  double psf = (datafile.PhaseEncodeVector(collections.IndexAt(c,a)))(2) * datafile.ReadOutTime(collections.IndexAt(c,a));
          if (row) psf = (datafile.PhaseEncodeVector(collections.IndexAt(c,a)))(1) * datafile.ReadOutTime(collections.IndexAt(c,a));
          if (!a) K[sl][i] = dv.GetSparseK_Matrix(psf);
          else K[sl][i] &= dv.GetSparseK_Matrix(psf);
	}
        KtK[sl][i] = K[sl][i].TransMultSelf();
      }
    }
    // Loop over all directions and gather the data
    for (int d=0; d<ovol.tsize(); d++) {
      y[d].resize(scans[0].zsize());
      yty[d].resize(scans[0].zsize());
      for (int sl=0; sl<scans[0].zsize(); sl++) {
        y[d][sl].resize(fesz);
        yty[d][sl].resize(fesz);
        for (unsigned int i=0; i<fesz; i++) { 
          for (unsigned int a=0; a<collections.NScans(c); a++) {
            if (!a) y[d][sl][i] = extract_row_col((scans[collections.ScanAt(c,a)])[d],sl,i,row);
            else y[d][sl][i] &= extract_row_col((scans[collections.ScanAt(c,a)])[d],sl,i,row);
          }
          yty[d][sl][i] = (y[d][sl][i].t()*y[d][sl][i]).AsScalar();
	}
      }
    }

    // Set priors
    double l_0 = 1e-10; double tau_0 = 1e10;
    double k_0 = 1e-10; double theta_0 = 1e10;
    // Initialize all parameters to their priors
    double l = l_0; double tau = tau_0;
    double k = k_0; double theta = theta_0;
    std::vector<std::vector<SpMat<double> > >                      Lambda(scans[0].zsize());
    std::vector<std::vector<std::vector<NEWMAT::ColumnVector> > >  mu(scans[0].tsize());
    for (int d=0; d<scans[0].tsize(); d++) {
      mu[d].resize(scans[0].zsize());
      for (int sl=0; sl<scans[0].zsize(); sl++) {
        Lambda[sl].resize(fesz);
        mu[d][sl].resize(fesz);
      }
    }
    for (int sl=0; sl<scans[0].zsize(); sl++) {
      for (unsigned int i=0; i<fesz; i++) {
        Lambda[sl][i] = k*theta*KtK[sl][i] + l*tau*Lii;
	NEWMAT::CroutMatrix tmpLambda = Lambda[sl][i].AsNEWMAT();
        for (int d=0; d<scans[0].tsize(); d++) {
          mu[d][sl][i] = (k*theta) * (tmpLambda.i()*K[sl][i].TransMult(y[d][sl][i]));
        }
      }
    }

    // Do the VB estimation
    for (int iter=0; iter<20; iter++) {

      cout << "Relative lambda = " << (l*tau)/(k*theta) << ",   phi = " << k*theta << ",   lambda = " << l*tau << endl;
      cout << "Updating phi" << endl;
      // Update parameters for phi (noise precision)
      k = k_0 + n*N/2;
      theta = 1.0 / theta_0;
      for (int sl=0; sl<scans[0].zsize(); sl++) {
        for (unsigned int i=0; i<fesz; i++) {
          theta += 0.5 * double(scans[0].tsize()) * (Lambda[sl][i].AsNEWMAT().i()*KtK[sl][i].AsNEWMAT()).Trace();
          for (int d=0; d<scans[0].tsize(); d++) {
            theta += (0.5 * (mu[d][sl][i].t()*(KtK[sl][i]*mu[d][sl][i]) + yty[d][sl][i])).AsScalar();
            theta -= (mu[d][sl][i].t()*(K[sl][i].TransMult(y[d][sl][i]))).AsScalar();
          }
        }
      }
      theta = 1.0/theta;

      cout << "Updating lambda" << endl;
      // Update parameters for lambda (regularisation weight)
      l = l_0 + m*N/2;
      tau = tau_update(mu,Lambda,Lijk,tau_0);
      /*
      tau = 1.0 / tau_0;
      for (int d=0; d<scans[0].tsize(); d++) {
        tau += tau_update_contrib(mu[d],Lambda,Lijk);
      }
      tau = 1.0 / tau;
      */

      cout << "Updating images" << endl;
      // Update image and variance of image
      for (int sl=0; sl<scans[0].zsize(); sl++) {
        for (unsigned int i=0; i<fesz; i++) {
          Lambda[sl][i] = k*theta*KtK[sl][i] + l*tau*Lijk[0][0];
	  NEWMAT::CroutMatrix tmpLambda = Lambda[sl][i].AsNEWMAT();
          for (int d=0; d<scans[0].tsize(); d++) {
            mu[d][sl][i] = mu_update(mu[d],K[sl][i],y[d][sl][i],tmpLambda,Lijk,k*theta,l*tau,sl,i);
          }
        }    
      }
    }
    // Add to output
    for (int d=0; d<scans[0].tsize(); d++) {
      for (int sl=0; sl<scans[0].zsize(); sl++) {
        add_to_slice(ovol[d],mu[d][sl],sl,row);
      }
    }
      
  } // End collections

  return(ovol);
}

NEWIMAGE::volume4D<float> vb_resample_3D(const TopupDatafileReader&                       datafile, 
                                         const TopupFileReader&                           topupfile, 
                                         const std::vector<unsigned int>&                 inindices, 
                                         NEWIMAGE::volume<char>&                          mask,
                                         std::vector<NEWIMAGE::volume4D<float> >&         scans)
{
  NEWIMAGE::volume<float>   field = topupfile.FieldAsVolume();
  TopupCollections          collections(inindices,datafile);
  NEWIMAGE::volume4D<float> ovol = scans[0];
  ovol = 0.0;

  cout << "NCollections = " << collections.NCollections() << endl;
  cout << "NScans = " << collections.NScans(0) << endl;

  for (unsigned int c=0; c<collections.NCollections(); c++) {
    // Find PE direction for this collection
    bool row = false;
    unsigned int pesz = scans[0].ysize();
    unsigned int fesz = scans[0].xsize();
    if (datafile.PhaseEncodeVector(collections.IndexAt(c,0))(1)) { // If PE in x-dir for this collection
      row = true;
      pesz = scans[0].xsize();
      fesz = scans[0].ysize();
    }
    int m=pesz;
    int n=2*m;
    int N=fesz*scans[0][0].zsize();
    std::vector<std::vector<MISCMATHS::SpMat<double> > >   K(scans[0].zsize());    // ALL K-matrices
    std::vector<std::vector<MISCMATHS::SpMat<double> > >   KtK(scans[0].zsize());  // ALL KtK-matrices
    std::vector<std::vector<NEWMAT::ColumnVector> >        y(scans[0].zsize());    // All data for one diffusion gradient
    std::vector<std::vector<double> >                      yty(scans[0].zsize());  // ALl yty for one diffusion gradient
    std::vector<double> sf(2,0.0);
    if (row) { sf[0] = scans[0].xdim()/scans[0].ydim(); sf[1] = scans[0].xdim()/scans[0].zdim(); }
    else { sf[0] = scans[0].ydim()/scans[0].xdim(); sf[1] = scans[0].ydim()/scans[0].zdim(); }   
    std::vector<std::vector<SpMat<double> > >  Lijk = GetLijk(pesz,sf);               // Regularisation sub-matrices
    MISCMATHS::SpMat<double>                   Lii = GetLii(pesz);                    // For initialization of mu
    // Pre-compute K- and KtK-matrices
    for (int sl=0; sl<scans[0].zsize(); sl++) {
      K[sl].resize(fesz);
      KtK[sl].resize(fesz);
      y[sl].resize(fesz);
      yty[sl].resize(fesz);
      for (unsigned int i=0; i<fesz; i++) {
	DispVec         dv(pesz);
	if (row) dv.SetFromRow(field,sl,i);
	else dv.SetFromColumn(field,sl,i);
        for (unsigned int a=0; a<collections.NScans(c); a++) { // Loop over acquisitions in collection
  	  double psf = (datafile.PhaseEncodeVector(collections.IndexAt(c,a)))(2) * datafile.ReadOutTime(collections.IndexAt(c,a));
          if (row) psf = (datafile.PhaseEncodeVector(collections.IndexAt(c,a)))(1) * datafile.ReadOutTime(collections.IndexAt(c,a));
          if (!a) K[sl][i] = dv.GetSparseK_Matrix(psf);
          else K[sl][i] &= dv.GetSparseK_Matrix(psf);
	}
        KtK[sl][i] = K[sl][i].TransMultSelf();
      }
    }
    // Loop over all directions
    for (int d=0; d<ovol.tsize(); d++) {

      // Gather all the data for this direction
      for (int sl=0; sl<scans[0].zsize(); sl++) {
        for (unsigned int i=0; i<fesz; i++) { 
          for (unsigned int a=0; a<collections.NScans(c); a++) {
            if (!a) y[sl][i] = extract_row_col((scans[collections.ScanAt(c,a)])[d],sl,i,row);
            else y[sl][i] &= extract_row_col((scans[collections.ScanAt(c,a)])[d],sl,i,row);
          }
          yty[sl][i] = (y[sl][i].t()*y[sl][i]).AsScalar();
	}
      }

      // Set priors
      double l_0 = 1e-10; double tau_0 = 1e10;
      double k_0 = 1e-10; double theta_0 = 1e10;
      // Initialize all parameters to their priors
      double l = l_0; double tau = tau_0;
      double k = k_0; double theta = theta_0;
      std::vector<std::vector<SpMat<double> > >        Lambda(scans[0].zsize());
      std::vector<std::vector<NEWMAT::ColumnVector> >  mu(scans[0].zsize());
      for (int sl=0; sl<scans[0].zsize(); sl++) {
        Lambda[sl].resize(fesz);
        mu[sl].resize(fesz);
        for (unsigned int i=0; i<fesz; i++) {
          Lambda[sl][i] = k*theta*KtK[sl][i] + l*tau*Lii;
          mu[sl][i] = (k*theta)*Lambda[sl][i].AsNEWMAT().i()*K[sl][i].TransMult(y[sl][i]);
        }
      }

      // Do the VB estimation
      for (int iter=0; iter<20; iter++) {

        cout << "Relative lambda = " << (l*tau)/(k*theta) << ",   phi = " << k*theta << ",   lambda = " << l*tau << endl;
        cout << "Updating phi" << endl;
        // Update parameters for phi (noise precision)
        k = k_0 + n*N/2;
        theta = 1.0 / theta_0;
        for (int sl=0; sl<scans[0][0].zsize(); sl++) {
          for (unsigned int i=0; i<fesz; i++) {
            theta += (0.5 * (mu[sl][i].t()*(KtK[sl][i]*mu[sl][i]) + yty[sl][i] + (Lambda[sl][i].AsNEWMAT().i()*KtK[sl][i].AsNEWMAT()).Trace())).AsScalar();
            theta -= (mu[sl][i].t()*(K[sl][i].TransMult(y[sl][i]))).AsScalar();
          }
        }
        theta = 1.0/theta;

        cout << "Updating lambda" << endl;
        // Update parameters for lambda (regularisation weight)
        l = l_0 + m*N/2;
        tau = tau_update(mu,Lambda,Lijk,tau_0);

        cout << "Updating image" << endl;
        // Update image and variance of image
        for (int sl=0; sl<scans[0][0].zsize(); sl++) {
          for (unsigned int i=0; i<fesz; i++) {
            Lambda[sl][i] = k*theta*KtK[sl][i] + l*tau*Lijk[0][0];
	    NEWMAT::CroutMatrix tmpLambda = Lambda[sl][i].AsNEWMAT();
            mu[sl][i] = mu_update(mu,K[sl][i],y[sl][i],tmpLambda,Lijk,k*theta,l*tau,sl,i);
          }
        }    
      }
      // Add to output
      for (int sl=0; sl<scans[0].zsize(); sl++) {
        add_to_slice(ovol[d],mu[sl],sl,row);
      }
      
    } // End directions
  } // End collections

  /*
  cout << "Computing K-matrices and stuff" << endl;
  // Pre-compute K-matrices and re-package data
  for (int k=0; k<scans[0][0].zsize(); k++) {
    NEWMAT::ColumnVector     tmpy(2*scans[0][0].ysize());

    K[k].resize(scans[0][0].xsize());
    KtK[k].resize(scans[0][0].xsize());
    y[k].resize(scans[0][0].xsize(),tmpy);
    yty[k].resize(scans[0][0].xsize(),0.0);
    for (int i=0; i<scans[0][0].xsize(); i++) {
      DispVec    dv(scans[0][0].ysize());
      dv.SetFromColumn(field,k,i);
      double sf = (datafile.PhaseEncodeVector(collections.IndexAt(0,0)))(2) * datafile.ReadOutTime(collections.IndexAt(0,0));
      K[k][i] = dv.GetSparseK_Matrix(sf);
      y[k][i].SubMatrix(1,scans[0][0].ysize(),1,1) = extract_row_col((scans[collections.ScanAt(0,0)])[0],k,i,false);
      sf = (datafile.PhaseEncodeVector(collections.IndexAt(0,1)))(2) * datafile.ReadOutTime(collections.IndexAt(0,1));
      K[k][i] &= dv.GetSparseK_Matrix(sf);
      KtK[k][i] = K[k][i].TransMultSelf();
      y[k][i].SubMatrix(scans[0][0].ysize()+1,2*scans[0][0].ysize(),1,1) = extract_row_col((scans[collections.ScanAt(0,1)])[0],k,i,false);
      yty[k][i] = (y[k][i].t()*y[k][i]).AsScalar();
    }
  }

  // Get blocks from regularisation matrix
  cout << "Computing regularisation matrix" << endl;
  std::vector<double>                       sf(2,scans[0][0].ydim()/scans[0][0].xdim());
  sf[1] = scans[0][0].ydim()/scans[0][0].zdim();
  std::vector<std::vector<SpMat<double> > > Lijk = GetLijk(scans[0][0].ysize(),sf);
  SpMat<double> Lii = GetLii(scans[0][0].ysize()); // For initialization of mu (images)

  // Initialize constants
  int m=scans[0][0].ysize();
  int n=2*m;
  int N=scans[0][0].xsize()*scans[0][0].zsize();
  // Set priors
  double l_0 = 1e-10; double tau_0 = 1e10;
  double k_0 = 1e-10; double theta_0 = 1e10;
  // Initialize all parameters to their priors
  double ll = l_0; double tau = tau_0;
  double kk = k_0; double theta = theta_0;
  std::vector<std::vector<SpMat<double> > >        Lambda(scans[0][0].zsize());
  std::vector<std::vector<NEWMAT::ColumnVector> >  mu(scans[0][0].zsize());
  for (int k=0; k<scans[0][0].zsize(); k++) {
    Lambda[k].resize(scans[0][0].xsize());
    mu[k].resize(scans[0][0].xsize());
    for (int i=0; i<scans[0][0].xsize(); i++) {
      Lambda[k][i] = kk*theta*KtK[k][i] + ll*tau*Lii;
      mu[k][i] = (kk*theta)*Lambda[k][i].AsNEWMAT().i()*K[k][i].TransMult(y[k][i]);
    }
  }

  for (int iter=0; iter<20; iter++) {

    cout << "Relative lambda = " << (ll*tau)/(kk*theta) << ",   phi = " << kk*theta << ",   lambda = " << ll*tau << endl;
    cout << "Updating phi" << endl;
    // Update parameters for phi (noise precision)
    kk = k_0 + n*N/2;
    theta = 1.0 / theta_0;
    for (int k=0; k<scans[0][0].zsize(); k++) {
      for (int i=0; i<scans[0][0].xsize(); i++) {
        theta += (0.5 * (mu[k][i].t()*(KtK[k][i]*mu[k][i]) + yty[k][i] + (Lambda[k][i].AsNEWMAT().i()*KtK[k][i].AsNEWMAT()).Trace())).AsScalar();
        theta -= (mu[k][i].t()*(K[k][i].TransMult(y[k][i]))).AsScalar();
      }
    }
    theta = 1.0/theta;

    cout << "Updating lambda" << endl;
    // Update parameters for lambda (regularisation weight)
    ll = l_0 + m*N/2;
    tau = tau_update(mu,Lambda,Lijk,tau_0);

    cout << "Updating image" << endl;
    // Update image and variance of image
    for (int k=0; k<scans[0][0].zsize(); k++) {
      for (int i=0; i<scans[0][0].xsize(); i++) {
        Lambda[k][i] = kk*theta*KtK[k][i] + ll*tau*Lijk[0][0];
        mu[k][i] = mu_update(mu,K[k][i],y[k][i],Lambda[k][i],Lijk,kk*theta,ll*tau,k,i);
      }
    }    
  }
  
  for (int k=0; k<ovol.zsize(); k++) {
    for (int i=0; i<ovol.xsize(); i++) {
      for (int j=0; j<ovol.ysize(); j++) {
	ovol[0](i,j,k) = mu[k][i](j+1);
      }
    }
  }
  */

  return(ovol);
}

double tau_update_contrib(const std::vector<std::vector<NEWMAT::ColumnVector> >&        mu,
                          const std::vector<std::vector<MISCMATHS::SpMat<double> > >&   Lambda,
		          const std::vector<std::vector<MISCMATHS::SpMat<double> > >&   Lijk)
{
  int nsl = mu.size();      // Number of slices
  int nrow = mu[0].size();  // Number of frequency encodes
  double tauc = 0.0;
  for (int k=0; k<nsl; k++) {
    for (int i=0; i<nrow; i++) {
      tauc += (0.5*((Lambda[k][i].AsNEWMAT().i()*Lijk[0][0].AsNEWMAT()).Trace() + mu[k][i].t()*(Lijk[0][0]*mu[k][i]))).AsScalar();
      for (int ii=max(0,i-2); ii<=min(nrow-1,i+2); ii++) {
        if (ii != i) tauc += (0.5*mu[k][i].t()*(Lijk[0][abs(ii-i)]*mu[k][ii])).AsScalar();
      }
      for (int kk=max(0,k-2); kk<=min(nsl-1,k+2); kk++) {
        if (kk != k) tauc += (0.5*mu[k][i].t()*(Lijk[abs(kk-k)][0]*mu[kk][i])).AsScalar();
      }
      if (k && i) tauc += (0.5*mu[k][i].t()*(Lijk[1][1]*mu[k-1][i-1])).AsScalar();
      if (k && i<(nrow-1)) tauc += (0.5*mu[k][i].t()*(Lijk[1][1]*mu[k-1][i+1])).AsScalar();
      if (k<(nsl-1) && i<(nrow-1)) tauc += (0.5*mu[k][i].t()*(Lijk[1][1]*mu[k+1][i+1])).AsScalar();
      if (k<(nsl-1) && i) tauc += (0.5*mu[k][i].t()*(Lijk[1][1]*mu[k+1][i-1])).AsScalar();
    }
  }

  return(tauc);
}
double tau_update(const std::vector<std::vector<std::vector<NEWMAT::ColumnVector> > >&  mu,
                  const std::vector<std::vector<MISCMATHS::SpMat<double> > >&           Lambda,
		  const std::vector<std::vector<MISCMATHS::SpMat<double> > >&           Lijk,
                  double                                                                tau_0)
{
  int nd = mu.size();          // Number of directions
  int nsl = mu[0].size();      // Number of slices
  int nrow = mu[0][0].size();  // Number of frequency encodes
  double tau = 1.0 / tau_0;
  for (int k=0; k<nsl; k++) {
    for (int i=0; i<nrow; i++) {
      tau += 0.5 * double(nd) * (Lambda[k][i].AsNEWMAT().i()*Lijk[0][0].AsNEWMAT()).Trace();
      for (int d=0; d<nd; d++) {
	tau += 0.5 * (mu[d][k][i].t()*(Lijk[0][0]*mu[d][k][i])).AsScalar();
	for (int ii=max(0,i-2); ii<=min(nrow-1,i+2); ii++) {
	  if (ii != i) tau += (0.5*mu[d][k][i].t()*(Lijk[0][abs(ii-i)]*mu[d][k][ii])).AsScalar();
	}
	for (int kk=max(0,k-2); kk<=min(nsl-1,k+2); kk++) {
	  if (kk != k) tau += (0.5*mu[d][k][i].t()*(Lijk[abs(kk-k)][0]*mu[d][kk][i])).AsScalar();
	}
	if (k && i) tau += (0.5*mu[d][k][i].t()*(Lijk[1][1]*mu[d][k-1][i-1])).AsScalar();
	if (k && i<(nrow-1)) tau += (0.5*mu[d][k][i].t()*(Lijk[1][1]*mu[d][k-1][i+1])).AsScalar();
	if (k<(nsl-1) && i<(nrow-1)) tau += (0.5*mu[d][k][i].t()*(Lijk[1][1]*mu[d][k+1][i+1])).AsScalar();
	if (k<(nsl-1) && i) tau += (0.5*mu[d][k][i].t()*(Lijk[1][1]*mu[d][k+1][i-1])).AsScalar();
      }
    }
  }

  return(1.0 / tau);
}

double tau_update(const std::vector<std::vector<NEWMAT::ColumnVector> >&        mu,
                  const std::vector<std::vector<MISCMATHS::SpMat<double> > >&   Lambda,
		  const std::vector<std::vector<MISCMATHS::SpMat<double> > >&   Lijk,
                  double                                                        tau_0)
{
  int nsl = mu.size();
  int nrow = mu[0].size();
  double tau = 1.0 / tau_0;
  for (int k=0; k<nsl; k++) {
    for (int i=0; i<nrow; i++) {
      tau += (0.5*((Lambda[k][i].AsNEWMAT().i()*Lijk[0][0].AsNEWMAT()).Trace() + mu[k][i].t()*(Lijk[0][0]*mu[k][i]))).AsScalar();
      for (int ii=max(0,i-2); ii<=min(nrow-1,i+2); ii++) {
        if (ii != i) tau += (0.5*mu[k][i].t()*(Lijk[0][abs(ii-i)]*mu[k][ii])).AsScalar();
      }
      for (int kk=max(0,k-2); kk<=min(nsl-1,k+2); kk++) {
        if (kk != k) tau += (0.5*mu[k][i].t()*(Lijk[abs(kk-k)][0]*mu[kk][i])).AsScalar();
      }
      if (k && i) tau += (0.5*mu[k][i].t()*(Lijk[1][1]*mu[k-1][i-1])).AsScalar();
      if (k && i<(nrow-1)) tau += (0.5*mu[k][i].t()*(Lijk[1][1]*mu[k-1][i+1])).AsScalar();
      if (k<(nsl-1) && i<(nrow-1)) tau += (0.5*mu[k][i].t()*(Lijk[1][1]*mu[k+1][i+1])).AsScalar();
      if (k<(nsl-1) && i) tau += (0.5*mu[k][i].t()*(Lijk[1][1]*mu[k+1][i-1])).AsScalar();
    }
  }
  tau = 1.0/tau;

  return(tau);
}
                    
NEWMAT::ColumnVector mu_update(const std::vector<std::vector<NEWMAT::ColumnVector> >&        mu,
                               const MISCMATHS::SpMat<double>&                               K,
		               const NEWMAT::ColumnVector&                                   y,
                               const NEWMAT::CroutMatrix&                                    Lambda,
                               const std::vector<std::vector<MISCMATHS::SpMat<double> > >&   Lijk,
                               double                                                        phi,
                               double                                                        lambda,
                               int                                                           sl,
                               int                                                           row)
{
  int nsl = mu.size();
  int nrow = mu[0].size();
  NEWMAT::ColumnVector  newmu = phi*K.TransMult(y);

  for (int i=max(0,row-2); i<=min(nrow-1,row+2); i++) {
    if (i != row) newmu -= lambda*Lijk[0][abs(i-row)]*mu[sl][i];
  }
  for (int i=max(0,sl-2); i<=min(nsl-1,sl+2); i++) {
    if (i != sl) newmu -= lambda*Lijk[abs(i-sl)][0]*mu[i][row];
  }
  if (sl && row) newmu -= lambda*Lijk[1][1]*mu[sl-1][row-1];
  if (sl && row<(nrow-1)) newmu -= lambda*Lijk[1][1]*mu[sl-1][row+1];
  if (sl<(nsl-1) && row<(nrow-1)) newmu -= lambda*Lijk[1][1]*mu[sl+1][row+1];
  if (sl<(nsl-1) && row) newmu -= lambda*Lijk[1][1]*mu[sl+1][row-1];

  newmu = Lambda.i()*newmu;
  
  return(newmu);
}
                               
NEWIMAGE::volume4D<float> vb_resample_2D(const TopupDatafileReader&                       datafile, 
                                         const TopupFileReader&                           topupfile, 
                                         const std::vector<unsigned int>&                 inindices, 
                                         NEWIMAGE::volume<char>&                          mask,
                                         std::vector<NEWIMAGE::volume4D<float> >&         scans)
{
  NEWIMAGE::volume<float>   field = topupfile.FieldAsVolume();
  TopupCollections          collections(inindices,datafile);
  NEWIMAGE::volume4D<float> ovol = scans[0];
  ovol = 0.0;

  cout << "NCollections = " << collections.NCollections() << endl;
  cout << "NScans = " << collections.NScans(0) << endl;

  for (unsigned int c=0; c<collections.NCollections(); c++) {
    // Find PE direction for this collection
    bool row = false;
    unsigned int pesz = scans[0].ysize();
    unsigned int fesz = scans[0].xsize();
    if (datafile.PhaseEncodeVector(collections.IndexAt(c,0))(1)) { // If PE in x-dir for this collection
      row = true;
      pesz = scans[0].xsize();
      fesz = scans[0].ysize();
    }
    std::vector<MISCMATHS::SpMat<double> > K(fesz);    // All K-matrices for one slice
    std::vector<MISCMATHS::SpMat<double> > KtK(fesz);  // All KtK-matrices for one slice
    std::vector<NEWMAT::ColumnVector> y(fesz);         // All data for one slice and direction
    std::vector<double>               yty(fesz,0.0);   // Well..
    double sf = scans[0][0].ydim()/scans[0][0].xdim();
    if (row) sf = 1.0/sf;
    std::vector<SpMat<double> >  Lij = GetLij(pesz,sf);// Regularisation sub-matrices
    MISCMATHS::SpMat<double>     Lii = GetLii(pesz);   // For initialization of mu
    int m=pesz;                                        // Number of elements along PE
    int N=fesz;                                        // Number if frequency encodes in one slice
    int n=collections.NScans(c)*m;                     // Number of data-points for one "column"
    for (int sl=0; sl<scans[0].zsize(); sl++) {
      cout << "slice = " << sl << endl;
      // Make all the K-matrices for the current slice
      for (unsigned int i=0; i<fesz; i++) { 
	DispVec         dv(pesz);
	if (row) dv.SetFromRow(field,sl,i);
	else dv.SetFromColumn(field,sl,i);
        for (unsigned int a=0; a<collections.NScans(c); a++) { // Loop over acquisitions in collection
  	  double psf = (datafile.PhaseEncodeVector(collections.IndexAt(c,a)))(2) * datafile.ReadOutTime(collections.IndexAt(c,a));
          if (row) psf = (datafile.PhaseEncodeVector(collections.IndexAt(c,a)))(1) * datafile.ReadOutTime(collections.IndexAt(c,a));
          if (!a) K[i] = dv.GetSparseK_Matrix(psf);
          else K[i] &= dv.GetSparseK_Matrix(psf);
	}
        KtK[i] = K[i].TransMultSelf();
      }
      for (int d=0; d<ovol.tsize(); d++) {
        // Gather all data for this slice and direction
        for (unsigned int i=0; i<fesz; i++) { 
          for (unsigned int a=0; a<collections.NScans(c); a++) {
            if (!a) y[i] = extract_row_col((scans[collections.ScanAt(c,a)])[d],sl,i,row);
	    else y[i] &= extract_row_col((scans[collections.ScanAt(c,a)])[d],sl,i,row);
	  }
          yty[i] = (y[i].t()*y[i]).AsScalar();
	}
        // Estimate phi, lambda and images for this slice and direction
        // Set priors
        double l_0 = 1e-10; double tau_0 = 1e10;
        double k_0 = 1e-10; double theta_0 = 1e10;
        // Initialize all parameters to their priors
        double l = l_0; double tau = tau_0;
        double k = k_0; double theta = theta_0;
        std::vector<SpMat<double> >        Lambda(N);
        std::vector<NEWMAT::ColumnVector>  mu(N);
        for (int i=0; i<N; i++) {
	  Lambda[i] = k*theta*KtK[i] + l*tau*Lii;
	  mu[i] = (k*theta)*Lambda[i].AsNEWMAT().i()*K[i].TransMult(y[i]);
        }
        for (int iter=0; iter<50; iter++) {
          cout << "relative lambda = " << (l*tau)/(k*theta) << ",   phi = " << k*theta << ",   lambda = " << l*tau << endl;
	  // Update parameters for phi (noise precision)
          k = k_0 + n*N/2;
	  theta = 1.0 / theta_0;
	  for (int i=0; i<N; i++) {
	    theta += (0.5 * (mu[i].t()*(KtK[i]*mu[i]) + yty[i] + (Lambda[i].AsNEWMAT().i()*KtK[i].AsNEWMAT()).Trace())).AsScalar();
	    theta -= (mu[i].t()*(K[i].TransMult(y[i]))).AsScalar();
	  }
	  theta = 1.0/theta;
	  // Update parameters for lambda (regularisation weight)
	  l = l_0 + m*N/2;
	  tau = 1.0 / tau_0;
	  for (int i=0; i<N; i++) {
	    tau += (0.5*((Lambda[i].AsNEWMAT().i()*Lij[0].AsNEWMAT()).Trace() + mu[i].t()*(Lij[0]*mu[i]))).AsScalar();
	    for (int j=max(0,i-2); j<=min(N-1,i+2); j++) {
	      if (j != i) {
	        tau += (0.5*mu[i].t()*(Lij[abs(j-i)]*mu[j])).AsScalar();
	      }
	    }
	  }
	  tau = 1.0/tau;
	  // Update image and variance of image
	  for (int i=0; i<N; i++) {
	    Lambda[i] = k*theta*KtK[i] + l*tau*Lij[0];
	    mu[i] = (k*theta)*K[i].TransMult(y[i]);
	    for (int j=max(0,i-2); j<=min(N-1,i+2); j++) {
	      if (j != i) mu[i] -= (l*tau)*Lij[abs(j-i)]*mu[j];
	    }
	    mu[i] = Lambda[i].AsNEWMAT().i()*mu[i];
	  }    
        }
        cout << "slice = " << sl << ",   relative lambda = " << (l*tau)/(k*theta) << ",   phi = " << k*theta << ",   lambda = " << l*tau << endl;
        add_to_slice(ovol[d],mu,sl,row);
      } // End direction
    } // End slice
  } // End collection	

  return(ovol);
}

/////////////////////////////////////////////////////////////////////
//
// Returns the diagonal and off-diagonal (those that contain non-zero
// values) blocks of a matrix A'*A where A is a matrix implementing
// a 3D 2nd derivative. If the full 3D volume is lxmxn, where the 
// first direction is of size sz (i.e. sz=l) and represents the 
// direction of phase-encoding then each "block" is an lxl matrix.
// It returns six such blocks, where the three in the first vector
// corresponds to the diagonal block and the two off-diagonal blocks
// when moving in the 2nd direction of the volume. The two in
// the second vector corresponds to moving one step in the third
// direction and the one in the third vector corresponds to moving
// two steps in the third direction.
// The scale-factors in sf corresponds to a scaling to the voxel-sizes
// in the 2nd and 3rd direction so that the scale of the derivative
// is the same in all directions (i.e. per-voxel-in-first-direction).
// If your volume has voxel-sizes vxs then sf[0]=vxs[0]/vxs[1] and
// sf[1]=vxs[0]/vxs[2].
//
/////////////////////////////////////////////////////////////////////

std::vector<std::vector<MISCMATHS::SpMat<double> > > GetLijk(int sz, const std::vector<double> sf)
{
  SpMat<double> tmpL(sz,sz);
  std::vector<std::vector<MISCMATHS::SpMat<double> > >  Lijk(3);
  Lijk[0].resize(3,tmpL); Lijk[1].resize(2,tmpL); Lijk[2].resize(1,tmpL);

  { // Lii
    double A11 = Sqr(2+2*sf[0]+2*sf[1]) + 2*Sqr(-1.0) + 2*Sqr(-sf[0]) + 2*Sqr(-sf[1]);
    double A12 = 2*(-1*(2+2*sf[0]*sf[1]));
    double A13 = 1.0;
    Lijk[0][0].Set(1,sz-1,A13); Lijk[0][0].Set(1,sz,A12); 
    Lijk[0][0].Set(1,1,A11); Lijk[0][0].Set(1,2,A12); Lijk[0][0].Set(1,3,A13);
    Lijk[0][0].Set(2,sz,A13); Lijk[0][0].Set(2,1,A12); 
    Lijk[0][0].Set(2,2,A11); Lijk[0][0].Set(2,3,A12); Lijk[0][0].Set(2,4,A13);
    for (int i=3; i<=sz-2; i++) {
      Lijk[0][0].Set(i,i-2,A13); Lijk[0][0].Set(i,i-1,A12); 
      Lijk[0][0].Set(i,i,A11); Lijk[0][0].Set(i,i+1,A12); Lijk[0][0].Set(i,i+2,A13);
    }
    Lijk[0][0].Set(sz-1,sz-3,A13); Lijk[0][0].Set(sz-1,sz-2,A12); 
    Lijk[0][0].Set(sz-1,sz-1,A11); Lijk[0][0].Set(sz-1,sz,A12); Lijk[0][0].Set(sz-1,1,A13);
    Lijk[0][0].Set(sz,sz-2,A13); Lijk[0][0].Set(sz,sz-1,A12); 
    Lijk[0][0].Set(sz,sz,A11); Lijk[0][0].Set(sz,1,A12); Lijk[0][0].Set(sz,2,A13);
  }
  { // Lij: One off in the first direction
    double A11 = 2*((-sf[0])*(2+2*sf[0]+2*sf[1]));
    double A12 = 2*(-1)*(-sf[0]);
    Lijk[0][1].Set(1,sz,A12); Lijk[0][1].Set(1,1,A11); Lijk[0][1].Set(1,2,A12);
    for (int i=2; i<=sz-1; i++) {
      Lijk[0][1].Set(i,i-1,A12); Lijk[0][1].Set(i,i,A11); Lijk[0][1].Set(i,i+1,A12);
    }
    Lijk[0][1].Set(sz,sz-1,A12); Lijk[0][1].Set(sz,sz,A11); Lijk[0][1].Set(sz,1,A12);
  }
  { // Lij: Two off in the first direction
    double A11 = Sqr(-sf[0]);
    for (int i=1; i<=sz; i++) {
      Lijk[0][2].Set(i,i,A11);
    }
  }
  { // Lij: One off in the second direction
    double A11 = 2*((-sf[1])*(2+2*sf[0]+2*sf[1]));
    double A12 = 1*(-1)*(-sf[1]);
    Lijk[1][0].Set(1,sz,A12); Lijk[1][0].Set(1,1,A11); Lijk[1][0].Set(1,2,A12);
    for (int i=2; i<=sz-1; i++) {
      Lijk[1][0].Set(i,i-1,A12); Lijk[1][0].Set(i,i,A11); Lijk[1][0].Set(i,i+1,A12);
    }
    Lijk[1][0].Set(sz,sz-1,A12); Lijk[1][0].Set(sz,sz,A11); Lijk[1][0].Set(sz,1,A12);
  }
  { // Lij: One off in BOTH directions
    double A11 = 2*(-sf[0])*(-sf[1]);
    for (int i=1; i<=sz; i++) {
      Lijk[1][1].Set(i,i,A11);
    }
  }    
  { // Lij: Two off in the second direction
    double A11 = Sqr(-sf[1]);
    for (int i=1; i<=sz; i++) {
      Lijk[2][0].Set(i,i,A11);
    }
  }    
  return(Lijk);
} 

/////////////////////////////////////////////////////////////////////
//
// Returns the diagonal and off-diagonal (those that contain non-zero
// values) blocks of a matrix A'*A where A is a matrix implementing
// a 2D 2nd derivative. If the full 2D slice is mxn, where the 
// first direction is of size sz (i.e. sz=m) and represents the 
// direction of phase-encoding then each "block" is an mxm matrix.
// It returns three such blocks in a vector corresponding to the
// diagonal block and to moving one and two steps in the 2n direction.
//
// The scale-factor in sf corresponds to a scaling to the voxel-size
// in the 2nd direction so that the scale of the derivative is the 
// same in both directions (i.e. per-voxel-in-first-direction).
// If your slice has voxel-sizes vxs then sf=vxs[0]/vxs[1].
//
/////////////////////////////////////////////////////////////////////

std::vector<SpMat<double> > GetLij(int sz, double sf)
{
  SpMat<double>               tmp(sz,sz);
  std::vector<SpMat<double> > Lij(3,tmp);

  { // Lii
    double A11 = Sqr(2+2*sf) + 2*Sqr(-1.0) + 2*Sqr(-sf);
    double A12 = 2 * (-1) * (2+2*sf);
    double A13 = 1;
    Lij[0].Set(1,sz-1,A13); Lij[0].Set(1,sz,A12); Lij[0].Set(1,1,A11);
    Lij[0].Set(1,2,A12); Lij[0].Set(1,3,A13); 
    Lij[0].Set(2,sz,A13); Lij[0].Set(2,1,A12); Lij[0].Set(2,2,A11);
    Lij[0].Set(2,3,A12); Lij[0].Set(2,4,A13); 
    for (int i=3; i<=sz-2; i++) {
      Lij[0].Set(i,i-2,A13); Lij[0].Set(i,i-1,A12); Lij[0].Set(i,i,A11);
      Lij[0].Set(i,i+1,A12); Lij[0].Set(i,i+2,A13); 
    }
    Lij[0].Set(sz-1,sz-3,A13); Lij[0].Set(sz-1,sz-2,A12); Lij[0].Set(sz-1,sz-1,A11);
    Lij[0].Set(sz-1,sz,A12); Lij[0].Set(sz-1,1,A13); 
    Lij[0].Set(sz,sz-2,A13); Lij[0].Set(sz,sz-1,A12); Lij[0].Set(sz,sz,A11);
    Lij[0].Set(sz,1,A12); Lij[0].Set(sz,2,A13); 
  }
  { // Lij: One off
    double A11 = 2 * (-sf) * (2+2*sf);
    double A12 = 2 * (-1) * (-sf);
    Lij[1].Set(1,sz,A12); Lij[1].Set(1,1,A11); Lij[1].Set(1,2,A12);
    for (int i=2; i<=sz-1; i++) {
      Lij[1].Set(i,i-1,A12); Lij[1].Set(i,i,A11); Lij[1].Set(i,i+1,A12);
    }
    Lij[1].Set(sz,sz-1,A12); Lij[1].Set(sz,sz,A11); Lij[1].Set(sz,1,A12);
  }
  { // Lij: Two off
    double A11 = Sqr(-sf);
    for (int i=1; i<=sz; i++) {
      Lij[2].Set(i,i,A11);
    }
  }
  return(Lij);
}

MISCMATHS::SpMat<double> GetLii(int sz)
{
  SpMat<double>  Lii(sz,sz);

  double A11 = 6;
  double A12 = -4;
  double A13 = 1;
  Lii.Set(1,sz-1,A13); Lii.Set(1,sz,A12); Lii.Set(1,1,A11);
  Lii.Set(1,2,A12); Lii.Set(1,3,A13); 
  Lii.Set(2,sz,A13); Lii.Set(2,1,A12); Lii.Set(2,2,A11);
  Lii.Set(2,3,A12); Lii.Set(2,4,A13); 
  for (int i=3; i<=sz-2; i++) {
    Lii.Set(i,i-2,A13); Lii.Set(i,i-1,A12); Lii.Set(i,i,A11);
    Lii.Set(i,i+1,A12); Lii.Set(i,i+2,A13); 
  }
  Lii.Set(sz-1,sz-3,A13); Lii.Set(sz-1,sz-2,A12); Lii.Set(sz-1,sz-1,A11);
  Lii.Set(sz-1,sz,A12); Lii.Set(sz-1,1,A13); 
  Lii.Set(sz,sz-2,A13); Lii.Set(sz,sz-1,A12); Lii.Set(sz,sz,A11);
  Lii.Set(sz,1,A12); Lii.Set(sz,2,A13);

  return(Lii); 
}

/*
MISCMATHS::SpMat<double> GetLii(int sz)
{
  SpMat<double> Lii(sz,sz);
  Lii.Set(1,sz-1,1); Lii.Set(1,sz,-8); Lii.Set(1,1,20); Lii.Set(1,2,-8); Lii.Set(1,3,1);
  Lii.Set(2,sz,1); Lii.Set(2,1,-8); Lii.Set(2,2,20); Lii.Set(2,3,-8); Lii.Set(2,4,1);
  for (int i=3; i<=sz-2; i++) {
    Lii.Set(i,i-2,1); Lii.Set(i,i-1,-8); Lii.Set(i,i,20); Lii.Set(i,i+1,-8); Lii.Set(i,i+2,1);
  }
  Lii.Set(sz-1,sz-3,1); Lii.Set(sz-1,sz-2,-8); Lii.Set(sz-1,sz-1,20); Lii.Set(sz-1,sz,-8); Lii.Set(sz-1,1,1);
  Lii.Set(sz,sz-2,1); Lii.Set(sz,sz-1,-8); Lii.Set(sz,sz,20); Lii.Set(sz,1,-8); Lii.Set(sz,2,1);

  return(Lii);    
}

std::vector<SpMat<double> > GetLij(int sz)
{
  SpMat<double>               tmp(sz,sz);
  std::vector<SpMat<double> > Lij(2,tmp);

  Lij[0].Set(1,sz,2); Lij[0].Set(1,1,-8); Lij[0].Set(1,2,2);
  for (int i=2; i<=sz-1; i++) {
    Lij[0].Set(i,i-1,2); Lij[0].Set(i,i,-8); Lij[0].Set(i,i+1,2);
  }
  Lij[0].Set(sz,sz-1,2); Lij[0].Set(sz,sz,-8); Lij[0].Set(sz,1,2);

  for (int i=1; i<=sz; i++) Lij[1].Set(i,i,1);

  return(Lij);
}

MISCMATHS::SpMat<double> GetS_Matrix(const std::vector<int>&   isz,
                                     const std::vector<bool>&  wrap)
{
  unsigned int mn=1;
  for (unsigned int i=0; i<isz.size(); i++) mn *= isz[i];
  SpMat<double> S(mn,mn);

  // Do first direction
  int nline = mn/isz[0];
  for (int line = 1; line<=nline; line++) {
    int offs = (line-1)*isz[0];
    // First row
    if (wrap[0]) {
      S.Set(offs+1,offs+isz[0],-1);
      S.Set(offs+1,offs+1,2);
      S.Set(offs+1,offs+2,-1);
    }
    else {
      S.Set(offs+1,offs+1,-1);
      S.Set(offs+1,offs+2,1);
    }
    // Rows 2 to N-1
    for (int i=2; i<isz[0]; i++) {
      S.Set(offs+i,offs+i-1,-1);
      S.Set(offs+i,offs+i,2);
      S.Set(offs+i,offs+i+1,-1);
    }
    // Last row
    if (wrap[0]) {
      S.Set(offs+isz[0],offs+isz[0]-1,-1);
      S.Set(offs+isz[0],offs+isz[0],2);
      S.Set(offs+isz[0],offs+1,-1);
    }
    else {
      S.Set(offs+isz[0],offs+isz[0]-1,-1);
      S.Set(offs+isz[0],offs+isz[0],1);
    }
  }

  // Do second direction (if there is one)
  if (isz.size() > 1) {
    int nline = mn/isz[0];
    // First line
    for (int i=1; i<=isz[0]; i++) {
      if (wrap[1]) {
        S.AddTo(i,i+(isz[1]-1)*isz[0],-1);
        S.AddTo(i,i,2);
        S.AddTo(i,i+isz[0],-1);        
      }
      else {
        S.AddTo(i,i,-1);
        S.AddTo(i,i+isz[0],1);        
      }
    }
    // Lines 2 to N-1
    for (int line=2; line<nline; line++) {
      int offs = (line-1)*isz[0];
      for (int i=1; i<=isz[0]; i++) {
        S.AddTo(offs+i,offs+i-isz[0],-1);
        S.AddTo(offs+i,offs+i,2);
        S.AddTo(offs+i,offs+i+isz[0],-1);
      }
    }
    // Last line
    int offs = (isz[1]-1)*isz[0];
    for (int i=1; i<=isz[0]; i++) {
      if (wrap[1]) {
        S.AddTo(offs+i,offs+i-1,-1);
        S.AddTo(offs+i,offs+i,2);
        S.AddTo(offs+i,i,-1);        
      }
      else {
        S.AddTo(offs+i,offs+i-1,-1);
        S.AddTo(offs+i,offs+i,1);        
      }
    }
  }
  
  return(S);
}
*/

NEWIMAGE::volume4D<float> jac_resample(const TopupDatafileReader&                       datafile, 
                                       const TopupFileReader&                           topupfile, 
                                       const std::vector<unsigned int>&                 inindices, 
                                       NEWIMAGE::volume<char>&                          mask,
                                       std::vector<NEWIMAGE::volume4D<float> >&         scans,
                                       NEWIMAGE::interpolation                          interp)
{
  NEWIMAGE::volume<float>  tmp = scans[0][0]; 
  NEWIMAGE::volume4D<float> df(tmp.xsize(),tmp.ysize(),tmp.zsize(),3);
  copybasicproperties(tmp,df[0]); copybasicproperties(tmp,df[1]); copybasicproperties(tmp,df[2]);
  NEWIMAGE::volume4D<float> ovol = scans[0];
  ovol = 0.0;
  for (unsigned int i=0; i<scans.size(); i++) {        // Loop over acquistions
    scans[i].setinterpolationmethod(interp);
    if (interp == NEWIMAGE::spline) scans[i].setsplineorder(3);
    scans[i].setextrapolationmethod(NEWIMAGE::periodic);
    if ((datafile.PhaseEncodeVector(inindices[i]))(1) && !(datafile.PhaseEncodeVector(inindices[i]))(2)) scans[i].setextrapolationvalidity(true,false,false);
    else if (!(datafile.PhaseEncodeVector(inindices[i]))(1) && (datafile.PhaseEncodeVector(inindices[i]))(2)) scans[i].setextrapolationvalidity(false,true,false);
    else scans[i].setextrapolationvalidity(false,false,false);
    df[0] = topupfile.FieldAsVolume();
    df[1] = float((datafile.PhaseEncodeVector(inindices[i]))(2) * datafile.ReadOutTime(inindices[i]) * ovol.ydim()) * df[0];
    df[0] *= ((datafile.PhaseEncodeVector(inindices[i]))(1) * datafile.ReadOutTime(inindices[i]) * ovol.xdim());
    df[2] = 0.0;
    BASISFIELD::splinefield xcomp = topupfile.Field();
    BASISFIELD::splinefield ycomp = xcomp;
    BASISFIELD::splinefield zcomp = xcomp;
    xcomp.ScaleField((datafile.PhaseEncodeVector(inindices[i]))(1) * datafile.ReadOutTime(inindices[i]) * ovol.xdim());
    ycomp.ScaleField((datafile.PhaseEncodeVector(inindices[i]))(2) * datafile.ReadOutTime(inindices[i]) * ovol.ydim());
    zcomp.ScaleField(0.0);
    NEWIMAGE::volume<float>  jac = tmp;
    NEWIMAGE::volume<char>   tmpmask = mask;
    jac = 0.0;
    NEWIMAGE::deffield2jacobian(xcomp,ycomp,zcomp,jac);
    NEWMAT::Matrix  M = topupfile.MoveMatrix(inindices[i]);
    for (int j=0; j<scans[0].tsize(); j++) {           // Loop over diffusion gradient direction    
      tmp = 0.0;
      tmpmask = 0;
      general_transform((scans[i])[j],M,df,tmp,tmpmask);
      mask *= tmpmask;
      ovol[j] += tmp*jac;  
    }      
  }
  ovol /= scans.size();
  return(ovol);
}

NEWIMAGE::volume4D<float> lsr_resample(const TopupDatafileReader&                       datafile, 
                                       const TopupFileReader&                           topupfile, 
                                       const std::vector<unsigned int>&                 inindices, 
                                       NEWIMAGE::volume<char>&                          mask,
                                       std::vector<NEWIMAGE::volume4D<float> >&         scans)
{
  NEWIMAGE::volume<float>      field = topupfile.FieldAsVolume();      // Get off-resonance field
  TopupCollections             collections(inindices,datafile);       // Divide data into "Collections"
  std::vector<NEWMAT::Matrix>  mis_map(collections.NCollections());    // Indicator of "missing" data
  NEWIMAGE::volume4D<float>    ovol = scans[0];                        // N.B. that scans is a vector
  ovol = 0.0;

  for (unsigned int c=0; c<collections.NCollections(); c++) {          // Loop over collections
    bool row = false;
    unsigned int sz = scans[0].ysize();
    if ((datafile.PhaseEncodeVector(collections.IndexAt(c,0)))(1)) {   // If this collection does PE in x-direction
      row = true;
      sz = scans[0].xsize();
    }
    mis_map[c].ReSize(scans[0].zsize(),(row) ? scans[0].ysize() : scans[0].xsize());
    mis_map[c] = 0;
    DispVec               dv(sz);
    NEWMAT::Matrix        K(collections.NScans(c)*sz,sz);             // Vertical concat of K-matrices
    NEWMAT::Matrix        Y(collections.NScans(c)*sz,ovol.tsize());   // All data-vectors for a given row
    NEWMAT::Matrix        StS=dv.GetS_Matrix(false);                  // Laplacian for regularisation
    StS = StS.t()*StS;
    for (int k=0; k<scans[0].zsize(); k++) {
      // cout << "k = " << k << endl;
      for (int ij=0; ij<((row) ? scans[0].ysize() : scans[0].xsize()); ij++) {
        if (row_col_is_alright(mask,k,ij,row)) {  // If there is not "too much" missing data in this row/column
          // Create K-matrix for this row/column (common to all diffusion directions)
          if (row) dv.SetFromRow(field,k,ij);
          else dv.SetFromColumn(field,k,ij);
          for (unsigned int a=0; a<collections.NScans(c); a++) { // Loop over acquisitions in collection
            double sf = (datafile.PhaseEncodeVector(collections.IndexAt(c,a)))(2) * datafile.ReadOutTime(collections.IndexAt(c,a));
            if (row) sf = (datafile.PhaseEncodeVector(collections.IndexAt(c,a)))(1) * datafile.ReadOutTime(collections.IndexAt(c,a));
            K.Rows(a*sz+1,(a+1)*sz) = dv.GetK_Matrix(sf);
          }
          // Prepare for multiple matrix-solve
          NEWMAT::Matrix       KtK = K.t()*K + 0.01*StS;
          NEWMAT::CroutMatrix  XtX = KtK;
          // Fill a matrix with data from all acquisitions and directions
          for (int d=0; d<ovol.tsize(); d++) {                              // Loop over diffusion directions
            for (unsigned int a=0; a<collections.NScans(c); a++) {  // Copy data from each acquisition
              // cout << "Indexing the scan array with " << (collected_indices[c])[a]-1 << ", direction = " << d << endl;
              Y.SubMatrix(a*sz+1,(a+1)*sz,d+1,d+1) = extract_row_col((scans[collections.ScanAt(c,a)])[d],k,ij,row);
            }
          }
          NEWMAT::Matrix KtY = K.t()*Y;
          NEWMAT::Matrix B = XtX.i() * KtY;
          add_to_rows_cols(ovol,B,k,ij,row);  // Write into output for all directions
        }
        else {  // Indicate missing data
          mis_map[c](k+1,ij+1) = 1;
        }
      }
    }
  }
  // Mask out rows/columns with "too much" missing data
  for (unsigned int c=0; c<collections.NCollections(); c++) {             // Loop over collections
    bool row = false;
    if ((datafile.PhaseEncodeVector(collections.IndexAt(c,0)))(1)) row = true;  // If this collection does PE in x-direction
    zero_out_rows_cols(ovol,mis_map[c],row);  
  }

  return(ovol);    
}

void resample_using_movement_parameters(const TopupDatafileReader&                      datafile,
                                        const TopupFileReader&                          topupfile,
                                        const std::vector<unsigned int>&                inindices,
                                        NEWIMAGE::volume<char>&                         mask,
                                        std::vector<NEWIMAGE::volume4D<float> >&        scans,
                                        NEWIMAGE::interpolation                         interp)
{
  NEWIMAGE::volume<char>  tmpmask = mask;
  NEWIMAGE::volume<float> tmp = scans[0][0];
  tmp = 0.0;

  for (unsigned int i=0; i<scans.size(); i++) {        // Loop over acquistions
    scans[i].setinterpolationmethod(interp);
    if (interp == NEWIMAGE::spline) scans[i].setsplineorder(3);
    scans[i].setextrapolationmethod(NEWIMAGE::periodic);
    if ((datafile.PhaseEncodeVector(inindices[i]))(1) && !(datafile.PhaseEncodeVector(inindices[i]))(2)) scans[i].setextrapolationvalidity(true,false,false);
    else if (!(datafile.PhaseEncodeVector(inindices[i]))(1) && (datafile.PhaseEncodeVector(inindices[i]))(2)) scans[i].setextrapolationvalidity(false,true,false);
    else scans[i].setextrapolationvalidity(false,false,false);
    NEWMAT::Matrix  M = topupfile.MoveMatrix(inindices[i]);
    for (int j=0; j<scans[0].tsize(); j++) {           // Loop over diffusion gradient direction    
      tmp = 0.0;
      tmpmask = 0;
      affine_transform((scans[i])[j],M,tmp,tmpmask);
      mask *= tmpmask;
      (scans[i])[j] = tmp;  // Overwrite original data  
    }      
  }
  return;
}

std::vector<unsigned int> parse_commaseparated_numbers(const std::string& list)
{
  std::vector<std::string> str_list = parse_commaseparated_list(list);
  std::vector<unsigned int> number_list(str_list.size(),0);
  for (unsigned int i=0; i<str_list.size(); i++) {
    number_list[i] = atoi(str_list[i].c_str());
  }
  
  return(number_list);
}

std::vector<std::string> parse_commaseparated_list(const std::string&  list)
{
  std::vector<std::string> ostr;

  size_t cur_pos = 0;
  size_t new_pos = 0;
  unsigned int n=0;
  while ((new_pos = list.find_first_of(',',cur_pos)) != string::npos) {
    ostr.resize(++n);
    ostr[n-1] = list.substr(cur_pos,new_pos-cur_pos);
    cur_pos = new_pos+1;        
  }
  ostr.resize(++n);
  ostr[n-1] = list.substr(cur_pos,string::npos);

  return(ostr);
} 

// Makes sure that data are suited for least-squares restoration and 
// divides it into groups depending on phase-encode vectors

TopupCollections::TopupCollections(const std::vector<unsigned int>&  indx,
                                   const TopupDatafileReader&        dfile)
{
  // First crude division
  std::vector<std::vector<unsigned int> >  tmp(2);
  for (unsigned int i=0; i<indx.size(); i++) {
    if ((dfile.PhaseEncodeVector(indx[i]))(1) && !(dfile.PhaseEncodeVector(indx[i]))(2)) tmp[0].push_back(i);
    else if (!(dfile.PhaseEncodeVector(indx[i]))(1) && (dfile.PhaseEncodeVector(indx[i]))(2)) tmp[1].push_back(i);
    else throw ApplyTopupException("Invalid phase-encode vector for least-squares restoration");
  }

  // Now ensure uniqueness
  bool unique[2] = { false, false };  
  for (unsigned int i=0; i<2; i++) {
    for (unsigned int j=1; j<tmp[i].size(); j++) {
      if (dfile.PhaseEncodeVector(indx[tmp[i][j]])(i+1)*dfile.ReadOutTime(indx[tmp[i][j]]) != dfile.PhaseEncodeVector(indx[tmp[i][0]])(i+1)*dfile.ReadOutTime(indx[tmp[i][0]])) {
        unique[i] = true;
      }
    }
  }
  if ((tmp[0].size() && !unique[0]) || (tmp[1].size() && !unique[1])) throw ApplyTopupException("Invalid combination of phase-encode vectors for least-squares restoration");

  // Finally repackage
  if (tmp[0].size() && tmp[1].size()) { scans.resize(2); indxs.resize(2); }
  else { scans.resize(1); indxs.resize(1); }
  for (unsigned int i=0, j=0; i<2; i++) {
    if (tmp[i].size()) scans[j++]=tmp[i];
  }
  for (unsigned int i=0, ii=0; i<2; i++) {
    if (tmp[i].size()) {
      scans[ii].resize(tmp[i].size());
      indxs[ii].resize(tmp[i].size());
      for (unsigned int j=0; j<tmp[i].size(); j++) {
        scans[ii][j] = tmp[i][j];
        indxs[ii][j] = indx[tmp[i][j]];
      }
      ii++;
    }
  }
}

unsigned int TopupCollections::IndexAt(unsigned int c, unsigned int i)
{
  if (c>=NCollections() || i>=NIndices(c)) throw ApplyTopupException("TopupCollections::IndexAt: Indices out of range");
  return(indxs[c][i]);
}

unsigned int TopupCollections::ScanAt(unsigned int c, unsigned int i)
{
  if (c>=NCollections() || i>=NIndices(c)) throw ApplyTopupException("TopupCollections::ScanAt: Indices out of range");
  return(scans[c][i]);
}

bool row_col_is_alright(const NEWIMAGE::volume<char>&   mask,
                        unsigned int                    k,
                        unsigned int                    ij,
                        bool                            row)
{
  if (row) { for (int i=0; i<mask.xsize(); i++) if (!mask(i,ij,k)) return(false); }
  else { for (int j=0; j<mask.ysize(); j++) if (!mask(ij,j,k)) return(false); }   
  return(true);
}

NEWMAT::ReturnMatrix extract_row_col(const NEWIMAGE::volume<float>&  vol,
                                     unsigned int                    k,
                                     unsigned int                    ij,
                                     bool                            row)
{
  NEWMAT::ColumnVector  v;
  if (row) v.ReSize(vol.xsize());
  else v.ReSize(vol.ysize());  
  if (row) { for (int i=0; i<vol.xsize(); i++) v(i+1) = vol(i,ij,k); }
  else { for (int j=0; j<vol.ysize(); j++) v(j+1) = vol(ij,j,k); }

  v.Release();
  return(v);  
}

void add_to_slice(NEWIMAGE::volume<float>&                  vol,
                  const std::vector<NEWMAT::ColumnVector>&  mu,
                  unsigned int                              sl,
                  bool                                      row)
{
  for (unsigned int i=0; i<mu.size(); i++) {
    for (int j=0; j<mu[0].Nrows(); j++) {
      if (row) vol(j,i,sl) += mu[i](j+1);
      else vol(i,j,sl) += mu[i](j+1);
    }
  }
}

void add_to_rows_cols(NEWIMAGE::volume4D<float>&  vol,
                      const NEWMAT::Matrix&       B,
                      unsigned int                k,
                      unsigned int                ij,
                      bool                        row)
{
  if (vol.tsize() != B.Ncols()) throw ApplyTopupException("add_to_rows_cols: Mismatch between vol and B");
  for (int d=0; d<vol.tsize(); d++) {
    if (row) { for (int i=0; i<vol.xsize(); i++) vol(i,ij,k,d) += B(i+1,d+1); }
    else { for (int j=0; j<vol.ysize(); j++) vol(ij,j,k,d) += B(j+1,d+1); }
  }  
}

void zero_out_rows_cols(NEWIMAGE::volume4D<float>&   vol,
                        const NEWMAT::Matrix&        map,
                        bool                         row)
{
  for (int d=0; d<vol.tsize(); d++) {
    for (int k=0; k<vol.zsize(); k++) {
      if (row) {
        for (int j=0; j<vol.ysize(); j++) {
          if (map(k+1,j+1)) {
	    for (int i=0; i<vol.xsize(); i++) vol(i,j,k,d) = 0.0;
  	  }
        }
      }
      else {
        for (int i=0; i<vol.xsize(); i++) {
	  if (map(k+1,i+1)) {
	    for (int j=0; j<vol.ysize(); j++) vol(i,j,k,d) = 0.0;
	  }
        }
      }
    }
  }
}

} // End namespace TOPUP
