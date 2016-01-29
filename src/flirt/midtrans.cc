/*  midtrans.cc

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2003 University of Oxford  */

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

#define _GNU_SOURCE 1
#define POSIX_SOURCE 1

#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/options.h"
#include <vector>

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;
using namespace std;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="midtrans \nCopyright(c) 2010, University of Oxford (Mark Jenkinson)";
string examples="midtrans [options] transform1 transform2 ... transformN\n  e.g. midtrans -o temp2mid.mat A2temp.mat B2temp.mat C2temp.mat\n       midtrans -o C2mid.mat A2C.mat B2C.mat ident.mat";

// Each (global) object below specificies as option and can be accessed
//  anywhere in this file (since they are global).  The order of the
//  arguments needed is: name(s) of option, default value, help message,
//       whether it is compulsory, whether it requires arguments
// Note that they must also be included in the main() function or they
//  will not be active.

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<bool> debug(string("--debug"), false, 
		     string("switch on debugging output"), 
		     false, no_argument);
Option<string> templte(string("--template"), string(""),
		  string("input filename for template image (needed for fix origin)"),
		  false, requires_argument);
Option<string> separateouts(string("--separate"), string(""),
		  string("basename for output of separate matrices (final name includes a number; e.g. img2mid0001.mat)"),
		  false, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		  string("output filename for matrix"),
		  false, requires_argument);
int nonoptarg;

////////////////////////////////////////////////////////////////////////////

// Local functions

Matrix approx_nth_root(const vector<Matrix>& matlist) 
{
  int NN = matlist.size();
  // initialise using the affine decomposition (arthimetic average + geometric mean)
  ColumnVector centre(3), params(12), totaltrans(3);
  centre=0.0;

  Matrix rootmat(3,3), prodmat(4,4), prodmat3(3,3);
  prodmat = IdentityMatrix(4);
  for (int n=0; n<NN; n++) {
    prodmat *= matlist[n];
  }
  // store and then zero the translational component
  totaltrans(1)=prodmat(1,4);  prodmat(1,4)=0.0;
  totaltrans(2)=prodmat(2,4);  prodmat(2,4)=0.0;
  totaltrans(3)=prodmat(3,4);  prodmat(3,4)=0.0;
  prodmat3 = prodmat.SubMatrix(1,3,1,3);
  if (verbose.value()) { cout << "Product matrix = " << endl << prodmat << endl; }

  // calculate the affine decomposition and approximate root
  decompose_aff(params, prodmat, centre, rotmat2euler);
  if (verbose.value()) { cout << "Affine params = " << params.t() << endl; }
  for (int n=1; n<=12; n++) {
    if ((n>=7) && (n<=9)) { 
      params(n) = exp(log(params(n))/NN);   // NN-th root
    } else {
      params(n) /= NN;
    }
  }
  compose_aff(params, 12, centre, rootmat, construct_rotmat_euler);

  // iterative solution process (no convergence guarantees though, so just hoping for the best!)
  Matrix err, apow, amat, bmat;
  float errval=1.0, errthresh=1e-8;
  int nit=0, max_it=20;
  amat = rootmat.SubMatrix(1,3,1,3);
  while ((errval>errthresh) && (++nit<max_it)) {
    if (verbose.value()) { cout << "Iteration number " << nit << endl; }
    apow = IdentityMatrix(3);
    for (int n=1; n<=NN; n++) { apow *= amat; }
    err = apow - prodmat3;
    errval = err.SumSquare()/9.0;
    bmat = (err*apow.i())*(-1.0/NN);
    amat += bmat;
  }
  if (errval>errthresh) {
    // Put alternative method (e.g. direct numerical optimisation) in here
    cerr << "Error:: failed to converge to tolerance" << endl;
    cerr << "Approximation error = " << errval << endl;
    exit(EXIT_FAILURE);
  }

  // calculate A^(n-1) + A^(n-2) + ... + A + I
  Matrix matsum;
  ColumnVector trans;
  matsum=IdentityMatrix(3);
  for (int n=1; n<=NN-1; n++) {
    matsum = amat*matsum + IdentityMatrix(3);
  }
  trans = matsum.i() * totaltrans;

  // put solution (amat) + translation back into rootmat
  rootmat.SubMatrix(1,3,1,3) = amat;
  rootmat(1,4)=trans(1);
  rootmat(2,4)=trans(2);
  rootmat(3,4)=trans(3);

  Matrix midtrans(4,4);
  midtrans=rootmat.i();  
return midtrans;
}


Matrix avg_affine_with_procrustes(const vector<Matrix>& matlist, const ColumnVector& temp_orig,
				  float rmax)
{
  // Based on minimisation of squared distance after factoring out average and individual poses
  // Uses min E = -Trace[A.P] is solved for a 6 dof A (orthogonal) by A=V.U' where P=U.D.V' (SVD)
  int NN = matlist.size();
  Matrix M44;
  vector<ColumnVector> xavlist;
  ColumnVector v0001(4), x_av(3);

  if (verbose.value()) { cout << "Initialising averages" << endl; }

  // Initially calculate the centre voxel locations for each native image by projecting template's centre
  v0001 << temp_orig(1) << temp_orig(2) << temp_orig(3) << 1; // Centre/Origin/COV of the template (in FLIRT coords)
  // calculate centre of individual space (taking template "centre" as v0001)
  for (int n=0; n<NN; n++) {
    M44 = matlist[n];
    x_av = M44.i() * v0001;
    x_av = x_av.Rows(1,3);
    xavlist.push_back(x_av);
  }

  // STAGE 1 (estimate best rigid body average correction for template pose)
  // Equations are: min E = -Trace[A(<P.x.x'> + <t.x'> - <P.x><x'> - <t><x'>)]
  //           and: s = <x> - A(<Px> + <t>)
  // Solve for A (3x3) and s (3x1) given sets of P (3x3) and t (3x1)
  // <P.x.x'> = (1/N) * \sum_i P_i <x_i.x_i'> 
  // <x_i.x_i'> = (1/5)*(r^2)*I + <x_i>.<x_i>'
  if (verbose.value()) { cout << "Stage 1" << endl; }
  ColumnVector s(3), x_i(3), t_i(3), t_av(3), Px_av(3);
  Matrix A(3,3), P_i(3,3), tx_av(3,3), Pxx_av(3,3), M2, U, V;
  DiagonalMatrix D;

  s=0; x_i=0; x_av=0; t_i=0; t_av=0; Px_av=0;
  A=0; P_i=0; tx_av=0; Pxx_av=0;

  for (int n=0; n<NN; n++) { 
    x_i = xavlist[n];
    P_i = matlist[n].SubMatrix(1,3,1,3);
    t_i = matlist[n].SubMatrix(1,3,4,4);
    x_av += x_i;
    t_av += t_i;
    tx_av += t_i*x_i.t();
    Px_av += P_i * x_i;
    Pxx_av += P_i * ((rmax*rmax/5)*IdentityMatrix(3) + x_i*x_i.t());
    if (debug.value()) { cout << "P_"<< n+1 << " = " << endl << P_i << endl; }
  }
  x_av /= (float) NN;
  t_av /= (float) NN;
  tx_av /= (float) NN;
  Px_av /= (float) NN;
  Pxx_av /= (float) NN;
  
  if (debug.value()) { cout << "Pxx_av = " << endl << Pxx_av << endl; }

  M2 = Pxx_av + tx_av - Px_av*x_av.t() - t_av*x_av.t();
  SVD(M2,D,U,V);
  A = V*U.t();
  s = x_av - A*(Px_av + t_av);
  if (verbose.value()) { cout << "A = " << endl << A << endl << "s = " << s.t() << endl; }

  // STAGE 2 (estimate a separate rigid body approx to the template for each image)
  // Equations are: min E = -Trace[R(<xx'>-<x><x'>)P'] = -Trace[R.P']   (as <xx'>-<x><x'>=0.2*r^2*I)
  //           and: u = (P-R)<x> + t
  // Solve for R (3x3) and u (3x1) given sets of P (3x3) and t (3x1), where P="AP" from above (4x4 sense)
  if (verbose.value()) { cout << "Stage 2" << endl; }
  vector<Matrix> matlist2;
  Matrix R, P;
  ColumnVector u, t;
  for (int n=0; n<NN; n++) {
    M44 = matlist[n];
    P = A * M44.SubMatrix(1,3,1,3);
    t = A * M44.SubMatrix(1,3,4,4) + s;
    SVD(P.t(), D, U, V);
    R = V*U.t();
    u = (P-R)*xavlist[n] + t;
    M44.SubMatrix(1,3,1,3) = R;
    M44.SubMatrix(1,3,4,4) = u;
    matlist2.push_back(M44);
    if (verbose.value()) { cout << "Subject number " << n << endl; }
    if (verbose.value()) { cout << "R = " << endl << R << endl << "u = " << u.t() << endl; }
  }

  // STAGE 3 (estimate the best affine average correction to the template using the best inidividual rigid body estimates from stage 2)
  // Equations are: M = -<ba'><aa'>^(-1)  ;   w = -M<Px> - M<t> + <Rx> + <u>
  //  where a_i = P_i.x_i-<Px> + t_i - <t> ; b_i = -(R_i.x_i -<Rx>) - (u_i-<u>)
  //  and P_i|t_i stands for "A.P_i" (4x4 sense) from after stage 1, and R_i|u_i is the result of stage 2
  // <aa'> = <Pxx'P'> - <Px><Px>' + <Pxt'> + <tx'P'> - <t><Px>' - <Px><t>' + <tt'> - <t><t>'
  // <ba'> = <Rxx'P'> - <Rx><Px>' + <Rxt'> - <Rx><t'> - <ux'P'> + <u><Px>' - <ut'> + <u><t>'
  // <Pxx'P'> = (1/5)*(r^2)*<PP'> + (1/N)*\sum_i P_i.<x_i>.<x_i>'.P_i'
  // <Rxx'P'> = (1/5)*(r^2)*<RP'> + (1/N)*\sum_i R_i.<x_i>.<x_i>'.P_i'
  
  if (verbose.value()) { cout << "Stage 3" << endl; }
  ColumnVector u_i(3), u_av(3), Rx_av(3);
  Matrix tt_av(3,3), ut_av(3,3), R_i(3,3), Pxt_av(3,3), Rxt_av(3,3), uxP_av(3,3), PP_av(3,3), RP_av(3,3), PxxP(3,3), RxxP(3,3);

  u_i=0; x_i=0; x_av=0; u_av=0; t_i=0; t_av=0; Px_av=0; Rx_av=0;
  P_i=0; tt_av=0; ut_av=0; R_i=0; Pxt_av=0; Rxt_av=0; uxP_av=0; PP_av=0; RP_av=0; PxxP=0; RxxP=0;

  for (int n=0; n<NN; n++) { 
    u_i = matlist2[n].SubMatrix(1,3,4,4);
    x_i = xavlist[n];
    x_av += x_i;
    u_av += u_i; 
    P_i = A * matlist[n].SubMatrix(1,3,1,3);
    t_i = A * matlist[n].SubMatrix(1,3,4,4) + s;
    t_av += t_i;
    tt_av += t_i*t_i.t();
    ut_av += u_i*t_i.t();
    Px_av += P_i * x_i;
    R_i = matlist2[n].SubMatrix(1,3,1,3);
    Rx_av += R_i * x_i;
    Pxt_av += P_i * x_i * t_i.t();
    Rxt_av += R_i * x_i * t_i.t();
    uxP_av += u_i * x_i.t() * P_i.t();
    PP_av += P_i * P_i.t();
    RP_av += R_i * P_i.t();
    PxxP += P_i * x_i * x_i.t() * P_i.t();
    RxxP += R_i * x_i * x_i.t() * P_i.t();
  }
  x_av /= (float) NN;
  u_av /= (float) NN;
  t_av /= (float) NN;
  tt_av /= (float) NN;
  ut_av /= (float) NN;
  Px_av /= (float) NN;
  Rx_av /= (float) NN;
  Pxt_av /= (float) NN;
  Rxt_av /= (float) NN;
  uxP_av /= (float) NN;
  PP_av /= (float) NN;
  RP_av /= (float) NN;
  PxxP = (rmax*rmax/5)*PP_av + PxxP/((float) NN);
  RxxP = (rmax*rmax/5)*RP_av + RxxP/((float) NN);

  Matrix aat, ba, M;
  aat = PxxP - Px_av*Px_av.t() + Pxt_av + Pxt_av.t() - t_av*Px_av.t() - Px_av*t_av.t() + tt_av - t_av*t_av.t();
  ba = -1*(RxxP - Rx_av*Px_av.t() + Rxt_av - Rx_av*t_av.t() - uxP_av + u_av*Px_av.t() - ut_av + u_av*t_av.t());

  ColumnVector w;
  M = -ba*(aat).i();
  w = -M*Px_av - M*t_av + Rx_av + u_av;
  if (verbose.value()) { cout << "M = " << endl << M << endl << "w = " << w.t() << endl; }

  if (verbose.value()) { cout << "Final calculations" << endl; }
  Matrix midtrans;
  // midtrans = M*A (as this is how the template is transformed)
  midtrans=IdentityMatrix(4);
  midtrans.SubMatrix(1,3,1,3) = M*A;
  midtrans.SubMatrix(1,3,4,4) = w + M*s;

  // Resulting matrix to get existing Image -> Old Template transform to Image -> New Template is:
  //   mat(Old Template -> New Template) * mat(Image -> Old Template)
  // That is, pre-multiply old matrices by current midtrans

  return midtrans;
}



int do_work(int argc, char* argv[]) 
{
  //int ntrans = argc - nonoptarg;

  // read in matrices
  vector<Matrix> matlist;
  for (int n=nonoptarg; n<argc; n++) {
    Matrix trans(4,4);
    trans = read_ascii_matrix(string(argv[n]));
    if (fabs(trans.Determinant()) < 1e-6) {
      cerr << "Could not read matrix " << argv[n] << endl;
      exit(EXIT_FAILURE);
    } 
    matlist.push_back(trans);
  }

  int NN = matlist.size();
  if (verbose.value()) { cout << "Number of matrices = " << NN << endl; }

  ColumnVector temp_orig(4);
  if (templte.set()) {
    volume<float> templ;
    read_volume(templ,templte.value());
    temp_orig = templ.cog("scaled_mm");
  } else {
    cout << "Assuming MNI origin as no template image specified" << endl;
    temp_orig(1)=90; temp_orig(2)=110; temp_orig(3)=90;  // from MNI
    if (debug.value()) { temp_orig(1)=0; temp_orig(2)=0; temp_orig(3)=0; }  // DEBUG ONLY!
  }
  if (verbose.value()) { cout << "Centre of template = " << temp_orig.t() << endl; }

  // Do the calculation
  Matrix midtrans;
  //midtrans = approx_nth_root(matlist);
  float rmax=80;  // size, in mm, of ROI to integrate error term over (template dependent - e.g. rats?)
  midtrans = avg_affine_with_procrustes(matlist,temp_orig,rmax);

  // show/save output matrix
  if (outname.set()) {
    write_ascii_matrix(midtrans,outname.value());
  }
  if (verbose.value() || outname.unset()) {
    cout << midtrans << endl;
  }

  if (separateouts.set()) {
    for (int n=1; n<=NN; n++) {
      write_ascii_matrix(midtrans * matlist[n-1],separateouts.value()+num2str(n,4)+".mat");
    }
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(outname);
    options.add(templte);
    options.add(separateouts);
    options.add(debug);
    options.add(verbose);
    options.add(help);
    
    nonoptarg = options.parse_command_line(argc, argv, 0, true);

    // line below stops the program if the help was requested or 
    //  a compulsory option was not set
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
    
    if ((argc - nonoptarg)<2) {
      cerr << "Must specify at least 2 transforms to find the mid-transform"
  	<< endl;
      options.usage();
      exit(EXIT_FAILURE);
    }

    // Call the local functions
    
    return do_work(argc,argv);
    
  } catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } catch(Exception &e) {
    cerr << e.what() << endl;
  } catch (...) {
    cerr << "Fatal error" << endl;
  } 
}

