/*   asl_file.cc file manipulator for multi-TI ASL data

      Michael Chappell - FMIRB Image Analysis Group

      Copyright (C) 2009 University of Oxford */

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

#include <iostream>
#include <math.h>
#include <string>
#include "newmatap.h"
#include "newmatio.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "utils/tracer_plus.h"
#include "stdlib.h"

#include "readoptions.h"
#include "asl_functions.h"

using namespace Utilities;
using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace MISCMATHS;

using namespace OXASL;

int main(int argc, char *argv[])
{
  try {

    //parse command line (puts these into the log file)
    ReadOptions& opts = ReadOptions::getInstance();
    opts.parse_command_line(argc,argv);

   //deal with input data type options
    bool isblocked=false; //indicates if data is in blocks of repeats rather than TIs
    bool ispairs=false; //indicates if data contains adjacent pairs of measurments
    bool isdiff=false; //indicates if we have differenced data
    bool tagfirst=true; //indicates that tag comes first in tag-control pairs
    //ispairs=opts.ispairs.value();
    
    // block format: isblocked indicates if data are in *blocks of repeats*
    string ibf=opts.inblockform.value();
    if      (ibf.compare("rpt")==0) isblocked=true;
    else if (ibf.compare("tis")==0) isblocked=false;
    else    throw Exception("Unrecognised input block format");
    
    string iaf=opts.inaslform.value();
    if      (iaf.compare("diff")==0) isdiff=true;
    else if (iaf.compare("tc")==0)   isdiff=false;
    else if (iaf.compare("ct")==0)  { isdiff=false; tagfirst=false;}
    else    throw Exception("Unrecognised input asl form");
    
    ispairs=!isdiff; //convienient to have ispairs as a bool

    bool outblocked;
    string obf=opts.outblockform.value();
    if      (obf.compare("rpt")==0) outblocked=true;
    else if (obf.compare("tis")==0) outblocked=false;
    else if (obf.compare("notset")==0) outblocked=isblocked; //default is to set output same as input
    else    throw Exception("Unrecognised output block format");

    /*
    bool outdiff;
    string oaf=opts.outaslform.value();
    if     (oaf.compare("diff")==0) outdiff=true;
    else if(oaf.compare("tc")==0) outdiff=false;
    else throw Exception("Unrecognised output asl form");
    */

    bool outpairs=ispairs; // outpairs indicates wehter the data we are processing for output is in the form of pairs - by deafult if input is in pairs then output the pairs

    //load data
    volume4D<float> data;
    read_volume4D(data,opts.datafile.value());

    // load mask
    // if a mask is not supplied then default to processing whole volume
    volume<float> mask(data.xsize(),data.ysize(),data.zsize());
    mask.setdims(data.xdim(),data.ydim(),data.zdim());
    mask=1;
    if (opts.maskfile.set()) read_volume(mask,opts.maskfile.value());

    Matrix datamtx;
    datamtx = data.matrix(mask);
    int nvox=datamtx.Ncols();

    cout << "Number of voxels is:" << nvox << endl;

    //establish number of repeats
    int ntis=opts.ntis.value();
    int nmeas=data.tsize()/ntis; //number of measurements at each TI
    int nrpts;
    if (ispairs) nrpts=nmeas/2; // number of repeats, note that is data is pairs then we have half as many true repeats
    else nrpts=nmeas;
    cout << "Number of repeats in data is:" << nrpts << endl;

    int ndata;
    if (ispairs) ndata=ntis*nrpts*2;
    else ndata = nmeas*ntis;
    if (ndata<data.tsize()) {
      //some leftover samples in the data, produce a warning
      cout << "Warning: spare measurements found at end of data!" << endl
	   << "    Number of measurements discarded from end is: " << data.tsize()-ndata << endl;

      // discard
      datamtx = datamtx.Rows(1,ndata);
    }


    //get data into 'standard' format
    vector<Matrix> asldata;
    data2stdform(datamtx,asldata,ntis,isblocked,ispairs);

    //deal with the splitting of pairs
    vector<Matrix> asldataodd;
    vector<Matrix> asldataeven;
    // int idx;
    
    if (opts.splitpairs.value() && opts.tcdiff.value()) {
	// doesn't make sense to try and do both splitpairs and tcdifference
	throw Exception("Cannot both split the pairs and difference them");
      }


  if (opts.splitpairs.value()) {
    // need to split the data here if plit pairs has been requested
      separatepairs(asldata,asldataodd,asldataeven);
    }

  //tag control difference
  if (opts.tcdiff.value()) {
    separatepairs(asldata,asldataodd,asldataeven); //split pairs ready for differencing
    //overwrite asldata with differenced data
    for (int ti=0; ti<ntis; ti++) {
      asldata[ti] = asldataeven[ti] - asldataodd[ti];
      if (!tagfirst) asldata[ti] *= -1.0; //if control image is first then the sign will be wrong here
    }
	outpairs=false;
  }

  // surround tag-control difference
  if (opts.surrtcdiff.value()) {
    //separatepairs(asldata,asldataodd,asldataeven); //split pairs ready for differencing
    //overwrite asldata with differenced data
    for (int ti=0; ti<ntis; ti++) {
      Matrix tempindata;
      tempindata=asldata[ti];
      Matrix tempdata(tempindata.Nrows()-1,tempindata.Ncols());
      for (int aq=1; aq<tempindata.Nrows(); aq++) {
	tempdata.Row(aq) = tempindata.Row(aq) - tempindata.Row(aq+1);
      }
      asldata[ti] = tempdata;
      if (!tagfirst) asldata[ti] *= -1.0; //if control image is first then the sign will be wrong here
    }
	outpairs=false;
  }


    vector<Matrix> asldataout; //the data we are about to use for output purposes
    int nout=1; //number of output cycles we need to go through
    string fsub="";
    
    if (opts.splitpairs.value()) {
      nout=2;
      outpairs=false;
    }

    // OUTPUT: Main output section - anything that is compatible with the splitting of pairs (splitpairs) can go in here
    for (int o=1; o<=nout; o++) {
      if(!opts.splitpairs.value()) asldataout=asldata;
      else if (o==1) { asldataout=asldataodd; fsub="_odd"; cout << "Dealing with odd members of pairs"<<endl;}
      else           { asldataout=asldataeven; fsub="_even"; cout << "Dealing with even members of pairs"<<endl;}

      //output data
      if (opts.out.set()) {
	Matrix outmtx;
	volume4D<float> dataout;
	stdform2data(asldataout,outmtx,outblocked,outpairs);
	dataout.setmatrix(outmtx,mask);
	save_volume4D(dataout,opts.out.value()+fsub);
      }

      //take mean at each TI
      if (opts.meanout.set()) {
	timeanout(asldataout,mask,opts.meanout.value()+fsub,outpairs);
	/*cout << "Outputting ASL data mean at each TI" << endl;
	  Matrix meanti(ntis,nvox);
	  for (int n=0; n<ntis; n++) 
	  {
	  meanti.Row(n+1)=mean(asldata[n],1);
	  }
	  volume4D<float> meanout;
	  meanout.setmatrix(meanti,mask);
	  save_volume4D(meanout,opts.meanout.value());*/
      }
    
      //split data into separate file for each TI
      if (opts.splitout.set()) {
	splitout(asldataout,mask,opts.splitout.value()+fsub);
	/*cout << "Splitting ASL data into files for each TI" << endl;
	  volume4D<float> blockout;
	  for (int n=0; n<ntis; n++) 
	  {
	  char cstr [5];
	  if (n<10) sprintf(cstr,"00%d",n);
	  else if (n<100) sprintf(cstr,"0%d",n);
	  else if (n<1000) sprintf(cstr,"%d",n);
	  else throw Exception("More than 1000 measurements in this ASL data file, sorry cannot handle this operation");
	  string tino(cstr);
	  
	  blockout.setmatrix(asldata[n],mask);
	  save_volume4D(blockout,opts.splitout.value()+tino);
	  }*/
      }

    //do epochwise output
    if (opts.epochout.set()) {
      int eplen=opts.epochlen.value();
      int epol=opts.epochover.value();

      bool tiunit; //indicates if we want epochs over TIs rather than over repeats
      string epunit=opts.epochunit.value();
      if      (epunit.compare("rpt")==0) tiunit=false;
      else if (epunit.compare("tis")==0) tiunit=true;
      else    throw Exception("Unrecognised epoch unit");


      //if (ispairs && !opts.tcdiff.value() ) throw Exception("Epoch ceration is currently incompatible with data that contains pairs if differencing has not been performed");
      /* NOT CORRECT - the pairs are stores within the TI
      if (outpairs) {
	//if we are dealing with pairs then we have two measurements for every repeat
	// epoch parameters are (meant to be) in number of repeats
	// hence multiply by 2 in this scenario
	eplen*=2;
	epol*=2;
	}*/

      epochout(asldataout,mask,opts.epochout.value()+fsub,eplen,epol,outpairs,tiunit);
      /*cout << "Creating ASL data epochs" << endl;
      int eplen = opts.epochlen.value();
      int epol = opts.epochover.value();
      if (epol>=eplen) throw Exception("The epoch overlap may not exceed or equal the length of the epoch");
      int epadv = eplen-epol;

      //NOTE that ispairs indicates that in asldata the measurements come in pairs that we collected at the same time
      //int ncollect=1; //number of measurements to collect for a single repeat-TI
      //if (ispairs) ncollect=2;

      int e=1;
      Matrix epoch_temp(ntis,nvox);
      volume4D<float> epochout;
      while (e*epadv<=nmeas)
	{
	  //Matrix ti_temp(eplen,nvox);
	  //go through the TIs
	  for (int ti=1; ti<=ntis; ti++) {
	    epoch_temp.Row(ti) = mean(asldata[ti-1].Rows((e-1)*epadv+1,e*epadv),1);
	  }

	  char cstr [5];
	  if (e<10) sprintf(cstr,"00%d",e);
	  else if (e<100) sprintf(cstr,"0%d",e);
	  else if (e<1000) sprintf(cstr,"%d",e);
	  else throw Exception("More than 1000 epochs in this ASL data file, sorry cannot handle this operation");
	  string epno(cstr);
	  
	  epochout.setmatrix(epoch_temp,mask);
	  save_volume4D(epochout,opts.epochout.value()+epno);

	  e++;
	}

      int unused=nmeas-(e-1)*epadv;
      cout << unused << endl;
      assert(unused>=0);
      if (unused>0) cout << "Number of measurements from end of data discarded: " << unused << endl;*/
      }
  }

// OUTPUT: supplementary output section, things that are not compatible with the splitting of pairs into separate outputs

    // do deconvolution
    if (opts.deconvout.set()) {
     
      //load aif
      volume4D<float> aif;
      read_volume4D(aif,opts.aif.value());
      Matrix aifmtx;
      aifmtx = aif.matrix(mask);

      deconvout(asldata,mask,aifmtx,opts.deconvout.value());
    }
  }

catch(Exception& e)
    {
      cerr << endl << e.what() << endl;
      return 1;
    }
  catch(X_OptionError& e)
    {
      cerr << endl << e.what() << endl;
      return 1;
    }

  cout << "Done." << endl;

  return 0;


}
