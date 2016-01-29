/*    Copyright (C) 2012 University of Oxford  */

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
#include "utils/options.h"
#include "warpfns/warpfns.h"
#include "warpfns/fnirt_file_reader.h"
#include "newimage/newimageall.h"
#include "csv.h"
#include "stdlib.h"
#include "string.h"
#include "miscmaths/miscmaths.h"

using namespace Utilities;
using namespace std;
using namespace NEWIMAGE;
using namespace MISCMATHS;


string title="surf2surf - conversions between surface formats and/or conventions";
string examples="Usage: surf2surf -i <inputSurface> -o <outputSurface> [options]";


Option<string> surfin(string("-i,--surfin"),string(""),
		      string("input surface"),
		      true,requires_argument);
Option<string> surfout(string("-o,--surfout"),string(""),
		       string("output surface"),
		       true,requires_argument);
Option<string> convin(string("--convin"),string("caret"),
		      string("input convention [default=caret] - only used if output convention is different"),
		      false,requires_argument);
Option<string> convout(string("--convout"),string("caret"),
		       string("output convention [default=same as input]"),
		       false,requires_argument);
Option<string> volin(string("--volin"),string(""),
		     string("\tinput ref volume - Must set this if changing conventions"),
		     false,requires_argument);
Option<string> volout(string("--volout"),string(""),
		      string("output ref volume [default=same as input]"),
		      false,requires_argument);
Option<string> xfm(string("--xfm"),"",
		   string("\tin-to-out ascii matrix or out-to-in warpfield [default=identity]"),
		   false,requires_argument);
Option<string> otype(string("--outputtype"),"GIFTI_BIN_GZ",
		     string("output type: ASCII, VTK, GIFTI_ASCII, GIFTI_BIN, GIFTI_BIN_GZ (default)"),
		     false,requires_argument);
Option<string> vals(string("--values"),"",
		     string("set output scalar values (e.g. --values=mysurface.func.gii or --values=1)"),
		     false,requires_argument);




int main(int argc,char *argv[]){

  OptionParser options(title,examples);
  
  
  options.add(surfin);
  options.add(surfout);
  options.add(convin);
  options.add(convout);
  options.add(volin);
  options.add(volout);
  options.add(xfm);
  options.add(otype);
  options.add(vals);

  
  options.parse_command_line(argc,argv);
  if(!options.check_compulsory_arguments(true)){
    options.usage();
    return(1);
  }
  
  // check options
  if(convin.value()!=convout.value()){
    if(volin.value()==""){
      cerr<<"Please specify input reference volume"<<endl;
      return(1);
    }
  }

  volume<short int> refvolin,refvolout;
  CSV csv;
  if(volin.set()){
    read_volume(refvolin,volin.value());
    csv.reinitialize(refvolin);
    csv.set_convention(convin.value());
  }
  if(volout.set()){
    read_volume(refvolout,volout.value());
  }
  else{
    refvolout.reinitialize(refvolin.xsize(),refvolin.ysize(),refvolin.zsize());
    refvolout=refvolin;
  }

  csv.load_rois(surfin.value());
  //csv.get_mesh(0).print("grot.txt");

  if(convin.value()!=convout.value()){  
    bool isWarp=false;
    volume4D<float> vox2vox_warp;
    Matrix          vox2vox;
    ColumnVector old_dims(3);
    old_dims << refvolin.xdim() << refvolin.ydim() << refvolin.zdim();
    if(xfm.set()){
      if(fsl_imageexists(xfm.value())){
	isWarp=true;
	FnirtFileReader ffr(xfm.value());
	vox2vox_warp=ffr.FieldAsNewimageVolume4D(true);
      }
      else{
	vox2vox=read_ascii_matrix(xfm.value());
      }
    }
    else
      vox2vox=IdentityMatrix(4);
    
    csv.set_refvol(refvolout);
    if(convout.set()){
      if(!isWarp)
	csv.switch_convention(convout.value(),vox2vox,old_dims);
      else
	csv.switch_convention(convout.value(),vox2vox_warp,refvolin,refvolout);
    }
    else{
      if(!isWarp)
	csv.switch_convention(convin.value(),vox2vox,old_dims);
      else
	csv.switch_convention(convin.value(),vox2vox_warp,refvolin,refvolout);
    }
  }

  if(vals.value()!=""){
    if(meshExists(vals.value())){
      CsvMesh m;
      m.load(vals.value());
      //OUT(m._pvalues.size());
      //OUT(csv.get_mesh(0).nvertices());
      if((int)m._pvalues.size()!=csv.get_mesh(0).nvertices()){
	cerr<<"values mesh does not have the correct number of vertices"<<endl;
	exit(1);
      }
      vector<float> v = m.getValuesAsVectors();
      ColumnVector vv(v.size());
      for(unsigned int i=0;i<v.size();i++)
	vv(i+1)=v[i];
      csv.get_mesh(0).set_pvalues(vv);
    }
    else{
      csv.get_mesh(0).set_pvalues((float)atoi(vals.value().c_str()));
    }
  }

  //csv.save_roi(0,surfout.value());
  int o;
  if(otype.value()=="ASCII")
    o=CSV_ASCII;
  else if(otype.value()=="VTK")
    o=CSV_VTK;
  else if(otype.value()=="GIFTI_ASCII")
    o=GIFTI_ENCODING_ASCII;
  else if(otype.value()=="GIFTI_BIN")
    o=GIFTI_ENCODING_B64BIN;
  else if(otype.value()=="GIFTI_BIN_GZ")
    o=GIFTI_ENCODING_B64GZ;
  else{
    cerr<<"Unknown format "<<otype.value()<<endl;
    exit(1);
  }
  csv.get_mesh(0).save(surfout.value(),o);


  return 0;
}



