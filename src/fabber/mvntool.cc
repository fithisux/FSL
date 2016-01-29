/* mvntool.cc - Tool for adding parameters to an MVN

   Michael Chappell, FMRIB Analysis Group

   Copyright (C) 2007 University of Oxford  */
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
using namespace std;

#ifdef __FABBER_LIBRARYONLY // Skip entire file if making fabber_library
int main() {cout << "MVNTOOL not built; compiled with __FABBER_LIBRARYONLY option." << endl; return 2;}
#else

#include <exception>
#include <stdexcept>
#include <map>
#include <string>
#include "dist_mvn.h"
#include "easyoptions.h"
#include "newimage/newimageall.h"

using namespace Utilities;
using namespace MISCMATHS;
using namespace NEWIMAGE;

/* Function declarations */
void Usage(const string& errorString = "");

int main(int argc, char** argv)
{
	try
	  {
	    cout << "FABBER: MVNtool" << endl;

	    EasyOptions args(argc, argv);

	    if (args.ReadBool("help"))
	      {
		Usage();
		return 0;
	      }

	    EasyLog::StartLogUsingStream(cout);

	    /* parse command line arguments*/
	    bool verbose=args.ReadBool("v");

	    string infile;
	    string outfile;
	    infile = args.Read("input");
	    string maskfile;
	    maskfile = args.Read("mask");
	    int param=0;
	    int cparam=0;
	    outfile = args.ReadWithDefault("output",infile);

	    bool ins; bool write;

	    double val;	double var;
	    string valimfile; string varimfile;
	    bool bval=false;bool bvar=false;bool cvar=false;
	    /* Choose what we want to do to/with the parameter - default is to read */
	    ins = args.ReadBool("new"); //insert a new parameter
	    write = args.ReadBool("write"); //overwrite an existing parameter

// determine which parameter we are reading/writing etc
	    string plistfile = args.ReadWithDefault("param-list","");
	    if (plistfile=="")
	      { //use must have specified a parameter number
		param = convertTo<int>(args.Read("param"));

		//deal with --cvar option
		cparam = convertTo<int>(args.ReadWithDefault("cvar","0"));
	      }
	    else
	      {
		string paramname = args.Read("param");
		// read in parameter names
		ifstream paramFile((plistfile).c_str());
		vector<string> paramNames;

		//NOTE covariance parameter not setup for parameter lists
	       
		string currparam;
		while (paramFile.good())
		  {
		    getline(paramFile,currparam);
		    paramNames.push_back(currparam);
		  }
		paramNames.pop_back(); //remove final empty line assocaited with eof
		
		//user will have specified a parameter name in the list
		string nplistfile = args.ReadWithDefault("new-param-list","");
		if (nplistfile=="")
		  { //simply deal with the parameter list
		    for (unsigned int i=0; i<paramNames.size(); i++)  
		      {
			if (paramname.compare(paramNames[i])==0) param=i+1;
		      }
		    if (param==0) throw Exception("Cannot find specfied parameter name in list");

		  }
		else
		  { //deal with current and new parameter lists
		    // if we are here we must be inserting a new parameter
		    ins=true;

		    //load new parameter list
		    ifstream newparamFile((nplistfile).c_str());
		    vector<string> newparamNames;
		    while (newparamFile.good())
		      {
			getline(newparamFile,currparam);
			newparamNames.push_back(currparam);
		      }
		    newparamNames.pop_back(); //remove final empty line assocaited with eof

		    // check parameter name is not in old list
		    for (unsigned int i=0; i<paramNames.size(); i++)  
		      {
			if (paramname.compare(paramNames[i])==0) throw Exception("Parameter name found in parameter list for this MVN, cannot insert an identical parameter");
		      }
		    //find parameter in new list
		    int newparam=0;
		    for (unsigned int i=0; i<newparamNames.size(); i++)  
		      {
			if (paramname.compare(newparamNames[i])==0) newparam=i+1;
		      }
		    if (newparam==0) throw Exception("Cannot find specfied parameter name in new parameter name list");
		    //cout << newparam<<endl;

		    if (newparam==1)
		      {
			param=1;
		      }
		    else
		      {
			//get name of previous parameter in list
			string preparamname = newparamNames[newparam-1-1];
			//cout << preparamname << endl;
			// find this in the old list
			for (unsigned int i=0; i<paramNames.size(); i++)  
			  {
			    if (preparamname.compare(paramNames[i])==0) param=i+1;
			  }
			if (param==0) throw Exception("Cannot complete this operation since the new list contains other parameters not present in the old list.\n You may be able to get around this by adding new parameters in a different order:\n new parameters that are sequential in the new list should be added starting with the earliest appearing.");
			// param currently hold location of previous parameter to the one we want to insert, so increment
			param++;
			//cout << param <<endl;
		      }

		    //need to write out the updated parameter list
		    vector<string> outParams;
		    for (int i=0; i<param-1; i++) { outParams.push_back(paramNames[i]); }
		    outParams.push_back(paramname);
		    for (unsigned int i=param-1; i<paramNames.size(); i++) { outParams.push_back(paramNames[i]); }

		    string paramlistout=args.ReadWithDefault("out-param-file","");
		    ofstream outparamFile((paramlistout).c_str());
		    if (paramlistout=="") {
		      for (unsigned int i=0; i < outParams.size(); i++)
			{
			  outparamFile << outParams[i] << endl;
			}
		    }
		  }
	      }


	    // Deal with --val and --var
	    if (ins | write) {
	      if (ins & write) { throw Invalid_option("Cannot insert and write at same time - choose either --new or --write"); }
		

		valimfile = args.ReadWithDefault("valim","");
		varimfile = args.ReadWithDefault("varim","");

		val = convertTo<double>(args.ReadWithDefault("val","-1e-6"));
		var = convertTo<double>(args.ReadWithDefault("var","-1e-6"));
	      }
	    else {
	      if (outfile==infile) throw Invalid_option("Output filename has not been specified");

	      bval = args.ReadBool("val");
	      bvar = args.ReadBool("var");
	      if (cparam>0) cvar=true;

	      if (bval & bvar) { throw Invalid_option("Cannot output value and variance at same time - choose either --val or --val"); }
	      if (bval & cvar) { throw Invalid_option("Cannot output value and covariance at same time"); }
	      if (bvar & cvar) { throw Invalid_option("Cannot output variance and covariance at same time"); }
	      if (!bval & !bvar & !cvar) { throw Invalid_option("Please select whether you want to extract the value (--val) or the variance (--var) or a covaraince (--cvar)"); }
	    }



	    /* Read in MVN file */
	    volume<float> mask;
            read_volume(mask,maskfile);
	    mask.binarise(1e-16,mask.max()+1,exclusive);

	    if (verbose) cout << "Read file" << endl;
	    vector<MVNDist*> vmvnin;
	    MVNDist::Load(vmvnin,infile,mask);

	    
	    if (ins | write) {
		/* section to deal with writing to or inserting into an MVN*/

	      // deal with the values we are going to insert
	      ColumnVector inmean(vmvnin.size());
	      ColumnVector invar(vmvnin.size());
	      if (valimfile=="") {
		inmean=val;
	      }
	      else {
		volume4D<float> valimvol;
		read_volume4D(valimvol,valimfile);
		inmean = (valimvol.matrix(mask)).AsColumn();
	      }
	      if (varimfile=="") {
		invar=var;
	      }
	      else {
		volume4D<float> varimvol;
		read_volume4D(varimvol,varimfile);
		invar = (varimvol.matrix(mask)).AsColumn();
	      }


		vector<MVNDist*> vmvnout(vmvnin);

		int oldsize;
		MVNDist mvnin;
		mvnin = *vmvnin[1];
		oldsize = mvnin.GetSize();

		if (ins)
		  { 
		    if (verbose) cout << "Inserting new parameter" << endl;
		    if (param > oldsize+1) throw Invalid_option("Cannot insert parameter here, not enough parameters in existing MVN");
		  }
		else
		  { if (param > oldsize) throw Invalid_option("Cannot edit this parameter, not enough parameters in existing MVN");}
		SymmetricMatrix mvncov;
		MVNDist mvnnew(1);
		mvnnew.means = 0;
		SymmetricMatrix covnew(1);
		covnew(1,1) = 0;
		mvnnew.SetCovariance(covnew);
		MVNDist mvnout;
		
		/* Loop over each enrty in mvnin - each voxel! */
		for (unsigned v=0;v < vmvnin.size(); v++)
		  {
		    mvnin = *vmvnin[v];
		    
		    if (ins)
		      { /* insert new parameter */
			//cout << "Add new parameter" << endl;
			
			if (param-1 >= 1)
			  {
			    MVNDist mvn1(param -1);
			    mvn1.CopyFromSubmatrix(mvnin,1,param -1,0);
			    mvnout=MVNDist(mvn1,mvnnew);
			  }
			else
			  {
			    mvnout = mvnnew;
			  }
			
			if (oldsize >= param)
			  {
			    MVNDist mvn2(oldsize+1-param);
			    mvn2.CopyFromSubmatrix(mvnin,param,oldsize,0);
			    mvnout=MVNDist(mvnout,mvn2);
			  }
   
			mvncov = mvnout.GetCovariance();
		      }
		    else
		      { 
			mvnout = mvnin;
			mvncov = mvnin.GetCovariance();
		      }
		    
		    /* Set parameters mean value and variance*/
		    //cout << "Set parameter mean and variance" << endl;
		    mvnout.means(param) = inmean(v+1);
		    mvncov(param,param) = invar(v+1);
		    mvnout.SetCovariance(mvncov);
		    
		    vmvnout[v] = new MVNDist(mvnout);
		  }

		/* Save MVN to output */
		if (verbose) cout << "Save file" << endl;
		MVNDist::Save(vmvnout,outfile,mask);
	      }

	    else {
	      /* Section to deal with reading a parameter out to an image */
	      int nVoxels = vmvnin.size();
	      Matrix image;
	      image.ReSize(1,nVoxels);

	      if (bval) {
		if (verbose) cout << "Extracting value for parameter:" << param << endl;
		for (int vox = 1; vox <= nVoxels; vox++)
		  {
		    image(1,vox) = vmvnin[vox-1]->means(param);
		  }
	      }
	      else if (bvar) {
		if (verbose) cout << "Extracting variance for parameter:" << param << endl;
		for (int vox = 1; vox <= nVoxels; vox++)
		  {
		    image(1,vox) = vmvnin[vox-1]->GetCovariance()(param,param);
		  }
	      }
	      else if (cvar) {
		if (verbose) cout << "Extracting co-variance for parameter " << param << "with parameter" << cparam << endl;
		for (int vox =1; vox <= nVoxels; vox++)
		  {
		    image(1,vox) = vmvnin[vox-1]->GetCovariance()(param,cparam);
		  }
	      }

	      if (verbose) cout << "Writing output file" << endl;

              volume4D<float> output(mask.xsize(),mask.ysize(),mask.zsize(),1);
	      copybasicproperties(mask,output);
	      output.setmatrix(image,mask);
	      output.setDisplayMaximumMinimum(output.max(),output.min());
	      output.set_intent(NIFTI_INTENT_NONE,0,0,0);
	      save_volume4D(output,outfile);
	    }
	      
	    if (verbose) cout << "Done." << endl;
	    
	    return 0;
	  }
	catch (const Invalid_option& e)
	  {
	    cout << Exception::what() << endl;
	    Usage();
	  }
	catch (Exception)
	  {
	    cout << Exception::what() << endl;
	  }
	catch (...)
	  {
	    cout << "There was an error!" << endl;
	  }

	return 1;
}

void Usage(const string& errorString)
{
  cout << "\nUsage: mvntool <arguments>\n"
       << "Arguments are mandatory unless they appear in [brackets].\n\n";

  cout << " --help : Prints this information." << endl
       << " --input=<MVNfile> : Name of input MVN file." << endl
       << " --param=<n> : Number of parameter to read/replace/insert." << endl << endl
       << " Extract behaviour (default):" << endl
       << " --output=<NIFIfile> : Name of file for output" << endl
       << "   [--val] : Write parameter value to file." << endl
       << "   [--var] : Write parameter variance to file." << endl << endl
       << " Writing behaviour:" << endl
       << "   [--write] : Overwrite an existing parameter" << endl
       << "   [--new] : Insert a new parameter." << endl
       << "[--output]=<MVNfile> : Name of file for output, overwrites input if not set" << endl
       << " --valim=<NIFITfile> : Image to write for mean of parameter." << endl
       << " --varim=<NIFITfile> : Image to write for variance of parameter." << endl
       << " --val=<mean_value>  : Mean value for parameter to be written." << endl
       << " --var=<variance>    : Variance of parameter to be written." << endl
       << endl;
}


#endif //__FABBER_LIBRARYONLY
