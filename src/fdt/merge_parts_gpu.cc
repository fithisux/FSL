/*  merge_parts_gpu.cc

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

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

#include "xfibresoptions.h"
#include "newmat.h"
#include "newimage/newimageall.h"

using namespace Xfibres;


void join_Parts(NEWIMAGE::volume<float> mask, string name_in, string name_out, string subjdir, int nvox, int nsamples, int nParts, float max, float min){
		int size_parts = nvox/nParts;
		int last_part = nvox - ((nParts-1)*size_parts);

		//if(mean) nsamples=1;

		Matrix result(nsamples,0);
		Matrix part;

		for(int i=0;i<(nParts-1);i++){
			part.ReSize(nsamples,size_parts);
			std::ostringstream num;
			num << i;
			std::string part_number;
			part_number.assign(num.str());
			std::string aux;
			while(part_number.size()<4){
				aux = "0" + part_number;
				part_number.assign(aux);
			}
			std::string file_name;
			file_name.assign(subjdir);
			file_name += ".bedpostX/diff_parts/data_part_";
			file_name += part_number; 
			file_name += "/"; 
			file_name += name_in; 
			file_name += "J"; 

			ifstream in;
			in.open (file_name.data(), ios::in | ios::binary);
			in.read((char*)&part(1,1), size_parts*nsamples*sizeof(Real));
			in.close();

			result |= part;
		}
		part.ReSize(nsamples,last_part);
		std::ostringstream num;
		num << nParts-1;;
		std::string part_number;
		part_number.assign(num.str());
		std::string aux;
		while(part_number.size()<4){
			aux = "0" + part_number;
			part_number.assign(aux);
		}
		std::string file_name;
		file_name.assign(subjdir);
		file_name += ".bedpostX/diff_parts/data_part_";
		file_name += part_number; 
		file_name += "/"; 
		file_name += name_in;
		file_name += "J";  

		ifstream in;
		in.open (file_name.data(), ios::in | ios::binary);
		in.read((char*)&part(1,1), last_part*nsamples*sizeof(Real));
		in.close();
		result |= part;

		NEWIMAGE::volume4D<float> tmp;
      		tmp.setmatrix(result,mask);
		if(max==-10) max=tmp.max();
		if(min==-10) min=tmp.min(); 
		tmp.setDisplayMaximumMinimum(max,min);
     		save_volume4D(tmp,subjdir+".bedpostX/"+name_out);
}

//////////////////////////////////////////////////////////
//       MERGE THE OUTPUTS FILES OF BEDPOSTX
//////////////////////////////////////////////////////////
//parameters:
// argc - 3 nvox
// argc - 2 nParts
// argc - 1 subjdir

int main(int argc, char *argv[])
{
	// Setup logging:
    	Log& logger = LogSingleton::getInstance();
    	xfibresOptions& opts = xfibresOptions::getInstance();
	opts.parse_command_line(argc-3,argv,logger);

    	NEWIMAGE::volume<float> mask;
    	read_volume(mask,opts.maskfile.value());

	cout << opts.maskfile.value() << endl;
    	
	int nvox = atoi(argv[argc-3]);
	int nParts = atoi(argv[argc-2]);
	string subjdir = argv[argc-1];

	int nsamples = opts.njumps.value()/opts.sampleevery.value();
	
	//////////////////////////////////////////////////////////////
	////////// JOIN Results of the Parts //////////////////////
	//////////////////////////////////////////////////////////////

	if(opts.modelnum.value()==1){
		join_Parts(mask,"mean_dsamples","mean_dsamples",subjdir, nvox, 1, nParts, -10, 0);
	}else if(opts.modelnum.value()>=2){
		join_Parts(mask,"mean_dsamples","mean_dsamples",subjdir, nvox, 1, nParts, -10, 0);
		join_Parts(mask,"mean_d_stdsamples","mean_d_stdsamples",subjdir, nvox, 1, nParts, -10, 0);
		//join_Parts(mask,"dsamples","dsamples",subjdir, nvox, nsamples, nParts, -10, 0);
		//join_Parts(mask,"d_stdsamples","d_stdsamples",subjdir, nvox, nsamples, nParts, -10, 0);
		if(opts.modelnum.value()==3){
			join_Parts(mask,"mean_Rsamples","mean_Rsamples",subjdir, nvox, 1, nParts, 1, 0);
		}	
	}	
	if (opts.f0.value()){
      		join_Parts(mask,"mean_f0samples","mean_f0samples",subjdir, nvox, 1, nParts, 1, 0);
		//join_Parts(mask,"f0samples","f0samples",subjdir, nvox, nsamples, nParts, 1, 0);
    	}
    	if (opts.rician.value()){
		join_Parts(mask,"mean_tausamples","mean_tausamples",subjdir, nvox, 1, nParts, -10, 0);
    	}
	join_Parts(mask,"mean_S0samples","mean_S0samples",subjdir, nvox, 1, nParts, -10, 0);

	for(int f=0;f<opts.nfibres.value();f++){
		join_Parts(mask,"th"+num2str(f+1)+"samples","merged_th"+num2str(f+1)+"samples",subjdir, nvox, nsamples, nParts, -10, -10);
		join_Parts(mask,"ph"+num2str(f+1)+"samples","merged_ph"+num2str(f+1)+"samples",subjdir, nvox, nsamples, nParts, -10, -10);
		join_Parts(mask,"f"+num2str(f+1)+"samples","merged_f"+num2str(f+1)+"samples",subjdir, nvox, nsamples, nParts, 1, 0);

		//join_Parts(mask,"mean_f"+num2str(f+1)+"samples",subjdir, nvox, 1, nParts, 1, 0);
		//join_Parts(mask,"dyads"+num2str(f+1),subjdir, nvox, nsamples, nParts, 1, -1);
	}	
		
  	return 0;
}

