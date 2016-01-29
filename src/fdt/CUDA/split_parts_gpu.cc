/*  split_parts_gpu.cc

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

#include "newimage/newimageall.h"
#include <string>
#include <stdio.h>
#include <stdlib.h>

void save_part(Matrix data, string name, int idpart){

	int nvox = data.Ncols();
	int ndir = data.Nrows();

	string file_name;
	file_name = name+"_"+num2str(idpart);

	ofstream out;
	out.open(file_name.data(), ios::out | ios::binary);
	out.write((char*)&data(1,1),nvox*ndir*sizeof(Real));
	out.close();
}


// parameters:
// 1. data.nii.gz
// 2. mask.nii.gz
// 3. grad_dev.nii.gz
// 4. gflag
// 5. nparts
// 6. output directory


int main(int argc, char *argv[]){

	std::string data_str = argv[1];
	std::string mask_str = argv[2];
	std::string grad_str = argv[3];
	int gflag = atoi(argv[4]);
	int nparts = atoi(argv[5]);
	std::string out_dir = argv[6];

	//printf("%s\n%s\n%s\n%i\n%i\n%s\n",data_str.data(),mask_str.data(),grad_str.data(),gflag,nparts,out_dir.data());
	NEWIMAGE::volume4D<float> data;
	NEWIMAGE::volume<float> mask;
    	read_volume4D(data,data_str);
    	read_volume(mask,mask_str);
	Matrix datam;
    	datam=data.matrix(mask); 

	int nvoxels=datam.Ncols();
	int ndirections=datam.Nrows();
	
	NEWIMAGE::volume4D<float> grad; 
	Matrix gradm;
	int dirs_grad=0;
	if(gflag){
		read_volume4D(grad,grad_str);
      		gradm=grad.matrix(mask);
		dirs_grad = gradm.Nrows();
	}
	

	int size_part=nvoxels/nparts;

	Matrix data_part;
	Matrix grad_part;
	string out_data;
	string out_grad;
	out_data.append(out_dir);
	out_grad.append(out_dir);
	out_data.append("/data");
	out_grad.append("/grad_dev");
	for(int i=0;i<(nparts-1);i++){
		data_part = datam.SubMatrix(1,ndirections,i*size_part+1,(i+1)*size_part);
		save_part(data_part,out_data,i);
		if(gflag){
			grad_part = gradm.SubMatrix(1,dirs_grad,i*size_part+1,(i+1)*size_part);
			save_part(grad_part,out_grad,i);
		}
	}
	// last part
	data_part = datam.SubMatrix(1,ndirections,(nparts-1)*size_part+1,nvoxels);
	save_part(data_part,out_data,(nparts-1));
	if(gflag){
		grad_part = gradm.SubMatrix(1,dirs_grad,(nparts-1)*size_part+1,nvoxels);
		save_part(grad_part,out_grad,(nparts-1));
	}
}
