/*  solver_mult_inverse.cu

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

#include "options.h"

//X = A.i() * B . Used in Levenberg-Marquardt
//MATRIX INVERSE AS NEWMAT LU SOLVER
//implemented in NEWMAT:newmat7.cpp GeneralSolvI.
__device__ void solver(	//INPUT
			float *A, 
			float *P,
			int length,
			//TO USE
			float *C,
			float *el,
			int *indx,	
			//OUTPUT
			float *B)  
{  
	//double C[NPARAMS*NPARAMS];

	for(int i=0;i<length;i++){
		for(int j=0;j<length;j++){
			C[i*length+j]=A[i*length+j];
		}
	}
	
 	bool d=true; 
  	//int indx[NPARAMS];

   	float* akk = C;   
	float big = fabs(*akk); 
	int mu = 0; 
	float* ai = akk; 
	int k;

	for (k = 1; k<length; k++){
      		ai += length; 
		const float trybig = fabs(*ai);
      		if (big < trybig){ 
			big = trybig; 
			mu = k; 
		}
   	}

   	if(length) for (k = 0;;){

		indx[k] = mu;
		if (mu != k){
         		float* a1 = C + length*k; 
			float* a2 = C + length*mu; 
			d = !d;
         		int j = length;
         		while (j--){ 
				const float temp = *a1; 
				*a1++ = *a2; 
				*a2++ = temp; 
			}
      		}

      		float diag = *akk; 
		big = 0; 
		mu = k + 1;
      		if (diag != 0){
         		ai = akk; 
			int i = length - k - 1;
         		while (i--){
            			ai += length; 
				float* al = ai; 
				float mult = *al / diag; 
				*al = mult;
            			int l = length - k - 1; 
				float* aj = akk;
				if (l-- != 0){
				
					float aux=al[1]-(mult* *(++aj));
					*(++al) = aux;
					//*(++al) = __dadd_rn (*al,-mult* *(++aj)); //FAIL in cuda 4.2 compiler
					
               				const float trybig = fabs(*al);
               				if (big < trybig){ 
						big = trybig; 
						mu = length - i - 1; 
					}
               				while (l--){ 
						float aux= al[1]-(mult* *(++aj));
						*(++al) = aux;
						//*(++al) = __dadd_rn (*al,-mult* *(++aj)); //FAIL in cuda 4.2 compiler
					}
           			 }
         		}
      		}
      		if (++k == length) break;      
      		akk += length + 1;
   	}


//////////////////////////////

	//double el[NPARAMS];

	for(int e=0;e<length;e++){
		el[e]=P[e];		
    	}
		
   	int j;
	int ii = length; 
	int ip;    
	float temp;
	int i;
     
	for (i=0; i<length; i++){
 		ip = indx[i]; 
		temp = el[ip]; 
		el[ip] = el[i];
		el[i] = temp;
      		if (temp != 0.0) { ii = i; break; }
   	}
	
  	float* bi; 
	float* ai2;
   	i = ii + 1;

  	if (i < length){
      		bi = el + ii; 
		ai2 = C + ii + i * length;
      		for (;;){
         		int ip = indx[i]; 
			float sum = el[ip]; 
			el[ip] = el[i];
         		float* aij = ai2; 
			float* bj = bi; 
			j = i - ii;
         		while (j--){ 
				sum -=  *aij++* *bj++; 
			}
         		el[i] = sum;
         		if (++i == length) break;
         		ai2 += length;
      		}
   	}

   	ai2 = C + length*length;

   	for (i = length - 1; i >= 0; i--){
      		float* bj = el+i; 
		ai2 -= length; 
		float* ajx = ai2+i;
      		float sum = *bj; 
		float diag = *ajx;
      		j = length - i; 
		while(--j){ 
			sum -= *(++ajx)* *(++bj);  
		}
      		el[i] = sum / diag;
			
   	}
	for(int e=0;e<length;e++){
		B[e]=el[e];
    	}
}

