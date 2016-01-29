/*************************************
Code by: Tejas Pendse; 31/8/12

Usage: "resample -i <input obj/phantom> -o <output> -m <MR param file> [optional/diagnostic params]"
Optional params:
-e,--te : Scanner Echo time
-r,--tr : Scanner Repetition time

## Theory calc currently only works for slice select gradient along Z! ##
*************************************/

#include "utils/options.h"
#include "newimage/newimageall.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include <time.h>

using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;

#define pi 3.1416

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="resample --Tejas Pendse 2012";
string examples="resample -i <input image> [-options] -o <output>";

/**Required**/
Option<string> opt_input(string("-i,--input"), "", 
		     string("Input Image (4D phantom for theoretical calculations)"), 
		     true, requires_argument);

Option<string> opt_output(string("-o,--output"), "", 
		     string("Output Image"), 
		     true, requires_argument);

/**Scanner options**/

Option<float> opt_te(string("-e,--te"), 0,
		     string("        Echo Time (TE) in seconds \t\t[For example: T1-weighted images for 3T TE=0.01 s]"), 
		     true, requires_argument);

Option<float> opt_tr(string("-r,--tr"), 0, 
		     string("        Repetition Time (TR) in seconds \t[For example: T1-weighted images for 3T TR=0.7 s]"), 
		     true, requires_argument);

Option<string> opt_mrpar(string("-m,--mrpar"), "", 
		     string("MRpar File"), 
		     true, requires_argument);

/**Optional**/
	/**Image options**/

Option<int> opt_nx(string("-x,--nx"), 0, 
		     string("        Number of Voxels along X (default: phantom)"), 
		     false, requires_argument);

Option<int> opt_ny(string("-y,--ny"), 0, 
		     string("        Number of Voxels along Y (default: phantom)"), 
		     false, requires_argument);

Option<int> opt_nz(string("-z,--nz"), 0, 
		     string("        Number of Voxels along Z (default: phantom)"), 
		     false, requires_argument);

Option<float> opt_dx(string("-a,--dx"), 0.0, 
		     string("        Size of voxels along X(default: phantom)"), 
		     false, requires_argument);

Option<float> opt_dy(string("-b,--dy"), 0.0, 
		     string("        Size of voxels along Y(default: phantom)"), 
		     false, requires_argument);

Option<float> opt_dz(string("-c,--dz"), 0.0, 
		     string("        Size of voxels along Z i.e., number of slices (default: phantom)"), 
		     false, requires_argument);

Option<float> opt_zstart(string("--zstart"), 0.0,
		     string("Starting position of the volume in mm (default = 0mm)"),
		     false,requires_argument);

Option<float> opt_sigma(string("--sigma"), 0.0, 
		     string("\tAdd noise with given sigma (default: 0 i.e., no noise)"), 
		     false, requires_argument);

Option<bool> opt_save_og(string("-s,--save"), false, 
		     string("Save original non-resample output image"), 
		     false, no_argument);


	/**Diagnostic**/
Option<bool> help(string("-?,--help"), false, 
		     string("Display help message"), 
		     false, no_argument);

Option<bool> verbose(string("-v,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);

int nonoptarg;
/******
'do_work' : (Step 5-6 can/should be done in k-space !? Not centered, so nope ?!)
1. Read input files/images
2. If MRpar file is given, compute theory value
	2.1. Compute theory volume
	2.2. Save volume if defined
3. Calculate output image size in mm = dimensions*voxel-size
	3.1. If output size is smaller than phantom size, isotropic resample by appropriate scale, i.e. w.r.t. largest dimesion diff
4. Pad zeroes according to difference between phantom and output mm sizes (right padded more (?why-check?))
5. Resample each slice - same scale since Nx=Ny
6. Resample along Z (!Trying to emulate possum!)
	6.1. Resample by scale of voxel size
	6.2. Select number of slices
7. Save output

******/
int do_work(int argc, char* argv[])
{
	//Read inputs
	bool compute_theory = true;
		
	//Read inputs: Scanner options
	double te = opt_te.value();
	double tr = opt_tr.value();

	//Check if noise is to be added
	bool addNoise = false;
	float sigma = opt_sigma.value();
	if(sigma != 0)
	{
		if(sigma > 0)
			addNoise = true;
		else
			cerr << "Error: Sigma can't be negative!" << endl;
	}

	//MRpar matrix file
	Matrix mrpar;

	//Read inputs: Sampling dimensions
	RowVector imdims(3);
	RowVector voxdims(3);
	imdims(1) = opt_nx.value();
	imdims(2) = opt_ny.value();
	imdims(3) = opt_nz.value();
	voxdims(1) = opt_dx.value();
	voxdims(2) = opt_dy.value();
	voxdims(3) = opt_dz.value();
	int zstart = (int)floor(opt_zstart.value() / voxdims(3));

	//Check inputs
	for( int i=1; i<=3; i++)
	{
		if(imdims(i) < 0)
		{
			cerr << "Error: Image dimensions must be positive!";
			cerr << "--IMG-DIM " << i << "=" << imdims(i) << endl;
			exit(0);
		}
		if(voxdims(i) < 0)
		{
			cerr << "Error: Voxel dimensions must be positive!\t";
			cerr << "--VOXEL-DIM " << i << "=" << voxdims(i) << endl;
			exit(0);
		}
	}

	//Read inputs: input/phantom image spec
	RowVector ph_imdims(3);
	RowVector ph_voxdims(3);

	if(verbose.value()) cout << "Reading input volume: " << opt_input.value() << endl;

	//Read inputs: input/phantom
	volume4D<double> phantom;
	volume<double> input;

	if(opt_mrpar.value() == "")
	{
		compute_theory = false;
		cout << "No MRPar file! Only resampling input 3D image..." << endl;

		if(verbose.value()) 
			if(te != 0 || tr != 0) 
				cerr << "Warning! No MRpar file given, TE/TR will not be used!" << endl;

		read_volume(input,opt_input.value());
		//Image size
		ph_imdims(1) = input.xsize();
		ph_imdims(2) = input.ysize();
		ph_imdims(3) = input.zsize();
		//Voxel size
		ph_voxdims(1) = input.xdim();
		ph_voxdims(2) = input.ydim();
		ph_voxdims(3) = input.zdim();
	}
	else
	{
		if(verbose.value())	cout << "Input object: 4D phantom" << endl;
		read_volume4D(phantom,opt_input.value());

		//Read the MRpar file
		mrpar = read_ascii_matrix(opt_mrpar.value());

		bool thereIsAnError = false;
		if(mrpar.Ncols() < 3)
		{
			cerr << "Error: MRpar file doesn't have enough number of columns!" << endl;
			thereIsAnError++;
		}

		if(mrpar.Nrows() < phantom.tsize())
		{
			cerr  << "Error: MRpar doesn't have enough number of rows!" << endl;
			thereIsAnError++;
		}

		//If theory calculations are needed, check if TE/TR are non-zero
		if(te == 0 || tr == 0)
		{
			cerr << "Error: Theory calculations need non-zero TE/TR!" << endl;
			thereIsAnError++;
		}

		if(thereIsAnError)	exit(EXIT_FAILURE);

		//Image size
		ph_imdims(1) = phantom.xsize();
		ph_imdims(2) = phantom.ysize();
		ph_imdims(3) = phantom.zsize();
		//Voxel size
		ph_voxdims(1) = phantom.xdim();
		ph_voxdims(2) = phantom.ydim();
		ph_voxdims(3) = phantom.zdim();
	}

	//Check if volume fits
	if(zstart > ph_imdims(3)*ph_voxdims(3))
	{
		cerr << "Error: Dimensions Z-dimensions don't fit in the object" << endl;
		exit(0);
	}

	//Check if any input dimension is 0 and replace by phantom dimension
	for( int i=1; i<=3; i++)
	{
		//If this is true, then voxels are given
		if(imdims(i) == 0)
			imdims(i) = ph_imdims(i);
		if(voxdims(i) == 0)
		{
			voxdims(i) = ph_imdims(i)*ph_voxdims(i)/imdims(i);
			if(verbose.value())
			{
				//unnecessary really
				cout << "Output voxel ";
				if(i==1) cout << "X-size calculated as:\t" << voxdims(i) << endl;
				else if(i==2) cout << "Y-size calculated as:\t" << voxdims(i) << endl;
				else if(i==3) cout << "Z-size calculated as:\t" << voxdims(i) << endl;
			}
		}
	}

	volume<double> output;
	//All initialisations done

/****************************
Calculate theory volume
****************************/

	if(compute_theory)
	{
		output.reinitialize((int)ph_imdims(1),(int)ph_imdims(2),(int)ph_imdims(3));	output=0;
		if(verbose.value()) cout << "Computing theoretical volume...\nCalculating signal values..." << endl;

		//Signal stores NMR signal from each tissue
		RowVector signal(phantom.tsize());
		for (int i=1; i<=phantom.tsize(); i++)
		{
			double t1 = mrpar(i,1);
			double t2 = mrpar(i,2);
			//Theory equation
			//signal(i) = (1-exp(-tr/mrpar(i,1))) * exp(-te/mrpar(i,2));	//Old equations: specific for T1
			//Changed 14.12.12 - new equation generalised
			signal(i) = mrpar(i,3)* (1 - 2*exp( (-tr + te/2)/t1 ) + exp( -tr/t1 )) * exp( -te/t2 );
            if(verbose.value()) cout << "Contrast values are: " << signal(i)<<endl;
		}

		if(verbose.value()) cout << "Calculating output values..." << endl;

    //Loop through all phantom voxels	
    
        //if(verbose.value()) cout << "Tissue: " << t << endl;
        for( int z=phantom.minz(); z <= phantom.maxz(); z++)
            for( int y=phantom.miny(); y <= phantom.maxy(); y++)
                for( int x=phantom.minx(); x <= phantom.maxx(); x++)
                    for( int t=0; t < phantom.tsize(); t++)
                    {
                    output(x,y,z) += phantom(x,y,z,t) * signal(t+1);
                    }

		//if any noise is to be added, add it to this volume
		//Code taken from possum:systemnoise
		if(addNoise)
		{
			if(verbose.value())	cout << "Adding noise..." << endl;
			//Was trying to mimic what 'possumX_postproc.sh' & 'systemnoise' do. But for the moment, just sigma as the option
			//float p98 = output.percentile(0.98);
			//float p02 = output.percentile(0.02);
			//float thresh = 0.1*p98 + 0.9*p02;
			//output.threshold(thresh);
			//float medint = output.percentile(0.5);
			//float sigma = medint/snr;
			//sigma = sigma/output.xsize();

			Matrix noise(output.xsize(),output.ysize());
			time_t loc;
			time(&loc);
			int seed_val=loc%105634;
			for(int z=output.minz(); z<=output.maxz(); z++)
			{
				srand(seed_val+z);
				noise = normrnd(output.xsize(),output.ysize(),0,sigma);
				for(int y=output.miny(); y<=output.maxy(); y++)
					for(int x=output.minx(); x<=output.maxx(); x++)
						output(x,y,z) += noise(x-output.minx()+1,y-output.miny()+1);
			}
		}

		//Output here is the original theory volume, i.e. not resampled
		if(opt_save_og.value())
		{
			if(verbose.value()) cout << "Saving original volume as '" << opt_output.value() <<"_full.nii.gz'..." << endl;
			save_volume(output,opt_output.value()+"_full");
		}
	}
	else
	{
		if(verbose.value()) cout << "No theory calculations..." << endl;
		//Expected output from theory initialized as the input
		output = input;
	}

/**
1. Check if any dimensions smaller than phimdims*phvoxdims
2. Pad to next largest factor of 2 and subsample by 2

**/
	int excess;
	int lgrdim=0;
	bool isSmaller = false;
	RowVector size_diff(3); size_diff=0;
	RowVector lfpad(3);	lfpad=0;
	RowVector rtpad(3);	rtpad=0;
	RowVector limits(3);	limits=0;

	//Since output is to be resampled, output mm > phantom mm
	//Check if any intended dimensions are smaller than phantom, else resize phantom.
	for(int i=1;i<=2;i++)
	{
		size_diff(i) = ph_imdims(i)*ph_voxdims(i) - imdims(i)*voxdims(i);
		if(size_diff(i) > 0)	isSmaller++;
	}

	if(isSmaller)
	{
		if(verbose.value())	cout << "..." << endl;
		//Find the largest dimension
		if(size_diff(1) > size_diff(2))
			lgrdim = 1;
		else
			lgrdim = 2;

		//Resample image by scale that elliminates largest difference
		float scale = (ph_imdims(lgrdim)*ph_voxdims(lgrdim))/(imdims(lgrdim)*voxdims(lgrdim));
		float stepx = scale / output.xdim();
		float stepy = scale / output.ydim();
		int sx,sy;
		sx = (int) Max(1.0, ( ((float) (output.maxx() - output.minx() + 1.0)) / stepx));
		sy = (int) Max(1.0, ( ((float) (output.maxy() - output.miny() + 1.0)) / stepy));
		volume<double> output_res(sx,sy,output.zsize());
		float fx,fy;
		int x1,y1;
		for(int z=output.minz(); z<=output.maxz(); z++)
			for( fy=0.0, y1=0; y1<sy; y1++, fy += stepy )
				for( fx=0.0, x1=0; x1<sx; x1++, fx += stepx )
					output_res(x1,y1,z) = output.interpolate(fx,fy,z);
		
		output.reinitialize(output_res);

		//New std output:voxels dims are 1
		output.setxdim(1);	output.setydim(1);	output.setzdim(1);
		ph_imdims(1) = output.xsize();	ph_imdims(2) = output.ysize();	ph_imdims(3) = output.zsize();
	}

	//For Z-dir : rtpad pads to bottom; lfpad pads above
	for(int i=1;i<=3;i++)
	{
		//Calculate the extra voxels to be padded
		excess = (int)abs(imdims(i)*voxdims(i) - ph_imdims(i)*ph_voxdims(i));
		limits(i) = ph_imdims(i) + excess;
		if(excess > 0)
		{
			//For odd excess : pad 1 voxel more on the 'right'
			lfpad(i) = floor(excess/2);
			if(excess % 2 != 0)
				rtpad(i) = floor(excess/2)+1;
			else
				rtpad(i) = excess/2;
		}
	}

	//output2 is padded version of prev output
	volume<double> output2((int)(lfpad(1)+output.xsize()+rtpad(1)),(int)(lfpad(2)+output.ysize()+rtpad(2)),output.zsize());	output2=0;

	//Write smaller sized output to padded output2. Offset from origin by lfpad values
	for( int z=output.minz(); z <= output.maxz(); z++)
		for( int y=output.miny(); y <= output.maxy(); y++)
			for( int x=output.minx(); x <= output.maxx(); x++)
				output2((int)(x+lfpad(1)),(int)(y+lfpad(2)),z) = output(x,y,z);

	if(verbose.value())	cout << "Resampling image..." << endl;

	float scalex = voxdims(1)/ph_voxdims(1);
	float scaley = voxdims(2)/ph_voxdims(2);
	float stepx = scalex / output2.xdim();
	float stepy = scaley / output2.ydim();
	int sx = (int) Max(1.0, ( ((float) (output2.maxx() - output2.minx() + 1.0)) / stepx));
	int sy = (int) Max(1.0, ( ((float) (output2.maxy() - output2.miny() + 1.0)) / stepy));

	volume<double> out_res(sx,sy,output2.zsize());

	float fx,fy;
	int x1,y1;

	for( int z=output2.minz(); z<=output2.maxz(); z++ )
		for( fy=0.0, y1=0; y1<sy; y1++, fy += stepy )
			for( fx=0.0, x1=0; x1<sx; x1++, fx += stepx )
				out_res(x1,y1,z) = output2.interpolate(fx,fy,z);

/****************************
Z-Direction RESAMPLING
****************************/

	float scalez = voxdims(3)/ph_voxdims(3);
	float stepz = scalez / out_res.zdim();
	int sz = (int) Max(1.0, ( ((float) (out_res.maxz() - out_res.minz() + 1.0)) / stepz));

	volume<double> out_res_z(out_res.xsize(),out_res.ysize(),sz);
	float fz;
	int z1;

	for( fz=0.0, z1=0; z1<sz; z1++, fz += stepz )
		for( int y=out_res.miny(); y<=out_res.maxy(); y++ )
			for( int x=out_res.minx(); x<=out_res.maxx(); x++ )
				out_res_z(x,y,z1) = out_res.interpolate(x,y,fz);
	
	out_res.reinitialize(out_res_z);

	//Select number of slices
	out_res_z.reinitialize(out_res.xsize(),out_res.ysize(),(int)imdims(3));
	out_res.setROIlimits(out_res.minx(),out_res.miny(),zstart,out_res.maxx(),out_res.maxy(),zstart+(int)imdims(3)-1);
	out_res.activateROI();	out_res_z = out_res.ROI();	out_res.deactivateROI();

	volume<double> output_abs;
	output_abs = out_res_z;

/************************
Finalising

### output_abs : final image to be saved
*************/
	if(verbose.value())	cout << "Setting image properties..." << endl;
	output_abs.setxdim(voxdims(1));	output_abs.setydim(voxdims(2));	output_abs.setzdim(voxdims(3));
	
	if(verbose.value())	cout << "Saving output file(s) '" << opt_output.value() <<".nii.gz'..." << endl;
	save_volume(output_abs,opt_output.value());	//add a "_abs" or "_phase" later, maybe?

//option to calc phase? (from fft maybe)

	return 0;
}

int main(int argc,char *argv[])
{
	Tracer tr("main");
	OptionParser options(title, examples);
	try
	{
		// must include all wanted options here (the order determines how
		//  the help message is printed)
 
		options.add(opt_input);
		options.add(opt_output);

		//Scanner options
		options.add(opt_te);
		options.add(opt_tr);
		options.add(opt_mrpar);

		//Resample options
		options.add(opt_nx);
		options.add(opt_ny);
		options.add(opt_nz);
		options.add(opt_dx);
		options.add(opt_dy);
		options.add(opt_dz);
		options.add(opt_zstart);

		options.add(opt_sigma);
		options.add(opt_save_og);
		options.add(verbose);
		options.add(help);

		nonoptarg = options.parse_command_line(argc, argv);
		
		// line below stops the program if the help was requested or 
		//  a compulsory option was not set
		if( (help.value()) || (!options.check_compulsory_arguments(true)) )
		{
			options.usage();
			exit(EXIT_FAILURE);
		}
	}
	catch(X_OptionError& e)
	{
		options.usage();
		cerr << endl << e.what() << endl;
		exit(EXIT_FAILURE);
	}
	catch(std::exception &e)
	{
		cerr << e.what() << endl;
	}

	// Call the local functions
	return do_work(argc,argv);
}
