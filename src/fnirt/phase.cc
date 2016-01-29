/*  phase.cc

    Matt Mellor, Medical Vision Lab

    Copyright (C) 2004 University of Oxford  */

#include "phase.h"
#include "newimage/imfft.h"
#include <iostream>

/////////////////////////////////////////////////
//-------------- FILTER DEFAULTS --------------//
/////////////////////////////////////////////////

#define DEFAULT_FILTER "range"

#define DEFAULT_PRM_1 4.5
#define DEFAULT_PRM_2 0.05

#define DEFAULT_RANGE_PRM_1 4.5
#define DEFAULT_RANGE_PRM_2 0.05

#define DEFAULT_GAUDIFF_PRM_1 4
#define DEFAULT_GAUDIFF_PRM_2 8



/////////////////////////////////////////////////
//-------------- CONSTRUCTORS -----------------//
/////////////////////////////////////////////////

namespace MONOGENIC {
// initialising image only
Monogenic::Monogenic(const volume<float>& img) :
     image(img)
{
  // initialise the trasnform and the filter
  transform = complexvolume(image,image*0);

  // set up the filter
  f_type = DEFAULT_FILTER;
  f_parm1 = DEFAULT_PRM_1;
  f_parm2 = DEFAULT_PRM_2;
	
  // perform the transform
  transform_est();
}

// initialising image and filter type
Monogenic::Monogenic(const volume<float>& img, string& filter_type) :
  image(img), f_type(filter_type)
{
  // initialise the trasnform and the filter
  transform = complexvolume(image,image*0);

  // set up the filter
  if( f_type == "range" ){
	f_parm1 = DEFAULT_RANGE_PRM_1;
	f_parm2 = DEFAULT_RANGE_PRM_2;
  }else if( f_type == "gaudiff" ){
        f_parm1 = DEFAULT_GAUDIFF_PRM_1;
        f_parm2 = DEFAULT_GAUDIFF_PRM_2;
  }else{
        cerr << "unknown filter type in Monogenic(const volume<float>& img, string& filter_type)\n";
        f_type = DEFAULT_FILTER;
        f_parm1 = DEFAULT_PRM_1;
        f_parm2 = DEFAULT_PRM_2;
  }
 

  // perform the transform
  transform_est();
}

// initialising image, filter type and parameters
Monogenic::Monogenic(const volume<float>& img, string& filter_type, float parameter_1, float parameter_2) :
  image(img), f_parm1(parameter_1), f_parm2(parameter_2), f_type(filter_type) 
{
  // initialise the trasnform and the filter
  transform = complexvolume(image,image*0);

	
  // perform the transform
  transform_est();
}


Monogenic::~Monogenic(){
  transform.destroy();
}

/////////////////////////////////////////////////
//-------------- FILTERING --------------------//
/////////////////////////////////////////////////

// Note: all filters assume that 'transform' has
// been Fourier transformed already. This saves
// on FTing when applying multiple filters.


// apply the appropriate filter
int Monogenic::apply_filter(){

  if ( f_type == "range" ){
    range(f_parm1,f_parm2);
  }else if ( f_type == "gaudiff" ){
    gau_diff(f_parm1,f_parm2);
  }else{
    cerr << "Unknown filter type in Monogenic::apply_filter()\n";
    return -1;
  }

  return 1;
}


// 'Range' filtering

int Monogenic::range(float alpha, float beta){

  double xo,yo,zo;
  int width_x,width_y,width_z;
  double r;
  double filtsum=0,tempsum=0;
  double filt_val, temp_val;
  complexvolume filter(this->xsize(),this->ysize(),this->zsize());
  complexvolume tempfilt(this->xsize(),this->ysize(),this->zsize());
 
  width_x = transform.xsize();
  width_y = transform.ysize();
  width_z = transform.zsize();

  float num_den = 0; // numbre by which we have to divide sum_filter to find the value in filter(0, 0, 0)

  for(int x=0; x<width_x/2; x++){
  for(int y=0; y<width_y/2; y++){
  for(int z=0; z<width_z/2; z++){

    
    xo = x; // we don't shift by 0.5 -- Ben
    yo = y; // we don't shift by 0.5 -- Ben
    zo = z; // we don't shift by 0.5 -- Ben
/*
    xo = x + 0.5;
    yo = y + 0.5;
    zo = z + 0.5;
*/
    if ((x != 0) || (y != 0) || (z != 0)) {
      r = sqrt(xo*xo + yo*yo + zo*zo);

      filt_val = 1/pow(r,alpha+beta);
      temp_val = 1/pow(r,alpha-beta);
    
      filter(x,y,z) = filt_val;
      filter(width_x-x, y, z) = filt_val;
      filter(x, width_y-y, z) = filt_val;
      filter(x, y, width_z-z) = filt_val;
      filter(width_x-x, y, width_z-z) = filt_val;
      filter(width_x-x, width_y-y, z) = filt_val;
      filter(x, width_y-y, width_z-z) = filt_val;
      filter(width_x-x, width_y-y, width_z-z) = filt_val;

      tempfilt(x,y,z) = temp_val;
      tempfilt(width_x-x, y, z) = temp_val;
      tempfilt(x, width_y-y, z) = temp_val;
      tempfilt(x, y, width_z-z) = temp_val;
      tempfilt(width_x-x, y, width_z-z) = temp_val;
      tempfilt(width_x-x, width_y-y, z) = temp_val;
      tempfilt(x, width_y-y, width_z-z) = temp_val;
      tempfilt(width_x-x, width_y-y, width_z-z) = temp_val;

      filtsum = 8 * filt_val; // filt_val is 8 times in filter volume.
      tempsum += 8 * tempfilt.re(x,y,z);
      num_den += 8;
    }
    
  }}}

  // we treat filter(0, 0, 0)
    filter(0, 0, 0) = filtsum / num_den;
    filtsum += filter.re(0, 0, 0);
    tempsum += tempfilt.re(0, 0, 0);
  
  // this should be done inside the previous loop really
/*
  for(int x=0; x<width_x; x++){
  for(int y=0; y<width_y; y++){
  for(int z=0; z<width_z; z++){

    filtsum += filter.re(x,y,z);
    tempsum += tempfilt.re(x,y,z);
  }}}
*/
 
  // normalise the filters to guarantee zero DC
  filter = filter/filtsum;
  tempfilt = tempfilt/tempsum;

  // the filter is actually the difference of the two
  filter -= tempfilt;
  
  // take the FT and apply to the transform
  fft3(filter);
  transform = transform*filter;

  return 1;
}


// difference of Gausian filtering
// This would be alot more efficient if the 
// filters were implemented directly in the 
// Fourier domain.
int Monogenic::gau_diff(float sigma1, float sigma2){

  double xo,yo,zo;
  int width_x,width_y,width_z;
  double r_sq;
  double filtsum=0,tempsum=0;
  double filt_val, temp_val;
  complexvolume filter(this->xsize(),this->ysize(),this->zsize());
  complexvolume tempfilt(this->xsize(),this->ysize(),this->zsize());
 
  width_x = transform.xsize();
  width_y = transform.ysize();
  width_z = transform.zsize();

  for(int x=0; x<width_x/2; x++){
  for(int y=0; y<width_y/2; y++){
  for(int z=0; z<width_z/2; z++){

    xo = x;
    yo = y;
    zo = z;

    r_sq = xo*xo + yo*yo + zo*zo;

    // I added a cast in the pow function -- Ben
    filt_val = exp(-r_sq/(pow((double) sigma1, (double) 2)));
    temp_val = exp(-r_sq/(pow((double) sigma2, (double) 2)));
    
    filter(x,y,z) = filt_val;
    filter(width_x-x, y, z) = filt_val;
    filter(x, width_y-y, z) = filt_val;
    filter(x, y, width_z-z) = filt_val;
    filter(width_x-x, y, width_z-z) = filt_val;
    filter(width_x-x, width_y-y, z) = filt_val;
    filter(x, width_y-y, width_z-z) = filt_val;
    filter(width_x-x, width_y-y, width_z-z) = filt_val;

    tempfilt(x,y,z) = temp_val;
    tempfilt(width_x-x, y, z) = temp_val;
    tempfilt(x, width_y-y, z) = temp_val;
    tempfilt(x, y, width_z-z) = temp_val;
    tempfilt(width_x-x, y, width_z-z) = temp_val;
    tempfilt(width_x-x, width_y-y, z) = temp_val;
    tempfilt(x, width_y-y, width_z-z) = temp_val;
    tempfilt(width_x-x, width_y-y, width_z-z) = temp_val;

  }}}


  // this should be done inside the previous loop really

  for(int x=0; x<width_x; x++){
  for(int y=0; y<width_y; y++){
  for(int z=0; z<width_z; z++){

    filtsum += filter.re(x,y,z);
    tempsum += tempfilt.re(x,y,z);

  }}}
 
  // normalise the filters to guarantee zero DC
  filter = filter/filtsum;
  tempfilt = tempfilt/tempsum;

  // the filter is actually the difference of the two
  filter -= tempfilt;
  
  // take the FT and apply to the transform
  fft3(filter);
  transform = transform*filter;

  return 1;
}

//////////////////////////////////////////////////
//-------- COMPUTE MONOGENIC SIGNAL ------------//
//////////////////////////////////////////////////

// compute the monogenic signal, assuming that the complex 
// volume 'transform' presently contains the image.
// The image is first filtered with a bandpass.

int Monogenic::transform_est(void){

  int xo,yo,zo;
  int rx,ry,rz;
  complexpoint h1l,h2l,h3l,tl; 
  float r,odd;
  complexvolume h1(transform.xsize(),transform.ysize(),transform.zsize());
  complexvolume h2(transform.xsize(),transform.ysize(),transform.zsize());
  complexvolume h3(transform.xsize(),transform.ysize(),transform.zsize());

  // compute FT
  fft3(transform);

  // apply filter
  apply_filter();
 
  // calculate the three oriented monogenic components
  
  rx = transform.xsize()/2;
  ry = transform.ysize()/2;
  rz = transform.zsize()/2;

  for(int x=0; x<transform.xsize(); x++){
  for(int y=0; y<transform.ysize(); y++){
  for(int z=0; z<transform.zsize(); z++){

    xo = x-rx;
    yo = y-ry;
    zo = z-rz;

    // Be careful I (Ben) added transtyping in the next line ("sqrt"
    // function).
    r = (float) sqrt((double) (xo*xo + yo*yo + zo*zo));
    
    tl = transform(x,y,z); 

    h1(x,y,z) = tl*float(xo)/(r+0.001);
    h2(x,y,z) = tl*float(yo)/(r+0.001);
    h3(x,y,z) = tl*float(zo)/(r+0.001);

  }}}
   
  ifft3(h1);
  ifft3(h2);
  ifft3(h3);
  ifft3(transform);
 
  for(int x=0; x<transform.xsize(); x++){
  for(int y=0; y<transform.ysize(); y++){
  for(int z=0; z<transform.zsize(); z++){

    h1l = h1(x,y,z); 
    h2l = h2(x,y,z); 
    h3l = h3(x,y,z);
 
    odd = sqrt(h1l.im()*h1l.im() + h2l.im()*h2l.im() + h3l.im()*h3l.im());
 
    transform(x,y,z) = complexpoint(transform.re(x,y,z),odd);

  }}} 

  return 1;

}



///////////////////////////////////////////
//-----------SIZE FUNCTIONS--------------//
///////////////////////////////////////////

int Monogenic::xsize(){
  return transform.xsize();
}

int Monogenic::ysize(){
  return transform.ysize();
}

int Monogenic::zsize(){
  return transform.zsize();
}



///////////////////////////////////////////
//----------ACCESS FUNCTIONS-------------//
///////////////////////////////////////////

volume<float> Monogenic::energy() const{
  return transform.abs();
}

volume<float> Monogenic::phase() const{
  return transform.phase();
}

volume<float>& Monogenic::re(){
  return transform.re();
}

volume<float>& Monogenic::im(){
  return transform.im();
}

complexref Monogenic::operator()(int x,int y, int z){
 return transform(x,y,z);
}
}
