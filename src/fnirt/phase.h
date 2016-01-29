
/*  phase.h

    Matt Mellor, Medical Vision Laboratory

    Copyright (C) 2004 University of Oxford  */

 
// Monogenic - a class for computing the monogenic signal
// 

#if !defined(__phase_h)
#define __phase_h

#include "newimage/newimageall.h"

 using namespace NEWIMAGE;

 namespace MONOGENIC {

class Monogenic{
 

 public:
    // Publicly available calls
    Monogenic(const volume<float>& img);
    Monogenic(const volume<float>& img, string& filter_type);
    Monogenic(const volume<float>& img, string& filter_type, float parameter_1, float parameter_2);
    ~Monogenic();
 
    complexref operator()(int x,int y, int z);
     

    int xsize();
    int ysize();
    int zsize();

    volume<float> energy() const;
    volume<float> phase() const;
    volume<float>& re();
    volume<float>& im();
  

 private:
    // private calls
    int gau_diff(float sigma1, float sigma2);
    int range(float alpha, float beta);
    int transform_est(void);
    int apply_filter();

    // Private data
    const volume<float> &image;
    complexvolume transform;
    float f_parm1, f_parm2;
    string f_type;
 
};

   //////////////////////////////////////////////////////////////////////////


 }

#endif //__phase_h

