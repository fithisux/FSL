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

#if !defined (MISCPLOT_H)
#define MISCPLOT_H

#include <string>
#include "newmatap.h"
#include "miscmaths/histogram.h" 
#include "gd.h"

#ifndef FALSE
#define FALSE false
#endif
#ifndef TRUE
#define TRUE true
#endif


namespace MISCPLOT{

	extern unsigned long sc_init[64];
  class miscplot
  {
    public:
      
   		//constructor
      miscplot(){
				req_xsize=0; req_ysize=0;explabel=string("");bp_colctr=0;
				null_shift=0.0;minmaxscale=0.0;spacing=0;ymin=0.0;ymax=0.0;
				bp_whiskerlength = 1.5; bp_notched = TRUE; 
				histogram_bins = 0;scat_ctr=0;gridswapdefault=false;
				Ylabel_fmt = string("");
				for(int ctr=0; ctr<64; ++ctr)
					sc[ctr]=MISCPLOT::sc_init[ctr];
      };

      //destructor
      ~miscplot() {
				GDCglobals_reset();
      }    
			
			static void Timeseries(const NEWMAT::Matrix& mat, string filename, 
				string title, float tr = 0.0, int ysize = 150, int width = 4, 
				int prec = 2, bool sci = false)
				{	miscplot plot;
					plot.timeseries(mat, filename, title, tr, ysize, width, prec, sci);
				}

      void timeseries(const NEWMAT::Matrix& mat, string filename, 
				string title, float tr = 0.0, int ysize = 150, int width = 4, 
				int prec = 2, bool sci = false);

      void histogram(const NEWMAT::Matrix& mat, string filename, string title);

      void boxplot(string filename, string title);
      void boxplot(const NEWMAT::Matrix& mat, string filename, string title);
      void boxplot(const NEWMAT::ColumnVector& vec, string filename, 
				string title);

      void gmmfit(const NEWMAT::Matrix& mat, Matrix& mu, Matrix& sig, 
		  Matrix& pi, string filename, string title, bool mtype = false, 
		  float offset = 0.0, float detailfactor = 0.0); 
    
      inline void ggmfit(const NEWMAT::Matrix& mat, Matrix& mu, 
				Matrix& sig, Matrix& pi, string filename, string title, 
				float offset = 0.0, float detailfactor = 0.0){
					this->gmmfit(mat, mu, sig, pi, filename, title, true, 
					offset, detailfactor);
				}

    // plot a mixture of K gaussians
    void gmmfit(const NEWMAT::Matrix& mat,const NEWMAT::ColumnVector& mu,const NEWMAT::ColumnVector& var,const NEWMAT::ColumnVector& pi,
		string filename,string title,
		bool mtype = false,float offset = 0.0, float detailfactor = 0.0);
    
      inline void add_label(string txt){ 
				labels.push_back(txt);}
      inline void remove_labels(int i) {
				for (int j=1;j<=i;j++) labels.pop_back();}
			inline void clear_labels(){
				labels.clear();}
      inline void add_xlabel(string txt){ 
				xlabels.push_back(txt);}
      inline void remove_xlabel(){ 
				xlabels.pop_back();}
			inline void clear_xlabel(){
					xlabels.clear();}
      inline void add_ylabel(string txt){ 
				ylabels.push_back(txt);}
			inline void clear_ylabel(){
				ylabels.clear();}
      void setscatter(Matrix &data, int width=15);
      void deletescatter();
			void GDCglobals_reset();
      void add_bpdata(const NEWMAT::Matrix& mat);
      void add_bpdata(const NEWMAT::ColumnVector& vec);
      inline void set_xysize(int xsize, int ysize){ 
				req_xsize=xsize; req_ysize=ysize;}
      inline void set_nullshift(double val){
				null_shift=val;};
      inline void set_minmaxscale(float val){
				minmaxscale = val;};
      inline void set_yrange(float val1, float val2){
				ymin = val1; ymax=val2;};
      inline void set_spacing(int val){
				spacing = val;};
      inline void set_bpnotches(bool val){
				bp_notched  = val;};
      inline void set_bpwhiskerlength(float val){
				bp_whiskerlength = val;};
      inline void set_histogram_bins(int val){
				histogram_bins = val;};
			inline void grid_swapdefault(){
				gridswapdefault = true;};
			inline void set_Ylabel_fmt(string what){
				Ylabel_fmt = what;}
			inline void col_replace(int which, unsigned long what){
				if(which>=0 && which<64)
					sc[which] = what;
			}

    private:
			unsigned long sc[64];
      vector<string> labels;
      vector<string> xlabels;
      vector<string> ylabels;
      vector<float> bp_min;
      vector<float> bp_max;
      vector<float> bp_median;
      vector<float> bp_medhi;
      vector<float> bp_medlo;
      vector<float> bp_wishi;
      vector<float> bp_wislo;
      vector<float> bp_iqr;
      vector<float> bp_q1;
      vector<float> bp_q3;
      vector<int> bp_outlierindex;
      vector<float> bp_outliervalue;
      int scat_ctr;

      string explabel;
      int req_xsize,req_ysize;
      double null_shift;
      double minmaxscale;
      float ymin,ymax;
			string Ylabel_fmt;
			
      int spacing;
      int bp_colctr;
      int histogram_bins;

      bool bp_notched;
			bool gridswapdefault;
      float bp_whiskerlength, bp_maxall, bp_minall;

      gdImagePtr outim;
      void add_legend(void* ptr, unsigned long cmap[], bool inside=false);
  };
}
#endif
