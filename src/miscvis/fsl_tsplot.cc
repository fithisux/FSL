/* {{{ Copyright etc. */

/*  fsl_tsplot - 

    Christian Beckmann, FMRIB Image Analysis Group

    Copyright (C) 2006-2007 University of Oxford  */

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

/* }}} */
/* {{{ defines, includes and typedefs */
 
#include "libvis/miscplot.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"
#include "utils/options.h"
#include <vector>
#include <cfloat> 
#include <limits>
#include <limits.h>
using namespace MISCPLOT;
using namespace MISCMATHS;
using namespace Utilities;
using namespace std;

// The two strings below specify the title and example usage that is
//  printed out as the help or usage message

string title="fsl_tsplot \nCopyright(c) 2007, University of Oxford (Christian F. Beckmann)";
string examples="fsl_tsplot [options] ";

// Each (global) object below specificies as option and can be accessed
//  anywhere in this file (since they are global).  The order of the
//  arguments needed is: name(s) of option, default value, help message,
//       whether it is compulsory, whether it requires arguments
// Note that they must also be included in the main() function or they
//  will not be active.

//Option<bool> verbose(string("-v,--verbose"), false, 
//		     string("switch on diagnostic messages"), 
//		     false, no_argument);
Option<bool> help(string("--help"), false,
		  string("display this message"),
		  false, no_argument);
Option< std::vector<string> >inname(string("-i,--in"), std::vector<string>(),
		  string("        comma-separated list of input file names (ASCII text matrix, one column per timecourse)"),
		  true, requires_argument);
Option<string> outname(string("-o,--out"), string(""),
		  string("output filename for the PNG file"),
		  true, requires_argument);
Option<float>  ymax(string("--ymax"),float(0.0),
		    string("        maximum y-value"),
		    false, requires_argument);
Option<float>  ymin(string("--ymin"),float(0.0),
		    string("        minimum y-value"),
		    false, requires_argument);  
Option<string> ptitle(string("-t,--title"), string(""),
		  string("plot title"),
		  false, requires_argument);
Option<string> xtitle(string("-x,--xlabel"), string(""),
		  string("X-axis label"),
		  false, requires_argument);
Option<string> ytitle(string("-y,--ylabel"), string(""),
		  string("Y-axis label"),
		  false, requires_argument);
Option<int> ysize(string("-h,--height"), 150,
		  string("plot height in pixels (default 150)"),
		  false, requires_argument);
Option<int> xsize(string("-w,--width"), 600,
		  string("plot width in pixels (default 600)"),
		  false, requires_argument);
Option<float> TR(string("-u,--unit"), 0.0,
		  string("scaling units for x-axis (default 1...length of infile"),
		  false, requires_argument);
Option<int> width(string("--numsize"), 4,
		  string("size of individual x-labels (number of digits)"),
		  false, requires_argument);
Option<int> prec(string("--precision"), 2,
		  string("precision of x-axis labels"),
		  false, requires_argument);
Option<bool> sci(string("--sci"), false,
		  string("        switch on sientific notation"),
		  false, no_argument);
Option<int> start(string("--start"), 0,
		  string("        Position of first column to plot"),
		  false, requires_argument);
Option<int> finish(string("--finish"), 0,
		  string("Position of final column to plot"),
		  false, requires_argument);
Option<string> labelname(string("-l,--legend"), string(""),
		  string("file name of ASCII text file, one row per legend entry"),
		  false, requires_argument);
Option< std::vector<string> >labels(string("-a,--labels"), std::vector<string>(),
		  string("comma-separated list of labels"),
		  false, requires_argument);

int nonoptarg;

////////////////////////////////////////////////////////////////////////////

// Local functions
int do_work(int argc, char* argv[]) 
{
  Matrix in;
  in=read_ascii_matrix(inname.value().at(0));
  for (int i = 1; i < (int)inname.value().size(); i++){
    Matrix tmp;
    tmp = read_ascii_matrix(inname.value().at(i));

    //increase tmp or in depending on what is wrong
    int diffrows=in.Nrows()-tmp.Nrows();
    if (diffrows>0)
    {
      Matrix diff(diffrows,tmp.Ncols());
      diff=std::numeric_limits<float>::quiet_NaN();
      tmp &=diff;
    }
    else if (diffrows<0)
    {
      Matrix diff(-diffrows,in.Ncols());
      diff=std::numeric_limits<float>::quiet_NaN();
      in &=diff;
    }

    in |= tmp;
  }
  if (start.value() && finish.value())
    in=in.SubMatrix(1,in.Nrows(),start.value(),finish.value());

  miscplot newplot;

  newplot.set_minmaxscale(1.001);
  if (labelname.value().size()>0)
    {
      ifstream fs(labelname.value().c_str());
      if (!fs) {
	cerr << "Could not open file " << labelname.value() << endl;
      }
      
      int ctr=1;
      string cline;
      while(!fs.eof())
	{
	  getline(fs, cline);
	  if (cline.size()>0 && ctr <= in.Ncols()){
	    newplot.add_label(cline);
	    ctr++;
	  }
	} 
      fs.close();
    }
  for (int i = 0; i < (int)labels.value().size(); i++) newplot.add_label(labels.value().at(i));
  if((xsize.value()>0.0)&&(ysize.value()>0.0))
    newplot.set_xysize(xsize.value(),ysize.value());

  if (ymin.value()==ymax.value()) newplot.set_yrange(in.Minimum(),in.Maximum());
  else newplot.set_yrange(ymin.value(),ymax.value());
  newplot.add_xlabel(xtitle.value());
  newplot.add_ylabel(ytitle.value());  
  newplot.timeseries(in.t(),outname.value(),ptitle.value(),TR.value(),
		     ysize.value(),width.value(),prec.value(),sci.value());

  return 0;
}

////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{

  Tracer tr("main");
  OptionParser options(title, examples);

  try {
    // must include all wanted options here (the order determines how
    //  the help message is printed)
    options.add(inname);
    options.add(outname);
    options.add(ptitle);
    options.add(labelname);
    options.add(labels);
    options.add(ymin);
    options.add(ymax);
    options.add(xtitle);
    options.add(ytitle);
    options.add(ysize);
    options.add(xsize);
    options.add(TR);
    options.add(prec);
    options.add(sci);
    options.add(start);
    options.add(finish);    
    options.parse_command_line(argc, argv);

    // line below stops the program if the help was requested or 
    //  a compulsory option was not set
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) )
      {
	options.usage();
	exit(EXIT_FAILURE);
      }
    
  }  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
                                                     exit(EXIT_FAILURE);
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  } 

  // Call the local functions
  return do_work(argc,argv);
}

