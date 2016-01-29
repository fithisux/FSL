/*  Copyright (C) 1999-2004 University of Oxford  */
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

#include <vector>
#include <string>
#include <fstream>

using namespace std;

namespace Utilities {
  bool string_to_T(pair<float,float> &p, const string& s)
  {
    string str(s), delin(",");
    vector<float> vf(0);
    if(str.find(":")!=string::npos)
      delin = ":";
    str=str+delin;
    vf.clear();
    while(str.size()) {
      float v = atof(str.substr(0,str.find(delin)).c_str());
      vf.push_back(v);
      str = str.substr(str.find(delin)+1,str.length()-str.find(delin)-1);
    }
    bool retval(false);
    if(vf.size() == 2) {
      p.first = vf[0];
      p.second = vf[1];
      retval = true;
    }
    return true;
  }
}

#include "options.h"

using namespace Utilities;

Option<bool> verbose(string("-V,--verbose"), false, 
		     string("switch on diagnostic messages"), 
		     false, no_argument);
Option<bool> debugging(string("-D"), false, 
		     string("switch on debugging mode"), 
		     false, no_argument);
Option<bool> help(string("-h,--help"), false,
		  string("display this message"),
		  false, no_argument);
Option<float> dof(string("-d,--dof"), 100.0,
		  string("number of degrees of freedom"),
		  true, requires_argument);
Option<string> mask(string("-m,--mask"), string("mask"),
		    string("brain mask volume"),
		    true, requires_argument);
Option<string> resid(string("-r,--res"), string("res4d"),
		     string("4d `residual-of-fit' image"),
		     true, requires_argument);
Option<string> config_file(string("-c,--config"), string(""),
			   string("Specify a config file to read the default settings from."),
			   false, requires_argument);
Option<int> segid(string("-s,--shared-seg-id"), -1,
		  "shared memory segment ID",
		  false, requires_argument);
HiddenOption<bool> noint(string("-n,--no-scientific-integrity"), false,
			 string("You complete putz"),
			 false, no_argument);
Option<string> zopt(string("--zopt"), string("Whoo!"),
		     string("string input"),
		     false, optional_argument);
Option< pair<float,float> > popt(string("-P,--popt"), std::make_pair(0.0, 0.0),
				 string("X,Y location"),
				 false, requires_argument);

FmribOption< std::vector<string> > 
strseq(string("-I"), std::vector<string>(),
       string("A coma seperated include path"),
       false, requires_argument);

string title = 
"opttst (Version 2.0)\n\n\
Copyright(c) 2000-2007, University of Oxford\n\
Author: Dave Flitney";

string examples =
"opttst --dof=<number> --mask=<filename> --res=<filename>\n\
opttst -d <number> -m <filename> -r <filename>\n\
opttst --verbose\n";

int main(int argc, char **argv) {

  OptionParser options(title, examples);

  try {

    options.add(verbose);
    options.add(debugging);
    options.add(help);
    options.add(config_file);
    options.add(segid);
    options.add(dof);
    options.add(mask);
    options.add(resid);
    options.add(noint);
    options.add(strseq);
    options.add(zopt);
    options.add(popt);

    for(unsigned int a = options.parse_command_line(argc, argv); 
	a < argc; ) {
      // Should be image names followed by optional image options

      string imagename(argv[a]);
      
      // Possible sub-options as follows:
      Option<string> lutname(string("-l,--lut"), string("Unset"), 
			     string("Lookup table name. One of: GreyScale; RedYellow; BlueLightblue; Red; Green; Blue, etc."), 
			     false, requires_argument);
      Option< std::pair<float,float> >
	ibricon(string("-b,--bricon"), std::pair<float,float>(), string("Initial bricon range, e.g., -1:2.5"), false, requires_argument);

      OptionParser imageOptions("", "image -l GreyScale -b 2.3,6");

      imageOptions.add(lutname);
      imageOptions.add(ibricon);

      //      ++a;
      a += imageOptions.parse_command_line(argc - a, &(argv[a])) ;

      cout << "imagename = " << imagename << endl;
      cout << "lutname   = " << lutname.value() << endl;
      cout << "ibricon   = " << ibricon.value().first << ", " << ibricon.value().second<< endl;
    }

    if(config_file.set())
      options.parse_config_file(config_file.value());

    if(help.value() || 
       !options.check_compulsory_arguments(true))
      options.usage();


    if(verbose.value()) {
      cout << "verbose = " << verbose.value() << endl;
      cout << "help = " << help.value() << endl;
      cout << "segid = " << segid.value() << endl;
      cout << "dof = " << dof.value() << endl;
      dof.set_T(50);
      cout << "dof.set_T = " << dof.value() << endl;
      cout << "mask = " << mask.value() << endl;
      cout << "resid = " << resid.value() << endl;
      cout << "noint = " << noint.value() << endl;

      if(zopt.set()) cout << "zopt = " << zopt.value() << endl;
      if(config_file.set()) cout << "config_file = " << config_file.value() << endl;

      for(int i =0; i < (int)strseq.value().size(); i++)
	cout << strseq.value().at(i) << endl;

      cout << "popt   = " << popt.value().first << ", " << popt.value().second<< endl;
      cout << endl << endl;
    } else {
     for(int i =0; i < (int)strseq.value().size(); i++)
	cout << strseq.value().at(i) << endl;      
    }

    std::ofstream of("saved_config");

    of << options << endl;

//     cerr << verbose << endl;
//     cerr << debugging << endl;
//     cerr << help << endl;
//     cerr << segid << endl;
//     cerr << dof << endl;
//     cerr << mask << endl;
//     cerr << resid << endl;
//     cerr << noint << endl;
//     cerr << zopt << endl;
//     cerr << strseq << endl;
//     cerr << popt << endl;

  } catch(X_OptionError& e) {
    options.usage();
    cerr << endl << "Exception:: " << e.what() << endl;
  } catch(std::exception &e) {
    cerr << e.what() << endl;
  }    
}
