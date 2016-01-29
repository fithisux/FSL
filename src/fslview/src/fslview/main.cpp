/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if defined(WIN32) 
#pragma warning (disable:4786)
#endif
 
//#include "qwindowsstyle.h"

#include "version.h"
#include "application.h"
#include "splashscreen.h"

#include <qapplication.h>

namespace Utilities {

  bool string_to_T(std::pair<float,float> &p, const string& s) {
    string str(s), delin(",");
    std::vector<float> vf(0);
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

bool string_to_T(std::vector<string> &sl, const string& s) {
  string str(s), delin(",");
  str = str + delin;
  sl.clear();
  while(str.size()) {
    string ss = str.substr(0,str.find(delin));
    sl.push_back(ss);
    str = str.substr(str.find(delin)+1,str.length()-str.find(delin)-1);
  }
  return (sl.size() > 0);
}

int main( int argc, char **argv )
{
  using namespace Utilities;

  Option<bool> verbose(string("-V,--verbose"), false, 
		       string("switch on diagnostic messages"), 
		       false, no_argument);
  Option<bool> help(string("-h,--help"), false,
		    string("display this message"),
		    false, no_argument);
  Option< std::vector<string> > mode(string("-m,--mode"), std::vector<string>(),
				     string("Initial viewer mode. Comma separated list of: 3d; single, ortho; lightbox"), false,
				     requires_argument);
  
  string title("fslview ("+string(Version)+"."+string(Release)+")\n\nCopyright(c) 2005, University of Oxford\nDave Flitney");
  string usage("fslview [-m 3d|ortho|lightbox] <baseimage> [-l lutname] [-b low,hi]"
	       "\n\t[ <overlay> [-l lutname] [-b low,hi] ] ..."
	       "\nfslview -m ortho,lightbox filtered_func_data thresh_zstat1 -t 0.5 thresh_zstat2 -l \"Cool\" -t 0.5");
  
  OptionParser options(title, usage);
  options.add(verbose);
  options.add(help);
  options.add(mode);

  OptionParser imageOptions("Per-image options", "image [-l GreyScale] [-t 0.1] [-b 2.3,6]");
  
  try {
    QApplication::setColorSpec( QApplication::CustomColor );
    QApplication app(argc,argv);			

   //OverlayOptionList overlays;
    ApplicationOptions appOpts;

    for(unsigned int pos = options.parse_command_line(qApp->argc(), qApp->argv());
	int(pos) < qApp->argc(); ) {
      Option<string> lutname(string("-l,--lut"), string("Unset"),
			     string("Lookup table name. As per GUI, one of: Greyscale;"
				    "\n\t\t\t\"Red-Yellow\"; \"Blue-Lightblue\"; Red; Green;"
				    "\n\t\t\tBlue; Yellow; Pink; Hot; Cool; Copper, etc."), 
			     false, requires_argument);
      Option<float> transparency(string("-t,--trans"), float(0.0), 
				 string("Initial transparency, e.g., 0.2"), 
				 false, requires_argument);
      Option< std::pair<float,float> > ibricon(string("-b,--bricon"), std::pair<float,float>(0.0,0.0), 
					       string("Initial bricon range, e.g., 2.3,6"), 
					       false, requires_argument);
  
  
      imageOptions.add(lutname);
      imageOptions.add(ibricon);
      imageOptions.add(transparency);

      // Should be an image name followed by image sub options
      string filename(qApp->argv()[pos]);

      pos += imageOptions.parse_command_line(qApp->argc() - pos, &(qApp->argv()[pos]));

      if(!imageOptions.check_compulsory_arguments())
	imageOptions.usage();

      appOpts.push_back(OverlayOption(filename, lutname, transparency, ibricon));
    }

    if(mode.set()) {
      appOpts.setModes(mode.value());
    }

    // Dummy - unused - options needed to make help text work properly
    Option<string> lutname(string("-l,--lut"), string("Unset"),
			   string("Lookup table name. As per GUI, one of: Greyscale;"
				  "\n\t\t\t\"Red-Yellow\"; \"Blue-Lightblue\"; Red; Green;"
				  "\n\t\t\tBlue; Yellow; Pink; Hot; Cool; Copper, etc."), 
			   false, requires_argument);
    Option<float> transparency(string("-t,--trans"), float(0.0), 
			       string("Initial transparency, e.g., 0.2"), 
			       false, requires_argument);
    Option< std::pair<float,float> > ibricon(string("-b,--bricon"), std::pair<float,float>(0.0,0.0), 
					     string("Initial bricon range, e.g., 2.3,6"), 
					     false, requires_argument);
  
  
    imageOptions.add(lutname);
    imageOptions.add(ibricon);
    imageOptions.add(transparency);
    
    if(help.value() || !options.check_compulsory_arguments()) {
      options.usage();
      imageOptions.brief_usage();
    } else {
      app.connect( &app, SIGNAL(lastWindowClosed()), &app, SLOT(quit()) );
  
      SplashScreen *s = new SplashScreen(0, appOpts);
      s->show();

      return app.exec();
    }

  } catch(X_OptionError& e) {
    //     options.usage();
    cerr << e.what() << endl;
  } catch(std::exception &e) {
    cerr << e.what() << endl;
    options.usage();
    imageOptions.brief_usage();
  } catch (...) {
    cerr << "Unhandled exception!" << endl;
  }

  return -1;
}
