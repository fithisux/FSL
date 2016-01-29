
/*  Copyright (C) 1999-2014 University of Oxford  */

/*  CCOPYRIGHT */

#include "options.h"
#include <fstream>

namespace Utilities {

  using namespace std;

  BaseOption * OptionParser::find_matching_option(const string& optstr)
  {
    for(Options::iterator o = options_.begin();
	o != options_.end();
	++o)
      if((*o)->matches(optstr))
	return *o;

    return 0;
  }

  unsigned int OptionParser::parse_option(const string& optstr, const string& valstr,
					  char *argv[], int valpos, int argc)
    throw(X_OptionError)
  {
    BaseOption * theOption = 0;

    if((theOption = find_matching_option(optstr)) == 0)
      throw X_OptionError(optstr, "Option doesn't exist");

    if(theOption->unset() || (overWriteMode_==Allow)) 
      {
	if(theOption->has_arg()) {
	  if(valstr.length() > 0) {
	    if(theOption->set_value(valstr,argv,valpos,argc))
	      return 1 + theOption->nrequired();
	    else {
	      string errstr = "Couldn't set_value! valstr=\"" + valstr;
	      for (int nn=valpos+1; nn<=valpos + theOption->nrequired(); nn++) {
		if (nn<argc)  errstr += " " + string(argv[nn]);
	      }
	      throw X_OptionError(optstr, errstr + "\""); 
	    }
	  } else if(!theOption->optional()) {
	    throw X_OptionError(optstr, "Missing non-optional argument");
	  }
	}
	if(theOption->optional()) 
	  theOption->use_default_value();
	else
	  theOption->set_value(string());
	return 1;
      } 
    else 
      {
	if( overWriteMode_!= Ignore)
	  throw X_OptionError(optstr, "Option already set");
	else
	  return 1;
      }

    throw X_OptionError(optstr);
    return 0;
  }


  unsigned int OptionParser::parse_long_option(const string& str)
  {
    string optstr(str);
    string valstr;

    string::size_type pos = 0;
    if((pos = str.find("=", 0)) != string::npos) {
      optstr = str.substr(0, pos);
      valstr = str.substr(pos + 1, str.length() - pos + 1);
    }

    parse_option(optstr, valstr, 0,0,0);

    return 1;
  }

  unsigned int OptionParser::parse_config_file(const string& filename)
  {
    ifstream cf(filename.c_str());

    if(cf.fail())
      throw X_OptionError(filename, "Couldn't open the file");
    
    OverwriteMode oldMode=overWriteMode_;
    overWriteMode_=Ignore;

    string optstr; char buffer[1024];
    while (cf >> optstr) {
      if(optstr[0] == '#')
	cf.getline(buffer, 1024);	     // Read and discard the rest of this line
      else if(optstr.substr(0,2) == "--")
	parse_long_option(optstr); // Parse a long option
      else {
	cf.getline(buffer, 1024);
	parse_option(optstr, string(buffer), 0, 0, 0);
      }
    }
    overWriteMode_=oldMode;
    return 1;
  }
 
  unsigned int OptionParser::parse_command_line(unsigned int argc, 
						char **argv, int skip, bool silentFail) 
  {
    unsigned int optpos = 1 + skip;
    unsigned int valpos = 1 + skip;

    while(optpos < argc) {

      unsigned int increments = 0;
      
      string optstr(argv[optpos]), valstr;
   
      if(optstr[0] != '-') {	// End of parsable options
	
	if (!silentFail) // Throwing an error if an unrecognised token occurs in the command line
	  throw X_OptionError(optstr, " is an unrecognised token");
	else // if silentFail is TRUE then ignores an unrecognised token and continues
	  break;
      }

      if(optstr[1] == '-') {	// Parse a long opt

	increments = parse_long_option(optstr);
	optpos += increments;

      } else {

	valpos = optpos + 1;

	for(unsigned int i = 1; i < optstr.length(); ++i)
	  {
	    string suboptstr = "-" + optstr.substr(i, 1);
	    
	    if (valpos<argc) valstr=string(argv[valpos]); else valstr=string();
	    increments = parse_option(suboptstr, valstr, argv, valpos, argc);
	    
	    valpos += increments - 1;
	  }
	
	optpos = valpos;
      }
    } 
    return optpos;		// User should process any remaining args
  }

  std::ostream& operator<<(std::ostream& os, const OptionParser p) 
  {
    for(OptionParser::Options::const_iterator o = p.options_.begin();
	o != p.options_.end(); ++o)
      os << *(*o) << std::endl;

    return os;
  }
}
