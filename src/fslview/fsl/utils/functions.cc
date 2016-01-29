/*  Copyright (C) 1999-2004 University of Oxford  */

/*  CCOPYRIGHT */

#include "options.h"

namespace Utilities {
  
  using namespace std;

  template<> string Option<bool>::config_key() const
  {
    if(set()) {
      string key(long_form());
      if( key == "" )
	key = short_form();
      
      return key;
    } else
      return "";
  } 
  
  template<> string Option<bool>::value_string() const { return ""; }

  template<> bool Option<bool>::set_value(const string& s)
  {
    if(s.length() == 0)
      {
	value_ = !default_;
	unset_=false;
      }
    else if (s == "true")
      {
	value_ = true;
	unset_=false;
      }
    else if (s == "false")
      {
	value_ = false;
	unset_=false;
      }
    return !unset_;
  }

  template<> ostream& Option<bool>::print(ostream& os) const 
  {
    os << "# " << help_text() << endl;
    if(set())
      os << config_key().substr(0, config_key().find("=")); 

    return os;
  }
  
  ostream& operator<<(ostream& os, const BaseOption& o)
  {
    return o.print(os);
  }

  bool string_to_T(bool& b, const string& s) {
    b = false;
    return false;
  }

  bool string_to_T(string& d, const string& s) {
    d = s;
    return true;
  }

  bool string_to_T(int& i, const string& s) {
    char *endptr = 0; const char *str = s.c_str();
    i = strtol(str, &endptr, 0);
    if(*endptr == str[s.length()])
      return true;
    else
      return false;
  }

  bool string_to_T(float& v, const string& s) {
    char *endptr = 0; const char *str = s.c_str();
    v = strtod(str, &endptr);
    if(*endptr == str[s.length()])
      return true;
    else
      return false;
  }

  bool string_to_T(vector<int>& vi, const string& s) {
    string str(s), delin(",");
    if(str.find(":")!=string::npos)
      delin = ":";
    str=str+delin;
    vi.clear();
    while(str.size()) {
      int v = atoi(str.substr(0,str.find(delin)).c_str());
      vi.push_back(v);
      str = str.substr(str.find(delin)+1,str.length()-str.find(delin)-1);
    }
    return true;
  }

  bool string_to_T(vector<float>& vi, const string& s) {
    string str(s), delin(",");
    if(str.find(":")!=string::npos)
      delin = ":";
    str=str+delin;
    vi.clear();
    while(str.size()) {
      float v = atof(str.substr(0,str.find(delin)).c_str());
      vi.push_back(v);
      str = str.substr(str.find(delin)+1,str.length()-str.find(delin)-1);
    }
    return true;
  }

  bool string_to_T(vector<string>& vi, const string& s) {
    string str(s), delin(",");
    if(str.find(":")!=string::npos)
      delin = ":";
    str=str+delin;
    vi.clear();
    while(str.size()) {
      string v = str.substr(0,str.find(delin));
      vi.push_back(v);
      str = str.substr(str.find(delin)+1,str.length()-str.find(delin)-1);
    }
    return true;
  }

//   ostream& operator<<(ostream &os, const BaseOption& o) {
//     string test=o.help_text();
//     if ((test.length()>=1) && (test[0]=='~')) {
//       test[0]=' ';
//       return os << "\t" << o.key() << test;
//     } else {
//       return os << "\t" << o.key() << "\t" << o.help_text();
//     }
//   }

  void BaseOption::usage(ostream& os) const {
    string test(help_text());
     if ((test.length()>=1) && (test[0]=='~')) {
       test[0]=' ';
       os << "\t" << key() << test;
     } else {
       os << "\t" << key() << "\t" << help_text();
     }
   }

  bool is_short_form(const string& s)
  {
    return (s.substr(0,2) != "--");
  }


  /*
    @return first short-form key (if any)
  */
  const string BaseOption::short_form() const
  {
    string::size_type pos(0), np;
    
    while( (np = key_.find(",", pos)) != string::npos ) {
      string candidate(key_.substr(pos, np - pos));
      if( is_short_form(candidate) )
	return candidate;
      else
	pos = np + 1;
    }
    string candidate(key_.substr(pos, np - pos));
    if( is_short_form(candidate) )
      return candidate;
    else
      return "";
  }

  /*
    @return first long-form key (if any)
  */
  const string BaseOption::long_form() const
  {
    string::size_type pos(0), np;
    
    while( (np = key_.find(",", pos)) != string::npos ) {
      string candidate(key_.substr(pos, np - pos));

      if( !is_short_form(candidate) )
	return candidate;
      else
	pos = np + 1;
    }
    string candidate(key_.substr(pos, np - pos));
    if( !is_short_form(candidate) )
      return candidate;
    else
      return "";
  }
}
