/*  easyoptions.h - FABBER's options-handling class

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

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

#pragma once
#include <map>
#include <string>
#include "easylog.h"

class EasyOptions {
 public:
    // Create an instance from command-line options
    EasyOptions(int argc, char** argv);
        // break down options into key=value clauses

    EasyOptions(const map<string,string>& src, 
                const map<string,const Matrix*>* dataIn, 
                map<string,Matrix>* dataOut) 
        : args(src) // copy key=value pairs from src; use key="" for no-argument options
        { 
          assert(inMatrices==NULL); 
          assert(outMatrices==NULL); 
          inMatrices=dataIn; 
          outMatrices=dataOut; 
          assert(outMatrices->empty());
	  args[""] = "fabber_library"; // This would normally hold argv[0].
        }

    // Below: option-reading values.  Once they are called, 
    // the corresponding key=value pair is removed... this is deliberate, to
    // ensure that eacy argument is used exactly once.  If you want to use an
    // argument in more than one place, you're probably doing something wrong...
    // this should make code more maintainable by making it easy to see where
    // options have an effect on code.
    
    string Read(const string& key)
        { return Read(key, "Missing mandatory option: --" + key + "\n"); }
        // Get a mandatory option.
        // throws if option is missing or has no =value clause
    string Read(const string& key, const string& msg);
        // supply msg if you want a customized (helpful) error message.
        
    bool ReadBool(const string& key);
        // true if present, false if absent, throws error if --key=value.
        
    string ReadWithDefault(const string& key, const string& def);
        // Like default, but if --key=value is omitted it just treats
        // as if --key=default.
        
    void CheckEmpty();
        // throws if there are any options left

    ~EasyOptions() { }; 
        // throwing an exception in a destructor would be a bad idea!
        // Also, note that args may not be empty, if an exception has
        // been thrown.

    friend ostream& operator<<(ostream& out, const EasyOptions& opts);

private:
    map<string,string> args;
        // all remaining unused options.
    void AddKey(const string& key); // internal helper function

    // These are only used in fabber_library mode:
    static const map<string,const Matrix*>* inMatrices;
    static map<string,Matrix>* outMatrices;
public:
    static bool UsingMatrixIO() { assert(!inMatrices == !outMatrices); return (inMatrices!=NULL); }
    static const Matrix& InMatrix(const string& filename);  // simpler syntax
    static Matrix& OutMatrix(const string& filename);
};
 
// Helper function:
Matrix read_vest_fabber(const string& filename);

typedef EasyOptions ArgsType;

// Use NEWMAT's exception handling base class
class Invalid_option : public Runtime_error
{
public:
  static unsigned long Select;          // for identifying exception
  Invalid_option(const char* c) : Runtime_error(c) { return; }
  Invalid_option(const string& s) : Runtime_error(s.c_str()) { return; }
};

// Convert a string into almost anything, using operator>>. See the C++ FAQ-Lite: 
// http://www.parashift.com/c++-faq-lite/misc-technical-issues.html#faq-39.3
template<typename T>
inline T convertTo(const std::string& s)
{
  T x;
  istringstream i(s);
  char c;
  if (!(i>>x) || (i.get(c)))
    throw Invalid_option(
        "convertTo failed... couldn't convert string '" + s + "' into a value.\n");
  return x;
}


