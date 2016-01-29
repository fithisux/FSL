/*  defaultschedule.h

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2008 University of Oxford  */

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

// Set default schedule file
//  Written by Mark Jenkinson  11/10/99
//  Modified by Mark Jenkinson    08/08

#if !defined(__defaultschedule_h)
#define __defaultschedule_h

#include <vector>
#include <string>

void setdefaultschedule(std::vector<string>& comms)
{
  comms.clear();
  comms.push_back("# 8mm scale");
  comms.push_back("setscale 8");
  comms.push_back("setoption smoothing 8");
  comms.push_back("clear S");
  comms.push_back("clear P");
  comms.push_back("search");

  comms.push_back("# 4mm scale");
  comms.push_back("setscale 4");
  comms.push_back("setoption smoothing 4");
  comms.push_back("clear U");
  comms.push_back("clear UA ");
  comms.push_back("clear UB");
  comms.push_back("clear US");
  comms.push_back("clear UP");

  comms.push_back("# remeasure costs at this scale");
  comms.push_back("measurecost 7 S 0 0 0 0 0 0 rel");
  comms.push_back("copy U US");
  comms.push_back("clear U");
  comms.push_back("measurecost 7 P 0 0 0 0 0 0 rel");
  comms.push_back("copy U UP");
  comms.push_back("dualsort US UP");

  comms.push_back("# optimise best 3 candidates (pre and post 8mm optimisations)");
  comms.push_back("clear U");
  comms.push_back("optimise 7 US:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 UP:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("# also try the identity transform as a starting point at this resolution");
  comms.push_back("clear UQ");
  comms.push_back("setrow UQ  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1");
  comms.push_back("optimise 7 UQ  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("sort U");
  comms.push_back("copy U UA");

  comms.push_back("# select best 4 optimised solutions and try perturbations of these");
  comms.push_back("clear U");
  comms.push_back("copy UA:1-4 U");
  comms.push_back("optimise 7 UA:1-4  1.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 UA:1-4 -1.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 UA:1-4  0.0   1.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 UA:1-4  0.0  -1.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 UA:1-4  0.0   0.0   1.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 UA:1-4  0.0   0.0  -1.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 7 UA:1-4  0.0   0.0   0.0   0.0   0.0   0.0   0.1  abs 4");
  comms.push_back("optimise 7 UA:1-4  0.0   0.0   0.0   0.0   0.0   0.0  -0.1  abs 4");
  comms.push_back("optimise 7 UA:1-4  0.0   0.0   0.0   0.0   0.0   0.0   0.2  abs 4");
  comms.push_back("optimise 7 UA:1-4  0.0   0.0   0.0   0.0   0.0   0.0  -0.2  abs 4");
  comms.push_back("sort U");
  comms.push_back("copy U UB");

  comms.push_back("# 2mm scale");
  comms.push_back("setscale 2");
  comms.push_back("setoption smoothing 2");
  comms.push_back("clear U");
  comms.push_back("clear UC");
  comms.push_back("clear UD");
  comms.push_back("clear UE");
  comms.push_back("clear UF");

  comms.push_back("# remeasure costs at this scale");
  comms.push_back("measurecost 7 UB 0 0 0 0 0 0 rel");
  comms.push_back("sort U");
  comms.push_back("copy U UC");

  comms.push_back("clear U");
  comms.push_back("optimise 7  UC:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("copy U UD");
  comms.push_back("setoption boundguess 1");
  comms.push_back("if MAXDOF > 7");
  comms.push_back(" clear U");
  comms.push_back("if MAXDOF > 7");
  comms.push_back(" optimise 9  UD:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1");
  comms.push_back("copy U UE");
  comms.push_back("if MAXDOF > 9");
  comms.push_back(" clear U");
  comms.push_back("if MAXDOF > 9");
  comms.push_back(" optimise 12 UE:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 2");
  comms.push_back("sort U");
  comms.push_back("copy U UF");

  comms.push_back("# 1mm scale");
  comms.push_back("setscale 1");
  comms.push_back("setoption smoothing 1");
  comms.push_back("setoption boundguess 1");
  comms.push_back("clear U");
  comms.push_back("# also try the qsform as a starting point at this resolution");
  comms.push_back("setrowqsform UF");
  comms.push_back("optimise 12 UF:1-2  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1");
  comms.push_back("# in addition, try qsform as the final transformation, not just an initialisation");
  comms.push_back("clear UG");
  comms.push_back("setrowqsform UG");
  comms.push_back("measurecost 12 UG:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1");
  comms.push_back("sort U");
}


void set2Ddefaultschedule(std::vector<string>& comms)
{
  comms.clear();
  comms.push_back("# 8mm scale");
  comms.push_back("setscale 8");
  comms.push_back("setoption smoothing 8");
  comms.push_back("setoption paramsubset 3  0 0 1 0 0 0 0 0 0 0 0 0  0 0 0 1 0 0 0 0 0 0 0 0  0 0 0 0 1 0 0 0 0 0 0 0");
  comms.push_back("clear U");
  comms.push_back("clear UA");
  comms.push_back("setrow UA 1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1");
  comms.push_back("optimise 12 UA:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4 ");
  comms.push_back("");
  comms.push_back("# 4mm scale");
  comms.push_back("setscale 4");
  comms.push_back("setoption smoothing 4");
  comms.push_back("setoption paramsubset 3  0 0 1 0 0 0 0 0 0 0 0 0  0 0 0 1 0 0 0 0 0 0 0 0  0 0 0 0 1 0 0 0 0 0 0 0");
  comms.push_back("clear UB");
  comms.push_back("clear UL");
  comms.push_back("clear UM");
  comms.push_back("# remeasure costs at this scale");
  comms.push_back("clear U");
  comms.push_back("measurecost 12 UA 0 0 0 0 0 0 rel");
  comms.push_back("sort U");
  comms.push_back("copy U UL");
  comms.push_back("# optimise best 3 candidates");
  comms.push_back("clear U");
  comms.push_back("optimise 12 UL:1-3  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("# also try the identity transform as a starting point at this resolution");
  comms.push_back("clear UQ");
  comms.push_back("setrow UQ  1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1");
  comms.push_back("optimise 7 UQ  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("sort U");
  comms.push_back("copy U UM");
  comms.push_back("# select best 4 optimised solutions and try perturbations of these");
  comms.push_back("clear U");
  comms.push_back("copy UM:1-4 U");
  comms.push_back("optimise 12 UM:1-4  0.0   0.0   1.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("optimise 12 UM:1-4  0.0   0.0  -1.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("sort U");
  comms.push_back("clear UB");
  comms.push_back("copy U UB");
  comms.push_back("");
  comms.push_back("# 2mm scale");
  comms.push_back("setscale 2");
  comms.push_back("setoption smoothing 2");
  comms.push_back("setoption paramsubset 3  0 0 1 0 0 0 0 0 0 0 0 0  0 0 0 1 0 0 0 0 0 0 0 0  0 0 0 0 1 0 0 0 0 0 0 0");
  comms.push_back("clear U");
  comms.push_back("clear UC");
  comms.push_back("clear UD");
  comms.push_back("clear UE");
  comms.push_back("clear UF");
  comms.push_back("# remeasure costs at this scale");
  comms.push_back("measurecost 12 UB 0 0 0 0 0 0 rel");
  comms.push_back("sort U");
  comms.push_back("copy U UC");
  comms.push_back("clear U");
  comms.push_back("optimise 12  UC:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 4");
  comms.push_back("copy U UD");
  comms.push_back("setoption boundguess 1");
  comms.push_back("if MAXDOF > 7");
  comms.push_back(" clear U");
  comms.push_back("if MAXDOF > 7");
  comms.push_back(" optimise 9  UD:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1");
  comms.push_back("copy U UE");
  comms.push_back("if MAXDOF > 9");
  comms.push_back(" clear U");
  comms.push_back("if MAXDOF > 9");
  comms.push_back(" optimise 12 UE:1  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 2");
  comms.push_back("sort U");
  comms.push_back("copy U UF");
  comms.push_back("");
  comms.push_back("# 1mm scale");
  comms.push_back("setscale 1");
  comms.push_back("setoption smoothing 1");
  comms.push_back("setoption boundguess 1");
  comms.push_back("setoption paramsubset 3  0 0 1 0 0 0 0 0 0 0 0 0  0 0 0 1 0 0 0 0 0 0 0 0  0 0 0 0 1 0 0 0 0 0 0 0");
  comms.push_back("clear U");
  comms.push_back("# also try the qsform as a starting point at this resolution");
  comms.push_back("setrowqsform UF");
  comms.push_back("optimise 12 UF:1-2  0.0   0.0   0.0   0.0   0.0   0.0   0.0  rel 1");
  comms.push_back("sort U");
}

#endif
