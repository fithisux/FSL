#define WANT_STREAM
#define WANT_MATH

#include "include.h"
#include "boolean.h"
#include "extreal.h"
#include "newran.h"

#ifdef use_namespace
using namespace NEWRAN;
#endif


void Histogram(Random* rx, int n)          // draw histogram with n obsv
{
   int i,j; int count[20];
   Real* a = new Real[n];
   if (!a) { std::cout << "\nNo memory for Histogram\n"; return; }
   for (i = 0; i < n; i++) a[i] = rx->Next();
   Real amax = a[0]; Real amin = a[0]; Real mean = a[0]; Real sd = 0;
   for (i = 1; i < n; i++)
   {
      if (amin > a[i]) amin = a[i]; else if (amax < a[i]) amax = a[i];
      mean += a[i];
   }
   mean /= n;
   for (i = 0; i < 20; i++) count[i]=0;
   for (i = 0; i < n; i++)
   {
     Real rat= (amax != amin) ? (a[i] - amin)/(amax - amin) : 1.0;
     j = (int)( 19.999 * rat ); count[j]++;
     Real diff = a[i] - mean; sd += diff*diff;
   }
   sd = sqrt(sd/(n-1));
   j = 0;
   for (i = 0; i < 20; i++) { if (j < count[i]) j = count[i]; }
   if (j > 70) { for (i = 0; i < 20; i++) count[i] = (int)((70L*count[i])/j); }
   std::cout << "\n";
   for (i = 0; i < 20; i++)
     { std::cout << "\n|"; for (j = 1; j < count[i]; j = j+1) std::cout << "*"; }
   std::cout << "\n" << rx->Name() << "\n";
   std::cout << "p. mean = " << std::setw(9) << std::setprecision(5) << rx->Mean()
	     << ", p. var = " << std::setw(9) << std::setprecision(5)
           << rx->Variance() << "\n";
   std::cout << "s. mean = " << std::setw(9) << std::setprecision(5) << mean
	     << ", s. var = " << std::setw(9) << std::setprecision(5) << sd*sd
	     << ", max = " << std::setw(9) << std::setprecision(5) << amax
	     << ", min = " << std::setw(9) << std::setprecision(5) << amin
        << "\n";
   std::cout << std::flush;
   delete a;
}
