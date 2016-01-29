#define WANT_STREAM

#include "include.h"
#include "newran.h"

#ifdef use_namespace
using namespace NEWRAN;
#endif

void test1(int);
void test2(int);
void test3(int);
void test4(int);

main()
{

   Random::Set(0.46875);

   Real* s1; Real* s2; Real* s3; Real* s4;
   cout << "\nBegin test\n";   // Forces cout to allocate memory at beginning
   { s1 = new Real[8000]; delete [] s1; }
   { s3 = new Real; delete s3;}
   {

      Real* A = new Real[3750];

      long n = 200000;
      long n_large = 1000000;

      test1(n);
      test2(n);
      test3(n_large);
      test4(n);

      cout << "\nEnd of tests\n";

      delete [] A;
   }

   { s2 = new Real[8000]; delete [] s2; }
   cout << "\n(The following memory checks are probably not valid with all\n";
   cout << "compilers - see documentation)\n";
   cout << "\nChecking for lost memory: "
      << (unsigned long)s1 << " " << (unsigned long)s2 << " ";
   if (s1 != s2) cout << " - error\n"; else cout << " - ok\n\n";
   { s4 = new Real; delete s4;}
   cout << "\nChecking for lost memory: "
      << (unsigned long)s3 << " " << (unsigned long)s4 << " ";
   if (s3 != s4) cout << " - error\n"; else cout << " - ok\n\n";

   return 0;
}

