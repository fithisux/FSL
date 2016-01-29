// newran.cpp -----------------------------------------------------------

// NEWRAN02

#define WANT_STREAM
#define WANT_MATH
#include "include.h"

#include "newran.h"
//#include "mother.h"

#ifdef use_namespace
namespace NEWRAN {
#endif

//********* classes which are used internally only **********************

//*********** chi-square-1 random number generator **********************

class ChiSq1 : public Normal              // generate non-central chi-square
					  // rv with 1 degree of freedom
{
   Real deltasq;                          // non-centrality parameter
   Real delta;

public:
   ChiSq1(Real);                          // non-centrality parameter
   ExtReal Mean() const { return 1.0+deltasq; }
   ExtReal Variance() const { return 2.0+4.0*deltasq; }
   Real Next();
};

//*********** Poisson random number generator, larger mu ****************

class Poisson1 : public AsymGen           // generate poisson rv, large mu
{
   Real mu, ln_mu;
public:
   Poisson1(Real);                        // constructor (Real=mean)
   Real Density(Real) const;              // Poisson density function
   Real Next() { return floor(AsymGen::Next()); }
   ExtReal Mean() const { return mu; }
   ExtReal Variance() const { return mu; }
};

//*********** Poisson random number generator, mu <= 10 ****************

class Poisson2 : public Random            // generate poisson rv, large mu
{
   DiscreteGen* dg;
public:
   Poisson2(Real);                        // constructor (Real=mean)
   ~Poisson2();
   Real Next() { return dg->Next(); }
   ExtReal Mean() const { return dg->Mean(); }
   ExtReal Variance() const { return dg->Variance(); }
};

//********** Gamma random number generator, alpha <= 1 *****************

class Gamma1 : public PosGen              // generate gamma rv
					  // shape parameter <= 1
{
   Real ln_gam, ralpha, alpha;
public:
   Gamma1(Real);                          // constructor (Real=shape)
   Real Density(Real) const;              // gamma density function
   Real Next();                           // carries out power transform
   ExtReal Mean() const { return alpha; }
   ExtReal Variance() const { return alpha; }
};

//********** Gamma random number generator, alpha > 1 ******************

class Gamma2 : public AsymGen             // generate gamma rv
					  // shape parameter > 1
{
   Real alpha, ln_gam;
public:
   Gamma2(Real);                          // constructor (Real=shape)
   Real Density(Real) const;              // gamma density function
   ExtReal Mean() const { return alpha; }
   ExtReal Variance() const { return alpha; }
};

//*********** Binomial random number generator, n > 40 *****************

class Binomial1 : public AsymGen           // generate binomial rv, larger n
{
   Real p, q, ln_p, ln_q, ln_n_fac; int n;
public:
   Binomial1(int nx, Real px);
   Real Density(Real) const;
   Real Next() { return floor(AsymGen::Next()); }
   ExtReal Mean() const { return p * n; }
   ExtReal Variance() const { return p * q * n; }
};

//******* Binomial random number generator, n < 40 or n*p <= 8 *************

class Binomial2 : public Random            // generate binomial rv, smaller n
{
   DiscreteGen* dg;
public:
   Binomial2(int nx, Real px);
   ~Binomial2();
   Real Next() { return dg->Next(); }
   ExtReal Mean() const { return dg->Mean(); }
   ExtReal Variance() const { return dg->Variance(); }
};

//************************ static variables ***************************

double Random::seed;
//unsigned long Random::iseed;                // for Mother
Real Random::Buffer[128];
Real Normal::Nxi;
Real* Normal::Nsx;
Real* Normal::Nsfx;
long Normal::count=0;

//**************************** utilities ******************************

inline Real square(Real x) { return x*x; }
inline ExtReal square(const ExtReal& x) { return x*x; }

static void ErrorNoSpace() { Throw(Bad_alloc("Newran: out of space")); }

//************************* end of definitions ************************


Real Random::Raw()                           // get new uniform random number
{
   // m = 2147483647 = 2^31 - 1; a = 16807;
   // 127773 = m div a; 2836 = m mod a
   long iseed = (long)seed;
   long hi = iseed / 127773L;                 // integer division
   long lo = iseed - hi * 127773L;            // modulo
   iseed = 16807 * lo - 2836 * hi;
   if (iseed <= 0) iseed += 2147483647L;
   seed = (double)iseed; return seed*4.656612875e-10;
}

Real Random::Density(Real) const
{ Throw(Logic_error("density function not defined")); return 0.0; }

#ifdef _MSC_VER
static void DoNothing(int) {}
#endif

Real Random::Next()                          // get new mixed random number
{
   if (!seed)
      Throw(Logic_error("Random number generator not initialised"));
   int i = (int)(Raw()*128);               // 0 <= i < 128
#ifdef _MSC_VER
   DoNothing(i); DoNothing(i);
#endif
   Real f = Buffer[i]; Buffer[i] = Raw();  // Microsoft release gets this wrong
   return f;

   // return Mother(&iseed);
}

double Random::Get()                  // get random number seed
{ return seed/2147483648L; }

void Random::Set(double s)            // set random number seed
                                      // s must be between 0 and 1
{
   if (s>=1.0 || s<=0.0)
      Throw(Logic_error("Newran: seed out of range"));
   //iseed = 2147483648L * s;         // for Mother
   seed = (long)(s*2147483648L);
   for (int i = 0; i<128; i++) Buffer[i] = Raw();
}


PosGen::PosGen()                             // Constructor
{
   #ifdef MONITOR
      cout << "constructing PosGen\n";
   #endif
   NotReady=true;
}

PosGen::~PosGen()
{
   if (!NotReady)
   {
      #ifdef MONITOR
	 cout << "freeing PosGen arrays\n";
      #endif
      delete [] sx; delete [] sfx;
   }
   #ifdef MONITOR
      cout << "destructing PosGen\n";
   #endif
}

void PosGen::Build(bool sym)                 // set up arrays
{
   #ifdef MONITOR
      cout << "building PosGen arrays\n";
   #endif
   int i;
   NotReady=false;
   sx=new Real[60]; sfx=new Real[60];
   if (!sx || !sfx) ErrorNoSpace();
   Real sxi=0.0; Real inc = sym ? 0.01 : 0.02;
   for (i=0; i<60; i++)
   {
      sx[i]=sxi; Real f1=Density(sxi); sfx[i]=f1;
      if (f1<=0.0) goto L20;
      sxi+=inc/f1;
   }
   Throw(Runtime_error("Newran: area too large"));
L20:
   if (i<50) Throw(Runtime_error("Newran: area too small"));
   xi = sym ? 2*i : i;
   return;
}

Real PosGen::Next()
{
   Real ak,y; int ir;
   if (NotReady) Build(false);
   do
   {
      Real r1=Random::Next();
      ir = (int)(r1*xi); Real sxi=sx[ir];
      ak=sxi+(sx[ir+1]-sxi)*Random::Next();
      y=sfx[ir]*Random::Next();
   }
   while ( y>=sfx[ir+1] && y>=Density(ak) );
   return ak;
}

Real SymGen::Next()
{
   Real s,ak,y; int ir;
   if (NotReady) Build(true);
   do
   {
      s=1.0;
      Real r1=Random::Next();
      if (r1>0.5) { s=-1.0; r1=1.0-r1; }
      ir = (int)(r1*xi); Real sxi=sx[ir];
      ak=sxi+(sx[ir+1]-sxi)*Random::Next();
      y=sfx[ir]*Random::Next();
   }
   while ( y>=sfx[ir+1] && y>=Density(ak) );
   return s*ak;
}

AsymGen::AsymGen(Real modex)                 // Constructor
{
   #ifdef MONITOR
      cout << "constructing AsymGen\n";
   #endif
   mode=modex; NotReady=true;
}

void AsymGen::Build()                        // set up arrays
{
   #ifdef MONITOR
      cout << "building AsymGen arrays\n";
   #endif
   int i;
   NotReady=false;
   sx=new Real[121]; sfx=new Real[121];
   if (!sx || !sfx)  ErrorNoSpace();
   Real sxi=mode;
   for (i=0; i<120; i++)
   {
      sx[i]=sxi; Real f1=Density(sxi); sfx[i]=f1;
      if (f1<=0.0) goto L20;
      sxi+=0.01/f1;
   }
   Throw(Runtime_error("Newran: area too large (a)"));
L20:
   ic=i-1; sx[120]=sxi; sfx[120]=0.0;
   sxi=mode;
   for (; i<120; i++)
   {
      sx[i]=sxi; Real f1=Density(sxi); sfx[i]=f1;
      if (f1<=0.0) goto L30;
      sxi-=0.01/f1;
   }
   Throw(Runtime_error("Newran: area too large (b)"));
L30:
   if (i<100)  Throw(Runtime_error("Newran: area too small"));
   xi=i;
   return;
}

Real AsymGen::Next()
{
   Real ak,y; int ir1;
   if (NotReady) Build();
   do
   {
      Real r1=Random::Next();
      int ir=(int)(r1*xi); Real sxi=sx[ir];
      ir1 = (ir==ic) ? 120 : ir+1;
      ak=sxi+(sx[ir1]-sxi)*Random::Next();
      y=sfx[ir]*Random::Next();
   }
   while ( y>=sfx[ir1] && y>=Density(ak) );
   return ak;
}

AsymGen::~AsymGen()
{
   if (!NotReady)
   {
      #ifdef MONITOR
	 cout << "freeing AsymGen arrays\n";
      #endif
      delete [] sx; delete [] sfx;
   }
   #ifdef MONITOR
      cout << "destructing AsymGen\n";
   #endif
}

PosGenX::PosGenX(PDF fx) { f=fx; }

SymGenX::SymGenX(PDF fx) { f=fx; }

AsymGenX::AsymGenX(PDF fx, Real mx) : AsymGen(mx) { f=fx; }


Normal::Normal()
{
   if (count) { NotReady=false; xi=Nxi; sx=Nsx; sfx=Nsfx; }
   else { Build(true); Nxi=xi; Nsx=sx; Nsfx=sfx; }
   count++;
}

Normal::~Normal()
{
   count--;
   if (count) NotReady=true;                     // disable freeing arrays
}

Real Normal::Density(Real x) const               // normal density
{ return (fabs(x)>8.0) ? 0 : 0.398942280 * exp(-x*x / 2); }

ChiSq1::ChiSq1(Real d) : Normal()                // chisquare with 1 df
{ deltasq=d; delta=sqrt(d); }

Real ChiSq1::Next()
{ Real z=Normal::Next()+delta; return z*z; }

ChiSq::ChiSq(int df, Real noncen)
{
  if (df<=0 || noncen<0.0) Throw(Logic_error("Newran: illegal parameters"));
  else if (df==1) { version=0; c1=new ChiSq1(noncen); }
  else if (noncen==0.0)
  {
     if (df==2) { version=1; c1=new Exponential(); }
     else { version=2; c1=new Gamma2(0.5*df); }
  }
  else if (df==2) { version=3; c1=new ChiSq1(noncen/2.0); }
  else if (df==3) { version=4; c1=new Exponential(); c2=new ChiSq1(noncen); }
  else { version=5; c1=new Gamma2(0.5*(df-1)); c2=new ChiSq1(noncen); }
  if (!c1 || (version>3 && !c2)) ErrorNoSpace();
  mean=df+noncen; var=2*df+4.0*noncen;
}

ChiSq::~ChiSq() { delete c1; if (version>3) delete c2; }

Real ChiSq::Next()
{
   switch(version)
   {
   case 0: return c1->Next();
   case 1: case 2: return c1->Next()*2.0;
   case 3: return c1->Next() + c1->Next();
   case 4: case 5: Real s1 = c1->Next()*2.0; Real s2 = c2->Next();
	   return s1 + s2; // this is to make it work with Microsoft VC5
   }
   return 0;
}

Pareto::Pareto(Real shape) : Shape(shape)
{
   if (Shape <= 0) Throw(Logic_error("Newran: illegal parameter"));
   RS = -1.0 / Shape;
}

Real Pareto::Next()
{ return pow(Random::Next(), RS); }

ExtReal Pareto::Mean() const
{
   if (Shape > 1) return Shape/(Shape-1.0);
   else return PlusInfinity;
}

ExtReal Pareto::Variance() const
{
   if (Shape > 2) return Shape/(square(Shape-1.0))/(Shape-2.0);
   else return PlusInfinity;
}

Real Cauchy::Density(Real x) const               // Cauchy density function
{ return (fabs(x)>1.0e15) ? 0 : 0.31830988618 / (1.0+x*x); }

Poisson1::Poisson1(Real mux) : AsymGen(mux)      // Constructor
{ mu=mux; ln_mu=log(mu); }

Real Poisson1::Density(Real x) const             // Poisson density function
{
   if (x < 0.0) return 0.0;
   double ix = floor(x);                         // use integer part
   double l = ln_mu * ix - mu - ln_gamma(1.0 + ix);
   return  (l < -40.0) ? 0.0 : exp(l);
}

Binomial1::Binomial1(int nx, Real px)
   : AsymGen((nx + 1) * px), p(px), q(1.0 - px), n(nx)
      { ln_p = log(p); ln_q = log(q); ln_n_fac = ln_gamma(n+1); }

Real Binomial1::Density(Real x) const            // Binomial density function
{
   double ix = floor(x);                         // use integer part
   if (ix < 0.0 || ix > n) return 0.0;
   double l = ln_n_fac - ln_gamma(ix+1) - ln_gamma(n-ix+1)
      + ix * ln_p + (n-ix) * ln_q;
   return  (l < -40.0) ? 0.0 : exp(l);
}

Poisson2::Poisson2(Real mux)
{
   Real probs[40];
   probs[0]=exp(-mux);
   for (int i=1; i<40; i++) probs[i]=probs[i-1]*mux/i;
   dg=new DiscreteGen(40,probs);
   if (!dg) ErrorNoSpace();
}

Poisson2::~Poisson2() { delete dg; }

Binomial2::Binomial2(int nx, Real px)
{
   Real qx = 1.0 - px;
   Real probs[40];
   int k = (int)(nx * px);
   probs[k] = exp(ln_gamma(nx+1) - ln_gamma(k+1) - ln_gamma(nx-k+1)
      + k * log(px) + (nx-k) * log(qx));
   int i;
   int m = (nx >= 40) ? 39 : nx;
   for (i=k+1; i<=m; i++) probs[i]=probs[i-1] * px * (nx-i+1) / qx / i;
   for (i=k-1; i>=0; i--) probs[i]=probs[i+1] * qx * (i+1) / px / (nx-i);
   dg = new DiscreteGen(m + 1, probs);
   if (!dg) ErrorNoSpace();
}

Binomial2::~Binomial2() { delete dg; }

Real Exponential::Density(Real x) const          // Negative exponential
{ return  (x > 40.0 || x < 0.0) ? 0.0 : exp(-x); }

Poisson::Poisson(Real mu)
{
   if (mu <= 8.0) method = new Poisson2(mu);
   else method = new Poisson1(mu);
   if (!method) ErrorNoSpace();
}

Binomial::Binomial(int nx, Real px)
{
   if (nx < 40 || nx * px <= 8.0) method = new Binomial2(nx, px);
   else method = new Binomial1(nx, px);
   if (!method) ErrorNoSpace();
}

NegativeBinomial::NegativeBinomial(Real NX, Real PX)
   : AsymGen(0.0), N(NX), P(PX), Q(1.0 + PX)
{
   p = 1.0 / Q;  ln_q = log(1.0 - p);
   c = N * log(p) - ln_gamma(N);  mode = (N - 1) * P;
   if (mode < 1.0) mode = 0.0;
}

Real NegativeBinomial::Next() { return floor(AsymGen::Next()); }

Real NegativeBinomial::Density(Real x) const
{
   if (x < 0.0) return 0.0;
   Real ix = floor(x);
   Real l = c + ln_gamma(ix + N) - ln_gamma(ix + 1) + ix * ln_q;
   return  (l < -40.0) ? 0.0 : exp(l);
}

Gamma1::Gamma1(Real alphax)                      // constructor (Real=shape)
{ ralpha=1.0/alphax; ln_gam=ln_gamma(alphax+1.0); alpha=alphax; }

Real Gamma1::Density(Real x) const               // density function for
{                                                // transformed gamma
   Real l = - pow(x,ralpha) - ln_gam;
   return  (l < -40.0) ? 0.0 : exp(l);
}

Real Gamma1::Next()                               // transform variable
{ return pow(PosGen::Next(),ralpha); }

Gamma2::Gamma2(Real alphax) : AsymGen(alphax-1.0) // constructor (Real=shape)
{ alpha=alphax; ln_gam=ln_gamma(alpha); }

Real Gamma2::Density(Real x) const                // gamma density function
{
   if (x<=0.0) return 0.0;
   double l = (alpha-1.0)*log(x) - x - ln_gam;
   return  (l < -40.0) ? 0.0 : exp(l);
}

Gamma::Gamma(Real alpha)                         // general gamma generator
{
   if (alpha<1.0) method = new Gamma1(alpha);
   else if (alpha==1.0) method = new Exponential();
   else method = new Gamma2(alpha);
   if (!method)  ErrorNoSpace();
}

DiscreteGen::DiscreteGen(int n1, Real* prob)     // discrete generator
						 // values on 0,...,n1-1
{
   #ifdef MONITOR
      cout << "constructing DiscreteGen\n";
   #endif
   Gen(n1, prob); val=0;
   mean=0.0; var=0.0;
   { for (int i=0; i<n; i++) mean = mean + i*prob[i]; }
   { for (int i=0; i<n; i++) var = var + square(i-mean) * prob[i]; }
}

DiscreteGen::DiscreteGen(int n1, Real* prob, Real* val1)
                                                 // discrete generator
                                                 // values on *val
{
   #ifdef MONITOR
      cout << "constructing DiscreteGen\n";
   #endif
   Gen(n1, prob); val = new Real[n1];
   if (!val)  ErrorNoSpace();
   for (int i=0; i<n1; i++) val[i]=val1[i];
   mean=0.0; var=0.0;
   { for (int i=0; i<n; i++) mean = mean + val[i]*prob[i]; }
   { for (int i=0; i<n; i++) var = var + square(val[i]-mean)*prob[i]; }
}


void DiscreteGen::Gen(int n1, Real* prob)
{
   n=n1;                                         // number of values
   p=new Real[n]; ialt=new int[n];
   if (!p || !ialt)  ErrorNoSpace();
   Real rn = 1.0/n; Real px = 0; int i;
   for (i=0; i<n; i++) { p[i]=0.0; ialt[i]=-1; }
   for (i=0; i<n; i++)
   {
      Real pmin=1.0; Real pmax=-1.0; int jmin=-1; int jmax=-1;
      for (int j=0; j<n; j++)
      {
         if (ialt[j]<0)
         {
            px=prob[j]-p[j];
            if (pmax<=px) { pmax=px; jmax=j; }
            if (pmin>=px) { pmin=px; jmin=j; }
         }
      }
      if ((jmax<0) || (jmin<0)) Throw(Runtime_error("Newran: method fails"));
      ialt[jmin]=jmax; px=rn-pmin; p[jmax]+=px; px*=n; p[jmin]=px;
      if ((px>1.00001)||(px<-.00001))
         Throw(Runtime_error("Newran: probs don't add to 1 (a)"));
   }
   if (px>0.00001) Throw(Runtime_error("Newran: probs don't add to 1 (b)"));
}

DiscreteGen::~DiscreteGen()
{
   delete [] p; delete [] ialt; delete [] val;
   #ifdef MONITOR
      cout << "destructing DiscreteGen\n";
   #endif
}

Real DiscreteGen::Next()                  // Next discrete random variable
{
   int i = (int)(n*Random::Next()); if (Random::Next()<p[i]) i=ialt[i];
   return val ? val[i] : (Real)i;
}

Real ln_gamma(Real xx)
{
   // log gamma function adapted from numerical recipes in C

   if (xx<1.0)                           // Use reflection formula
   {
      double piz = 3.14159265359 * (1.0-xx);
      return log(piz/sin(piz))-ln_gamma(2.0-xx);
   }
   else
   {
      static double cof[6]={76.18009173,-86.50532033,24.01409822,
         -1.231739516,0.120858003e-2,-0.536382e-5};

      double x=xx-1.0; double tmp=x+5.5;
      tmp -= (x+0.5)*log(tmp); double ser=1.0;
      for (int j=0;j<=5;j++) { x += 1.0; ser += cof[j]/x; }
      return -tmp+log(2.50662827465*ser);
   }
}


Real NegatedRandom::Next() { return - rv->Next(); }

ExtReal NegatedRandom::Mean() const { return - rv->Mean(); }

ExtReal NegatedRandom::Variance() const { return rv->Variance(); }

Real ScaledRandom::Next() { return rv->Next() * s; }

ExtReal ScaledRandom::Mean() const { return rv->Mean() * s; }

ExtReal ScaledRandom::Variance() const { return rv->Variance() * (s*s); }

Real ShiftedRandom::Next() { return rv->Next() + s; }

ExtReal ShiftedRandom::Mean() const { return rv->Mean() + s; }

ExtReal ShiftedRandom::Variance() const { return rv->Variance(); }

Real ReverseShiftedRandom::Next() { return s - rv->Next(); }

ExtReal ReverseShiftedRandom::Mean() const { return - rv->Mean() + s; }

ExtReal ReverseShiftedRandom::Variance() const { return rv->Variance(); }

Real ReciprocalRandom::Next() { return s / rv->Next(); }

ExtReal RepeatedRandom::Mean() const { return rv->Mean() * (Real)n; }

ExtReal RepeatedRandom::Variance() const { return rv->Variance() * (Real)n; }

RepeatedRandom& Random::operator()(int n)
{
   RepeatedRandom* r = new RepeatedRandom(*this, n);
   if (!r) ErrorNoSpace(); return *r;
}

NegatedRandom& operator-(Random& rv)
{
   NegatedRandom* r = new NegatedRandom(rv);
   if (!r) ErrorNoSpace(); return *r;
}

ShiftedRandom& operator+(Random& rv, Real s)
{
   ShiftedRandom* r = new ShiftedRandom(rv, s);
   if (!r) ErrorNoSpace(); return *r;
}

ShiftedRandom& operator-(Random& rv, Real s)
{
   ShiftedRandom* r = new ShiftedRandom(rv, -s);
   if (!r) ErrorNoSpace(); return *r;
}

ScaledRandom& operator*(Random& rv, Real s)
{
   ScaledRandom* r = new ScaledRandom(rv, s);
   if (!r) ErrorNoSpace(); return *r;
}

ShiftedRandom& operator+(Real s, Random& rv)
{
   ShiftedRandom* r = new ShiftedRandom(rv, s);
   if (!r) ErrorNoSpace(); return *r;
}

ReverseShiftedRandom& operator-(Real s, Random& rv)
{
   ReverseShiftedRandom* r = new ReverseShiftedRandom(rv, s);
   if (!r) ErrorNoSpace(); return *r;
}

ScaledRandom& operator*(Real s, Random& rv)
{
   ScaledRandom* r = new ScaledRandom(rv, s);
   if (!r) ErrorNoSpace(); return *r;
}

ScaledRandom& operator/(Random& rv, Real s)
{
   ScaledRandom* r = new ScaledRandom(rv, 1.0/s);
   if (!r) ErrorNoSpace(); return *r;
}

ReciprocalRandom& operator/(Real s, Random& rv)
{
   ReciprocalRandom* r = new ReciprocalRandom(rv, s);
   if (!r) ErrorNoSpace(); return *r;
}

AddedRandom& operator+(Random& rv1, Random& rv2)
{
   AddedRandom* r = new AddedRandom(rv1, rv2);
   if (!r) ErrorNoSpace(); return *r;
}

MultipliedRandom& operator*(Random& rv1, Random& rv2)
{
   MultipliedRandom* r = new MultipliedRandom(rv1, rv2);
   if (!r) ErrorNoSpace(); return *r;
}

SubtractedRandom& operator-(Random& rv1, Random& rv2)
{
   SubtractedRandom* r = new SubtractedRandom(rv1, rv2);
   if (!r) ErrorNoSpace(); return *r;
}

DividedRandom& operator/(Random& rv1, Random& rv2)
{
   DividedRandom* r = new DividedRandom(rv1, rv2);
   if (!r) ErrorNoSpace(); return *r;
}

Real AddedRandom::Next() { return rv1->Next() + rv2->Next() ; }

ExtReal AddedRandom::Mean() const { return rv1->Mean() + rv2->Mean() ; }

ExtReal AddedRandom::Variance() const
   { return rv1->Variance() + rv2->Variance() ; }

Real SubtractedRandom::Next() { return rv1->Next() - rv2->Next() ; }

ExtReal SubtractedRandom::Mean() const { return rv1->Mean() - rv2->Mean() ; }

ExtReal SubtractedRandom::Variance() const
   { return rv1->Variance() + rv2->Variance() ; }

Real MultipliedRandom::Next() { return rv1->Next() * rv2->Next() ; }

ExtReal MultipliedRandom::Mean() const { return rv1->Mean() * rv2->Mean() ; }

ExtReal MultipliedRandom::Variance() const
{
   ExtReal e1 = square(rv1->Mean()); ExtReal e2 = square(rv2->Mean());
   ExtReal v1 = rv1->Variance(); ExtReal v2 = rv2->Variance();
   ExtReal r=v1*v2+v1*e2+e1*v2;
   return r;
}

Real DividedRandom::Next() { return rv1->Next() / rv2->Next() ; }

void Random::load(int*,Real*,Random**)
   { Throw(Logic_error("Newran: illegal combination")); }

void SelectedRandom::load(int* i, Real* probs, Random** rvx)
{
   probs[*i]=p; rvx[*i]=rv; (*i)++;
   delete this;
}

Real SelectedRandom::Next()
   { Throw(Logic_error("Newran: Next not defined")); return 0.0; }

Real AddedSelectedRandom::Next()
   { Throw(Logic_error("Newran: Next not defined")); return 0.0; }

Real RepeatedRandom::Next()
   { Real sum=0.0; for (int i=0; i<n; i++) sum+=rv->Next(); return sum; }

MixedRandom::MixedRandom(int nx, Real* probs, Random** rvx)
{
   n = nx;
   rv = new Random*[n]; if (!rv) ErrorNoSpace();
   for (int i=0; i<n; i++) rv[i]=rvx[i];
   Build(probs);
}

MixedRandom::MixedRandom(AddedSelectedRandom& sr)
{
   n = sr.nelems();                       // number of terms;
   Real* probs = new Real[n]; rv = new Random*[n];
   if (!probs || !rv) ErrorNoSpace();
   int i=0; sr.load(&i,probs,rv);
   Build(probs); delete [] probs;
}

void MixedRandom::Build(Real* probs)
{
   int i;
   dg=new DiscreteGen(n,probs);
   if (!dg) ErrorNoSpace();
   mean=0.0; var=0.0;
   for (i=0; i<n; i++) mean = mean + (rv[i])->Mean()*probs[i];
   for (i=0; i<n; i++)
   {
      ExtReal sigsq=(rv[i])->Variance();
      ExtReal mudif=(rv[i])->Mean()-mean;
      var = var + (sigsq+square(mudif))*probs[i];
   }

}

MixedRandom::~MixedRandom()
{
   for (int i=0; i<n; i++) rv[i]->tDelete();
   delete [] rv; delete dg;
}

Real MixedRandom::Next()
   { int i = (int)(dg->Next()); return (rv[i])->Next(); }

int AddedSelectedRandom::nelems() const
   { return rv1->nelems() + rv2->nelems(); }

void AddedSelectedRandom::load(int* i, Real* probs, Random** rvx)
{
   rv1->load(i, probs, rvx); rv2->load(i, probs, rvx);
   delete this;
}

AddedSelectedRandom& operator+(SelectedRandom& rv1,
   SelectedRandom& rv2)
{
   AddedSelectedRandom* r = new AddedSelectedRandom(rv1, rv2);
   if (!r) ErrorNoSpace(); return *r;
}

AddedSelectedRandom& operator+(AddedSelectedRandom& rv1,
   SelectedRandom& rv2)
{
   AddedSelectedRandom* r = new AddedSelectedRandom(rv1, rv2);
   if (!r) ErrorNoSpace(); return *r;
}

AddedSelectedRandom& operator+(SelectedRandom& rv1,
   AddedSelectedRandom& rv2)
{
   AddedSelectedRandom* r = new AddedSelectedRandom(rv1, rv2);
   if (!r) ErrorNoSpace(); return *r;
}

AddedSelectedRandom& operator+(AddedSelectedRandom& rv1,
   AddedSelectedRandom& rv2)
{
   AddedSelectedRandom* r = new AddedSelectedRandom(rv1, rv2);
   if (!r) ErrorNoSpace(); return *r;
}

SelectedRandom& Random::operator()(double p)
{
   SelectedRandom* r = new SelectedRandom(*this, p);
   if (!r) ErrorNoSpace(); return *r;
}




// Identification routines for each class - may not work on all compilers?

const char* Random::Name()            { return "Random";           }
const char* Uniform::Name()           { return "Uniform";          }
const char* Constant::Name()          { return "Constant";         }
const char* PosGen::Name()            { return "PosGen";           }
const char* SymGen::Name()            { return "SymGen";           }
const char* AsymGen::Name()           { return "AsymGen";          }
const char* PosGenX::Name()           { return "PosGenX";          }
const char* SymGenX::Name()           { return "SymGenX";          }
const char* AsymGenX::Name()          { return "AsymGenX";         }
const char* Normal::Name()            { return "Normal";           }
const char* ChiSq::Name()             { return "ChiSq";            }
const char* Cauchy::Name()            { return "Cauchy";           }
const char* Exponential::Name()       { return "Exponential";      }
const char* Poisson::Name()           { return "Poisson";          }
const char* Binomial::Name()          { return "Binomial";         }
const char* NegativeBinomial::Name()  { return "NegativeBinomial"; }
const char* Gamma::Name()             { return "Gamma";            }
const char* Pareto::Name()            { return "Pareto";           }
const char* DiscreteGen::Name()       { return "DiscreteGen";      }
const char* SumRandom::Name()         { return "SumRandom";        }
const char* MixedRandom::Name()       { return "MixedRandom";      }


// ********************** permutation generator ***************************

void RandomPermutation::Next(int N, int M, int p[], int start)
{
   // N = size of urn; M = number of draws
   if (N < M) Throw(Logic_error("Newran: N < M in RandomPermutation"));
   int i;
   int* q = new int [N];
   if (!q) ErrorNoSpace();
   for (i = 0; i < N; i++) q[i] = start + i;
   for (i = 0; i < M; i++)
   {
      int k = i + (int)(U.Next() * (N - i));       // uniform on i ... N-1
      p[i] = q[k]; q[k] = q[i];                    // swap i and k terms
                                                   // but put i-th term into p
   }
   delete [] q;
}

void RandomCombination::SortAscending(int n, int gm[])
{
   // from numerical recipies in C - Shell sort

   const double aln2i = 1.442695022; const double tiny = 1.0e-5;
   int m = n; int lognb2 = (int)(aln2i * log((double)n) + tiny);
   while (lognb2--)
   {
      m >>= 1;
      for (int j = m; j<n; j++)
      {
         int* gmj = gm+j; int i = j-m; int* gmi = gmj-m; int t = *gmj;
         while (i>=0 && *gmi>t)  { *gmj = *gmi; gmj = gmi; gmi -= m; i -= m; }
         *gmj = t;
      }
   }
}

#ifdef use_namespace
}
#endif


