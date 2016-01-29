#if !defined(__possumfns_h)
#define __possumfns_h

#include "newmat.h"
#include "newimage/newimageall.h"
#include "newimage/costfns.h"

using namespace NEWMAT;
using namespace NEWIMAGE;

class PMatrix {
private:
  double* dmat;
  float * fmat;
  int rows;
  int cols;
  void initialize(int nrows, int ncols);
public:
  PMatrix();
  PMatrix(int nrows, int ncols);
  ~PMatrix();
  void ReSize(int nrows, int ncols);
  void destroy();
  PMatrix& operator=(const PMatrix& src);
  PMatrix& operator=(double val);
  double& time(int r);
  double time(int r) const;
  float& operator()(int r, int c);
  float& at(int r, int c);
  float operator()(int r, int c) const;
  float at(int r, int c) const;
  bool bounds_check(int r, int c) const;
  int Nrows() const { return rows; }
  int Ncols() const { return cols; }
};

int write_binary_matrix(const PMatrix& mat, const string& filename);
int read_binary_matrix(PMatrix& mres, const string& filename);
int read_binary_matrix(PMatrix& mres, const string& filename,
                       int row1, int row2, int col1, int col2);

void voxel4(const double x,const double y,const double z, 
            const RowVector& tissue,const PMatrix& H,const int nreadp,const int v,
            const double xdim,const double ydim,const double zdim,
	    const double* b0time, const double* b0xtime,const double* b0ytime,const double* b0ztime,
            const double* b0timecourse,const int Nb0,
	    const double b0, const double b0x,const double b0y,const double b0z, 
            const double* timecourse,const double* activation,const int Nact,
	    const string outputname, const double* table_slcprof, const double dslcp, const double dslcp_first, const int Nslc,
            const double den,const double RFtrans, const int opt_test,
            const int nospeedup,
            const int save_kcoord,
            double* sreal, double* simag);

void voxel3(const double x,const double y,const double z, 
	    const RowVector& tissue,const PMatrix& H,int const segA, const int nrf,const int nreadp,const int nonzero,const int v,
            const double xdim,const double ydim,const double zdim,
            const double b1,const double b2,const double b3,const double b4,const double b5,const double b6,const double b7,const double b8,const double b9,
            const double bx1,const double bx2,const double bx3,const double bx4,const double bx5,const double bx6,const double bx7,const double bx8,const double bx9,
            const double by1,const double by2,const double by3,const double by4,const double by5,const double by6,const double by7,const double by8,const double by9,
            const double bz1,const double bz2,const double bz3,const double bz4,const double bz5,const double bz6,const double bz7,const double bz8,const double bz9, 
            const double* timecourse,const double* activation,const int Nact, const string outputname,
            const double* table_slcprof,  const double dslcp, const double dslcp_first, const int Nslc,
	    const double den, const double RFtrans, const int opt_test,
            const int nospeedup,
            const int save_kcoord,const bool rfavg,
            double* sreal, double* simag);

void voxel2(const double x,const double y,const double z, 
	    const RowVector& tissue,const PMatrix& H,const int nrf,const int nreadp, const int v,
            const double xdim,const double ydim,const double zdim,
            const double b0, const double b0gxx,const double b0gyy,const double b0gzz,                    
            const double* timecourse,const double* activation,const int Nact, const string outputname,
            const double* table_slcprof, const double dslcp, const double dslcp_first, const int Nslc,
            const double den, const double RFtrans, const int opt_test,
            const int nospeedup,
            const int save_kcoord,
            double* sreal, double* simag);

void voxel1(const double x,const double y,const double z, 
            const RowVector& tissue,const PMatrix& H,const int nreadp,const int v,
            const double xdim,const double ydim,const double zdim,
            const double b0, const double b0gxx,const double b0gyy,const double b0gzz, 
            const double* timecourse,const double* activation,const int Nact,
	    const string outputname, const double* table_slcprof, const double dslcp, const double dslcp_first, const int Nslc,
            const double den,const double RFtrans, const int opt_test,
            const int nospeedup,
            const int save_kcoord,
            double* sreal, double* simag);


int calc_gradientsROI(volume<double>& b, volume<double>& b0gx, volume<double>& b0gy, volume<double>& b0gz, 
                      const int myid, const int Nxx, const int numprocs);

int calc_gradients4DROI(volume4D<double>& b, volume4D<double>& b0gx, volume4D<double>& b0gy, volume4D<double>& b0gz, 
                        const int myid, const int Nxx, const int numprocs);

void sorter(PMatrix& pulse, const Matrix& motion, const string pulsefile,const int opt_test);

Matrix sorter_old(const Matrix& pulse, const Matrix& motion);

double mj_sinc(const double x);

double round_ivana(const double x, const int n);


#endif

//#ifndef gammabar
//#define 42.58
//#endif
