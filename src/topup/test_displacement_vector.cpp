#include <cstring>
#include "miscmaths/miscmaths.h"
#include "newimage/newimageall.h"
#include "topup_file_io.h"
#include "displacement_vector.h"

using namespace TOPUP;

int main(int argc, char *argv[])
{
  // Read a field

  TopupFileReader          topupfile("/Users/jesper/data/topup/FromAdam/test_out");
  NEWIMAGE::volume<float>  field = topupfile.FieldAsVolume();
  
  // Pick out a suitable column

  DispVec  v;
  v.SetFromColumn(field,26,72);
  
  // Print it out along with two inverses and K matrices
  NEWMAT::ColumnVector  cv = v.GetDisplacements();
  write_ascii_matrix(cv,"Original_field.txt");
  cv = v.GetInverseDisplacements(0.0876);
  write_ascii_matrix(cv,"Inverse.txt");
  cv = v.GetInverseDisplacements(-0.0876);
  write_ascii_matrix(cv,"Inverse_of_neg.txt");
  NEWMAT::Matrix M; //  = v.GetK_Matrix(0.0876);
  // write_ascii_matrix(M,"K_matrix.txt");
  M = v.GetK_Matrix(-0.0876);
  write_ascii_matrix(M,"K_matrix_of_neg.txt");
  
  // All seems to work fine.
 
  exit(EXIT_SUCCESS); 
}
