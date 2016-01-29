All the templates in this directory are slightly modified versions of those in
the IML++ library (math.nist.gov/iml++/) kindly made available by Prof Roldan
Pozo. A complete description of the routines can be found at the web-address
above.

The reason for the small modifications is tiny differences between the API
assumed by IML++, and that of NEWMAT. I had a choice between creating a 
wrapper for NEWMAT, or make changes in the IML++ templates, and I choose
the latter.

The changes are:

1. IML++ assumes zero-offset random access to matrices and vectors. So
e.g. the first element of vector x is x(0). This is in contrast with
NEWMAT which has a one-offset access. And to avoid confusion I also
opted for one-offset access when writing the sparse matrix class SpMat.
So, I have added 1 to all vector/matrix indicies.

2. IML++ assumes the existence of a dot(x,y) global function for 
calculating the dot-product of vectors x and y. And NEWMAT provides
a DotProduct(x,y) global function. So, references to dot has been
changed to DotProduct.

3. IML++ assumes the existence of a global norm(x) function for 
calculating the L2-norm of a vector x. And NEWMAT provides
an x.FrobeniusNorm() function. 

Jesper Andersson FMRIB Image Analysis Group
 
