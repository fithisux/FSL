/* GDCHART 0.11.2dev  ARRAY_ALLOC.C  Mar 2002 */
/* Copyright Bruce Verderaime 2002 */
/* Contributed in part or whole by Don Fasen and Darren Kent */

#if !defined (__GNUC__)  || defined (__STRICT_ANSI__)

#include <stdio.h>

#include "array_alloc.h"

/*F***************************************************************************
*
* Function Name:         allocate_array_float
*
* Function Prototype:    float	** allocate_array_float (int nr, int nc)
*
* Allocates memory for a two-dimensional array of floats
* referenced as arr[r][c]
*
***************************************************************************
*
*   access array as arr[i1=0-nr][i2=0-nc] = *(&arr[0][0]+nr*i1+i2)
*                       row      col  makes fewest pointer array elements
*                                     assuming Ncol>>Nrow
*  Ex: NC=2
*           _______                  __________
*    **arr |arr[0] |------>arr[0]-->|arr[0][0] |
*          |_______|                |__________|
*          |arr[1] |-----           |arr[0][1] |
*          |_______|     |          |__________|
*          |  .    |      --------->|arr[1][0] |
*          |  .    |                |__________|
*          |_______|                |arr[1][1] |
*          |arr[NR-1]
*          |       |
*         pointer array             data array
*
******************************************************************************
*
* Parameters :
*
* Name               Type           Use       Description
* ----               ----           ---       -----------
* nr                 int          Input     Number of rows
* nc                 int          Input     Number of columns
*
******************************************************************************
*
* Function Return:   float** with allocated memory all set up as float[][].
*
******************************************************************************
*
* Author:  Don Fasen
* Date  :  12/18/98
* Modified: Darren Kent 20/03/2002
*
****************************************************************************F*/
float **
allocate_array_float( int	nr,		/* number of pointer array elements */
					  int	nc )	/* larger number for min memory usage */
{
	int		i;
	float	**arr;

	/* Try to allocate memory for the pointer array, 'nr' pointers-to-float: */
	arr = (float**)malloc( nr * sizeof(float*) );
	/* no test for success ala glib, i.e., if you're outa memory, ... */
#ifdef DBUG
	if( arr == NULL )
		{
		perror( "Unable to malloc pointer array in allocate_array_float()\n" );
		exit( 2 );
		}
#endif
	
	/* Now allocate memory for the data array, 'nr * nc' floats: */
	arr[0] = (float*)malloc( nr * nc * sizeof(float) );
#ifdef DBUG
	if( arr == NULL )
		{
		perror( "Unable to malloc data array in allocate_array_float()\n" );
		exit( 2 );
		}
#endif
	/* fill in the array pointers to point to successive columns. */
	/* Don't need to do the first (zeroth) one, because it was set up by the data malloc. */
	/* Note that adding 1 to (float*)a.arr[0] adds sizeof(float*) automatically. */

	for (i = 1; i < nr; i++) 
		arr[i] = arr[0] + i*nc;

	return arr;

} /* End of FUNCTION allocate_array_float */

/***********************************************************************************
*
* Same Implementation for a Char array.
*
************************************************************************************/
char **
allocate_array_char( int	nr,
					 int	nc )
{
	int		i;
	char	**arr;

	/* Try to allocate memory for the pointer array, 'nr' pointers-to-chars: */
	arr = (char**)malloc( nr * sizeof(char*) );
	/* no test for success ala glib, i.e., if you're outa memory, ... */
#ifdef DBUG
	if( arr == NULL )
		{
		perror( "Unable to malloc pointer array in allocate_array_char()\n" );
		exit( 2 );
		}
#endif
	
	/* Now allocate memory for the data array, 'nr * nc' char: */
	arr[0] = (char*)malloc( nr * nc * sizeof(char) );
#ifdef DBUG
	if( arr == NULL )
		{
		perror( "Unable to malloc data array in allocate_array_char()\n" );
		exit( 2 );
		}
#endif
	/* fill in the array pointers to point to successive columns. */
	/* Don't need to do the first (zeroth) one, because it was set up by the data malloc. */
	/* Note that adding 1 to (char*)a.arr[0] adds sizeof(char*) automatically. */

	for (i = 1; i < nr; i++) 
		arr[i] = arr[0] + i*nc;

	return arr;

} /* End of FUNCTION allocate_array_char */

/***********************************************************************************
*
*  Same Implementation for a Int array.
*
************************************************************************************/
int **
allocate_array_int( int	nr,
					int	nc )

{
	int		i;
	int		**arr;

	/* Try to allocate memory for the pointer array, 'nr' pointers-to-ints: */
	arr = (int**)malloc( nr * sizeof(int*) );
	/* no test for success ala glib, i.e., if you're outa memory, ... */
#ifdef DBUG
	if( arr == NULL )
		{
		perror( "Unable to malloc pointer array in allocate_array_int()\n" );
		exit( 2 );
		}
#endif
	
	/* Now allocate memory for the data array, 'nr * nc' int: */
	arr[0] = (int*)malloc( nr * nc * sizeof(int) );
#ifdef DBUG
	if( arr == NULL )
		{
		perror( "Unable to malloc data array in allocate_array_int()\n" );
		exit( 2 );
		}
#endif
	/* fill in the array pointers to point to successive columns. */
	/* Don't need to do the first (zeroth) one, because it was set up by the data malloc. */
	/* Note that adding 1 to (int*)a.arr[0] adds sizeof(int*) automatically. */

	for (i = 1; i < nr; i++) 
		arr[i] = arr[0] + i*nc;

	return arr;

} /* End of FUNCTION allocate_array_int */

#endif /* !defined (__GNUC__)  || defined (__STRICT_ANSI__) */
