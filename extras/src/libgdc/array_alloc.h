#ifndef _ARRAY_ALLOC_H
#define _ARRAY_ALLOC_H

/* ----- no dynamic auto array, need to mallocate space ----- */
#if !defined (__GNUC__)  || defined (__STRICT_ANSI__)

#include <stdlib.h>

float	**allocate_array_float( int nrs, int ncs );
char	**allocate_array_char( int nrs, int ncs );
int		**allocate_array_int( int nrs, int ncs );

#define CREATE_ARRAY1( arr, type, nels )		type *arr = (type*)malloc( (nels) * sizeof(type) )
												/* array2 can't do ptr (*) types */
#define CREATE_ARRAY2( arr, type, nrs, ncs )	type **arr = allocate_array_##type( nrs, ncs )

#define FREE_ARRAY1( arr )						if( arr ) { free( (void*)arr ); arr = NULL; } else
#define FREE_ARRAY2( arr )						if( arr )						\
													{							\
													if( arr[0] )				\
														free( (void*)arr[0] );	\
													free( (void*)arr );			\
													arr = NULL;					\
													}							\
												else

/* ----- avoid alloc, let the stack do the work ----- */
#else

#define CREATE_ARRAY1( arr, type, nels )		type arr[ nels ]
#define CREATE_ARRAY2( arr, type, nrs, ncs )	type arr[ nrs ][ ncs ]
#define FREE_ARRAY1( arr )
#define FREE_ARRAY2( arr )

#endif

#endif /*_ARRAY_ALLOC_H */
