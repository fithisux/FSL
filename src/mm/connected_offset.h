/*  connected_offset.h

    Mark Woolrich, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  COPYRIGHT  */

#if !defined(connected_offset_h)
#define connected_offset_h

#include "connected_offset.h"

namespace Mm{
    
  class Connected_Offset
    {	
    public:
      Connected_Offset(int px, int py, int pz, int pind, int ppartner_ind)
	: x(px),y(py),z(pz),ind(pind),partner_ind(ppartner_ind)
	{
	}
      
      int x;
      int y;
      int z;

      // index for connected offset
      // for use in storing neighbourhoods as vectors.
      int ind;

      // index for connected offset that would link back to this one
      // in the neighbourhood of the voxel that this connected offset
      // corresponds to.
      // for use in storing neighbourhoods as vectors.
      int partner_ind;
    };
}
#endif
