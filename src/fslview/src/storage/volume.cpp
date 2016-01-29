/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "volume.hpp"

Volume& Volume::operator=(const Volume& rhs)
{
  if(this == &rhs)
    return *this;

  for(int z = 0; z < rhs.inqZ(); ++z)
    for(int y = 0; y < rhs.inqY(); ++y)
      for(int x = 0; x < rhs.inqX(); ++x)
	setValue(x, y, z, rhs.value(x, y, z));

  return *this;
}

template class VolumeStore<char>;
template class VolumeStore<unsigned char>;
template class VolumeStore<short>;
template class VolumeStore<unsigned short>;
template class VolumeStore<int>;
template class VolumeStore<unsigned int>;
template class VolumeStore<long long>;
template class VolumeStore<unsigned long long>;
template class VolumeStore<float>;
template class VolumeStore<double>;
template class VolumeStore<long double>;
