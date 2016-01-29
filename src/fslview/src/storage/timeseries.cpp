/*  FSLView - 2D/3D Interactive Image Viewer

    James Saunders, David Flitney and Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include "timeseries.hpp"

template class TimeSeriesStore<char>;
template class TimeSeriesStore<unsigned char>;
template class TimeSeriesStore<short>;
template class TimeSeriesStore<unsigned short>;
template class TimeSeriesStore<int>;
template class TimeSeriesStore<unsigned int>;
template class TimeSeriesStore<long long>;
template class TimeSeriesStore<unsigned long long>;
template class TimeSeriesStore<float>;
template class TimeSeriesStore<double>;
template class TimeSeriesStore<long double>;
