#if !defined(HISTOGRAMFNS_H)
#define HISTOGRAMFNS_H

#include <valarray>

template <class T>
int find_histogram(const std::valarray<T>& vol, std::valarray<int>& hist, unsigned int bins, 
		   T& min, T& max);

template <class T>
void find_thresholds(const std::valarray<T>& vol, std::valarray<int>& hist, unsigned int bins, 
		     T& minval, T& maxval);

#endif
