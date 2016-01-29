#include "histogramfns.h"

template <class T>
int find_histogram(const std::valarray<T>& vol, std::valarray<int>& hist, unsigned int bins, 
		   T& min, T& max)
{
  // size and zero the histogram
  hist.resize(bins); hist = 0;

  if(min == max) { min = vol.min(); max = vol.max(); }

  int validsize(-1);

  if(min != max) {
    double fA = bins / double(max - min);
    double fB = (bins * -min) / double(max - min);
    
    validsize = 0;

    for(unsigned int i = 0; i < vol.size(); ++i) {
      unsigned int idx = unsigned(fA * vol[i] + fB);
      ++hist[ std::max(unsigned(0), std::min(idx, bins - 1)) ];
      ++validsize;
    }      
  }

  return validsize;
}

template <class T>
void find_thresholds(const std::valarray<T>& vol, std::valarray<int>& hist, unsigned int bins, 
		     T& minval, T& maxval)
{
  const unsigned int max_passes(10);

  unsigned int pass(1);
  unsigned int lowest_bin(0), highest_bin(bins-1), bottom_bin(0), top_bin(0);
  T min(vol.min()), max(vol.max());
  T thresh2(0), thresh98(0);

  while((pass == 1) ||
	(double(thresh98 - thresh2) < (double(max - min) / 10))) {

    if(pass > 1) {
      // increase range slightly from the 2-98% range found
      bottom_bin = std::max(int(bottom_bin) - 1,             0);
      top_bin    = std::min(int(top_bin)    + 1, int(bins) - 1);

      double fA = ((max - min) / double(bins));
      T tmpmin = min + ( bottom_bin * fA);
      max      = min + ((top_bin+1) * fA);
      min = tmpmin;
    }

    if(pass == max_passes) { min = max = max_passes; } // give up and revert to full range ...

    int validsize = find_histogram(vol, hist, bins, min, max);

    if(validsize < 1) {
      minval = thresh2 = min; maxval = thresh98 = max;
      return;
    }

    if(pass == max_passes) { // ... _but_ ignore end bins
      validsize -= hist[lowest_bin] + hist[highest_bin];
      ++lowest_bin; --highest_bin;
    }

    if (validsize < 0) {
      // ie zero range
      thresh2=thresh98=min;
      break;
    }

    double fA = ((max-min)/double(bins));

    int count;

    for(count=0, bottom_bin=lowest_bin; count<validsize/50; ++bottom_bin)
      count += hist[bottom_bin];
    --bottom_bin;
    thresh2 =  min + (bottom_bin*fA);

    for(count=0, top_bin=highest_bin; count<validsize/50; --top_bin)
      count += hist[top_bin]; 
    ++top_bin;
    thresh98 = min + ((top_bin+1)*fA);

    if (pass==max_passes)
      break;

    ++pass;
  }

  find_histogram(vol, hist, bins, thresh2, thresh98);

  minval =  thresh2;
  maxval = thresh98;
}
