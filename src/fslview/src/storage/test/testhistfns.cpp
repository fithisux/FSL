#include "../histogramfns.hpp"
#include "../image.h"
#include <iostream>

template <typename T>
void printValarray(const std::valarray<T>& va, const unsigned int n)
{
  for(unsigned int i = 0; i < va.size() / n; ++i) {
    for(unsigned int j = 0; j < n; ++j)
      std::cout << va[i * n + j] << ' ';
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

int main(int argc, char * argv[]) 
{
  try{
    Image::Handle         image = Image::load(argv[1]);
    ImageInfo::Handle imageInfo = image->getInfo();
    Volume::Handle       volume = image->getVolume(0);

    std::cout << "BriCon range = " << image->getVolume(0)->inqMin() << " to " << image->getVolume(0)->inqMax() << std::endl;

    unsigned int voxels = volume->inqX() * volume->inqY() * volume->inqZ();
  
    std::valarray<float> data(voxels);

    for(unsigned int i = 0; i < voxels; ++i)
      data[i] = volume->value(i);

    // Now calculate the threshold limits and associated histograms. 
    float min(0), max(0);
    std::valarray<int> histogram;

//     std::cout << "mean = " << data.sum() / float(data.size()) << std::endl;
    std::cout << "min = " << data.min() << ", max = " << data.max() << " - ";
    find_thresholds(data, histogram, 100, min, max);
    std::cout << "2% threshold = " << min << ", 98% threshold = " << max << std::endl;

    printValarray(histogram, 10);

    return 1;

  } catch(...) {
    std::cout << "Oops. It's all gone wrong!" << std::endl; 
    return 0;
  }
}
