#include "histogramfns.hpp"

template<>
int find_histogram(const std::valarray<unsigned char>&, std::valarray<int>&, int,
		   unsigned char&, unsigned char&);
template<>
int find_histogram(const std::valarray<short>&, std::valarray<int>&, int,
		   short&, short&);
template<>
int find_histogram(const std::valarray<int>&, std::valarray<int>&, int,
		   int&, int&);
template<>
int find_histogram(std::valarray<float> const&, std::valarray<int>&, int,
		   float&, float&);
template<>
int find_histogram(const std::valarray<double>&, std::valarray<int>&, int,
		   double&, double&);

template <>
void find_thresholds(const std::valarray<unsigned char>&, std::valarray<int>&, int, 
		     unsigned char&, unsigned char&);
template <>
void find_thresholds(const std::valarray<short>&, std::valarray<int>&, int, 
		     short&, short&);
template <>
void find_thresholds(const std::valarray<int>&, std::valarray<int>&, int, 
		     int&, int&);
template <>
void find_thresholds<float>(const std::valarray<float>&, std::valarray<int>&, int, 
		     float&, float&);
template <>
void find_thresholds(const std::valarray<double>&, std::valarray<int>&, int, 
		     double&, double&);

void meter(float& f)
{
  f = 1;
}
