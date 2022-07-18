#pragma once

#include <algorithm>
#include "sphere_helpers.cc"

// indexing function that roughly groups together nearby points
double indexOf(const double3 p) {
  double elevation = asin(p.z);
  double azimuth = atan2(p.y, p.x);
  int nStrips = 8;
  double stripWidth = (pi/2) / nStrips;
  stripWidth += 0.00001; // avoid roundoff shenanigans
  
  return (floor(elevation / stripWidth) + (azimuth + pi) / (2*pi)) / nStrips;
}

bool indexCompare(const double3 left, const double3 right) {
  return indexOf(left) < indexOf(right);
}

size_t findInList(const double3 p, const std::vector<double3>& centers) {
  double pIndex = indexOf(p);
  size_t minIdx = 0;
  size_t maxIdx = centers.size();
  while (maxIdx - minIdx > 1) {
    size_t checkIdx = (maxIdx + minIdx) / 2;
    if (pIndex < indexOf(centers[checkIdx])) {
      maxIdx = checkIdx;
    } else {
      minIdx = checkIdx;
    }
  }
  return minIdx;
}

// estimate the proportion of the upper hemisphere which is within an angle of phi
// of at least one of the unit vectors in centers
double monteCarloSingleIndex(const std::vector<double3>& centers, const double phi) {
  std::vector<double3> c(centers);
  std::sort(c.begin(), c.end(), indexCompare);
  int nSamples = 10000;
  int sampleHits = 0;
  double thresh = cos(phi);
  std::vector<double3> points = sampleSphere(nSamples);
  //std::sort(points.begin(), points.end(), indexCompare);
  for (double3 p : points) {
    // find the closest center after p (measured along the index)
    size_t startingIndex = findInList(p, c);
    if (dot(p, c[startingIndex]) > thresh) {
      ++sampleHits;
      continue;
    }
    // iterate outward from the closest center
    for (int step = 0; step < c.size(); ++step) {
      if ((startingIndex - step >= 0 && dot(p, c[startingIndex - step]) > thresh) ||
          (startingIndex + step < c.size() && dot(p, c[startingIndex + step]) > thresh)) {
        ++sampleHits;
        break;
      }
    }
  }
  return sampleHits / (double) nSamples;
}
