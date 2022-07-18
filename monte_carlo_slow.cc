#pragma once

#include "sphere_helpers.cc"

// estimate the proportion of the upper hemisphere which is within an angle of phi
// of at least one of the unit vectors in centers
double monteCarloSlow(const std::vector<double3>& centers, const double phi) {
  int nSamples = 10000;
  int sampleHits = 0;
  double thresh = cos(phi);
  std::vector<double3> points = sampleSphere(nSamples);
  for (double3 p : points) {
    for (double3 s : centers) {
      if(dot(p, s) > thresh) {
        ++sampleHits;
        break;
      }
    }
  }
  return sampleHits / (double) nSamples;
}