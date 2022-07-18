#pragma once

#include <functional>
#include <algorithm>

#include "sphere_helpers.cc"
#include "integrator.cc"
#include "interval_union.cc"

// One of the two ways we solve this problem with numerical integration.
// In this approach, we integrate over coelevation (pi/2 - elevation, or the
// separation angle from the north pole) from 0 to pi/2. The integrand is a
// function representing how much of the sphere at that coelevation is contained
// in one of the circles. We evaluate this function for a given coelevation by
// slicing the sphere at that coelevation and searching through all the circles 
// that intersect with that slice, recording which azimuth intervals of the slice
// are inside of circles.

// gets the center with the greatest z-coordinate less than or equal to the given z value,
// or the center with the least z-coordinate if that's still greater than
// the given z value.
// note: sortedCenters is sorted by DECREASING z-coordinate.
size_t getNextCenter(const std::vector<double3>& sortedCenters, const double z) {
  size_t minIdx = 0;
  size_t maxIdx = sortedCenters.size();
  while (maxIdx - minIdx > 1) {
    size_t checkIdx = (maxIdx + minIdx) / 2;
    if (sortedCenters[checkIdx].z < z) {
      maxIdx = checkIdx;
    } else {
      minIdx = checkIdx;
    }
  }
  return minIdx;
}

// the integrand that we integrate over: the length of the latitude line 
double integrandA(double coel, const std::vector<double3>& sortedCenters, const double phi) {
  const double minCoel = coel - phi;
  const double maxCoel = coel + phi;
  // i borked this and don't wanna debug it rn
  //size_t firstCenter = getNextCenter(sortedCenters, minCoel);
  //size_t lastCenter = getNextCenter(sortedCenters, maxCoel);
  size_t firstCenter = 0;
  size_t lastCenter = sortedCenters.size()-1;
  //printf("Checking centers %d to %d\n", firstCenter, lastCenter);
  IntervalUnion intervals;
  
  for (size_t centerIdx = firstCenter; centerIdx <= lastCenter; ++centerIdx) {
    double3 s = sortedCenters[centerIdx];
    double sAz = atan2(s.y, s.x);
    sAz += (sAz < 0) ? 2*pi : 0;
    double sCoel = acos(s.z);
    // the difference in coelevations is larger than the radius of the circle, so there are no intersections
    if (abs(coel - sCoel) >= phi) {
      continue;
    }
    // the sum of coelevations is smaller than the radius of the circle, so the circle contains the entire
    // latitude line
    if (coel + sCoel <= phi) {
      //printf("    IntegrandA eval at %f is %f\n", coel, sin(coel) * 2*pi);
      return sin(coel) * 2*pi;
    }
    // finding the intersection of the latitude line at coel and the circle around s w/ radius phi:
    // consider the point on the latitude line closest to s, call it p. Then, rotate the coordinate
    // system about the z-axis such that s and p lie in the x-z plane; their coordinates are 
    // (sin(sCoel), 0, cos(sCoel)) and (sin(coel), 0, cos(coel)) respectively. Now, consider rotating p
    // about the z-axis by an angle theta: its coordinates become p_theta =
    // (sin(coel)*cos(theta), sin(coel)*sin(theta), cos(coel)). We're looking for the angle theta such
    // that the angle between s and p_theta is phi, meaning dot(s, p_theta) = cos(phi). Thus, we need
    // to solve for theta in: sin(sCoel)*sin(coel)*cos(theta) + cos(sCoel)*cos(coel) = cos(phi).
    
    double cosTheta = (cos(phi) - cos(sCoel) * cos(coel)) / (sin(sCoel) * sin(coel));
    
    // There is guaranteed to be either 0 or 1 solution; we get 0 solutions when the circle doesn't contain
    // any of the latitude line (in which cosTheta will be > 1), or the circle contains all of the
    // latitude line (in which case cosTheta will be < -1). These cases should be covered by the checks above,
    // but floating point arithmetic is evil and I've been burned by this exact situation before.
    if (cosTheta >= 1) {
      continue;
    }
    if (cosTheta <= -1) {
      //printf("    IntegrandA eval at %f is %f\n", coel, sin(coel) * 2*pi);
      return sin(coel) * 2*pi;
    }
    double theta = acos(cosTheta);
    
    double minAz = sAz - theta;
    double maxAz = sAz + theta;
    
    intervals.addInterval(minAz, maxAz);
  }
  //printf("    IntegrandA eval at %f is %f\n", coel, sin(coel) * intervals.getLength());
  return sin(coel) * intervals.getLength();
  //return sin(coel);
}

// estimate the proportion of the upper hemisphere which is within an angle of phi
// of at least one of the unit vectors in centers
double integralA(const std::vector<double3>& centers, const double phi) {
  using namespace std::placeholders;
  std::vector<double3> sortedCenters(centers);
  // sort by DECREASING z-coordinate to get increasing coelevation
  std::sort(sortedCenters.begin(), sortedCenters.end(),
    [](double3 left, double3 right) {
      return left.z >= right.z;
    });
  
  Integrator theIntegrator;
  auto integrand = std::bind(integrandA, _1, sortedCenters, phi);
  double result = theIntegrator.integrate(integrand, 0, pi/2) / (2*pi);
  return result;
}