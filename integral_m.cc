#pragma once

#include <functional>

#include "sphere_helpers.cc"
#include "integrator.cc"
#include "interval_union.cc"
#include "intersection_finder.cc"

// precondition: the intersection exists (halfAngle < phi), halfAngle > 0
double angleToIntersection(const double halfAngle, const double phi) {
  // same derivation as in integral_a.cc; by choosing coordinates such that the center of one circle is
  // at the north pole, we see that this is a special case where phi = coel, and 2*halfAngle is sCoel.
  double cosTheta = cos(phi) * (1 - cos(2*halfAngle)) / ( sin(2*halfAngle) * sin(phi) );
  // cos(phi) > 1 iff there's no intersection: when halfAngle = phi, using sin(2x) = 2sin(x)cos(x) and 
  // cos(2x) = cos^2(x) - sin^2(x), we get 2*cos(phi)*sin^2(phi) / (2*cos(phi)*sin^2(phi)) = 1
  // here, we assume this code is only called when the intersection exists
  return acos(cosTheta);
}

double getContributionToContour(const double contourRadius, const IntersectionStructure& intersections, const std::vector<double3>& centers, const double phi, const size_t centerIdx) {
  auto centerIntersectionsEntry = intersections.find(centerIdx);
  IntervalUnion intersectionBlocks;
  if (centerIntersectionsEntry != intersections.end()) {
    auto centerIntersections = centerIntersectionsEntry->second;
    for (IntersectionEntry ixn: centerIntersections) {
      if (ixn.halfAngle >= contourRadius) {
        continue;
      }
      double deltaAngle = angleToIntersection(ixn.halfAngle, contourRadius);
      intersectionBlocks.addInterval(ixn.direction - deltaAngle, ixn.direction + deltaAngle);
    }
  }
  
  // check intersection with equator
  double el = asin(centers[centerIdx].z);
  if (el < contourRadius) {
    double deltaAngle = (el > 0.000001) ? angleToIntersection(el, contourRadius) : pi/2;
    // 0 points north (except at the north pole where it points to (1,0,0),
    // so pi points south, which is where the interval is centered.
    //printf("    hit equator! center %d has elevation %f at contourRadius %f, covering an angle of %f\n", centerIdx, el, contourRadius, 2*deltaAngle);
    intersectionBlocks.addInterval(pi - deltaAngle, pi + deltaAngle);
  }
  return 2*pi - intersectionBlocks.getLength();
}

double integrandM(const double contourRadius, const IntersectionStructure& intersections, const std::vector<double3>& centers, const double phi) {
  double result = 0;
  for (size_t centerIdx = 0; centerIdx < centers.size(); ++centerIdx) {
    result += getContributionToContour(contourRadius, intersections, centers, phi, centerIdx);
  }
  return result * sin(contourRadius);
}

double integralM(const std::vector<double3>& centers, const double phi) {
  using namespace std::placeholders;
  IntersectionStructure intersections = findIntersections(centers, phi);
  
  Integrator theIntegrator;
  auto integrand = std::bind(integrandM, _1, intersections, centers, phi);
  double result = theIntegrator.integrate(integrand, 0, phi) / (2*pi);
  return result;
}