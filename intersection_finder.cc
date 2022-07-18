#pragma once

#include <unordered_map>
#include <vector>

#include "sphere_helpers.cc"

struct IntersectionEntry {
  size_t otherIdx;
  double halfAngle;
  double direction;
  
  bool operator<(const IntersectionEntry& other) const {
    return halfAngle < other.halfAngle;
  }
};

IntersectionEntry getIntersectionEntry(const std::vector<double3>& centers, const size_t myIdx, const size_t otherIdx) {
  const double3& me = centers[myIdx];
  const double3& other = centers[otherIdx];
  
  double cosDist = dot(me, other);
  double halfAngle = 0.5 * acos(cosDist);
  double3 mainDirection = (me.z < 0.99999) ? double3{0,0,1} : double3{1,0,0};
  double3 mainDirectionTangent = mainDirection.minusProjectionUnit(me);
  double3 otherDirection = normalize(other - me);
  double3 otherDirectionTangent = otherDirection.minusProjectionUnit(me);
  double direction = acos(dot(mainDirectionTangent, otherDirectionTangent));
  
  // acos measures just the angle between the vectors (between 0 and pi), but we
  // need the angle specifically measured counterclockwise. if the angle was measured
  // clockwise, we actually want 2pi minus that angle.
  double3 clockwiseDirection = cross(me, mainDirectionTangent);
  if (dot(clockwiseDirection, otherDirectionTangent) < 0) {
    direction = 2*pi - direction;
  }
  return IntersectionEntry{otherIdx, halfAngle, direction};
}

using IntersectionStructure = std::unordered_map<size_t, std::vector<IntersectionEntry>>;

IntersectionStructure findIntersections(const std::vector<double3>& centers, double phi) {
  IntersectionStructure result;
  double thresh = cos(2*phi);
  for (size_t i = 0; i < centers.size(); ++i) {
    for (size_t j = i+1; j < centers.size(); ++j) {
      double cosDist = dot(centers[i], centers[j]);
      if (cosDist > thresh) {
        IntersectionEntry entryForI = getIntersectionEntry(centers, i, j);
        IntersectionEntry entryForJ = getIntersectionEntry(centers, j, i);
        //printf("    Centers %d and %d intersect at half-angle %f with directions %f and %f\n", i, j, entryForI.halfAngle, entryForI.direction, entryForJ.direction);
        
        auto iIntersections = result.find(i);
        if (iIntersections == result.end()) {
          std::vector<IntersectionEntry> entriesForI {entryForI};
          result.emplace(i, entriesForI);
        } else {
          iIntersections->second.push_back(entryForI);
        }
        
        auto jIntersections = result.find(j);
        if (jIntersections == result.end()) {
          std::vector<IntersectionEntry> entriesForJ {entryForJ};
          result.emplace(j, entriesForJ);
        } else {
          jIntersections->second.push_back(entryForJ);
        }
      }
    }
  }
  return result;
}
