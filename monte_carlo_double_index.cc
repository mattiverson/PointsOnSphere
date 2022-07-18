#pragma once

#include <unordered_map>

#include "sphere_helpers.cc"

// Monte Carlo approximation, with a small optimization.
// We draw a grid over the sphere at the start, and group
// the circle centers together into their cells. When a monte
// carlo point searches for circle centers it might be near, it
// starts by searching in its cell, then checks the surrounding cells,
// so it doesn't have to check every center.

// A bin/cell, which we use to organize nearby points.
struct Bin {
  int x, y;
  size_t binLength;
  bool inBounds;
  
  void checkInBounds() {
    bool outOfBounds = x < 0;
    outOfBounds |= y < 0;
    outOfBounds |= x >= binLength;
    outOfBounds |= y >= binLength;
    
    inBounds = !outOfBounds;
  }
  
  Bin(int x_, int y_, size_t binLength_):
    x{x_}
  , y{y_}
  , binLength{binLength_}  {
    checkInBounds();
  }
  
  Bin operator+(const Bin& other) const {
    return Bin(x + other.x, y + other.y, binLength);
  }
  
  bool operator==(const Bin& other) const {
    return x == other.x && y == other.y;
  }
};

// a hashing function for a bin, since we use bins as keys in a hash map.
struct BinHash {
  size_t operator() (const Bin& bin) const {
    return std::hash<int>()(bin.x) ^ std::hash<int>()(bin.y);
  }
};

using BinPoints = std::unordered_map<Bin, std::vector<size_t>, BinHash>;

// the order to traverse nearby cells.
const std::vector<Bin> binDeltas {Bin(0,0,0), Bin(1,0,0), Bin(0,1,0), Bin(-1,0,0), Bin(0,-1,0), Bin(1,-1,0), Bin(1,1,0), Bin(-1,1,0), Bin(-1,-1,0)};

// protect against roundoff error: if phi is pi/n for some
// integer n, a point at (-1,0,0) would end up in an off-by-one bin when
// cos(n * (pi/n)) = -0.99999999999. Technically other parts of the
// code would also protect against this, but I don't want to create
// an awful-to-find floating-point bug if I change that later.
inline double adjustedPhi(const double phi) {
  return phi + 0.001;
}

size_t indexSize(const double phi) {
  // to make sure that one circle never touches more than 3 slices
  // of the index, the width of an index slice must be at least the
  // radius of the circle.
  return ceil(pi / adjustedPhi(phi));
}

// calculate where the boundaries of each bin should be. This is
// based on the radius of the circle so we can guarantee that points
// will only touch circles in their own bins or the 8 surrounding bins.
std::vector<double> makeBinThresholds(const double phi) {
  const size_t nThresholds = indexSize(phi);
  std::vector<double> thresholds;
  thresholds.reserve(nThresholds);
  for (size_t i = 0; i < nThresholds; ++i) {
    //printf("bin threshold: %.2f\n", -cos((pi * i)/nThresholds));
    thresholds.push_back(-cos(i*adjustedPhi(phi)));
  }
  return thresholds;
}

// finds which bin a coordinate falls into.
size_t getBin(const std::vector<double>& binThresholds, const double x) {
  size_t minIdx = 0;
  size_t maxIdx = binThresholds.size();
  while (maxIdx - minIdx > 1) {
    size_t checkIdx = (maxIdx + minIdx) / 2;
    if (x < binThresholds[checkIdx]) {
      maxIdx = checkIdx;
    } else {
      minIdx = checkIdx;
    }
  }
  return minIdx;
}

// populates the given BinPoints data structure. Each non-empty bin will map to a vector of
// size_t's, corresponding to which elements of points are in that bin.
void fillBinPoints(const std::vector<double3>& points, BinPoints& binPoints, const double phi) {
  const std::vector<double> binThresholds = makeBinThresholds(phi);
  size_t nBins = binThresholds.size();
  for (size_t i = 0; i < points.size(); ++i) {
    Bin bin(getBin(binThresholds, points[i].x), getBin(binThresholds, points[i].y), nBins);
    auto entry = binPoints.find(bin);
    if (entry == binPoints.end()) {
      std::vector<size_t> indicesInBin {i};
      binPoints.emplace(bin, indicesInBin);
      //printf("Added new bin (%d, %d) and added point no. %d\n", bin.x, bin.y, i);
    } else {
      entry->second.push_back(i);
    }
  }
}

// estimate the proportion of the upper hemisphere which is within an angle of phi
// of at least one of the unit vectors in centers
double monteCarloDoubleIndex(const std::vector<double3>& centers, const double phi) {
  int nSamples = 10000;
  int sampleHits = 0;
  double thresh = cos(phi);
  std::vector<double3> points = sampleSphere(nSamples);
  
  const std::vector<double> binThresholds = makeBinThresholds(phi);
  int nBins = binThresholds.size();
  //printf("There are %d bins.\n", nBins);
  BinPoints centersByBin;
  fillBinPoints(centers, centersByBin, phi);
  for (double3 p : points) {
    Bin pointBin(getBin(binThresholds, p.x), getBin(binThresholds, p.y), nBins);
    // search only the bins close to the point
    for (Bin delta: binDeltas) {
      auto entry = centersByBin.find(pointBin + delta);
      
      // if the bin isn't in the map, no centers are in that bin
      if (entry == centersByBin.end()) {
        //printf("bin (%d, %d) empty\n", pointBin.x + delta.x, pointBin.y + delta.y);
        continue;
      }
      //printf("bin (%d, %d) contains %d points\n", entry->first.x, entry->first.y, entry->second.size());
      for (size_t centerIdx: entry->second) {
        double3 s = centers[centerIdx];
        if(dot(p, s) > thresh) {
          //printf("hit\n");
          ++sampleHits;
          goto endOfPointLoop;
        }
      }
    }
    endOfPointLoop: continue;
  }
    //for (double3 s : centers) {
    //  if(dot(p, s) > thresh) {
    //    ++sampleHits;
    //    break;
    //  }
    //}
  return sampleHits / (double) nSamples;
}

