#pragma once

struct Interval {
  double min;
  double max;
  
  bool operator<(const Interval& other) const {
    return min < other.min;
  }
  
  // Precondition: intersects(this, b) == true
  void merge(const Interval& b) {
    min = (min < b.min) ? min : b.min;
    max = (max > b.max) ? max : b.max;
  }
  
};

bool intersects(const Interval& a, const Interval& b) {
  return (a.min <= b.min && b.min <= a.max) ||
         (a.min <= b.max && b.max <= a.max) ||
         (b.min <= a.min && a.min <= b.max) ||
         (b.min <= a.max && a.max <= b.max);
}

class IntervalUnion {
private:
  std::vector<Interval> intervals;
  
public:
  void addInterval(double a, double b) {
    addInterval(Interval{a, b});
  }
  
  void addInterval(Interval x) {
    //printf("Adding interval (%f, %f)\n", x.min, x.max);
    // this will happen when the interval contains 0 or 2*pi
    if (x.min < 0) {
      addInterval(2*pi + x.min, 2*pi);
      addInterval(0, x.max);
      return;
    }
    if (x.max > 2*pi) {
      addInterval(x.min, 2*pi);
      addInterval(0, x.max - 2*pi);
      return;
    }
    std::vector<size_t> intervalsToRemove;
    for (auto it = intervals.begin(); it != intervals.end(); ++it) {
      //Interval i = intervals[idx];
      if (intersects(x, *it)) {
        x.merge(*it);
        --it;
        intervals.erase(it+1);
      }
    }
    
    //for (auto idxIt = intervalsToRemove.rbegin(); idxIt != intervalsToRemove.rend(); ++idxIt) {
    //  intervals.erase(idxIt);
    //}
    intervals.push_back(x);
    //printf("Added interval (%f, %f), nIntervals %d, new length %f\n", x.min, x.max, intervals.size(), getLength());
  }
  
  double getLength() {
    double len = 0;
    for (Interval i: intervals) {
      len += (i.max - i.min);
    }
    return len;
  }
};