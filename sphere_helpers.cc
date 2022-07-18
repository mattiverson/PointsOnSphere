#pragma once

#include <random>
#include <math.h>
#include <vector>

double pi = 3.14159265358979323846;

struct double3;
double3 operator*(const double a, const double3& d);
double dot(const double3& a, const double3& b);
double3 normalize(const double3& x);

// a 3-dimensional vector.
struct double3 {
  double x, y, z;
  
  // vector-vector addition.
  double3 operator+(const double3& other) const {
    return double3{x + other.x, y + other.y, z + other.z};
  }
  
  // vector-vector subtraction.
  double3 operator-(const double3& other) const {
    return double3{x - other.x, y - other.y, z - other.z};
  }
  
  // vector-scalar division.
  double3 operator/(const double a) const {
    return double3{x / a, y / a, z / a};
  }
  
  // precondition: this and other are unit-length.
  // returns this vector minus its projection onto other.
  double3 minusProjection(const double3& other) {
    double proj = dot(*this, other);
    return *this - (proj*other);
  }
  
  // precondition: this and other are unit-length.
  // returns this vector minus its projection onto other, normalized.
  double3 minusProjectionUnit(const double3& other) {
    return normalize(minusProjection(other));
  }
};

// scalar-vector multiplication.
double3 operator*(const double a, const double3& d) {
  return double3{a * d.x, a * d.y, a * d.z};
}

// dot product.
double dot(const double3& a, const double3& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

// cross product.
double3 cross(const double3& a, const double3& b) {
  return double3{a.y * b.z - a.z * b.y,
                 a.z * b.x - a.x * b.z,
                 a.x * b.y - a.y * b.x};
}

// normalize: returns a unit vector parallel to this.
double3 normalize(const double3& x) {
  double len = sqrt(dot(x, x));
  return x / len;
}


std::random_device randomDevice;
std::default_random_engine engine(randomDevice());
std::uniform_real_distribution<> uniformSample(0, 1);

// Sample a random point uniformly from the unit sphere.
double3 sampleSphere() {
  double u = uniformSample(engine);
  double v = uniformSample(engine);
  
  // uniformly sampling (lat, lon) would bias toward the north pole; this biases the latitude sample toward the equator
  // in the correct way to cancel out that effect and produce a uniform sample over the sphere.
  v = acos(2*v-1);
  
  double x = cos(2*pi*u) * sin(v);
  double y = sin(2*pi*u) * sin(v);
  double z = cos(v);
  
  // mirror into the upper hemisphere, since that's all we're concerned about in this problem.
  z = (z < 0) ? -z : z;
  return double3{x, y, z};
}

// Sample n random points uniformly from the unit sphere.
std::vector<double3> sampleSphere(int n) {
  std::vector<double3> results;
  results.reserve(n);
  for (int i = 0; i < n; ++i) {
    results.push_back(sampleSphere());
  }
  return results;
}