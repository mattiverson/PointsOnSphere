#pragma once

#include<vector>
#include<math.h>
#include "sphere_helpers.cc"

using Alg = double(const std::vector<double3>&, const double);
using Test = void(Alg);

//-------------------------------------------------------------------------
// Correctness Tests
// ------------------------------------------------------------------------

// single circle, centered at the north pole, with a 90 degree radius: covers the entire hemisphere.
void testSinglePointFullCoverage(Alg algorithm) {
  double3 pole{0, 0, 1};
  std::vector<double3> spherePoints = {pole};
  double phi = pi/2;
  
  // in a plane, area of a circle is the integral from 0 to R of 2*pi*r dr => pi*R^2
  // on a sphere, area of a circle is the integral from 0 to R of 2*pi*sin(r) dr => 2*pi*(1-cos(R))
  // in this case, we're interested in the area as a fraction of a hemisphere's area (which is 2*pi)
  // as extra verification, using the 2-term taylor series cos(x) = 1 - 1/2 * x^2 gives us back pi*R^2
  double expected = 1;
  double actual = algorithm(spherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testSinglePointFullCoverage\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

// single circle, centered on the equator, with a 90 degree radius: covers half the hemisphere.
void testSinglePointHalfCoverage(Alg algorithm) {
  double isq2 = 1/sqrt(2);
  double3 pole{isq2, isq2, 0};
  std::vector<double3> spherePoints = {pole};
  double phi = pi/2;
  
  // in a plane, area of a circle is the integral from 0 to R of 2*pi*r dr => pi*R^2
  // on a sphere, area of a circle is the integral from 0 to R of 2*pi*sin(r) dr => 2*pi*(1-cos(R))
  // in this case, we're interested in the area as a fraction of a hemisphere's area (which is 2*pi)
  // as extra verification, using the 2-term taylor series cos(x) = 1 - 1/2 * x^2 gives us back pi*R^2
  double expected = 0.5;
  double actual = algorithm(spherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testSinglePointHalfCoverage\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

// single circle, centered at the north pole, with a 30 degree radius.
void testSingleBigPointAtPole(Alg algorithm) {
  double3 pole{0, 0, 1};
  std::vector<double3> spherePoints = {pole};
  double phi = 30.0 * pi/180;
  
  // in a plane, area of a circle is the integral from 0 to R of 2*pi*r dr => pi*R^2
  // on a sphere, area of a circle is the integral from 0 to R of 2*pi*sin(r) dr => 2*pi*(1-cos(R))
  // in this case, we're interested in the area as a fraction of a hemisphere's area (which is 2*pi)
  // as extra verification, using the 2-term taylor series cos(x) = 1 - 1/2 * x^2 gives us back pi*R^2
  double expected = (1 - cos(phi));
  double actual = algorithm(spherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testSingleBigPointAtPole\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

// single circle, centered at the north pole, with a 3 degree radius.
void testSingleSmallPointAtPole(Alg algorithm) {
  double3 pole{0, 0, 1};
  std::vector<double3> spherePoints = {pole};
  double phi = 3.0 * pi/180;
  
  double expected = (1 - cos(phi));
  double actual = algorithm(spherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testSingleSmallPointAtPole\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

// single circle, centered at 45 degrees north, with a 30 degree radius.
void testSingleBigPointAtMidLatitude(Alg algorithm) {
  double isq2 = 1/sqrt(2);
  double3 pole{0, isq2, isq2};
  std::vector<double3> spherePoints = {pole};
  double phi = 30.0 * pi/180;
  
  double expected = (1 - cos(phi));
  double actual = algorithm(spherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testSingleBigPointAtMidLat\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

// single circle, centered at 45 degrees north, with a 3 degree radius.
void testSingleSmallPointAtMidLatitude(Alg algorithm) {
  double isq2 = 1/sqrt(2);
  double3 pole{isq2, 0, isq2};
  std::vector<double3> spherePoints = {pole};
  double phi = 3.0 * pi/180;
  
  double expected = (1 - cos(phi));
  double actual = algorithm(spherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testSingleSmallPointAtMidLat\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

// single circle, centered on the equator, with a 30 degree radius.
void testSingleBigPointAtEquator(Alg algorithm) {
  double3 equator{1, 0, 0};
  std::vector<double3> spherePoints = {equator};
  double phi = 30.0 * pi/180;
  
  // since we're only interested in the upper hemisphere,
  // the expected area is half the circle's area
  double expected = (1 - cos(phi)) / 2.0;
  double actual = algorithm(spherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testSingleBigPointAtEquator\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

// single circle, centered on the equator, with a 3 degree radius.
void testSingleSmallPointAtEquator(Alg algorithm) {
  double3 equator{0, 1, 0};
  std::vector<double3> spherePoints = {equator};
  double phi = 3.0 * pi/180;
  
  double expected = (1 - cos(phi)) / 2.0;
  double actual = algorithm(spherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testSingleSmallPointAtEquator\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

// two circles, touching but not overlapping at the north pole, 30 degree radius.
void testTwoDisjointPoints(Alg algorithm) {
  double sq32 = sqrt(3)/2;
  double3 p1{0.5, 0, sq32};
  double3 p2{-0.5, 0, sq32};
  std::vector<double3> spherePoints = {p1, p2};
  double phi = 30.0 * pi/180;
  
  double expected = 2*(1 - cos(phi));
  double actual = algorithm(spherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testTwoDisjointPoints\t\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

// two circles, partially overlapping around the north pole, 30 degree radius.
void testTwoIntersectingPoints(Alg algorithm) {
  double sq32 = sqrt(3)/2;
  double3 p1{0.5, 0, sq32};
  double3 p2{0, 0, 1};
  std::vector<double3> spherePoints = {p1, p2};
  double phi = 30.0 * pi/180;
  
  double expected = 0.21426941566577; // based on integralM results
  double actual = algorithm(spherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testTwoIntersectingPoints\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

// two circles in the same place, completely overlapping around the north pole, 30 degree radius.
void testTwoDuplicatePoints(Alg algorithm) {
  double sq32 = sqrt(3)/2;
  double3 p1{0.5, 0, sq32};
  double3 p2{0.5, 0, sq32};
  std::vector<double3> spherePoints = {p1, p2};
  double phi = 30.0 * pi/180;
  
  double expected = (1 - cos(phi));
  double actual = algorithm(spherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testTwoDuplicatePoints\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}
// 5 points, sampled at random from around the sphere, 30 degree radius.
std::vector<double3> fiveRandomSpherePoints = sampleSphere(5);
void testFiveRandomPoints(Alg algorithm) {
  double phi = 30.0 * pi/180;
  
  double expected = 0.40976470016892; // based on integralM results
  double actual = algorithm(fiveRandomSpherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testFiveRandomPoints\t\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

// 5 points chosen around the sphere, 4 equally spaced around the equator and one at the north pole, 30 degree radius.
void testFivePoints(Alg algorithm) {
  double phi = 30.0 * pi/180;
  double3 p1{1, 0, 0};
  double3 p2{-1, 0, 0};
  double3 p3{0, 1, 0};
  double3 p4{0, -1, 0};
  double3 p5{0, 0, 1};
  std::vector<double3> spherePoints = {p1, p2, p3, p4, p5};
  double expected = 3*(1 - cos(phi));
  double actual = algorithm(spherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testFivePoints\t\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

void runAllAccuracyTests(Alg algorithm) {
  testSinglePointFullCoverage(algorithm);
  testSinglePointHalfCoverage(algorithm);
  
  testSingleBigPointAtPole(algorithm);
  testSingleSmallPointAtPole(algorithm);
  testSingleBigPointAtMidLatitude(algorithm);
  testSingleSmallPointAtMidLatitude(algorithm);
  testSingleBigPointAtEquator(algorithm);
  testSingleSmallPointAtEquator(algorithm);
  
  testTwoDisjointPoints(algorithm);
  testTwoIntersectingPoints(algorithm);
  //testTwoDuplicatePoints(algorithm); removed the requirement that algorithms handle duplicates correctly
  
  testFiveRandomPoints(algorithm);
  testFivePoints(algorithm);
}

//-------------------------------------------------------------------------
// Stress Tests
// ------------------------------------------------------------------------

// stress test with 100 random points.
static std::vector<double3> hundredRandomSpherePoints = sampleSphere(100);
void testHundredRandomPoints(Alg algorithm) {
  double phi = 10.0 * pi/180;
  
  double expected = 0.78453789923011; // based on IntegralM
  double actual = algorithm(hundredRandomSpherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testHundredRandomPoints\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

// stress test with 1000 random points.
static std::vector<double3> thousandRandomSpherePoints = sampleSphere(1000);
void testThousandRandomPoints(Alg algorithm) {
  double phi = 3.0 * pi/180;
  
  double expected = 0.74175418923266; // based on IntegralM
  double actual = algorithm(thousandRandomSpherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testOneThousandRandomPoints\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

// stress test with 10,000 random points.
static std::vector<double3> tenKRandomSpherePoints = sampleSphere(10000);
void testTenKRandomPoints(Alg algorithm) {
  double phi = 1.0 * pi/180;
  
  double expected = 0.7807;
  double actual = algorithm(tenKRandomSpherePoints, phi);
  double delta = abs(expected - actual);
  double relDelta = (delta / expected) * 100;
  bool pass = delta < 0.01;
  printf("  testTenKRandomSpherePoints\t%s: expected %.3e, actual %.3e, delta %.3e (%.1f\%)\n", (pass ? "PASSED" : "FAILED"), expected, actual, delta, relDelta);
}

void runStressTest(Alg algorithm) {
  testHundredRandomPoints(algorithm);
  testThousandRandomPoints(algorithm);
  //testTenKRandomPoints(algorithm);
}