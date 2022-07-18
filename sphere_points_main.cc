#include <iostream>
#include <chrono>

#include "sphere_helpers.cc"
#include "monte_carlo_slow.cc"
#include "monte_carlo_double_index.cc"
#include "test_cases.cc"
#include "integrator.cc"
#include "integral_a.cc"
#include "integral_m.cc"

using std::chrono::steady_clock;

// run the test suites for the MonteCarloSlow algorithm.
void testMonteCarloSlow() {
  printf("Testing monteCarloSlow:\n");
  runAllAccuracyTests(monteCarloSlow);
  
  printf("\nStarting stress tests:\n");
  steady_clock::time_point start = steady_clock::now();
  runStressTest(monteCarloSlow);
  steady_clock::time_point stop = steady_clock::now();
  int durationMicros = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  printf("monteCarloSlow ran in %06d us\n\n\n\n\n", durationMicros);
}

// run the test suites for the MonteCarloDoubleIndex algorithm.
void testMonteCarloDoubleIndex() {
  printf("Testing monteCarloDoubleIndex:\n");
  runAllAccuracyTests(monteCarloDoubleIndex);
  
  printf("\nStarting stress tests:\n");
  steady_clock::time_point start = steady_clock::now();
  runStressTest(monteCarloDoubleIndex);
  steady_clock::time_point stop = steady_clock::now();
  int durationMicros = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  printf("monteCarloDoubleIndex ran in %06d us\n\n\n\n\n", durationMicros);
}

// run the test suites for the IntegralA algorithm.
void testIntegralA() {
  printf("Testing IntegralA:\n");
  runAllAccuracyTests(integralA);
  
  printf("\nStarting stress tests:\n");
  steady_clock::time_point start = steady_clock::now();
  runStressTest(integralA);
  steady_clock::time_point stop = steady_clock::now();
  int durationMicros = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  printf("integralA ran in %06d us\n\n\n\n\n", durationMicros);
}

// run the test suites for the IntegralM algorithm.
void testIntegralM() {
  printf("Testing IntegralM:\n");
  runAllAccuracyTests(integralM);
  
  printf("\nStarting stress tests:\n");
  steady_clock::time_point start = steady_clock::now();
  runStressTest(integralM);
  steady_clock::time_point stop = steady_clock::now();
  int durationMicros = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  printf("integralM ran in %06d us\n\n\n\n\n", durationMicros);
}

// this was used for debugging, not currently in use
void testIntegrandA() {
  double phi = 30.0 * pi/180;
  
  double3 p1{1, 0, 0};
  double3 p2{-1, 0, 0};
  double3 p3{0, 1, 0};
  double3 p4{0, -1, 0};
  double3 p5{0, 0, 1};
  std::vector<double3> badPoints = {p1, p2, p3, p4, p5};
  
  std::vector<double3> goodPoints = {p1};
  printf("bad integrand: %f\n", integrandA(pi/2, badPoints, phi));
  printf("good integrand: %f\n", integrandA(pi/2, goodPoints, phi));
}

// runs the tests.
int main() {
  integralTests();
  
  testMonteCarloSlow();
  testMonteCarloDoubleIndex();
  testIntegralA();
  testIntegralM();
  
  return 0;
}