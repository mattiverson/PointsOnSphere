#pragma once

#include <unordered_map>

  // G7-K15 baybeeeeeeeeeeeeeeeeeeeeeee
  static constexpr double nodes[] = {0.004272314439594, 0.025446043828621, 0.067567788320115, 0.129234407200303, 0.206956382266154, 0.297077424311301, 0.396107522496051, 0.500000000000000, 0.603892477503949, 0.702922575688699, 0.793043617733846, 0.870765592799697, 0.932432211679885, 0.974553956171379, 0.995727685560406};
  static constexpr double kronWeights[] = {0.011467661005265, 0.031546046314989, 0.052395005161125, 0.070326629857763, 0.084502363319634, 0.095175289032393, 0.102216470037649, 0.104741070542364, 0.102216470037649, 0.095175289032393, 0.084502363319634, 0.070326629857763, 0.052395005161125, 0.031546046314989, 0.011467661005265};
  static constexpr double gaussWeights[] = {0, 0.064742483084435, 0, 0.139852695744638, 0, 0.190915025252559, 0, 0.208979591836735, 0, 0.190915025252559, 0, 0.139852695744638, 0, 0.064742483084435, 0};
  
  static constexpr size_t integralLength = 15;

class Integrator {
private:
  // store function evaluations to reuse them, since computing them
  // is very expensive.
  // this does add a ton of overhead, so this might actually be more of a loss than a gain.
  std::unordered_map<double, double> pastEvaluations;

public:  
  double tolerance = 0.0001;
  double minInterval = 0.01;
  
  template<typename Integrand>
  double integrate(Integrand integrand, double min, double max) {
    //double pointsToEval[integralLength];
    //double evaluations[integralLength];
    double kronIntegral = 0;
    double gaussIntegral = 0;
    for (size_t i = 0; i < integralLength; ++i) {
      double x = min + (max - min) * nodes[i];
      double y = 0;
      x = min + (max - min) * nodes[i];
      if (pastEvaluations.find(x) != pastEvaluations.end()) {
        y = pastEvaluations.find(x)->second;
      } else {
        y = integrand(x);
        // store evaluations at points that we might need later
        if (i == 0 || i == 7 || i == 14) {
          pastEvaluations.emplace(x, y);
        }
      }
      kronIntegral += y * kronWeights[i];
      gaussIntegral += y * gaussWeights[i];
    }
    kronIntegral *= (max - min);
    gaussIntegral *= (max - min);
    //printf("  Integral on (%f, %f) is %f to %f\n", min, max, gaussIntegral, kronIntegral);
    if (abs(kronIntegral-gaussIntegral) < tolerance * (max - min) || (max - min) < minInterval) {
      //printf("  Accepting integral on (%f, %f) as %f\n", min, max, kronIntegral);
      return kronIntegral;
    } else {
      //printf("  We must go deeper: splitting (%f, %f)\n", min, max);
      double mid = (max + min) / 2;
      return integrate(integrand, min, mid) + integrate(integrand, mid, max);
    }
  }
  
};

//-------------------------------------------------------------------------
// Integrator Tests
// ------------------------------------------------------------------------

double constantFunction(double x) {
  return 1;
}

double linearFunction(double x) {
  return x;
}

double sineFunction(double x) {
  return sin(x);
}

double cosineFunction(double x) {
  return cos(x);
}

double reciprocalFunction(double x) {
  return 1/x;
}

void integralTests() {
  printf("Testing numerical integration:\n");
  Integrator a;
  // integral from 0 to 1 of [1 dx] = 1
  printf("integral from 0 to 1 of [1 dx] (should = 1): %.15f\n", a.integrate(constantFunction, 0, 1));
  Integrator b;
  // integral from 0 to 5 of [x dx] = 1/2 * 5^2 = 12.5
  printf("integral from 0 to 5 of [x dx] (should = 12.5): %.15f\n", b.integrate(linearFunction, 0, 5));
  Integrator c;
  // integral from 0 to pi/3 of [sin(x) dx] = -cos(pi/3) + cos(0) = 1 - sqrt(3)/2 = 0.133974596215561
  printf("integral from 0 to pi/3 of [sin(x) dx] (should = 0.5): %.15f\n", c.integrate(sineFunction, 0, pi/3));
  Integrator d;
  // integral from 0 to pi/3 of [cos(x) dx] = sin(pi/3) = sqrt(3)/2 = 0.866025403784439
  printf("integral from 0 to pi/3 of [cos(x) dx] (should = sqrt(3)/2 = 0.866025403784439): %.15f\n", d.integrate(cosineFunction, 0, pi/3));
  Integrator e;
  // integral from 0.01 to 1 of [1/x dx] = log(1) - log(0.01) = log(100) = 4.60517018598809
  printf("integral from 0.01 to 1 of [1/x dx] (should = -log(0.01) = 4.60517018598809): %.15f\n", e.integrate(reciprocalFunction, 0.01, 1));
  printf("End of numerical integration tests\n\n");
  
}