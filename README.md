# Sphere Coverage

This repo solves the problem:

## Given a radius and a set of points on the unit sphere, what proportion of the upper hemisphere would be covered by circles of the given radius centered at each of the given points?

There are 3 different approaches used:

- Monte Carlo
- IntegralA
- IntegralM

Monte Carlo is self-explanatory. There is a completely naive implementation, as well as one that uses spacial indexing for a performance boost.

IntegralA is an approach based on numerical integration. 
In this approach, we integrate over coelevation (pi/2 - elevation, or the
separation angle from the north pole) from 0 to pi/2. The integrand is a
function representing how much of the sphere at that coelevation is contained
in one of the circles. We evaluate this function for a given coelevation by
slicing the sphere at that coelevation and searching through all the circles 
that intersect with that slice, recording which azimuth intervals of the slice
are inside of circles.

IntegralM is another approach based on numerical integration.
In this approach, we integrate over distance from 0 to the given radius. The
integrand is a function representing the measure of the set of points that are
the given distance from the closest center.

#Performance

The algorithms all have tunable parameters that can control speed v.s. accuracy.
For Monte Carlo, the number of samples can be configured; for the numerical integration-based
solutions, the minimum interval size and the error tolerance per interval length can be changed.

With the current parameter choices, the Monte Carlo solutions complete the stress tests in about 100ms
(or 50ms with indexing) and have relative errors on the order of 1% for most test cases. In the test
cases where a very small proportion of the hemisphere is covered, its relative errors swell to nearly 50%.

The numerical integration solutions run slower but have much higher accuracy: IntegralA takes 1200ms to complete
the stress tests, but its relative error is on the order of 0.01% for most tests. The one exception is the small-circle-at-pole test case, where its error is 0.4%, but this is a worst-case scenario for the integrator. Also, it should be noted that IntegralA currently considers every circle in every slice, even though it could easily filter down to only the centers that intersect a given slice, so a large performance improvement is still on the table.

IntegralM takes 522ms to complete the stress tests, and its relative error is at least on par with IntegralA in every test, on the order of 0.01%. For most of the simpler tests, IntegralM achieves accuracy near the double-precision limit, roughly 0.000000000000001%.

The full results can be seen in output.txt.