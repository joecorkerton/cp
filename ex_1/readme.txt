Used Monte-carlo simulation with gsl random number generator to numerically integrate an 8 dimensional integral. N random points were evaluated, and for each N this was repeated 25 times. An average of these 25 repeats was taken to give the average value and error in the numerical integration for a particular N. Iterated over N from 10 to 1E7 to investigate how the error varied with N.

As expected, found standard deviation = N^-1/2. There is an exact analytic solution to the integral, taking the value 537.1873411... For N = 1E7, calculated the value 537.190 with an error 0.014.
