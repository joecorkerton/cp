// Monte-Carlo & Numerical Integration 

#include <iostream>
#include <math.h>
#include <gsl/gsl_rng.h>       
#include <time.h>       
using namespace std;

const int d = 8; //# of dimensions
const double s = atan(1)/2; //limits
const double vol  = pow (s, d); //integration volume

double integrand (double x) {
	
	return sin(x); 	
}

double estimate (int N, gsl_rng *rng) {

	double sum = 0; //calculate f(r)
	double sum_x = 0; 
	
	for (int i = 0; i < N ; i++) {
		sum_x = 0;
		for (int j = 0; j < d; j++) {
			sum_x += (double)gsl_rng_uniform(rng) * s; //generate random position
		}
		sum += integrand(sum_x); //integrand is a function of the sum of positions	 
	}
	
	return sum * vol * 1E6 / N;	
}

void loop (int N, double &mean, double &s_dev, gsl_rng *rng) {
	
	int nt = 25; //# of times to estimate integral
	double sum = 0;
	double sum2 = 0;
	double tmp = 0;
	
	
	for (int i = 0; i < nt; i++) {
		tmp = estimate(N, rng);
		sum += tmp;
		sum2 += tmp * tmp;
	}
	
	mean = sum / nt;
	s_dev = sqrt(sum2 / nt - mean*mean);
}

int main() {
		
	double N = 1E8;
	double mean = 0;
	double s_dev = 0;
	double i = 10;
	cout.precision(14);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, time(NULL));
	
	while (i < N) {
		loop((int)i, mean, s_dev, rng);
		cout <<  (int)i << "," << mean << "," << s_dev << endl;
		i = i * 1.1;
	}
	gsl_rng_free(rng);
	
	return 0;
}

