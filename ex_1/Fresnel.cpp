// Fresnel integrals using GSL integration routine

#include <iostream>
#include <cmath>
#include <gsl/gsl_integration.h>
using namespace std;

const double pi = 4*atan(1);
gsl_function F;

double c (double x, void *params) {
	
	double alpha = *(double *) params; //recovers alpha
	return cos(pi*alpha*x*x/2);
}

double s (double x, void *params) {
	
	double alpha = *(double *) params;
	return sin(pi*alpha*x*x/2);
}

void integration (double &result, double &error, gsl_integration_workspace *w, double u) {
	
	double alpha = 1.0; //parameter in integrand
	
	F.params = &alpha;
	
	gsl_integration_qag (&F, 0, u, 0, 1e-7, 1000, 4, w, &result, &error);
	//compute integration within a relative error of 1e7
}

void cornu () {
	
	gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
	double u = 10; //upper limit on integral
	double result_c, result_s, error;
	double i = 0.01;
	
	while (i < u) {
		F.function = &c;
		integration(result_c, error, w, i);
		F.function = &s;
		integration(result_s, error, w, i);
		cout << result_s << "	" << result_c << endl;
		i += 0.01;
	}
	
	gsl_integration_workspace_free (w);
}

int main() {
	
	cout.precision(18);
	cornu();
	
	return 0;
}
