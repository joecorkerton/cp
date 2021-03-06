//Solving ODE of a driven pendulum
//d^2 theta/dt^2 = -g/l sin(theta) - q d theta/dt + F sin(omega_d t)
//we set l = g and omega_d = 2/3 rad s^-1

#include <iostream>
#include <cmath>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
using namespace std;

const double pi = 4 * atan(1);

class pendulum {
	
	public:
		//constructor
		pendulum(double a, double b, double c) {
		
			q = a;
			F = b;
			theta0 = c;
		}
		double get_q() {
		
			return q;
		}
		double get_F() {
	
			return F;
		}
		double get_theta0() {
		
			return theta0;
		}
	
	private:
		double q;
		double F;
		double theta0;
};

int coupled_eqn (double t, const double y[], double f[], void *params) {
	
	//working in transformed variables y[0] = theta, y[1] = d theta/dt
	pendulum p = *(pendulum *)params;
	f[0] = y[1];
	f[1] = - sin(y[0]) - p.get_q() * y[1] + p.get_F() * sin(2*t/3);
	
	return GSL_SUCCESS;	
}

double energy (double theta, double v, double mass) {
	
	//mass per unit length^2 to reduce rounding errors
	double E = mass * v * v / 2 + mass * (1 - cos(theta));
	
	return E;
}

double solve(double q, double F, double theta0, double T) {
	
	pendulum p = pendulum (q, F, theta0);
	const int n_eqns = 2;
	const double mass = 1/(9.81*9.81); //mass per unit length^2 with l=g
	double y[n_eqns] = {theta0, 0};
	double t = 0.0;
	double t_max = 2 * pi * T;
	double h = 1E-3; //step size
	double E, period, theta_tmp, t_tmp = 0, sum_T = 0; //Energy, period, temporary value to store previous value of theta and t, and sum over periods
	int i = 0;

	//Create a stepping function
	gsl_odeiv2_step *gsl_step = gsl_odeiv2_step_alloc (gsl_odeiv2_step_rkf45, n_eqns);
	//Adaptive step control, equally weight theta and d theta/dt
	gsl_odeiv2_control *gsl_control = gsl_odeiv2_control_standard_new (0, 1E-7, 1, 1); //relative error within 1E-7
	//Create evolution function
	gsl_odeiv2_evolve *gsl_evolve = gsl_odeiv2_evolve_alloc (n_eqns);
	//setup the system, 2nd argument is jacobian which isn't needed
	gsl_odeiv2_system gsl_sys = {coupled_eqn, NULL, n_eqns, &p};
	
	while (t < t_max) {
		
		theta_tmp = y[0];
		E = energy(y[0], y[1], mass);
		cout << t << " " << y[0] << " " << y[1] << " " << E << endl;
		int status = gsl_odeiv2_evolve_apply (gsl_evolve, gsl_control, gsl_step, &gsl_sys, &t, t_max, &h, y);
		
		if ((y[0] * theta_tmp) < 0) { //if theta changes sign
			
			sum_T += t - t_tmp;
			i += 1;
			t_tmp = t;
			
		}
		if (status != GSL_SUCCESS) break;
	}
	gsl_odeiv2_evolve_free(gsl_evolve);
	gsl_odeiv2_step_free(gsl_step);
	period = sum_T * 2 / i;
	
	return period;
}

int main() {
	
	cout.precision(10);
	double theta0 = 0.20001;
	
	/*while (theta0 < pi) {
		double T = solve(0.5, 0.5, theta0, 1000);
		cout << theta0 << " " << T << endl;
		theta0 += 0.02;
	}*/
	double T = solve(0.5, 1.2, theta0, 1000);
	
	return 0;
}
