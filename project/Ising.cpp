//Use Metropolis algorithm to investigate ferromagnetism using the Ising model.
//Periodic boundary conditions used so there are no edge effects

#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>
using namespace std;

//set up constants
int const N = 200; //Size of the lattice (N x N)
double const mu = 1; //Magnetic moment
double const H = 0; //External magnetic field
double const J = 1; //Exchange energy, set to 1 for simplicity

signed char metropolis (signed char n1, signed char n2, signed char n3, signed char n4, signed char s, gsl_rng *rng, double T) { 
	
	//signed chars are 4 nearest neighbours + current value of spin
	double dE; //change in energy if spin is reversed
	
	//Find change in energy from Ising model, E = - J sum <ij> {Si Sj} - mu * H sum i {Si}
	dE = J * (n1 + n2 + n3 + n4) * (2 * s) + mu * H * (2 * s);
	
	//flip spin if dE < 0 or exp( - dE / (k_B * T)) > p, for p randomly sampled in [0,1]
	if (dE < 0) 
		return (-1 * s);
	
	else {
		double p = exp(-dE / T);
		
		if (p > gsl_rng_uniform(rng))
			return (-1 * s);
		
		else 
			return s;
	}
}

void nearest_neighbour (int n[4], int i) { //find position of the nearest neighbours for a given index in the lattice

	if (i == 0)  { //if on the top left corner of the lattice	
		n[0] = i + 1; n[1] =  N - 1; n[2] = N * (N - 1); n[3] =  N;
	}		
	
	else if (i == (N * (N - 1))) {//if on bottom left corner
		n[0] = i + 1; n[1] = N * (N - 1); n[2] = i - N; n[3] = 0;
	}		
	
	else if ((i % N) == 0) { //if on the left edge & not in a corner
		n[0] = i + 1; n[1] = i + (N - 1); n[2] = i - N; n[3] = i + N;
	}		
	
	else if (i == (N - 1)) { //if on top right corner
		n[0] = N; n[1] = i - 1; n[2] = N * N - 1; n[3] = i + N;
	}		
	
	else if (i == (N * N - 1)) { //if on bottom right corner
		n[0] = N * (N - 1); n[1] = i - 1; n[2] = i - N; n[3] = N - 1;
	}		
	
	else if ((i % N) == (N - 1)) { //if on right edge & not in corner
		n[0] = i - (N - 1); n[1] = i - 1; n[2] = i - N; n[3] = i + N;
	}	
	
	else if (i < N) { //if on top edge & not in a corner
		n[0] = i + 1; n[1] = i - 1; n[2] = i + N * (N - 1); n[3] = i + N;
	}		
	
	else if (i > N * (N - 1)) { //if on bottom edge & not in a corner
		n[0] = i + 1; n[1] = i - 1; n[2] = i - N; n[3] = i % N;
	}		
	
	else {
		n[0] = i + 1; n[1] = i - 1; n[2] = i - N; n[3] = i + N;
	}
}

double energy (signed char spins[N * N], double &var) { //calculate sum E and sum E^2. E = - J sum <ij> {Si Sj} - mu * H sum i {Si}

	double E = 0, E2 = 0, tmp;
	int n[4] = {}; //array to pass values of nearest neighbours
	
	for (int i = 0; i < N * N; i++) {
		
		nearest_neighbour(n, i);
		tmp = - J * spins[i] * (spins[n[0]] + spins[n[1]] + spins[n[2]] + spins[n[3]]) - mu * H * spins[i];
		E += tmp;
		E2 += tmp * tmp;
	}
	
	tmp = E / (2 * N * N);
	var = (E2 / (4 * N * N) - tmp * tmp) / (N * N - 1); //1/2 factors to avoid double counting
	
	return tmp;
}

double magnetisation (signed char spins[N * N]) { //Magnetisation is given by the average spin value; M = 1/N sum i {Si}

	double M = 0;

	for (int i = 0; i < N * N; i++) {
		
		M += (double) spins[i];
	}
	
	return (M / (N * N));
}

void iterate (signed char spins[N * N], gsl_rng *rng, double T) {

	int n[4] = {}; //array to pass values of nearest neighbours

	//Iterate over the whole lattice, applying Metropolis algorithm to determine whether to flip spin or not
	for (int i = 0; i < N * N; i++) {
	
		nearest_neighbour(n, i);
		spins[i] = metropolis(spins[n[0]], spins[n[1]], spins[n[2]], spins[n[3]], spins[i], rng, T);
	}
}

void output(signed char spins[N * N]) { //outputs the current state of the lattice in a format to produce a 2d heat map
	
	for (int i = 0; i < N * N; i++) {
		
		cout << (int) (i / N) << " " << i % N << " " << (int)spins[i] << endl;
		
		if (i % N == (N - 1))
			cout << endl;
	}
}

void equilibrium(signed char spins[N * N], gsl_rng *rng, double j, double T) { //iterate 2 * j times to reach equilibrium, then a further j times to calculate properties of the system

	int n = 0; //counter for the number of iterations
	double m = 0, e = 0; //running total of magnetisation/energy
	double var = 0; //running total of the variance used to calculate heat capacity
	double tmp;

	while (n < (int)(2 * j)) {
	
		iterate(spins, rng, T);
		n++;
	}
	
	n = 0;
	/*
	while (n < (int)(j)) {
		
		m += magnetisation(spins);
		e += energy(spins, tmp);
		var += tmp;
		iterate(spins, rng, T);
		n++;
	}
	
	//output average magnetisation/energy/heat capacity per site
	//average heat capacity per site, C = sigma^2 / (T^2) in units k_B ^ 2 * K / J
	cout << N << " " << T << " " << abs(m / j) << " " << (e / j) << " " << (var / (j * T * T)) << endl;
	*/
	
	output(spins);
}

int main () {
	
	//Spins represented by signed char, 1 = up, -1 = down
	signed char spins[N * N];
	
	double Tmax = 4.1E0; // maximum temperature, units J / k_B
	double Tmin = 3E0; //minimum temperature to start loop from
	double dT = 4.0; //step size
	double T = Tmin; //Current temperature
	

	//initialise and seed rng
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, time(NULL));
	
	double j = 1.010E3; //number of iterations
	
	//calculate properties of the system after equilibrium is reached for a range of temperatures
	while (T < Tmax) {
	
		//Initialize lattice spins to all point up
		for (int i = 0; i < N * N; i++)
			spins[i] = 1;
	
		equilibrium(spins, rng, j, T);
		T += dT;
	}
	
	gsl_rng_free(rng);
	
	//output running time
	//cout << (clock()/ CLOCKS_PER_SEC) << endl;
	
	return 0;
}
