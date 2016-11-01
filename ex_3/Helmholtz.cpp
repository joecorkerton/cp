//Calculate magnetic field from Helmholtz coils of wire with current I
//Current is 1/(mu0) to get simple solutions
#include <iostream>
#include <cmath>

using namespace std;

double const pi = 4 * atan(1);
int const N = 1E3; //# of elements to break up wire into
double const radius = 1.0; //Radius of loops

class vector3 {

	public:
		double x, y, z;
		vector3 (double a = 0, double b = 0, double c = 0) {
			
			x = a;
			y = b;
			z = c;
		}
		
		double length () {
			
			return sqrt(x * x + y * y + z * z);
		}
		
		vector3 operator % (const vector3 &v) { //compute cross product with overloaded % operator
			
			return vector3 (y * v.z - z * v.y,
					z * v.x - x * v.z,
					x * v.y - y * v.x);
		}
		
		vector3 operator - (const vector3 &v) { //compute difference between two vectors
		
			return vector3 (x - v.x,
					y - v.y,
					z - v.z);
		}
		
		vector3 operator + (const vector3 &v) {
		
			return vector3 (x + v.x,
					y + v.y,
					z + v.z);
		}
		
		vector3 operator * (const double a) {
		
			return vector3 (x * a,
					y * a,
					z * a);
		}
		
};

void initialise_loop (vector3 loop[N], vector3 centre) { 
	//assume wire perpendicular to x axis
	
	double step = 2 * pi / (N + 1);
	
	for (int i = 0; i <= N; i++) {
		
		loop[i].x = centre.x;
		loop[i].y = radius * cos (i * step) + centre.y;
		loop[i].z = radius * sin (i * step) + centre.z;
	}
}

void initialise_loop_reverse (vector3 loop[N], vector3 centre) { 
	//assume wire perpendicular to x axis
	
	double step = 2 * pi / (N + 1);
	
	for (int i = 0; i <= N; i++) {
		
		loop[i].x = centre.x;
		loop[i].y = radius * cos (i * step) + centre.y;
		loop[i].z = - radius * sin (i * step) + centre.z;
	}
}

double B_axial (double x) {
	//calculates analytic result for field on axis of a single coil

	return (radius * radius / (2 * pow((radius * radius + x * x), 1.5)));
}

vector3 bfield (vector3 loop[N], vector3 r) {
	//Applies biot savart law to calculate field from a line element. Sums over whole loop to give total field.
	
	vector3 B, dl, R;
	
	if ((abs(abs(r.y) - abs(radius)) < 1E-7) && (abs(abs(r.x) - abs(loop[0].x)) < 1E-7)) //avoid problems calculating field ontop of the loop
		return B;
	
	for (int i = 0; i <= N; i++) {
		
		if (i == N) 
			dl = loop[0] - loop[N];
			
		else
			dl = loop[i + 1] - loop[i];
		
		R = r - loop[i];
		B = B + ((dl % R) * (pow(R.length(), -3)/(4 * pi))); //dB = dl x R / (4 pi R ^ 3)
	}
	
	return B * 0.2;
}

void plotxy (vector3 loop1[N], vector3 loop2[N], double max_x, double max_y) {
	
	double step = 0.05; //step size for x,y co-ordinates
	vector3 r = vector3(-max_x, -max_y, 0);
	vector3 B;
	cout.precision(8);
	//double Btheory;
	
	while (r.x <= max_x) {
	
		while (r.y <= max_y) {
		
			B = bfield (loop1, r);
			B = B + bfield (loop2, r);
			cout << r.x << " " << r.y << " " << B.x << " " << B.y << " " << B.z << " " << B.length() << endl;
			r.y += step;
		}
		
		//used to calculate axial result
		/* B = bfield (loop, r);
		Btheory = B_axial(r.x);
		cout << r.x << " " << B.x << " " << Btheory << " " << Btheory - B.x << endl; */
		
		r.x += step;
		r.y = -max_y;
		cout << endl;
	}
}

int main() {

	double const max_x = 2;
	double const max_y = 2; //setting +- x & y values to plot data from

	//Initialise loops
	vector3 centre1 = vector3 (-0.5, 0, 0);
	vector3 centre2 = vector3 (0.5, 0, 0);
	vector3 loop1[N];
	vector3 loop2[N];
	
	initialise_loop(loop1, centre1);
	initialise_loop_reverse(loop2, centre2);
	plotxy(loop1, loop2, max_x, max_y);
	
	return 0;
}
