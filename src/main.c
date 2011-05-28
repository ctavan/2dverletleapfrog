// Verlet/LeapFrog
// Patric Zimmermann, Dortje Schirok, Christoph Tavan
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>

#define PI 3.141592654

/**
 * Structure definitions
 */
struct TVector {
	double x;
	double y;

	TVector& operator=(TVector const& vector2) {
		if(this != &vector2) {
			x = vector2.x;
			y = vector2.y;
		}
		return *this;
	}

	TVector& operator=(double const& value) {
		x = value;
		y = value;
		return *this;
	}

	TVector& operator+=(TVector const& vector) {
		x += vector.x;
		y += vector.y;
		return *this;
	}

	TVector& operator-=(TVector const& vector) {
		x -= vector.x;
		y -= vector.y;
		return *this;
	}

	TVector& operator/=(double const& value) {
		x /= value;
		y /= value;
		return *this;
	}

	TVector& operator*=(double const& value) {
		x *= value;
		y *= value;
		return *this;
	}

	TVector operator*(double const& value) const {
		TVector temp(*this);
		temp.x *= value;
		temp.y *= value;
		return temp;
	}

	TVector operator/(double const& value) const {
		TVector temp(*this);
		temp.x /= value;
		temp.y /= value;
		return temp;
	}

	TVector operator+(TVector const& vector) const {
		TVector temp(*this);
		temp.x += vector.x;
		temp.y += vector.y;
		return temp;
	}

	TVector operator-(TVector const& vector) const {
		TVector temp(*this);
		temp.x -= vector.x;
		temp.y -= vector.y;
		return temp;
	}
};

struct TPartDist {				// for particle distances
	double dx;
	double dy;
	double r2;
};

struct TAverager
{
	double sum;
	double ssum;
	int number;
	
	void init()
	{
		sum = 0.0;
		ssum = 0.0;
		number = 0;
	}

	void add(double const& value)
	{
		sum += value;
		// ssum += value*value;
		number++;
	}

	double average()
	{
		return sum/(double)number;
	}
};

/**
 * Global variables
 */
const bool debug = true;			// Enable/disable debug output

TVector force();

// System parameters
TVector r;			// Particle position r(t)
TVector r_prev;		// r(t-tau) (needed for verlet only)
TVector r_next;		// r(t+tau)
TVector v;			// Particle velocity

// TVector vsum;				// Velocity of center of mass
// double vsum2 = 0;			// Mean squared velocity = 2 E_kin
// double epot = 0;			// Total potential energy
// double virial = 0;			// Total virial

double dt = 0;				// Length of one timestep
double t = 0;				// Current time
int nt, nequil, nproduct;	// Number of timesteps
double tmax = 0;			// Total simulation time

FILE* outTrajectory;		// Outputfile for trajectories

int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		printf("Usage: %s <dt> <nt>\n", argv[0]);
		printf("\t<dt>       Length of one timestep\n");
		printf("\t<nt>       Number of timesteps to be performed\n");
		exit(EXIT_SUCCESS);
	}

	outTrajectory = fopen("outTrajectory.txt", "w+");

	dt = atof(argv[1]);
	nt = atof(argv[2]);
	tmax = nt*dt;

	printf("Timestep length:\t%g\n", dt);
	printf("Total steps:\t\t%d\n", nt);
	printf("Total time:\t\t%g\n", tmax);

	r.x = 1.0;
	r.y = 0.0;

	double r2 = (r.x*r.x+r.y*r.y);

	r_next.x = sqrt(r2-(dt*dt/sqrt(r2)));
	r_next.y = dt/sqrt(sqrt(r2));

	printf("Startposition r(t-tau)\t= (%e,%e)\n", r.x, r.y);
	printf("Startposition r(t)\t= (%e,%e)\n", r_next.x, r_next.y);


	fprintf(outTrajectory, "#t\tn\tr_x\t\tr_y\t\tv_x\t\tv_y\n");
	fprintf(outTrajectory, "%e\t%e\t%e\n", t, r.x, r.y);
	fprintf(outTrajectory, "%e\t%e\t%e\n", t, r_next.x, r_next.y);
	for (int n = 0; n <= nt; n++)
	{
		// Current time
		t = dt*n;
		// if(debug) printf("t:\t%6.3f\t\n", t);

		TVector f = force();

		// if(debug) printf("Steplength:\t%e\t\n", sqrt((r_prev.x-r.x)*(r_prev.x-r.x)+(r_prev.y-r.y)*(r_prev.y-r.y)));

		r_prev = r;
		r = r_next;
		r_next = r*2.0-r_prev+f*dt*dt;
		fprintf(outTrajectory, "%e\t%e\t%e\n", t, r_next.x, r_next.y);
	}

	exit(EXIT_SUCCESS);
}

TVector force() {
	double r2 = (r.x*r.x+r.y*r.y);
	double r3 = r2*sqrt(r2);

	return (r/r3)*-1.0;
}