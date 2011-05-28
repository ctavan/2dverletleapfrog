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

struct TAverager
{
	double sum;
	int number;
	
	void init()
	{
		sum = 0.0;
		number = 0;
	}

	void add(double const& value)
	{
		sum += value;
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
double e_kin(TVector v);
double e_pot(TVector v);

// System parameters
TVector r;			// Particle position r(t)
TVector r_prev;		// r(t-tau) (needed for verlet only)
TVector r_next;		// r(t+tau)
TVector v;			// Particle velocity
TVector v_next;

double ekin = 0;			// Total kinetic energy
double epot = 0;			// Total potential energy
double etot = 0;			// Total energy
TAverager avg_epot;			// Averager for the total energy
TAverager avg_ekin;			// Averager for the total energy
TAverager avg_etot;			// Averager for the total energy

double dt = 0;				// Length of one timestep
double t = 0;				// Current time
int nt, nequil, nproduct;	// Number of timesteps
double tmax = 0;			// Total simulation time

FILE* outTrajectory;		// Outputfile for trajectories
FILE* outEnergies;			// Outputfile for energies

int main(int argc, char *argv[])
{
	if (argc < 4)
	{
		printf("Usage: %s <dt> <nt> <algo>\n", argv[0]);
		printf("\t<dt>       Length of one timestep\n");
		printf("\t<nt>       Number of timesteps to be performed\n");
		printf("\t<algo>     1 for verlet / 2 for leap frog\n");
		exit(EXIT_SUCCESS);
	}

	dt = atof(argv[1]);
	nt = atof(argv[2]);
	tmax = nt*dt;

	int algo = atoi(argv[3]);

	char filename[50];
	sprintf(filename, "outTrajectory_%s_%g.txt", (algo == 1) ? "verlet" : "leapfrog", dt);
	outTrajectory = fopen(filename, "w+");
	sprintf(filename, "outEnergies_%s_%g.txt", (algo == 1) ? "verlet" : "leapfrog", dt);
	outEnergies = fopen(filename, "w+");

	printf("Timestep length:\t%g\n", dt);
	printf("Total steps:\t\t%d\n", nt);
	printf("Total time:\t\t%g\n", tmax);
	printf("Algorithm:\t\t%s\n", (algo == 1) ? "verlet" : "leapfrog");


	// Init averagers
	avg_epot.init();
	avg_ekin.init();
	avg_etot.init();

	fprintf(outEnergies, "#t\t\tE_kin\t\t<E_kin>\t\tE_pot\t\t<E_pot>\t\tE_tot\t\t<E_tot>\n");


	// Init Particle position & velocity
	r.x = 1.0;
	r.y = 0.0;

	double r2 = (r.x*r.x+r.y*r.y);

	r_next.x = sqrt(r2-(dt*dt/sqrt(r2)));
	r_next.y = dt/sqrt(sqrt(r2));

	v_next.x = 0.0;
	v_next.y = 1.0/sqrt(sqrt(r2));

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

		if (algo == 1)
		{
			r_prev = r;
			r = r_next;
			r_next = r*2.0-r_prev+f*dt*dt;
			v_next = (r_next-r_prev)/(2.0*dt);
		}
		else
		{
			r = r_next;
			v = v_next;
			v_next = v + f*dt;
			r_next = r + v*dt;
		}
		ekin = e_kin(v_next);
		avg_ekin.add(ekin);
		epot = e_pot(r_next);
		avg_epot.add(epot);
		etot = ekin+epot;
		avg_etot.add(etot);
		fprintf(outEnergies, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", t, ekin, avg_ekin.average(), epot, avg_epot.average(), etot, avg_etot.average());
		fprintf(outTrajectory, "%e\t%e\t%e\t%e\t%e\n", t, r_next.x, r_next.y, v_next.x, v_next.y);
	}

	printf("Averages:\n");
	printf("Kinetic energy:\t\t%g\n", avg_ekin.average());
	printf("Potential energy:\t%g\n", avg_epot.average());
	printf("Total energy:\t\t%g\n", avg_etot.average());
	printf("Virial theorem <E_pot> = - 2 <E_kin>.\n");
	printf("  Difference:\t\t%g\n", avg_epot.average()+2.0*avg_ekin.average());
	printf("  Relative Diff.:\t%.3f %%\n", fabs((avg_epot.average()+2.0*avg_ekin.average())/avg_epot.average()));

	fclose(outTrajectory); outTrajectory = NULL;
	fclose(outEnergies); outEnergies = NULL;

	exit(EXIT_SUCCESS);
}

TVector force()
{
	double r2 = (r.x*r.x+r.y*r.y);
	double r3 = r2*sqrt(r2);

	return (r/r3)*-1.0;
}

double e_kin(TVector v)
{
	return 0.5*(v.x*v.x+v.y*v.y);
}

double e_pot(TVector r)
{
	return -1.0/sqrt(r.x*r.x+r.y*r.y);
}