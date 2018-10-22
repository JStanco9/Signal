// John Stanco 5.11.18

/**
 * 		Implements methods for sampling normal and uniform distributions
 **/


#include "signal.h"
#include <random>

double pi = std::acos(-1);
std::random_device rseed;
std::mt19937_64 generator(rseed());
std::normal_distribution<double> normal(0, 1);
std::uniform_real_distribution<double> UNI01(0, 1);
std::uniform_real_distribution<double> UNI11(-1, 1);

#define max_signed_int generator.max()
#define half_signed max_signed_int/2
#define double_signed max_signed_int*2

#define rand_int32 (signed)(generator() - half_signed)
#define UNI (double)generator() / double_signed / 2
#define UNI1 (double)(rand_int32) / double_signed
//#define UNI UNI01(generator)
//#define UNI1 UNI11(generator)


static double table[ 257 ][ 2 ];
static double r;
static double ks[ 256 ];
static double ws[ 256 ];
static double ys[ 256 ];


double boxMuller(){
	double s = 2;
	double u = 0;
	double v = 0;
	while(s > 1){
		u = UNI;
		v = UNI;
		s = u * u + v * v;
	}
	double y = u * sqrt(-2 * log(s) / s);
	return (UNI1 > 0) ? y : -y;
}


double invGaussian(double x){
	return sqrt(-2 * log(x));
}


void compute_zig_table(){
	double A = .00492867323399;
	table[255][1] = .00126028593;
	for(size_t i = 255; i > 0; i--){
		table[i][0] = invGaussian(table[i][1]);
		table[i - 1][1] = table[i][1] + A / table[i][0];
	}
	table[255][0] = A / table[255][1];
	r = table[255][0];
}


void compute_ks(){
	for(size_t i = 1; i < 256; i++){
		ks[i] = table[i - 1][0] / table[i][0] * double_signed;
	}
	ks[0] = r / table[255][0] * double_signed;
}


void compute_ws(){
	for(size_t i = 1; i < 256; i++){
		ws[i] = table[i][0] / double_signed;
	}
	ws[0] = table[255][0] / double_signed;
}


void compute_ys(){
	for(size_t i = 0; i < 256; i++){
		ys[i] = table[i][1];
	}
}


int zigSet(){
	compute_zig_table();
	compute_ys();
	compute_ws();
	compute_ks();
	return 1;
}


static int set = zigSet();


double ziggurat() {
	double x, y; int j; unsigned short i;
	while(true){
		j = rand_int32;
		i = j&255;
		if(abs(j) < ks[i]){
			return j * ws[i];
		}
		while(i == 0){
			x = -log(UNI) / x;
			y = -log(UNI);
			if(y + y > x * x){
				return (j > 0)? x + r : -x - r;
			}
		}
		x = j * ws[i];
		y = UNI * (ys[i] - ys[i - 1]);
		if(y < exp(-.5 * x * x) - ys[i]){
			return x;
		}
	}
}


double min(const double x1, const double x2){
	return(x1 > x2)? x2 : x1;
}


double randDouble(const double fMin, const double fMax){
    double f = UNI;
    return fMin + f * (fMax - fMin);
}


int randInt(const int a, const int b){
	int interval = b - a;
	if(interval < 0){
		throw "|  function: rand_int  |  file: myrand.cpp  |  error:  arguments must be specified in ascending order  |";
	} else if(interval == 0){
		return 0;
	}
	return a + generator() % (interval + 1);
}


int randInt(const int a){
	if(a == 0){ return 0; }
	return (a > 0)? randInt(0, a - 1) : randInt(a + 1, 0);
}


int flip(const double x){
	if(0 > x || x > 1){
		throw "|  function: flip  |  file: myrand.cpp  |  error:  input value 'x' must be between 0 and 1  |";
	}
	if(UNI < x){ return 1; }
	return 0;
}


double randNorm() {
	return ziggurat();
}


double randNorm( double mean, double var ) {
	return ziggurat() * var + mean;
}
