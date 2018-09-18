//John Stanco 9.14.18


#include </usr/local/include/signal.h>
#include </usr/local/include/fftw3.h>

#include <iostream>
#include <vector>
#include <ctime>


void line( cx_double f[], size_t N ) {
	double k = ( 2 * pi ) / N;
	for( size_t i = 0; i < N; ++i ) {
		f[i] = cx_double{ cos( k * i ), sin( k * i ) };
	}
}


void line( std::vector<cx_double> &f ) {
	double k = ( 2 * pi ) / f.size();
	for( size_t i = 0; i < f.size(); ++i ) {
		f[i] = cx_double{ cos( k * i ), sin( k * i ) };
	}
}


void print( cx_double in[], size_t N ) {
	for( size_t i = 0; i < N; ++i ) {
		std::cout << "Re: " << in[i].real() << "\t\t";
		std::cout << "Im: " << in[i].imag() << "\n";
	}
}


void print( const std::vector<cx_double> &f ) {
	for( size_t i = 0; i < f.size(); ++i ) {
		std::cout << "Re: " << f[i].real() << "\t\t";
		std::cout << "Im: " << f[i].imag() << "\n";
	}
}


void fft_static( cx_double f[], size_t N ) {
	fftw_complex in[N];

	size_t n_bytes = N * sizeof( fftw_complex );

	fftw_plan p = fftw_plan_dft_1d( N, in, in, FFTW_FORWARD, FFTW_ESTIMATE );

	memcpy( in, f, n_bytes );
	fftw_execute( p );
	memcpy( f, in, n_bytes );

	fftw_destroy_plan( p );
}


void fft_dynamic( cx_double f[], size_t N ) {
	fftw_complex *in;

	size_t n_bytes = N * sizeof( fftw_complex );
	in = ( fftw_complex* )malloc( n_bytes );

	fftw_plan p = fftw_plan_dft_1d( N, in, in, FFTW_FORWARD, FFTW_ESTIMATE );

	memcpy( in, f, n_bytes );
	fftw_execute( p );
	memcpy( f, in, n_bytes );

	free( in );
	fftw_destroy_plan( p );
}


void fft( std::vector<cx_double> &f ) {
	fftw_complex *in;
	size_t N = f.size();
	size_t n_bytes = N * sizeof( fftw_complex );
	in = ( fftw_complex* )malloc( n_bytes );

	fftw_plan p = fftw_plan_dft_1d( N, in, in, FFTW_FORWARD, FFTW_ESTIMATE );

	memcpy( in, f.data(), n_bytes );
	fftw_execute( p );
	memcpy( f.data(), in, n_bytes );

	free( in );
	fftw_destroy_plan( p );
}


double avg_time_fft_static() {
	size_t N;
	size_t iter = 1000;
	clock_t T = 0;

	for( size_t i = 0; i < iter; ++i ) {
		N = randInt( 1 << 4, 1 << 6 );
		cx_double *f = ( cx_double* )malloc( N * sizeof( cx_double ) );
		line( f, N );
		//print( f, N );
		clock_t t = clock();
		fft_static( f, N );
		t = clock() - t;
		free( f );
		T += t;
	}

	//print( f, N );

	T /= iter;

	return (double)T / CLOCKS_PER_SEC;
}


double avg_time_fft_dynamic() {
	size_t N;
	size_t iter = 1000;
	clock_t T = 0;

	for( size_t i = 0; i < iter; ++i ) {
		N = randInt( 1 << 4, 1 << 6 );
		cx_double *f = ( cx_double* )malloc( N * sizeof( cx_double ) );
		line( f, N );
		//print( f, N );
		clock_t t = clock();
		fft_dynamic( f, N );
		t = clock() - t;
		free( f );
		T += t;
	}

	//print( f, N );

	T /= iter;

	return (double)T / CLOCKS_PER_SEC;
}



double avg_time_fft_vector() {
	size_t N;
	size_t iter = 1000;
	clock_t T = 0;

	for( size_t i = 0; i < iter; ++i ) {
		//N = randInt( 1 << 4, 1 << 6 );
		N = 64;
		std::vector<cx_double> f( N );
		line( f );
		//print( f );
		clock_t t = clock();
		fft( f );
		t = clock() - t;
		T += t;
	}

	//print( f );

	T /= iter;
	return (double)T / CLOCKS_PER_SEC;
}


void test_fft() {

	size_t N = 1 << 6;
	Signal::TSignal<cx_double> f{ N };

	

}




int main() {

	test_fft();
	//std::cout << avg_time_fft_dynamic() << std::endl;

	return 0;
}
