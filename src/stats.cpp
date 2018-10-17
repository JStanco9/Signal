//John Stanco 9.12.18

#include "signal.h"

#include <random>

/**
 * 		Implements statistical analysis interface for TSignal<double>
 **/

double Signal::mean( const Signal::TSignal<double>& x ) {
	double mu = 0;
	for( size_t i = 0; i < x.length(); i++ ){
		mu += x( i );
	}
	return mu / x.length();
}


double Signal::var( const Signal::TSignal<double>& x ) {
	if( x.length() < 2 ) { return 0; }
	double sigma = 0;
	double mu = mean( x );
	for( size_t i = 0; i < x.length(); i++ ) {
		sigma += pow( x( i ) - mu, 2 );
	}
	return sigma / ( x.length() - 1 );
}


Signal::TSignal<double> Signal::deleted_averages( const TSignal<double> &x ) {
	double mu = mean( x );
	size_t M = x.length();
	size_t M_1 = M - 1;
	double Mxmu = M * mu;

	Signal::TSignal<double> mu_del( M );
	for( size_t i = 0; i < M; i++ ) {
		mu_del( i ) = ( Mxmu - x( i ) ) / ( M_1 );
	}
	return mu_del;
}


/// Divides signal up into blocks and returns mean of each block
Signal::TSignal<double> Signal::blocked_means( const Signal::TSignal<double> &x, const size_t blocksize ) {
	if( x.length() % blocksize ){
		throw "Block size must divide the total number of samples";
	}
	size_t nblocks = x.length() / blocksize;
	size_t jmax;
	Signal::TSignal<double> mu_b( nblocks );
	for( size_t i = 0; i < nblocks; i++ ) {
		jmax = ( i+1 ) * blocksize;
		mu_b( i )=0;
		for( size_t j = i * blocksize; j < jmax; j++ ) {
			mu_b( i ) += x( j );
		}
	}
	return mu_b/blocksize;
}


double Signal::skewness( const Signal::TSignal<double> &x ) {
	double s = 0;
	double v = 0;
	double mu = mean( x );
	double dx; double dv;
	size_t M=x.length();
	for( size_t i = 0; i < M; i++ ) {
		dx = x( i ) - mu;
		dv = dx * dx;
		v += dv;
		s += dv * dx;
	}
	v /= ( M - 1 );
	return s / ( ( M - 1 ) *v*v*v );
}


double Signal::kurtosis( const Signal::TSignal<double> &x ) {
	double v=0;
	double k=0;
	double mu = mean( x );
	double dx;
	size_t M = x.length();
	for( size_t i=0; i<M; i++ ) {
		dx = x( i ) - mu;
		dx *= dx;
		v += dx;
		dx *= dx;
		k += dx;
	}
	v /= ( M - 1 );
	return k / ( ( M - 1 ) *v*v*v*v ) - 3;
}


double Signal::blocked_var( const Signal::TSignal<double> &x, const size_t blocksize ) {
	return var( blocked_means( x, blocksize ) );
}


double Signal::blocked_mean( const Signal::TSignal<double> &x, const size_t blocksize ) {
	return mean( blocked_means( x, blocksize ) );
}


double Signal::jackknife_mean( const Signal::TSignal<double> &x, const size_t blocksize ) {
	return mean( deleted_averages( blocked_means(x,blocksize ) ) );
}


double Signal::jackknife_var( const Signal::TSignal<double> &x, const size_t blocksize ) {
	return var( deleted_averages( blocked_means( x, blocksize ) ) );
}


double Signal::bootstrap_mean( const Signal::TSignal<double> &x, const size_t blocksize ) {
	const auto& mu_b = blocked_means( x,blocksize );
	size_t nblocks = x.length() / blocksize;
	Signal::TSignal<double> bootstrap_ensemble( nblocks );
	/// Smaller ensemble of statistically uncorrelated means
	for( size_t i = 0; i < nblocks; i++ ) {
		bootstrap_ensemble( i ) = mu_b( randInt( nblocks ) );
	}
	return mean( bootstrap_ensemble );
}


double Signal::bootstrap_var( const Signal::TSignal<double> &x, const size_t blocksize ){
	const auto& mu_b = blocked_means( x, blocksize );
	size_t nblocks = x.length() / blocksize;
	Signal::TSignal<double> bootstrap_ensemble( nblocks );
	/// Smaller ensemble of statistically uncorrelated means
	for( size_t i = 0; i < nblocks; i++ ) {
		bootstrap_ensemble( i ) = mu_b( randInt(nblocks ) );
	}
	return var( bootstrap_ensemble );
}
