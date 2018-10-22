//created by John Stanco 8.10.18

/*

 Interface to implement statistical averaging procedures
 on time-series signals from interface Signal<T> in namespace Signal.

 Interface and implementation can be found in files signal.h and signal.cpp

 */

#include "signal.hpp"


#ifndef SIGNAL_STATS_H
#define SIGNAL_STATS_H


/// Implements blocking, jacknife, bootstrapping analysis of real-valued random variables
template<typename Container>
double mean( const Container& X ) {
	double mu = 0;
	for( auto& x : X ){
		mu += x;
	}
	return mu / X.size();
}


template<typename Container>
double var( const Container& X ) {
	if( X.size() < 2 ) { return 0; }
	double sigma = 0;
	double mu = mean( X );
	for( auto& x : X ) {
		sigma += pow( x - mu, 2 );
	}
	return sigma / ( X.size() - 1 );
}


template<typename Container>
TSignal<double> deleted_averages( const Container &X ) {
	double mu = mean( X );
	size_t M_1 = X.size() - 1;
	double Mxmu = X.size() * mu;

	auto mu_del = Container( X.size() );
	for( size_t i = 0; i < X.size(); ++i ) {
		mu_del[i] = ( Mxmu - X[i] ) / ( M_1 );
	}
	return mu_del;
}


/// Divides signal up into blocks and returns mean of each block
template<typename Container>
Container blocked_means( const Container &X, const size_t blocksize ) {
	if( X.size() % blocksize ){
		throw "Block size must divide the total number of samples";
	}
	size_t nblocks = X.size() / blocksize;
	size_t jmax;
	auto mu_b = Container( nblocks );
	for( size_t i = 0; i < nblocks; ++i ) {
		jmax = ( i+1 ) * blocksize;
		mu_b[i]=0;
		for( size_t j = i * blocksize; j < jmax; ++j ) {
			mu_b[i] += X[j];
		}
	}
	return mu_b/blocksize;
}


template<typename Container>
double skewness( const Container &X ) {
	double s = 0;
	double v = 0;
	double mu = mean( X );
	double dx;
	for( auto& x : X ) {
		dx = x - mu;
		dx *= dx;
		v += dx;
		dx *= dx;
		s += dx;
		// v += pow( x[i] - mu, 2 );
		// s += pow( x[i] - mu, 3 );
	}
	v /= ( X.size() - 1 );
	return s / ( ( X.size() - 1 ) *v*v*v );
}


template<typename Container>
double kurtosis( const Container &X ) {
	double v=0;
	double k=0;
	double mu = mean( X );
	double dx;
	for( auto& x : X ) {
		dx = x - mu;
		dx *= dx;
		v += dx;
		dx *= dx;
		k += dx;
		// v += pow( x[i] - mu, 2 );
		// k += pow( x[i] - mu, 3 );
	}
	v /= ( X.size() - 1 );
	return k / ( ( X.size() - 1 ) *v*v*v*v ) - 3;
}


template<typename Container>
double blocked_var( const Container &X, const size_t blocksize ) {
	return var( blocked_means( X, blocksize ) );
}


template<typename Container>
double blocked_mean( const Container &X, const size_t blocksize ) {
	return mean( blocked_means( X, blocksize ) );
}


template<typename Container>
double jackknife_mean( const Container &X, const size_t blocksize ) {
	return mean( deleted_averages( blocked_means( X, blocksize ) ) );
}


template<typename Container>
double jackknife_var( const Container &X, const size_t blocksize ) {
	return var( deleted_averages( blocked_means( X, blocksize ) ) );
}


template<typename Container>
double bootstrap_mean( const Container &X, const size_t blocksize ) {
	const auto& mu_b = blocked_means( X, blocksize );
	size_t nblocks = X.size() / blocksize;
	auto bootstrap_ensemble = Container( nblocks );
	/// Smaller ensemble of statistically uncorrelated means
	for( size_t i = 0; i < nblocks; ++i ) {
		bootstrap_ensemble[i] = mu_b[randInt( nblocks )];
	}
	return mean( bootstrap_ensemble );
}


template<typename Container>
double bootstrap_var( const Container &X, const size_t blocksize ){
	const auto& mu_b = blocked_means( X, blocksize );
	size_t nblocks = X.size() / blocksize;
	auto bootstrap_ensemble = Container( nblocks );
	/// Smaller ensemble of statistically uncorrelated means
	for( size_t i = 0; i < nblocks; ++i ) {
		bootstrap_ensemble[i] = mu_b[randInt( nblocks )];
	}
	return var( bootstrap_ensemble );
}


template<typename Container>
void normalize_std_inpl( Container& X ) {
	double mu = mean( X );
	double stdev = sqrt( var( X ) );
	for( auto& x : X ) {
		x = ( x - mu ) / stdev;
	}
}


template<typename Container>
Container normalize_std( const Container& X ) {
	auto norm = Container( X.size() );
	double mu = mean( X );
	double stdev = sqrt( var( X ) );
	for( size_t i = 0; i < X.size(); ++i ) {
		norm[i] = ( X[i] - mu ) / stdev;
	}
	return norm;
}


#endif /* SIGNAL_STATS_H */
