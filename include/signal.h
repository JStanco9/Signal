//Created by John Stanco on 2/5/18


/**
 * 		Interface for the class TSignal<T>:
 * 		- provides storage time-series signals/random variables.
 * 		- Supports computation of Fast Signal-Processing routines.
 * 			FFT, IFFT, Cross-Corellation, Auto-Correlation, and Convolution
 * 			through implementation of FFTW3 library for complex and real-valued signals.
 * 		- Supports Statistical Analysis of signals including
 * 			blocking, jackknife, and bootstrap analysis for real signals.
 * 			includes calculation of statistical moments up to 4.
 **/


#include <complex>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fftw3.h>


#ifndef SIGNAL_H
#define SIGNAL_H
#ifndef SIGNAL_TSIGNAL_H
#define SIGNAL_TSIGNAL_H


typedef std::complex<double> cx_double;
extern double pi;

namespace Signal{

	template<class T, class A> class TSignal;
	template<class T, class A> void swap( TSignal<T, A> &a, TSignal<T, A> &b );

	template<class T, class A = std::allocator<T>>
	class TSignal {
	private:
		T *dat;
		size_t len;
		size_t cap;
		A alloc;

		friend void swap<T, A>( TSignal& a, TSignal& b );
		void expand();
	public:
		//TSignal();
		TSignal( const size_t n = 0 );
		TSignal( const T arr[] );
		TSignal( const TSignal &other );
		TSignal& operator=( TSignal rhs );
		~TSignal();

		size_t size() const;
		/// accessing underlying array
		T* data();
		T* data() const;

		void push_back( const T &item );

		T& operator()(const size_t index);
		T operator()(const size_t index) const;
		T& operator[](const size_t index);
		T operator[](const size_t index) const;
		bool operator==( const TSignal &rhs ) const;
		bool operator!=( const TSignal &rhs ) const;

		class iterator {
			friend class TSignal;
			T *_elem;
		public:
			iterator( T *elem );
			iterator& operator++();
			//iterator& operator++(int);
			T& operator*();
			bool operator==( const iterator &rhs );
			bool operator!=( const iterator &rhs );
		};

		class const_iterator {
			friend class TSignal;
			T *_elem;
		public:
			const_iterator( T *elem );
			const_iterator& operator++();
			//const_iterator& operator++(int);
			T& operator*();
			bool operator==( const const_iterator &rhs );
			bool operator!=( const const_iterator &rhs );
		};

		iterator begin();
		const_iterator begin() const;
		iterator back();
		const_iterator back() const;
		iterator end();
		const_iterator end() const;

	};


	class TSignalOutOfRange : public std::out_of_range {
	public:
		TSignalOutOfRange() : std::out_of_range{ "Index out range." } {}
	};


	template<class T, class A>
	void swap( TSignal<T, A> &a, TSignal<T, A> &b ) {
		std::swap( a.dat, b.dat );
		std::swap( a.alloc, b.alloc );
		std::swap( a.len, b.len );
		std::swap( a.cap, b.cap );
	}

	/*
	template<class T, class A>
	TSignal<T, A>::TSignal() : alloc{ A{} }, len{ 0 }, cap{ 0 } {
		dat = alloc.allocate( 1 );
	}
	*/

	template<class T, class A>
	TSignal<T, A>::TSignal( const size_t n ) : alloc{ A{} }, len{ n }, cap{ n } {
		dat = alloc.allocate( n + 1 );
	}


	template<class T, class A>
	TSignal<T, A>::TSignal( const TSignal<T, A>& other ) :
	alloc{ other.alloc }, len{ other.len }, cap{ other.cap } {
		dat = alloc.allocate( other.cap + 1 );
		memcpy( dat, other.dat, (cap + 1) * sizeof( T ));
	}


	template<class T, class A>
	TSignal<T, A>& TSignal<T, A>::operator=( TSignal rhs ) {
	 swap( *this, rhs );
	 return *this;
	}


	template<class T, class A>
	TSignal<T, A>::~TSignal() {
		alloc.deallocate( dat, cap + 1 );
	}


	template<class T, class A>
	size_t TSignal<T, A>::size() const { return len; }


	template<class T, class A>
	T* TSignal<T, A>::data() { return dat; }


	template<class T, class A>
	T* TSignal<T, A>::data() const { return dat; }


	template<class T, class A>
	void TSignal<T, A>::expand() {
		int newCap = ( cap == 0 )? cap+1 : cap * 2;
		auto *newDat = alloc.allocate( newCap + 1 );
		memmove( newDat, dat, len * sizeof( T ) );
		auto *tmp = dat;
		dat = newDat;
		alloc.deallocate( tmp, cap + 1 );
		cap = newCap;
	}


	template<class T, class A>
	void TSignal<T, A>::push_back( const T &item ) {
		if( len == cap ) { expand(); }
		dat[len++] = item;
	}


	template<class T, class A>
	T& TSignal<T, A>::operator()( const size_t index ) {
		if( index >= len ) { throw TSignalOutOfRange(); }
		return dat[index];
	}


	template<class T, class A>
	T TSignal<T, A>::operator()( const size_t index ) const {
		if( index >= len ) { throw TSignalOutOfRange(); }
		return dat[index];
	}


	template<class T, class A>
	T& TSignal<T, A>::operator[](const size_t index) {
		return this->operator()( index );
	}


	template<class T, class A>
	T TSignal<T, A>::operator[](const size_t index) const {
		return this->operator()( index );
	}


	template<class T, class A>
	bool TSignal<T, A>::operator==( const TSignal &rhs ) const {
		if( len != rhs.len ) { return false; }
		auto itl = begin();
		auto itr = rhs.begin();
		for( ; itl != end() && itr != rhs.end(); ++itl, ++itr ) {
			if( *itl != *itr ) { return false; }
		}
		return true;
	}


	template<class T, class A>
	bool TSignal<T, A>::operator!=( const TSignal &rhs ) const {
		if( len != rhs.len ) { return true; }
		auto itl = begin();
		auto itr = rhs.begin();
		for( ; itl != end() && itr != rhs.end(); ++itl, ++itr ) {
			if( *itl != *itr ) { return true; }
		}
		return false;
	}


	template<class T, class A>
	TSignal<T, A>::iterator::iterator( T *elem ) : _elem{ elem } {}


	template<class T, class A>
	class TSignal<T, A>::iterator& TSignal<T, A>::iterator::operator++() {
		++_elem;
		return *this;
	}


	template<class T, class A>
	T& TSignal<T, A>::iterator::operator*() {
		return *_elem;
	}


	template<class T, class A>
	bool TSignal<T, A>::iterator::operator==( const iterator& rhs ) {
		return _elem = rhs._elem;
	}


	template<class T, class A>
	bool TSignal<T, A>::iterator::operator!=( const iterator& rhs ) {
		return _elem != rhs._elem;
	}


	template<class T, class A>
	TSignal<T, A>::const_iterator::const_iterator( T *elem ) : _elem{ elem } {}


	template<class T, class A>
	class TSignal<T, A>::const_iterator& TSignal<T, A>::const_iterator::operator++() {

		++_elem;
		return *this;
	}


	template<class T, class A>
	T& TSignal<T, A>::const_iterator::operator*() {
		return *_elem;
	}


	template<class T, class A>
	bool TSignal<T, A>::const_iterator::operator==( const const_iterator& rhs ) {
		return _elem = rhs._elem;
	}


	template<class T, class A>
	bool TSignal<T, A>::const_iterator::operator!=( const const_iterator& rhs ) {
		return _elem != rhs._elem;
	}


	template<class T, class A>
	class TSignal<T, A>::iterator TSignal<T, A>::begin() {
		return iterator{ &dat[0] };
	}


	template<class T, class A>
	class TSignal<T, A>::const_iterator TSignal<T, A>::begin() const {
		return const_iterator{ &dat[0] };
	}


	template<class T, class A>
	class TSignal<T, A>::iterator TSignal<T, A>::back() {
		return iterator{ &dat[len - 1] };
	}


	template<class T, class A>
	class TSignal<T, A>::const_iterator TSignal<T, A>::back() const {
		return const_iterator{ &dat[len - 1] };
	}


	template<class T, class A>
	class TSignal<T, A>::iterator TSignal<T, A>::end() {
		/// dat array always contains 1 extra element
		return iterator{ &dat[len] };
	}


	template<class T, class A>
	class TSignal<T, A>::const_iterator TSignal<T, A>::end() const {
		/// dat array always contains 1 extra element
		return const_iterator{ &dat[len] };
	}

	typedef TSignal<cx_double> cx_signal;
	typedef TSignal<double> signal;

}


template<template<class, class> class Container, class V, class A, class Stream>
void print( const Container<V, A> &s, Stream& os ){
	for( auto& x : s ) {
		os << x << "\n";
	}
}


template<template<class, class> class Container, class A, class Stream>
void print( const Container<cx_double, A> &s, Stream& os ) {
	for( auto& z : s ) {
		os << std::left << std::setw(24) << std::setprecision(16) << z.real() << std::left << z.imag() << "\n";
	}
}


template<template<class, class> class Container, class V, class A>
void print( const Container<V, A> &s, std::basic_ostream<char> &os = std::cout ) {
	for( auto& x : s ) {
		os << x << "\n";
	}
}


template<template<class, class> class Container, class A>
void print( const Container<cx_double, A> &s, std::basic_ostream<char> &os = std::cout ) {
	for( auto& z : s ) {
		os << std::left << std::setw(24) << std::setprecision(16) << z.real() << std::left << z.imag() << "\n";
	}
}


template<template<class, class> class Container, class A, class Stream>
void printReal( const Container<cx_double, A> &s, Stream &os ) {
	if( !os.is_open() ) { throw "Attempted to print to unopen file stream"; }
	for( auto& z : s ) {
		os << std::left << std::setprecision(16) << z.real() << "\n";
	}
}


template<template<class, class> class Container, class A>
void printReal( const Container<cx_double, A> &s, std::basic_ostream<char> &os = std::cout ) {
	for( auto& z : s ) {
		os << std::left << std::setprecision(16) << z.real() << "\n";
	}
}


template<template<class, class> class Container, class A, class Stream>
void printImag( const Container<cx_double, A> &s, Stream &os ) {
	if( !os.is_open() ) { throw "Attempted to print to unopen file stream"; }
	for( auto& z : s ) {
		os << std::left << std::setprecision(16) << z.imag() << "\n";
	}
}


template<template<class, class> class Container, class A>
void printImag( const Container<cx_double, A> &s, std::basic_ostream<char> &os = std::cout ) {
	for( auto& z : s ) {
		os << std::left << std::setprecision(16) << z.imag() << "\n";
	}
}


#endif /* SIGNAL_TSIGNAL_H */
#ifndef MYRAND_H
#define MYRAND_H

int flip( double x = .5 );
int randInt( int, int );
int randInt( int );
double min( double, double );
double randDouble( double, double );
double plusMinus( double );
double randNorm();
double randNorm( double, double );

#endif /* MYRAND_H */
#ifndef SIGNAL_STATS_H
#define SIGNAL_STATS_H


/// Implements blocking, jacknife, bootstrapping analysis of real-valued random variables
template<class Container>
double mean( const Container& X ) {
	double mu = 0;
	for( auto& x : X ){
		mu += x;
	}
	return mu / X.size();
}


template<class Container>
double var( const Container& X ) {
	if( X.size() < 2 ) { return 0; }
	double sigma = 0;
	double mu = mean( X );
	for( auto& x : X ) {
		sigma += pow( x - mu, 2 );
	}
	return sigma / ( X.size() - 1 );
}


template<class Container>
Container deleted_averages( const Container &X ) {
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
template<class Container>
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


template<class Container>
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


template<class Container>
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


template<class Container>
double blocked_var( const Container &X, const size_t blocksize ) {
	return var( blocked_means( X, blocksize ) );
}


template<class Container>
double blocked_mean( const Container &X, const size_t blocksize ) {
	return mean( blocked_means( X, blocksize ) );
}


template<class Container>
double jackknife_mean( const Container &X, const size_t blocksize ) {
	return mean( deleted_averages( blocked_means( X, blocksize ) ) );
}


template<class Container>
double jackknife_var( const Container &X, const size_t blocksize ) {
	return var( deleted_averages( blocked_means( X, blocksize ) ) );
}


template<class Container>
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


template<class Container>
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


template<class Container>
void normalize_std_inpl( Container& X ) {
	double mu = mean( X );
	double stdev = sqrt( var( X ) );
	for( auto& x : X ) {
		x = ( x - mu ) / stdev;
	}
}


template<class Container>
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
#ifndef SIGNAL_FFT_H
#define SIGNAL_FFT_H


class FFT_base {
	FFT_base( const FFT_base& other ) = delete;
	FFT_base& operator=( FFT_base other ) = delete;
public:
	template<template<class, class> class Container, class T, class A>
	Container<T, A>& shift( Container<T, A>& );
};


class FFT : public FFT_base {
	FFT( const FFT &other ) = delete;
	FFT& operator=( FFT other ) = delete;
public:
	template<template<class, class> class Container, class A>
	static Container<cx_double, A> c2c( Container<cx_double, A> &f );
	template<template<class, class> class Container, class A>
	static Container<double, A> r2r( Container<double, A> &f );
	template<template<class, class> class Container, class A>
	static Container<cx_double, A> r2c( Container<double, A> &f );
	template<template<class, class> class Container, class A>
	static void inPlace( Container<cx_double, A>& );
	template<template<class, class> class Container, class A>
	static void inPlace( Container<double, A>& );
};


class IFFT : public FFT_base {
	IFFT( const IFFT &other ) = delete;
	IFFT& operator=( IFFT other ) = delete;
public:
	template<template<class, class> class Container, class A>
	static Container<cx_double, A> c2c( Container<cx_double, A> &f );
	template<template<class, class> class Container, class A>
	static Container<double, A> r2r( Container<double, A> &f );
	template<template<class, class> class Container, class A>
	static Container<double, A> c2r( Container<cx_double, A> &f );
	template<template<class, class> class Container, class A>
	static void inPlace( Container<cx_double, A>& );
	template<template<class, class> class Container, class A>
	static void inPlace( Container<double, A>& );
};


template<template<class, class> class Container, class T, class A>
Container<T, A>& FFT_base::shift( Container<T, A> &f ) {
	auto left = f.size() / 2;
 	auto right = f.size() - left;
	auto n_bytes = sizeof( T );

	auto *tmp = ( T* )malloc( n_bytes*left );
	memmove( tmp,f.data()+right,left*n_bytes );
	memmove( f.data()+left,f.data(),right*n_bytes );
	memmove( f.data(),tmp,left*n_bytes );
	free( tmp ); //might not be necessary
	return f;
}


template<template<class, class> class Container, class A>
Container<cx_double, A> FFT::c2c( Container<cx_double, A> &f ) {
	auto s = Container<cx_double, A>( f.size() );
	auto arrsize = f.size()*sizeof( fftw_complex );
	auto *in = ( fftw_complex* )fftw_malloc( arrsize );
	auto p = fftw_plan_dft_1d( f.size(), in, in, FFTW_FORWARD, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( s.data(), in, arrsize );
	fftw_free( in );
	return s;
}


template<template<class, class> class Container, class A>
Container<double, A> FFT::r2r( Container<double, A> &f ) {
	auto s = Container<double, A>( f.size() );
	auto arrsize = f.size()*sizeof( double );
	auto *in = ( double* )fftw_malloc( arrsize );
	auto p = fftw_plan_r2r_1d( f.size() , in, in, FFTW_REDFT10, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( s.data(), in, arrsize );
	fftw_free( in );
	return s;
}


template<template<class, class> class Container, class A>
void FFT::inPlace( Container<double, A> &f ) {
	auto arrsize = f.size()*sizeof( double );
	auto *in = ( double* )fftw_malloc( arrsize );
	auto p = fftw_plan_r2r_1d( f.size(), in, in, FFTW_REDFT10, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( f.data(), in, arrsize );
	fftw_free( in );
}


template<template<class, class> class Container, class A>
void FFT::inPlace( Container<cx_double, A> &f ) {
	auto arrsize = f.size()*sizeof( fftw_complex );
	auto *in = ( fftw_complex* )fftw_malloc( arrsize );
	auto p = fftw_plan_dft_1d( f.size(), in, in, FFTW_FORWARD, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( f.data(), in, arrsize );
	fftw_free(in);
}


template<template<class, class> class Container, class A>
Container<cx_double, A> FFT::r2c( Container<double, A> &f ) {
	auto s = Container<cx_double, A>( f.size() );
	auto in_arrsize = f.size()*sizeof( double );
	auto out_arrsize = ( f.size()/2+1 )*sizeof( fftw_complex );
	auto *in = ( double* )fftw_malloc( in_arrsize );
	auto *out = ( fftw_complex* )fftw_malloc( out_arrsize );
	auto *rslt = ( fftw_complex* )fftw_malloc( f.size()*sizeof( fftw_complex ) );
	auto p = fftw_plan_dft_r2c_1d( f.size(), in, out, FFTW_ESTIMATE );

	memcpy( in, f.data(), in_arrsize );
	fftw_execute(p);
	/// FFTW only computes 1st half, use hermiticity for 2nd half | TODO: CHECK
	memcpy( s.data(), out, out_arrsize );
	for( size_t i = f.size()/2+1; i < f.size(); ++i ) {
		s( s.size() - i ) = conj( s[i] );
	}
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	fftw_free(rslt);
	return s;
}


template<template<class, class> class Container, class A>
Container<cx_double, A> IFFT::c2c( Container<cx_double, A> &f ) {
	auto s = Container<cx_double, A>( f.size() );
	auto arrsize = f.size()*sizeof( fftw_complex );
	auto *in = ( fftw_complex* )fftw_malloc( arrsize );
	auto p = fftw_plan_dft_1d( f.size(), in, in, FFTW_BACKWARD, FFTW_ESTIMATE );

	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( s.data(), in, arrsize );
	fftw_free( in );
	return s;
}


template<template<class, class> class Container, class A>
Container<double, A> IFFT::r2r( Container<double, A> &f ) {
	auto s = Container<double, A>( f.size() );
	auto arrsize = f.size()*sizeof( double );
	auto *in = ( double* )fftw_malloc( arrsize );
	auto p = fftw_plan_r2r_1d( f.size() , in, in, FFTW_REDFT10, FFTW_ESTIMATE );

	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( s.data(), in, arrsize );
	fftw_free( in );
	return s;
}


template<template<class, class> class Container, class A>
Container<double, A> IFFT::c2r( Container<cx_double, A> &f ) {
	auto s = Container<double, A>( f.size() );
	auto in_arrsize = ( f.size()/2+1 )*sizeof( fftw_complex );
	auto out_arrsize = f.size()*sizeof( double );
	auto *in = ( fftw_complex* )fftw_malloc( in_arrsize );
	auto *out = ( double* )fftw_malloc( out_arrsize );
	auto p = fftw_plan_dft_c2r_1d( f.size(), in, out, FFTW_ESTIMATE );

	memcpy( in, f.data(), in_arrsize );
	fftw_execute(p);
	memcpy( s.data(), out, out_arrsize );

	fftw_destroy_plan( p );
	fftw_free( in );
	fftw_free( out );
	return s;
}


template<template<class, class> class Container, class A>
void IFFT::inPlace( Container<double, A> &f ) {
	auto arrsize = f.size()*sizeof( double );
	auto *in = ( double* )fftw_malloc( arrsize );
	auto p = fftw_plan_r2r_1d( f.size(), in, in, FFTW_REDFT10, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( f.data(), in, arrsize );
	fftw_free( in );
}


template<template<class, class> class Container, class A>
void IFFT::inPlace( Container<cx_double, A> &f ) {
	auto arrsize = f.size()*sizeof( fftw_complex );
  auto *in = ( fftw_complex* )fftw_malloc( arrsize );
	auto p = fftw_plan_dft_1d( f.size(), in, in, FFTW_BACKWARD, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( f.data(), in, arrsize );
	fftw_free(in);
}


inline void multiply( fftw_complex l, fftw_complex r, fftw_complex rslt ) {
	rslt[0] = l[0]*r[0] - l[1]*r[1];
	rslt[1] = l[0]*r[1] + l[1]*r[0];
}


inline void mult_arr( fftw_complex *l, fftw_complex *r, fftw_complex *rslt, const size_t len ) {
	for( size_t i = 0; i < len; ++i ) {
		multiply(l[i],r[i],rslt[i]);
	}
}


inline void conj( fftw_complex z ){ z[1] *= -1; }


inline void conj_arr( fftw_complex *f, const size_t len ) {
	for( size_t i = 0; i < len; ++i ) { conj( f[i] ); }
}


/// DCT-II - input data must be real and symmetric
template<template<class, class> class Container, class A>
Container<double, A> fft( const Container<double, A> &f ){
	auto s = Container<double, A>( f.size() );
	auto arrsize = f.size()*sizeof( double );
	auto *in = ( double* )fftw_malloc( arrsize );
	auto p = fftw_plan_r2r_1d( f.size(), in, in, FFTW_REDFT10, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute(p);
	memcpy( s.data(), in, arrsize );
	fftw_destroy_plan( p );
	fftw_free( in );
	return s;
}


template<template<class, class> class Container, class A>
Container<double, A> ifft( const Container<double, A> &f ) {
	auto s = Container<double, A>( f.size() );
	auto arrsize = f.size()*sizeof( double );
	auto *in = ( double* )fftw_malloc( arrsize );
	auto p = fftw_plan_r2r_1d( f.size(), in, in, FFTW_REDFT01, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute(p);
	memcpy( s.data(), in, arrsize );
	fftw_destroy_plan( p );
	fftw_free( in );
	return s;
}


/// IDCT-II - input data must be real and symmetric
template<template<class, class> class Container, class A>
Container<cx_double, A> fft_r2c( const Container<double, A> &f ) {
	auto s = Container<cx_double, A>( f.size() );
	auto in_arrsize = f.size()*sizeof( double );
	auto out_arrsize = (f.size()/2+1)*sizeof( fftw_complex );
	auto *in = ( double* )fftw_malloc(in_arrsize);
	auto *out = ( fftw_complex* )fftw_malloc( out_arrsize );
	auto *rslt = ( fftw_complex* )fftw_malloc( f.size()*sizeof( fftw_complex ) );
	auto p = fftw_plan_dft_r2c_1d( f.size(), in, out, FFTW_ESTIMATE );

	memcpy( in, f.data(), in_arrsize );
	fftw_execute(p);
	memcpy( s.data(), out, out_arrsize );

	for( size_t i = f.size() / 2 + 1; i < f.size(); ++i ) {
		s( s.size() - i ) = conj( s( i ) );
	}

	fftw_destroy_plan( p );
	fftw_free( in );
	fftw_free( out );
	fftw_free( rslt );
	return s;
}


template<template<class, class> class Container, class A>
Container<double, A> ifft_c2r( const Container<cx_double, A> &f ) {
	auto s = Container<double, A>( f.size() );
	auto in_arrsize = ( f.size() / 2 + 1 )*sizeof( fftw_complex );
	auto out_arrsize = f.size()*sizeof( double );
	auto *in = ( fftw_complex* )fftw_malloc( in_arrsize );
	auto *out = ( double* )fftw_malloc( out_arrsize );
	auto p = fftw_plan_dft_c2r_1d( f.size(), in, out, FFTW_ESTIMATE );

	memcpy( in, f.data(), in_arrsize );
	fftw_execute( p );
	memcpy( s.data(), out, out_arrsize );
	fftw_destroy_plan( p );
	fftw_free( in );
	fftw_free( out );
	return s;
}


template<template<class, class> class Container, class A>
Container<double, A> convolve( const Container<double, A> &l, const Container<double, A> &r ) {
	if( l.size() != r.size() ) { throw "Signals must be of same size"; }

	auto n = l.size();
	auto N = n / 2 + 1;
	auto in_arrsize = n*sizeof( double );
	auto out_arrsize = N*sizeof( fftw_complex );

	auto s = Container<double, A>( n );

	auto *in = ( double* )fftw_malloc( in_arrsize );
	auto *lout = ( fftw_complex* )fftw_malloc( out_arrsize );
	auto *rout = ( fftw_complex* )fftw_malloc( out_arrsize );
	auto *rslt = ( fftw_complex* )fftw_malloc( out_arrsize );

	auto pl = fftw_plan_dft_r2c_1d( n, in, lout, FFTW_ESTIMATE );
	auto pr = fftw_plan_dft_r2c_1d( n, in, rout, FFTW_ESTIMATE );
	auto pc = fftw_plan_dft_c2r_1d( n, rout, in, FFTW_ESTIMATE );

	memcpy( in, l.data(), in_arrsize );
	fftw_execute( pl );
	memcpy( in, r.data(), in_arrsize );
	fftw_execute( pr );
	mult_arr( lout, rout, rout, N );
	fftw_execute( pc );
	memcpy( s.data(), in, in_arrsize );
	fftw_free( in );
	fftw_free( lout );
	fftw_free( rout );
	fftw_free( rslt );
	fftw_destroy_plan( pl );
	fftw_destroy_plan( pr );
	fftw_destroy_plan( pc );
	return s;
}


template<template<class, class> class Container, class A>
Container<double, A> crosscorr( const Container<double, A> &l, const Container<double, A> &r ){
	return convolve( l, r );
}


template<template<class, class> class Container, class A>
Container<double, A> autocorr( const Container<double, A> &f ){
	return crosscorr( f, f );
}


template<template<class, class> class Container, class A>
Container<cx_double, A> fft( const Container<cx_double, A> &f ) {
	Container<cx_double, A> s( f.size() );
	auto arrsize = sizeof( fftw_complex )*f.size();
	auto *in = ( fftw_complex* )fftw_malloc( arrsize );
	auto p = fftw_plan_dft_1d( f.size(), in, in, FFTW_FORWARD, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( s.data(), in, arrsize );
	fftw_destroy_plan( p );
	fftw_free( in );
	return s;
}


template<template<class, class> class Container, class A>
Container<cx_double, A> ifft( const Container<cx_double, A> &f ) {
	Container<cx_double, A> s( f.size() );
	auto arrsize = f.size()*sizeof( fftw_complex );
	auto *in = ( fftw_complex* )fftw_malloc( arrsize );
	auto p = fftw_plan_dft_1d( f.size(), in, in, FFTW_BACKWARD, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( s.data(), in, arrsize );
	fftw_destroy_plan( p );
	fftw_free( in );
	return s;
}


template<template<class, class> class Container, class A>
Container<cx_double, A> convolve( const Container<cx_double, A> &l, const Container<cx_double, A> &r ) {
  if( l.size() != r.size() ) { throw "Signals must be of same size"; }

	auto N = l.size();
	auto arrsize = N*sizeof( fftw_complex );
	Container<cx_double, A> s( N );

	auto *in = ( fftw_complex* )fftw_malloc( arrsize );
	auto *lout = ( fftw_complex* )fftw_malloc( arrsize );
	auto *rout = ( fftw_complex* )fftw_malloc( arrsize );

	auto pl = fftw_plan_dft_1d( N, in, lout, FFTW_FORWARD, FFTW_ESTIMATE );
	auto pr = fftw_plan_dft_1d( N, in, rout, FFTW_FORWARD, FFTW_ESTIMATE );
	auto pc = fftw_plan_dft_1d( N, in, lout, FFTW_BACKWARD, FFTW_ESTIMATE );

	memcpy( in, l.data(), arrsize );
	fftw_execute( pl );
	memcpy( in, r.data(), arrsize );
	fftw_execute( pr );
	//cblas_zsbmv('U',N,0,FFTW_ONE,...)
	mult_arr( lout, rout, in, N );
	fftw_execute( pc );
	memcpy( s.data(), lout, arrsize );
	fftw_free( in );
	fftw_free( lout );
	fftw_free( rout );
	fftw_destroy_plan( pl );
	fftw_destroy_plan( pr );
	fftw_destroy_plan( pc );
	return s;
}


template<template<class, class> class Container, class A>
Container<cx_double, A> crosscorr( const Container<cx_double, A> &l, const Container<cx_double, A> &r ) {
  if( l.size() != r.size() ) { throw "Signals must be of same size"; }

	auto N = l.size();
	auto arrsize = N*sizeof( fftw_complex );
	Container<cx_double, A> s( N );

	auto *in = ( fftw_complex* )fftw_malloc( arrsize );
	auto *lout = ( fftw_complex* )fftw_malloc( arrsize );
	auto *rout = ( fftw_complex* )fftw_malloc( arrsize );

	auto pl = fftw_plan_dft_1d( N, in, lout, FFTW_FORWARD, FFTW_ESTIMATE );
	auto pr = fftw_plan_dft_1d( N, in, rout, FFTW_FORWARD, FFTW_ESTIMATE );
	auto pc = fftw_plan_dft_1d( N, in, lout, FFTW_BACKWARD, FFTW_ESTIMATE );

	memcpy( in, l.data(), arrsize );
	conj_arr( in, N );
	fftw_execute( pl );
	memcpy( in, r.data(), arrsize );
	fftw_execute( pr );
  //cblas_zsbmv('U',N,0,FFTW_ONE,...) // can be used for fast element-wise vector multiplication
	mult_arr( lout, rout, in, N );
	fftw_execute( pc );
	memcpy( s.data(), lout, arrsize );
	fftw_free( in );
	fftw_free( lout );
	fftw_free( rout );
	fftw_destroy_plan( pl );
	fftw_destroy_plan( pr );
	fftw_destroy_plan( pc );
	return s;
}


template<template<class, class> class Container, class A>
Container<cx_double, A> autocorr( const Container<cx_double, A> &f ) {
	return crosscorr( f, f );
}


template<template<class, class> class Container, class T, class A>
Container<T, A>& fftshift( Container<T, A>& f ) {
	size_t left = f.size() / 2;
	size_t right = f.size() - left;
	size_t n_bytes = sizeof( T );

	T *tmp = ( T* )malloc( n_bytes*left );
	memmove( tmp,f.data() + right, left*n_bytes );
	memmove( f.data() + left, f.data(), right*n_bytes );
	memmove( f.data(), tmp, left*n_bytes );
	return f;
}

#endif /* SIGNAL_FFT_H */
#ifndef SIGNAL_FUNCTOR_H
#define SIGNAL_FUNCTOR_H


template<class U, class V>
class UnaryFunctor {
public:
  virtual V operator()( const U& ) = 0;
  virtual ~UnaryFunctor() {}
};


template<class U, class V>
class BinaryFunctor {
public:
  virtual V operator()( const U&, const U& ) = 0;
  virtual ~BinaryFunctor() {}
};


template<class U, class V>
class UnaryFunctionWrap : public UnaryFunctor<U, V> {
  //wraps function pointer
  V( *f )( const U& );
public:
  UnaryFunctionWrap( V( *function )( const U& ) );
  V operator()( const U &x );
};


template<class U, class V>
UnaryFunctionWrap<U, V>::UnaryFunctionWrap( V( *function )( const U& ) ) : f{ function } {}


template<class U, class V>
V UnaryFunctionWrap<U, V>::operator()( const U &x ) { return f( x ); }


template<class U, class V>
class BinaryFunctionWrap : public BinaryFunctor<U, V> {
  //wraps function pointer
  V( *f )( const U&, const U& );
public:
  BinaryFunctionWrap( V( *function )( const U&, const U& ) );
  V operator()( const U &x, const U &y );
};


template<class U, class V>
BinaryFunctionWrap<U, V>::BinaryFunctionWrap( V( *function )( const U&, const U& ) ) : f{ function } {}


template<class U, class V>
V BinaryFunctionWrap<U, V>::operator()( const U &x, const U &y ) { return f( x, y ); }


template<class U, class V>
Signal::TSignal<V> map( const Signal::TSignal<U> &s, UnaryFunctor<U, V> &f ) {
  Signal::TSignal<V> image{};
  for( auto& x : s ) { image.push_back( f( x ) ); }
  return image;
}


template<class U>
void apply( Signal::TSignal<U> &s, UnaryFunctor<U, U> &f ) {
  for( auto& x : s ) { x = f( x ); }
}


template<class U>
U reduce( const Signal::TSignal<U> &s, BinaryFunctor<U, U> &f, const U& id ) {
  auto acc = id;
  for( auto& x : s ) { acc = f( acc, x ); }
  return acc;
}

#endif /* SIGNAL_FUNCTOR_H */
#endif /* SIGNAL__H */
