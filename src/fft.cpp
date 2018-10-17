//John Stanco 9.12.18


#include "signal.h"

#include </usr/local/include/fftw3.h>

/**
 * 		Implementation of Signal Processing/FFT interface for TSignal<double> and TSignal<cx_double>
 * 		Makefile links with to installed FFTW3 library usr/local/lib/libfftw3.a
 **/


template<class T>
Signal::TSignal<T>& Signal::FFT_base::shift( Signal::TSignal<T> &f ) {
	auto left = f.length() / 2;
 	auto right = f.length() - left;
	auto n_bytes = sizeof( T );

	auto *tmp = ( T* )malloc( n_bytes*left );
	memmove( tmp,f.data()+right,left*n_bytes );
	memmove( f.data()+left,f.data(),right*n_bytes );
	memmove( f.data(),tmp,left*n_bytes );
	return f;
}


Signal::TSignal<cx_double> Signal::FFT::c2c( Signal::TSignal<cx_double> &f ) {
	TSignal<cx_double> s{ f.length() };
	auto arrsize = f.length()*sizeof( fftw_complex );
	auto *in = ( fftw_complex* )fftw_malloc( arrsize );
	auto p = fftw_plan_dft_1d( f.length(), in, in, FFTW_FORWARD, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( s.data(), in, arrsize );
	fftw_free( in );
	return s;
}


Signal::TSignal<double> Signal::FFT::r2r( Signal::TSignal<double> &f ) {
	TSignal<double> s{ f.length() };
	auto arrsize = f.length()*sizeof( double );
	auto *in = ( double* )fftw_malloc( arrsize );
	auto p = fftw_plan_r2r_1d( f.length() , in, in, FFTW_REDFT10, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( s.data(), in, arrsize );
	fftw_free( in );
	return s;
}


void Signal::FFT::inPlace( Signal::TSignal<double> &f ) {
	auto arrsize = f.length()*sizeof( double );
	auto *in = ( double* )fftw_malloc( arrsize );
	auto p = fftw_plan_r2r_1d( f.length(), in, in, FFTW_REDFT10, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( f.data(), in, arrsize );
	fftw_free( in );
}


void Signal::FFT::inPlace( Signal::TSignal<cx_double> &f ) {
	auto arrsize = f.length()*sizeof( fftw_complex );
	auto *in = ( fftw_complex* )fftw_malloc( arrsize );
	auto p = fftw_plan_dft_1d( f.length(), in, in, FFTW_FORWARD, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( f.data(), in, arrsize );
	fftw_free(in);
}


Signal::TSignal<cx_double> Signal::FFT::r2c( Signal::TSignal<double> &f ) {
	Signal::TSignal<cx_double> s{ f.length() };
	auto in_arrsize = f.length()*sizeof( double );
	auto out_arrsize = ( f.length()/2+1 )*sizeof( fftw_complex );
	auto *in = ( double* )fftw_malloc( in_arrsize );
	auto *out = ( fftw_complex* )fftw_malloc( out_arrsize );
	auto *rslt = ( fftw_complex* )fftw_malloc( f.length()*sizeof( fftw_complex ) );
	auto p = fftw_plan_dft_r2c_1d( f.length(), in, out, FFTW_ESTIMATE );

	memcpy( in, f.data(), in_arrsize );
	fftw_execute(p);
	/// FFTW only computes 1st half, use hermiticity for 2nd half | TODO: CHECK
	memcpy( s.data(), out, out_arrsize );

	for( size_t i = f.length()/2+1; i < f.length(); ++i ) {
		s( s.length() - i ) = conj( s( i ) );
	}

	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	fftw_free(rslt);
	return s;
}


Signal::TSignal<cx_double> Signal::IFFT::c2c( Signal::TSignal<cx_double> &f ) {
	Signal::TSignal<cx_double> s{ f.length() };
	auto arrsize = f.length()*sizeof( fftw_complex );
	auto *in = ( fftw_complex* )fftw_malloc( arrsize );
	auto p = fftw_plan_dft_1d( f.length(), in, in, FFTW_BACKWARD, FFTW_ESTIMATE );

	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( s.data(), in, arrsize );
	fftw_free( in );
	return s;
}


Signal::TSignal<double> Signal::IFFT::r2r( Signal::TSignal<double> &f ) {
	Signal::TSignal<double> s{ f.length() };
	auto arrsize = f.length()*sizeof( double );
	auto *in = ( double* )fftw_malloc( arrsize );
	auto p = fftw_plan_r2r_1d( f.length() , in, in, FFTW_REDFT10, FFTW_ESTIMATE );

	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( s.data(), in, arrsize );
	fftw_free( in );
	return s;
}


Signal::TSignal<double> Signal::IFFT::c2r( Signal::TSignal<cx_double> &f ) {
	Signal::TSignal<double> s{ f.length() };
	auto in_arrsize = ( f.length()/2+1 )*sizeof( fftw_complex );
	auto out_arrsize = f.length()*sizeof( double );
	auto *in = ( fftw_complex* )fftw_malloc( in_arrsize );
	auto *out = ( double* )fftw_malloc( out_arrsize );
	auto p = fftw_plan_dft_c2r_1d( f.length(), in, out, FFTW_ESTIMATE );

	memcpy( in, f.data(), in_arrsize );
	fftw_execute(p);
	memcpy( s.data(), out, out_arrsize );

	fftw_destroy_plan( p );
	fftw_free( in );
	fftw_free( out );
	return s;
}


void Signal::IFFT::inPlace( Signal::TSignal<double> &f ) {
	auto arrsize = f.length()*sizeof( double );
	auto *in = ( double* )fftw_malloc( arrsize );
	auto p = fftw_plan_r2r_1d( f.length(), in, in, FFTW_REDFT10, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( f.data(), in, arrsize );
	fftw_free( in );
}


void Signal::IFFT::inPlace( Signal::TSignal<cx_double> &f ) {
	auto arrsize = f.length()*sizeof( fftw_complex );
  auto *in = ( fftw_complex* )fftw_malloc( arrsize );
	auto p = fftw_plan_dft_1d( f.length(), in, in, FFTW_BACKWARD, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( f.data(), in, arrsize );
	fftw_free(in);
}


void assign( fftw_complex z1, const fftw_complex z2 ) {
	z1[0] = z2[0];
	z1[1] = z2[1];
}


void assign_conj( fftw_complex z1, const fftw_complex z2 ) {
	z1[0] = z2[0];
	z1[1] = -z2[1];
}


void multiply( fftw_complex l, fftw_complex r, fftw_complex rslt ) {
	rslt[0] = l[0]*r[0] - l[1]*r[1];
	rslt[1] = l[0]*r[1] + l[1]*r[0];
}


void mult_arr( fftw_complex *l, fftw_complex *r, fftw_complex *rslt, const size_t len ) {
	for( size_t i = 0; i < len; ++i ) {
		multiply(l[i],r[i],rslt[i]);
	}
}


void conj( fftw_complex z ){ z[1] *= -1; }


void conj_arr( fftw_complex *f, const size_t len ) {
	for( size_t i = 0; i < len; ++i ) { conj( f[i] ); }
}


/// DCT-II - input data must be real and symmetric
Signal::TSignal<double> Signal::fft( const Signal::TSignal<double> &f ){
	Signal::TSignal<double> s{ f.length() };
	auto arrsize = f.length()*sizeof( double );
	auto *in = ( double* )fftw_malloc( arrsize );
	auto p = fftw_plan_r2r_1d( f.length(), in, in, FFTW_REDFT10, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute(p);
	memcpy( s.data(), in, arrsize );
	fftw_destroy_plan( p );
	fftw_free( in );
	return s;
}


Signal::TSignal<double> Signal::ifft( const Signal::TSignal<double> &f ) {
	Signal::TSignal<double> s{ f.length() };
	auto arrsize = f.length()*sizeof( double );
	auto *in = ( double* )fftw_malloc( arrsize );
	auto p = fftw_plan_r2r_1d( f.length(), in, in, FFTW_REDFT01, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute(p);
	memcpy( s.data(), in, arrsize );
	fftw_destroy_plan( p );
	fftw_free( in );
	return s;
}


/// IDCT-II - input data must be real and symmetric
Signal::TSignal<cx_double> Signal::fft_r2c( const Signal::TSignal<double> &f ) {
	Signal::TSignal<cx_double> s{ f.length() };
	auto in_arrsize = f.length()*sizeof( double );
	auto out_arrsize = (f.length()/2+1)*sizeof( fftw_complex );
	auto *in = ( double* )fftw_malloc(in_arrsize);
	auto *out = ( fftw_complex* )fftw_malloc( out_arrsize );
	auto *rslt = ( fftw_complex* )fftw_malloc( f.length()*sizeof( fftw_complex ) );
	auto p = fftw_plan_dft_r2c_1d( f.length(), in, out, FFTW_ESTIMATE );

	memcpy( in, f.data(), in_arrsize );
	fftw_execute(p);
	memcpy( s.data(), out, out_arrsize );

	for( size_t i = f.length() / 2 + 1; i < f.length(); ++i ) {
		s( s.length() - i ) = conj( s( i ) );
	}

	fftw_destroy_plan( p );
	fftw_free( in );
	fftw_free( out );
	fftw_free( rslt );
	return s;
}


Signal::TSignal<double> Signal::ifft_c2r( const Signal::TSignal<cx_double> &f ) {
	Signal::TSignal<double> s{ f.length() };
	auto in_arrsize = ( f.length() / 2 + 1 )*sizeof( fftw_complex );
	auto out_arrsize = f.length()*sizeof( double );
	auto *in = ( fftw_complex* )fftw_malloc( in_arrsize );
	auto *out = ( double* )fftw_malloc( out_arrsize );
	auto p = fftw_plan_dft_c2r_1d( f.length(), in, out, FFTW_ESTIMATE );

	memcpy( in, f.data(), in_arrsize );
	fftw_execute( p );
	memcpy( s.data(), out, out_arrsize );
	fftw_destroy_plan( p );
	fftw_free( in );
	fftw_free( out );
	return s;
}


Signal::TSignal<double> Signal::convolve( const Signal::TSignal<double> &l, const Signal::TSignal<double> &r ) {
	if( l.length() != r.length() ) { throw "Signals must be of same length"; }

	auto n = l.length();
	auto N = n / 2 + 1;
	auto in_arrsize = n*sizeof( double );
	auto out_arrsize = N*sizeof( fftw_complex );

	TSignal<double> s{ n };

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


Signal::TSignal<double> Signal::crosscorr( const Signal::TSignal<double> &l, const Signal::TSignal<double> &r ){
	return Signal::convolve( l, r );
}


Signal::TSignal<double> Signal::autocorr( const Signal::TSignal<double> &f ){
	return Signal::crosscorr( f, f );
}


Signal::TSignal<cx_double> Signal::fft( const TSignal<cx_double> &f ) {
	Signal::TSignal<cx_double> s{ f.length() };
	auto arrsize = sizeof( fftw_complex )*f.length();
	auto *in = ( fftw_complex* )fftw_malloc( arrsize );
	auto p = fftw_plan_dft_1d( f.length(), in, in, FFTW_FORWARD, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( s.data(), in, arrsize );
	fftw_destroy_plan( p );
	fftw_free( in );
	return s;
}


Signal::TSignal<cx_double> Signal::ifft( const TSignal<cx_double> &f ) {
	Signal::TSignal<cx_double> s{ f.length() };
	auto arrsize = f.length()*sizeof( fftw_complex );
	auto *in = ( fftw_complex* )fftw_malloc( arrsize );
	auto p = fftw_plan_dft_1d( f.length(), in, in, FFTW_BACKWARD, FFTW_ESTIMATE );
	memcpy( in, f.data(), arrsize );
	fftw_execute( p );
	memcpy( s.data(), in, arrsize );
	fftw_destroy_plan( p );
	fftw_free( in );
	return s;
}


Signal::TSignal<cx_double> Signal::convolve( const TSignal<cx_double> &l, const TSignal<cx_double> &r ) {
  if( l.length() != r.length() ) { throw "Signals must be of same length"; }

	auto N = l.length();
	auto arrsize = N*sizeof( fftw_complex );
	Signal::TSignal<cx_double> s{ N };

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


Signal::TSignal<cx_double> Signal::crosscorr( const Signal::TSignal<cx_double> &l, const Signal::TSignal<cx_double> &r ) {
  if( l.length() != r.length() ) { throw "Signals must be of same length"; }

	auto N = l.length();
	auto arrsize = N*sizeof( fftw_complex );
	Signal::TSignal<cx_double> s{ N };

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


Signal::TSignal<cx_double> Signal::autocorr( const Signal::TSignal<cx_double> &f ) {
	return crosscorr( f, f );
}


template<class T>
Signal::TSignal<T>& Signal::fftshift( Signal::TSignal<T>& f ) {
	size_t left = f.length() / 2;
	size_t right = f.length() - left;
	size_t n_bytes = sizeof( T );

	T *tmp = ( T* )malloc( n_bytes*left );
	memmove( tmp,f.data() + right, left*n_bytes );
	memmove( f.data() + left, f.data(), right*n_bytes );
	memmove( f.data(), tmp, left*n_bytes );
	return f;
}
