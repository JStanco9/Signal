// John Stanco 9.14.18

/*

	Interface for implementing FFT, mapping both complex or real signals to
	their fourier transform - Implementation uses fftw3 library.

*/

#include "signal.h"

#ifndef SIGNAL_FFT_H
#define SIGNAL_FFT_H

class FFT_base {
	FFT_base( const FFT_base& other ) = delete;
	FFT_base& operator=( FFT_base other ) = delete;
public:
	template<class T>
	Signal::TSignal<T>& shift( Signal::TSignal<T>& );
};


class FFT : public FFT_base {
	FFT( const FFT &other ) = delete;
	FFT& operator=( FFT other ) = delete;
public:
	static Signal::TSignal<cx_double> c2c( Signal::TSignal<cx_double> &f );
	static Signal::TSignal<double> r2r( Signal::TSignal<double> &f );
	static Signal::TSignal<cx_double> r2c( Signal::TSignal<double> &f );
	static void in_place( Signal::TSignal<cx_double>& );
	static void in_place( Signal::TSignal<double>& );
};


class IFFT : public FFT_base {
	IFFT( const IFFT &other ) = delete;
	IFFT& operator=( IFFT other ) = delete;
public:
	static Signal::TSignal<cx_double> c2c( Signal::TSignal<cx_double> &f );
	static Signal::TSignal<double> r2r( Signal::TSignal<double> &f );
	static Signal::TSignal<double> c2r( Signal::TSignal<cx_double> &f );
	static void in_place( Signal::TSignal<cx_double>& );
	static void in_place( Signal::TSignal<double>& );
};


namespace Signal {

	TSignal<double> convolve( const TSignal<double> &l, const TSignal<double> &r );
	TSignal<double> crosscorr( const TSignal<double> &l, const TSignal<double> &r );
	TSignal<double> autocorr( const TSignal<double> &f );
	TSignal<double>& fftshift( TSignal<double>& f );
	TSignal<double> fft( const TSignal<double> &f );
	TSignal<double> ifft( const TSignal<double> &f );
	TSignal<cx_double> fft_r2c( const TSignal<double> &f );
	TSignal<double> ifft_c2r( const TSignal<cx_double> &f );


	TSignal<cx_double> convolve( const TSignal<cx_double> &l, const TSignal<cx_double> &r );
	TSignal<cx_double> crosscorr( const TSignal<cx_double> &l, const TSignal<cx_double> &r );
	TSignal<cx_double> autocorr( const TSignal<cx_double> &f );
	TSignal<cx_double>& fftshift( TSignal<cx_double>& f );
	TSignal<cx_double> fft( const TSignal<cx_double> &f );
	TSignal<cx_double> ifft( const TSignal<cx_double> &f );

}

#endif /* SIGNAL_FFT_H */
