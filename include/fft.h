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
	template<template<typename> class Container, typename T>
	Container<T>& shift( Container<T>& );
};


class FFT : public FFT_base {
	FFT( const FFT &other ) = delete;
	FFT& operator=( FFT other ) = delete;
public:
	template<template<typename, typename> class Container, typename A>
	static Container<cx_double, A> c2c( Container<cx_double, A> &f );
	template<template<typename, typename> class Container, typename A>
	static Container<double, A> r2r( Container<double, A> &f );
	template<template<typename, typename> class Container, typename A>
	static Container<cx_double, A> r2c( Container<double, A> &f );
	template<template<typename, typename> class Container, typename A>
	static void inPlace( Container<cx_double, A>& );
	template<template<typename, typename> class Container, typename A>
	static void inPlace( Container<double, A>& );
};


class IFFT : public FFT_base {
	IFFT( const IFFT &other ) = delete;
	IFFT& operator=( IFFT other ) = delete;
public:
	template<template<typename, typename> class Container, typename A>
	static Container<cx_double, A> c2c( Container<cx_double, A> &f );
	template<template<typename, typename> class Container, typename A>
	static Container<double, A> r2r( Container<double, A> &f );
	template<template<typename, typename> class Container, typename A>
	static Container<double, A> c2r( Container<cx_double, A> &f );
	template<template<typename, typename> class Container, typename A>
	static void inPlace( Container<cx_double, A>& );
	template<template<typename, typename> class Container, typename A>
	static void inPlace( Container<double, A>& );
};


template<template<typename, typename> class Container, typename A>
Container<double, A> convolve( const Container<double, A> &l, const Container<double, A> &r );
template<template<typename, typename> class Container, typename A>
Container<double, A> crosscorr( const Container<double, A> &l, const Container<double, A> &r );
template<template<typename, typename> class Container, typename A>
Container<double, A> autocorr( const Container<double, A> &f );
template<template<typename, typename> class Container, typename A>
Container<double, A>& fftshift( Container<double, A>& f );
template<template<typename, typename> class Container, typename A>
Container<double, A> fft( const Container<double, A> &f );
template<template<typename, typename> class Container, typename A>
Container<double, A> ifft( const Container<double, A> &f );
template<template<typename, typename> class Container, typename A>
Container<cx_double, A> fft_r2c( const Container<double, A> &f );
template<template<typename, typename> class Container, typename A>
Container<double, A> ifft_c2r( const Container<cx_double, A> &f );


template<template<typename, typename> class Container, typename A>
Container<cx_double, A> convolve( const Container<cx_double, A> &l, const Container<cx_double, A> &r );
template<template<typename, typename> class Container, typename A>
Container<cx_double, A> crosscorr( const Container<cx_double, A> &l, const Container<cx_double, A> &r );
template<template<typename, typename> class Container, typename A>
Container<cx_double, A> autocorr( const Container<cx_double, A> &f );
template<template<typename, typename> class Container, typename A>
Container<cx_double, A>& fftshift( Container<cx_double, A>& f );
template<template<typename, typename> class Container, typename A>
Container<cx_double, A> fft( const Container<cx_double, A> &f );
template<template<typename, typename> class Container, typename A>
Container<cx_double, A> ifft( const Container<cx_double, A> &f );

#endif /* SIGNAL_FFT_H */
