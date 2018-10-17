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
namespace Signal{
	double mean( const TSignal<double> &);
	double var( const TSignal<double> &);
	double jackknife_var( const TSignal<double> &, const size_t );
	double jackknife_mean( const TSignal<double> &, const size_t );
	double bootstrap_var( const TSignal<double> &, const size_t );
	double bootstrap_mean( const TSignal<double> &, const size_t );
	double blocked_var( const TSignal<double> &, const size_t );
	double blocked_mean( const TSignal<double> &, const size_t );
	TSignal<double> blocked_means( const TSignal<double> &, const size_t );
	TSignal<double> deleted_averages( const TSignal<double> & );
	double skewness( const TSignal<double> & );
	double kurtosis( const TSignal<double> & );
	void normalize_std_inpl( TSignal<double>& );
	TSignal<double> normalize_std( const TSignal<double>& );
}

#endif /* SIGNAL_STATS_H */