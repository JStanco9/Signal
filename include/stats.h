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
	template<typename Container> double mean( const Container &);
	template<typename Container> double var( const Container& );
	template<typename Container> double jackknife_var( const Container&, const size_t );
	template<typename Container> double jackknife_mean( const Container&, const size_t );
	template<typename Container> double bootstrap_var( const Container&, const size_t );
	template<typename Container> double bootstrap_mean( const Container&, const size_t );
	template<typename Container> double blocked_var( const Container&, const size_t );
	template<typename Container> double blocked_mean( const Container&, const size_t );
	template<typename Container> Container blocked_means( const Container&, const size_t );
	template<typename Container> Container deleted_averages( const Container& );
	template<typename Container> double skewness( const Container& );
	template<typename Container> double kurtosis( const Container& );
	template<typename Container> void normalize_std_inpl( Container& );
	template<typename Container> Container normalize_std( const Container& );
}

#endif /* SIGNAL_STATS_H */
