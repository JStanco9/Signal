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


#ifndef SIGNAL_H
#define SIGNAL_H
#ifndef SIGNAL_TSIGNAL_H
#define SIGNAL_TSIGNAL_H


typedef std::complex<double> cx_double;
extern double pi;

namespace Signal{

	template<class T> class TSignal;
	template<class T> void swap( TSignal<T> &a, TSignal<T> &b );

	template<class T>
	class TSignal {
	private:
		T *dat;
		size_t len;
		size_t cap;

		friend void swap<T>( TSignal& a, TSignal& b );
		void expand();
	public:
		TSignal();
		TSignal( const size_t n );
		TSignal( const T arr[] );
		TSignal( const TSignal &other );
		TSignal& operator=( TSignal rhs );
		~TSignal();

		size_t length() const;
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


	template<class T>
	void swap( TSignal<T> &a, TSignal<T> &b ) {
		std::swap( a.dat, b.dat );
		std::swap( a.len, b.len );
		std::swap( a.cap, b.cap );
	}


	template<class T>
	TSignal<T>::TSignal() : dat{ new T[1] }, len{ 0 }, cap{ 0 } {}


	template<class T>
	TSignal<T>::TSignal( const size_t n ) : dat{ new T[n + 1]() }, len{ n }, cap{ n } {}


	template <class T>
	TSignal<T>::TSignal( const TSignal<T>& other ) :
	dat{ new T[other.len + 1]() }, len{ other.len }, cap{ other.cap } {
		memcpy( dat, other.dat, (cap + 1) * sizeof( T ));
	}


	template<class T>
	TSignal<T>& TSignal<T>::operator=( TSignal rhs ) {
	 swap( *this, rhs );
	 return *this;
	}


	template <class T>
	TSignal<T>::~TSignal() {
		delete[] dat;
	}


	template<class T>
	size_t TSignal<T>::length() const { return len; }


	template<class T>
	T* TSignal<T>::data() { return dat; }


	template<class T>
	T* TSignal<T>::data() const { return dat; }


	template<class T>
	void TSignal<T>::expand() {
		( cap == 0 )? ++cap : cap *= 2;
		T *tmp = ( T* )realloc( dat, ( cap + 1 ) * sizeof( T ) );
		if( tmp ) { dat = tmp; }
		//should not need to free original data if realloc works.
	}


	template<class T>
	void TSignal<T>::push_back( const T &item ) {
		if( len == cap ) { expand(); }
		dat[len++] = item;
	}


	template<class T>
	T& TSignal<T>::operator()( const size_t index ) {
		if( index >= len ) { throw TSignalOutOfRange(); }
		return dat[index];
	}


	template<class T>
	T TSignal<T>::operator()( const size_t index ) const {
		if( index >= len ) { throw TSignalOutOfRange(); }
		return dat[index];
	}


	template<class T>
	T& TSignal<T>::operator[](const size_t index) {
		return this->operator()( index );
	}


	template<class T>
	T TSignal<T>::operator[](const size_t index) const {
		return this->operator()( index );
	}


	template<class T>
	bool TSignal<T>::operator==( const TSignal &rhs ) const {
		if( len != rhs.len ) { return false; }
		auto itl = begin();
		auto itr = rhs.begin();
		for( ; itl != end() && itr != rhs.end(); ++itl, ++itr ) {
			if( *itl != *itr ) { return false; }
		}
		return true;
	}


	template<class T>
	bool TSignal<T>::operator!=( const TSignal &rhs ) const {
		if( len != rhs.len ) { return true; }
		auto itl = begin();
		auto itr = rhs.begin();
		for( ; itl != end() && itr != rhs.end(); ++itl, ++itr ) {
			if( *itl != *itr ) { return true; }
		}
		return false;
	}


	template<class T>
	TSignal<T>::iterator::iterator( T *elem ) : _elem{ elem } {}


	template<class T>
	typename TSignal<T>::iterator& TSignal<T>::iterator::operator++() {
		++_elem;
		return *this;
	}


	template<class T>
	T& TSignal<T>::iterator::operator*() {
		return *_elem;
	}


	template<class T>
	bool TSignal<T>::iterator::operator==( const iterator& rhs ) {
		return _elem = rhs._elem;
	}


	template<class T>
	bool TSignal<T>::iterator::operator!=( const iterator& rhs ) {
		return _elem != rhs._elem;
	}


	template<class T>
	TSignal<T>::const_iterator::const_iterator( T *elem ) : _elem{ elem } {}


	template<class T>
	typename TSignal<T>::const_iterator& TSignal<T>::const_iterator::operator++() {

		++_elem;
		return *this;
	}


	template<class T>
	T& TSignal<T>::const_iterator::operator*() {
		return *_elem;
	}


	template<class T>
	bool TSignal<T>::const_iterator::operator==( const const_iterator& rhs ) {
		return _elem = rhs._elem;
	}


	template<class T>
	bool TSignal<T>::const_iterator::operator!=( const const_iterator& rhs ) {
		return _elem != rhs._elem;
	}


	template<class T>
	typename TSignal<T>::iterator TSignal<T>::begin() {
		return iterator{ &dat[0] };
	}


	template<class T>
	typename TSignal<T>::const_iterator TSignal<T>::begin() const {
		return const_iterator{ &dat[0] };
	}


	template<class T>
	typename TSignal<T>::iterator TSignal<T>::back() {
		return iterator{ &dat[len - 1] };
	}


	template<class T>
	typename TSignal<T>::const_iterator TSignal<T>::back() const {
		return const_iterator{ &dat[len - 1] };
	}


	template<class T>
	typename TSignal<T>::iterator TSignal<T>::end() {
		/// dat array always contains 1 extra element
		return iterator{ &dat[len] };
	}


	template<class T>
	typename TSignal<T>::const_iterator TSignal<T>::end() const {
		/// dat array always contains 1 extra element
		return const_iterator{ &dat[len] };
	}


	template<>
	class TSignal<double> {
	private:
		double *dat;
		size_t len;
		size_t cap;

		friend void swap<double>( TSignal& a, TSignal& b );
	public:
		TSignal();
		TSignal( const size_t n );
		TSignal( const TSignal &other );
		TSignal( const double arr[], const size_t n );
		TSignal& operator=( TSignal other );
		~TSignal();

		size_t length() const;
		double* data();
		double* data() const;

		double& operator()( const size_t index );
		double operator()( const size_t index ) const;
		double& operator[]( const size_t index );
		double operator[]( const size_t index ) const;
		bool operator==( const TSignal &other ) const;
		bool operator!=( const TSignal &other ) const;
		TSignal operator*( const TSignal &other ) const;
		TSignal operator/( const double x ) const;

		class iterator {
			friend class TSignal;
			double *_elem;
		public:
			iterator( double *elem );
			iterator& operator++();
			//iterator& operator++(int);
			double& operator*();
			bool operator==( const iterator &rhs );
			bool operator!=( const iterator &rhs );
		};

		class const_iterator {
			friend class TSignal;
			double *_elem;
		public:
			const_iterator( double *elem );
			const_iterator& operator++();
			//const_iterator& operator++(int);
			double& operator*();
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


	template<>
	class TSignal<cx_double> {
	private:
		cx_double *dat;
		size_t len;
		size_t cap;

		friend void swap<cx_double>( TSignal& a, TSignal& b );
	public:
		TSignal();
		TSignal( const size_t n );
		TSignal( const cx_double arr[], const size_t n );
		TSignal& operator=( TSignal other );
		TSignal( const TSignal &other );
		~TSignal();

		size_t length() const;
		cx_double* data();
		cx_double* data() const;

		cx_double& operator()( const size_t index );
		cx_double operator()( const size_t index ) const;
		cx_double& operator[]( const size_t index );
		cx_double operator[]( const size_t index ) const;
		bool operator==( const TSignal &other ) const;
		bool operator!=( const TSignal &other ) const;
		TSignal operator*( const TSignal &other );

		class iterator {
			friend class TSignal;
			cx_double *_elem;
		public:
			iterator( cx_double *elem );
			iterator& operator++();
			//iterator& operator++(int);
			cx_double& operator*();
			bool operator==( const iterator &rhs );
			bool operator!=( const iterator &rhs );
		};

		class const_iterator {
			friend class TSignal;
			cx_double *_elem;
		public:
			const_iterator( cx_double *elem );
			const_iterator& operator++();
			//const_iterator& operator++(int);
			cx_double& operator*();
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

	void print( const TSignal<cx_double> & );
	void print_real( const TSignal<cx_double> &, const char* );
	void print_imag( const TSignal<cx_double> &, const char* );
	void print_real( const TSignal<cx_double> &, const std::string& );
	void print_imag( const TSignal<cx_double> &, const std::string& );
	void print( const TSignal<double> & );
	void print( const TSignal<double> &, const char* );
	void print( const TSignal<double> &, const std::string& );

	typedef TSignal<cx_double> cx_signal;
	typedef TSignal<double> signal;

}

#endif /* SIGNAL_TSIGNAL_H */
#ifndef SIGNAL_STATS_H
#define SIGNAL_STATS_H


/// Implements blocking, jacknife, bootstrapping analysis of real-valued random variables
namespace Signal {
	double mean( const TSignal<double> & );
	double var( const TSignal<double> & );
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
#ifndef MYRAND_H
#define MYRAND_H

int flip( const double x = .5 );
int randInt( const int, const int );
int randInt( const int );
double min( const double, const double );
double randDouble( const double, const double );
double plusMinus( const double );
double randNorm();
double randNorm( const double, const double );

#endif /* MYRAND_H */

#ifndef SIGNAL_FFT_H
#define SIGNAL_FFT_H

namespace Signal {

	class FFT_base {
		FFT_base( const FFT_base& other ) = delete;
		FFT_base& operator=( FFT_base other ) = delete;
	public:
		template<class T>
		Signal::TSignal<T>& shift( Signal::TSignal<T>& );
	};


	class FFT : public FFT_base {
	public:
		static Signal::TSignal<cx_double> c2c( Signal::TSignal<cx_double> &f );
		static Signal::TSignal<double> r2r( Signal::TSignal<double> &f );
		static Signal::TSignal<cx_double> r2c( Signal::TSignal<double> &f );
		static void inPlace( Signal::TSignal<cx_double>& );
		static void inPlace( Signal::TSignal<double>& );
	};


	class IFFT : public FFT_base {
		IFFT( const IFFT &other ) = delete;
		IFFT& operator=( IFFT other ) = delete;
	public:
		static Signal::TSignal<cx_double> c2c( Signal::TSignal<cx_double> &f );
		static Signal::TSignal<double> r2r( Signal::TSignal<double> &f );
		static Signal::TSignal<double> c2r( Signal::TSignal<cx_double> &f );
		static void inPlace( Signal::TSignal<cx_double>& );
		static void inPlace( Signal::TSignal<double>& );
	};

	TSignal<double> convolve( const TSignal<double> &l, const TSignal<double> &r );
	TSignal<double> crosscorr( const TSignal<double> &l, const TSignal<double> &r );
	TSignal<double> autocorr( const TSignal<double> &f );
	TSignal<double> fft( const TSignal<double> &f );
	TSignal<double> ifft( const TSignal<double> &f );
	TSignal<cx_double> fft_r2c( const TSignal<double> &f );
	TSignal<double> ifft_c2r( const TSignal<cx_double> &f );


	TSignal<cx_double> convolve( const TSignal<cx_double> &l, const TSignal<cx_double> &r );
	TSignal<cx_double> crosscorr( const TSignal<cx_double> &l, const TSignal<cx_double> &r );
	TSignal<cx_double> autocorr( const TSignal<cx_double> &f );
	TSignal<cx_double> fft( const TSignal<cx_double> &f );
	TSignal<cx_double> ifft( const TSignal<cx_double> &f );

	template<class T> TSignal<T>& fftshift( TSignal<T> &f );

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
