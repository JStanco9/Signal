//created by John Stanco 7.27.18

#include "signal.h"

#include <iostream>

//TODO - Make all print functions use streams instead of FILE*/printf ( Make code entirely C++ vs C )

namespace Signal{

	/*
	void swap( TSignal<double>& a, TSignal<double>& b ) {
		std::swap(a.dat, b.dat);
		std::swap(a.len, b.len);
		std::swap(a.cap, b.cap);
	}
	*/


	TSignal<double>::TSignal() : dat{ new double[1]() }, len{ 0 }, cap{ 0 } {}


	TSignal<double>::TSignal( const size_t n ) : dat{ new double[n + 1] }, len{ n }, cap{ n } {}


	TSignal<double>::TSignal( const double arr[], const size_t n ) :
	dat{ new double[n + 1]() }, len{ n }, cap{ n } {
		memcpy( dat, arr, ( cap + 1 )*sizeof( double ) );
	}


	TSignal<double>::TSignal( const TSignal& other ) :
	dat{ new double[other.cap + 1]() }, len{ other.len }, cap{ other.cap } {
		memcpy( dat, other.dat, ( cap + 1 )*sizeof( double ) );
	}


	TSignal<double>& TSignal<double>::operator=( TSignal other ) { swap( *this, other ); return *this; }


	TSignal<double>::~TSignal() { delete[] dat; }


	size_t TSignal<double>::length() const { return len; }


	double* TSignal<double>::data() { return dat; }


	double* TSignal<double>::data() const { return dat; }


	double& TSignal<double>::operator()( const size_t index ) {
		if( index >= len ) {
			throw TSignalOutOfRange();
		}
		return dat[index];
	}


	double TSignal<double>::operator()( const size_t index ) const {
		if( index >= len ) {
			throw TSignalOutOfRange();
		}
		return dat[index];
	}


	bool TSignal<double>::operator==( const TSignal<double> &other ) const {
		if( len != other.len ) { return false; }
		for( size_t i = 0; i < len; i++ ) {
			if( dat[ i ] != other( i ) ) { return false; }
		}
		return true;
	}


	bool TSignal<double>::operator!=( const TSignal<double> &other ) const {
		if( len != other.len ) { return true; }
		for( size_t i = 0; i < len; i++ ) {
			if( dat[ i ] != other( i ) ) { return true; }
		}
		return false;
	}


	TSignal<double> TSignal<double>::operator*( const TSignal<double> &other ) const {
		if( len != other.len ) { throw "TSignals must be of same size"; }
		TSignal tmp{ len };
		for( size_t i = 0; i < len; ++i ) {
			tmp.dat[i] = dat[i]*other.dat[i];
		}
		return tmp;
	}


	TSignal<double> TSignal<double>::operator/( const double x ) const {
		TSignal<double> s{ *this };
		for( size_t i = 0; i < len; ++i ) { s.dat[i] /= x; }
		return s;
	}


	TSignal<double>::iterator::iterator( double *elem ) : _elem{ elem } {}


	TSignal<double>::iterator& TSignal<double>::iterator::operator++() {
		++_elem;
		return *this;
	}


	double& TSignal<double>::iterator::operator*() {
		return *_elem;
	}


	bool TSignal<double>::iterator::operator==( const iterator& rhs ) {
		return _elem = rhs._elem;
	}


	bool TSignal<double>::iterator::operator!=( const iterator& rhs ) {
		return _elem != rhs._elem;
	}


	TSignal<double>::const_iterator::const_iterator( double *elem ) : _elem{ elem } {}


	TSignal<double>::const_iterator& TSignal<double>::const_iterator::operator++() {
		++_elem;
		return *this;
	}


	double& TSignal<double>::const_iterator::operator*() {
		return *_elem;
	}


	bool TSignal<double>::const_iterator::operator==( const const_iterator& rhs ) {
		return _elem = rhs._elem;
	}


	bool TSignal<double>::const_iterator::operator!=( const const_iterator& rhs ) {
		return _elem != rhs._elem;
	}


	TSignal<double>::iterator TSignal<double>::begin() {
		return iterator{ &dat[0] };
	}


	TSignal<double>::const_iterator TSignal<double>::begin() const {
		return const_iterator{ &dat[0] };
	}


	TSignal<double>::iterator TSignal<double>::back() {
		return iterator{ &dat[len - 1] };
	}


	TSignal<double>::const_iterator TSignal<double>::back() const {
		return const_iterator{ &dat[len - 1] };
	}


	TSignal<double>::iterator TSignal<double>::end() {
		//data array always contains 1 extra element
		return iterator{ &dat[len] };
	}


	TSignal<double>::const_iterator TSignal<double>::end() const {
		//data array always contains 1 extra element
		return const_iterator{ &dat[len] };
	}


	/////////////////////// COMPLEX Signal ////////////////////////

	/*
	void swap( TSignal<cx_double>& a, TSignal<cx_double>& b ) {
		std::swap(a.dat, b.dat);
		std::swap(a.len, b.len);
		std::swap(a.cap, b.cap);
	}
	*/

	TSignal<cx_double>::TSignal() : dat{ new cx_double[1]() }, len{ 0 }, cap{ 0 } {}


	TSignal<cx_double>::TSignal( const size_t n ) :
	dat{ new cx_double[n] }, len{ n }, cap{ n } {}


	TSignal<cx_double>::TSignal( const cx_double arr[], const size_t n ) :
		dat{ new cx_double[n + 1] }, len{ n }, cap{ n } {
		memcpy( dat, arr, ( cap + 1 )*sizeof( cx_double ) );
	}


	TSignal<cx_double>::TSignal( const TSignal &other ) :
	dat{ new cx_double[other.cap + 1] }, len{ other.len }, cap{ other.cap } {
		memcpy( dat, other.dat, ( cap + 1 )*sizeof( cx_double ) );
	}


	TSignal<cx_double>& TSignal<cx_double>::operator=( TSignal other ) {
		swap( *this, other );
		return *this;
	}


	TSignal<cx_double>::~TSignal() { delete[] dat; }


	size_t TSignal<cx_double>::length() const { return len; }


	cx_double* TSignal<cx_double>::data() { return dat; }


	cx_double* TSignal<cx_double>::data() const { return dat; }


	cx_double& TSignal<cx_double>::operator()( const size_t index ) {
	 if(index>=len){ throw "Index out of bounds."; }
	 return dat[index];
	}


	cx_double TSignal<cx_double>::operator()( const size_t index ) const {
	 if( index>=len ) { throw "Index out of bounds."; }
	 return dat[index];
	}


	bool TSignal<cx_double>::operator==( const TSignal<cx_double> &other ) const {
		if( len != other.len ) { return false; }
		for( size_t i = 0; i < len; i++ ) {
			if( dat[ i ] != other( i ) ) { return false; }
		}
		return true;
	}


	bool TSignal<cx_double>::operator!=( const TSignal<cx_double> &other ) const {
		if( len != other.len ) { return true; }
		for( size_t i = 0; i < len; i++ ) {
			if( dat[i] != other( i ) ) { return true; }
		}
		return false;
	}


	TSignal<cx_double> TSignal<cx_double>::operator*( const TSignal<cx_double> &other ) {
		if( len != other.len ){ throw "TSignals must be of same size"; }
		TSignal<cx_double> tmp( len );
		for( size_t i = 0; i<len; ++i ) { tmp.dat[i] = dat[i]*other.dat[i]; }
		return tmp;
	}

	TSignal<cx_double>::iterator::iterator( cx_double *elem ) : _elem{ elem } {}


	TSignal<cx_double>::iterator& TSignal<cx_double>::iterator::operator++() {
	  ++_elem;
	  return *this;
	}


	cx_double& TSignal<cx_double>::iterator::operator*() {
	  return *_elem;
	}


	bool TSignal<cx_double>::iterator::operator==( const iterator& rhs ) {
	  return _elem = rhs._elem;
	}


	bool TSignal<cx_double>::iterator::operator!=( const iterator& rhs ) {
	  return _elem != rhs._elem;
	}


	TSignal<cx_double>::const_iterator::const_iterator( cx_double *elem ) : _elem{ elem } {}


	TSignal<cx_double>::const_iterator& TSignal<cx_double>::const_iterator::operator++() {
	  ++_elem;
	  return *this;
	}


	cx_double& TSignal<cx_double>::const_iterator::operator*() {
	  return *_elem;
	}


	bool TSignal<cx_double>::const_iterator::operator==( const const_iterator& rhs ) {
	  return _elem = rhs._elem;
	}


	bool TSignal<cx_double>::const_iterator::operator!=( const const_iterator& rhs ) {
	  return _elem != rhs._elem;
	}


	TSignal<cx_double>::iterator TSignal<cx_double>::begin() {
	  return iterator{ &dat[0] };
	}


	TSignal<cx_double>::const_iterator TSignal<cx_double>::begin() const {
	  return const_iterator{ &dat[0] };
	}


	TSignal<cx_double>::iterator TSignal<cx_double>::back() {
	  return iterator{ &dat[len - 1] };
	}


	TSignal<cx_double>::const_iterator TSignal<cx_double>::back() const {
	  return const_iterator{ &dat[len - 1] };
	}


	TSignal<cx_double>::iterator TSignal<cx_double>::end() {
	  //data array always contains 1 extra element
	  return iterator{ &dat[len] };
	}


	TSignal<cx_double>::const_iterator TSignal<cx_double>::end() const {
	  //data array always contains 1 extra element
	  return const_iterator{ &dat[len] };
	}

}

void Signal::print(const Signal::TSignal<double> &f){
	for( size_t i = 0; i < f.length(); ++i ){
		std::cout << f( i ) << "\n";
	}
}


void Signal::print(const Signal::TSignal<double> &f, const char* filename){
	FILE *pFile = fopen(filename, "w");
	if(pFile){ for( size_t i = 0; i < f.length(); ++i ) {
			fprintf(pFile,"%lf\n", f(i)); }
	}
}

void Signal::print(const Signal::TSignal<double> &f, std::string const& filename){
	FILE *pFile = fopen(filename.c_str(), "w");
	if(pFile){ for(size_t i=0;i<f.length();i++){
			fprintf(pFile,"%lf\n", f(i)); }
	}
}


void Signal::print(const Signal::TSignal<cx_double> &f){
	for(size_t i=0;i<f.length();i++){
		std::cout << "re: " <<  f( i ).real() << "\t" << "im: " << f( i ).imag() << "\n";
	}
}

void Signal::print_real(const Signal::TSignal<cx_double> &f, const char* filename){
	FILE *pFile = fopen(filename, "w");
	if(pFile){ for(size_t i=0;i<f.length();i++){
			fprintf(pFile,"%lf\n", f(i).real()); }
	}
}

void Signal::print_imag(const Signal::TSignal<cx_double> &f, const char* filename){
	FILE *pFile = fopen(filename, "w");
	if(pFile){ for(size_t i=0;i<f.length();i++){
			fprintf(pFile,"%lf\n", f(i).imag()); }
	}
}