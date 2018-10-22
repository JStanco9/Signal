//John Stanco 9.18.18


#include <complex>
#include <memory>


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

	template<typename Container, typename value_type, typename Stream>
	void print( const Container<value_type> &s, Stream& ofs ){
		for( auto& x : s ) {
			ofs << x << "\n";
		}
	}


	template<typename Container, typename Stream>
	void print( const Container<cx_double, A> &s, Stream& ofs ) {
		for( auto& z : s ) {
			ofs << std::left << std::setw(24) << std::setprecision(16) << z.real() << std::left << z.imag() << "\n";
		}
	}


	template<typename Container>
	void print( const Container<cx_double, A> &s, std::basic_ostream<char> &ofs = stds::cout ) {
		for( auto& z : s ) {
			ofs << std::left << std::setw(24) << std::setprecision(16) << z.real() << std::left << z.imag() << "\n";
		}
	}


	template<typename Container, typename value_type>
	void print( const Container<value_type> &s, std::basic_ostream<char> ofs = std::cout ) {
		for( auto& x : s ) {
			ofs << x << "\n";
		}
	}


	template<typename Container, typename Stream>
	void printReal( const Container<cx_double, A> &s, Stream &os ) {
		if( !os.is_open() ) {
			throw "Attempted to print to unopen file stream";
		}
		for( auto& z : s ) {
			os << std::left << std::setprecision(16) << z.real() << "\n";
		}
	}


	template<typename Container>
	void printReal( const Container<cx_double, A> &s, std::basic_ostream<char> &os = std::cout ) {
		if( !os.is_open() ) {
			throw "Attempted to print to unopen file stream";
		}
		for( auto& z : s ) {
			os << std::left << std::setprecision(16) << z.real() << "\n";
		}
	}


	template<typename Container, typename Stream>
	void printImag( const Container<cx_double, A> &s, Stream &os ) {
		if( !os.is_open() ) {
			throw "Attempted to print to unopen file stream";
		}
		for( auto& z : s ) {
			os << std::left << std::setprecision(16) << z.imag() << "\n";
		}
	}


	template<typename Container>
	void printImag( const Container<cx_double, A> &s, std::basic_ostream<char> &os = std::cout ) {
		if( !os.is_open() ) {
			throw "Attempted to print to unopen file stream";
		}
		for( auto& z : s ) {
			os << std::left << std::setprecision(16) << z.imag() << "\n";
		}
	}

#endif /* SIGNAL_TSIGNAL_H */
