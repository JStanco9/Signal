//John Stanco 9.18.18

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

}

#endif /* SIGNAL_TSIGNAL_H */
