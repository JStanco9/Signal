//John Stanco 9.17.18


#include "signal.h"

/*

    Interface for using map/apply/reduce with TSignal<T> class
    One could choose to make map/apply/reduce member functions
    of the TSignal<T> class.  However, in order to provide some
    degree of uniformity with using std::vector<T> as another
    potential data type, I instead define functions that take in
    a TSignal<T> as a parameter.

    This allows one to provide an identical interface w.r.t
    std::vector. ( probably down the line if one sees a need )

    Main Idea:

    One could generalize maps from {U} -> {V} as functors
    of type <U, V>.  This would allow one to define a map operation
    that maps a TSignal<U> and a functor F from {U} -> {V} and return a
    TSignal<V>, which is the image of F in {V}.

*/

#ifndef SIGNAL_FUNCTOR_H
#define SIGNAL_FUNCTOR_H


template<class U, class V>
class UnaryFunctor {
public:
  virtual V operator()( const U& ) = 0;
};


template<class U, class V>
class BinaryFunctor {
public:
  virtual V operator()( const U&, const U& ) = 0;
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
