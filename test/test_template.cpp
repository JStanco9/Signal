//John Stanco 9.12.18

#include <iostream>
#include <cmath>


template<class T, bool isNumber>
class A {
	T obj;
public:
	A() {}
	A( const T& x ) : obj{ x } {}
	void print() { std::cout << obj << "\n"; }
	T& self() { return obj; }
	T self() const { return obj; }
};

template<class T>
class A<T, true> {
	T obj;
public:
	A() {}
	A( const T& x ) : obj{ x } {}
	void print() { std::cout << obj << "\n"; }

	template<class U>
	A operator*( const A<U, true>& other ) { return A{ other.self() * obj }; }
	T& self() { return obj; }
	T self() const  { return obj; }
};

template<class T>
using Number = A<T, true>;

template<class T>
using Object = A<T, false>;


class Dog {
public:
	void bark() { std::cout << "Woof\n"; }
};


class B {
public:
	template<class U, class V>
	V operator()( const U& );
};


template<>
double B::operator()( const double& x ) { return x * x; }


class C {
	static double x;
public:
	C() { x = 2; }
	static double addx( double y ) { return x + y; }
};


struct Constants {
public:
	constexpr static double pi 		= 3.141592653589793;
	constexpr static double e			= 2.718281828459045;
	constexpr static double h 		= 6.62607004081e-34;
	constexpr static double hbar 	= h / ( 2 * pi );
	//constexpr static double kb 		= 
	constexpr static double c 		= 2.99792458e8;
};



int main() {

	Number<double> b{ 4. };
	Number<int> a{ 3 };
	(b * a).print();

	Dog Rover;

	Rover.bark();

	Object<Dog> wrappedDog( Rover );
	wrappedDog.self().bark();



	B f;
	double c = f.operator()<double, double>( 3.0 );

	C cee;

	std::cout << C::addx( 3.0 ) << "\n";


	double r = 14.7;
	std::cout << Constant::pi*r*r << std::endl;

	std::cout << 0.3 - 0.2 << std::endl;


	return 0;
}