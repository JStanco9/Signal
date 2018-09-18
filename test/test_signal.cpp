//John Stanco 9.10.18

#include "/usr/local/include/signal.h"


#include <iostream>


// TODO : write unit testing framework ( https://accu.org/index.php/journals/368 )
//				Develop into its own library -> that one can include in all test makefiles

//				Use this testing framework to write tests for all src file features separately



size_t n_tests = 0;
size_t n_passed = 0;
size_t n_failed = 0;


void runTest( bool ( *test )( void ) ) {
	test()? n_passed++ : n_failed++;
	n_tests++;
}


bool testSignalDoubleEqualsEmpty() {
	Signal::TSignal<double> s{};
	return s == s;
}


bool testSignalDoubleEqualsFilled() {
	Signal::TSignal<double> s{ 5 };
	size_t i = 0;
	for( auto& x : s ) {
		x = i++;
	}
	return s == s;
}


bool testSignalDoubleCopyEmpty() {
	Signal::TSignal<double> s{};
	Signal::TSignal<double>test{ s };
	return s == test;
}


bool testSignalDoubleCopyFilled() {
	Signal::TSignal<double> s{ 5 };
	size_t i = 0;
	for( auto& x : s ) {
		x = i++;
	}
	Signal::TSignal<double> test{ s };
	return s == test;
}


bool testSignalDoubleAssignEmpty() {
	Signal::TSignal<double> s{};
	Signal::TSignal<double> test{ 5 };
	size_t i = 0;
	for( auto& x : s ) {
		x = i++;
	}
	test = s;
	return s == test;
}


bool testSignalDoubleAssignFilled() {
	Signal::TSignal<double> s{ 5 };
	Signal::TSignal<double> test{};
	size_t i = 0;
	for( auto& x : s ) {
		x = i++;
	}
	test = s;
	return s == test;
}


bool testSignalDoubleArrayCtr() {
	double arr[6] = { 1., 2., 3., 4., 5., 6. };
	Signal::TSignal<double> s{ arr, 6 };
	Signal::TSignal<double> test{ 6 };
	size_t i = 1;
	for( auto& x : test ) {
		x = i++;
	}
	return s == test;
}


void testSignalDouble() {
	runTest( &testSignalDoubleEqualsEmpty );
	runTest( &testSignalDoubleEqualsFilled );
	runTest( &testSignalDoubleCopyEmpty );
	runTest( &testSignalDoubleCopyFilled );
	runTest( &testSignalDoubleAssignEmpty );
	runTest( &testSignalDoubleAssignFilled );
	runTest( &testSignalDoubleArrayCtr );
}


bool testSignalComplexDoubleEqualsEmpty() {
	Signal::TSignal<cx_double> s{};
	return s == s;
}


bool testSignalComplexDoubleEqualsFilled() {
	Signal::TSignal<cx_double> s{ 5 };
	size_t i = 0;
	for( auto& x : s ) {
		x = cx_double( i, 3 - i ); i++;
	}
	return s == s;
}


bool testSignalComplexDoubleCopyEmpty() {
	Signal::TSignal<cx_double> s{};
	Signal::TSignal<cx_double>test{ s };
	return s == test;
}


bool testSignalComplexDoubleCopyFilled() {
	Signal::TSignal<cx_double> s{ 5 };
	size_t i = 0;
	for( auto& x : s ) {
		x = cx_double( i, 3 - i ); i++;
	}
	Signal::TSignal<cx_double> test{ s };
	return s == test;
}


bool testSignalComplexDoubleAssignEmpty() {
	Signal::TSignal<cx_double> s{};
	Signal::TSignal<cx_double> test{ 5 };
	size_t i = 0;
	for( auto& x : s ) {
		x = cx_double( i, 3 - i ); i++;
	}
	test = s;
	return s == test;
}


bool testSignalComplexDoubleAssignFilled() {
	Signal::TSignal<cx_double> s{ 5 };
	Signal::TSignal<cx_double> test{};
	size_t i = 0;
	for( auto& x : s ) {
		x = cx_double( i, 3 - i ); i++;
	}
	test = s;
	return s == test;
}


bool testSignalComplexDoubleArrayCtr() {
	cx_double arr[6];
	Signal::TSignal<cx_double> test{ 6 };
	for( size_t i = 0; i < 6; ++i ) {
		arr[i] = cx_double( i + 1, 2 - i );
		test( i ) = cx_double( i + 1, 2 - i );
	}
	Signal::TSignal<cx_double> s{ arr, 6 };
	return s == test;
}


void testSignalComplexDouble() {
	runTest( &testSignalComplexDoubleEqualsEmpty );
	runTest( &testSignalComplexDoubleEqualsFilled );
	runTest( &testSignalComplexDoubleCopyEmpty );
	runTest( &testSignalComplexDoubleCopyFilled );
	runTest( &testSignalComplexDoubleAssignEmpty );
	runTest( &testSignalComplexDoubleAssignFilled );
	runTest( &testSignalComplexDoubleArrayCtr );
}


bool testTSignalIterator() {
	Signal::TSignal<int> s{ 5 };
	Signal::TSignal<int> test{ 5 };
	for( size_t i = 0; i < 5; ++i ) {
		test( i ) = 4;
	}
	for( auto& x : s ) {
		x = 4;
	}
	return s == test;
}


bool testTSignalPush_back() {
	Signal::TSignal<int> s{};
	Signal::TSignal<int> test{ 5 };

	for( size_t i = 0; i < 5; ++i ) {
		test( i ) = i;
		s.push_back( i );
	}

	return s == test;
}


bool testTSignalCopyCtr() {
	Signal::TSignal<int> s{};

	for( size_t i = 0; i < 5; ++i ) {
		s.push_back( i );
	}

	Signal::TSignal<int> test{ s };
	return s == test;
}


bool testTSignalAssignment() {
	Signal::TSignal<int> s{};
	Signal::TSignal<int> test{};

	for( size_t i = 0; i < 5; ++i ) {
		s.push_back( i );
	}

	test = s;
	return s == test;
}

int add2( const int &x ) {
	return x + 2;
}

bool testTSignalFunctorMap() {
	UnaryFunctionWrap<int, int> f( add2 );
	Signal::TSignal<int> s{ 6 };
	Signal::TSignal<int> test{ 6 };

	for( size_t i = 0; i < 6; ++i ) {
		s( i ) = i;
		test( i ) = i + 2;
	}

	auto image = map( s, f );
	return image == test;
}


void testTSignal() {
	runTest( &testTSignalIterator );
	runTest( &testTSignalPush_back );
	runTest( &testTSignalCopyCtr );
	runTest( &testTSignalAssignment );
	runTest( &testTSignalFunctorMap );
}


void testSignal() {

	testSignalDouble();
	testSignalComplexDouble();
	testTSignal();

	std::cout << "Passed: " << n_passed << "/" << n_tests << std::endl;
	std::cout << "Failed: " << n_failed << "/" << n_tests << std::endl;

}


int main() {

	testSignal();
	return 0;
}
