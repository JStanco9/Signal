//John Stanco 9.10.18

#include <tsignal.h>
#include <iostream>
#include <unittest.h>
#include <vector>


void testSignalDoubleEqualsEmpty() {
	Signal::TSignal<double> s{};
	ASSERT_TRUE( s == s );;
}


void testSignalDoubleEqualsFilled() {
	Signal::TSignal<double> s{ 5 };
	size_t i = 0;
	for( auto& x : s ) {
		x = i++;
	}
	ASSERT_TRUE( s == s );;
}


void testSignalDoubleCopyEmpty() {
	Signal::TSignal<double> s{};
	Signal::TSignal<double>test{ s };
	ASSERT_TRUE( s == test );
}


void testSignalDoubleCopyFilled() {
	Signal::TSignal<double> s{ 5 };
	size_t i = 0;
	for( auto& x : s ) {
		x = i++;
	}
	Signal::TSignal<double> test{ s };
	ASSERT_TRUE( s == test );
}


void testSignalDoubleAssignEmpty() {
	Signal::TSignal<double> s{};
	Signal::TSignal<double> test{ 5 };
	size_t i = 0;
	for( auto& x : s ) {
		x = i++;
	}
	test = s;
	ASSERT_TRUE( s == test );
}


void testSignalDoubleAssignFilled() {
	Signal::TSignal<double> s{ 5 };
	Signal::TSignal<double> test{};
	size_t i = 0;
	for( auto& x : s ) {
		x = i++;
	}
	test = s;
	ASSERT_TRUE( s == test );
}


void testSignalDouble() {
	testSignalDoubleEqualsEmpty();
	testSignalDoubleEqualsFilled();
	testSignalDoubleCopyEmpty();
	testSignalDoubleCopyFilled();
	testSignalDoubleAssignEmpty();
	testSignalDoubleAssignFilled();
}


void testSignalComplexDoubleEqualsEmpty() {
	Signal::TSignal<cx_double> s{};
	ASSERT_TRUE( s == s );;
}


void testSignalComplexDoubleEqualsFilled() {
	Signal::TSignal<cx_double> s{ 5 };
	size_t i = 0;
	for( auto& x : s ) {
		x = cx_double( i, 3 - i ); i++;
	}
	ASSERT_TRUE( s == s );
}


void testSignalComplexDoubleCopyEmpty() {
	Signal::TSignal<cx_double> s{};
	Signal::TSignal<cx_double>test{ s };
	ASSERT_TRUE( s == test );
}


void testSignalComplexDoubleCopyFilled() {
	Signal::TSignal<cx_double> s{ 5 };
	size_t i = 0;
	for( auto& x : s ) {
		x = cx_double( i, 3 - i ); i++;
	}
	Signal::TSignal<cx_double> test{ s };
	ASSERT_TRUE( s == test );
}


void testSignalComplexDoubleAssignEmpty() {
	Signal::TSignal<cx_double> s{};
	Signal::TSignal<cx_double> test{ 5 };
	size_t i = 0;
	for( auto& x : s ) {
		x = cx_double( i, 3 - i ); i++;
	}
	test = s;
	ASSERT_TRUE( s == test );
}


void testSignalComplexDoubleAssignFilled() {
	Signal::TSignal<cx_double> s{ 5 };
	Signal::TSignal<cx_double> test{};
	size_t i = 0;
	for( auto& x : s ) {
		x = cx_double( i, 3 - i ); i++;
	}
	test = s;
	ASSERT_TRUE( s == test );
}


void testSignalComplexDouble() {
	testSignalComplexDoubleEqualsEmpty();
	testSignalComplexDoubleEqualsFilled();
	testSignalComplexDoubleCopyEmpty();
	testSignalComplexDoubleCopyFilled();
	testSignalComplexDoubleAssignEmpty();
	testSignalComplexDoubleAssignFilled();
}


void testTSignalIterator() {
	Signal::TSignal<int> s{ 5 };
	Signal::TSignal<int> test{ 5 };
	for( size_t i = 0; i < 5; ++i ) {
		test( i ) = 4;
	}
	for( auto& x : s ) {
		x = 4;
	}
	ASSERT_TRUE( s == test );
}


void testTSignalPush_back() {
	Signal::TSignal<int> s{};
	Signal::TSignal<int> test{ 5 };

	for( size_t i = 0; i < 5; ++i ) {
		test( i ) = i;
		s.push_back( i );
	}
	ASSERT_TRUE( s == test );
}


void testTSignalCopyCtr() {
	Signal::TSignal<int> s{};

	for( size_t i = 0; i < 5; ++i ) {
		s.push_back( i );
	}

	Signal::TSignal<int> test{ s };
	ASSERT_TRUE( s == test );
}


void testTSignalAssignment() {
	Signal::TSignal<int> s{};
	Signal::TSignal<int> test{};

	for( size_t i = 0; i < 5; ++i ) {
		s.push_back( i );
	}

	test = s;
	ASSERT_TRUE( s == test );
}

int add2( const int &x ) {
	return x + 2;
}

void testTSignalFunctorMap() {
	UnaryFunctionWrap<int, int> f( add2 );
	Signal::TSignal<int> s{ 6 };
	Signal::TSignal<int> test{ 6 };

	for( size_t i = 0; i < 6; ++i ) {
		s( i ) = i;
		test( i ) = i + 2;
	}

	auto image = map( s, f );
	ASSERT_TRUE( image == test );
}


void testTSignal() {
	testTSignalIterator();
	testTSignalPush_back();
	testTSignalCopyCtr();
	testTSignalAssignment();
	testTSignalFunctorMap();
}


void testFFT() {
	size_t N = 1 << 4;
	Signal::TSignal<cx_double> f( 2 * N );
	for( size_t i = 0; i < f.size(); ++i ) {
		f[i] = { cos( i * pi / N ), sin( i * pi / N ) };
	}
	print( f );
	std::cout << "\n";
	auto F = fft( f );
	print( fftshift( F ) );

}


void testSignal() {
	testSignalDouble();
	testSignalComplexDouble();
	testTSignal();
}


int main() {
	try{
		testSignal();
	} catch( const char* e ) {
		std::cout << e << std::endl;
		exit( EXIT_FAILURE );
	}

	return 0;
}
