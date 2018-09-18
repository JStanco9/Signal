//John Stanco 9.12.18

#include <iostream>
#include <random>

std::random_device rseed;
std::mt19937_64 generator(rseed());



template<class T>
T random() {
	T e = ( T )generator();
	return e;
}



int main() {

	std::cout << random<int>() << std::endl;
	std::cout << random<size_t>() << std::endl;
	std::cout << random<float>() << std::endl;
	std::cout << random<double>() << std::endl;


	std::cout << sizeof(int) << std::endl;
	std::cout << sizeof(size_t) << std::endl; 
	std::cout << sizeof(float) << std::endl;
	std::cout << sizeof(double) << std::endl;

	return 0;
}