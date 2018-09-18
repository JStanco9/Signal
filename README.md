# signal

Template-based container interface that works very similarly to std::vector, with added functionality for statistics and signal processing.

 - tsignal.h/cpp | The container itself, TSignal<T>, which behaves much like std::vector.
 - stats.h/cpp   | Interface for statistical analysis of TSignal<double> class.
 - fft.h/cpp     | Interface for Fourier Transforms of TSignal<double> and TSignal<cx_double> ( implements fftw3 library ).
 - functor.h     | Interface for apply/map/reduce with two base classes UnaryFunctor<U, V> and BinaryFunctor<U, V>.
 - myrand.h/cpp  | Interface for random number generation.
 
 - signal.h      | Contains all header files in one, providing all functionality when included.

To make : make - make install ( installs library libsignal.a to /usr/local/lib and installs header signal.h to /usr/local/include )
