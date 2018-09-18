# signal

Template-based container interface that works very similarly to std::vector, with added functionality of equipping contaners with binary operations that classify a container as existing in a vector space.  This includes multiplication by a scalar, and dot-product operations.  

 - tsignal.h/cpp | The container itself, TSignal<T>, which behaves much like std::vector.
 - stats.h/cpp  | Interface for statistical analysis of TSignal<double> class.
 - fft.h/cpp    | Interface for Fourier Transforms of TSignal<double> and TSignal<cx_double> ( implements fftw3 library ) 
 - myrand.h/cpp | Interface for random number generation

Note - All headers are combined in the signal.h header file, which provides all functionality when included.
