CXX= g++
CPPFLAGS= -Wall -lgsl -O3 -lgslcblas
LDFLAGS= -Wall -lgsl -O3 -lgslcblas
LDFLAGS_CYGWIN= -Wall -O3 -Wl,--stack,1004194304 -lgsl -lgslcblas

OBJS= fitujeme.o read_data.o misc.o write_data.o reddening.o filter_zeropoints.o lbol.o deriv.o mni.o bolometric.o priors.o ln85.o mni_dependence.o physical_parameters.o
FIT_ONE= fit_single.o
FIT_EXPTIME= fit_exptime.o

clean:
	rm -f *.o


fit_single: $(OBJS) $(FIT_ONE)
	$(CXX) -o fit_single.exe $(OBJS) $(FIT_ONE) $(LDFLAGS) cmpfit-1.2/mpfit.o


.cpp.o:
	$(CXX) -c $(CPPFLAGS) $<	
