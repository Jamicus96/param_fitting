RAT_CONFIG = `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux
VPATH = ./antinuFit/:./antinuFit/models/:./scripts/:./temp_files/
INCLUDE = ${RAT_CONFIG} -I./antinuFit -I./antinuFit/models -I./scripts
COMPILE = g++ -g -std=c++1y
# error: undefined reference to '__dso_handle':
# "If I’m not mistaken, it is related to a combination of complex C++ object destruction of static objects and the nostdlib compiler option.
# In an embedded system, you likely don’t need the destruction of static objects. So try this compiler option: -fno-use-cxa-atexit"

all: make_PDFs re_combine_fits fit_params clean


# make_PDFs
make_PDFs: make_PDFs.cpp
	${COMPILE} scripts/make_PDFs.cpp -o scripts/make_PDFs.exe ${INCLUDE}

# re_combine_fits
re_combine_fits: E_systematics.o model_alphaN.o model_geoNu.o model_Reactor.o fitter.o re_combine_fits.o 
	${COMPILE} E_systematics.o model_alphaN.o model_geoNu.o model_Reactor.o fitter.o re_combine_fits.o -o scripts/re_combine_fits.exe ${INCLUDE}

re_combine_fits.o: re_combine_fits.cpp fitter.hpp fitVars.hpp model.hpp model_alphaN.hpp model_geoNu.hpp model_Reactor.hpp utils.hpp
	${COMPILE} scripts/re_combine_fits.cpp -c ${INCLUDE}

# fit_params
fit_params: E_systematics.o model_alphaN.o model_geoNu.o model_Reactor.o fitter.o fit_params.o 
	${COMPILE} E_systematics.o model_alphaN.o model_geoNu.o model_Reactor.o fitter.o fit_params.o -o scripts/fit_params.exe ${INCLUDE}

fit_params.o: fit_params.cpp fitter.hpp fitVars.hpp model.hpp model_alphaN.hpp model_geoNu.hpp model_Reactor.hpp utils.hpp
	${COMPILE} scripts/fit_params.cpp -c ${INCLUDE}

# antinuFit
fitter.o: fitter.cpp fitter.hpp fitVars.hpp model.hpp
	${COMPILE} antinuFit/fitter.cpp -c ${INCLUDE}

E_systematics.o: E_systematics.cpp E_systematics.hpp fitter.hpp fitVars.hpp
	${COMPILE} antinuFit/E_systematics.cpp -c ${INCLUDE}

model_alphaN.o: model_alphaN.cpp model_alphaN.hpp model.hpp fitter.hpp E_systematics.hpp fitVars.hpp
	${COMPILE} antinuFit/models/model_alphaN.cpp -c ${INCLUDE}

model_geoNu.o: model_geoNu.cpp model_geoNu.hpp model.hpp fitter.hpp E_systematics.hpp fitVars.hpp
	${COMPILE} antinuFit/models/model_geoNu.cpp -c ${INCLUDE}

model_Reactor.o: model_Reactor.cpp model_Reactor.hpp model.hpp fitter.hpp E_systematics.hpp fitVars.hpp
	${COMPILE} antinuFit/models/model_Reactor.cpp -c ${INCLUDE}



# Cleaning
clean:
	rm *.o

delete :
	rm scripts/*.exe
