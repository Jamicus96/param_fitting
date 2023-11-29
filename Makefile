RAT_CONFIG = `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux
VPATH = ./antinuFit/:./antinuFit/models/:./scripts/:./temp_files/
INCLUDE = ${RAT_CONFIG} -I./antinuFit -I./antinuFit/models -I./scripts
COMPILE = g++ -g -std=c++1y
# error: undefined reference to '__dso_handle':
# "If I’m not mistaken, it is related to a combination of complex C++ object destruction of static objects and the nostdlib compiler option.
# In an embedded system, you likely don’t need the destruction of static objects. So try this compiler option: -fno-use-cxa-atexit"

all: fit_params clean



make_PDFs: make_PDFs.cpp
	${COMPILE} scripts/make_PDFs.cpp -o scripts/make_PDFs.exe ${INCLUDE}

re_combine_fits: re_combine_fits.cpp
	${COMPILE} scripts/re_combine_fits.cpp -o scripts/re_combine_fits.exe ${INCLUDE}


fit_params: fitVar.o model.o model_alphaN.o model_geoNu.o model_Reactor.o fitter.o fit_params.o 
	${COMPILE} fitVar.o model.o model_alphaN.o model_geoNu.o model_Reactor.o fitter.o fit_params.o -o scripts/fit_params.exe ${INCLUDE}

fit_params.o: fit_params.cpp fitter.hpp fitVar.hpp model.hpp model_alphaN.hpp model_geoNu.hpp model_Reactor.hpp
	${COMPILE} scripts/fit_params.cpp -c ${INCLUDE}

fitter.o: fitter.cpp fitter.hpp fitVar.hpp model.hpp
	${COMPILE} antinuFit/fitter.cpp -c ${INCLUDE}

fitVar.o: fitVar.cpp fitVar.hpp
	${COMPILE} antinuFit/fitVar.cpp -c ${INCLUDE}

model.o: model.cpp model.hpp fitVar.hpp
	${COMPILE} antinuFit/model.cpp -c ${INCLUDE}

model_alphaN.o: model_alphaN.cpp model_alphaN.hpp model.hpp fitVar.hpp
	${COMPILE} antinuFit/models/model_alphaN.cpp -c ${INCLUDE}

model_geoNu.o: model_geoNu.cpp model_geoNu.hpp model.hpp fitVar.hpp
	${COMPILE} antinuFit/models/model_geoNu.cpp -c ${INCLUDE}

model_Reactor.o: model_Reactor.cpp model_Reactor.hpp model.hpp fitVar.hpp
	${COMPILE} antinuFit/models/model_Reactor.cpp -c ${INCLUDE}



clean :
	rm *.o

delete :
	rm scripts/*.exe
