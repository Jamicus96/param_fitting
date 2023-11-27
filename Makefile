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
	${COMPILE} -nostartfiles fitVar.o model.o model_alphaN.o model_geoNu.o model_Reactor.o fitter.o fit_params.o -o scripts/fit_params.exe ${INCLUDE}

fit_params.o: fit_params.cpp fitter.hpp fitVar.hpp model.hpp model_alphaN.hpp model_geoNu.hpp model_Reactor.hpp
	${COMPILE} scripts/fit_params.cpp -o fit_params.o ${INCLUDE}

fitter.o: fitter.hpp fitVar.hpp model.hpp
	${COMPILE} -nostartfiles antinuFit/fitter.hpp -o fitter.o ${INCLUDE}

fitVar.o: fitVar.hpp
	${COMPILE} -nostartfiles antinuFit/fitVar.hpp -o fitVar.o ${INCLUDE}

model.o: model.hpp fitVar.hpp
	${COMPILE} -nostartfiles antinuFit/model.hpp -o model.o ${INCLUDE}

model_alphaN.o: model_alphaN.hpp model.hpp fitVar.hpp
	${COMPILE} -nostartfiles antinuFit/models/model_alphaN.hpp -o model_alphaN.o ${INCLUDE}

model_geoNu.o: model_geoNu.hpp model.hpp fitVar.hpp
	${COMPILE} -nostartfiles antinuFit/models/model_geoNu.hpp -o model_geoNu.o ${INCLUDE}

model_Reactor.o: model_Reactor.hpp model.hpp fitVar.hpp
	${COMPILE} -nostartfiles antinuFit/models/model_Reactor.hpp -o model_Reactor.o ${INCLUDE}



clean :
	rm *.o

delete :
	rm scripts/*.exe
