RAT_CONFIG = `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux -lRooFitCore -lRooFit
VPATH = ./antinuFit/:./antinuFit/models/:./scripts/:./temp_files/
INCLUDE = ${RAT_CONFIG} -I./antinuFit -I./antinuFit/models -I./scripts
COMPILE = g++ -g -std=c++1y -shared -fPIC


all: make_PDFs re_combine_fits fit_params clean



make_PDFs: make_PDFs.cpp
	${COMPILE} scripts/make_PDFs.cpp -o scripts/make_PDFs.exe ${INCLUDE}

re_combine_fits: re_combine_fits.cpp
	${COMPILE} scripts/re_combine_fits.cpp -o scripts/re_combine_fits.exe ${INCLUDE}


fit_params: fitVar.o model.o model_alphaN.o model_geoNu.o model_Reactor.o fitter.o fit_params.o 
	${COMPILE} fit_params.o fitter.o fitVar.o model.o -o scripts/fit_params.exe ${INCLUDE}

fit_params.o: fit_params.cpp fitter.hpp fitVar.hpp model.hpp model_alphaN.hpp model_geoNu.hpp model_Reactor.hpp
	${COMPILE} scripts/fit_params.cpp -o fit_params.o ${INCLUDE}

fitter.o: fitter.cpp fitter.hpp fitVar.hpp model.hpp
	${COMPILE} antinuFit/fitter.cpp -o fitter.o ${INCLUDE}

fitVar.o: fitVar.cpp fitVar.hpp
	${COMPILE} antinuFit/fitVar.cpp -o fitVar.o ${INCLUDE}

model.o: model.hpp fitVar.hpp
	${COMPILE} antinuFit/model.hpp -o model.o ${INCLUDE}

model_alphaN.o: model_alphaN.hpp model.hpp fitVar.hpp
	${COMPILE} antinuFit/models/model_alphaN.hpp -o model_alphaN.o ${INCLUDE}

model_geoNu.o: model_geoNu.hpp model.hpp fitVar.hpp
	${COMPILE} antinuFit/models/model_geoNu.hpp -o model_geoNu.o ${INCLUDE}

model_Reactor.o: model_Reactor.cpp model_Reactor.hpp model.hpp fitVar.hpp
	${COMPILE} antinuFit/models/model_Reactor.cpp -o model_Reactor.o ${INCLUDE}



clean :
	rm *.o

delete :
	rm scripts/*.exe
