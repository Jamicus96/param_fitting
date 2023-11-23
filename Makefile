RAT_CONFIG = `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux -lRooFitCore -lRooFit
VPATH=./antinuFit/:./antinuFit/models/:./scripts/

all: make_PDFs re_combine_fits fit_params


make_PDFs: make_PDFs.cpp
	g++ -g -std=c++1y make_PDFs.cpp -o make_PDFs.exe ${RAT_CONFIG}

re_combine_fits: re_combine_fits.cpp
	g++ -g -std=c++1y re_combine_fits.cpp -o re_combine_fits.exe ${RAT_CONFIG}


fit_params: fit_params.o fitter.o fitVar.o model.o
	g++ -g -std=c++1y fit_params.o fitter.o fitVar.o model.o -o fit_params.exe ${RAT_CONFIG}

fit_params.o: fit_params.cpp fitter.hpp fitVar.hpp model.hpp model_alphaN.hpp model_geoNu.hpp model_Reactor.hpp
	g++ -g -std=c++1y -Wall fit_params.cpp ${RAT_CONFIG}

fitter.o: fitter.cpp fitter.hpp fitVar.hpp model.hpp
	g++ -g -std=c++1y -W -Wall fitter.cpp ${RAT_CONFIG}

fitVar.o: fitVar.cpp fitVar.hpp
	g++ -g -std=c++1y -W -Wall fitVar.cpp ${RAT_CONFIG}

model.o: model.cpp model.hpp fitVar.hpp
	g++ -g -std=c++1y -W -Wall model.cpp ${RAT_CONFIG}

model_alphaN.o: model_alphaN.hpp model.hpp fitVar.hpp
	g++ -g -std=c++1y -W -Wall model_alphaN.hpp ${RAT_CONFIG}

model_geoNu.o: model_geoNu.hpp model.hpp fitVar.hpp
	g++ -g -std=c++1y -W -Wall model_geoNu.hpp ${RAT_CONFIG}

model_Reactor.o: model_Reactor.cpp model_Reactor.hpp model.hpp fitVar.hpp
	g++ -g -std=c++1y -W -Wall model_Reactor.cpp ${RAT_CONFIG}

clean :
	rm *.o vec
