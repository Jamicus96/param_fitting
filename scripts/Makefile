RAT_CONFIG = `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux -lRooFitCore -lRooFit

all: make_PDFs.exe fit_params.exe re_combine_fits.exe

make_PDFs.exe: make_PDFs.cpp
	g++ -g -std=c++1y make_PDFs.cpp -o make_PDFs.exe ${RAT_CONFIG}

fit_params.exe: fit_params.cpp
	g++ -g -std=c++1y fit_params.cpp fitting_utils.hpp fitting_utils.cpp -o fit_params.exe ${RAT_CONFIG}

re_combine_fits.exe: re_combine_fits.cpp
	g++ -g -std=c++1y re_combine_fits.cpp -o re_combine_fits.exe ${RAT_CONFIG}

clean :
	rm *.exe
