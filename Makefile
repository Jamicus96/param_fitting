RAT_CONFIG = `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux
VPATH = ./antinuFit/:./antinuFit/models/:./scripts/:./temp_files/
INCLUDE = ${RAT_CONFIG} -I./antinuFit -I./antinuFit/models -I./scripts
COMPILE = g++ -g -std=c++1y
# error: undefined reference to '__dso_handle':
# "If I’m not mistaken, it is related to a combination of complex C++ object destruction of static objects and the nostdlib compiler option.
# In an embedded system, you likely don’t need the destruction of static objects. So try this compiler option: -fno-use-cxa-atexit"

all: make_PDFs cut_data re_combine_fits fit_params clean


# Make PDFs
make_PDFs: make_PDFs.cpp
	${COMPILE} scripts/make_PDFs.cpp -o scripts/make_PDFs.exe ${INCLUDE}

# Cut Data
cut_data: cut_data.cpp
	${COMPILE} scripts/cut_data.cpp -o scripts/cut_data.exe ${INCLUDE}

# re_combine_fits
re_combine_fits: re_combine_fits.o 
	${COMPILE} re_combine_fits.o -o scripts/re_combine_fits.exe ${INCLUDE}

re_combine_fits.o: re_combine_fits.cpp fitter.hpp fitVars.hpp model_alphaN.hpp model_geoNu.hpp model_Reactor.hpp utils.hpp
	${COMPILE} scripts/re_combine_fits.cpp -c ${INCLUDE}

# fit_params
fit_params: fit_params.o 
	${COMPILE} fit_params.o -o scripts/fit_params.exe ${INCLUDE}

fit_params.o: fit_params.cpp fitter.hpp fitVars.hpp model_alphaN.hpp model_geoNu.hpp model_Reactor.hpp utils.hpp
	${COMPILE} scripts/fit_params.cpp -c ${INCLUDE}



# Cleaning
clean:
	rm *.o

delete :
	rm scripts/*.exe
