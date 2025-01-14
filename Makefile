RAT_CONFIG = `root-config --cflags --libs` -I${RATROOT}/include/libpq -I${RATROOT}/include -I${RATROOT}/include/external -L${RATROOT}/lib -lRATEvent_Linux
VPATH = ./antinuFit/:./antinuFit/models/:./fitting/:./cutting/:./testing/
INCLUDE = ${RAT_CONFIG} -I./antinuFit -I./antinuFit/models -I./fitting -I./cutting
COMPILE = g++ -g -std=c++1y
# error: undefined reference to '__dso_handle':
# "If I’m not mistaken, it is related to a combination of complex C++ object destruction of static objects and the nostdlib compiler option.
# In an embedded system, you likely don’t need the destruction of static objects. So try this compiler option: -fno-use-cxa-atexit"

all: cutting fitting tests clean


# Cutting and PDF making
cutting: make_PDFs cut_data cut_ratds get_accidentals_and_vetos reScaleReactorIBD clean

make_PDFs: make_PDFs.cpp 
	${COMPILE} cutting/make_PDFs.cpp -o cutting/make_PDFs.exe ${INCLUDE}

cut_data: cut_data.o 
	${COMPILE} cut_data.o -o cutting/cut_data.exe ${INCLUDE}
cut_data.o: cut_data.cpp cutting_utils.hpp
	${COMPILE} cutting/cut_data.cpp -c ${INCLUDE}

cut_ratds: cut_ratds.o 
	${COMPILE} cut_ratds.o -o cutting/cut_ratds.exe ${INCLUDE}
cut_ratds.o: cut_ratds.cpp cutting_utils.hpp
	${COMPILE} cutting/cut_ratds.cpp -c ${INCLUDE}

get_accidentals_and_vetos: get_accidentals_and_vetos.o 
	${COMPILE} get_accidentals_and_vetos.o -o cutting/get_accidentals_and_vetos.exe ${INCLUDE}
get_accidentals_and_vetos.o: get_accidentals_and_vetos.cpp cutting_utils.hpp
	${COMPILE} cutting/get_accidentals_and_vetos.cpp -c ${INCLUDE}

reScaleReactorIBD: reScaleReactorIBD.cpp 
	${COMPILE} cutting/reScaleReactorIBD.cpp -o cutting/reScaleReactorIBD.exe ${INCLUDE}

# Fitting
fitting: fit_params re_combine_fits clean

re_combine_fits: re_combine_fits.o 
	${COMPILE} re_combine_fits.o -o fitting/re_combine_fits.exe ${INCLUDE}
re_combine_fits.o: re_combine_fits.cpp fitter.hpp fitVars.hpp E_systematics.hpp model_alphaN.hpp model_geoNu.hpp model_Reactor.hpp fitting_utils.hpp
	${COMPILE} fitting/re_combine_fits.cpp -c ${INCLUDE}

fit_params: fit_params.o 
	${COMPILE} fit_params.o -o fitting/fit_params.exe ${INCLUDE}
fit_params.o: fit_params.cpp fitter.hpp fitVars.hpp E_systematics.hpp model_alphaN.hpp model_geoNu.hpp model_Reactor.hpp fitting_utils.hpp cutting_utils.hpp
	${COMPILE} fitting/fit_params.cpp -c ${INCLUDE}


# Testing
tests: test_Esys clean

test_Esys: test_Esys.o 
	${COMPILE} test_Esys.o -o testing/test_Esys.exe ${INCLUDE}
test_Esys.o: test_Esys.cpp fitVars.hpp E_systematics.hpp
	${COMPILE} testing/test_Esys.cpp -c ${INCLUDE}


# Cleaning
clean:
	rm *.o

delete :
	rm fitting/*.exe cutting/*.exe testing/*.exe
