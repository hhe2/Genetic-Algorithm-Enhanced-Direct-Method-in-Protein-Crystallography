a.out: HIO_boundary_hist_match.o
	mpiicpc -o a.out HIO_boundary_hist_match.o -mcmodel=large -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -O
HIO_boundary_hist_match.o: HIO_boundary_hist_match.cpp ../crysPara.hpp ../../header/*.hpp
	mpiicpc -c HIO_boundary_hist_match.cpp -mcmodel=large -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -O
