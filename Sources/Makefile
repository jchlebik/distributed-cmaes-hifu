###########################################################################################
# @file         Makefile
#
# @author       Jakub Chlebik \n
#               Faculty of Information Technology \n
#               Brno University of Technology \n
#               xchleb07@stud.fit.vutbr.com
#
# @brief        Makefile to build MPI and sequential version of CMAES optimzer.
#
# @version      0.1
#
# @date         2020-11-26 (created) \n
############################################################################################

USER:=$(shell whoami)
CC=g++
ICC=icpc
MPICC=mpiicc
CFLAGS=-std=c++11 -O3

LIBS=-mkl -ffast-math -xHost

LIBS_FFTW=-lfftw3f_omp -lfftw3f -lm -ffast-math

LIBS_OMP_I=-qopenmp
LIBS_OMP_G=-fopenmp
#-liomp5 -lpthread -L<icc_dir>/lib

USR_INCPATH=/home/$(USER)/libs/include
USR_LIBPATH=/home/$(USER)/libs/lib

LOGHOOKHEADERSPATH=./LoggingFunc/
HELPERSHEADERSPATH=./Helpers/
FITNESSHEADERSPATH=./FitnessFunc/
HIFUFITNESSPATH=$(FITNESSHEADERSPATH)HIFUScore/
BENCHFITNESSPATH=$(FITNESSHEADERSPATH)Benchmarks/

CMAESPPHEADERPATH=./ThirdParty/CMA-ESpp/cma-es/

HIFUDIFFUSIONSOURCEFILES=$(HIFUFITNESSPATH)KWaveDiffusionSolver/*.cc
HIFUSCORESOURCEFILES=$(HIFUFITNESSPATH)*.cc

HIFUSCOREPATH=$(HIFUFITNESSPATH)heatDiffusionScore.cc
BENCHSCOREPATH=$(BENCHFITNESSPATH)AckleyScore.cc

CMAESPPLOHGOOKSPATH=LoggingFuncs/LoggingHooks_CMAESpp.cc

.PHONY: hifu_test hifu_bcast hifu_pp hifu_circle hifu_seq bench_bcast bench_pp bench_circle bench_seq

all:	bench hifu

bench:  bench_bcast bench_pp bench_circle bench_seq

hifu: hifu_bcast hifu_pp hifu_circle hifu_seq

perf: perf_bcast perf_pp

perf_bcast: $(HIFUSCOREPATH) $(HIFUDIFFUSIONSOURCEFILES) $(HIFUSCORESOURCEFILES) $(CMAESPPLOHGOOKSPATH) Cmaes_MPI_alltoall.cc
	$(MPICC) $(CFLAGS) $^ $(LIBS) $(LIBS_OMP_I) -g -parallel-source-info=2 -o perfBcast

perf_pp: $(HIFUSCOREPATH) $(HIFUDIFFUSIONSOURCEFILES) $(HIFUSCORESOURCEFILES) $(CMAESPPLOHGOOKSPATH) Cmaes_MPI_PP.cc
	$(MPICC) $(CFLAGS) $^ $(LIBS) $(LIBS_OMP_I) -g -parallel-source-info=2 -o perfPp

hifu_test:	$(HIFUDIFFUSIONSOURCEFILES) $(HIFUSCORESOURCEFILES) $(HIFUSCOREPATH) hifuTest.cc
	$(ICC) -std=c++11 -g $^ $(LIBS) $(LIBS_OMP_I) -qopt-report=5 -qopt-report-phase=vec,par -vecabi=cmdtarget -o hifuTest

hifu_bcast: $(HIFUSCOREPATH) $(HIFUDIFFUSIONSOURCEFILES) $(HIFUSCORESOURCEFILES) $(CMAESPPLOHGOOKSPATH) Cmaes_MPI_alltoall.cc
	$(MPICC) $(CFLAGS) $^ $(LIBS) $(LIBS_OMP_I) -o hifuBcast

hifu_pp: $(HIFUSCOREPATH) $(HIFUDIFFUSIONSOURCEFILES) $(HIFUSCORESOURCEFILES) $(CMAESPPLOHGOOKSPATH) Cmaes_MPI_PP.cc
	$(MPICC) $(CFLAGS) $^ $(LIBS) $(LIBS_OMP_I) -o hifuPp

hifu_circle: $(HIFUSCOREPATH) $(HIFUDIFFUSIONSOURCEFILES) $(HIFUSCORESOURCEFILES) $(CMAESPPLOHGOOKSPATH) Cmaes_MPI_circle.cc
	$(MPICC) $(CFLAGS) $^ $(LIBS) $(LIBS_OMP_I) -o hifuCircle

hifu_seq:	$(HIFUSCOREPATH) $(HIFUDIFFUSIONSOURCEFILES) $(HIFUSCORESOURCEFILES) $(CMAESPPLOHGOOKSPATH) Cmaes_SEQ.cc
	$(ICC) $(CFLAGS) $^ $(LIBS) $(LIBS_OMP_I) -o hifuSeq

bench_bcast: $(BENCHSCOREPATH) $(CMAESPPLOHGOOKSPATH) Cmaes_MPI_alltoall.cc
	$(MPICC) $(CFLAGS) $^ $(LIBS) -o benchBcast

bench_pp:	$(BENCHSCOREPATH) $(CMAESPPLOHGOOKSPATH) Cmaes_MPI_PP.cc
	$(MPICC) $(CFLAGS) $^ $(LIBS) -o benchPp

bench_circle: $(BENCHSCOREPATH) $(CMAESPPLOHGOOKSPATH) Cmaes_MPI_circle.cc
	$(MPICC) $(CFLAGS) $^ $(LIBS) -o benchCircle

bench_seq: $(BENCHSCOREPATH) $(CMAESPPLOHGOOKSPATH) Cmaes_SEQ.cc
	$(ICC) $(CFLAGS) $^ $(LIBS) -o benchSeq

clean:
	rm -f benchSeq
	rm -f benchCircle
	rm -f benchPp
	rm -f benchBcast
	rm -f hifuSeq
	rm -f hifuCircle
	rm -f hifuPp
	rm -f hifuBcast
	rm -f hifuTest
