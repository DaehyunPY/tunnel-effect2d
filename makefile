# q3 square tunneling

TARGET = a.out
OBJECTS = sde.o const.o main.o
COMMON_MOD = sde.f08 const.f08 main.f08

F08 = /usr/local/bin/gfortran
FFLAGS = -fimplicit-none -fbounds-check

INCLUDE = -I/usr/local/include
LIB = -L/usr/local/lib
LAPACK = -llapack -lblas
FFW = -lfftw3 -lm





all: clean ${TARGET}

${TARGET}: ${OBJECTS}
	${F08} -o $@ ${OBJECTS}

${OBJECTS}: ${COMMON_MOD}
	${F08} -c ${COMMON_MOD} ${FFLAGS}





clean: ;rm -f *.o *.mod
