# q3 square tunneling

TARGET = a.out
OBJECTS = const.o globals.o subfuncs.o subprogs.o main.o
COMMON_MOD = const.f90 globals.f90 subfuncs.f90 subprogs.f90 main.f90

FORTRAN = /usr/local/bin/gfortran 
FFLAGS = -fimplicit-none -fbounds-check 
# FORTRAN = /opt/intel/bin/ifort
# FORTRAN = ifort



# LIB = -L/usr/local/lib
# INCLUDE = -I/usr/local/include
# LAPACK = -lscalapack -lmpi -lblas

# LIB = -L/usr/local/opt/lapack/lib
# INCLUDE = -I/usr/local/opt/lapack/include
# LIB = -L/usr/lib
# INCLUDE = -I/usr/include
LAPACK = -llapack

# LIB = -L/usr/lib
# INCLUDE = -I/usr/include
# LAPACK = -lopenblas

# FFW = -lfftw3 -lm





all: clean ${TARGET}

${TARGET}: ${OBJECTS}
	# ${FORTRAN} -o $@ ${OBJECTS} 
	${FORTRAN} -o $@ ${OBJECTS} ${LAPACK}
	# ${FORTRAN} -o $@ ${OBJECTS} ${INCLUDE} ${LIB} ${LAPACK}

${OBJECTS}: ${COMMON_MOD}
	${FORTRAN} -c ${COMMON_MOD} 
	# ${FORTRAN} -c ${COMMON_MOD} ${FFLAGS} 





clean: ;rm -f *.o *.mod
