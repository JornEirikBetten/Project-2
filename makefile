# Comment lines
# General makefile for c++ - choose PROG = name of given program
# Here we define compiler option, libraries and the target 
CPPflags=c++ -O3
# Here we define the library functions needed 
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG=jacobisolver
${PROG} : 	main.o tridiagmat.o jacobi.o
		${CPPflags} main.o tridiagmat.o jacobi.o ${LIB} -o ${PROG}
main.o : 				main.cpp
					${CPPflags} -c main.cpp
tridiagmat.o : 		tridiagmat.cpp
					${CPPflags} -c tridiagmat.cpp
jacobi.o : 				jacobi.cpp
					${CPPflags} -c jacobi.cpp


		
