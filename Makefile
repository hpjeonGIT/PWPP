#
# Makefile for plane wave basis pseudopotential DFT method
#
# Byoungseon Jeon
# Dept. Applied Science, UC.Davis
# Oct.30. 2008
#
# suffix definition for compiling
.SUFFIXES: .o .cc
#
# select compiler
CXX = g++
#
# select compiling option - add optimization options
FLAGS = -Wall  -O3 #-g
LIB=   -lm  -L/usr/local/lib -lfftw3 -framework vecLib
INC= -I/usr/local/include
#
# Object files
OBJTS = main.o density.o hamiltn.o solver.o
TARGET = pwpp
#
# generation of executable
${TARGET}:${OBJTS}
	${CXX} ${FLAGS} -o ${TARGET} ${OBJTS} ${LIB}
#
# generation of object files
.cc.o:
	${CXX} ${FLAGS} -c $< ${INC}
#
# clean object files and executable
clean:
	rm -rf *.o  *.cc~ core ${TARGET}
