# Makefile for PA2
CC = cc

default: heat2dcol

heat2d: heat2d.c svalidate.c svalidate.h heat2dhelper.c heat2dhelper.h
	mpicc heat2d.c svalidate.c heat2dhelper.c -lm -o heat2d

heat2dcol: heat2dcol.c svalidate.c svalidate.h heat2dhelper.c heat2dhelper.h
	mpicc heat2dcol.c svalidate.c heat2dhelper.c -lm -o heat2dcol

heat2drow: heat2drow.c svalidate.c svalidate.h heat2dhelper.c heat2dhelper.h
	mpicc heat2drow.c svalidate.c heat2dhelper.c -lm -o heat2drow

clean:
	- /bin/rm heat2d
	- /bin/rm heat2drow
	- /bin/rm heat2dcol
