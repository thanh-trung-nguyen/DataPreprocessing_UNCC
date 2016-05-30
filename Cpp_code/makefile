CC = g++

OBJ1 = main.o MeasuredData.o InputOutput.o VecMatCal.o
OBJ2 = waveeq1d_forwardproblem.o WaveEq1D.o Grid1D.o InputOutput.o VecMatCal.o

data_preprocessing:  $(OBJ1)
	$(CC) -o data_preprocessing $(OBJ1)

waveeq1d:  $(OBJ2)
	$(CC) -o waveeq1d $(OBJ2)



all: $(OBJ1)

сlean:
	rm $(OBJ1)
