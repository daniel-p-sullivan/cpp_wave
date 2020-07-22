OBJS = CFD_CPP.o FileWriter.o DomainFunctions.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

CFD_CPP_2D : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o CFD_CPP_2D

CFD_CPP.o : CFD_CPP.cpp main.h Configuration.h DomainFunctions.h FileWriter.h
	$(CC) $(CFLAGS) CFD_CPP.cpp

FileWriter.o : FileWriter.cpp Configuration.h main.h
	$(CC) $(CFLAGS) FileWriter.cpp

DomainFunctions.o : DomainFunctions.cpp main.h DomainFunctions.h Configuration.h
	$(CC) $(CFLAGS) DomainFunctions.cpp

clean :
	\rm *.o  CFD_CPP_2D
