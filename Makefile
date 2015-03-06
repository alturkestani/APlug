# Tariq Alturkestani
# Adaptive Parallelism for Scientific Applications
# makefile for APlug stand alone driver 

SRC = main.cpp 
LIB = APlug.cpp 
INC = APlug.h 

main : ${SRC} ${LIB} ${INC}
	g++ -o APlug -Wall -Wextra -O3 ${SRC} ${LIB} 
	
clean:
	rm -f APlug *.o

