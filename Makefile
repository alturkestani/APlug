# Tariq Alturkestani
# Adaptive Parallelism for Scientific Applications
	# makefile for APlug stand alone driver 

SRC = main.cpp 
LIB = APlug.cpp 
INC = APlug.h 

main : ${SRC} ${LIB} ${INC}
		g++ -o APlug -Wall -Wextra -O3 ${SRC} ${LIB}

GUI_SRC = gui/main_gui.cpp

gui : ${GUI_SRC} ${LIB} ${INC}
	g++ -o APlugGUI -std=c++11 -Wall -O3 ${GUI_SRC} ${LIB} \
	$(shell pkg-config --cflags --libs Qt5Widgets)

.PHONY: gui
	
clean:
		rm -f APlug APlugGUI *.o

