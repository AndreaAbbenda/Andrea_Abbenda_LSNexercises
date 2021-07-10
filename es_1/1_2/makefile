LIBS:=`root-config --libs`
INCS:=`root-config --cflags`


main: main.o random.o
	g++ -no-pie main.o random.o -o main ${INCS} ${LIBS}
main.o: main.C random.h
	g++ -c main.C -o main.o ${INCS}
random.o: random.cpp random.h
	g++ -c random.cpp -o random.o
