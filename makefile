.DEFAULT_GOAL := all

src=$(wildcard *.cpp)
obj=$(src:.cpp=.o)

all: main.bin clean

run: all
	./main.bin

clean:
	rm *.o

%.o: %.cpp
	g++ -c $< -o $@

main.bin: $(obj)
	g++ $(obj) -o main.bin -std=c++11 -lsfml-graphics -lsfml-window -lsfml-system

