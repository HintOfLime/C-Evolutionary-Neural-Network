.DEFAULT_GOAL := all

src=$(wildcard *.cpp)
obj=$(src:.cpp=.o)

all: raytracer.bin clean

run: all
	./main.bin

clean:
	rm *.o

%.o: %.cpp
	g++ -c $< -o $@

raytracer.bin: $(obj)
	g++ $(obj) -o main.bin -std=c++11 -lsfml-graphics -lsfml-window -lsfml-system

