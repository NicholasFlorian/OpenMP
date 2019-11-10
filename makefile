
main: main.cc
	g++ main.cc -fopenmp -o main -lm 

clean:
	rm main
