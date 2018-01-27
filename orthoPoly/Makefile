CXX=g++ -std=c++11 -fopenmp -I$(PWD) -O3 -DNDEBUG

orthoPoly1d.o: orthoPoly1D.H orthoPoly1D.C
	$(CXX) -c orthoPoly1D.C -o orthoPoly1D.o

orthoPoly: orthoPoly3D.H orthoPoly3D.C orthoPoly1D.o
	$(CXX) -c orthoPoly3D.C -o orthoPoly3D.o
	$(CXX) orthoPoly3D.o orthoPoly1D.o -o orthoPoly

test: tmp.C
	$(CXX) tmp.C -o tmp

all:
	make orthoPoly1d.o
	make orthoPoly
	make test

clean:
	rm -f orthoPoly tmp *.o
