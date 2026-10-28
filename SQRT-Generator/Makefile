CXX      = g++
CXXFLAGS = -O2 -std=c++17 -fopenmp

ROOTS_6PT = /home/dschiebe/Documents/PhD-Files/FullHexaBoxDE/222-228-Doublebox_7MI/Roots.txt
GRAM_6PT  = /home/dschiebe/Documents/PhD-Files/FullHexaBoxDE/229-231-PentaBox_3MI/Reconstruction/SixPointGram.txt

all: gen_points gen_6pt

gen_points: gen_points.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

gen_6pt: gen_6pt.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

# Quick test targets
test-simple:
	./gen_points 5 42 --roots

test-6pt:
	./gen_6pt --roots $(ROOTS_6PT) --gram $(GRAM_6PT) --npoints 5 --sqrt

clean:
	rm -f gen_points gen_6pt
