CXX    = g++
prefix = ../htslib

ScaffoldIwormContigs:
	$(CXX) -I${HOME}/include -L${HOME}/lib -I$(prefix) -L$(prefix) ScaffoldIwormContigs.cpp error_checker.cpp -lhts -lz -o scaffold_iworm_contigs

clean:
	rm scaffold_iworm_contigs


