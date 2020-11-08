

cwd = $(shell pwd)

sift_bam_max_cov: sift_bam_max_cov.cpp htslib/version.h
	g++ -std=c++11 -o _sift_bam_max_cov sift_bam_max_cov.cpp -Wall -O2 -L./htslib/build/lib/ -I./htslib/build/include -lhts


htslib/version.h :
	./build_htslib.sh


