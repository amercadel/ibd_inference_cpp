#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <map>
#include <cassert>
#include <cstring>
#include <algorithm>



struct VCFReader {
	std::ifstream vcf;
	int G, M, p1 = 0; // p1 points to first site in gap
	std::vector<std::vector<int>> gap; // stores haplotype data in the gap
	std::vector<std::string> ID, fixedFields;

	VCFReader(std::string file, int _G, int _M, std::ofstream& out) {
		vcf = std::ifstream(file);
		this->G = std::max(_G, 1);
		this->M = _M; // gap size of 0 is equivalent to a gap size of 1 when handling the VCF file and gap
		gap = std::vector<std::vector<int>>(this->G, std::vector<int>(this->M));
		ID = std::vector<std::string>(this->M);
		fixedFields = std::vector<std::string>(this->G);
		preprocess(out);
		initGap();
	}

	VCFReader(std::string file, int _G, int _M) {
		vcf = std::ifstream(file);
		this->G = std::max(_G, 1);
		this->M = _M; // gap size of 0 is equivalent to a gap size of 1 when handling the VCF file and gap
		gap = std::vector<std::vector<int>>(this->G, std::vector<int>(this->M));
		ID = std::vector<std::string>(this->M);
		fixedFields = std::vector<std::string>(G);
		preprocess();
		initGap();
	}
	
	// gets haplotype IDs and moves input stream pointer to start of raw data
	void preprocess(std::ofstream& out);

	void preprocess();

	// initializes sliding window for the gap
	void initGap();

	// reads the next site in the VCF file
	void nextSite();

	int getAllele(int g, int idx);

	std::string getFixedField(int g);

	void editGap(int g, int idx, int allele);

	void close();
};

struct biPBWT {
	int M, site, N, L, W, G; ;
	double rho;
	std::vector<int> pre, div, backwardPre, backwardDiv;
	std::vector<int> a, b, d, e;
	std::vector<int> block, rBlock;

	biPBWT() {}
	biPBWT(int _M) {
		this->M = _M;
		site = 0;
		pre = std::vector<int>(this->M);
		iota(pre.begin(), pre.end(), 0);
		div = std::vector<int>(this->M);
		backwardPre = std::vector<int>(this->M), backwardDiv = std::vector<int>(this->M);
		a = std::vector<int>(this->M), b = std::vector<int>(this->M), d = std::vector<int>(this->M), e = std::vector<int>(this->M);
		block = std::vector<int>(this->M), rBlock = std::vector<int>(this->M);
	}

	void nextSite(VCFReader& vcf, std::ifstream& backward);

	void countingSort(std::vector<std::vector<int>>& v, int idx);

	// gets blocks at the current site and processes them
	void getBlocks(VCFReader& orig, VCFReader& edit);

	void processBlock(std::vector<std::vector<int>>& link, int start, int end, VCFReader& orig, VCFReader& edit);
};

int runPBWT(std::string input, std::string writeTo, std::string genetic_map, int checkpoint, int L, int W, int G, double rho);