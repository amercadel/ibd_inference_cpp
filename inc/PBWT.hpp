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

using namespace std;

struct VCFReader {
	ifstream vcf;
	int G, M, p1 = 0; // p1 points to first site in gap
	vector<vector<int> > gap; // stores haplotype data in the gap
	vector<string> ID, fixedFields;

	VCFReader(string file, int _G, int _M, ofstream& out) {
		vcf = ifstream(file);
		G = max(_G, 1), M = _M; // gap size of 0 is equivalent to a gap size of 1 when handling the VCF file and gap
		gap = vector<vector<int> >(G, vector<int>(M));
		ID = vector<string>(M);
		fixedFields = vector<string>(G);
		preprocess(out);
		initGap();
	}

	VCFReader(string file, int _G, int _M) {
		vcf = ifstream(file);
		G = max(_G, 1), M = _M; // gap size of 0 is equivalent to a gap size of 1 when handling the VCF file and gap
		gap = vector<vector<int> >(G, vector<int>(M));
		ID = vector<string>(M);
		fixedFields = vector<string>(G);
		preprocess();
		initGap();
	}
	
	// gets haplotype IDs and moves input stream pointer to start of raw data
	void preprocess(ofstream& out);

	void preprocess();

	// initializes sliding window for the gap
	void initGap();

	// reads the next site in the VCF file
	void nextSite();

	int getAllele(int g, int idx) {
		return gap[(p1 + g) % G][idx];
	}

	string getFixedField(int g);

	void editGap(int g, int idx, int allele);

	void close();
};


struct biPBWT {
	int M, site, N, G, L, W;
    float rho;
	vector<int> pre, div, backwardPre, backwardDiv;
	vector<int> a, b, d, e;
	vector<int> block, rBlock;

	biPBWT() {}
	biPBWT(int _M, int N, int G, int L, int W, float rho) {
		this->M = _M;
        this->rho = rho;
        this->N = N;
        this->G = G;
        this->L = L;
        this->W = W;
		site = 0;
		pre = vector<int>(M);
		iota(pre.begin(), pre.end(), 0);
		div = vector<int>(M);
		backwardPre = vector<int>(M), backwardDiv = vector<int>(M);
		a = vector<int>(M), b = vector<int>(M), d = vector<int>(M), e = vector<int>(M);
		block = vector<int>(M), rBlock = vector<int>(M);
	}

	void nextSite(VCFReader& vcf, ifstream& backward);

	void countingSort(vector<vector<int> >& v, int idx);

	// gets blocks at the current site and processes them
	void getBlocks(VCFReader& orig, VCFReader& edit);

	void processBlock(vector<vector<int> >& link, int start, int end, VCFReader& orig, VCFReader& edit);
};


int runPBWT(string input_vcf, string writeTo, string genetic_map_file, int checkpoint, int L, int W, int G, float rho);