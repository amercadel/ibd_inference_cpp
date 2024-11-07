#include <PBWT.hpp>

using namespace std;

	
	// gets haplotype IDs and moves input stream pointer to start of raw data
void VCFReader::preprocess(ofstream& out) { 
	// skip meta-info lines and get header line
	string header;
	while (getline(vcf, header)) {
		out << header << '\n';
		if ((int)header.size() < 2 || header[0] != '#' || header[1] != '#') break;
	}
	// input sample IDs
	stringstream ss(header);
	for (int i = 0; i < 9; ++i) getline(ss, ID[0], '\t'); // skip fixed columns, assumes 9 columns (FORMAT column) 
	for (int i = 0; i < this->M / 2; ++i) {
		getline(ss, ID[2 * i], '\t');
		ID[2 * i + 1] = ID[2 * i] + "-1";
		ID[2 * i] += "-0";
	}
}

void VCFReader::preprocess() { 
	// skip meta-info lines and get header line
	string header;
	while (getline(vcf, header)) {
		if ((int)header.size() < 2 || header[0] != '#' || header[1] != '#') break;
	}

	// input sample IDs
	stringstream ss(header);
	for (int i = 0; i < 9; ++i) getline(ss, ID[0], '\t'); // skip fixed columns, assumes 9 columns (FORMAT column) 
	for (int i = 0; i < this->M / 2; ++i) {
		getline(ss, ID[2 * i], '\t');
		ID[2 * i + 1] = ID[2 * i] + "-1";
		ID[2 * i] += "-0";
	}
}

	// initializes sliding window for the gap
void VCFReader::initGap() {
	for (int i = 0; i < this->G; ++i) nextSite();
}

	// reads the next site in the VCF file
void VCFReader::nextSite() {
	fixedFields[p1].clear();
	char s[2 * this->M + 5000]; // assumes fixed fields take up less than 5000 characters
	vcf.getline(s, 2 * this->M + 5000);
	// skip fixed fields
	int offset = 0, tabs = 0; // offset = position in "s" of first sequence - points to the first character after 9 tabs
	while (tabs < 9) {
		if (s[offset] == '\t') ++tabs;
		if (tabs < 9) fixedFields[p1] += s[offset];
		++offset;
	}

	for (int i = 0; i < this->M; ++i) {
		assert(s[offset + (i / 2) * 4 + 1] == '|'); // sanity check
		gap[this->p1][i] = (s[offset + 2 * i] == '0' ? 0 : 1);
	}
	this->p1 = (this->p1 + 1) % G;
}

int VCFReader::getAllele(int g, int idx) {
	return gap[(this->p1 + g) % this->G][idx];
}

string VCFReader::getFixedField(int g) {
	return fixedFields[(this->p1 + g) % this->G];
}

void VCFReader::editGap(int g, int idx, int allele) {
	gap[(this->p1 + g) % this->G][idx] = allele;
}

void VCFReader::close() {vcf.close();}



void biPBWT::nextSite(VCFReader& vcf, ifstream& backward) {
	// update forward PBWT
	int u = 0, v = 0, p = site + 1, q = site + 1;
	for (int i = 0; i < M; ++i) {
		int id = pre[i];
		if (div[i] > p) p = div[i];
		if (div[i] > q) q = div[i];
		if (vcf.getAllele(0, id) == 0) {
			a[u] = id;
			d[u] = p;
			++u;
			p = 0;
		}
		else {
			b[v] = id;
			e[v] = q;
			++v;
			q = 0;
		}
	}
	for (int i = 0; i < u; ++i) {
		pre[i] = a[i];
		div[i] = d[i];
	}
	for (int i = 0; i < v; ++i) {
		pre[u + i] = b[i];
		div[u + i] = e[i];
	}
	++site;

	// update reverse PBWT
	int rsite = (this->N - 1) - this->site - this->G; // index of the corresponding reverse site
	backward.seekg((long long)rsite * this->M * 8);

	// initialize backwardPre, backwardDiv, and rBlock
	int id = 0;
	for (int i = 0; i < this->M; ++i) {
		backward.read((char*)&backwardPre[i], sizeof backwardPre[i]);
		int rDiv; backward.read((char*)&rDiv, sizeof rDiv);
		backwardDiv[i] = rDiv;
		rDiv = (this->N - 1) - rDiv; // get forward index for position comparision

		if (rDiv < this->site + (this->G - 1) + this->L) ++id;
		rBlock[backwardPre[i]] = id;
	}

	// initialize block
	id = 0;
	for (int i = 0; i < this->M; ++i) {
		if (div[i] > site - this->L) ++id;
		block[pre[i]] = id;
	}

}

void biPBWT::countingSort(vector<vector<int>>& v, int idx) {
	vector<vector<vector<int>>> table(this->M + 1);
	for (int i = 0; i < this->M; ++i) {
		table[v[i][idx]].push_back(v[i]);
	}
	int p = 0;
	for (int i = 0; i <= this->M; ++i) {
		for (int j = 0; j < (int)table[i].size(); ++j) {
			v[p++] = table[i][j];
		}
	}
}

// gets blocks at the current site and processes them
void biPBWT::getBlocks(VCFReader& orig, VCFReader& edit) {
	// Algorithm 2 - block matching
	vector<vector<int>> link(this->M, vector<int>(3)); // [sample ID, forward block ID, reverse block ID]
	for (int i = 0; i < this->M; ++i) {
		link[i][0] = i, link[i][1] = block[i], link[i][2] = rBlock[i];
	}
	// radix sort
	countingSort(link, 2);
	countingSort(link, 1);

	int start = 0;
	for (int i = 1; i < M; ++i) {
		if (link[i][1] != link[i - 1][1] || link[i][2] != link[i - 1][2]) {
			processBlock(link, start, i, orig, edit);
			start = i;	
		}
	}
	processBlock(link, start, this->M, orig, edit);
}

void biPBWT::processBlock(vector<vector<int>>& link, int start, int end, VCFReader& orig, VCFReader& edit) { // [start, end)
	int blockSize = end - start;
	if (blockSize < W) return; // width too small
	
	vector<int> zero(G), one(G);
	for (int j = start; j < end; ++j) {
		int id = link[j][0];
		for (int k = 0; k < G; ++k) {
			if (orig.getAllele(k, id) == 0) ++zero[k];
			else ++one[k];
		}
	}

	// error correct
	for (int k = 0; k < this->G; ++k) {
		if (min(zero[k], one[k]) <= blockSize * this->rho) {
			int correctAllele = zero[k] < one[k] ? 1 : 0;
			for (int j = start; j < end; ++j) {
				int id = link[j][0];
				edit.editGap(k, id, correctAllele);
			}
		}
	}
}



int runPBWT(string input, string writeTo, string genetic_map, int checkpoint, int L, int W, int G, double rho) {

	ifstream backward(writeTo + ".rpbwt"), meta(writeTo + ".meta"), gMap(genetic_map);
	ofstream out(writeTo + ".smooth.vcf");
	int M, N;

	// retrieve M and N from meta file
	meta >> M >> N;

	VCFReader orig(input, G, M, out), edit(input, G, M); // orignal vcf and editted vcf

	biPBWT bipbwt(M);
	int site;
	for (site = 0; site + G < N; ++site) {
		out << edit.getFixedField(0);
		for (int i = 0; i < M / 2; ++i) {
			out << '\t' << edit.getAllele(0, 2 * i) << '|' << edit.getAllele(0, 2 * i + 1);
		}
		out << '\n';
		bipbwt.nextSite(orig, backward);

		orig.nextSite();
		edit.nextSite();

		bipbwt.getBlocks(orig, edit);

		if (site % checkpoint == 0) cout << "Checkpoint " << site << endl;
	}
	for (int j = 0; j < G; ++j) {
		out << edit.getFixedField(0);
		for (int i = 0; i < M / 2; ++i) {
			out << '\t' << edit.getAllele(j, 2 * i) << '|' << edit.getAllele(j, 2 * i + 1);
		}
		out << '\n';
		++site;
	}

	orig.close();
	edit.close();
	meta.close();
	gMap.close();
	out.close();

	return 0;
}