// * History
// * ...
// * Dec 8, 2018. decompression
// * ... *//
# include <cstdio>
# include <iostream>
# include <fstream>
# include <vector>
# include <map>
# include <queue>
# include <string>
# include <cstring>
# include <chrono>
# include <thread>
# include <mutex>
# include <cstdlib>
# include "libbsc/bsc.h"
using namespace std;

int L;
int ispe, isorder;
size_t max_rid;
const int nthreads = 24;
string folder;

struct READS_t {
	string str;
	size_t prid; //reads id of parent
	READS_t() {}
	READS_t(char *_str, const size_t &_prid): str(_str), prid(_prid) {}
};

char complement[256];

const char invert_code_rule[4] = {'A', 'C', 'G', 'T'}; //decoding rule //A - 0; C - 1; G - 2; T - 3;

class CStopWatch {
	private:
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	std::chrono::steady_clock::time_point t_temp;
	int running;
	std::vector<double> elapsed;

	public:
	CStopWatch();
	~CStopWatch();
	void start();
	double stop();
	void resume();
	double totalTime();
};

CStopWatch stopwatch;

void init() {
	// ht[1] = new HashMap;
	memset(complement, 125, 256);
	complement['A'] = 'T';
	complement['a'] = 'T';
	complement['C'] = 'G';
	complement['c'] = 'G';
	complement['G'] = 'C';
	complement['g'] = 'C';
	complement['T'] = 'A';
	complement['t'] = 'A';
	complement['N'] = 'N';
	complement['n'] = 'N';
}

inline void reverseComplement(char* start) {
	char* left = start; // sequence starts
	char* right = start + L - 1;
	while (right > left) {
		char tmp = complement[(uint8_t)*left];
		*left = complement[(uint8_t)*right];
		*right = tmp;
		++left;
		--right;
	}
	if (left == right)
		*left = complement[(uint8_t)*left];
}

inline void decodeStr(string pstr, char *enstr, char *dstr) {
	bool debug = false;
	// if (strcmp(enstr, "\n") == 0) debug = true;
	// if (debug) {
	// 	fprintf(stderr, "special pstr: %s\n", pstr.c_str());
	// 	exit(0);
	// }
	// if (enstr[0] == ' ' && enstr[1] == '\n') { // equal
	if (enstr[0] == '\n') { // equal
		strcpy(dstr, pstr.c_str());
	} else {
		// if (strcmp(enstr, "6GAGAT32T3C7T30A TCATATACAGG\n") == 0) debug = true;
		// if (strcmp(enstr, "88CTCGAGACAGAGT\n") == 0) debug = true;
		// if (strcmp(enstr, "AAAA 16G5G6G20A5A2T4C19A9C\n") == 0) debug = true;
		// if (strcmp(enstr, " CTACCATATTCGCCACAGCCCCA\n") == 0) debug = true;
		// if (strcmp(enstr, "TG TGAAGGCGG\n") == 0) debug = true;
		// if (strcmp(enstr, "TGACACGAAACAAGGAAGT \n") == 0) debug = true;
		// if (strcmp(enstr, "C AT\n") == 0) debug = true;
		// if (strcmp(enstr, "TG TGAAGGCGG\n") == 0) debug = true;
		// check space ' '
		int spidx = -1;
		int len = strlen(enstr) - 1;
		for (int i = 0; i < len; ++i) {
			if (enstr[i] == ' ') {
				spidx = i;
				break;
			}
		}
		int dstridx = 0;
		// fprintf(stderr, "--spidx: %d\n", spidx);

		/*if (spidx == 0) {
			for (int pstridx = len - 1; pstr[pstridx]; ++pstridx) {
				dstr[dstridx++] = pstr[pstridx];
			}
			++spidx;
			for (; spidx < len; ++spidx) {
				dstr[dstridx++] = enstr[spidx];
			}
			dstr[dstridx++] = '\0';
			return;
		} else 
		if (spidx == len - 1) {
			for (int i = 0; i < spidx; ++i) {
				dstr[dstridx++] = enstr[i];
			}
			for (int i = 0; dstridx < L; ++i) {
				dstr[dstridx++] = pstr[i];
			}
			dstr[dstridx++] = '\0';
			return;
		}*/

		// fprintf(stderr, "pstr: %s\n", pstr.c_str());
		bool shiftRight = false;
		if (spidx == 0) {
			shiftRight = true;
		} else 
		if (spidx > -1) {
			for (int i = 0; i < spidx; ++i) {
				if (enstr[i] >= '0' && enstr[i] <= '9') {
					shiftRight = true;
					break;
				}
			}
		}

		if (spidx > -1 && spidx != len - 1 && !shiftRight) {
			bool flag = true;
			for (int i = spidx + 1; i < len; ++i) {
				if (enstr[i] >= '0' && enstr[i] <= '9') {
					flag = false;
					break;
				}
			}
			if (flag) {
				shiftRight = true;
			}
		}

		int pstridx = 0, enstridx = 0, maxenstridx = len;
		if (spidx > -1) {
			if (shiftRight) {
				maxenstridx = spidx;
				pstridx = len - spidx - 1;
			} else {
				for (enstridx = 0; enstridx < spidx; ++enstridx) {
					dstr[dstridx++] = enstr[enstridx];
				}
				++ enstridx;
				// if (spidx == len - 1) {
				// 	pstridx = spidx - 1;
				// }
			}
		}

		if (debug) {
			fprintf(stderr, "pstr: %s\n", pstr.c_str());
			fprintf(stderr, "spidx: %d\n", spidx);
			fprintf(stderr, "shiftRight: %d\n", shiftRight);
			fprintf(stderr, "len: %d\n", len);
			fprintf(stderr, "maxenstridx: %d\n", maxenstridx);
			fprintf(stderr, "dstridx: %d\n", dstridx);
			fprintf(stderr, "pstridx: %d\n", pstridx);
		}

		int eq_char_num = 0;
		for (; enstridx < maxenstridx; ++enstridx) {
			if (enstr[enstridx] >= 'A' && enstr[enstridx] <= 'Z') {
				if (eq_char_num > 0) {
					for (int j = 0; j < eq_char_num; ++j) {
						dstr[dstridx++] = pstr[pstridx++];
					}
					eq_char_num = 0;
				}
				dstr[dstridx++] = enstr[enstridx];
				if (debug) {
					fprintf(stderr, "dstr[%d]: %c\n", dstridx-1, dstr[dstridx-1]);
					fprintf(stderr, "enstr[%d]: %c\n", enstridx, enstr[enstridx]);
				}
				pstridx++;
			} else {
				eq_char_num = eq_char_num * 10 + enstr[enstridx] - '0';
			}
		}
		dstr[dstridx] = '\0';
		if (debug) {
			fprintf(stderr, "%s\n", dstr);
			fprintf(stderr, "maxenstridx: %d\n", maxenstridx);
		}
		for (; pstridx < L && dstridx < L; ++pstridx) {
			dstr[dstridx++] = pstr[pstridx];
		}
		if (shiftRight) {
			for (int i = maxenstridx + 1; i < len && dstridx < L; ++i) {
				dstr[dstridx++] = enstr[i];
			}
		} /*else {
			for (; dstridx < L; ++dstridx) {
				dstr[dstridx++] = pstr[pstridx++];
			}
		}*/
		dstr[L] = '\0';
		if (debug) {
			fprintf(stderr, "%s\n", dstr);
			// exit(1);
		}
	}
}

struct BITPOOL {
	int n, a[8];
};

inline int getDir(std::ifstream& fpdir, BITPOOL& curdirpool) {
	if (curdirpool.n >= 8) {
		uint8_t dirbin;
		fpdir.read((char*)&dirbin, sizeof(uint8_t));
		for (int i = 0; i < 8; ++i) {
		// for (int i = 7; i >= 0; --i) {
			curdirpool.a[i] = dirbin&1;
			dirbin >>= 1;
		}
		curdirpool.n = 0;
	}
	int res = curdirpool.a[curdirpool.n];
	++curdirpool.n;
	return res;
}

struct DNAPOOL {
	int n, a[4];
};

inline int getDNAcode(std::ifstream& fpif, DNAPOOL& curdnapool) {
	// fprintf(stderr, "yyyyyy\n");
	if (curdnapool.n >= 4) {
		uint8_t bin;
		// fprintf(stderr, "zzzz\n");
		fpif.read((char*)&bin, sizeof(uint8_t));
		// fprintf(stderr, "%d\n", bin);
		if (fpif.eof()) return -1;

		for (int i = 0; i < 4; ++i) {
			curdnapool.a[i] = bin&3;
			bin >>= 2;
		}
		curdnapool.n = 0;
	}
	int res = curdnapool.a[curdnapool.n];
	++curdnapool.n;
	return res;
}

inline bool getReads(std::ifstream& fp, DNAPOOL& curdnapool, char *r) {
	int code, i;
	for (i = 0; i < L; ++i) {
		code = getDNAcode(fp, curdnapool);
		if (code < 0) break;
		r[i] = invert_code_rule[code];
	}
	r[i] = '\0';
	if (i == L) return true;
	return false;
}

bool checkStrOrDig(char *str) {
	int len = strlen(str);
	// if (str[0] == ' ' && str[1] == '\n') { // ' \n'
	if (str[0] == '\n') { // ' \n'
		return true;
	}
	for (int i = 0; i < len - 1; ++i) {
		if (str[i] >= 'A' && str[i] <= 'Z') {
			return true;
		}
	}
	return false;
}

int mainxxx(int argc, char const *argv[]) {
	string parent = "AAAGGCCCCAGTTTGGCAGACCGACAGCGTGAATACCTTTTAGACATGATCCCTCCCCGGTCTATATCGCAGTCCATCAGTGGACAGAAATAACGCCTGTT";
	char enstr[1<<10], dstr[1<<10];
	strcpy(enstr, "C AT0\n");
	decodeStr(parent, enstr, dstr);
	fprintf(stderr, "%s\n", dstr);
	return 0;
}

static const char alphanum[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
inline std::string generateString(const std::string &chr, int length = 5) {
	srand(time(0));
	string res = chr + "_";
	for (int i = 0; i < length; ++i) {
		res += alphanum[rand() % 62];
	}
	return res;
}


void decompressSingle(string result) {
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string rootstrfn = folder + "rootstr.txt";
	mstcom::bsc::BSC_decompress((rootstrfn+".bsc").c_str(), rootstrfn.c_str());
	string orderfn;

	std::ifstream fporder;
	if (isorder) {
		orderfn = folder + string("order.bin");
		mstcom::bsc::BSC_decompress((orderfn+".bsc").c_str(), orderfn.c_str());
		fporder.open(orderfn, std::ios::binary);
	}
	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	string singletonfn = folder + "singleton.bin"; 
	mstcom::bsc::BSC_decompress((singletonfn+".bsc").c_str(), singletonfn.c_str());
	string readsnfn = folder + "readsN.txt";
	mstcom::bsc::BSC_decompress((readsnfn+".bsc").c_str(), readsnfn.c_str());

	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	BITPOOL curdirpool;
	curdirpool.n = 8;
	int dir;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");
	FILE *fproot = fopen(rootstrfn.c_str(), "r");
	// sprintf(name, "rootstr.bin");
	// std::ifstream fproot(name, std::ios::binary);

	FILE *fpdecstr = fopen(result.c_str(), "w");
	char str[1<<10], enstr[1<<10], dstr[1<<10];;

	// int num = 112;

	// DNAPOOL rootdnapool;
	// rootdnapool.n = 4;
	if (isorder) {
		cout << "max_rid: " << max_rid << endl;
		uint32_t rid;
		char **readsstr = new char*[max_rid];
		bool *flag = new bool[max_rid];
		memset(flag, false, sizeof(bool)*max_rid);

		while (fscanf(fproot, "%s", str) != EOF) {
		// while (getReads(fproot, rootdnapool, str)) {
			vector<READS_t> reads;
			READS_t r(str, 0);
			reads.push_back(r);
			int idx = 0;
			// fprintf(fpdecstr, "%s\n", str);
			fporder.read((char*)&rid, sizeof(uint32_t));
			readsstr[rid] = strdup(str);
			flag[rid] = true;

			// while(fscanf(fpenstr, "%s", enstr) != EOF) {
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					// fprintf(stderr, "xx enstr: %s\n", enstr);
					//
					decodeStr(reads[idx].str, enstr, dstr);
					//
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					// fprintf(fpdecstr, "%s\n", dstr);
					fporder.read((char*)&rid, sizeof(uint32_t));
					readsstr[rid] = strdup(dstr);
					flag[rid] = true;
					// --num;
					// fprintf(stderr, "dstr: %s\n", dstr);
					// if (num <= 0) exit(0);

					READS_t r(dstr, idx);
					reads.push_back(r);
					idx = reads.size() - 1;

				} else { // back
					int back_step = atoi(enstr);
					while (back_step > 0) {
						idx = reads[idx].prid;
						--back_step;
					}
				}
			}
			reads.clear();
		}

		fclose(fpenstr);
		fclose(fproot);
		// fproot.close();
		fpdir.close();
		
		FILE *fpreadsn = fopen(readsnfn.c_str(), "r");
		while (fscanf(fpreadsn, "%s", str) != EOF) {
			// fprintf(fpdecstr, "%s\n", str);
			fporder.read((char*)&rid, sizeof(uint32_t));
			readsstr[rid] = strdup(str);
			flag[rid] = true;
		}
		fclose(fpreadsn);

		ifstream fpsingleton(singletonfn, std::ios::binary);
		DNAPOOL curdnapool;
		curdnapool.n = 4;
		int code;
		char *seq = (char*)alloca((L + 3) * sizeof(char));
		// bool flag = true;
		for (size_t id = 0; id < max_rid; ++id) {
			if (!flag[id]) {
				for (int i = 0; i < L; ++i) {
					code = getDNAcode(fpsingleton, curdnapool);
					if (code < 0) {
						break;
					}
					seq[i] = invert_code_rule[code];
				}
				seq[L] = '\0';
				// fporder.read((char*)&rid, sizeof(uint32_t));
				readsstr[id] = strdup(seq);
				// flag[rid] = true;
				// fprintf(fpdecstr, "%s\n", seq);
			}
		}
		fporder.close();
		fpsingleton.close();
		for (size_t id = 0; id < max_rid; ++id) {
			fprintf(fpdecstr, "%s\n", readsstr[id]);
		}

	} else {
		while (fscanf(fproot, "%s", str) != EOF) {
			vector<READS_t> reads;
			READS_t r(str, 0);
			reads.push_back(r);
			int idx = 0;
			fprintf(fpdecstr, "%s\n", str);

			// while(fscanf(fpenstr, "%s", enstr) != EOF) {
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					// fprintf(stderr, "xx enstr: %s\n", enstr);
					//
					decodeStr(reads[idx].str, enstr, dstr);
					//
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					fprintf(fpdecstr, "%s\n", dstr);
					// --num;
					// fprintf(stderr, "dstr: %s\n", dstr);
					// if (num <= 0) exit(0);

					READS_t r(dstr, idx);
					reads.push_back(r);
					idx = reads.size() - 1;

				} else { // back
					int back_step = atoi(enstr);
					while (back_step > 0) {
						idx = reads[idx].prid;
						--back_step;
					}
				}
			}
			reads.clear();
		}

		fclose(fpenstr);
		fclose(fproot);
		// fproot.close();
		fpdir.close();
		ifstream fpsingleton(singletonfn, std::ios::binary);
		DNAPOOL curdnapool;
		curdnapool.n = 4;
		int code;
		char *seq = (char*)alloca((L + 3) * sizeof(char));
		bool flag = true;
		while (flag) {
			for (int i = 0; flag && i < L; ++i) {
				code = getDNAcode(fpsingleton, curdnapool);
				if (code < 0) {
					flag = false;
					break;
				}
				seq[i] = invert_code_rule[code];
			}
			if (!flag) break;
			seq[L] = '\0';
			fprintf(fpdecstr, "%s\n", seq);
		}
		fpsingleton.close();

		FILE *fpreadsn = fopen(readsnfn.c_str(), "r");
		while (fscanf(fpreadsn, "%s", str) != EOF) {
			fprintf(fpdecstr, "%s\n", str);
		}
		fclose(fpreadsn);
	}

	fclose(fpdecstr);
}

void decompressPE(string result) {
	fprintf(stderr, "decompressPE()...\n");
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string rootstrfn = folder + "rootstr.txt";
	mstcom::bsc::BSC_decompress((rootstrfn+".bsc").c_str(), rootstrfn.c_str());

	string orderfn = folder + string("order.bin");
	std::ifstream fporder;
	mstcom::bsc::BSC_decompress((orderfn+".bsc").c_str(), orderfn.c_str());
	fporder.open(orderfn, std::ios::binary);

	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string singletonfn = folder + "singleton.bin"; 
	mstcom::bsc::BSC_decompress((singletonfn+".bsc").c_str(), singletonfn.c_str());
	
	string leftreadsnfn = folder + "readsN.txt.left";
	mstcom::bsc::BSC_decompress((leftreadsnfn+".bsc").c_str(), leftreadsnfn.c_str());

	string rightreadsnfn = folder + "readsN.txt.right";
	mstcom::bsc::BSC_decompress((rightreadsnfn+".bsc").c_str(), rightreadsnfn.c_str());

	BITPOOL curdirpool;
	curdirpool.n = 8;
	int dir, file;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");
	FILE *fproot = fopen(rootstrfn.c_str(), "r");
	// sprintf(name, "rootstr.bin");
	// std::ifstream fproot(name, std::ios::binary);

	FILE *fpdecstr1 = fopen((result+"_1").c_str(), "w");
	FILE *fpdecstr2 = fopen((result+"_2").c_str(), "w");
	char str[1<<10], enstr[1<<10], dstr[1<<10];;

	// int num = 112;

	// DNAPOOL rootdnapool;
	// rootdnapool.n = 4;
	if (isorder) {
		// cout << "max_rid: " << max_rid << endl;
		uint32_t rid;
		char **readsstr = new char*[max_rid];
		bool *flag = new bool[max_rid];
		memset(flag, false, sizeof(bool)*max_rid);

		// root string
		while (fscanf(fproot, "%s", str) != EOF) {
			vector<READS_t> reads;
			READS_t r(str, 0);
			reads.push_back(r);
			int idx = 0;
			// fprintf(fpdecstr, "%s\n", str);
			fporder.read((char*)&rid, sizeof(uint32_t));
			file = getDir(fpdir, curdirpool);
			if (file) {
				// fprintf(stderr, "%lu\n", rid + (max_rid >> 1));
				readsstr[rid + (max_rid >> 1)] = strdup(str);
				flag[rid + (max_rid >> 1)] = true;
			} else {
				// fprintf(stderr, "%lu\n", rid);
				readsstr[rid] = strdup(str);
				flag[rid] = true;
			}
			// fprintf(stderr, "1111\n");

			// while(fscanf(fpenstr, "%s", enstr) != EOF) {
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					// fprintf(stderr, "xx enstr: %s\n", enstr);
					//
					decodeStr(reads[idx].str, enstr, dstr);
					//

					file = getDir(fpdir, curdirpool);
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					// fprintf(fpdecstr, "%s\n", dstr);
					fporder.read((char*)&rid, sizeof(uint32_t));
					if (file) {
						// fprintf(stderr, "%lu --\n", rid + (max_rid >> 1));
						readsstr[rid + (max_rid >> 1)] = strdup(dstr);
						flag[rid + (max_rid >> 1)] = true;
					} else {
						// fprintf(stderr, "%lu\n", rid);
						readsstr[rid] = strdup(dstr);
						flag[rid] = true;
					}
					// --num;
					// fprintf(stderr, "dstr: %s\n", dstr);
					// if (num <= 0) exit(0);

					READS_t r(dstr, idx);
					reads.push_back(r);
					idx = reads.size() - 1;

				} else { // back
					int back_step = atoi(enstr);
					while (back_step > 0) {
						idx = reads[idx].prid;
						--back_step;
					}
				}
			}
			// fprintf(stderr, "222\n");
			reads.clear();
		}

		fclose(fpenstr);
		fclose(fproot);
		// fproot.close();
		fpdir.close();
		
		// fprintf(stderr, "xxxxxxxx\n");

		FILE *fpreadsnleft = fopen(leftreadsnfn.c_str(), "r");
		while (fscanf(fpreadsnleft, "%s", str) != EOF) {
			// fprintf(fpdecstr, "%s\n", str);
			fporder.read((char*)&rid, sizeof(uint32_t));
			readsstr[rid] = strdup(str);
			flag[rid] = true;
		}
		fclose(fpreadsnleft);

		// fprintf(stderr, "yyy\n");
		// --
		FILE *fpreadsnright = fopen(rightreadsnfn.c_str(), "r");
		while (fscanf(fpreadsnright, "%s", str) != EOF) {
			// fprintf(fpdecstr, "%s\n", str);
			fporder.read((char*)&rid, sizeof(uint32_t));
			readsstr[rid + (max_rid>>1)] = strdup(str);
			flag[rid + (max_rid>>1)] = true;
		}
		fclose(fpreadsnright);
		fporder.close();

		// fprintf(stderr, "zzz\n");
		ifstream fpsingleton(singletonfn, std::ios::binary);
		DNAPOOL curdnapool;
		curdnapool.n = 4;
		int code;
		char *seq = (char*)alloca((L + 3) * sizeof(char));
		// bool flag = true;
		for (size_t id = 0; id < max_rid; ++id) {
			if (!flag[id]) {
				for (int i = 0; i < L; ++i) {
					code = getDNAcode(fpsingleton, curdnapool);
					if (code < 0) {
						break;
					}
					seq[i] = invert_code_rule[code];
				}
				seq[L] = '\0';
				// fporder.read((char*)&rid, sizeof(uint32_t));
				readsstr[id] = strdup(seq);
				// flag[rid] = true;
				// fprintf(fpdecstr, "%s\n", seq);
			}
		}
		fpsingleton.close();
		// fprintf(stderr, "www\n");

		for (size_t id = 0; id < (max_rid>>1); ++id) {
			fprintf(fpdecstr1, "%s\n", readsstr[id]);
		}

		for (size_t id = (max_rid>>1); id < max_rid; ++id) {
			fprintf(fpdecstr2, "%s\n", readsstr[id]);
		}
	} else {
		size_t leftnum = 0;
		uint32_t rid;
		char **readsstr = new char*[max_rid>>1];
		// bool *flag = new bool[max_rid];
		// memset(flag, false, sizeof(bool)*max_rid);

		// fprintf(stderr, "fdafdafdaf\n");
		// root string
		while (fscanf(fproot, "%s", str) != EOF) {
			vector<READS_t> reads;
			READS_t r(str, 0);
			reads.push_back(r);
			int idx = 0;
			// fprintf(fpdecstr, "%s\n", str);
			file = getDir(fpdir, curdirpool);
			if (file) {
				fporder.read((char*)&rid, sizeof(uint32_t));
				// fprintf(stderr, "%lu\n", rid);
				readsstr[rid] = strdup(str);
			} else {
				fprintf(fpdecstr1, "%s\n", str);
				++leftnum;
			}
			// fprintf(stderr, "1111\n");

			// while(fscanf(fpenstr, "%s", enstr) != EOF) {
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					// fprintf(stderr, "xx enstr: %s\n", enstr);
					//
					decodeStr(reads[idx].str, enstr, dstr);
					//
					file = getDir(fpdir, curdirpool);
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					// fprintf(fpdecstr, "%s\n", dstr);
					if (file) {
						fporder.read((char*)&rid, sizeof(uint32_t));
						// fprintf(stderr, "%lu --\n", rid);
						readsstr[rid] = strdup(dstr);
					} else {
						fprintf(fpdecstr1, "%s\n", dstr);
						++leftnum;
					}
					// --num;
					// fprintf(stderr, "dstr: %s\n", dstr);
					// if (num <= 0) exit(0);

					READS_t r(dstr, idx);
					reads.push_back(r);
					idx = reads.size() - 1;
				} else { // back
					int back_step = atoi(enstr);
					while (back_step > 0) {
						idx = reads[idx].prid;
						--back_step;
					}
				}
			}
			// fprintf(stderr, "222\n");
			reads.clear();
		}

		// fprintf(stderr, "xxxxxxxx\n");
		fclose(fpenstr);
		// fprintf(stderr, "xxxxxxxx\n");
		fclose(fproot);
		// fproot.close();
		fpdir.close();
		
		// fprintf(stderr, "xxxxxxxx\n");

		FILE *fpreadsnleft = fopen(leftreadsnfn.c_str(), "r");
		while (fscanf(fpreadsnleft, "%s", str) != EOF) {
			fprintf(fpdecstr1, "%s\n", str);
			++leftnum;
			// fprintf(fpdecstr1, "%s\n", dstr);
			// fporder.read((char*)&rid, sizeof(uint32_t));
			// readsstr[rid] = strdup(str);
			// flag[rid] = true;
		}
		fclose(fpreadsnleft);

		// fprintf(stderr, "yyy\n");
		// --
		FILE *fpreadsnright = fopen(rightreadsnfn.c_str(), "r");
		while (fscanf(fpreadsnright, "%s", str) != EOF) {
			// fprintf(fpdecstr, "%s\n", str);
			fporder.read((char*)&rid, sizeof(uint32_t));
			readsstr[rid] = strdup(str);
		}
		fclose(fpreadsnright);

		// fprintf(stderr, "zzz\n");
		ifstream fpsingleton(singletonfn, std::ios::binary);
		DNAPOOL curdnapool;
		curdnapool.n = 4;
		int code;
		char *seq = (char*)alloca((L + 3) * sizeof(char));
		bool flag = true;
		while (flag) {
			if (leftnum >= (max_rid >> 1)) break;
			for (int i = 0; i < L; ++i) {
				code = getDNAcode(fpsingleton, curdnapool);
				if (code < 0) {
					flag = false;
					break;
				}
				seq[i] = invert_code_rule[code];
			}
			seq[L] = '\0';
			if (!flag) break;
			// readsstr[id] = strdup(seq);
			fprintf(fpdecstr1, "%s\n", seq);
			++leftnum;
		}

		flag = true;
		while (flag) {
			for (int i = 0; flag && i < L; ++i) {
				code = getDNAcode(fpsingleton, curdnapool);
				if (code < 0) {
					flag = false;
					break;
				}
				seq[i] = invert_code_rule[code];
			}
			seq[L] = '\0';
			if (!flag) break;

			fporder.read((char*)&rid, sizeof(uint32_t));
			readsstr[rid] = strdup(seq);
			// readsstr[id] = strdup(seq);
			// fprintf(fpdecstr, "%s\n", seq);
			// ++leftnum;
		}

		fporder.close();
		fpsingleton.close();
		// fprintf(stderr, "www\n");

		for (size_t id = 0; id < (max_rid >> 1); ++id) {
			fprintf(fpdecstr2, "%s\n", readsstr[id]);
		}
	}

	fclose(fpdecstr1);
	fclose(fpdecstr2);
}

inline void build(size_t *&tr, size_t &M, size_t n) {
    for (M=1; M<=n+1; M<<=1); 
    tr = new size_t[M + n + 1];
	memset(tr, 0, sizeof(size_t)*(M + n + 1));
}

inline void up(size_t *tr, size_t x) {
	// cout << "up(x): " << x << endl;
    tr[x]=tr[x<<1]+tr[x<<1|1];
}

inline void update(size_t *tr, size_t M, size_t x, int y) {
    for(tr[x+=M]+=y,x>>=1;x;x>>=1)
        up(tr, x);
}

inline size_t query(size_t *tr, size_t M, size_t s, size_t t) {
    size_t ans=0;
    s=s+M-1;t=t+M+1;
    for(;s^t^1;s>>=1,t>>=1) {
        if(~s&1) {
        	ans+=tr[s^1];
        }
        if(t&1) {
        	ans+=tr[t^1];
        }
    }
    return ans;
}

inline size_t binarySearch(size_t *tr, size_t M, size_t s, size_t max_rid, bool *flag, size_t dist) { // find a 
	size_t m, num, e = max_rid, oe = max_rid, tempe, rdist;
	while (e > s) {
		num = query(tr, M, s, e);
		tempe = e;
		if (e - s + 1 >= num) {
			rdist = e - s + 1 - num;
			if (dist < rdist) {
				// oe = e;
				e = s + ((e - s) >> 1) - 1;
				if (e == s) {
					++ e;
				} else 
				if (e < s) {
					++ e;
					num = query(tr, M, s, e);
					if (dist == e - s + 1 - num) {
						break;
					}
				}
			} else 
			if (dist > rdist) {
				e += ((oe - e) >> 1);
			} else {
				while (e >= s && flag[e - 1]) {
					-- e;
				} 
				num = query(tr, M, s, e);
				if (!(e - s + 1 >= num && e - s + 1 - num == dist)) {
					++ e;
				}
				break;
			}
			if (tempe > e) oe = tempe;
		} else {
			e += ((oe - e) >> 1);
		}
	}

	while (e <= max_rid && flag[e - 1]) {
		++ e;
	}
	return e;
}

void decompressPEX(string result) {
	fprintf(stderr, "decompressPE()...\n");
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string rootstrfn = folder + "rootstr.txt";
	mstcom::bsc::BSC_decompress((rootstrfn+".bsc").c_str(), rootstrfn.c_str());

	string orderfn = folder + string("order.bin");
	std::ifstream fporder;
	if (isorder) {
		mstcom::bsc::BSC_decompress((orderfn+".bsc").c_str(), orderfn.c_str());
		fporder.open(orderfn, std::ios::binary);
	}
	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string singletonfn = folder + "singleton.bin"; 
	mstcom::bsc::BSC_decompress((singletonfn+".bsc").c_str(), singletonfn.c_str());
	
	string readsnfn = folder + "readsN.txt";
	mstcom::bsc::BSC_decompress((readsnfn+".bsc").c_str(), readsnfn.c_str());

	string distfn = folder + "dist.bin";
	mstcom::bsc::BSC_decompress((distfn+".bsc").c_str(), distfn.c_str());

	BITPOOL curdirpool;
	curdirpool.n = 8;
	int dir, file;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");
	FILE *fproot = fopen(rootstrfn.c_str(), "r");
	// sprintf(name, "rootstr.bin");
	// std::ifstream fproot(name, std::ios::binary);

	FILE *fpdecstr1 = fopen((result+"_1").c_str(), "w");
	FILE *fpdecstr2 = fopen((result+"_2").c_str(), "w");
	char str[1<<10], enstr[1<<10], dstr[1<<10];;

	// int num = 112;

	if (isorder) {
		cout << "max_rid: " << max_rid << endl;
		uint32_t rid = 0;
		char **readsstr = new char*[max_rid];
		
		// root string
		while (fscanf(fproot, "%s", str) != EOF) {
			vector<READS_t> reads;
			READS_t r(str, 0);
			reads.push_back(r);
			int idx = 0;
			
			readsstr[rid ++] = strdup(str);
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					decodeStr(reads[idx].str, enstr, dstr);
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					
					readsstr[rid++] = strdup(dstr);

					READS_t r(dstr, idx);
					reads.push_back(r);
					idx = reads.size() - 1;

				} else { // back
					int back_step = atoi(enstr);
					while (back_step > 0) {
						idx = reads[idx].prid;
						--back_step;
					}
				}
			}
			reads.clear();
		}

		fclose(fpenstr);
		fclose(fproot);
		// fproot.close();
		fpdir.close();
		
		FILE *fpreadsn = fopen(readsnfn.c_str(), "r");
		while (fscanf(fpreadsn, "%s", str) != EOF) {
			readsstr[rid ++] = strdup(str);
		}
		fclose(fpreadsn);

		ifstream fpsingleton(singletonfn, std::ios::binary);
		DNAPOOL curdnapool;
		curdnapool.n = 4;
		int code;
		char *seq = (char*)alloca((L + 3) * sizeof(char));
		// bool flag = true;
		while (rid < max_rid) {
			for (int i = 0; i < L; ++i) {
				code = getDNAcode(fpsingleton, curdnapool);
				if (code < 0) {
					break;
				}
				seq[i] = invert_code_rule[code];
			}
			seq[L] = '\0';
			readsstr[rid++] = strdup(seq);
		}
		fpsingleton.close();

		size_t *ids = new size_t[max_rid + 3];

		size_t *tr, M;
		build(tr, M, max_rid);

		bool *flag = new bool[max_rid + 3];
		memset(flag, false, sizeof(bool)*max_rid); 

		size_t half = max_rid >> 1, pid;
		std::ifstream fpdist(distfn.c_str(), std::ios::binary);
		
		uint32_t tempid;
		int32_t dist;

		for (size_t num = 0; num < max_rid; ++num) {
			if (!flag[num]) {
				fporder.read((char*)&tempid, sizeof(uint32_t));
				// fscanf(fpdist, "%d", &dist);
				fpdist.read((char*)&dist, sizeof(int32_t));

				ids[num] = tempid;
				if (dist < 0) {
					dist = 0 - dist;
					ids[num] += half;
				}

				pid = binarySearch(tr, M, num + 1 + 1, max_rid, flag, dist);

				if (ids[num] < half) {
					ids[pid - 1] = half + ids[num];
				} else {
					ids[pid - 1] = ids[num] - half;
				}

				update(tr, M, pid, 1);
				flag[pid - 1] = true;
			}
		}

		fpdist.close();
		fporder.close();
		delete[] tr;
		delete[] flag;

		size_t *od = new size_t[max_rid + 3];
		for (size_t id = 0; id < max_rid; ++id) {
			od[ids[id]] = id;
		}
		delete[] ids;

		for (size_t id = 0; id < half; ++id) {
			fprintf(fpdecstr1, "%s\n", readsstr[od[id]]);
		}

		for (size_t id = half; id < max_rid; ++id) {
			fprintf(fpdecstr2, "%s\n", readsstr[od[id]]);
		}
	} else {
		size_t leftnum = 0;
		uint32_t rid = 0;
		char **readsstr = new char*[max_rid + 1];

		while (fscanf(fproot, "%s", str) != EOF) {
			vector<READS_t> reads;
			READS_t r(str, 0);
			reads.push_back(r);
			int idx = 0;
			
			readsstr[rid ++] = strdup(str);
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					decodeStr(reads[idx].str, enstr, dstr);
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					readsstr[rid ++] = strdup(dstr);
					
					READS_t r(dstr, idx);
					reads.push_back(r);
					idx = reads.size() - 1;
				} else { // back
					int back_step = atoi(enstr);
					while (back_step > 0) {
						idx = reads[idx].prid;
						--back_step;
					}
				}
			}
			reads.clear();
		}

		fclose(fpenstr);
		fclose(fproot);
		fpdir.close();
		
		FILE *fpreadsn = fopen(readsnfn.c_str(), "r");
		while (fscanf(fpreadsn, "%s", str) != EOF) {
			readsstr[rid ++] = strdup(str);
		}
		fclose(fpreadsn);

		ifstream fpsingleton(singletonfn, std::ios::binary);
		DNAPOOL curdnapool;
		curdnapool.n = 4;
		int code;
		char *seq = (char*)alloca((L + 3) * sizeof(char));
		while (rid < max_rid) {
			for (int i = 0; i < L; ++i) {
				code = getDNAcode(fpsingleton, curdnapool);
				if (code < 0) {
					break;
				}
				seq[i] = invert_code_rule[code];
			}
			seq[L] = '\0';
			readsstr[rid ++] = strdup(seq);
		}
		fpsingleton.close();

		// cout << "xxxx\n";
		///
		size_t *ids = new size_t[max_rid + 3];

		size_t *tr, M;
		build(tr, M, max_rid);

		bool *flag = new bool[max_rid + 3];
		memset(flag, false, sizeof(bool)*max_rid); 

		size_t half = max_rid >> 1, pid;
		std::ifstream fpdist(distfn.c_str(), std::ios::binary);
		
		uint32_t tempid;
		int32_t dist;
		size_t id = 0;

		for (size_t num = 0; num < max_rid; ++num) {
			if (!flag[num]) {
				fpdist.read((char*)&dist, sizeof(int32_t));

				ids[num] = id ++;
				if (dist < 0) {
					dist = 0 - dist;
					ids[num] += half;
				}

				pid = binarySearch(tr, M, num + 1 + 1, max_rid, flag, dist);

				if (ids[num] < half) {
					ids[pid - 1] = half + ids[num];
				} else {
					ids[pid - 1] = ids[num] - half;
				}

				update(tr, M, pid, 1);
				flag[pid - 1] = true;
			}
		}

		fpdist.close();
		fporder.close();
		delete[] tr;
		delete[] flag;

		size_t *od = new size_t[max_rid + 3];
		for (size_t id = 0; id < max_rid; ++id) {
			od[ids[id]] = id;
		}
		delete[] ids;

		for (size_t id = 0; id < half; ++id) {
			fprintf(fpdecstr1, "%s\n", readsstr[od[id]]);
		}

		for (size_t id = half; id < max_rid; ++id) {
			fprintf(fpdecstr2, "%s\n", readsstr[od[id]]);
		}
	}

	fclose(fpdecstr1);
	fclose(fpdecstr2);
}

int main(int argc, char const *argv[]) {
	stopwatch.start();

	folder = generateString("mstcomdec", 10); //creat a temp folder for current input

	string cmd = "mkdir -p " + folder;
	fprintf(stderr, "%s\n", folder.c_str());
	system(cmd.c_str());
	folder += "/";

	init();

	cmd = "tar -xf " + string(argv[1]) + " -C " + folder;
	// tar -xf $filename -C $decomp
	system(cmd.c_str());

	string parfn = folder + "par.txt";
	cout << parfn << endl;
	FILE *fp = fopen(parfn.c_str(), "r");
	fscanf(fp, "%d%d%d", &L, &ispe, &isorder);
	fprintf(stderr, "L: %d;\nispe: %d\nisorder: %d\n", L, ispe, isorder);
	if (isorder || ispe) {
		fscanf(fp, "%lu", &max_rid);
		fprintf(stderr, "max_rid: %lu\n", max_rid);
	}
	fclose(fp);

	if (ispe) {
		// decompressPE(argv[2]);
		decompressPEX(argv[2]);
	} else {
		decompressSingle(argv[2]);
	}

	cmd = "rm -rf " + folder;
	system(cmd.c_str());

	cout << "Time of decompressioin = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	return 0;
}

CStopWatch::CStopWatch(){
	running = 0;
}


CStopWatch::~CStopWatch(){
}


void CStopWatch::start(){
	running = 1;
	elapsed.clear();
	t1 = std::chrono::steady_clock::now();
}


double CStopWatch::stop(){
	t2  = std::chrono::steady_clock::now();
	running = 0;
	double el = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0;
	elapsed.push_back(el);
	return el;
}


void CStopWatch::resume(){
	running = 1;
	t1 = std::chrono::steady_clock::now();
}


double CStopWatch::totalTime(){
	double sum = 0;
	for (std::vector<double>::iterator it = elapsed.begin(); it != elapsed.end(); ++it) {
		sum += *it;
	}
	return sum;
}
// * History
// * ...
// * Dec 8, 2018. decompression
// * ... *//
# include <cstdio>
# include <iostream>
# include <fstream>
# include <vector>
# include <map>
# include <queue>
# include <string>
# include <cstring>
# include <chrono>
# include <thread>
# include <mutex>
# include <cstdlib>
# include "libbsc/bsc.h"
# include "mstcom.h"
using namespace std;

int L;
int ispe, isorder;
uint32_t max_rid;
const int nthreads = 24;
string folder;

struct READS_t {
	string str;
	uint32_t rid, prid; //reads id of parent
	READS_t() {}
	// READS_t(char *_str, const uint32_t &_prid): str(_str), prid(_prid) {}
	READS_t(char *_str, const uint32_t &_rid, const uint32_t &_prid): str(_str), rid(_rid), prid(_prid) {}
};

char complement[256];

const char invert_code_rule[4] = {'A', 'C', 'G', 'T'}; //decoding rule //A - 0; C - 1; G - 2; T - 3;

class CStopWatch {
	private:
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	std::chrono::steady_clock::time_point t_temp;
	int running;
	std::vector<double> elapsed;

	public:
	CStopWatch();
	~CStopWatch();
	void start();
	double stop();
	void resume();
	double totalTime();
};

CStopWatch stopwatch;

void init() {
	memset(complement, 125, 256);
	complement['A'] = 'T';
	complement['a'] = 'T';
	complement['C'] = 'G';
	complement['c'] = 'G';
	complement['G'] = 'C';
	complement['g'] = 'C';
	complement['T'] = 'A';
	complement['t'] = 'A';
	complement['N'] = 'N';
	complement['n'] = 'N';
}

inline void reverseComplement(char* start) {
	char* left = start; // sequence starts
	char* right = start + L - 1;
	while (right > left) {
		char tmp = complement[(uint8_t)*left];
		*left = complement[(uint8_t)*right];
		*right = tmp;
		++left;
		--right;
	}
	if (left == right)
		*left = complement[(uint8_t)*left];
}

inline void decodeStr(string pstr, char *enstr, char *dstr) {
	bool debug = false;
	if (enstr[0] == '\n') { // equal
		strcpy(dstr, pstr.c_str());
	} else {
		// check space ' '
		int spidx = -1;
		int len = strlen(enstr) - 1;
		for (int i = 0; i < len; ++i) {
			if (enstr[i] == ' ') {
				spidx = i;
				break;
			}
		}
		int dstridx = 0;

		bool shiftRight = false;
		if (spidx == 0) {
			shiftRight = true;
		} else 
		if (spidx > -1) {
			for (int i = 0; i < spidx; ++i) {
				if (enstr[i] >= '0' && enstr[i] <= '9') {
					shiftRight = true;
					break;
				}
			}
		}

		if (spidx > -1 && spidx != len - 1 && !shiftRight) {
			bool flag = true;
			for (int i = spidx + 1; i < len; ++i) {
				if (enstr[i] >= '0' && enstr[i] <= '9') {
					flag = false;
					break;
				}
			}
			if (flag) {
				shiftRight = true;
			}
		}

		int pstridx = 0, enstridx = 0, maxenstridx = len;
		if (spidx > -1) {
			if (shiftRight) {
				maxenstridx = spidx;
				pstridx = len - spidx - 1;
			} else {
				for (enstridx = 0; enstridx < spidx; ++enstridx) {
					dstr[dstridx++] = enstr[enstridx];
				}
				++ enstridx;
				// if (spidx == len - 1) {
				// 	pstridx = spidx - 1;
				// }
			}
		}

		if (debug) {
			fprintf(stderr, "pstr: %s\n", pstr.c_str());
			fprintf(stderr, "spidx: %d\n", spidx);
			fprintf(stderr, "shiftRight: %d\n", shiftRight);
			fprintf(stderr, "len: %d\n", len);
			fprintf(stderr, "maxenstridx: %d\n", maxenstridx);
			fprintf(stderr, "dstridx: %d\n", dstridx);
			fprintf(stderr, "pstridx: %d\n", pstridx);
		}

		int eq_char_num = 0;
		for (; enstridx < maxenstridx; ++enstridx) {
			if (enstr[enstridx] >= 'A' && enstr[enstridx] <= 'Z') {
				if (eq_char_num > 0) {
					for (int j = 0; j < eq_char_num; ++j) {
						dstr[dstridx++] = pstr[pstridx++];
					}
					eq_char_num = 0;
				}
				dstr[dstridx++] = enstr[enstridx];
				if (debug) {
					fprintf(stderr, "dstr[%d]: %c\n", dstridx-1, dstr[dstridx-1]);
					fprintf(stderr, "enstr[%d]: %c\n", enstridx, enstr[enstridx]);
				}
				pstridx++;
			} else {
				eq_char_num = eq_char_num * 10 + enstr[enstridx] - '0';
			}
		}
		dstr[dstridx] = '\0';
		if (debug) {
			fprintf(stderr, "%s\n", dstr);
			fprintf(stderr, "maxenstridx: %d\n", maxenstridx);
		}
		for (; pstridx < L && dstridx < L; ++pstridx) {
			dstr[dstridx++] = pstr[pstridx];
		}
		if (shiftRight) {
			for (int i = maxenstridx + 1; i < len && dstridx < L; ++i) {
				dstr[dstridx++] = enstr[i];
			}
		} 
		dstr[L] = '\0';
		if (debug) {
			fprintf(stderr, "%s\n", dstr);
			// exit(1);
		}
	}
}

struct BITPOOL {
	int n, a[8];
};

inline int getDir(std::ifstream& fpdir, BITPOOL& curdirpool) {
	if (curdirpool.n >= 8) {
		uint8_t dirbin;
		fpdir.read((char*)&dirbin, sizeof(uint8_t));
		for (int i = 0; i < 8; ++i) {
		// for (int i = 7; i >= 0; --i) {
			curdirpool.a[i] = dirbin&1;
			dirbin >>= 1;
		}
		curdirpool.n = 0;
	}
	int res = curdirpool.a[curdirpool.n];
	++curdirpool.n;
	return res;
}

struct DNAPOOL {
	int n, a[4];
};

inline int getDNAcode(std::ifstream& fpif, DNAPOOL& curdnapool) {
	// fprintf(stderr, "yyyyyy\n");
	if (curdnapool.n >= 4) {
		uint8_t bin;
		// fprintf(stderr, "zzzz\n");
		fpif.read((char*)&bin, sizeof(uint8_t));
		// fprintf(stderr, "%d\n", bin);
		if (fpif.eof()) return -1;

		for (int i = 0; i < 4; ++i) {
			curdnapool.a[i] = bin&3;
			bin >>= 2;
		}
		curdnapool.n = 0;
	}
	int res = curdnapool.a[curdnapool.n];
	++curdnapool.n;
	return res;
}

inline bool getReads(std::ifstream& fp, DNAPOOL& curdnapool, char *r) {
	int code, i;
	for (i = 0; i < L; ++i) {
		code = getDNAcode(fp, curdnapool);
		if (code < 0) break;
		r[i] = invert_code_rule[code];
	}
	r[i] = '\0';
	if (i == L) return true;
	return false;
}

bool checkStrOrDig(char *str) {
	int len = strlen(str);
	// if (str[0] == ' ' && str[1] == '\n') { // ' \n'
	if (str[0] == '\n') { // ' \n'
		return true;
	}
	for (int i = 0; i < len - 1; ++i) {
		if (str[i] >= 'A' && str[i] <= 'Z') {
			return true;
		}
	}
	return false;
}

int mainxxx(int argc, char const *argv[]) {
	string parent = "AAAGGCCCCAGTTTGGCAGACCGACAGCGTGAATACCTTTTAGACATGATCCCTCCCCGGTCTATATCGCAGTCCATCAGTGGACAGAAATAACGCCTGTT";
	char enstr[1<<10], dstr[1<<10];
	strcpy(enstr, "C AT0\n");
	decodeStr(parent, enstr, dstr);
	fprintf(stderr, "%s\n", dstr);
	return 0;
}

static const char alphanum[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
inline std::string generateString(const std::string &chr, int length = 5) {
	srand(time(0));
	string res = chr + "_";
	for (int i = 0; i < length; ++i) {
		res += alphanum[rand() % 62];
	}
	return res;
}

#ifdef false
void decompressSingle(string result) {
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string rootstrfn = folder + "rootstr.txt";
	mstcom::bsc::BSC_decompress((rootstrfn+".bsc").c_str(), rootstrfn.c_str());
	string orderfn;

	std::ifstream fporder;
	if (isorder) {
		orderfn = folder + string("order.bin");
		mstcom::bsc::BSC_decompress((orderfn+".bsc").c_str(), orderfn.c_str());
		fporder.open(orderfn, std::ios::binary);
	}
	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	string singletonfn = folder + "singleton.bin"; 
	mstcom::bsc::BSC_decompress((singletonfn+".bsc").c_str(), singletonfn.c_str());
	string readsnfn = folder + "readsN.txt";
	mstcom::bsc::BSC_decompress((readsnfn+".bsc").c_str(), readsnfn.c_str());

	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	BITPOOL curdirpool;
	curdirpool.n = 8;
	int dir;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");
	FILE *fproot = fopen(rootstrfn.c_str(), "r");
	// sprintf(name, "rootstr.bin");
	// std::ifstream fproot(name, std::ios::binary);

	FILE *fpdecstr = fopen(result.c_str(), "w");
	char str[1<<10], enstr[1<<10], dstr[1<<10];;

	// int num = 112;

	// DNAPOOL rootdnapool;
	// rootdnapool.n = 4;
	if (isorder) {
		cout << "max_rid: " << max_rid << endl;
		uint32_t rid;
		char **readsstr = new char*[max_rid];
		bool *flag = new bool[max_rid];
		memset(flag, false, sizeof(bool)*max_rid);

		while (fscanf(fproot, "%s", str) != EOF) {
		// while (getReads(fproot, rootdnapool, str)) {
			vector<READS_t> reads;
			READS_t r(str, 0);
			reads.push_back(r);
			int idx = 0;
			// fprintf(fpdecstr, "%s\n", str);
			fporder.read((char*)&rid, sizeof(uint32_t));
			readsstr[rid] = strdup(str);
			flag[rid] = true;

			// while(fscanf(fpenstr, "%s", enstr) != EOF) {
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					// fprintf(stderr, "xx enstr: %s\n", enstr);
					//
					decodeStr(reads[idx].str, enstr, dstr);
					//
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					// fprintf(fpdecstr, "%s\n", dstr);
					fporder.read((char*)&rid, sizeof(uint32_t));
					readsstr[rid] = strdup(dstr);
					flag[rid] = true;
					// --num;
					// fprintf(stderr, "dstr: %s\n", dstr);
					// if (num <= 0) exit(0);

					READS_t r(dstr, idx);
					reads.push_back(r);
					idx = reads.size() - 1;

				} else { // back
					int back_step = atoi(enstr);
					while (back_step > 0) {
						idx = reads[idx].prid;
						--back_step;
					}
				}
			}
			reads.clear();
		}

		fclose(fpenstr);
		fclose(fproot);
		// fproot.close();
		fpdir.close();
		
		FILE *fpreadsn = fopen(readsnfn.c_str(), "r");
		while (fscanf(fpreadsn, "%s", str) != EOF) {
			// fprintf(fpdecstr, "%s\n", str);
			fporder.read((char*)&rid, sizeof(uint32_t));
			readsstr[rid] = strdup(str);
			flag[rid] = true;
		}
		fclose(fpreadsn);

		ifstream fpsingleton(singletonfn, std::ios::binary);
		DNAPOOL curdnapool;
		curdnapool.n = 4;
		int code;
		char *seq = (char*)alloca((L + 3) * sizeof(char));
		// bool flag = true;
		for (uint32_t id = 0; id < max_rid; ++id) {
			if (!flag[id]) {
				for (int i = 0; i < L; ++i) {
					code = getDNAcode(fpsingleton, curdnapool);
					if (code < 0) {
						break;
					}
					seq[i] = invert_code_rule[code];
				}
				seq[L] = '\0';
				// fporder.read((char*)&rid, sizeof(uint32_t));
				readsstr[id] = strdup(seq);
				// flag[rid] = true;
				// fprintf(fpdecstr, "%s\n", seq);
			}
		}
		fporder.close();
		fpsingleton.close();
		for (uint32_t id = 0; id < max_rid; ++id) {
			fprintf(fpdecstr, "%s\n", readsstr[id]);
		}

	} else {
		while (fscanf(fproot, "%s", str) != EOF) {
			vector<READS_t> reads;
			READS_t r(str, 0);
			reads.push_back(r);
			int idx = 0;
			fprintf(fpdecstr, "%s\n", str);

			// while(fscanf(fpenstr, "%s", enstr) != EOF) {
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					// fprintf(stderr, "xx enstr: %s\n", enstr);
					//
					decodeStr(reads[idx].str, enstr, dstr);
					//
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					fprintf(fpdecstr, "%s\n", dstr);
					// --num;
					// fprintf(stderr, "dstr: %s\n", dstr);
					// if (num <= 0) exit(0);

					READS_t r(dstr, idx);
					reads.push_back(r);
					idx = reads.size() - 1;

				} else { // back
					int back_step = atoi(enstr);
					while (back_step > 0) {
						idx = reads[idx].prid;
						--back_step;
					}
				}
			}
			reads.clear();
		}

		fclose(fpenstr);
		fclose(fproot);
		// fproot.close();
		fpdir.close();
		ifstream fpsingleton(singletonfn, std::ios::binary);
		DNAPOOL curdnapool;
		curdnapool.n = 4;
		int code;
		char *seq = (char*)alloca((L + 3) * sizeof(char));
		bool flag = true;
		while (flag) {
			for (int i = 0; flag && i < L; ++i) {
				code = getDNAcode(fpsingleton, curdnapool);
				if (code < 0) {
					flag = false;
					break;
				}
				seq[i] = invert_code_rule[code];
			}
			if (!flag) break;
			seq[L] = '\0';
			fprintf(fpdecstr, "%s\n", seq);
		}
		fpsingleton.close();

		FILE *fpreadsn = fopen(readsnfn.c_str(), "r");
		while (fscanf(fpreadsn, "%s", str) != EOF) {
			fprintf(fpdecstr, "%s\n", str);
		}
		fclose(fpreadsn);
	}

	fclose(fpdecstr);
}

void decompressPE(string result) {
	fprintf(stderr, "decompressPE()...\n");
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string rootstrfn = folder + "rootstr.txt";
	mstcom::bsc::BSC_decompress((rootstrfn+".bsc").c_str(), rootstrfn.c_str());

	string orderfn = folder + string("order.bin");
	std::ifstream fporder;
	mstcom::bsc::BSC_decompress((orderfn+".bsc").c_str(), orderfn.c_str());
	fporder.open(orderfn, std::ios::binary);

	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string singletonfn = folder + "singleton.bin"; 
	mstcom::bsc::BSC_decompress((singletonfn+".bsc").c_str(), singletonfn.c_str());
	
	string leftreadsnfn = folder + "readsN.txt.left";
	mstcom::bsc::BSC_decompress((leftreadsnfn+".bsc").c_str(), leftreadsnfn.c_str());

	string rightreadsnfn = folder + "readsN.txt.right";
	mstcom::bsc::BSC_decompress((rightreadsnfn+".bsc").c_str(), rightreadsnfn.c_str());

	BITPOOL curdirpool;
	curdirpool.n = 8;
	int dir, file;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");
	FILE *fproot = fopen(rootstrfn.c_str(), "r");
	// sprintf(name, "rootstr.bin");
	// std::ifstream fproot(name, std::ios::binary);

	FILE *fpdecstr1 = fopen((result+"_1").c_str(), "w");
	FILE *fpdecstr2 = fopen((result+"_2").c_str(), "w");
	char str[1<<10], enstr[1<<10], dstr[1<<10];;

	// int num = 112;

	// DNAPOOL rootdnapool;
	// rootdnapool.n = 4;
	if (isorder) {
		// cout << "max_rid: " << max_rid << endl;
		uint32_t rid;
		char **readsstr = new char*[max_rid];
		bool *flag = new bool[max_rid];
		memset(flag, false, sizeof(bool)*max_rid);

		// root string
		while (fscanf(fproot, "%s", str) != EOF) {
			vector<READS_t> reads;
			READS_t r(str, 0);
			reads.push_back(r);
			int idx = 0;
			// fprintf(fpdecstr, "%s\n", str);
			fporder.read((char*)&rid, sizeof(uint32_t));
			file = getDir(fpdir, curdirpool);
			if (file) {
				// fprintf(stderr, "%lu\n", rid + (max_rid >> 1));
				readsstr[rid + (max_rid >> 1)] = strdup(str);
				flag[rid + (max_rid >> 1)] = true;
			} else {
				// fprintf(stderr, "%lu\n", rid);
				readsstr[rid] = strdup(str);
				flag[rid] = true;
			}
			// fprintf(stderr, "1111\n");

			// while(fscanf(fpenstr, "%s", enstr) != EOF) {
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					// fprintf(stderr, "xx enstr: %s\n", enstr);
					//
					decodeStr(reads[idx].str, enstr, dstr);
					//

					file = getDir(fpdir, curdirpool);
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					// fprintf(fpdecstr, "%s\n", dstr);
					fporder.read((char*)&rid, sizeof(uint32_t));
					if (file) {
						// fprintf(stderr, "%lu --\n", rid + (max_rid >> 1));
						readsstr[rid + (max_rid >> 1)] = strdup(dstr);
						flag[rid + (max_rid >> 1)] = true;
					} else {
						// fprintf(stderr, "%lu\n", rid);
						readsstr[rid] = strdup(dstr);
						flag[rid] = true;
					}
					// --num;
					// fprintf(stderr, "dstr: %s\n", dstr);
					// if (num <= 0) exit(0);

					READS_t r(dstr, idx);
					reads.push_back(r);
					idx = reads.size() - 1;

				} else { // back
					int back_step = atoi(enstr);
					while (back_step > 0) {
						idx = reads[idx].prid;
						--back_step;
					}
				}
			}
			// fprintf(stderr, "222\n");
			reads.clear();
		}

		fclose(fpenstr);
		fclose(fproot);
		// fproot.close();
		fpdir.close();
		
		// fprintf(stderr, "xxxxxxxx\n");

		FILE *fpreadsnleft = fopen(leftreadsnfn.c_str(), "r");
		while (fscanf(fpreadsnleft, "%s", str) != EOF) {
			// fprintf(fpdecstr, "%s\n", str);
			fporder.read((char*)&rid, sizeof(uint32_t));
			readsstr[rid] = strdup(str);
			flag[rid] = true;
		}
		fclose(fpreadsnleft);

		// fprintf(stderr, "yyy\n");
		// --
		FILE *fpreadsnright = fopen(rightreadsnfn.c_str(), "r");
		while (fscanf(fpreadsnright, "%s", str) != EOF) {
			// fprintf(fpdecstr, "%s\n", str);
			fporder.read((char*)&rid, sizeof(uint32_t));
			readsstr[rid + (max_rid>>1)] = strdup(str);
			flag[rid + (max_rid>>1)] = true;
		}
		fclose(fpreadsnright);
		fporder.close();

		// fprintf(stderr, "zzz\n");
		ifstream fpsingleton(singletonfn, std::ios::binary);
		DNAPOOL curdnapool;
		curdnapool.n = 4;
		int code;
		char *seq = (char*)alloca((L + 3) * sizeof(char));
		// bool flag = true;
		for (size_t id = 0; id < max_rid; ++id) {
			if (!flag[id]) {
				for (int i = 0; i < L; ++i) {
					code = getDNAcode(fpsingleton, curdnapool);
					if (code < 0) {
						break;
					}
					seq[i] = invert_code_rule[code];
				}
				seq[L] = '\0';
				// fporder.read((char*)&rid, sizeof(uint32_t));
				readsstr[id] = strdup(seq);
				// flag[rid] = true;
				// fprintf(fpdecstr, "%s\n", seq);
			}
		}
		fpsingleton.close();
		// fprintf(stderr, "www\n");

		for (size_t id = 0; id < (max_rid>>1); ++id) {
			fprintf(fpdecstr1, "%s\n", readsstr[id]);
		}

		for (size_t id = (max_rid>>1); id < max_rid; ++id) {
			fprintf(fpdecstr2, "%s\n", readsstr[id]);
		}
	} else {
		size_t leftnum = 0;
		uint32_t rid;
		char **readsstr = new char*[max_rid>>1];
		// bool *flag = new bool[max_rid];
		// memset(flag, false, sizeof(bool)*max_rid);

		// fprintf(stderr, "fdafdafdaf\n");
		// root string
		while (fscanf(fproot, "%s", str) != EOF) {
			vector<READS_t> reads;
			READS_t r(str, 0);
			reads.push_back(r);
			int idx = 0;
			// fprintf(fpdecstr, "%s\n", str);
			file = getDir(fpdir, curdirpool);
			if (file) {
				fporder.read((char*)&rid, sizeof(uint32_t));
				// fprintf(stderr, "%lu\n", rid);
				readsstr[rid] = strdup(str);
			} else {
				fprintf(fpdecstr1, "%s\n", str);
				++leftnum;
			}
			// fprintf(stderr, "1111\n");

			// while(fscanf(fpenstr, "%s", enstr) != EOF) {
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					// fprintf(stderr, "xx enstr: %s\n", enstr);
					//
					decodeStr(reads[idx].str, enstr, dstr);
					//
					file = getDir(fpdir, curdirpool);
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					// fprintf(fpdecstr, "%s\n", dstr);
					if (file) {
						fporder.read((char*)&rid, sizeof(uint32_t));
						// fprintf(stderr, "%lu --\n", rid);
						readsstr[rid] = strdup(dstr);
					} else {
						fprintf(fpdecstr1, "%s\n", dstr);
						++leftnum;
					}
					// --num;
					// fprintf(stderr, "dstr: %s\n", dstr);
					// if (num <= 0) exit(0);

					READS_t r(dstr, idx);
					reads.push_back(r);
					idx = reads.size() - 1;
				} else { // back
					int back_step = atoi(enstr);
					while (back_step > 0) {
						idx = reads[idx].prid;
						--back_step;
					}
				}
			}
			// fprintf(stderr, "222\n");
			reads.clear();
		}

		// fprintf(stderr, "xxxxxxxx\n");
		fclose(fpenstr);
		// fprintf(stderr, "xxxxxxxx\n");
		fclose(fproot);
		// fproot.close();
		fpdir.close();
		
		// fprintf(stderr, "xxxxxxxx\n");

		FILE *fpreadsnleft = fopen(leftreadsnfn.c_str(), "r");
		while (fscanf(fpreadsnleft, "%s", str) != EOF) {
			fprintf(fpdecstr1, "%s\n", str);
			++leftnum;
			// fprintf(fpdecstr1, "%s\n", dstr);
			// fporder.read((char*)&rid, sizeof(uint32_t));
			// readsstr[rid] = strdup(str);
			// flag[rid] = true;
		}
		fclose(fpreadsnleft);

		// fprintf(stderr, "yyy\n");
		// --
		FILE *fpreadsnright = fopen(rightreadsnfn.c_str(), "r");
		while (fscanf(fpreadsnright, "%s", str) != EOF) {
			// fprintf(fpdecstr, "%s\n", str);
			fporder.read((char*)&rid, sizeof(uint32_t));
			readsstr[rid] = strdup(str);
		}
		fclose(fpreadsnright);

		// fprintf(stderr, "zzz\n");
		ifstream fpsingleton(singletonfn, std::ios::binary);
		DNAPOOL curdnapool;
		curdnapool.n = 4;
		int code;
		char *seq = (char*)alloca((L + 3) * sizeof(char));
		bool flag = true;
		while (flag) {
			if (leftnum >= (max_rid >> 1)) break;
			for (int i = 0; i < L; ++i) {
				code = getDNAcode(fpsingleton, curdnapool);
				if (code < 0) {
					flag = false;
					break;
				}
				seq[i] = invert_code_rule[code];
			}
			seq[L] = '\0';
			if (!flag) break;
			// readsstr[id] = strdup(seq);
			fprintf(fpdecstr1, "%s\n", seq);
			++leftnum;
		}

		flag = true;
		while (flag) {
			for (int i = 0; flag && i < L; ++i) {
				code = getDNAcode(fpsingleton, curdnapool);
				if (code < 0) {
					flag = false;
					break;
				}
				seq[i] = invert_code_rule[code];
			}
			seq[L] = '\0';
			if (!flag) break;

			fporder.read((char*)&rid, sizeof(uint32_t));
			readsstr[rid] = strdup(seq);
			// readsstr[id] = strdup(seq);
			// fprintf(fpdecstr, "%s\n", seq);
			// ++leftnum;
		}

		fporder.close();
		fpsingleton.close();
		// fprintf(stderr, "www\n");

		for (size_t id = 0; id < (max_rid >> 1); ++id) {
			fprintf(fpdecstr2, "%s\n", readsstr[id]);
		}
	}

	fclose(fpdecstr1);
	fclose(fpdecstr2);
}
#endif

inline void build(uint32_t *&tr, uint32_t &M, uint32_t n) {
    for (M=1; M<=n+1; M<<=1); 
    tr = new uint32_t[M + n + 1];
	memset(tr, 0, sizeof(uint32_t)*(M + n + 1));
}

inline void up(uint32_t *tr, uint32_t x) {
	// cout << "up(x): " << x << endl;
    tr[x]=tr[x<<1]+tr[x<<1|1];
}

inline void update(uint32_t *tr, uint32_t M, uint32_t x, int y) {
    for(tr[x+=M]+=y,x>>=1;x;x>>=1)
        up(tr, x);
}

inline uint32_t query(uint32_t *tr, uint32_t M, uint32_t s, uint32_t t) {
	if (s >= t) return 0;
    uint32_t ans=0;
    s=s+M-1;t=t+M+1;
    for(;s^t^1;s>>=1,t>>=1) {
        if(~s&1) {
        	ans+=tr[s^1];
        }
        if(t&1) {
        	ans+=tr[t^1];
        }
    }
    return ans;
}

inline uint32_t binarySearch(uint32_t *tr, uint32_t M, uint32_t s, uint32_t max_rid, bool *flag, uint32_t dist) { // find a 
	uint32_t m, num, e = max_rid, oe = max_rid, tempe, rdist;
	while (e > s) {
		num = query(tr, M, s, e);
		tempe = e;
		if (e - s + 1 >= num) {
			rdist = e - s + 1 - num;
			if (dist < rdist) {
				// oe = e;
				e = s + ((e - s) >> 1) - 1;
				if (e == s) {
					++ e;
				} else 
				if (e < s) {
					++ e;
					num = query(tr, M, s, e);
					if (dist == e - s + 1 - num) {
						break;
					}
				}
			} else 
			if (dist > rdist) {
				e += ((oe - e) >> 1);
			} else {
				while (e >= s && flag[e - 1]) {
					-- e;
				} 
				num = query(tr, M, s, e);
				if (!(e - s + 1 >= num && e - s + 1 - num == dist)) {
					++ e;
				}
				break;
			}
			if (tempe > e) oe = tempe;
		} else {
			e += ((oe - e) >> 1);
		}
	}

	while (e <= max_rid && flag[e - 1]) {
		++ e;
	}
	return e;
}

#ifdef false
void decompressPEX(string result) {
	fprintf(stderr, "decompressPE()...\n");
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string rootstrfn = folder + "rootstr.txt";
	mstcom::bsc::BSC_decompress((rootstrfn+".bsc").c_str(), rootstrfn.c_str());

	string orderfn = folder + string("order.bin");
	std::ifstream fporder;
	if (isorder) {
		mstcom::bsc::BSC_decompress((orderfn+".bsc").c_str(), orderfn.c_str());
		fporder.open(orderfn, std::ios::binary);
	}
	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string singletonfn = folder + "singleton.bin"; 
	mstcom::bsc::BSC_decompress((singletonfn+".bsc").c_str(), singletonfn.c_str());
	
	string readsnfn = folder + "readsN.txt";
	mstcom::bsc::BSC_decompress((readsnfn+".bsc").c_str(), readsnfn.c_str());

	string distfn = folder + "dist.bin";
	mstcom::bsc::BSC_decompress((distfn+".bsc").c_str(), distfn.c_str());

	BITPOOL curdirpool;
	curdirpool.n = 8;
	int dir, file;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");
	FILE *fproot = fopen(rootstrfn.c_str(), "r");
	// sprintf(name, "rootstr.bin");
	// std::ifstream fproot(name, std::ios::binary);

	FILE *fpdecstr1 = fopen((result+"_1").c_str(), "w");
	FILE *fpdecstr2 = fopen((result+"_2").c_str(), "w");
	char str[1<<10], enstr[1<<10], dstr[1<<10];;

	// int num = 112;

	if (isorder) {
		cout << "max_rid: " << max_rid << endl;
		uint32_t rid = 0;
		char **readsstr = new char*[max_rid];
		
		// root string
		while (fscanf(fproot, "%s", str) != EOF) {
			vector<READS_t> reads;
			READS_t r(str, 0);
			reads.push_back(r);
			int idx = 0;
			
			readsstr[rid ++] = strdup(str);
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					decodeStr(reads[idx].str, enstr, dstr);
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					
					readsstr[rid++] = strdup(dstr);

					READS_t r(dstr, idx);
					reads.push_back(r);
					idx = reads.size() - 1;

				} else { // back
					int back_step = atoi(enstr);
					while (back_step > 0) {
						idx = reads[idx].prid;
						--back_step;
					}
				}
			}
			reads.clear();
		}

		fclose(fpenstr);
		fclose(fproot);
		// fproot.close();
		fpdir.close();
		
		FILE *fpreadsn = fopen(readsnfn.c_str(), "r");
		while (fscanf(fpreadsn, "%s", str) != EOF) {
			readsstr[rid ++] = strdup(str);
		}
		fclose(fpreadsn);

		ifstream fpsingleton(singletonfn, std::ios::binary);
		DNAPOOL curdnapool;
		curdnapool.n = 4;
		int code;
		char *seq = (char*)alloca((L + 3) * sizeof(char));
		// bool flag = true;
		while (rid < max_rid) {
			for (int i = 0; i < L; ++i) {
				code = getDNAcode(fpsingleton, curdnapool);
				if (code < 0) {
					break;
				}
				seq[i] = invert_code_rule[code];
			}
			seq[L] = '\0';
			readsstr[rid++] = strdup(seq);
		}
		fpsingleton.close();

		uint32_t *ids = new uint32_t[max_rid + 3];

		uint32_t *tr, M;
		build(tr, M, max_rid);

		bool *flag = new bool[max_rid + 3];
		memset(flag, false, sizeof(bool)*max_rid); 

		uint32_t half = max_rid >> 1, pid;
		std::ifstream fpdist(distfn.c_str(), std::ios::binary);
		
		uint32_t tempid;
		int32_t dist;

		for (uint32_t num = 0; num < max_rid; ++num) {
			if (!flag[num]) {
				fporder.read((char*)&tempid, sizeof(uint32_t));
				// fscanf(fpdist, "%d", &dist);
				fpdist.read((char*)&dist, sizeof(int32_t));

				ids[num] = tempid;
				if (dist < 0) {
					dist = 0 - dist;
					ids[num] += half;
				}

				pid = binarySearch(tr, M, num + 1 + 1, max_rid, flag, dist);

				if (ids[num] < half) {
					ids[pid - 1] = half + ids[num];
				} else {
					ids[pid - 1] = ids[num] - half;
				}

				update(tr, M, pid, 1);
				flag[pid - 1] = true;
			}
		}

		fpdist.close();
		fporder.close();
		delete[] tr;
		delete[] flag;

		uint32_t *od = new uint32_t[max_rid + 3];
		for (uint32_t id = 0; id < max_rid; ++id) {
			od[ids[id]] = id;
		}
		delete[] ids;

		for (uint32_t id = 0; id < half; ++id) {
			fprintf(fpdecstr1, "%s\n", readsstr[od[id]]);
		}

		for (uint32_t id = half; id < max_rid; ++id) {
			fprintf(fpdecstr2, "%s\n", readsstr[od[id]]);
		}
	} else {
		uint32_t rid = 0;
		char **readsstr = new char*[max_rid + 1];

		while (fscanf(fproot, "%s", str) != EOF) {
			vector<READS_t> reads;
			READS_t r(str, 0);
			reads.push_back(r);
			int idx = 0;
			
			readsstr[rid ++] = strdup(str);
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					decodeStr(reads[idx].str, enstr, dstr);
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					readsstr[rid ++] = strdup(dstr);
					
					READS_t r(dstr, idx);
					reads.push_back(r);
					idx = reads.size() - 1;
				} else { // back
					int back_step = atoi(enstr);
					while (back_step > 0) {
						idx = reads[idx].prid;
						--back_step;
					}
				}
			}
			reads.clear();
		}

		fclose(fpenstr);
		fclose(fproot);
		fpdir.close();
		
		FILE *fpreadsn = fopen(readsnfn.c_str(), "r");
		while (fscanf(fpreadsn, "%s", str) != EOF) {
			readsstr[rid ++] = strdup(str);
		}
		fclose(fpreadsn);

		ifstream fpsingleton(singletonfn, std::ios::binary);
		DNAPOOL curdnapool;
		curdnapool.n = 4;
		int code;
		char *seq = (char*)alloca((L + 3) * sizeof(char));
		while (rid < max_rid) {
			for (int i = 0; i < L; ++i) {
				code = getDNAcode(fpsingleton, curdnapool);
				if (code < 0) {
					break;
				}
				seq[i] = invert_code_rule[code];
			}
			seq[L] = '\0';
			readsstr[rid ++] = strdup(seq);
		}
		fpsingleton.close();

		// cout << "xxxx\n";
		///
		uint32_t *ids = new uint32_t[max_rid + 3];

		uint32_t *tr, M;
		build(tr, M, max_rid);

		bool *flag = new bool[max_rid + 3];
		memset(flag, false, sizeof(bool)*max_rid); 

		uint32_t half = max_rid >> 1, pid;
		std::ifstream fpdist(distfn.c_str(), std::ios::binary);
		
		uint32_t tempid;
		int32_t dist;
		uint32_t id = 0;

		for (uint32_t num = 0; num < max_rid; ++num) {
			if (!flag[num]) {
				fpdist.read((char*)&dist, sizeof(int32_t));

				ids[num] = id ++;
				if (dist < 0) {
					dist = 0 - dist;
					ids[num] += half;
				}

				pid = binarySearch(tr, M, num + 1 + 1, max_rid, flag, dist);

				if (ids[num] < half) {
					ids[pid - 1] = half + ids[num];
				} else {
					ids[pid - 1] = ids[num] - half;
				}

				update(tr, M, pid, 1);
				flag[pid - 1] = true;
			}
		}

		fpdist.close();
		fporder.close();
		delete[] tr;
		delete[] flag;

		uint32_t *od = new uint32_t[max_rid + 3];
		for (uint32_t id = 0; id < max_rid; ++id) {
			od[ids[id]] = id;
		}
		delete[] ids;

		for (uint32_t id = 0; id < half; ++id) {
			fprintf(fpdecstr1, "%s\n", readsstr[od[id]]);
		}

		for (uint32_t id = half; id < max_rid; ++id) {
			fprintf(fpdecstr2, "%s\n", readsstr[od[id]]);
		}
	}

	fclose(fpdecstr1);
	fclose(fpdecstr2);
}
#endif

void decompressPEXX(string result) {
	fprintf(stderr, "decompressPE()...\n");
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string rootstrfn = folder + "rootstr.txt";
	mstcom::bsc::BSC_decompress((rootstrfn+".bsc").c_str(), rootstrfn.c_str());

	string orderfn = folder + string("order.bin");
	std::ifstream fporder;
	if (isorder) {
		mstcom::bsc::BSC_decompress((orderfn+".bsc").c_str(), orderfn.c_str());
		fporder.open(orderfn, std::ios::binary);
	}
	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string singletonfn = folder + "singleton.bin"; 
	mstcom::bsc::BSC_decompress((singletonfn+".bsc").c_str(), singletonfn.c_str());
	
	string readsnfn = folder + "readsN.txt";
	mstcom::bsc::BSC_decompress((readsnfn+".bsc").c_str(), readsnfn.c_str());

	string distfn = folder + "dist.bin";
	mstcom::bsc::BSC_decompress((distfn+".bsc").c_str(), distfn.c_str());

	BITPOOL curdirpool;
	curdirpool.n = 8;
	int dir, file;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");
	FILE *fproot = fopen(rootstrfn.c_str(), "r");
	// sprintf(name, "rootstr.bin");
	// std::ifstream fproot(name, std::ios::binary);

	FILE *fpdecstr1 = fopen((result+"_1").c_str(), "w");
	FILE *fpdecstr2 = fopen((result+"_2").c_str(), "w");
	char str[1<<10], enstr[1<<10], dstr[1<<10];;

	// int num = 112;
	uint32_v *idsv = new uint32_v[max_rid + 1];
	for (uint32_t i = 0; i < max_rid; ++i) {
		kv_init(idsv[i]);
	}
	// uint32_t *ids = new uint32_t[max_rid + 1];
	uint32_t idsid = 0;
	uint32_t *od = new uint32_t[max_rid + 3];

	if (isorder) {
		cout << "max_rid: " << max_rid << endl;
		uint32_t rid = 0;
		char **readsstr = new char*[max_rid];
		
		// root string
		// int input = 0;
		// cout << "aaaa\n";
		while (fscanf(fproot, "%s", str) != EOF) {
			// cout << str << endl;
			// fprintf(stderr, "%s\n", str);
			// cout << "BBBB\n";

			vector<READS_t> reads;
			READS_t r(str, rid, 0);
			reads.push_back(r);
			int idx = 0;
			
			kv_push(uint32_t, idsv[idsid], rid);
			od[rid] = idsid ++;

			readsstr[rid ++] = strdup(str);
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string

					decodeStr(reads[idx].str, enstr, dstr);
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					// rid 
					if (enstr[0] == '\n') {
					// if (strcmp(reads[idx].str.c_str(), dstr) == 0) {
						// cout << "xxxx\n";
						// exit (0);
						kv_push(uint32_t, idsv[od[reads[idx].rid]], rid);
						od[rid] = od[reads[idx].rid];
					} else {
						kv_push(uint32_t, idsv[idsid], rid);
						od[rid] = idsid ++;
					}


					READS_t r(dstr, rid, idx);
					reads.push_back(r);

					readsstr[rid ++] = strdup(dstr);
					idx = reads.size() - 1;

				} else { // back
					int back_step = atoi(enstr);
					while (back_step > 0) {
						idx = reads[idx].prid;
						--back_step;
					}
				}
			}
			reads.clear();
		}

		fclose(fpenstr);
		fclose(fproot);
		// fproot.close();
		fpdir.close();
		
		FILE *fpreadsn = fopen(readsnfn.c_str(), "r");
		while (fscanf(fpreadsn, "%s", str) != EOF) {
			kv_push(uint32_t, idsv[idsid], rid);
			od[rid] = idsid ++;

			readsstr[rid ++] = strdup(str);
		}
		fclose(fpreadsn);

		ifstream fpsingleton(singletonfn, std::ios::binary);
		DNAPOOL curdnapool;
		curdnapool.n = 4;
		int code;
		char *seq = (char*)alloca((L + 3) * sizeof(char));
		// bool flag = true;
		while (rid < max_rid) {
			for (int i = 0; i < L; ++i) {
				code = getDNAcode(fpsingleton, curdnapool);
				if (code < 0) {
					break;
				}
				seq[i] = invert_code_rule[code];
			}
			seq[L] = '\0';
			
			kv_push(uint32_t, idsv[idsid], rid);
			od[rid] = idsid ++;

			readsstr[rid++] = strdup(seq);
		}
		fpsingleton.close();

		uint32_t maxidsid = idsid;

		cout << "maxidsid: " << maxidsid << endl;
		FILE *fp2 = fopen("idx2.txt", "w");
		for (uint32_t i = 0; i < maxidsid; ++i) {
			for (uint32_t j = 0; j < idsv[i].n; ++j) {
				// fprintf(fp2, "%s\n", );
				fprintf(fp2, "%u ", idsv[i].a[j]);
			}
			fprintf(fp2, "\n");
		}
		fclose(fp2);
		// exit (0);


		uint32_t *ids = new uint32_t[max_rid + 3];

		uint32_t *tr, M;
		build(tr, M, max_rid);

		bool *flag = new bool[max_rid + 3];
		memset(flag, false, sizeof(bool)*max_rid); 
		bool *tvi = new bool[maxidsid + 3];
		memset(tvi, false, sizeof(bool)*maxidsid);

		uint32_t half = max_rid >> 1, pid;
		std::ifstream fpdist(distfn.c_str(), std::ios::binary);
		
		uint32_t tempid, trid, num, cnt = 0;
		int32_t dist;

		for (uint32_t i = 0; i < maxidsid; ++i) {
			for (uint32_t j = 0; j < idsv[i].n; ++j) {
				num = idsv[i].a[j];
				if (!flag[num]) {
					flag[num] = true;

					++ cnt;
					// cout << "cnt: " << cnt << endl;
					fporder.read((char*)&tempid, sizeof(uint32_t));
					// fscanf(fpdist, "%d", &dist);
					fpdist.read((char*)&dist, sizeof(int32_t));

					// cout << "i: " << i << endl;
					// cout << "tempid: " << tempid << endl;
					// cout << "dist: " << dist << endl;

					ids[num] = tempid;
					if (dist < 0) {
						dist = 0 - dist;
						ids[num] += half;
					}
					if (dist == 0) {
						trid = i + 1;
					} else {
					cout << "before binarySearch\n";
					
					trid = binarySearch(tr, M, i + 1 + 1, maxidsid, tvi, dist);
					}
					cout << "trid: " << trid << endl;

					pid = max_rid + 1;
					// get tid 
					if (trid == 348345) {
						cout << "idsv[trid - 1].n: " << idsv[trid - 1].n << endl;
					} 
					for (uint32_t k = 0; k < idsv[trid - 1].n; ++k) {
						// if (trid == 52798) cout << "idsv[trid - 1].a[k]: " << idsv[trid - 1].a[k] << endl;
						if (!flag[idsv[trid - 1].a[k]]) {
							pid = idsv[trid - 1].a[k];
							break;
						}
					}
					if (pid == max_rid + 1) {
						cout << "error\n";
						exit(0);
					}

					if (ids[num] < half) {
						ids[pid] = half + ids[num];
					} else {
						ids[pid] = ids[num] - half;
					}
					flag[pid] = true;

					if (!tvi[trid - 1]) {
						bool tmpflag = true;
						for (uint32_t k = 0; k < idsv[trid - 1].n; ++k) {
							if (!flag[idsv[trid - 1].a[k]]) {
								tmpflag = false;
								break;
							}
						}

						if (tmpflag) {
							update(tr, M, trid, 1);
							tvi[trid - 1] = true;
						}
					}
				}
			}
		}

		cout << "end \n";

#ifdef false
		for (uint32_t num = 0; num < max_rid; ++num) {
			if (!flag[num]) {
				fporder.read((char*)&tempid, sizeof(uint32_t));
				// fscanf(fpdist, "%d", &dist);
				fpdist.read((char*)&dist, sizeof(int32_t));

				ids[num] = tempid;
				if (dist < 0) {
					dist = 0 - dist;
					ids[num] += half;
				}

				pid = binarySearch(tr, M, num + 1 + 1, max_rid, flag, dist);

				if (ids[num] < half) {
					ids[pid - 1] = half + ids[num];
				} else {
					ids[pid - 1] = ids[num] - half;
				}

				update(tr, M, pid, 1);
				flag[pid - 1] = true;
			}
		}
#endif

		fpdist.close();
		fporder.close();
		delete[] tr;
		delete[] flag;
		delete[] tvi;

		// uint32_t *od = new uint32_t[max_rid + 3];
		for (uint32_t id = 0; id < max_rid; ++id) {
			od[ids[id]] = id;
		}
		delete[] ids;

		for (uint32_t id = 0; id < half; ++id) {
			fprintf(fpdecstr1, "%s\n", readsstr[od[id]]);
		}

		for (uint32_t id = half; id < max_rid; ++id) {
			fprintf(fpdecstr2, "%s\n", readsstr[od[id]]);
		}
	} 
#ifdef false
	else {
		uint32_t rid = 0;
		char **readsstr = new char*[max_rid + 1];

		while (fscanf(fproot, "%s", str) != EOF) {
			vector<READS_t> reads;
			READS_t r(str, 0);
			reads.push_back(r);
			int idx = 0;
			
			readsstr[rid ++] = strdup(str);
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					decodeStr(reads[idx].str, enstr, dstr);
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					readsstr[rid ++] = strdup(dstr);
					
					READS_t r(dstr, idx);
					reads.push_back(r);
					idx = reads.size() - 1;
				} else { // back
					int back_step = atoi(enstr);
					while (back_step > 0) {
						idx = reads[idx].prid;
						--back_step;
					}
				}
			}
			reads.clear();
		}

		fclose(fpenstr);
		fclose(fproot);
		fpdir.close();
		
		FILE *fpreadsn = fopen(readsnfn.c_str(), "r");
		while (fscanf(fpreadsn, "%s", str) != EOF) {
			readsstr[rid ++] = strdup(str);
		}
		fclose(fpreadsn);

		ifstream fpsingleton(singletonfn, std::ios::binary);
		DNAPOOL curdnapool;
		curdnapool.n = 4;
		int code;
		char *seq = (char*)alloca((L + 3) * sizeof(char));
		while (rid < max_rid) {
			for (int i = 0; i < L; ++i) {
				code = getDNAcode(fpsingleton, curdnapool);
				if (code < 0) {
					break;
				}
				seq[i] = invert_code_rule[code];
			}
			seq[L] = '\0';
			readsstr[rid ++] = strdup(seq);
		}
		fpsingleton.close();

		// cout << "xxxx\n";
		///
		uint32_t *ids = new uint32_t[max_rid + 3];

		uint32_t *tr, M;
		build(tr, M, max_rid);

		bool *flag = new bool[max_rid + 3];
		memset(flag, false, sizeof(bool)*max_rid); 

		uint32_t half = max_rid >> 1, pid;
		std::ifstream fpdist(distfn.c_str(), std::ios::binary);
		
		uint32_t tempid;
		int32_t dist;
		uint32_t id = 0;

		for (uint32_t num = 0; num < max_rid; ++num) {
			if (!flag[num]) {
				fpdist.read((char*)&dist, sizeof(int32_t));

				ids[num] = id ++;
				if (dist < 0) {
					dist = 0 - dist;
					ids[num] += half;
				}

				pid = binarySearch(tr, M, num + 1 + 1, max_rid, flag, dist);

				if (ids[num] < half) {
					ids[pid - 1] = half + ids[num];
				} else {
					ids[pid - 1] = ids[num] - half;
				}

				update(tr, M, pid, 1);
				flag[pid - 1] = true;
			}
		}

		fpdist.close();
		fporder.close();
		delete[] tr;
		delete[] flag;

		uint32_t *od = new uint32_t[max_rid + 3];
		for (uint32_t id = 0; id < max_rid; ++id) {
			od[ids[id]] = id;
		}
		delete[] ids;

		for (uint32_t id = 0; id < half; ++id) {
			fprintf(fpdecstr1, "%s\n", readsstr[od[id]]);
		}

		for (uint32_t id = half; id < max_rid; ++id) {
			fprintf(fpdecstr2, "%s\n", readsstr[od[id]]);
		}
	}
#endif
	fclose(fpdecstr1);
	fclose(fpdecstr2);
}

int main(int argc, char const *argv[]) {
	stopwatch.start();

	folder = generateString("mstcomdec", 10); //creat a temp folder for current input

	string cmd = "mkdir -p " + folder;
	fprintf(stderr, "%s\n", folder.c_str());
	system(cmd.c_str());
	folder += "/";
	cout << folder << endl;

	init();

	cmd = "tar -xf " + string(argv[1]) + " -C " + folder;
	// tar -xf $filename -C $decomp
	system(cmd.c_str());

	string parfn = folder + "par.txt";
	cout << parfn << endl;
	FILE *fp = fopen(parfn.c_str(), "r");
	fscanf(fp, "%d%d%d", &L, &ispe, &isorder);
	fprintf(stderr, "L: %d;\nispe: %d\nisorder: %d\n", L, ispe, isorder);
	if (isorder || ispe) {
		fscanf(fp, "%u", &max_rid);
		fprintf(stderr, "max_rid: %lu\n", max_rid);
	}
	fclose(fp);

	if (ispe) {
		// decompressPE(argv[2]);
		// decompressPEX(argv[2]);
		decompressPEXX(argv[2]);
	} else {
		decompressSingle(argv[2]);
	}

	cmd = "rm -rf " + folder;
	system(cmd.c_str());

	cout << "Time of decompressioin = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	return 0;
}

CStopWatch::CStopWatch(){
	running = 0;
}


CStopWatch::~CStopWatch(){
}


void CStopWatch::start(){
	running = 1;
	elapsed.clear();
	t1 = std::chrono::steady_clock::now();
}


double CStopWatch::stop(){
	t2  = std::chrono::steady_clock::now();
	running = 0;
	double el = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0;
	elapsed.push_back(el);
	return el;
}


void CStopWatch::resume(){
	running = 1;
	t1 = std::chrono::steady_clock::now();
}


double CStopWatch::totalTime(){
	double sum = 0;
	for (std::vector<double>::iterator it = elapsed.begin(); it != elapsed.end(); ++it) {
		sum += *it;
	}
	return sum;
}
