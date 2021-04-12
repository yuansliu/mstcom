# include "mstcom.h"
using namespace std;

class dec_class{
int L, L1, L2;
int ispe, isorder;
uint32_t max_ctg_length;
size_t max_rid;
const int nthreads = 24;
string folder, infile, outfile;
uint32_t MAXNO;

struct READS_t {
	string str;
	size_t prid, prechildrid, n; //reads id of parent
	READS_t(): prechildrid(0), n(0) {}
	READS_t(char *_str): str(_str), n(0) {}
	READS_t(char *_str, const size_t &_prid): str(_str), prid(_prid), prechildrid(0), n(0) {}
	READS_t(char *_str, const size_t &_prid, size_t _n): str(_str), prid(_prid), prechildrid(0), n(_n) {}
	void set(char *_str) { 
		str = _str;
	}
};

char complement[256];

const char invert_code_rule[5] = {'A', 'C', 'G', 'T', 'N'}; //decoding rule //A - 0; C - 1; G - 2; T - 3;

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

inline void reverseReads(char* start) {
	char* left = start; // sequence starts
	char* right = start + L - 1;
	while (right > left) {
		char tmp = *left;
		*left = *right;
		*right = tmp;
		++left;
		--right;
	}
}

inline void reverseComplement(char* start, int L) {
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

inline void reverseReads(char* start, int L) {
	char* left = start; // sequence starts
	char* right = start + L - 1;
	while (right > left) {
		char tmp = *left;
		*left = *right;
		*right = tmp;
		++left;
		--right;
	}
}

inline void decodeStr(string pstr, char *enstr, char *dstr) {
	bool debug = false;
	// debug = true;
	// if (strcmp(enstr, "68A3A3CACAGAAATCG3GAG3CGAGG") == 0) debug = true;
	// if (debug) {
	// 	fprintf(stderr, "special pstr: %s\n", pstr.c_str());
	// 	exit(0);
	// }
	// if (enstr[0] == ' ' && enstr[1] == '\n') { // equal
	if (enstr[0] == '\0') { // equal
		strcpy(dstr, pstr.c_str());
	} else {
		// if (strcmp(enstr, "6GAGAT32T3C7T30A TCATATACAGG\n") == 0) debug = true;
		// if (strcmp(enstr, "11G5T10T32CGGA12T5AATTTTT2T TTTTTTTT") == 0) debug = true;
		// if (strcmp(enstr, "AAAA 16G5G6G20A5A2T4C19A9C\n") == 0) debug = true;
		// if (strcmp(enstr, " CTACCATATTCGCCACAGCCCCA\n") == 0) debug = true;
		// if (strcmp(enstr, "TG TGAAGGCGG\n") == 0) debug = true;
		// if (strcmp(enstr, "TGACACGAAACAAGGAAGT \n") == 0) debug = true;
		// if (strcmp(enstr, "C AT\n") == 0) debug = true;
		// if (strcmp(enstr, "TG TGAAGGCGG\n") == 0) debug = true;
		// check space ' '
		int spidx = -1;
		// int len = strlen(enstr) - 1;
		int len = strlen(enstr);
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
			// AATTATGC
			//   TTATGCAA
		} else 
		if (spidx > -1) {
			// AATTATGC
			//   TTACGCAA
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

		if (spidx > 0 && spidx == len - 1) {
			shiftRight = false;
			// cout << pstr << endl;
			// cout << enstr << endl;
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
			dstr[dstridx] = '\0';
			cout << enstr << endl;
			cout << "dstr: " << dstr << endl;
			fprintf(stderr, "pstr: %s\n", pstr.c_str());
			fprintf(stderr, "spidx: %d\n", spidx);
			fprintf(stderr, "shiftRight: %d\n", shiftRight);
			fprintf(stderr, "len: %d\n", len);
			fprintf(stderr, "maxenstridx: %d\n", maxenstridx);
			fprintf(stderr, "dstridx: %d\n", dstridx);
			fprintf(stderr, "pstridx: %d\n", pstridx);
					cout << "L: " << L << endl;
		}

		int eq_char_num = 0;
		for (; enstridx < maxenstridx; ++enstridx) {
			if (debug) cout << "enstridx: " << enstridx << endl;
			if (enstr[enstridx] >= 'A' && enstr[enstridx] <= 'Z') {
				if (eq_char_num > 0) {
					if (debug) cout << "eq_char_num: " << eq_char_num << endl;
					for (int j = 0; j < eq_char_num; ++j) {
						dstr[dstridx++] = pstr[pstridx++];
					}
					eq_char_num = 0;
				}
				dstr[dstridx++] = enstr[enstridx];
				if (debug) {
					// fprintf(stderr, "dstr[%d]: %c\n", dstridx-1, dstr[dstridx-1]);
					// fprintf(stderr, "enstr[%d]: %c\n", enstridx, enstr[enstridx]);
				}
				pstridx++;
			} else {
				eq_char_num = eq_char_num * 10 + enstr[enstridx] - '0';
			}
			if (debug) {
				dstr[dstridx] = '\0';
				cout << "dstr: " << dstr << endl;
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

inline void decodeStrFromCtg(char *pstr, char *enstr, char *dstr, bool isleftshift) {
	// isleftshift == 1; leftshift
	bool debug = false;
	// if (strcmp(enstr, "83G GAA") == 0) debug = true;
	// check space ' '
	int spidx = -1;
	int len = strlen(enstr);
	for (int i = 0; i < len; ++i) {
		if (enstr[i] == ' ') {
			spidx = i;
			break;
		}
	}
	int dstridx = 0;
	// fprintf(stderr, "--spidx: %d\n", spidx);

	// fprintf(stderr, "pstr: %s\n", pstr.c_str());
	bool shiftRight = false;
	if (spidx > -1) {
		shiftRight = !isleftshift;
	}
	if (debug) {
		cout << "shiftRight: " << shiftRight << endl;
		cout << "spidx: " << spidx << endl;
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
		fprintf(stderr, "pstr: %s\n", pstr);
		fprintf(stderr, "spidx: %d\n", spidx);
		fprintf(stderr, "shiftRight: %d\n", shiftRight);
		fprintf(stderr, "len: %d\n", len);
		fprintf(stderr, "maxenstridx: %d\n", maxenstridx);
		fprintf(stderr, "dstridx: %d\n", dstridx);
		fprintf(stderr, "pstridx: %d\n", pstridx);
	}

	if (debug) cout << "enstridx: " << enstridx << endl;

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

inline bool getDir(std::ifstream& fpdir, BITPOOL& curdirpool, int &res) {
	if (curdirpool.n >= 8) {
		uint8_t dirbin;

		if (!fpdir.read((char*)&dirbin, sizeof(uint8_t))) {
			return false;
		}

		for (int i = 0; i < 8; ++i) {
		// for (int i = 7; i >= 0; --i) {
			curdirpool.a[i] = dirbin&1;
			dirbin >>= 1;
		}
		curdirpool.n = 0;
	}
	res = curdirpool.a[curdirpool.n];
	++curdirpool.n;
	return true;
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

struct BASESPOOL {
	int n, a[5];
};

inline int getBasesCode(std::ifstream& fpif, BASESPOOL& curbasespool) {
	// fprintf(stderr, "yyyyyy\n");
	if (curbasespool.n >= 5) {
		uint16_t bin;
		// fprintf(stderr, "zzzz\n");
		fpif.read((char*)&bin, sizeof(uint16_t));
		// fprintf(stderr, "%d\n", bin);
		if (fpif.eof()) return -1;

		for (int i = 0; i < 5; ++i) {
			curbasespool.a[i] = bin&7;
			bin >>= 3;
		}
		curbasespool.n = 0;
	}
	int res = curbasespool.a[curbasespool.n];
	++curbasespool.n;
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

bool checkStrOrDig_v0(char *str) {
	int len = strlen(str);
	// if (str[0] == ' ' && str[1] == '\n') { // ' \n'
	if (str[0] == '\n') { // ' \n'
		return true;
	}
	if (len >= L) return true;
	for (int i = 0; i < len - 1; ++i) {
		if ((str[i] >= 'A' && str[i] <= 'Z') || str[i] == '$') {
			return true; // string
		}
	}
	return false; // dig
}

bool checkStrOrDig(char *str) {
	int len = strlen(str);
	// if (str[0] == ' ' && str[1] == '\n') { // ' \n'
	if (str[0] == '\n') { // ' \n'
		return true;
	}
	if (len >= L) return true;
	for (int i = 0; i < len; ++i) {
		if ((str[i] >= 'A' && str[i] <= 'Z') || str[i] == '$') {
			return true; // string
		}
	}
	return false; // dig
}

void decompressFile(string file) {
	std::ifstream f(file + ".lzma", std::ios::binary);
	if (f.fail()) { // bsc
		cout << file + ".bsc" << endl;
		mstcom::bsc::BSC_decompress((file+".bsc").c_str(), file.c_str());
	} else {
		f.close();
		cout << file + ".lzma" << endl;
		mstcom::lzma::lzma_decompress((file+".lzma").c_str(), file.c_str());
	}
}

int debugcnt = 0, mindebugcnt = 0;
	
inline void decodeStr(string pstr, char *enstr, int curshiftdirection, char *dstr) {
	// bool debug = false;

	// if (debugcnt >= mindebugcnt && debugcnt < mindebugcnt + 10) {
		// debug = true;
	// }
	// ++ debugcnt;

	// if (debug) {
	// 	cout << "pstr: " << pstr << endl;
	// 	cout << "shift direction: " << curshiftdirection << endl;
	// 	cout << enstr << endl;
	// }

	string shiftpartsubstr = "";
	int len = strlen(enstr);
	int i = 0, j, k;
	// get shift part
	for (i = 0; i < len; ++i) {
		if (enstr[i] >= '0' && enstr[i] <= '9') break;
		shiftpartsubstr += enstr[i];
	}
	int shiftlen = shiftpartsubstr.length();

	// if (debug) {
		// cout << "shiftpartsubstr: " << shiftpartsubstr << endl;
	// }

	if (shiftlen != 0) {
		if (curshiftdirection == 0) {
			for (j = 0; j < shiftlen; ++j) {
				dstr[j] = shiftpartsubstr[j];
			}
			for (k = 0; j < L; ++j, ++k) {
				dstr[j] = pstr[k];
			}
		} else {
			for (j = 0, k = shiftlen; k < L; ++j, ++k) {
				dstr[j] = pstr[k];
			}
			for (k = 0; k < shiftlen; ++k, ++j) {
				dstr[j] = shiftpartsubstr[k];
			}
		}
	} else { // shiftlen == 0
		for (j = 0, k = 0; k < L; ++j, ++k) {
			dstr[j] = pstr[k];
		}
	}
	dstr[L] = '\0';

	// get mismatch positions and bases
	int v = 0, prepos = 0, pos;
	if (curshiftdirection == 0) {
		prepos += shiftlen;
	}
	for (; i < len; ++i) {
		if (enstr[i] >= '0' && enstr[i] <= '9') {
			v = v * 10 + enstr[i] - '0';
		} else { // a base
			pos = prepos + v;
			dstr[pos] = enstr[i];
			// if (debug) {
			// 	// cout << "pos: " << pos << "; base: " << enstr[i] << endl;
			// }
			prepos = pos + 1;
			v = 0;
		}
	}
	// if (debug) {
	// 	fprintf(stderr, "%s\n", dstr);
	// 	// exit(1);
	// }
}

inline void decodeStr(string pstr, char *enstr, int curshiftdirection, int filelistval, char *dstr) {

	// if (debugcnt >= mindebugcnt && debugcnt < mindebugcnt + 10) {
	// bool debug = false;
	// if (strcmp(enstr, "ACTCCGTTAAGAGGACGGAGTCAGAGGGTT3T1N14G5A6C11N0N8T0A5G6G") == 0
	// if (strcmp(enstr, "CAGATTTAATTTTGAGTTCTTGAGGAATACGAAA1A13T5G2G14A14G12T") == 0
	// if (strcmp(enstr, "TGT14G0T9T6T5G0T5T1C12C32T") == 0
	// if (strcmp(enstr, "TTT7A17T10A12G0T20A9T14C") == 0
	// if (strcmp(enstr, "AGGT27A0C21T1C42G") == 0
	// if (strcmp(enstr, "AA8A1A3A3A3A3A3A3A3A2C0A3A3A7A0G2A0G2A0G3G0A0T1G0A2G0A10A0G") == 0
	// 	// && strcmp(pstr.c_str(), "ACAAAAAAACTCTACCTCAAAACGAGAACAATGGGTCCGGCCTCACGTTACCGCACTAGAACCGAGTAACGTTGAAGACGGAGGACCCAAGTTCGTTAAG") == 0
	// 	) {
	// 	debug = true;
	// }
	// if (debug) {
	// 	// cout << "pstr: " << pstr << endl;
	// 	cout << "shift direction: " << curshiftdirection << endl;
	// 	cout << enstr << endl;
	// }
	// ++ debugcnt;
	int pl = pstr.length();
	int cl = filelistval == 0? L1: L2;

	if (pl == cl) {
		// decodeStr(pstr, enstr, curshiftdirection, dstr);
		string shiftpartsubstr = "";
		int len = strlen(enstr);
		int i = 0, j, k;
		// get shift part
		for (i = 0; i < len; ++i) {
			if (enstr[i] >= '0' && enstr[i] <= '9') break;
			shiftpartsubstr += enstr[i];
		}
		int shiftlen = shiftpartsubstr.length();

		// if (debug) {
			// cout << "shiftpartsubstr: " << shiftpartsubstr << endl;
		// }

		if (shiftlen != 0) {
			if (curshiftdirection == 0) {
				for (j = 0; j < shiftlen; ++j) {
					dstr[j] = shiftpartsubstr[j];
				}
				for (k = 0; j < pl; ++j, ++k) {
					dstr[j] = pstr[k];
				}
			} else {
				for (j = 0, k = shiftlen; k < pl; ++j, ++k) {
					dstr[j] = pstr[k];
				}
				for (k = 0; k < shiftlen; ++k, ++j) {
					dstr[j] = shiftpartsubstr[k];
				}
			}
		} else { // shiftlen == 0
			for (j = 0, k = 0; k < pl; ++j, ++k) {
				dstr[j] = pstr[k];
			}
		}
		dstr[pl] = '\0';

		// get mismatch positions and bases
		int v = 0, prepos = 0, pos;
		if (curshiftdirection == 0) {
			prepos += shiftlen;
		}
		for (; i < len; ++i) {
			if (enstr[i] >= '0' && enstr[i] <= '9') {
				v = v * 10 + enstr[i] - '0';
			} else { // a base
				pos = prepos + v;
				dstr[pos] = enstr[i];
				if (debug) {
					// cout << "pos: " << pos << "; base: " << enstr[i] << endl;
				}
				prepos = pos + 1;
				v = 0;
			}
		}
	} else {
		// pl != cl
		// cout << "111" << endl;
		string shiftpartsubstr = "";
		int len = strlen(enstr);
		int i = 0, j, k;
		// get shift part
		for (i = 0; i < len; ++i) {
			if (enstr[i] >= '0' && enstr[i] <= '9') break;
			shiftpartsubstr += enstr[i];
		}
		int shiftlen = shiftpartsubstr.length();

		// if (debug) {
		// 	cout << "shiftpartsubstr: " << shiftpartsubstr << endl;
		// }

		// always hold shiftlen > 0
		if (shiftlen != 0) {
			if (curshiftdirection == 0) {
				// cast: curshiftdirection == 0
				//       AATTGCATGC parent
				//    TCGAATTGCA
				for (j = 0; j < shiftlen; ++j) {
					dstr[j] = shiftpartsubstr[j];
				}
				for (k = 0; j < cl; ++j, ++k) {
					dstr[j] = pstr[k];
				}
			} else {
				// case curshiftdirection == 1
				//    AATTGCATGC parent
				//      TTGCATGCGAAAA  

				// if (debug) cout << "pl - (cl - shiftlen): " << pl - (cl - shiftlen) << endl;

				for (j = 0, k = pl - (cl - shiftlen); k < pl; ++j, ++k) {
					dstr[j] = pstr[k];
				}
				for (k = 0; k < shiftlen; ++k, ++j) {
					dstr[j] = shiftpartsubstr[k];
				}
				// shiftlen = pl - (cl - shiftlen);
			}
		} 
		dstr[cl] = '\0';
		// cout << "222" << endl;
		// get mismatch positions and bases
		int v = 0, prepos = 0, pos;
		if (curshiftdirection == 0) {
			prepos += shiftlen;
		}
		for (; i < len; ++i) {
			if (enstr[i] >= '0' && enstr[i] <= '9') {
				v = v * 10 + enstr[i] - '0';
			} else { // a base
				pos = prepos + v;
				dstr[pos] = enstr[i];
				if (debug) {
					// cout << "pos: " << pos << "; base: " << enstr[i] << endl;
				}
				prepos = pos + 1;
				v = 0;
			}
		}
		// if (debug) {
		// 	fprintf(stderr, "%s\n", dstr);
		// 	// exit(1);
		// }
		// cout << "333" << endl;
	}

	// if (debug) {
	// 	cout << "pstr: " << pstr << endl;
	// 	cout << "dstr: " << dstr << endl;
	// 	cout << "----" << endl;
	// }
}

inline void decodeStr_v1(string pstr, char *enstr, int curshiftdirection, char *dstr) {
	bool debug = false;
	if (debugcnt >= mindebugcnt && debugcnt < mindebugcnt + 10) {
		// debug = true;
	}
	// ++ debugcnt;
	// if (strcmp(enstr, "C33A") == 0) debug = true;
	if (debug) {
		// cout << "shift direction: " << curshiftdirection << endl;
		// cout << enstr << endl;
	}

	if (debug) {
		cout << "---------" << endl;
		cout << "pstr: " << pstr << endl;
		cout << "enstr: " << enstr << endl;
		cout << "---------" << endl;
	}
	int len = strlen(enstr);
	int i = 0, j = 0, k, t, shiftlen;
	// get mismatch positions and bases
	int v = 0, prepos = 0, pos;
	
	if (curshiftdirection) { // pure shift and shift offset < 0
		for (i = 0; i < len; ++i) {
			dstr[i] = enstr[i];
		}
		shiftlen = len;
		for (j = 0, k = shiftlen; k < L; ++j, ++k) {
			dstr[k] = pstr[j];
		}
	} else {
		if (isupper(enstr[0])) {
			t = 0;
			while (isupper(enstr[t + 1])) {
				++ t;
			}
			if (t == len - 1) {
				// pure shift and shift offset > 0
				shiftlen = len;
				for (j = 0, k = shiftlen; k < L; ++j, ++k) {
					dstr[j] = pstr[k];
				}
				for (k = 0; k < shiftlen; ++k, ++j) {
					dstr[j] = enstr[k];
				}
			} else 
			if (t >= 0) { // at least one base
				// shift with mismatch; shift offset < 0
				if (debug) {
					cout << "t: " << t << endl;
				}
				for (i = 0; i <= t; ++i) {
					dstr[i] = enstr[i];
				}
				shiftlen = t + 1;
				for (j = 0, k = shiftlen; k <= L; ++j, ++k) {
					// if (debug) {
						// cout << "j: " << j << "; k: " << k << endl;
					// }
					dstr[k] = pstr[j];
				}
				prepos += shiftlen;

				if (debug) {
					dstr[L] = '\0';
					cout << "dstr: " << dstr << endl;
					cout << "len: " << len << endl;
				}
				for (i = shiftlen; i < len; ++i) {
					if (enstr[i] >= '0' && enstr[i] <= '9') {
						v = v * 10 + enstr[i] - '0';
					} else { // a base
						pos = prepos + v;
						dstr[pos] = enstr[i];
						if (debug) {
							// cout << "pos: " << pos << "; base: " << enstr[i] << endl;
						}
						prepos = pos + 1;
						v = 0;
					}
				}
				if (debug) {
					cout << "dstr: " << dstr << endl;
				}
			} 
		} else { // t == 0; shift offset >= 0
			t = len - 1;
			while (isupper(enstr[t - 1])) {
				-- t;
			}
			++ t;
			if (t < len) {
				shiftlen = len - t;
			} else {
				shiftlen = 0;
			}
			// if (debug) {
			// 	cout << "t: " << t << endl;
			// 	cout << "shiftlen: " << shiftlen << endl;
			// }
			for (j = 0, k = shiftlen; k < L; ++j, ++k) {
				dstr[j] = pstr[k];
			}
			// if (debug) {
			// 	cout << "j: " << j << endl;
			// 	dstr[j] = '\0';
			// 	cout << "dstr: " << dstr << endl;
			// }
			for (k = t; k < len; ++k, ++j) {
				// if (debug) {
				// 	cout << "j: " << j << endl;
				// }
				dstr[j] = enstr[k];
			}
			// if (debug) {
			// 	dstr[j] = '\0';
			// 	cout << "dstr: " << dstr << endl;
			// }
			len = t;
			// if (debug) {
				// cout << "len: " << len << endl;
			// }
			for (; i < len; ++i) {
				if (enstr[i] >= '0' && enstr[i] <= '9') {
					v = v * 10 + enstr[i] - '0';
				} else { // a base
					pos = prepos + v;
					dstr[pos] = enstr[i];
					if (debug) {
						// cout << "pos: " << pos << "; base: " << enstr[i] << endl;
					}
					prepos = pos + 1;
					v = 0;
				}
			}
		}
	}
	dstr[L] = '\0';
	if (debug) {
		cout << "dstr len: " << strlen(dstr) << endl; 
		// fprintf(stderr, "%s\n", dstr);
		// exit(1);
	}
}

inline void decodeStr(string pstr, char *mismatchstr, char *shiftpartsubstr, int curshiftdirection, char *dstr) {
	bool debug = false;

	if (debugcnt >= mindebugcnt && debugcnt < mindebugcnt + 10) {
		// debug = true;
	}
	++ debugcnt;

	if (debug) {
		cout << "pstr: " << pstr << endl;
		cout << "shift direction: " << curshiftdirection << endl;
		cout << "mismatchstr: " << mismatchstr << endl;
		cout << "shiftstr: " << shiftpartsubstr << endl;
	}

	int len = strlen(mismatchstr);
	int i = 0, j, k;
	int shiftlen = strlen(shiftpartsubstr);

	if (debug) {
		// cout << "shiftpartsubstr: " << shiftpartsubstr << endl;
	}

	if (shiftlen != 0) {
		if (curshiftdirection == 0) {
			for (j = 0; j < shiftlen; ++j) {
				dstr[j] = shiftpartsubstr[j];
			}
			for (k = 0; j < L; ++j, ++k) {
				dstr[j] = pstr[k];
			}
		} else {
			for (j = 0, k = shiftlen; k < L; ++j, ++k) {
				dstr[j] = pstr[k];
			}
			for (k = 0; k < shiftlen; ++k, ++j) {
				dstr[j] = shiftpartsubstr[k];
			}
		}
	} else { // shiftlen == 0
		for (j = 0, k = 0; k < L; ++j, ++k) {
			dstr[j] = pstr[k];
		}
	}
	dstr[L] = '\0';

	// get mismatch positions and bases
	int v = 0, prepos = 0, pos;
	if (curshiftdirection == 0) {
		prepos += shiftlen;
	}
	for (i = 0; i < len; ++i) {
		if (isdigit(mismatchstr[i])) {
			v = v * 10 + mismatchstr[i] - '0';
		} else { // a base
			pos = prepos + v;
			dstr[pos] = mismatchstr[i];
			if (debug) {
				// cout << "pos: " << pos << "; base: " << enstr[i] << endl;
			}
			prepos = pos + 1;
			v = 0;
		}
	}
	if (debug) {
		fprintf(stderr, "%s\n", dstr);
		// exit(1);
	}
}

// void decompressSingleDFS_best_version_v0(string result) {
void decompressSingleDFS(string result) {
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string isleftbinfn = folder + "isleft.bin";
	decompressFile(isleftbinfn);
	std::ifstream fpisleft(isleftbinfn, std::ios::binary);

	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());

	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string dupfn = folder + "dup.bin";
	decompressFile(dupfn);
	std::ifstream fpdup(dupfn.c_str(), std::ios::binary);

	string isdupfn = folder + "isdup.bin";
	decompressFile(isdupfn);
	std::ifstream fpisdup(isdupfn.c_str(), std::ios::binary);

	BITPOOL curdirpool, curisduppool, curisleftpool;
	curdirpool.n = 8;
	curisduppool.n = 8;
	curisleftpool.n = 8;
	int dir;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");

	FILE *fpdecstr = fopen(result.c_str(), "w");
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));
	max_ctg_length = 1024;
	char *enstr = new char[max_ctg_length];
	char *dstr = (char*)alloca((L + 3) * sizeof(char));

	int isori, dn, idx;
	vector<READS_t> reads;
	// FILE *fpw = fopen("ddec.seq", "w");
	uint32_t dupno, rcdupno;

	int curshiftdirection;

	uint32_t printnum = 0, dupnum = 0;

	// while (fgets(enstr, max_ctg_length, fpenstr) != NULL) {
	while (fscanf(fpenstr, "%s", enstr) != EOF) {
		if (enstr[0] == '-') break;

		if (checkStrOrDig(enstr)) { // is encode string
			// cout << enstr << endl;
			uint32_t enstrlen = strlen(enstr);
			
			if (enstrlen < L) { // 
				// 0 shift > 0; 1 shift < 1;
				if (enstr[0] >= '0' && enstr[0] <= '9') { // no shift
					curshiftdirection = 0;
				} else {
					curshiftdirection = getDir(fpisleft, curisleftpool);
				}
				// 0 for left shift; 1 for right shift;

				decodeStr(reads[idx].str, enstr, curshiftdirection, dstr); //

				dir = getDir(fpdir, curdirpool);
				if (dir) {
					reverseComplement(dstr);
				}

				READS_t r(dstr, idx);
				reads.push_back(r);
				idx = reads.size() - 1;

				// fprintf(fpw, "%s\n", dstr);
				string tempstr = dstr;

				reverseReads(dstr);

				fprintf(fpdecstr, "%s\n", dstr);
				// ++ printnum;
				// if (strcmp(dstr, "83GGAAATGGCTGGGTCAAATGGTATTTCTAGTTCTAGATCCCTGAGGAATCGCCACACTGACTTCCACAATGGTTGAACTAGTTTACATTCCCACCAAGA") == 0) 
				// 	cout << "dstr in !isctgtree\n";

				int isdup = getDir(fpisdup, curisduppool);

				if (isdup) {
					reverseReads(dstr);
					strcpy(rcstr, dstr);
					reverseComplement(rcstr);

					reverseReads(dstr);
					reverseReads(rcstr);

					fpdup.read((char*)&dupno, sizeof(uint32_t));
					fpdup.read((char*)&rcdupno, sizeof(uint32_t));

					// for (uint32_t w = 0; w < dupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }
					// for (uint32_t w = 0; w < rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }

					for (int w = 0; w < dupno; ++w) {
						fprintf(fpdecstr, "%s\n", dstr);
					}
					for (int w = 0; w < rcdupno; ++w) {
						fprintf(fpdecstr, "%s\n", rcstr);
					}
				}

			} else { // root node // enstrlen == L
				reads.clear();

				READS_t r(enstr, 0);
				reads.push_back(r);
				idx = 0;

				// preshiftdirection = 0; // 0 for left shift; 1 for right shift;

				// fprintf(fpw, "%s\n", enstr);

				string tempstr = enstr;

				reverseReads(enstr);

				fprintf(fpdecstr, "%s\n", enstr);
				// ++ printnum;

				int isdup = getDir(fpisdup, curisduppool);

				if (isdup) {
					reverseReads(enstr);
					strcpy(rcstr, enstr);
					reverseComplement(rcstr);

					reverseReads(enstr);
					reverseReads(rcstr);

					fpdup.read((char*)&dupno, sizeof(uint32_t));
					fpdup.read((char*)&rcdupno, sizeof(uint32_t));

					// for (uint32_t w = 0; w < dupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }
					// for (uint32_t w = 0; w < rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }

					for (int w = 0; w < dupno; ++w) {
						fprintf(fpdecstr, "%s\n", enstr);
					}
					for (int w = 0; w < rcdupno; ++w) {
						fprintf(fpdecstr, "%s\n", rcstr);
					}
				}
			}
		} else {
			int back_step = atoi(enstr);
			while (back_step > 0) {
				idx = reads[idx].prid;
				--back_step;
			}
		}
	}

	delete[] enstr;

	// cout << "debugcnt: " << debugcnt << endl;
	// cout << "printnum: " << printnum << endl;
	// cout << "dupnum: " << dupnum << endl;
	// fclose(fpw);
	fpdir.close();
	fpisdup.close();
	fpdup.close();
	fpisleft.close();
	// FILE *fpencsg = fopen("decsg.txt", "w");
	while (fscanf(fpenstr, "%s", dstr) != EOF) {

		// fprintf(fpencsg, "%s\n", dstr);

		reverseReads(dstr);

		fprintf(fpdecstr, "%s\n", dstr);
	}
	fclose(fpenstr);
	// fclose(fpencsg);

	fclose(fpdecstr);
}

// reads singleton from a seperate file
void decompressSingleDFS_best_version_v0(string result) {
// void decompressSingleDFS(string result) {
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string isleftbinfn = folder + "isleft.bin";
	decompressFile(isleftbinfn);
	std::ifstream fpisleft(isleftbinfn, std::ios::binary);

	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	
	string singletonfn = folder + "singleton.bin";
	mstcom::bsc::BSC_decompress((singletonfn+".bsc").c_str(), singletonfn.c_str());

	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string dupfn = folder + "dup.bin";
	decompressFile(dupfn);
	std::ifstream fpdup(dupfn.c_str(), std::ios::binary);

	string isdupfn = folder + "isdup.bin";
	decompressFile(isdupfn);
	std::ifstream fpisdup(isdupfn.c_str(), std::ios::binary);

	BITPOOL curdirpool, curisduppool, curisleftpool;
	curdirpool.n = 8;
	curisduppool.n = 8;
	curisleftpool.n = 8;
	int dir;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");

	FILE *fpdecstr = fopen(result.c_str(), "w");
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));
	max_ctg_length = 1024;
	char *enstr = new char[max_ctg_length];
	char *dstr = (char*)alloca((L + 3) * sizeof(char));

	int isori, dn, idx;
	vector<READS_t> reads;
	// FILE *fpw = fopen("ddec.seq", "w");
	uint32_t dupno, rcdupno;

	int curshiftdirection;

	uint32_t printnum = 0, dupnum = 0;

	// while (fgets(enstr, max_ctg_length, fpenstr) != NULL) {
	while (fscanf(fpenstr, "%s", enstr) != EOF) {
		if (enstr[0] == '-') break;

		if (checkStrOrDig(enstr)) { // is encode string
			// cout << enstr << endl;
			uint32_t enstrlen = strlen(enstr);
			
			if (enstrlen < L) { // 
				// 0 shift > 0; 1 shift < 1;
				if (enstr[0] >= '0' && enstr[0] <= '9') { // no shift
					curshiftdirection = 0;
				} else {
					curshiftdirection = getDir(fpisleft, curisleftpool);
				}
				// 0 for left shift; 1 for right shift;

				decodeStr(reads[idx].str, enstr, curshiftdirection, dstr); //

				dir = getDir(fpdir, curdirpool);
				if (dir) {
					reverseComplement(dstr);
				}

				READS_t r(dstr, idx);
				reads.push_back(r);
				idx = reads.size() - 1;

				// fprintf(fpw, "%s\n", dstr);
				string tempstr = dstr;

				reverseReads(dstr);

				fprintf(fpdecstr, "%s\n", dstr);
				// ++ printnum;
				// if (strcmp(dstr, "83GGAAATGGCTGGGTCAAATGGTATTTCTAGTTCTAGATCCCTGAGGAATCGCCACACTGACTTCCACAATGGTTGAACTAGTTTACATTCCCACCAAGA") == 0) 
				// 	cout << "dstr in !isctgtree\n";

				int isdup = getDir(fpisdup, curisduppool);

				if (isdup) {
					reverseReads(dstr);
					strcpy(rcstr, dstr);
					reverseComplement(rcstr);

					reverseReads(dstr);
					reverseReads(rcstr);

					fpdup.read((char*)&dupno, sizeof(uint32_t));
					fpdup.read((char*)&rcdupno, sizeof(uint32_t));

					// for (uint32_t w = 0; w < dupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }
					// for (uint32_t w = 0; w < rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }

					for (int w = 0; w < dupno; ++w) {
						fprintf(fpdecstr, "%s\n", dstr);
					}
					for (int w = 0; w < rcdupno; ++w) {
						fprintf(fpdecstr, "%s\n", rcstr);
					}
				}

			} else { // root node // enstrlen == L
				reads.clear();

				READS_t r(enstr, 0);
				reads.push_back(r);
				idx = 0;

				// preshiftdirection = 0; // 0 for left shift; 1 for right shift;

				// fprintf(fpw, "%s\n", enstr);

				string tempstr = enstr;

				reverseReads(enstr);

				fprintf(fpdecstr, "%s\n", enstr);
				// ++ printnum;

				int isdup = getDir(fpisdup, curisduppool);

				if (isdup) {
					reverseReads(enstr);
					strcpy(rcstr, enstr);
					reverseComplement(rcstr);

					reverseReads(enstr);
					reverseReads(rcstr);

					fpdup.read((char*)&dupno, sizeof(uint32_t));
					fpdup.read((char*)&rcdupno, sizeof(uint32_t));

					// for (uint32_t w = 0; w < dupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }
					// for (uint32_t w = 0; w < rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }

					for (int w = 0; w < dupno; ++w) {
						fprintf(fpdecstr, "%s\n", enstr);
					}
					for (int w = 0; w < rcdupno; ++w) {
						fprintf(fpdecstr, "%s\n", rcstr);
					}
				}
			}
		} else {
			int back_step = atoi(enstr);
			while (back_step > 0) {
				idx = reads[idx].prid;
				--back_step;
			}
		}
	}
	delete[] enstr;

	// cout << "debugcnt: " << debugcnt << endl;
	// cout << "printnum: " << printnum << endl;
	// cout << "dupnum: " << dupnum << endl;
	// fclose(fpw);
	fpdir.close();
	fpisdup.close();
	fpdup.close();
	fpisleft.close();
	// FILE *fpencsg = fopen("decsg.txt", "w");
	while (fscanf(fpenstr, "%s", dstr) != EOF) {
		// cout << "dstr: " << dstr << endl;
		// fprintf(fpencsg, "%s\n", dstr);

		reverseReads(dstr);
		// cout << "dstr: " << dstr << endl;

		fprintf(fpdecstr, "%s\n", dstr);
	}
	fclose(fpenstr);
	// fclose(fpencsg);
	cout << "1111" << endl;

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
		reverseReads(seq);
		fprintf(fpdecstr, "%s\n", seq);
	}
	fpsingleton.close();

	fclose(fpdecstr);
}

// remove singleton reads ID by adding bit stream
// void decompressSingleOrder_work_well_v1_adding_bit(string result) {
void decompressSingleOrder(string result) {
	cout << "decompressSingleOrder(string result)...\n";
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string isleftbinfn = folder + "isleft.bin";
	decompressFile(isleftbinfn);
	std::ifstream fpisleft(isleftbinfn, std::ios::binary);

	string isrootbinfn = folder + "isroot.bin";
	decompressFile(isrootbinfn);
	std::ifstream fpisroot(isrootbinfn, std::ios::binary);

	string orderfn = folder + string("parent.bin");
	decompressFile(orderfn);
	// mstcom::lzma::lzma_decompress((orderfn+".lzma").c_str(), orderfn.c_str());
	std::ifstream fporder(orderfn, std::ios::binary);

	BITPOOL curdirpool, curisleftpool, curisrootpool;
	curdirpool.n = 8;
	curisleftpool.n = 8;
	curisrootpool.n = 8;
	int dir;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");

	char enstr[1<<10], dstr[1<<10], ori[1<<10], aft[1<<10], ctg[1<<10];

	// cout << "max_rid: " << max_rid << endl;
	uint32_t rid = 0, prid, trid;
	char **readsstr = new char*[max_rid];
	bool *isrc = new bool[max_rid];
	bool *direction = new bool[max_rid];
	memset(direction, false, sizeof(bool)*max_rid);
	uint32_t *parent = new uint32_t[max_rid];
	int curshiftdirection;

	uint32_t *newid = new uint32_t[max_rid], newidnum = 0;
	// FILE *fprid2newid = fopen("rid2newid.decomp.txt", "w");

	while (fgets(enstr, 1024, fpenstr) != NULL) {
		enstr[strlen(enstr) - 1] = '\0';
		readsstr[rid] = strdup(enstr);
		int len = strlen(enstr);
		if (len == L) { // root node or singleton reads
			bool isroot = getDir(fpisroot, curisrootpool);
			if (isroot) {
				newid[newidnum] = rid;
				newidnum ++;
				// fprintf(fprid2newid, "%u %u\n", rid, newidnum - 1);
			}
		} else 
		if (len > 0) {
		// if (enstr[0] != '\0') { // not duplicate reads
			newid[newidnum] = rid;
			newidnum ++;
			// fprintf(fprid2newid, "%u %u\n", rid, newidnum - 1);
		}
		++ rid; 
	}
	// fclose(fprid2newid);
	fclose(fpenstr);
	// cout << "newidnum: " << newidnum << endl;
	
	// FILE *fpout = fopen("decorder.txt", "w");
	// while (fgets(enstr, 1024, fpenstr) != NULL) {
		// enstr[strlen(enstr) - 1] = '\0';
		// readsstr[rid] = strdup(enstr);
	FILE *fppid = fopen("pid.decomp.txt", "w");

	for (rid = 0; rid < max_rid; ++rid) {
		parent[rid] = rid;
		int len = strlen(readsstr[rid]);
		// if (readsstr[rid][0] != '0') {
		if (len < L) {
			fporder.read((char*)&trid, sizeof(uint32_t));
			fprintf(fppid, "%u %u\n", newid[trid], trid);
			parent[rid] = newid[trid];
			isrc[rid] = getDir(fpdir, curdirpool);

			if (isupper(readsstr[rid][0])) {
				curshiftdirection = getDir(fpisleft, curisleftpool);
				if (curshiftdirection) direction[rid] = true;
			}
		} 
	}
	// fclose(fpout);
	fclose(fppid);
	fpdir.close();
	fporder.close();
	// cout << readsstr[174540] << endl;
	// cout << "xxxx\n";

	for (rid = 0; rid < max_rid; ++rid) {
		// if (parent[rid] == rid) {
		// 	cout << readsstr[rid] << endl;
		// }
		if (parent[rid] != rid) {
			uint32_v v;
			kv_init(v);

			prid = rid;
			bool debug = false;
			// cout << "prid: " << prid << endl;
			while (parent[prid] != prid) {
				// if (strcmp(readsen[prid], "CTCATAA2G4T2G9AA12G5T20A7ATTC TGATTTGTAGGGCATGACTT") == 0) debug = true;
				// if (prid == 132422) debug = true;
				// if (prid == 42673) debug = true;
				kv_push(uint32_t, v, prid);
				prid = parent[prid];
			}
			for (uint32_t n = v.n - 1; n>= 0; --n) {
				trid = v.a[n];

				curshiftdirection = 0;
				if (direction[trid]) curshiftdirection = 1;
				decodeStr(readsstr[parent[trid]], readsstr[trid], curshiftdirection, dstr); //

				if (isrc[trid]) {
					reverseComplement(dstr);
				}

				readsstr[trid] = strdup(dstr);
				parent[trid] = trid;

				if (n == 0) break;
			}
			kv_destroy(v);
		}
	}

	delete[] parent;
	delete[] direction;
	delete[] newid;
	// cout << "xxxxxxxyyyyyy\n";
	FILE *fpdecstr = fopen(result.c_str(), "w");

	// cout << "yyy" << endl;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		// cout << readsstr[id] << endl;
		reverseReads(readsstr[rid]);
		fprintf(fpdecstr, "%s\n", readsstr[rid]);
		free(readsstr[rid]);
	}
	fclose(fpdecstr);
	delete[] readsstr;
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

	string isleftbinfn = folder + "isleft.bin";
	decompressFile(isleftbinfn);
	std::ifstream fpisleft(isleftbinfn, std::ios::binary);

	string dupfn = folder + "dup.bin";
	decompressFile(dupfn);
	std::ifstream fpdup(dupfn.c_str(), std::ios::binary);

	string isdupfn = folder + "isdup.bin";
	decompressFile(isdupfn);
	std::ifstream fpisdup(isdupfn.c_str(), std::ios::binary);

	string filefn = folder + "file.bin";
	decompressFile(filefn);
	std::ifstream fpfile(filefn.c_str(), std::ios::binary);

	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string distfn = folder + "dist.bin";
	decompressFile(distfn);

	string smdistfn = folder + "smdist.bin";
	decompressFile(smdistfn);

	BITPOOL curdirpool, curisduppool, curisleftpool, curfilepool;
	curdirpool.n = 8;
	curisduppool.n = 8;
	curisleftpool.n = 8;
	curfilepool.n = 8;
	int dir, file;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");

	FILE *fpdecstr1 = fopen((result+"_1").c_str(), "w");
	FILE *fpdecstr2 = fopen((result+"_2").c_str(), "w");
	char str[1<<10], rcstr[1<<10], enstr[1<<10], dstr[1<<10], ori[1<<10], aft[1<<10], ctg[1<<10];

	// int num = 112;
	int dn, idx, curshiftdirection;
	uint32_t dupno, rcdupno;
	vector<READS_t> reads;
	size_t leftnum = 0;
	uint32_t rid = 0;
	char **readsstr = new char*[max_rid + 1];

	while (fscanf(fpenstr, "%s", enstr) != EOF) {
		if (enstr[0] == '-') break;

		if (checkStrOrDig(enstr)) {
			uint32_t enstrlen = strlen(enstr);
			
			if (enstrlen < L) { // 
				// 0 shift > 0; 1 shift < 1;
				if (enstr[0] >= '0' && enstr[0] <= '9') { // no shift
					curshiftdirection = 0;
				} else {
					curshiftdirection = getDir(fpisleft, curisleftpool);
				}
				// 0 for left shift; 1 for right shift;

				decodeStr(reads[idx].str, enstr, curshiftdirection, dstr); //

				dir = getDir(fpdir, curdirpool);
				if (dir) {
					reverseComplement(dstr);
				}

				READS_t r(dstr, idx);
				reads.push_back(r);
				idx = reads.size() - 1;

				// fprintf(fpw, "%s\n", dstr);
				string tempstr = dstr;

				reverseReads(dstr);

				// fprintf(fpdecstr, "%s\n", dstr);
				readsstr[rid ++] = strdup(dstr);
				// ++ printnum;
				// if (strcmp(dstr, "83GGAAATGGCTGGGTCAAATGGTATTTCTAGTTCTAGATCCCTGAGGAATCGCCACACTGACTTCCACAATGGTTGAACTAGTTTACATTCCCACCAAGA") == 0) 
				// 	cout << "dstr in !isctgtree\n";

				int isdup = getDir(fpisdup, curisduppool);

				if (isdup) {
					reverseReads(dstr);
					strcpy(rcstr, dstr);
					reverseComplement(rcstr);

					reverseReads(dstr);
					reverseReads(rcstr);

					fpdup.read((char*)&dupno, sizeof(uint32_t));
					fpdup.read((char*)&rcdupno, sizeof(uint32_t));

					// for (int w = 0; w < dupno + rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// }
					// for (uint32_t w = 0; w < dupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }
					// for (uint32_t w = 0; w < rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }

					for (int w = 0; w < dupno; ++w) {
						// fprintf(fpdecstr, "%s\n", dstr);
						readsstr[rid ++] = strdup(dstr);
					}
					for (int w = 0; w < rcdupno; ++w) {
						// fprintf(fpdecstr, "%s\n", rcstr);
						readsstr[rid ++] = strdup(rcstr);
					}
				}

			} else { // root node // enstrlen == L
				reads.clear();

				READS_t r(enstr, 0);
				reads.push_back(r);
				idx = 0;

				// preshiftdirection = 0; // 0 for left shift; 1 for right shift;
				// fprintf(fpw, "%s\n", enstr);

				string tempstr = enstr;

				reverseReads(enstr);

				// fprintf(fpdecstr, "%s\n", enstr);
				readsstr[rid ++] = strdup(enstr);
				// ++ printnum;

				bool debug = false;
				int isdup = getDir(fpisdup, curisduppool);

				if (isdup) {
					reverseReads(enstr);
					strcpy(rcstr, enstr);
					reverseComplement(rcstr);

					reverseReads(enstr);
					reverseReads(rcstr);

					fpdup.read((char*)&dupno, sizeof(uint32_t));
					fpdup.read((char*)&rcdupno, sizeof(uint32_t));

					// for (int w = 0; w < dupno + rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// }
					// for (uint32_t w = 0; w < dupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }
					// for (uint32_t w = 0; w < rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }

					for (int w = 0; w < dupno; ++w) {
						// fprintf(fpdecstr, "%s\n", enstr);
						readsstr[rid ++] = strdup(enstr);
					}
					for (int w = 0; w < rcdupno; ++w) {
						// fprintf(fpdecstr, "%s\n", rcstr);
						readsstr[rid ++] = strdup(rcstr);
					}
				}
			}
		} else {
			int back_step = atoi(enstr);
			while (back_step > 0) {
				idx = reads[idx].prid;
				--back_step;
			}
		}
	}

	
	while (fscanf(fpenstr, "%s", dstr) != EOF) {

		reverseReads(dstr);

		readsstr[rid ++] = strdup(dstr);
	}
	fclose(fpenstr);

	fpdir.close();
	fpisdup.close();
	fpdup.close();
	fpisleft.close();

	///
	size_t *ids = new size_t[max_rid + 3];

	size_t *tr, M;
	build(tr, M, max_rid);

	bool *flag = new bool[max_rid + 3];
	memset(flag, false, sizeof(bool)*max_rid); 

	size_t half = max_rid >> 1, pid;
	std::ifstream fpdist(distfn.c_str(), std::ios::binary);
	std::ifstream smfpdist(smdistfn.c_str(), std::ios::binary);
	
	uint32_t tempid;
	uint32_t dist;
	// int32_t dist;
	uint16_t disttemp;

	size_t id = 0;
	for (size_t num = 0; num < max_rid; ++num) {
		if (!flag[num]) {
			dir = getDir(fpfile, curfilepool);
			// fpdist.read((char*)&dist, sizeof(int32_t));
			// fpdist.read((char*)&dist, sizeof(uint32_t));
			// cout << dist << endl;
			file = getDir(fpfile, curfilepool);
			if (file) {
				fpdist.read((char*)&dist, sizeof(uint32_t));
			} else {
				smfpdist.read((char*)&disttemp, sizeof(uint16_t));
				dist = (uint32_t)disttemp;
			}
			
			ids[num] = id ++;
			if (dir) {
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

	// cout << "yyyy\n";
	fpdist.close();
	smfpdist.close();
	delete[] tr;
	delete[] flag;

	// FILE *fp = fopen("decids.txt", "w");
	// for (size_t id = 0; id < max_rid; ++id) {
	// 	fprintf(fp, "%u\n", ids[id]);
	// }
	// fclose(fp);

	size_t *od = new size_t[max_rid + 3];
	for (size_t id = 0; id < max_rid; ++id) {
		od[ids[id]] = id;
	}
	delete[] ids;

	// fp = fopen("decod.txt", "w");
	// for (size_t id = 0; id < max_rid; ++id) {
	// 	fprintf(fp, "%u\n", od[id]);
	// }
	// fclose(fp);

	for (size_t id = 0; id < half; ++id) {
		fprintf(fpdecstr1, "%s\n", readsstr[od[id]]);
	}

	for (size_t id = half; id < max_rid; ++id) {
		fprintf(fpdecstr2, "%s\n", readsstr[od[id]]);
	}

	fclose(fpdecstr1);
	fclose(fpdecstr2);
}

void decompressPEXL1L2(string result) {
	fprintf(stderr, "decompressPE()...\n");
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string isleftbinfn = folder + "isleft.bin";
	decompressFile(isleftbinfn);
	std::ifstream fpisleft(isleftbinfn, std::ios::binary);

	string dupfn = folder + "dup.bin";
	decompressFile(dupfn);
	std::ifstream fpdup(dupfn.c_str(), std::ios::binary);

	string isdupfn = folder + "isdup.bin";
	decompressFile(isdupfn);
	std::ifstream fpisdup(isdupfn.c_str(), std::ios::binary);

	string filefn = folder + "file.bin";
	decompressFile(filefn);
	std::ifstream fpfile(filefn.c_str(), std::ios::binary);

	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string distfn = folder + "dist.bin";
	decompressFile(distfn);

	string smdistfn = folder + "smdist.bin";
	decompressFile(smdistfn);

	BITPOOL curdirpool, curisduppool, curisleftpool, curfilepool;
	curdirpool.n = 8;
	curisduppool.n = 8;
	curisleftpool.n = 8;
	curfilepool.n = 8;
	int dir, file;
	int *filelist = new int[max_rid];

	///
	size_t *ids = new size_t[max_rid + 3];

	size_t *tr, M;
	build(tr, M, max_rid);

	bool *flag = new bool[max_rid + 3];
	memset(flag, false, sizeof(bool)*max_rid); 

	size_t half = max_rid >> 1, pid;
	std::ifstream fpdist(distfn.c_str(), std::ios::binary);
	std::ifstream smfpdist(smdistfn.c_str(), std::ios::binary);
	
	uint32_t tempid;
	uint32_t dist;
	// int32_t dist;
	uint16_t disttemp;

	size_t id = 0;
	for (size_t num = 0; num < max_rid; ++num) {
		if (!flag[num]) {
			dir = getDir(fpfile, curfilepool);
			filelist[num] = dir;
			file = getDir(fpfile, curfilepool);
			if (file) {
				fpdist.read((char*)&dist, sizeof(uint32_t));
			} else {
				smfpdist.read((char*)&disttemp, sizeof(uint16_t));
				dist = (uint32_t)disttemp;
			}
			
			ids[num] = id ++;
			if (dir) {
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
			filelist[pid - 1] = dir > 0?0:1;
		}
	}

	// cout << "yyyy\n";
	fpdist.close();
	smfpdist.close();
	delete[] tr;
	delete[] flag;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");

	FILE *fpdecstr1 = fopen((result+"_1").c_str(), "w");
	FILE *fpdecstr2 = fopen((result+"_2").c_str(), "w");
	char str[1<<10], rcstr[1<<10], enstr[1<<10], dstr[1<<10], ori[1<<10], aft[1<<10], ctg[1<<10];

	// int num = 112;
	int dn, idx, curshiftdirection;
	uint32_t dupno, rcdupno;
	vector<READS_t> reads;
	size_t leftnum = 0;
	uint32_t rid = 0;
	char **readsstr = new char*[max_rid + 1];

	while (fscanf(fpenstr, "%s", enstr) != EOF) {
		if (enstr[0] == '-') break;

		if (checkStrOrDig(enstr)) {
			uint32_t enstrlen = strlen(enstr);
			
			if (enstrlen < min(L1, L2)) { // 
				// 0 shift > 0; 1 shift < 1;
				if (enstr[0] >= '0' && enstr[0] <= '9') { // no shift
					curshiftdirection = 0;
				} else {
					curshiftdirection = getDir(fpisleft, curisleftpool);
				}
				// 0 for left shift; 1 for right shift;
				// cout << "1010101010" << endl;
				decodeStr(reads[idx].str, enstr, curshiftdirection, filelist[rid], dstr); //
				// cout << "12121212121212" << endl;

				int dstrlen = (filelist[rid] == 0? L1: L2);
				dir = getDir(fpdir, curdirpool);
				if (dir) {
					reverseComplement(dstr, dstrlen);
				}

				READS_t r(dstr, idx);
				reads.push_back(r);
				idx = reads.size() - 1;

				// fprintf(fpw, "%s\n", dstr);
				string tempstr = dstr;

				reverseReads(dstr, dstrlen);

				// fprintf(fpdecstr, "%s\n", dstr);
				readsstr[rid ++] = strdup(dstr);
				// ++ printnum;
				// if (strcmp(dstr, "83GGAAATGGCTGGGTCAAATGGTATTTCTAGTTCTAGATCCCTGAGGAATCGCCACACTGACTTCCACAATGGTTGAACTAGTTTACATTCCCACCAAGA") == 0) 
				// 	cout << "dstr in !isctgtree\n";

				int isdup = getDir(fpisdup, curisduppool);

				if (isdup) {
					reverseReads(dstr, dstrlen);
					strcpy(rcstr, dstr);
					reverseComplement(rcstr, dstrlen);

					reverseReads(dstr, dstrlen);
					reverseReads(rcstr, dstrlen);

					fpdup.read((char*)&dupno, sizeof(uint32_t));
					fpdup.read((char*)&rcdupno, sizeof(uint32_t));

					for (int w = 0; w < dupno; ++w) {
						// fprintf(fpdecstr, "%s\n", dstr);
						readsstr[rid ++] = strdup(dstr);
					}
					for (int w = 0; w < rcdupno; ++w) {
						// fprintf(fpdecstr, "%s\n", rcstr);
						readsstr[rid ++] = strdup(rcstr);
					}
				}

			} else {			
			// if (enstrlen == L1 || enstrlen == L2) { // root node // enstrlen == L1 || == L2
				reads.clear();

				READS_t r(enstr, 0);
				reads.push_back(r);
				idx = 0;

				// preshiftdirection = 0; // 0 for left shift; 1 for right shift;
				// fprintf(fpw, "%s\n", enstr);

				string tempstr = enstr;
				int enstrlen = strlen(enstr);
				reverseReads(enstr, enstrlen);

				// fprintf(fpdecstr, "%s\n", enstr);
				readsstr[rid ++] = strdup(enstr);
				// ++ printnum;

				bool debug = false;
				int isdup = getDir(fpisdup, curisduppool);

				if (isdup) {
					reverseReads(enstr, enstrlen);
					strcpy(rcstr, enstr);
					reverseComplement(rcstr, enstrlen);

					reverseReads(enstr, enstrlen);
					reverseReads(rcstr, enstrlen);

					fpdup.read((char*)&dupno, sizeof(uint32_t));
					fpdup.read((char*)&rcdupno, sizeof(uint32_t));

					for (int w = 0; w < dupno; ++w) {
						// fprintf(fpdecstr, "%s\n", enstr);
						readsstr[rid ++] = strdup(enstr);
					}
					for (int w = 0; w < rcdupno; ++w) {
						// fprintf(fpdecstr, "%s\n", rcstr);
						readsstr[rid ++] = strdup(rcstr);
					}
				}
			}
		} else {
			int back_step = atoi(enstr);
			while (back_step > 0) {
				idx = reads[idx].prid;
				--back_step;
			}
		}
	}
	// cout << "over decode..." << endl;
	
	while (fscanf(fpenstr, "%s", dstr) != EOF) {

		reverseReads(dstr, strlen(dstr));

		readsstr[rid ++] = strdup(dstr);
	}
	fclose(fpenstr);
	// cout << "1111" << endl;
	fpdir.close();
	fpisdup.close();
	fpdup.close();
	fpisleft.close();


	// FILE *fp = fopen("decids.txt", "w");
	// for (size_t id = 0; id < max_rid; ++id) {
	// 	fprintf(fp, "%u\n", ids[id]);
	// }
	// fclose(fp);

	// cout << "11112222" << endl;
	size_t *od = new size_t[max_rid + 3];
	for (size_t id = 0; id < max_rid; ++id) {
		od[ids[id]] = id;
	}
	delete[] ids;
	// cout << "111122223333333" << endl;

	// fp = fopen("decod.txt", "w");
	// for (size_t id = 0; id < max_rid; ++id) {
	// 	fprintf(fp, "%u\n", od[id]);
	// }
	// fclose(fp);

	// cout << "1111222233333331111" << endl;
	// fp = fopen("decstr.txt", "w");
	// for (size_t id = 0; id < max_rid; ++id) {
	// 	fprintf(fp, "%s\n", readsstr[id]);
	// }
	// fclose(fp);

	for (size_t id = 0; id < half; ++id) {
		fprintf(fpdecstr1, "%s\n", readsstr[od[id]]);
	}
	// cout << "1111222233333334444" << endl;

	for (size_t id = half; id < max_rid; ++id) {
		// cout << readsstr[od[id]] << endl;
		fprintf(fpdecstr2, "%s\n", readsstr[od[id]]);
	}
	// cout << "1111222233333334444555" << endl;

	fclose(fpdecstr1);
	fclose(fpdecstr2);
}

// void decompressPEOrder_work_well(string result) {
void decompressPEOrder(string result) {
	fprintf(stderr, "decompressPEOrder()...\n");
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string orderfn = folder + string("order.bin");
	decompressFile(orderfn);
	std::ifstream fporder(orderfn, std::ios::binary);

	string isleftbinfn = folder + "isleft.bin";
	decompressFile(isleftbinfn);
	std::ifstream fpisleft(isleftbinfn, std::ios::binary);

	string dupfn = folder + "dup.bin";
	decompressFile(dupfn);
	std::ifstream fpdup(dupfn.c_str(), std::ios::binary);

	string isdupfn = folder + "isdup.bin";
	decompressFile(isdupfn);
	std::ifstream fpisdup(isdupfn.c_str(), std::ios::binary);

	string filefn = folder + "file.bin";
	decompressFile(filefn);
	std::ifstream fpfile(filefn.c_str(), std::ios::binary);

	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string distfn = folder + "dist.bin";
	decompressFile(distfn);

	string smdistfn = folder + "smdist.bin";
	decompressFile(smdistfn);

	BITPOOL curdirpool, curisduppool, curisleftpool, curfilepool;
	curdirpool.n = 8;
	curisduppool.n = 8;
	curisleftpool.n = 8;
	curfilepool.n = 8;
	int dir, file;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");

	FILE *fpdecstr1 = fopen((result+"_1").c_str(), "w");
	FILE *fpdecstr2 = fopen((result+"_2").c_str(), "w");
	char str[1<<10], rcstr[1<<10], enstr[1<<10], dstr[1<<10], ori[1<<10], aft[1<<10], ctg[1<<10];

	// int num = 112;
	int dn, idx, curshiftdirection;
	uint32_t dupno, rcdupno;
	vector<READS_t> reads;
	size_t leftnum = 0;
	uint32_t rid = 0;
	char **readsstr = new char*[max_rid + 1];

	while (fscanf(fpenstr, "%s", enstr) != EOF) {
		if (enstr[0] == '-') break;

		if (checkStrOrDig(enstr)) {
			uint32_t enstrlen = strlen(enstr);
			
			if (enstrlen < L) { // 
				// 0 shift > 0; 1 shift < 1;
				if (enstr[0] >= '0' && enstr[0] <= '9') { // no shift
					curshiftdirection = 0;
				} else {
					curshiftdirection = getDir(fpisleft, curisleftpool);
				}
				// 0 for left shift; 1 for right shift;

				decodeStr(reads[idx].str, enstr, curshiftdirection, dstr); //

				dir = getDir(fpdir, curdirpool);
				if (dir) {
					reverseComplement(dstr);
				}

				READS_t r(dstr, idx);
				reads.push_back(r);
				idx = reads.size() - 1;

				// fprintf(fpw, "%s\n", dstr);
				string tempstr = dstr;

				reverseReads(dstr);

				// fprintf(fpdecstr, "%s\n", dstr);
				readsstr[rid ++] = strdup(dstr);
				// ++ printnum;
				// if (strcmp(dstr, "83GGAAATGGCTGGGTCAAATGGTATTTCTAGTTCTAGATCCCTGAGGAATCGCCACACTGACTTCCACAATGGTTGAACTAGTTTACATTCCCACCAAGA") == 0) 
				// 	cout << "dstr in !isctgtree\n";

				int isdup = getDir(fpisdup, curisduppool);

				if (isdup) {
					reverseReads(dstr);
					strcpy(rcstr, dstr);
					reverseComplement(rcstr);

					reverseReads(dstr);
					reverseReads(rcstr);

					fpdup.read((char*)&dupno, sizeof(uint32_t));
					fpdup.read((char*)&rcdupno, sizeof(uint32_t));

					// for (int w = 0; w < dupno + rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// }
					// for (uint32_t w = 0; w < dupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }
					// for (uint32_t w = 0; w < rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }

					for (int w = 0; w < dupno; ++w) {
						// fprintf(fpdecstr, "%s\n", dstr);
						readsstr[rid ++] = strdup(dstr);
					}
					for (int w = 0; w < rcdupno; ++w) {
						// fprintf(fpdecstr, "%s\n", rcstr);
						readsstr[rid ++] = strdup(rcstr);
					}
				}

			} else { // root node // enstrlen == L
				reads.clear();

				READS_t r(enstr, 0);
				reads.push_back(r);
				idx = 0;

				// preshiftdirection = 0; // 0 for left shift; 1 for right shift;
				// fprintf(fpw, "%s\n", enstr);

				string tempstr = enstr;

				reverseReads(enstr);

				// fprintf(fpdecstr, "%s\n", enstr);
				readsstr[rid ++] = strdup(enstr);
				// ++ printnum;

				bool debug = false;
				int isdup = getDir(fpisdup, curisduppool);

				if (isdup) {
					reverseReads(enstr);
					strcpy(rcstr, enstr);
					reverseComplement(rcstr);

					reverseReads(enstr);
					reverseReads(rcstr);

					fpdup.read((char*)&dupno, sizeof(uint32_t));
					fpdup.read((char*)&rcdupno, sizeof(uint32_t));

					// for (int w = 0; w < dupno + rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// }
					// for (uint32_t w = 0; w < dupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }
					// for (uint32_t w = 0; w < rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }

					for (int w = 0; w < dupno; ++w) {
						// fprintf(fpdecstr, "%s\n", enstr);
						readsstr[rid ++] = strdup(enstr);
					}
					for (int w = 0; w < rcdupno; ++w) {
						// fprintf(fpdecstr, "%s\n", rcstr);
						readsstr[rid ++] = strdup(rcstr);
					}
				}
			}
		} else {
			int back_step = atoi(enstr);
			while (back_step > 0) {
				idx = reads[idx].prid;
				--back_step;
			}
		}
	}
	
	while (fscanf(fpenstr, "%s", dstr) != EOF) {
		reverseReads(dstr);
		readsstr[rid ++] = strdup(dstr);
	}
	fclose(fpenstr);

	fpdir.close();
	fpisdup.close();
	fpdup.close();
	fpisleft.close();

	///
	size_t *ids = new size_t[max_rid + 3];

	size_t *tr, M;
	build(tr, M, max_rid);

	bool *flag = new bool[max_rid + 3];
	memset(flag, false, sizeof(bool)*max_rid); 

	size_t half = max_rid >> 1, pid;
	std::ifstream fpdist(distfn.c_str(), std::ios::binary);
	std::ifstream smfpdist(smdistfn.c_str(), std::ios::binary);
	
	uint32_t tempid;
	uint32_t dist;
	// int32_t dist;
	uint16_t disttemp;

	for (size_t num = 0; num < max_rid; ++num) {
		if (!flag[num]) {
			fporder.read((char*)&tempid, sizeof(uint32_t));

			dir = getDir(fpfile, curfilepool);
			file = getDir(fpfile, curfilepool);
			if (file) {
				fpdist.read((char*)&dist, sizeof(uint32_t));
			} else {
				smfpdist.read((char*)&disttemp, sizeof(uint16_t));
				dist = (uint32_t)disttemp;
			}
			++ dist;
			// fpdist.read((char*)&dist, sizeof(int32_t));

			ids[num] = tempid;
			if (dir) {
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
	fporder.close();
	// cout << "yyyy\n";
	fpdist.close();
	smfpdist.close();
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

	fclose(fpdecstr1);
	fclose(fpdecstr2);
}

void decompressPEOrderL1L2(string result) {
	fprintf(stderr, "decompressPEOrder()...\n");
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string orderfn = folder + string("order.bin");
	decompressFile(orderfn);
	std::ifstream fporder(orderfn, std::ios::binary);

	string isleftbinfn = folder + "isleft.bin";
	decompressFile(isleftbinfn);
	std::ifstream fpisleft(isleftbinfn, std::ios::binary);

	string dupfn = folder + "dup.bin";
	decompressFile(dupfn);
	std::ifstream fpdup(dupfn.c_str(), std::ios::binary);

	string isdupfn = folder + "isdup.bin";
	decompressFile(isdupfn);
	std::ifstream fpisdup(isdupfn.c_str(), std::ios::binary);

	string filefn = folder + "file.bin";
	decompressFile(filefn);
	std::ifstream fpfile(filefn.c_str(), std::ios::binary);

	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string distfn = folder + "dist.bin";
	decompressFile(distfn);

	string smdistfn = folder + "smdist.bin";
	decompressFile(smdistfn);

	BITPOOL curdirpool, curisduppool, curisleftpool, curfilepool;
	curdirpool.n = 8;
	curisduppool.n = 8;
	curisleftpool.n = 8;
	curfilepool.n = 8;
	int dir, file;
	int *filelist = new int[max_rid];

	///
	size_t *ids = new size_t[max_rid + 3];

	size_t *tr, M;
	build(tr, M, max_rid);

	bool *flag = new bool[max_rid + 3];
	memset(flag, false, sizeof(bool)*max_rid); 

	size_t half = max_rid >> 1, pid;
	std::ifstream fpdist(distfn.c_str(), std::ios::binary);
	std::ifstream smfpdist(smdistfn.c_str(), std::ios::binary);
	
	uint32_t tempid;
	uint32_t dist;
	// int32_t dist;
	uint16_t disttemp;

	for (size_t num = 0; num < max_rid; ++num) {
		if (!flag[num]) {
			fporder.read((char*)&tempid, sizeof(uint32_t));
			dir = getDir(fpfile, curfilepool);
			filelist[num] = dir;
			file = getDir(fpfile, curfilepool);
			if (file) {
				fpdist.read((char*)&dist, sizeof(uint32_t));
			} else {
				smfpdist.read((char*)&disttemp, sizeof(uint16_t));
				dist = (uint32_t)disttemp;
			}
			++ dist;
			// fpdist.read((char*)&dist, sizeof(int32_t));

			ids[num] = tempid;
			if (dir) {
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
			filelist[pid - 1] = dir > 0?0:1;
		}
	}
	fporder.close();
	// cout << "yyyy\n";
	fpdist.close();
	smfpdist.close();
	delete[] tr;
	delete[] flag;

	size_t *od = new size_t[max_rid + 3];
	for (size_t id = 0; id < max_rid; ++id) {
		od[ids[id]] = id;
	}
	delete[] ids;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");

	FILE *fpdecstr1 = fopen((result+"_1").c_str(), "w");
	FILE *fpdecstr2 = fopen((result+"_2").c_str(), "w");
	char str[1<<10], rcstr[1<<10], enstr[1<<10], dstr[1<<10], ori[1<<10], aft[1<<10], ctg[1<<10];

	// int num = 112;
	int dn, idx, curshiftdirection;
	uint32_t dupno, rcdupno;
	vector<READS_t> reads;
	size_t leftnum = 0;
	uint32_t rid = 0;
	char **readsstr = new char*[max_rid + 1];

	// FILE *fpdirinf = fopen("decdir.txt", "w");

	while (fscanf(fpenstr, "%s", enstr) != EOF) {
		if (enstr[0] == '-') break;

		if (checkStrOrDig(enstr)) {
			uint32_t enstrlen = strlen(enstr);
			
			if (enstrlen < min(L1, L2)) { // 
				// 0 shift < 0; 1 shift > 1;
				if (enstr[0] >= '0' && enstr[0] <= '9') { // no shift
					curshiftdirection = 0;
					// fprintf(fpdirinf, "n\n");
					// cout << enstr << endl;
					// exit(0);
				} else {
					curshiftdirection = getDir(fpisleft, curisleftpool);
					// fprintf(fpdirinf, "%d\n", curshiftdirection);
				}
				// 0 for left shift; 1 for right shift;

				// if (strcmp(enstr, "2G10C6A16C1G2T10C1G12G3G3A0C1G0A3A1G1A7A0C") == 0) {
				// 	cout << "curshiftdirection: " << curshiftdirection << endl;
				// }

				decodeStr(reads[idx].str, enstr, curshiftdirection, filelist[rid], dstr); //

				int dstrlen = (filelist[rid] == 0? L1: L2);
				dir = getDir(fpdir, curdirpool);
				if (dir) {
					reverseComplement(dstr, dstrlen);
				}

				READS_t r(dstr, idx);
				reads.push_back(r);
				idx = reads.size() - 1;

				// fprintf(fpw, "%s\n", dstr);
				string tempstr = dstr;

				reverseReads(dstr, dstrlen);

				// fprintf(fpdecstr, "%s\n", dstr);
				readsstr[rid ++] = strdup(dstr);
				// ++ printnum;
				// if (strcmp(enstr, "ACTCCGTTAAGAGGACGGAGTCAGAGGGTT3T1N14G5A6C11N0N8T0A5G6G") == 0) {
				// 	cout << dstr << endl;
				// 	exit(0);
				// }

				int isdup = getDir(fpisdup, curisduppool);

				if (isdup) {
					reverseReads(dstr, dstrlen);
					strcpy(rcstr, dstr);
					reverseComplement(rcstr, dstrlen);

					reverseReads(dstr, dstrlen);
					reverseReads(rcstr, dstrlen);

					fpdup.read((char*)&dupno, sizeof(uint32_t));
					fpdup.read((char*)&rcdupno, sizeof(uint32_t));

					// for (int w = 0; w < dupno + rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// }
					// for (uint32_t w = 0; w < dupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }
					// for (uint32_t w = 0; w < rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }

					for (int w = 0; w < dupno; ++w) {
						// fprintf(fpdecstr, "%s\n", dstr);
						readsstr[rid ++] = strdup(dstr);
					}
					for (int w = 0; w < rcdupno; ++w) {
						// fprintf(fpdecstr, "%s\n", rcstr);
						readsstr[rid ++] = strdup(rcstr);
					}
				}

			} else { // root node // enstrlen == L
				reads.clear();

				READS_t r(enstr, 0);
				reads.push_back(r);
				idx = 0;

				// preshiftdirection = 0; // 0 for left shift; 1 for right shift;
				// fprintf(fpw, "%s\n", enstr);

				string tempstr = enstr;
				int enstrlen = strlen(enstr);
				reverseReads(enstr, enstrlen);

				// fprintf(fpdecstr, "%s\n", enstr);
				readsstr[rid ++] = strdup(enstr);
				// ++ printnum;

				bool debug = false;
				int isdup = getDir(fpisdup, curisduppool);

				if (isdup) {
					reverseReads(enstr, enstrlen);
					strcpy(rcstr, enstr);
					reverseComplement(rcstr, enstrlen);

					reverseReads(enstr, enstrlen);
					reverseReads(rcstr, enstrlen);

					fpdup.read((char*)&dupno, sizeof(uint32_t));
					fpdup.read((char*)&rcdupno, sizeof(uint32_t));

					// for (int w = 0; w < dupno + rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// }
					// for (uint32_t w = 0; w < dupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }
					// for (uint32_t w = 0; w < rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", tempstr.c_str());
					// 	++ printnum;
					// 	++ dupnum;
					// }

					for (int w = 0; w < dupno; ++w) {
						// fprintf(fpdecstr, "%s\n", enstr);
						readsstr[rid ++] = strdup(enstr);
					}
					for (int w = 0; w < rcdupno; ++w) {
						// fprintf(fpdecstr, "%s\n", rcstr);
						readsstr[rid ++] = strdup(rcstr);
					}
				}
			}
		} else {
			int back_step = atoi(enstr);
			while (back_step > 0) {
				idx = reads[idx].prid;
				--back_step;
			}
		}
	}
	
	while (fscanf(fpenstr, "%s", dstr) != EOF) {
		// if (strcmp(dstr, "CTCTAAGAGAGAAAAAATAAGAGAGGAAAAGAGAGAAAAGAAAACAACCTCCTTTTACCTACTATTATCATTCATTTATTGTCGTTAAAAGAANNNNTCTAC") == 0) {
		// 	cout << "it is here" << endl;
		// }
		reverseReads(dstr, strlen(dstr));
		readsstr[rid ++] = strdup(dstr);
	}
	fclose(fpenstr);

	fpdir.close();
	fpisdup.close();
	fpdup.close();
	fpisleft.close();
	// fclose(fpdirinf);

	for (size_t id = 0; id < half; ++id) {
		fprintf(fpdecstr1, "%s\n", readsstr[od[id]]);
	}

	for (size_t id = half; id < max_rid; ++id) {
		fprintf(fpdecstr2, "%s\n", readsstr[od[id]]);
	}

	fclose(fpdecstr1);
	fclose(fpdecstr2);
}

void getPars(int argc, char* argv[]) {
	bool isinfile = false, isoutfile = false; 
	int oc;
	while ((oc = getopt(argc, argv, "i:o:")) >= 0) {
		switch (oc) {
			case 'i':
				infile = optarg;
				isinfile = true;
				break;
			case 'o':
				outfile = optarg;
				isoutfile = true;
				break;
			case 'h':
				show_usage(argv[0]);
				exit(0);
			case '?':
				std::cerr << "Error parameters.\n Please run 'mstcom -h'\n";
				exit(1);
				break;
		}
	}

	if (!isinfile || !isoutfile) {
		fprintf(stderr, "Required parameters (input or output file) are not provided!!\n\n");
		exit(1);
	}
	
	std::ifstream f;

	f.open(infile);
	if (f.fail()) {
		fprintf(stderr, "Input file '%s' does not exist.\n", infile.c_str());
		exit(1);
	}
	f.close();

	// f.open(outfile);
	// if (f.is_open()) {
	// 	fprintf(stderr, "The output file '%s' exist.\n", outfile.c_str());
	// 	exit(1);
	// }
	// f.close();

	folder = generateString("mstcom_dec", 10); //creat a temp folder for current input
	// string cmd = "mkdir -p " + folder;
	// system(cmd.c_str());
	mode_t mode = 0777;
	int nError = mkdir(folder.c_str(), mode);
	if (nError != 0) {
		fprintf(stderr, "Failed to creat the folder '%s'\n", folder.c_str());
		exit(EXIT_FAILURE);
	}
	folder += "/";
}

public:
	int dec_class_main(int argc, char *argv[]) {
		getPars(argc, argv);
		init();

		// decompressSingleTest(outfile.c_str());
		// return 0;
 
		string cmd = "tar -xf " + infile + " -C " + folder;
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
		L1 = L2 = L;
		if (ispe) {
			fscanf(fp, "%d%d", &L1, &L2);
			fprintf(stderr, "L1: %d; L2: %d\n", L1, L2);
		}
		fclose(fp);
		// return 0;

		if (ispe) {
			// decompressPE(argv[2]);
			if (isorder) {
				if (L1 == L2) decompressPEOrder(outfile.c_str());
				else decompressPEOrderL1L2(outfile.c_str());
			} else {
				if (L1 == L2) decompressPEX(outfile.c_str());
				else decompressPEXL1L2(outfile.c_str());
			}
		} else {
			if (isorder) {
				decompressSingleOrder(outfile.c_str());
			} else 
#ifdef DFSENCODING			
				decompressSingleDFS(outfile.c_str());
#endif

#ifndef DFSENCODING				
				decompressSingleV1(outfile.c_str());
#endif
		}

		cmd = "rm -rf " + folder;
		system(cmd.c_str());
		
		return 0;
	}

};

int decompress_main(int argc, char *argv[]) {
	dec_class dec;
	dec.dec_class_main(argc, argv);
	return 0;
}
