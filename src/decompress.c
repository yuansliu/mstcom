# include "mstcom.h"
using namespace std;

class dec_class{
int L;
int ispe, isorder;
size_t max_rid;
const int nthreads = 24;
string folder, infile, outfile;

struct READS_t {
	string str;
	size_t prid, prechildrid; //reads id of parent
	READS_t(): prechildrid(0) {}
	READS_t(char *_str, const size_t &_prid): str(_str), prid(_prid), prechildrid(0) {}
	void set(char *_str) { 
		str = _str;
	}
};

char complement[256];

const char invert_code_rule[4] = {'A', 'C', 'G', 'T'}; //decoding rule //A - 0; C - 1; G - 2; T - 3;

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
		// if (strcmp(enstr, "88CTCGAGACAGAGT\n") == 0) debug = true;
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
					cout << "L: " << L << endl;
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
		if ((str[i] >= 'A' && str[i] <= 'Z') || str[i] == '$') {
			return true; // string
		}
	}
	return false; // dig
}

void split2TwoEncodedStr(char *str, char *ori, char *aft) {
	char *s = str;
	*aft = '\0';
	for (; (*ori = *s) != '\0'; ++ori, ++s) {
		if (*s == '|') {
			*ori = '\0';
			++ s;
			for (; (*aft = *s) != '\0'; ++aft, ++s);
			break;
		}
	}
}

int getDup(char *str) {
	char *s = str;
	int res = 0;
	for (; *s != '\0'; ++s) {
		if(*s == '$') {
			if (*(s+1) >= '0' && *(s+1) <= '9') {
				res = atoi(s + 1);
			} else {
				res = 1;
			}
			*s = '\0';
			break;
		}
	}
	return res;
}

int mainxxx(int argc, char const *argv[]) {
	string parent = "AAAGGCCCCAGTTTGGCAGACCGACAGCGTGAATACCTTTTAGACATGATCCCTCCCCGGTCTATATCGCAGTCCATCAGTGGACAGAAATAACGCCTGTT";
	char enstr[1<<10], dstr[1<<10];
	strcpy(enstr, "C AT0\n");
	decodeStr(parent, enstr, dstr);
	fprintf(stderr, "%s\n", dstr);
	return 0;
}

void decompressSingle_(string result) {
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string rootstrfn = folder + "rootstr.txt";
	mstcom::bsc::BSC_decompress((rootstrfn+".bsc").c_str(), rootstrfn.c_str());
	string orderfn;

	std::ifstream fporder;
	if (isorder) {
		orderfn = folder + string("order.bin");
		// mstcom::bsc::BSC_decompress((orderfn+".bsc").c_str(), orderfn.c_str());
		mstcom::lzma::lzma_decompress((orderfn+".lzma").c_str(), orderfn.c_str());
		fporder.open(orderfn, std::ios::binary);
	}
	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	string singletonfn = folder + "singleton.bin"; 
	mstcom::bsc::BSC_decompress((singletonfn+".bsc").c_str(), singletonfn.c_str());
	// string readsnfn = folder + "readsN.txt";
	// mstcom::bsc::BSC_decompress((readsnfn+".bsc").c_str(), readsnfn.c_str());

	// exit(0);
	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	BITPOOL curdirpool;
	curdirpool.n = 8;
	int dir, changed;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");
	FILE *fproot = fopen(rootstrfn.c_str(), "r");
	// sprintf(name, "rootstr.bin");
	// std::ifstream fproot(name, std::ios::binary);

	FILE *fpdecstr = fopen(result.c_str(), "w");
	char str[1<<10], rcstr[1<<10], enstr[1<<10], dstr[1<<10], ori[1<<10], aft[1<<10], ctg[1<<10];

	// int num = 112;
	int dn;
	// DNAPOOL rootdnapool;
	// rootdnapool.n = 4;
	if (isorder) {
		cout << "max_rid: " << max_rid << endl;
		uint32_t rid, prid, prerootrid = 0;
		char **readsstr = new char*[max_rid];
		bool *flag = new bool[max_rid];
		memset(flag, false, sizeof(bool)*max_rid);

		while (fscanf(fproot, "%s", str) != EOF) {
			if (str[0] == '-') break;

			dn = getDup(str);

			split2TwoEncodedStr(str, ori, aft);
			vector<READS_t> reads;
			READS_t r(ori, 0);
			if (strlen(aft) > 0) {
				decodeStr(ori, aft, ctg);
				r.set(ctg);
			}
			reads.push_back(r);
			int idx = 0;

			fporder.read((char*)&rid, sizeof(uint32_t));
			rid += prerootrid;
			readsstr[rid] = strdup(ori);
			flag[rid] = true;

			prerootrid = rid;

			if (dn > 0) {
				strcpy(rcstr, ori);
				reverseComplement(rcstr);

				prid = rid;
				for (int w = 0; w < dn; ++w) {
					fporder.read((char*)&rid, sizeof(uint32_t));
					rid += prid;

					dir = getDir(fpdir, curdirpool);
					if (dir) { // fprintf(fpdecstr, "%s\n", rcstr);
						readsstr[rid] = strdup(rcstr);
					} else { // fprintf(fpdecstr, "%s\n", ori);
						readsstr[rid] = strdup(ori);
					}
					flag[rid] = true;
					prid = rid;
				}
			} 

			// while(fscanf(fpenstr, "%s", enstr) != EOF) {
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					// fprintf(stderr, "xx enstr: %s\n", enstr);
					dn = getDup(enstr);
					// cout << "000\n";
					//
					if (dn == 0) enstr[strlen(enstr) - 1] = '\0';
					split2TwoEncodedStr(enstr, ori, aft);
					// cout << "0001111\n";
					// cout << enstr << endl;
					// cout << ori << endl;
					// cout << aft << endl;
					// if (strlen(aft) > 0) exit(0);
					// cout << "---" << endl;

					decodeStr(reads[idx].str, ori, dstr);

					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}

					fporder.read((char*)&rid, sizeof(uint32_t));
					rid += reads[idx].prechildrid;
					readsstr[rid] = strdup(dstr);
					flag[rid] = true;

					reads[idx].prechildrid = rid;

					if (dn == 0) {
					} else 
					if (dn > 0) {
						// dir = getDir(fpdir, curdirpool);
						// if (dir) {
						// 	reverseComplement(dstr);
						// }
						
						// fprintf(fpdecstr, "%s\n", dstr);
						// cout << "00011112222\n";
						//
						strcpy(rcstr, dstr);
						reverseComplement(rcstr);
						// cout << "00011112222333\n";
						// cout << "111\n";

						// fporder.read((char*)&rid, sizeof(uint32_t));
						// readsstr[rid] = strdup(dstr);
						// }
						// flag[rid] = true;
						prid = rid;

						// if (rid == 0) cout << enstr << endl;

						for (int w = 0; w < dn; ++w) {
							fporder.read((char*)&rid, sizeof(uint32_t));
							rid += prid;

							dir = getDir(fpdir, curdirpool);
							if (dir) { // fprintf(fpdecstr, "%s\n", rcstr);
								readsstr[rid] = strdup(rcstr);
							} else { // fprintf(fpdecstr, "%s\n", dstr);
								readsstr[rid] = strdup(dstr);
							}
							flag[rid] = true;
							prid = rid;
						}
					}
					// cout << "222\n";
					// --num;
					// fprintf(stderr, "dstr: %s\n", dstr);
					// if (num <= 0) exit(0);

					READS_t r(dstr, idx);
					if (strlen(aft) > 0) { //
						decodeStr(dstr, aft, ctg);
						r.set(ctg);
					}
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
		
		// fproot.close();
		fpdir.close();

		/*char *seq = (char*)alloca((L + 3) * sizeof(char));
		FILE *fpsg = fopen(singletonfn.c_str(), "r");
		for (size_t id = 0; id < max_rid; ++id) {
			if (!flag[id]) {
				fscanf(fpsg, "%s", seq);
				readsstr[id] = strdup(seq);
			}
		}
		fclose(fpsg);*/

		uint32_t prerid = 0;
		while (fscanf(fproot, "%s", str) != EOF) {
			// fprintf(fpdecstr, "%s\n", str);
			fporder.read((char*)&rid, sizeof(uint32_t));
			rid += prerid;
			readsstr[rid] = strdup(str);
			flag[rid] = true;
			prerid = rid;
		}
		fclose(fproot);

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
			if (str[0] == '-') break;
			dn = getDup(str);

			split2TwoEncodedStr(str, ori, aft);
			vector<READS_t> reads;
			READS_t r(ori, 0);
			if (strlen(aft) > 0) {
				decodeStr(ori, aft, ctg);
				r.set(ctg);
			}
			reads.push_back(r);
			int idx = 0;
			fprintf(fpdecstr, "%s\n", ori);

			if (dn > 0) {
				strcpy(rcstr, ori);
				reverseComplement(rcstr);

				for (int w = 0; w < dn; ++w) {
					dir = getDir(fpdir, curdirpool);

					if (dir) fprintf(fpdecstr, "%s\n", rcstr);
					else fprintf(fpdecstr, "%s\n", ori);
				}
			}

			// while(fscanf(fpenstr, "%s", enstr) != EOF) {
			while (fgets(enstr, 1024, fpenstr) != NULL) {
				if (enstr[0] == '-') break;

				if (checkStrOrDig(enstr)) { // is encode string
					// fprintf(stderr, "xx enstr: %s\n", enstr);
					dn = getDup(enstr);
					//
					if (dn == 0) enstr[strlen(enstr) - 1] = '\0';
					split2TwoEncodedStr(enstr, ori, aft);
					// cout << enstr << endl;
					// cout << ori << endl;
					// cout << aft << endl;
					// if (strlen(aft) > 0) exit(0);
					// cout << "---" << endl;

					decodeStr(reads[idx].str, ori, dstr);
					//
					dir = getDir(fpdir, curdirpool);
					if (dir) {
						reverseComplement(dstr);
					}
					fprintf(fpdecstr, "%s\n", dstr);

					if (dn > 0) {
						strcpy(rcstr, dstr);
						reverseComplement(rcstr);

						for (int w = 0; w < dn; ++w) {
							dir = getDir(fpdir, curdirpool);

							if (dir) fprintf(fpdecstr, "%s\n", rcstr);
							else fprintf(fpdecstr, "%s\n", dstr);
						}
					}

					READS_t r(dstr, idx);

					if (strlen(aft) > 0) { //
						decodeStr(dstr, aft, ctg);
						// cout << r.str << endl;
						r.set(ctg);
						// cout << r.str << endl;
						// exit(0);
					}
					// --num;
					// fprintf(stderr, "dstr: %s\n", dstr);
					// if (num <= 0) exit(0);

					// READS_t r(ctg, idx);
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

		while (fscanf(fproot, "%s", str) != EOF) {
			fprintf(fpdecstr, "%s\n", str);
		}
		fclose(fproot);
	}

	fclose(fpdecstr);
}

void decompressSingle(string result) {
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	// string rootstrfn = folder + "rootstr.txt";
	// mstcom::bsc::BSC_decompress((rootstrfn+".bsc").c_str(), rootstrfn.c_str());
	// string orderfn;

	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	string singletonfn = folder + "singleton.bin"; 
	mstcom::bsc::BSC_decompress((singletonfn+".bsc").c_str(), singletonfn.c_str());

	// exit(0);
	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	BITPOOL curdirpool;
	curdirpool.n = 8;
	int dir, changed;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");

	FILE *fpdecstr = fopen(result.c_str(), "w");
	char str[1<<10], rcstr[1<<10], enstr[1<<10], dstr[1<<10], ori[1<<10], aft[1<<10], ctg[1<<10];

	// int num = 112;
	int isori, dn, idx;
	vector<READS_t> reads;

	while (fgets(enstr, 1024, fpenstr) != NULL) {
		if (enstr[0] == '-') break;

		if (checkStrOrDig(enstr)) { // is encode string
			dn = getDup(enstr);
			if (dn == 0) enstr[strlen(enstr) - 1] = '\0';

			isori = checkIsOriReads(enstr, L);

			if (isori == 0) {
				// 
				split2TwoEncodedStr(enstr, ori, aft);
				decodeStr(reads[idx].str, ori, dstr);
				//
				dir = getDir(fpdir, curdirpool);
				if (dir) {
					reverseComplement(dstr);
				}
				fprintf(fpdecstr, "%s\n", dstr);

				if (dn > 0) {
					strcpy(rcstr, dstr);
					reverseComplement(rcstr);

					for (int w = 0; w < dn; ++w) {
						dir = getDir(fpdir, curdirpool);

						if (dir) fprintf(fpdecstr, "%s\n", rcstr);
						else fprintf(fpdecstr, "%s\n", dstr);
					}
				}

				READS_t r(dstr, idx);

				if (strlen(aft) > 0) { //
					decodeStr(dstr, aft, ctg);
					r.set(ctg);
				}
				reads.push_back(r);
				idx = reads.size() - 1;
			} else { // root node
				reads.clear();

				split2TwoEncodedStr(enstr, ori, aft);
				READS_t r(ori, 0);
				if (strlen(aft) > 0) {
					decodeStr(ori, aft, ctg);
					r.set(ctg);
				}
				reads.push_back(r);
				idx = 0;

				fprintf(fpdecstr, "%s\n", ori);

				if (dn > 0) {
					strcpy(rcstr, ori);
					reverseComplement(rcstr);

					for (int w = 0; w < dn; ++w) {
						dir = getDir(fpdir, curdirpool);

						if (dir) fprintf(fpdecstr, "%s\n", rcstr);
						else fprintf(fpdecstr, "%s\n", ori);
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

	fpdir.close();
	while (fscanf(fpenstr, "%s", str) != EOF) {
		fprintf(fpdecstr, "%s\n", str);
	}
	fclose(fpenstr);

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

	fclose(fpdecstr);
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

void decompressSingleOrder(string result) {
	cout << "decompressSingleOrder(string result)...\n";
	string encodestrfn = folder + "encodestr.txt";
	mstcom::bsc::BSC_decompress((encodestrfn+".bsc").c_str(), encodestrfn.c_str());

	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string orderfn = folder + string("parent.bin");
	decompressFile(orderfn);
	// mstcom::lzma::lzma_decompress((orderfn+".lzma").c_str(), orderfn.c_str());
	std::ifstream fporder(orderfn, std::ios::binary);


	BITPOOL curdirpool;
	curdirpool.n = 8;
	int dir, changed;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");

	char enstr[1<<10], dstr[1<<10], ori[1<<10], aft[1<<10], ctg[1<<10];

	cout << "max_rid: " << max_rid << endl;
	uint32_t rid = 0, prid, trid;
	char **readsstr = new char*[max_rid];
	char **readsen = new char*[max_rid];
	bool *isrc = new bool[max_rid];
	uint32_t *parent = new uint32_t[max_rid];
	int isori;

	// FILE *fpout = fopen("decorder.txt", "w");
	while (fgets(enstr, 1024, fpenstr) != NULL) {
		enstr[strlen(enstr) - 1] = '\0';
		readsen[rid] = strdup(enstr);
		parent[rid] = rid;
		isori = checkIsOriReads(enstr, L);
		// if (rid == 384827) {
		// 	cout << enstr << endl;
		// 	cout << "isori: " << isori << endl;
		// }
		// exit(0);
		if (isori == 0) {
			fporder.read((char*)&trid, sizeof(uint32_t));
			// cout << "trid: " << trid << endl;
			// exit(0);
			parent[rid] = trid;
			isrc[rid] = getDir(fpdir, curdirpool);
		} else 
		if (isori == 1) {
			readsstr[rid] = strdup(readsen[rid]);
		} else { //isori == 2
			split2TwoEncodedStr(readsen[rid], ori, aft);
			if(strlen(aft) > 0) {
				decodeStr(ori, aft, ctg);
				if (strlen(readsen[rid]) > 0) free(readsen[rid]);
				readsen[rid] = strdup(ctg);
			}
			readsstr[rid] = strdup(ori);
		}
		// if (rid == 384827) exit(0);
		// if (parent[rid] != rid)
		// 	fprintf(fpout, "%s %u\n", enstr, parent[rid]);
		// else
		// 	fprintf(fpout, "%s\n", enstr);
		++ rid;

		// if (rid > 5) exit(0);
	}
	// fclose(fpout);
	fclose(fpenstr);
	fpdir.close();
	fporder.close();
	// cout << readsstr[174540] << endl;
	// cout << "xxxx\n";
// #ifdef false
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
				// if (debug) {
				// 	cout << "trid: " << trid << endl;
				// }
				split2TwoEncodedStr(readsen[trid], ori, aft);
				if (debug) {
					cout << "trid: " << trid << endl;
					cout << readsen[trid] << endl;
					cout << "parent[trid]: " << parent[trid] << endl;
					cout << "ref: " << readsen[parent[trid]] << endl;
					cout << "aft: " << aft << endl;
				}
				if (strlen(readsen[trid]) == 0) { // original duplicate reads
					decodeStr(readsstr[parent[trid]], ori, dstr);
				} else {
					if (readsen[trid][0] == '-') {
						ori[0] = '\0';
					} 
					decodeStr(readsen[parent[trid]], ori, dstr);
				}

				if (isrc[trid]) {
					reverseComplement(dstr);
				}

				if (debug) {
					cout << "ori: " << ori << endl;
					cout << "dstr: " << dstr << endl;
					cout << "------" << endl;
				}
				readsstr[trid] = strdup(dstr);
				parent[trid] = trid;

				if (strlen(readsen[trid]) > 0) free(readsen[trid]);
				
				// cout << "xxxx\n";
				if (strlen(aft) > 0) {
				// cout << "xxxx1111\n";
					decodeStr(dstr, aft, ctg);
				// cout << "xxxx11112222\n";
					readsen[trid] = strdup(ctg);
					if (debug) cout << "ctg: " << ctg << endl;
				} else {
					// cout << "yyy1111\n";
					readsen[trid] = strdup(dstr);
					// cout << "yyy1111222\n";
				}
				// cout << "n: " << n << endl;
				if (n == 0) break;
			}
			kv_destroy(v);
		}
	}
	// cout << "xxxxxxxyyyyyy\n";
// #endif
	FILE *fpdecstr = fopen(result.c_str(), "w");

	// cout << readsstr[0] << endl;
	// cout << "yyy" << endl;
	for (uint32_t id = 0; id < max_rid; ++id) {
		// cout << readsstr[id] << endl;
		fprintf(fpdecstr, "%s\n", readsstr[id]);
		// free(readsstr[id]);
		// free(readsen[id]);
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
		mstcom::lzma::lzma_decompress((orderfn+".lzma").c_str(), orderfn.c_str());
		fporder.open(orderfn, std::ios::binary);
	}
	string dirbinfn = folder + "dir.bin";
	mstcom::bsc::BSC_decompress((dirbinfn+".bsc").c_str(), dirbinfn.c_str());
	std::ifstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string singletonfn = folder + "singleton.bin"; 
	mstcom::bsc::BSC_decompress((singletonfn+".bsc").c_str(), singletonfn.c_str());
	
	string distfn = folder + "dist.bin";
	mstcom::lzma::lzma_decompress((distfn+".lzma").c_str(), distfn.c_str());

	string smdistfn = folder + "smdist.bin";
	mstcom::lzma::lzma_decompress((smdistfn+".lzma").c_str(), smdistfn.c_str());

	BITPOOL curdirpool;
	curdirpool.n = 8;
	int dir, file;

	FILE *fpenstr = fopen(encodestrfn.c_str(), "r");
	FILE *fproot = fopen(rootstrfn.c_str(), "r");
	// sprintf(name, "rootstr.bin");
	// std::ifstream fproot(name, std::ios::binary);

	FILE *fpdecstr1 = fopen((result+"_1").c_str(), "w");
	FILE *fpdecstr2 = fopen((result+"_2").c_str(), "w");
	char str[1<<10], rcstr[1<<10], enstr[1<<10], dstr[1<<10], ori[1<<10], aft[1<<10], ctg[1<<10];

	// int num = 112;
	int dn;

	size_t leftnum = 0;
	uint32_t rid = 0;
	char **readsstr = new char*[max_rid + 1];

	while (fscanf(fproot, "%s", str) != EOF) {
		if (str[0] == '-') break;
		dn = getDup(str);

		split2TwoEncodedStr(str, ori, aft);
		vector<READS_t> reads;
		READS_t r(ori, 0);
		if (strlen(aft) > 0) {
			decodeStr(ori, aft, ctg);
			r.set(ctg);
		}
		reads.push_back(r);
		int idx = 0;
		readsstr[rid ++] = strdup(ori);

		if (dn > 0) {
			strcpy(rcstr, ori);
			reverseComplement(rcstr);

			for (int w = 0; w < dn; ++w) {
				dir = getDir(fpdir, curdirpool);

				if (dir) readsstr[rid ++] = strdup(rcstr);
				else readsstr[rid ++] = strdup(ori);
			}
		}
		
		while (fgets(enstr, 1024, fpenstr) != NULL) {
			// cout << enstr << endl;
			if (enstr[0] == '-') break;

			if (checkStrOrDig(enstr)) { // is encode string
				dn = getDup(enstr);

				if (dn == 0) enstr[strlen(enstr) - 1] = '\0';
				split2TwoEncodedStr(enstr, ori, aft);

				decodeStr(reads[idx].str, ori, dstr);

				dir = getDir(fpdir, curdirpool);
				if (dir) {
					reverseComplement(dstr);
				}
				readsstr[rid ++] = strdup(dstr);
				// cout << dstr << endl;

				if (dn > 0) {
					strcpy(rcstr, dstr);
					reverseComplement(rcstr);

					for (int w = 0; w < dn; ++w) {
						dir = getDir(fpdir, curdirpool);

						if (dir) readsstr[rid ++] = strdup(rcstr);
						else readsstr[rid ++] = strdup(dstr);
					}
				}
				
				READS_t r(dstr, idx);
				if (strlen(aft) > 0) { //
					decodeStr(dstr, aft, ctg);
					r.set(ctg);
				}
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
		// if (aaaa > 5) 
		// exit(0);
	}

	fclose(fpenstr);
	
	while (fscanf(fproot, "%s", str) != EOF) {
		readsstr[rid ++] = strdup(str);
	}
	fclose(fproot);

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

	// FILE *fp = fopen("seq.txt", "w");
	// for (uint32_t tid = 0; tid < rid; ++tid) {
	// 	fprintf(fp, "%s\n", readsstr[tid]);
	// }
	// fclose(fp);

	// cout << "xxxx\n";
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

	if (isorder) {
		for (size_t num = 0; num < max_rid; ++num) {
			if (!flag[num]) {
				fporder.read((char*)&tempid, sizeof(uint32_t));

				dir = getDir(fpdir, curdirpool);
				file = getDir(fpdir, curdirpool);
				if (file) {
					fpdist.read((char*)&dist, sizeof(uint32_t));
				} else {
					smfpdist.read((char*)&disttemp, sizeof(uint16_t));
					dist = (uint32_t)disttemp;
				}
				// fpdist.read((char*)&dist, sizeof(int32_t));

				ids[num] = tempid;
				// if (dist < 0) {
					// dist = 0 - dist;
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
	} else {
		// cout << "xxxxx\n";
		size_t id = 0;
		for (size_t num = 0; num < max_rid; ++num) {
			if (!flag[num]) {
				dir = getDir(fpdir, curdirpool);
				// fpdist.read((char*)&dist, sizeof(int32_t));
				// fpdist.read((char*)&dist, sizeof(uint32_t));
				// cout << dist << endl;
				file = getDir(fpdir, curdirpool);
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
				// ids[num] = id ++;
				// if (dist < 0) {
				// 	dist = 0 - dist;
				// 	ids[num] += half;
				// }

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
	}
	// cout << "yyyy\n";

	fpdir.close();
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

	folder = generateString("mstcom", 10); //creat a temp folder for current input
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
		fclose(fp);
		// return 0;

		if (ispe) {
			// decompressPE(argv[2]);
			decompressPEX(outfile.c_str());
		} else {
			if (isorder) {
				decompressSingleOrder(outfile.c_str());
			} else 
				decompressSingle(outfile.c_str());
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
