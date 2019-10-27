// * History
// 改进计算minimizer的并行方法：
// 1. 先算minimizer并存入一个数组；使用原子运算统计bucket每行的个数；
// * ... *//
# include <cstdio>
# include <iostream>
# include <fstream>
# include <vector>
# include <queue>
# include <assert.h>
# include <chrono>
# include <thread>
# include <mutex>
# include <cstdlib>
# include <ctime>
# include <algorithm>
# include <map>
# include <unordered_map>
# include <unistd.h>
# include "mstcom.h"
// # include "util.h"
# include "libbsc/bsc.h"

using namespace std;

int nthreads = 24;
int max_dif_thr;
// const int K = 14;
int L;
size_t RN;
uint32_t max_rid;
bseq1_t *seq;
 
// input parameters
bool isorder;
bool ispe;
// string outputfn;

// --- for minimizer
int bsize = 20;
mm128_v *B;
uint64_v *Bs;

mm192_v *BL;

int startp[2], endp[2], sublen[2]; // the start and end positions of indexed substing 

// const int max_kmer = 45;
// const int min_kmer = 31;

/*const int max_kmer = 35;
const int min_kmer = 20;*/

// const int max_kmer = 39;
// const int min_kmer = 25;

// const int max_kmer = 29;
// const int min_kmer = 15;
int kmer, max_kmer, min_kmer;
int *kmervec, kmervecsize;

mutex *bmtx;

std::string folder;
void countTree();

struct READS_t {
	uint32_t prid; //reads id of parent
	uint32_t root, nextid; //root node of current tree 
	uint32_v crid;
	int dif;
	int16_t shift;
	bool isrc;
	READS_t(): nextid(0), shift(0), isrc(false) {}
	READS_t(int16_t _shift, bool _isrc): nextid(0), shift(_shift), isrc(_isrc) {}

	uint32_t getChildren() {
		if (nextid < crid.n) {
			nextid++;
			return crid.a[nextid - 1];
		}
		return max_rid + 1;
	}
};

READS_t *reads;
bool *isnextrnd;
uint32_t *prid2;
int16_t *shift2;
bool *isrc2;
int *min_dif2;
uint32_t *co;
bool *iscombined;

char complement[256];

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

inline void removeChild(size_t prid, size_t rid);
inline void cutCircle();
inline void connectWithPrid2();

static const char alphanum[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
inline std::string generateString(const std::string &chr, int length = 5) {
	srand(time(0));
	string res = chr + "_";
	for (int i = 0; i < length; ++i) {
		res += alphanum[rand() % 62];
	}
	return res;
}

char *cmd;

inline void init() {
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

string infile, infile1, outfile;

inline void displayHelp(const char* prog) {
	printf("mstcom v0.1, by Yuansheng Liu, April 2019.\n");
	printf("Usage: %s -i <.fastq> -o <output> [option parameters]\n", prog);
	printf("\t options:\n \t\t -f <.fastq> -t <threads>\n\n");
	// printf("-----------\n");
	printf("\t\t -i is fastq file\n");
	printf("\t\t -f is fastq file; only for paired-end reads\n");
	printf("\t\t -o is the output file\n");
	printf("\t\t -t is the number of threads\n");
	printf("\t\t -h print help message\n");

	printf("Example:\n\t\t");
	// printf("./bfmem -r H.all.fa -q M.all.fa -o hm-100.txt\n\n");
}

inline void getPars(int argc, char* argv[]) {
	isorder = false;
	ispe = false;
	bool isinfile = false, isoutfile = false; 
	int oc;
	while ((oc = getopt(argc, argv, "i:f:o:t:ph")) >= 0) {
		switch (oc) {
			case 'i':
				infile = optarg;
				isinfile = true;
				break;
			case 'f':
				infile1 = optarg;
				ispe = true;
				// cout << "f: ispe\n";
				break;
			case 'o':
				outfile = optarg;
				isoutfile = true;
				break;
			case 'p':
				isorder = true;
				break;
			case 't':
				nthreads = atoi(optarg);
				break;
			case 'h':
				displayHelp(argv[0]);
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

	if (ispe) {
	// fprintf(stderr, "ispe: %d\n", ispe);
		// fprintf(stderr, "%s\n", );
		f.open(infile1);
		if (f.fail()) {
			fprintf(stderr, "Input file '%s' does not exist.\n", infile1.c_str());
			exit(1);
		}
		f.close();
	}

	f.open(outfile);
	if (f.is_open()) {
		fprintf(stderr, "The output file '%s' exist.\n", outfile.c_str());
		exit(1);
	}
	f.close();

	folder = generateString("mstcom", 10); //creat a temp folder for current input

	sprintf(cmd, "mkdir -p %s", folder.c_str());
	fprintf(stderr, "%s\n", folder.c_str());
	system(cmd);
	folder += "/";
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

inline bool getReads(char const *fn) {
	bseq_file_t *fp;
	fp = bseq_open(fn);
	if (fp == 0) return false;
	seq = bseq_read(fp, &RN, L);
	bseq_close(fp);
	return true;
}

inline bool getRightReads(char const *fn) {
	bseq_file_t *fp;
	fp = bseq_open(fn);
	if (fp == 0) return false;
	fprintf(stdout, "Loading the second file %s ...\n", fn);
	size_t leftmaxrid = RN;
	bseq_read_second(seq, fp, &RN, L);
	bseq_close(fp); 
	if (leftmaxrid != (RN>>1)) {
		fprintf(stderr, "Two files contains different number of reads.\n");
		return false;
	}
	return true;
}

inline int getReadsLength(char const *fn) {
	bseq_file_t *fp;
	fp = bseq_open(fn);
	if (fp == 0) return false;
	int readslen = bseq_read_len(fp);
	bseq_close(fp);
	return readslen;
}

uint32_t rid_pthread;
const int mask = 1<<25 - 1;
mutex *rmtx;

const int intnum[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

bool debug = false;

inline int16_t diffstrlen(char *parent, char *child, const int16_t &b) {
	int eq_char_num = 0;
	int16_t dif = 0;
	// if (debug) {
	// 	cout << "b: " << b << endl;
	// 	cout << parent << endl;
	// 	cout << child << endl;
	// }
	// AAATGC --> parent
	//   ATGCAA --> child
	if (b >= 0) {
		dif = b;
		for (int i = b, j = 0; i < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				if (eq_char_num > 1) {
					dif += intnum[eq_char_num];
					// dif += 2;
					eq_char_num = 0;
				} else {
	                dif += eq_char_num;
	                eq_char_num = 0;
				}
				++ dif;
			} else ++ eq_char_num;
		}
	} else {
		dif = -b;
		for (int i = -b, j = 0; i < L; ++i, ++j) {
			if (child[i] != parent[j]) {
				if (eq_char_num > 1) {
					dif += intnum[eq_char_num];
					// dif += 2;
					eq_char_num = 0;
				} else {
	                dif += eq_char_num;
	                eq_char_num = 0;
				}
				++ dif;
			} else ++ eq_char_num;
		}
	}
	return dif;
}

inline int16_t diffstrlen(char *parent, char *_child, const int16_t &b, const bool &isrc) {
	char *child = (char*)alloca((L + 1) * sizeof(char));
	strcpy(child, _child);
	if (isrc) {
		reverseComplement(child);
	}
	int eq_char_num = 0;
	int16_t dif = 0;
	// if (debug) {
	// 	cout << "b: " << b << endl;
	// 	cout << parent << endl;
	// 	cout << child << endl;
	// }
	if (b >= 0) {
		dif = b;
		for (int i = b, j = 0; i < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				if (eq_char_num > 1) {
					dif += intnum[eq_char_num];
					// dif += 2;
					eq_char_num = 0;
				} else {
	                dif += eq_char_num;
	                eq_char_num = 0;
				}
				++ dif;
			} else ++ eq_char_num;
		}
	} else {
		dif = -b;
		for (int i = -b, j = 0; i < L; ++i, ++j) {
			if (child[i] != parent[j]) {
				if (eq_char_num > 1) {
					dif += intnum[eq_char_num];
					// dif += 2;
					eq_char_num = 0;
				} else {
	                dif += eq_char_num;
	                eq_char_num = 0;
				}
				++ dif;
			} else ++ eq_char_num;
		}
	}
	return dif;
}

inline void encode(char *parent, char *child, const int16_t &_shift, char *en_str) { //str1 is the short one string
	// char *en_str = (char*)alloca(((L<<2) + 1) * sizeof(char));
	// en_str = (char*)calloc(L<<1, sizeof(char));
	char *int_str = (char*)alloca(20 * sizeof(char));
	int16_t shift = _shift;

	int en_str_len = 0;
	int eq_char_num = 0;

	if (shift >= 0) {
		int i, j;
		for (i = shift, j = 0; i < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				if (eq_char_num > 1) {
					sprintf(int_str, "%d", eq_char_num);
					for (char *tk = int_str; *tk != '\0'; ++tk) {
						en_str[en_str_len++] = *tk;
					}
					eq_char_num = 0;
				} else {
					for (int i0 = j - eq_char_num; i0 < j; ++i0) {
	                    en_str[en_str_len++] = child[i0];
	                }
	                eq_char_num = 0;
				}
				en_str[en_str_len++] = child[j]; 
			} else ++ eq_char_num;
		}
		if (shift > 0)
		en_str[en_str_len++] = ' ';
		while (j < L) {
			en_str[en_str_len++] = child[j++];
		}
		en_str[en_str_len] = '\0';
	} else {
		int i, j;
		shift = 0 - shift;
		for (j = 0; j < shift; ++j) {
			en_str[en_str_len++] = child[j];
		}
		en_str[en_str_len++] = ' ';
		int spidx = en_str_len - 1;
		for (i = 0, j = shift; j < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				if (eq_char_num > 1) {
					sprintf(int_str, "%d", eq_char_num);
					for (char *tk = int_str; *tk != '\0'; ++tk) {
						en_str[en_str_len++] = *tk;
					}
					eq_char_num = 0;
				} else {
					for (int i0 = j - eq_char_num; i0 < j; ++i0) {
	                    en_str[en_str_len++] = child[i0];
	                }
	                eq_char_num = 0;
				}
				en_str[en_str_len++] = child[j]; 
			} else ++ eq_char_num;	
		}
		en_str[en_str_len] = '\0';

		// ???
		if (spidx > 0) {
			bool isdig = false;
			for (int i = 0; i < en_str_len; ++i) {
				if (en_str[i] >= '0' && en_str[i] <= '9') {
					isdig = true;
					break;
				}
			}
			if (!isdig) {
				en_str[en_str_len++] = '0';
				en_str[en_str_len] = '\0';
			}
		}
 	}
	// en_str[en_str_len] = '\0';
	// cout << en_str << endl;
	// cout << "dif: " << dif << "\n";
	// return dif;
	// return string(en_str);
	// return en_str_len;
}

inline void labelRoot() {
	size_t num_of_nodes = 0;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid) { //a root node 
			uint32_t noderid = rid;
			queue<uint32_t> q;
			q.push(noderid);
			while (!q.empty()) {
				noderid = q.front();
				reads[noderid].root = rid;
				q.pop();
				++ num_of_nodes;
				for (uint32_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}
		}
	}
	cout << "num_of_nodes: " << num_of_nodes << endl;
	// cout << "Time of labelRoot() = " << stopwatch.stop() << std::endl;
	// stopwatch.resume();
}

void calcMinimizersFun() { // before calling this method, MUST set rid_pthread = 0
	mm128_t minimizer;
	int mask = (1<<bsize) - 1, bucketidx, ncnt;
	while (1) {
		uint32_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;
		// cout << "rid: " << rid << endl;
		if (kmer == max_kmer) {
			isnextrnd[rid] = true;
			reads[rid].prid = rid; // itself
			kv_init(reads[rid].crid);
			// reads[rid].prid2 = rid;
			reads[rid].root = rid;

			ncnt = 0;
			for (int i = 0; seq[rid].seq[i]; ++i) {
				if (seq[rid].seq[i] == 'N') {
					++ ncnt;
				}
			}
			if (ncnt > L/2) {
				isnextrnd[rid] = false;
			}
		}

		if (isnextrnd[rid]) {
			mm_sketch(seq[rid].seq, L, kmer, rid, &minimizer);
			bucketidx = minimizer.x & mask;
			mm128_v *p = &B[bucketidx];
			bmtx[bucketidx].lock(); // the lock can be remove by first counting the number of this bucketidx
			kv_push(mm128_t, *p, minimizer);
			bmtx[bucketidx].unlock();
		}
	}
}

void calcLargeMinimizersFun() { // before calling this method, MUST set rid_pthread = 0
	mm192_t minimizer;
	int mask = (1<<bsize) - 1, bucketidx, ncnt;
	while (1) {
		uint32_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;
		// cout << "rid: " << rid << endl;
		if (kmer == max_kmer) {
			isnextrnd[rid] = true;
			reads[rid].prid = rid; // itself
			kv_init(reads[rid].crid);
			// reads[rid].prid2 = rid;
			reads[rid].root = rid;

			ncnt = 0;
			for (int i = 0; seq[rid].seq[i]; ++i) {
				if (seq[rid].seq[i] == 'N') {
					++ncnt;
				}
			}
			if (ncnt > L/2) {
				// reads[rid].isnextrnd = false;
				isnextrnd[rid] = false;
			}
		}

		if (isnextrnd[rid]) {
			mm_large_sketch(seq[rid].seq, L, kmer, rid, &minimizer);
			bucketidx = ((minimizer.x.x & mask) + (minimizer.x.y & mask)) & mask;
			mm192_v *p = &BL[bucketidx];
			bmtx[bucketidx].lock(); // the lock can be remove by first counting the number of this bucketidx
			kv_push(mm192_t, *p, minimizer);
			bmtx[bucketidx].unlock();
		}
	}
}

inline void calcMinimizers() {
	bmtx = new mutex[1 << bsize];
	rid_pthread = 0;
	std::vector<thread> threadVec;
	if (kmer <= 31) {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(calcMinimizersFun));
		} 
	} else {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(calcLargeMinimizersFun));
		}
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
	delete[] bmtx;
}

int cmpgroup(const void *a_, const void *b_) {
	uint64_t a = *(uint64_t*)a_, b = *(uint64_t*)b_;
	uint32_t rid_a = a>>32, rid_b = b>>32;
	int pos_a = (uint32_t)a>>1, pos_b = (uint32_t)b>>1;
	if (pos_a == pos_b) {
		return rid_a - rid_b;
	} 
	return pos_b - pos_a;// big to small
}

typedef std::unordered_map<std::string, std::vector<uint64_t> > StrMap;

inline void addStrMap(StrMap *smap, string bases, uint64_t y) {
	StrMap::iterator iter = smap->find(bases);
	if (iter != smap->end()) {
		iter->second.push_back(y);
	} else {
		vector<uint64_t> newvec;
		newvec.push_back(y);
		pair<StrMap::iterator, bool> res = smap->insert(make_pair(bases, newvec));
		if (!res.second) {
			res.first->second.push_back(y);
		}
	}
}

inline void linkReads2Tree(uint64_v *p) {
	uint64_t y = p->a[0];
	uint32_t prid = y >> 32, rid, trid;
	int ppos = (uint32_t)y >> 1, pos, tpos;
	int pdir = y & 1, dir, tdir;
	int16_t mindiff, difflen;
	int max_numberoffind;

	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));

	for (int k = 1; k < p->n; ++k) {
		max_numberoffind = 300;
		// if (k == 9) debug = true;
		y = p->a[k];
		rid = y >> 32;

		if (reads[rid].prid == rid || (reads[rid].prid != rid && abs(reads[rid].shift) > 1 && reads[rid].dif > 5)) {

			if (reads[rid].prid == rid) {
				mindiff = L;
			} else {
				mindiff = reads[rid].dif;
			}
			
			pos = (uint32_t)y>>1;
			dir = y & 1;
			strcpy(stra, seq[rid].seq);
			if (dir) {
				reverseComplement(stra);
			}
			int prenum = 0;
			if (k - 1 > 10000) max_numberoffind = 1000;

			for (int q = k - 1; q >= 0; --q) {
				tpos = (uint32_t) p->a[q] >> 1;
				tdir = p->a[q] & 1;
				// if (k - q >= 100 && tpos - pos >= 5) {
				if (k - q >= 100 && tpos - pos >= 15) { // tpos > pos ?
					break;
				} else { //tpos == pos; maybe equal
					trid = p->a[q] >> 32;

					// if (reads[trid].prid != rid) { /// ???? not sure
						strcpy(strb, seq[trid].seq);
						if (tdir) {
							reverseComplement(strb);
						}

						if(strcmp(stra, strb) == 0) {
							prid = trid;
							ppos = tpos;
							pdir = tdir;
							mindiff = 0;
							break;
						} else {
							difflen = diffstrlen(strb, stra, tpos - pos);
							if (difflen < mindiff) {
								prid = trid;
								ppos = tpos;
								pdir = tdir;

								mindiff = difflen;
							}
						}
					// }
				}
				if (debug) return;

				++prenum;
				// if (prenum >= 100) break;
				if (prenum >= max_numberoffind) break;
				// if (prenum >= 200) break;
				// if ((prenum >= 100 && abs(tpos - pos) >= 5) || prenum >= 10000) break;
				// if (prenum >= 10000) break;
				// if (prenum >= 300) break;
			}

			// fprintf(stderr, "mindiff: %d\n", mindiff);
			if (mindiff <= max_dif_thr && (reads[rid].prid == rid || (reads[rid].prid != rid && mindiff < reads[rid].dif))) {
				if (reads[rid].prid != rid) {
					rmtx[reads[rid].prid & mask].lock();
					removeChild(reads[rid].prid, rid);
					rmtx[reads[rid].prid & mask].unlock();
				}

				reads[rid].prid = prid;
				reads[rid].isrc = pdir ^ dir;
				reads[rid].shift = ppos - pos;
				reads[rid].dif = mindiff;

				if (pdir) {
					reads[rid].shift = pos - ppos;
				}
				
				// rmtx[prid & mask].lock();
				kv_push(uint32_t, reads[prid].crid, rid);
				// rmtx[prid & mask].unlock();
				if (mindiff == 0) {
					// reads[rid].isnextrnd = false;
					isnextrnd[rid] = false;
				}
			}
		}
	}
}

inline void linkReads2TreeWithIndex(uint64_v *p) {
	uint64_t y, ty;
	uint32_t prid, rid, trid;
	int ppos, pos, tpos;
	int pdir, dir, tdir;
	int16_t mindiff, difflen;

	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));

	StrMap *smap = new StrMap[2]; // smap[0], smap[1];
	string idxstr, curstr;
	int cstartp;

	y = p->a[0];
	strcpy(stra, seq[y >> 32].seq);
	if (y & 1) reverseComplement(stra);

	curstr = (string)(stra);
	idxstr = curstr.substr(startp[0], sublen[0]);
	addStrMap(&smap[0], idxstr, y);
	idxstr = curstr.substr(startp[1], sublen[1]);
	addStrMap(&smap[1], idxstr, y);

	for (int k = 1; k < p->n; ++k) {
		y = p->a[k];
		rid = y >> 32;
		pos = (uint32_t)y>>1;
		dir = y & 1;
		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra);
		}
		curstr = (string)(stra);

		if (reads[rid].prid == rid || (reads[rid].prid != rid && abs(reads[rid].shift) > 1 && reads[rid].dif > 5)) {

			if (reads[rid].prid == rid) {
				mindiff = L;
			} else {
				mindiff = reads[rid].dif;
			}
			
			int numsearched = 0;
			for (int q = k - 1; q >= 0; --q) {
				tpos = (uint32_t) p->a[q] >> 1;
				tdir = p->a[q] & 1;
				// if (k - q >= 100 && tpos - pos >= 5) {
				if (k - q >= 100 && tpos - pos >= 15) { // tpos > pos ?
					break;
				} else { //tpos == pos; maybe equal
					trid = p->a[q] >> 32;

					strcpy(strb, seq[trid].seq);
					if (tdir) {
						reverseComplement(strb);
					}

					if(strcmp(stra, strb) == 0) {
						prid = trid;
						ppos = tpos;
						pdir = tdir;
						mindiff = 0;
						break;
					} else {
						difflen = diffstrlen(strb, stra, tpos - pos);
						if (difflen < mindiff) {
							prid = trid;
							ppos = tpos;
							pdir = tdir;

							mindiff = difflen;
						}
					}
				}
				if (debug) return;

				++ numsearched;
				if (numsearched >= 100) break;
			}

			if (mindiff != 0 && k > 100) {
				for (int shift = 0; shift <= 15 && mindiff != 0; ++shift) {
					if (mindiff < shift) break;

					for (int mpi = 0; mpi < 2 && mindiff != 0; ++mpi) {
						if (mindiff < shift) break;

						cstartp = startp[mpi] - shift;
						if (cstartp >= 0) {
							idxstr = curstr.substr(cstartp, sublen[mpi]);

							StrMap::iterator iter = smap[mpi].find(idxstr);
							if (iter != smap[mpi].end()) {
								for (size_t j = 0; j < iter->second.size(); ++j) {
									ty = iter->second[j];
									//
									tpos = (uint32_t) ty >> 1;
									tdir = ty & 1;

									if (abs(tpos - pos) <= mindiff && shift == abs(tpos - pos)) {
										trid = ty >> 32;

										strcpy(strb, seq[trid].seq);
										if (tdir) {
											reverseComplement(strb);
										}

										if (strcmp(stra, strb) == 0) {
											prid = trid;
											ppos = tpos;
											pdir = tdir;
											mindiff = 0;
											break;
										} else {
											difflen = diffstrlen(strb, stra, tpos - pos);

											if (difflen < mindiff) {
												prid = trid;
												ppos = tpos;
												pdir = tdir;

												mindiff = difflen;

												if (mindiff < shift) break;
											}
										}
									}
								}
							}
						}
					}
				}
			}

			// fprintf(stderr, "mindiff: %d\n", mindiff);
			if (mindiff <= max_dif_thr && (reads[rid].prid == rid || (reads[rid].prid != rid && mindiff < reads[rid].dif))) {
				if (reads[rid].prid != rid) {
					rmtx[reads[rid].prid & mask].lock();
					removeChild(reads[rid].prid, rid);
					rmtx[reads[rid].prid & mask].unlock();
				}

				reads[rid].prid = prid;
				reads[rid].isrc = pdir ^ dir;
				reads[rid].shift = ppos - pos;
				reads[rid].dif = mindiff;

				if (pdir) {
					reads[rid].shift = pos - ppos;
				}
				
				// rmtx[prid & mask].lock();
				kv_push(uint32_t, reads[prid].crid, rid);
				// rmtx[prid & mask].unlock();
				if (mindiff == 0) {
					// reads[rid].isnextrnd = false;
					isnextrnd[rid] = false;
				}
			}
		}

		idxstr = curstr.substr(startp[0], sublen[0]);
		addStrMap(&smap[0], idxstr, y);
		idxstr = curstr.substr(startp[1], sublen[1]);
		addStrMap(&smap[1], idxstr, y);

	}

	delete[] smap;
}

// FILE *fptmp;
// mutex fptmpmtx;
mmrec_v rec, lrec;
mutex recmtx;

void sortBucketsFun() {
	uint32_t j, start_a, n;
	mmrec_t ttmmrec;

	while (1) {
		uint32_t bid = __sync_fetch_and_add(&rid_pthread, 1);
		if (bid >= (1ul << bsize)) break;

		mm128_v *b = &B[bid];
		// cout << "bb: " << bb << endl;
		if (b->n > 0) {
			// cout << "b->n: " << b->n << endl;
			radix_sort_128x(b->a, b->a + b->n);

			kv_resize(uint64_t, Bs[bid], b->n);
			for (j = 0; j < b->n; ++j) {
				Bs[bid].a[j] = 0;
			}
			Bs[bid].n = b->n;

			for (j = 1, n = 1, start_a = 0; j <= b->n; ++j) {
				if (j == b->n || b->a[j].x != b->a[j-1].x) {
					// fprintf(stderr, "j: %lu;\nstart_a: %lu\nn: %lu", j, start_a, n);
					assert(j - start_a == n);

					if (n >= 2) { // tree
						// record bid, start_a and n;
						ttmmrec.x = bid, ttmmrec.y = start_a, ttmmrec.z = n;
						recmtx.lock();
						kv_push(mmrec_t, rec, ttmmrec);
						recmtx.unlock();
					} 
					start_a = j, n = 1;
				} else {++n;}
			}
		}
	}
}

bool mm192cmp(const mm192_t &a, const mm192_t &b) {
	return a.x.y < b.x.y || (a.x.y == b.x.y && a.x.x < b.x.x);
}

bool ivcmpx(const uint64_t &a, const uint64_t &b) {
	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));

	// cout << (a>>32) << endl;
	strcmp(stra, seq[a >> 32].seq);
	// cout << stra << endl;
	if (a & 1) reverseComplement(stra);

	// cout << (b>>32) << endl;
	strcmp(strb, seq[b >> 32].seq);
	// cout << strb << endl;
	if (b & 1) reverseComplement(strb);

	return strcmp(stra, strb) >= 0;
}

int ivcmp(const void *a_, const void *b_) {
	uint64_t a = *(uint64_t*)a_, b = *(uint64_t*)b_;
	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));

	strcmp(stra, seq[(uint32_t)(a >> 32)].seq);
	if (a & 1) reverseComplement(stra);

	strcmp(strb, seq[(uint32_t)(b >> 32)].seq);
	if (b & 1) reverseComplement(strb);
	return strcmp(stra, strb);
}

void sortLargeBucketsFun() {
	uint32_t j, start_a, n;
	mmrec_t ttmmrec;

	while (1) {
		uint32_t bid = __sync_fetch_and_add(&rid_pthread, 1);
		if (bid >= (1ul << bsize)) break;

		mm192_v *b = &BL[bid];
		// cout << "bb: " << bb << endl;
		if (b->n > 0) {
			// cout << "b->n: " << b->n << endl;
			// radix_sort_128x(b->a, b->a + b->n);
			sort(b->a, b->a + b->n, mm192cmp);

			kv_resize(uint64_t, Bs[bid], b->n);
			for (j = 0; j < b->n; ++j) {
				Bs[bid].a[j] = 0;
			}
			Bs[bid].n = b->n;

			for (j = 1, n = 1, start_a = 0; j <= b->n; ++j) {
				if (j == b->n || b->a[j].x != b->a[j-1].x) {
					// fprintf(stderr, "j: %lu;\nstart_a: %lu\nn: %lu", j, start_a, n);
					assert(j - start_a == n);

					if (n >= 2) { // tree
						// record bid, start_a and n;
						ttmmrec.x = bid, ttmmrec.y = start_a, ttmmrec.z = n;
						recmtx.lock();
						kv_push(mmrec_t, rec, ttmmrec);
						recmtx.unlock();
					} 
					start_a = j, n = 1;
				} else {++n;}
			}
		}
	}
}

bool reccmp(const mmrec_t &a, const mmrec_t &b) {
	return a.z > b.z;
}

inline void sortBuckets() {
	kv_init(rec);
	kv_resize(mmrec_t, rec, 1<<16);

	rid_pthread = 0;
	std::vector<thread> threadVec;
	if (kmer <= 31) {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(sortBucketsFun));
		}
	} else {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(sortLargeBucketsFun));
		}
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();

	// sort rec
	sort (rec.a, rec.a + rec.n, reccmp);
	
	fprintf(stderr, "rec.n: %lu\n", rec.n);
	for (int i = 0; i < 5; ++i) {
		fprintf(stderr, "%lu\n", rec.a[i].z);
	}
}

void processBucketsFun(bool isfind) {
	mm128_t minimizer;
	int mask = (1<<bsize) - 1, bucketidx;
	uint32_t rid;
	uint32_t bid, j, start_a, n;
	uint64_t aa;

	char *temp_str = (char*)alloca((L + 1) * sizeof(char));
	// cout << "text xxx\n";
	uint64_v *iv = new uint64_v[L];
	uint64_v p;

	while (1) {
		uint32_t i = __sync_fetch_and_add(&rid_pthread, 1);
		if (i >= rec.n) break;

		bid = rec.a[i].x;
		start_a = rec.a[i].y;
		n = rec.a[i].z;

		// fprintf(stderr, "%lu %lu %lu\n", rec.a[i].x, rec.a[i].y, rec.a[i].z);
		// if (i > 10) exit(0);

		for (int k = 0; k < L; ++k) {
			kv_init(iv[k]);
		}

		if (kmer <= 31) {
			mm128_v *b = &B[bid];
			for (uint32_t k = 0; k < n; ++k) {
				aa = b->a[start_a + k].y;
				kv_push(uint64_t, iv[(uint32_t)aa >> 1], aa);
			}
		} else {
			mm192_v *b = &BL[bid];
			for (uint32_t k = 0; k < n; ++k) {
				aa = b->a[start_a + k].y;
				kv_push(uint64_t, iv[(uint32_t)aa >> 1], aa);
			}
		}
		// cout << "xxx\n";
		// cout << "n: " << n << endl;

		kv_init(p);
		kv_resize(uint64_t, p, n);
		// cout << "xxx12112121221\n";
		for (uint32_t k = L - 1; k > 0; --k) {
			// cout << "k: " << k << "; iv[k].n: " << iv[k].n << endl;
			if (iv[k].n > 0) {
				// sort iv[k].a, iv[k].a + iv[k].n
				// cout << "---\n" << endl;
				// for (int ss = 0; ss < iv[k].n; ++ss) {
				// 	cout << seq[(iv[k].a[ss] >> 32)].seq << endl;
				// }
				// cout << "---\n" << endl;

				// if (iv[k].n > 300) //sort (iv[k].a, iv[k].a + iv[k].n, ivcmp);
					// qsort(iv[k].a, iv[k].n, sizeof(uint64_t), ivcmp);

				// cout << "---\n" << endl;
				// for (int ss = 0; ss < iv[k].n; ++ss) {
				// 	cout << seq[(iv[k].a[ss] >> 32)].seq << endl;
				// }
				// cout << "---\n" << endl;
				// exit (0);

				for (int i1 = 0; i1 < iv[k].n; ++i1) {
					kv_push(uint64_t, p, iv[k].a[i1]);
				}
			}
		}
		// cout << "xxx111\n";
		for (int k = 0; k < L; ++k) {
			kv_destroy(iv[k]);
		}
		// cout << "xxx222\n";
		for (uint32_t k = 0; k < n; ++k) {
			Bs[bid].a[start_a + k] = p.a[k];
		}
		// cout << "xxx333\n";
		/*uint64_v *p = new uint64_v;
		kv_init(*p);
		kv_resize(uint64_t, *p, n);
		for (uint32_t k = 0; k < n; ++k) {
			kv_push(uint64_t, *p, b->a[start_a + k].y);
		}
		qsort(p->a, n, sizeof(uint64_t), cmpgroup);
		for (uint32_t k = 0; k < n; ++k) {
			Bs[bid].a[start_a + k] = p->a[k];
		}*/

		/*if (n > 5) {
			for (uint32_t k = 0; k < 5; ++k) {
				cout << ((uint32_t)p.a[k]>>1) << " " << seq[p.a[k]>>32].seq << endl;
			}
			exit (0);
		}*/

		// if (n > 4060000) {
		/*if (kmer < max_kmer) {
			cout << "n: " << n << endl;
			int pos;
			uint64_t y;
			char *str = (char*)alloca((L + 1) * sizeof(char));
			FILE *fp = fopen("seq.txt", "w");
			for (int i = 0; i < n; ++i) {
				y = p.a[i];
				rid = y >> 32;
				strcpy(str, seq[rid].seq);
				if (y&1) reverseComplement(str);
				pos = (uint32_t)y>>1;

				fprintf(fp, "%d %s\n", pos, str);
			}
			fclose(fp);
			exit (0);
		}*/

		if (isfind) {
			// if (n <= 100) {
			// 	linkReads2Tree(&p); // divide into two situation
			// } else {
			// }

			// if (n <= 5000) {
			if (n <= 10000) {
				linkReads2TreeWithIndex(&p);
			} else {
				recmtx.lock();
				kv_push(mmrec_t, lrec, rec.a[i]);
				recmtx.unlock();
			}
		}
		kv_destroy(p);
		
	}
	delete[] iv;
}

inline void processBuckets(bool isfind) {
	// fptmp = fopen("fptmp.txt", "w");

	// if (kmer == max_kmer) nthreads = 1; // for test
	kv_init(lrec);
	kv_resize(mmrec_t, lrec, 1<<8);

	rid_pthread = 0;
	std::vector<thread> threadVec;
	for (int i = 0; i < nthreads; ++i) {
		threadVec.push_back(std::thread(processBucketsFun, isfind));
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
	// fclose(fptmp);
	// kv_destroy(rec);

	
	// for test
	// for test
	// for test
	// for test
	// uint64_t y;
	// uint32_t rid;
	// uint32_t bid, j, start_a, n;
	// bid = rec.a[0].x;
	// start_a = rec.a[0].y;
	// n = rec.a[0].z;

	// uint64_v *b = &Bs[bid];
	// uint64_v *p = (uint64_v*) calloc(1, sizeof(uint64_v));

	// kv_init(*p);
	// kv_resize(uint64_t, *p, n);
	// for (uint32_t k = 0; k < n; ++k) {
	// 	kv_push(uint64_t, *p, b->a[start_a + k]);
	// }
	// char *stra = (char*)alloca((L + 1) * sizeof(char));
	// FILE *fp = fopen("tempseq.txt", "w");
	// int pos, dir;

	// for (uint32_t k = 0; k < n; ++k) {
	// 	y = p->a[k];
	// 	rid = y >> 32;

	// 	pos = (uint32_t)y>>1;
	// 	dir = y & 1;
	// 	strcpy(stra, seq[rid].seq);
	// 	if (dir) {
	// 		reverseComplement(stra);
	// 	}
	// 	fprintf(fp, "%d %s\n", pos, stra);

	// }
	// fclose(fp);
	
	// kv_destroy(*p);
	// free(p);
	// exit(0);
}

void processLargeBucketsFun(StrMap *smap, uint64_v *p) {
	uint64_t y, ty;
	uint32_t prid, rid, trid;
	int ppos, pos, tpos;
	int pdir, dir, tdir;
	int16_t shift;
	int16_t mindiff, oridiff, difflen;

	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));

	string idxstr, curstr;
	int cstartp;

	while (1) {
		uint32_t k = __sync_fetch_and_add(&rid_pthread, 1);
		if (k >= p->n) break;

		// k
		y = p->a[k];
		rid = y >> 32;
		pos = (uint32_t)y>>1;
		dir = y & 1;
		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra);
		}
		curstr = (string)(stra);

		if (reads[rid].prid == rid || (reads[rid].prid != rid && abs(reads[rid].shift) > 1 && reads[rid].dif > 5)) {
			if (reads[rid].prid == rid) {
				mindiff = L;
			} else {
				mindiff = reads[rid].dif;
			}
			
			int numsearched = 0;
			for (int q = k - 1; q >= 0; --q) {
				tpos = (uint32_t) p->a[q] >> 1;
				tdir = p->a[q] & 1;

				if (k - q >= 100 && tpos - pos >= 15) { // tpos > pos ?
					break;
				} else { //tpos == pos; maybe equal
					trid = p->a[q] >> 32;

					strcpy(strb, seq[trid].seq);
					if (tdir) {
						reverseComplement(strb);
					}

					if(strcmp(stra, strb) == 0) {
						prid = trid;
						ppos = tpos;
						pdir = tdir;
						mindiff = 0;
						break;
					} else {
						difflen = diffstrlen(strb, stra, tpos - pos);
						if (difflen < mindiff) {
							prid = trid;
							ppos = tpos;
							pdir = tdir;

							mindiff = difflen;
						}
					}
				}
				if (debug) return;

				++ numsearched;
				if (numsearched >= 100) break;
			}

			if (mindiff != 0 && k > 100) {
				for (int shift = 0; shift <= 15 && mindiff != 0; ++shift) {
					if (mindiff < shift) break;

					for (int mpi = 0; mpi < 2 && mindiff != 0; ++mpi) {
						if (mindiff < shift) break;

						cstartp = startp[mpi] - shift;
						if (cstartp >= 0) {
							idxstr = curstr.substr(cstartp, sublen[mpi]);

							StrMap::iterator iter = smap[mpi].find(idxstr);
							if (iter != smap[mpi].end()) {
								for (size_t j = 0; j < iter->second.size(); ++j) {
									ty = iter->second[j];
									//
									tpos = (uint32_t) ty >> 1;
									tdir = ty & 1;
									trid = ty >> 32;

									// cout << "111\n";
									if (co[trid] < k && abs(tpos - pos) <= mindiff && shift == abs(tpos - pos)) {

										strcpy(strb, seq[trid].seq);
										if (tdir) {
											reverseComplement(strb);
										}

										if (strcmp(stra, strb) == 0) {
											prid = trid;
											ppos = tpos;
											pdir = tdir;
											mindiff = 0;
											break;
										} else {
											difflen = diffstrlen(strb, stra, tpos - pos);

											if (difflen < mindiff) {
												prid = trid;
												ppos = tpos;
												pdir = tdir;

												mindiff = difflen;

												if (mindiff < shift) break;
											}
										}
									}
									// cout << "222\n";
								}
							}
						}
					}
				}
			}

			

			if (mindiff <= max_dif_thr && (reads[rid].prid == rid || (reads[rid].prid != rid && mindiff < reads[rid].dif))) {
				if (reads[rid].prid != rid) {
					rmtx[reads[rid].prid & mask].lock();
					removeChild(reads[rid].prid, rid);
					rmtx[reads[rid].prid & mask].unlock();
				}

				reads[rid].prid = prid;
				reads[rid].isrc = pdir ^ dir;
				reads[rid].shift = ppos - pos;
				reads[rid].dif = mindiff;

				if (pdir) {
					reads[rid].shift = pos - ppos;
				}
				
				rmtx[prid & mask].lock();
				kv_push(uint32_t, reads[prid].crid, rid);
				rmtx[prid & mask].unlock();

				if (mindiff == 0) {
					isnextrnd[rid] = false;
				}
			}
		}
	}
}

void processLargeBucketsRandomFun(StrMap *smap, uint64_v *p) {
	uint64_t y, ty;
	uint32_t prid, rid, trid;
	int ppos, pos, tpos;
	int pdir, dir, tdir;
	int16_t shift;
	int16_t mindiff, oridiff, difflen;

	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));

	string idxstr, curstr;
	int cstartp;

	while (1) {
		uint32_t k = __sync_fetch_and_add(&rid_pthread, 1);
		if (k >= p->n) break;
		// k
		y = p->a[k];
		rid = y >> 32;
		pos = (uint32_t)y>>1;
		dir = y & 1;
		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra);
		}
		curstr = (string)(stra);

		if (reads[rid].prid == rid || (reads[rid].prid != rid && abs(reads[rid].shift) > 1 && reads[rid].dif > 5)) {
			if (reads[rid].prid == rid) {
				mindiff = L;
			} else {
				mindiff = reads[rid].dif;
			}
			
			int numsearched = 0;
			for (int q = k - 1; q >= 0; --q) {
				tpos = (uint32_t) p->a[q] >> 1;
				tdir = p->a[q] & 1;

				if (k - q >= 100 && tpos - pos >= 15) { // tpos > pos ?
					break;
				} else { //tpos == pos; maybe equal
					trid = p->a[q] >> 32;

					strcpy(strb, seq[trid].seq);
					if (tdir) {
						reverseComplement(strb);
					}

					if(strcmp(stra, strb) == 0) {
						prid = trid;
						ppos = tpos;
						pdir = tdir;
						mindiff = 0;
						break;
					} else {
						difflen = diffstrlen(strb, stra, tpos - pos);
						if (difflen < mindiff) {
							prid = trid;
							ppos = tpos;
							pdir = tdir;

							mindiff = difflen;
						}
					}
				}
				if (debug) return;

				++ numsearched;
				if (numsearched >= 100) break;
			}

			if (mindiff != 0 && k > 100) {
				for (int shift = 0; shift <= 15 && mindiff != 0; ++shift) {
					if (mindiff < shift) break;

					for (int mpi = 0; mpi < 2 && mindiff != 0; ++mpi) {
						if (mindiff < shift) break;

						cstartp = startp[mpi] - shift;
						if (cstartp >= 0) {
							idxstr = curstr.substr(cstartp, sublen[mpi]);

							StrMap::iterator iter = smap[mpi].find(idxstr);
							if (iter != smap[mpi].end()) {
								for (size_t j = 0; j < iter->second.size(); ++j) {
									ty = iter->second[j];
									//
									tpos = (uint32_t) ty >> 1;
									tdir = ty & 1;
									trid = ty >> 32;

									// cout << "111\n";
									if (co[trid] != k && abs(tpos - pos) <= mindiff && shift == abs(tpos - pos)) {

										strcpy(strb, seq[trid].seq);
										if (tdir) {
											reverseComplement(strb);
										}

										if (strcmp(stra, strb) == 0) {
											prid = trid;
											ppos = tpos;
											pdir = tdir;
											mindiff = 0;
											break;
										} else {
											difflen = diffstrlen(strb, stra, tpos - pos);

											if (difflen < mindiff) {
												prid = trid;
												ppos = tpos;
												pdir = tdir;

												mindiff = difflen;

												if (mindiff < shift) break;
											}
										}
									}
									// cout << "222\n";
								}
							}
						}
					}
				}
			}

			

			if (mindiff <= max_dif_thr && (reads[rid].prid == rid || (reads[rid].prid != rid && mindiff < reads[rid].dif))) {
				if (reads[rid].prid != rid) {
					rmtx[reads[rid].prid & mask].lock();
					removeChild(reads[rid].prid, rid);
					rmtx[reads[rid].prid & mask].unlock();
				}

				reads[rid].prid = prid;
				reads[rid].isrc = pdir ^ dir;
				reads[rid].shift = ppos - pos;
				reads[rid].dif = mindiff;

				if (pdir) {
					reads[rid].shift = pos - ppos;
				}
				
				rmtx[prid & mask].lock();
				kv_push(uint32_t, reads[prid].crid, rid);
				rmtx[prid & mask].unlock();

				if (mindiff == 0) {
					isnextrnd[rid] = false;
				}
			}
		}
	}
}

inline void processLargeBuckets() {
	if (lrec.n == 0) return;
	stopwatch.resume();

	cout << "lrec.n: " << lrec.n << endl;

	char *stra = (char*)alloca((L + 1) * sizeof(char));

	std::vector<thread> threadVec;
	uint64_v p;
	uint64_t y;
	string idxstr, curstr;

	uint32_t n, start_a;
	for (size_t i = 0; i < lrec.n; ++i) {
		// cout << "processLargeBuckets() i: " << i << endl;

		uint64_v *b = &Bs[lrec.a[i].x];
		start_a = lrec.a[i].y;
		n = lrec.a[i].z;

		StrMap *smap = new StrMap[2];

		kv_init(p);
		kv_resize(uint64_t, p, n);
		for (uint32_t k = 0; k < n; ++k) {
			// creat hash table here?
			y = b->a[start_a + k];
			strcpy(stra, seq[y >> 32].seq);
			if (y & 1) reverseComplement(stra);

			co[y >> 32] = k;

			curstr = (string)(stra);
			idxstr = curstr.substr(startp[0], sublen[0]);
			addStrMap(&smap[0], idxstr, y);
			idxstr = curstr.substr(startp[1], sublen[1]);
			addStrMap(&smap[1], idxstr, y);

			kv_push(uint64_t, p, y);
		}
		// cout << "before processLargeBucketsFun() ...\n";

		// nthreads = 1;

		if (n >= 500000) {
			rid_pthread = rand()% (n - 20010);
			p.n = rid_pthread + 20000;
			// if (p->n > 500000 && k >= 20000) break;
			for (int i = 0; i < nthreads; ++i) {
				threadVec.push_back(std::thread(processLargeBucketsRandomFun, smap, &p));
			}
		} else {
			rid_pthread = 0;
			for (int i = 0; i < nthreads; ++i) {
				threadVec.push_back(std::thread(processLargeBucketsFun, smap, &p));
			}
		}
		std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
			thr.join();
		});
		threadVec.clear();

		// cout << "after ...\n";

		kv_destroy(p);
		smap[0].clear();
		smap[1].clear();
		delete[] smap;
	}
	kv_destroy(lrec);

	cout << "Time of processLargeBuckets() = " << stopwatch.stop() << std::endl;
	stopwatch.resume();
}

// int minicount;

void calcMinimizers2Fun() { // before calling this method, MUST set rid_pthread = 0
	mm128_t minimizer;
	int mask = (1<<bsize) - 1, bucketidx, ncnt;
	while (1) {
		uint32_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;
		// cout << "rid: " << rid << endl;
		if (kmer == max_kmer) {
			prid2[rid] = rid; // itself  //???
		}

		// if (reads[rid].isnextrnd) {
		if (isnextrnd[rid]) {
			// __sync_fetch_and_add(&minicount, 1);
			mm_sketch(seq[rid].seq, L, kmer, rid, &minimizer);
			bucketidx = minimizer.x & mask;
			mm128_v *p = &B[bucketidx];
			bmtx[bucketidx].lock();
			kv_push(mm128_t, *p, minimizer);
			bmtx[bucketidx].unlock();
		}
	}
}

void calcLargeMinimizers2Fun() { // before calling this method, MUST set rid_pthread = 0
	mm192_t minimizer;
	int mask = (1<<bsize) - 1, bucketidx, ncnt;
	while (1) {
		uint32_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;
		// cout << "rid: " << rid << endl;
		if (kmer == max_kmer) {
			prid2[rid] = rid; // itself  //???
		}

		// if (reads[rid].isnextrnd) {
		if (isnextrnd[rid]) {
			// __sync_fetch_and_add(&minicount, 1);
			mm_large_sketch(seq[rid].seq, L, kmer, rid, &minimizer);
			bucketidx = ((minimizer.x.x & mask) + (minimizer.x.y & mask)) & mask;
			mm192_v *p = &BL[bucketidx];
			bmtx[bucketidx].lock(); // the lock can be remove by first counting the number of this bucketidx
			kv_push(mm192_t, *p, minimizer);
			bmtx[bucketidx].unlock();
		}
	}
}

inline void calcMinimizers2() {
	// minicount = 0;
	bmtx = new mutex[1 << bsize];
	rid_pthread = 0;
	std::vector<thread> threadVec;
	if (kmer <= 31) {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(calcMinimizers2Fun));
		}
	} else {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(calcLargeMinimizers2Fun));
		}
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
	delete[] bmtx;
	// cout << "minicount: " << minicount << endl;
}

inline void findPrid2(uint64_v *p) {
	
	uint64_t y = p->a[0];
	uint32_t prid = y >> 32, rid, trid;
	int ppos = (uint32_t)y >> 1, pos, tpos;
	int pdir = y & 1, dir, tdir;
	int16_t shift;
	int max_numberoffind;
	// if (pdir) {
	// 	ppos = L - ppos + kmer - 2;
	// }
	int16_t mindiff, oridiff, difflen;

	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));

	for (int k = 1; k < p->n; ++k) {

		max_numberoffind = 600;

		// if (k == 9) debug = true;
		y = p->a[k];
		rid = y >> 32;

		if (reads[prid2[rid]].root == reads[rid].root || (reads[prid2[rid]].root != reads[rid].root && abs(shift2[rid]) > 1 && min_dif2[rid] > 5)) {
			// if (reads[rid].prid2 == rid) {
			if (prid2[rid] == rid) {
				mindiff = max_dif_thr + 1;
			} else {
				mindiff = min_dif2[rid];
			}
			
			pos = (uint32_t)y>>1;
			dir = y & 1;
			strcpy(stra, seq[rid].seq);
			if (dir) {
				reverseComplement(stra);
			}
			// if (dir) {
			// 	pos = L - pos + kmer - 2;
			// }
			if (k - 1 > 10000) max_numberoffind = 2000;

			int prenum = 0;
			for (int q = k - 1; q >= 0; --q) {
				tpos = (uint32_t) p->a[q] >> 1;
				tdir = p->a[q] & 1;
				// if (tdir) {
				// 	tpos = L - tpos + kmer - 2;
				// }

				if (abs(tpos - pos) <= mindiff) { //tpos == pos; maybe equal
					trid = p->a[q] >> 32;

					if (reads[trid].root != reads[rid].root) {
						strcpy(strb, seq[trid].seq);
						if (tdir) {
							reverseComplement(strb);
						}
						difflen = diffstrlen(strb, stra, tpos - pos);

						if (difflen < mindiff) {
							prid = trid;
							ppos = tpos;
							pdir = tdir;

							mindiff = difflen;
						}

						++ prenum; //
					}
				} else {
					break;
				}
				if (debug) return;

				// if (tpos != pos) ++prenum;
				// if (prenum >= 300) break;
				if (prenum >= max_numberoffind) break;
				// if (prenum >= 600) break;
				// if (prenum >= 300 && abs(tpos - pos) >= 10) break;
				// if (prenum >= 500) break;
				// if (prenum >= 20000) break;
			}

			// fprintf(stderr, "mindiff: %d\n", mindiff);
			// if (mindiff <= L/3) { //????????????? to be confirmed
			if (mindiff <= max_dif_thr && (reads[prid2[rid]].root == reads[rid].root || (reads[prid2[rid]].root != reads[rid].root && mindiff < min_dif2[rid]))) {
				 // && reads[prid].prid != rid

				prid2[rid] = prid;
				isrc2[rid] = pdir ^ dir;
				shift2[rid] = ppos - pos;
				min_dif2[rid] = mindiff;

				if (pdir) {
					shift2[rid] = pos - ppos;
				}
			}
		}
	}
}

inline void findPrid2WithIndex(uint64_v *p) {
	uint64_t y, ty;
	uint32_t prid, rid, trid;
	int ppos, pos, tpos;
	int pdir, dir, tdir;
	int16_t shift;
	int16_t mindiff, oridiff, difflen;

	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));

	StrMap *smap = new StrMap[2]; // smap[0], smap[1];
	string idxstr, curstr;
	int cstartp;

	y = p->a[0];
	strcpy(stra, seq[y >> 32].seq);
	if (y & 1) reverseComplement(stra);

	curstr = (string)(stra);
	idxstr = curstr.substr(startp[0], sublen[0]);
	addStrMap(&smap[0], idxstr, y);
	idxstr = curstr.substr(startp[1], sublen[1]);
	addStrMap(&smap[1], idxstr, y);

	for (int k = 1; k < p->n; ++k) {
		y = p->a[k];
		
		rid = y >> 32;
		pos = (uint32_t)y>>1;
		dir = y & 1;
		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra);
		}
		curstr = (string)(stra);

		// if (reads[prid2[rid]].root == reads[rid].root || (reads[prid2[rid]].root != reads[rid].root && abs(shift2[rid]) > 1 && min_dif2[rid] > 5)) {
		if (reads[prid2[rid]].root == reads[rid].root) { // now in the same tree
			// if (reads[rid].prid2 == rid) {
			mindiff = max_dif_thr + 1;
			
			// if (prid2[rid] == rid) {
			// } else {
			// 	mindiff = min_dif2[rid];
			// }
			
			int numsearched = 0;
			for (int q = k - 1; q >= 0; --q) {
				tpos = (uint32_t) p->a[q] >> 1;
				tdir = p->a[q] & 1;

				if (abs(tpos - pos) <= mindiff) { //tpos == pos; maybe equal
					trid = p->a[q] >> 32;

					if (reads[trid].root != reads[rid].root) {
						strcpy(strb, seq[trid].seq);
						if (tdir) {
							reverseComplement(strb);
						}
						difflen = diffstrlen(strb, stra, tpos - pos);

						if (difflen < mindiff) {
							prid = trid;
							ppos = tpos;
							pdir = tdir;

							mindiff = difflen;
						}

						++ numsearched;
					}
				} else {
					break;
				}
				if (debug) return;
				if (numsearched >= 300) break;
			}

			if (k > 300) { // search in hash table

				for (int shift = 0; shift <= 50; ++shift) { // large than 15 ???
					for (int mpi = 0; mpi < 2; ++mpi) {
						cstartp = startp[mpi] - shift;
						if (cstartp >= 0) {
							idxstr = curstr.substr(cstartp, sublen[mpi]);

							StrMap::iterator iter = smap[mpi].find(idxstr);
							if (iter != smap[mpi].end()) {
								for (size_t j = 0; j < iter->second.size(); ++j) {
									ty = iter->second[j];
									//
									tpos = (uint32_t) ty >> 1;
									tdir = ty & 1;

									if (abs(tpos - pos) <= mindiff && shift == abs(tpos - pos)) {
										trid = ty >> 32;

										if (reads[trid].root != reads[rid].root) {
											strcpy(strb, seq[trid].seq);
											if (tdir) {
												reverseComplement(strb);
											}
											difflen = diffstrlen(strb, stra, tpos - pos);

											if (difflen < mindiff) {
												prid = trid;
												ppos = tpos;
												pdir = tdir;

												mindiff = difflen;
											}
										}

									}
								}
							}
						}
					}
				}

			}

			// if (mindiff <= max_dif_thr && (reads[prid2[rid]].root == reads[rid].root || (reads[prid2[rid]].root != reads[rid].root && mindiff < min_dif2[rid]))) {
			if (mindiff <= max_dif_thr) {
				 // && reads[prid].prid != rid
				prid2[rid] = prid;
				isrc2[rid] = pdir ^ dir;
				shift2[rid] = ppos - pos;
				min_dif2[rid] = mindiff;

				if (pdir) {
					shift2[rid] = pos - ppos;
				}
			}
		}


		idxstr = curstr.substr(startp[0], sublen[0]);
		addStrMap(&smap[0], idxstr, y);
		idxstr = curstr.substr(startp[1], sublen[1]);
		addStrMap(&smap[1], idxstr, y);

	}

	smap[0].clear();
	smap[1].clear();
	delete[] smap;
}

void processBuckets2Fun() {
	mm128_t minimizer;
	int mask = (1<<bsize) - 1, bucketidx;
	uint32_t rid;
	uint32_t bid, j, start_a, n;

	char *temp_str = (char*)alloca((L + 1) * sizeof(char));
	// cout << "text xxx\n";

	while (1) {
		uint32_t i = __sync_fetch_and_add(&rid_pthread, 1);
		if (i >= rec.n) break;

		// if (rec.a[i].z > 50000) continue; // the cluster to large, not consider
		// if (rec.a[i].z > 10000) {
		// if (rec.a[i].z > 5000) {
		if (rec.a[i].z > 10000) {
			recmtx.lock();
			kv_push(mmrec_t, lrec, rec.a[i]);
			recmtx.unlock();
			continue;
		}

		bid = rec.a[i].x;
		start_a = rec.a[i].y;
		n = rec.a[i].z;

		// fprintf(stderr, "%lu: %lu %lu %lu\n", i, rec.a[i].x, rec.a[i].y, rec.a[i].z);
		// if (i > 10) exit(0);
		// fprintf(stderr, "%lu\n", i);

		// mm128_v *b = &B[bid];
		uint64_v *b = &Bs[bid];
		uint64_v *p = (uint64_v*) calloc(1, sizeof(uint64_v));

		kv_init(*p);
		kv_resize(uint64_t, *p, n);
		for (uint32_t k = 0; k < n; ++k) {
			// kv_push(uint64_t, *p, b->a[start_a + k].y);
			kv_push(uint64_t, *p, b->a[start_a + k]);
		}
		// qsort(p->a, n, sizeof(uint64_t), cmpgroup);
		// fprintf(stderr, "n: %lu; p.n: %lu\n", n, p->n);
		// if (n <= 300) {
		// 	findPrid2(p);
		// } else {
		// }

		findPrid2WithIndex(p);
		
		kv_destroy(*p);
		free(p);
		// fprintf(stderr, "%lu\n", i);
	}
}

inline void processBuckets2() {
	CStopWatch tstopwatch;
	tstopwatch.start();
	tstopwatch.resume();

	kv_init(lrec);
	kv_resize(mmrec_t, lrec, 1<<11);

	rid_pthread = 0;
	std::vector<thread> threadVec;
	// fprintf(stderr, "rec.n: %lu\n", rec.n);
	for (int i = 0; i < nthreads; ++i) {
		threadVec.push_back(std::thread(processBuckets2Fun));
	}
	// cout << "111xxxx22222\n";
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	// cout << "111xxxx222224444\n";
	threadVec.clear();
	// fclose(fptmp);
	kv_destroy(rec);

	cout << "Time of processBuckets2() = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();
}

void processLargeBuckets2Fun(StrMap *smap, uint64_v *p) {
	uint64_t y, ty;
	uint32_t prid, rid, trid;
	int ppos, pos, tpos;
	int pdir, dir, tdir;
	int16_t shift;
	int16_t mindiff, oridiff, difflen;

	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));

	string idxstr, curstr;
	int cstartp;

	while (1) {
		uint32_t k = __sync_fetch_and_add(&rid_pthread, 1);
		if (k >= p->n) break;
		// k
		y = p->a[k];
		rid = y >> 32;
		pos = (uint32_t)y>>1;
		dir = y & 1;
		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra);
		}
		curstr = (string)(stra);

		if (reads[prid2[rid]].root == reads[rid].root || (reads[prid2[rid]].root != reads[rid].root && abs(shift2[rid]) > 1 && min_dif2[rid] > 5)) {
			if (prid2[rid] == rid) {
				mindiff = max_dif_thr + 1;
			} else {
				mindiff = min_dif2[rid];
			}
			
			int numsearched = 0;

			for (int q = k - 1; q >= 0; --q) {
				tpos = (uint32_t) p->a[q] >> 1;
				tdir = p->a[q] & 1;

				if (abs(tpos - pos) <= mindiff) { //tpos == pos; maybe equal
					trid = p->a[q] >> 32;

					if (reads[trid].root != reads[rid].root) {
						strcpy(strb, seq[trid].seq);
						if (tdir) {
							reverseComplement(strb);
						}
						difflen = diffstrlen(strb, stra, tpos - pos);

						if (difflen < mindiff) {
							prid = trid;
							ppos = tpos;
							pdir = tdir;

							mindiff = difflen;
						}

						++ numsearched;
					}
				} else {
					break;
				}
				if (debug) return;

				if (numsearched >= 300) break;
			}

			if (k > 300) {
				for (int shift = 0; shift <= 50; ++shift) {
					if (mindiff < shift) break;

					for (int mpi = 0; mpi < 2; ++mpi) {
						if (mindiff < shift) break;

						cstartp = startp[mpi] - shift;
						if (cstartp >= 0) {
							idxstr = curstr.substr(cstartp, sublen[mpi]);

							StrMap::iterator iter = smap[mpi].find(idxstr);
							if (iter != smap[mpi].end()) {
								for (size_t j = 0; j < iter->second.size(); ++j) {
									ty = iter->second[j];
									//
									tpos = (uint32_t) ty >> 1;
									tdir = ty & 1;
									trid = ty >> 32;

									if (co[trid] < k && abs(tpos - pos) <= mindiff && shift == abs(tpos - pos)) {
										if (reads[trid].root != reads[rid].root) {
											strcpy(strb, seq[trid].seq);
											if (tdir) {
												reverseComplement(strb);
											}
											difflen = diffstrlen(strb, stra, tpos - pos);

											if (difflen < mindiff) {
												prid = trid;
												ppos = tpos;
												pdir = tdir;

												mindiff = difflen;

												if (mindiff < shift) break;
											}
										}

									}
								}
							}
						}
					}
				}
			}

			if (mindiff <= max_dif_thr && (reads[prid2[rid]].root == reads[rid].root || (reads[prid2[rid]].root != reads[rid].root && mindiff < min_dif2[rid]))) {
				prid2[rid] = prid;
				isrc2[rid] = pdir ^ dir;
				shift2[rid] = ppos - pos;
				min_dif2[rid] = mindiff;

				if (pdir) {
					shift2[rid] = pos - ppos;
				}
			}
		}
	}
}

void processLargeBuckets2RandomFun(StrMap *smap, uint64_v *p) {
	uint64_t y, ty;
	uint32_t prid, rid, trid;
	int ppos, pos, tpos;
	int pdir, dir, tdir;
	int16_t shift;
	int16_t mindiff, oridiff, difflen;

	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));

	string idxstr, curstr;
	int cstartp;

	while (1) {
		uint32_t k = __sync_fetch_and_add(&rid_pthread, 1);
		if (k >= p->n) break;

		// k
		y = p->a[k];
		rid = y >> 32;
		pos = (uint32_t)y>>1;
		dir = y & 1;
		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra);
		}
		curstr = (string)(stra);

		if (reads[prid2[rid]].root == reads[rid].root || (reads[prid2[rid]].root != reads[rid].root && abs(shift2[rid]) > 1 && min_dif2[rid] > 5)) {
			if (prid2[rid] == rid) {
				mindiff = max_dif_thr + 1;
			} else {
				mindiff = min_dif2[rid];
			}
			
			int numsearched = 0;

			for (int q = k - 1; q >= 0; --q) {
				tpos = (uint32_t) p->a[q] >> 1;
				tdir = p->a[q] & 1;

				if (abs(tpos - pos) <= mindiff) { //tpos == pos; maybe equal
					trid = p->a[q] >> 32;

					if (reads[trid].root != reads[rid].root) {
						strcpy(strb, seq[trid].seq);
						if (tdir) {
							reverseComplement(strb);
						}
						difflen = diffstrlen(strb, stra, tpos - pos);

						if (difflen < mindiff) {
							prid = trid;
							ppos = tpos;
							pdir = tdir;

							mindiff = difflen;
						}

						++ numsearched;
					}
				} else {
					break;
				}
				if (debug) return;

				if (numsearched >= 300) break;
			}

			if (k > 300) {
				for (int shift = 0; shift <= 50; ++shift) {
					if (mindiff < shift) break;

					for (int mpi = 0; mpi < 2; ++mpi) {
						if (mindiff < shift) break;

						cstartp = startp[mpi] - shift;
						if (cstartp >= 0) {
							idxstr = curstr.substr(cstartp, sublen[mpi]);

							StrMap::iterator iter = smap[mpi].find(idxstr);
							if (iter != smap[mpi].end()) {
								for (size_t j = 0; j < iter->second.size(); ++j) {
									ty = iter->second[j];
									//
									tpos = (uint32_t) ty >> 1;
									tdir = ty & 1;
									trid = ty >> 32;
									// if (co[trid] < k && abs(tpos - pos) <= mindiff && shift == abs(tpos - pos)) {
									if (co[trid] != k && abs(tpos - pos) <= mindiff && shift == abs(tpos - pos)) {

										if (reads[trid].root != reads[rid].root) {
											strcpy(strb, seq[trid].seq);
											if (tdir) {
												reverseComplement(strb);
											}
											difflen = diffstrlen(strb, stra, tpos - pos);

											if (difflen < mindiff) {
												prid = trid;
												ppos = tpos;
												pdir = tdir;

												mindiff = difflen;

												if (mindiff < shift) break;
											}
										}

									}
								}
							}
						}
					}
				}
			}

			if (mindiff <= max_dif_thr && (reads[prid2[rid]].root == reads[rid].root || (reads[prid2[rid]].root != reads[rid].root && mindiff < min_dif2[rid]))) {
				prid2[rid] = prid;
				isrc2[rid] = pdir ^ dir;
				shift2[rid] = ppos - pos;
				min_dif2[rid] = mindiff;

				if (pdir) {
					shift2[rid] = pos - ppos;
				}
			}
		}
	}
}

inline void processLargeBuckets2() {
	if (lrec.n == 0) return;
	CStopWatch tstopwatch;
	tstopwatch.start();
	tstopwatch.resume();

	cout << "lrec.n: " << lrec.n << endl;

	char *stra = (char*)alloca((L + 1) * sizeof(char));

	std::vector<thread> threadVec;
	uint64_v p;
	uint64_t y;
	string idxstr, curstr;

	uint32_t n, start_a;
	for (size_t i = 0; i < lrec.n; ++i) {
		// cout << "processLargeBuckets() i: " << i << endl;

		uint64_v *b = &Bs[lrec.a[i].x];
		start_a = lrec.a[i].y;
		n = lrec.a[i].z;

		StrMap *smap = new StrMap[2];

		kv_init(p);
		kv_resize(uint64_t, p, n);
		for (uint32_t k = 0; k < n; ++k) {
			y = b->a[start_a + k];
			strcpy(stra, seq[y >> 32].seq);
			if (y & 1) reverseComplement(stra);

			co[y >> 32] = k;

			curstr = (string)(stra);
			idxstr = curstr.substr(startp[0], sublen[0]);
			addStrMap(&smap[0], idxstr, y);
			idxstr = curstr.substr(startp[1], sublen[1]);
			addStrMap(&smap[1], idxstr, y);

			kv_push(uint64_t, p, y);
		}

		if (n >= 500000) {
			// if (p->n > 500000 && k >= 20000) break;
			rid_pthread = rand()% (n - 20010);
			p.n = rid_pthread + 20010;
			for (int i = 0; i < nthreads; ++i) {
				threadVec.push_back(std::thread(processLargeBuckets2RandomFun, smap, &p));
			}
		} else {
			rid_pthread = 0;
			for (int i = 0; i < nthreads; ++i) {
				threadVec.push_back(std::thread(processLargeBuckets2Fun, smap, &p));
			}
		}
		std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
			thr.join();
		});
		threadVec.clear();

		kv_destroy(p);
		smap[0].clear();
		smap[1].clear();
		delete[] smap;
	}
	kv_destroy(lrec);

	cout << "Time of processLargeBuckets2() = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();
}

int tree_num;

void countTreeFun() {
	while (1) {
		uint32_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;

		if (reads[rid].prid == rid && reads[rid].crid.n > 0) {
			__sync_fetch_and_add(&tree_num, 1);
		}
	}
}

inline void countTree() {
	tree_num = 0;
	rid_pthread = 0;
	std::vector<thread> threadVec;
	for (int i = 0; i < nthreads; ++i) {
		threadVec.push_back(std::thread(countTreeFun));
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
	fprintf(stderr, "trees number: %d\n", tree_num);
}

string fnv[65];

void outputORI();

inline void treeConstuction() {
	CStopWatch tstopwatch;
	tstopwatch.start();
	stopwatch.start();

	if (min_kmer <= 31) {
		B = (mm128_v*)calloc(1 << bsize, sizeof(mm128_v));
	}
	if (max_kmer > 31) {
		BL = (mm192_v*)calloc(1 << bsize, sizeof(mm192_v));
	}
	Bs = (uint64_v*)calloc(1 << bsize, sizeof(uint64_v));

	co = new uint32_t[max_rid];

	for (int kid = 0; kid < kmervecsize; ++kid) {
		kmer = kmervec[kid];

		cout << "kmer: " << kmer << endl;
		tstopwatch.resume();

		if (kmer <= 31) {
			for (int i = 0; i < (1 << bsize); ++i) {
				kv_init(B[i]);
				kv_init(Bs[i]);
				// kv_resize(mm128_t, B[i], 1024);
			}
		} else {
			for (int i = 0; i < (1 << bsize); ++i) {
				kv_init(BL[i]);
				kv_init(Bs[i]);
				// kv_resize(mm128_t, B[i], 1024);
			}
		}
		// cout << "kmer: " << kmer << endl;
		calcMinimizers();

		cout << "Time of calcMinimizers() = " << tstopwatch.stop() << std::endl;
		tstopwatch.resume();
		// cout << "after calcMinimizers()\n";
		sortBuckets();

		cout << "Time of sortBuckets() = " << tstopwatch.stop() << std::endl;
		tstopwatch.resume();
		// cout << "after sortBuckets()\n";
		processBuckets(true);

		cout << "Time of processBuckets() = " << tstopwatch.stop() << std::endl;
		tstopwatch.resume();
		// cout << "after processBuckets()\n";

		// processLargeBuckets();

		if (kmer < max_kmer) {
			// dump bucket
			fnv[kmer] = folder + to_string(kmer);
			// cout << "kmer: " << kmer << " " << fnv[kmer] << endl;
			bucket_dump(fnv[kmer], Bs);
		}
		kv_destroy(rec);

		// cout << "after processBuckets()\n";
		if (kmer == max_kmer) {
			countTree();

			/*uint32_t cnt = 0;
			for (uint32_t rd = 0; rd < max_rid; ++rd) {
				if (isnextrnd[rd]) {
					++ cnt;
				}
			}
			cout << "cnt: " << cnt << endl;
			exit (0);*/
		}

		if (kmer <= 31) { 
			for (int i = 0; i < (1 << bsize); ++i) {
				kv_destroy(B[i]);
				kv_destroy(Bs[i]);
			}
		} else {
			for (int i = 0; i < (1 << bsize); ++i) {
				kv_destroy(BL[i]);
				kv_destroy(Bs[i]);
			}	
		}
		// cout << "Time of initialTreeConstuction kmer =" << kmer << " = " << stopwatch.stop() << std::endl;
		// stopwatch.resume();
		-- kmer;
	}
	tstopwatch.resume();

	cutCircle();
	labelRoot();

	cout << "Time of cutCircle() && labelRoot() = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();
	// countTree();
	cout << "Time of step 1: " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	// outputORI();

	// exit(0);
	// for kmer = 31;
	kmer = max_kmer;
	if (kmer <= 31) {
		for (int i = 0; i < (1 << bsize); ++i) {
			kv_init(B[i]);
			kv_init(Bs[i]);
		}
	} else {
		for (int i = 0; i < (1 << bsize); ++i) {
			kv_init(BL[i]);
			kv_init(Bs[i]);
		}
	}
	prid2 = new uint32_t[max_rid];

	calcMinimizers2();
	sortBuckets();
	processBuckets(false);
	fnv[kmer] = folder + to_string(kmer);
	// cout << "kmer: " << kmer << " " << fnv[kmer] << endl;
	bucket_dump(fnv[kmer], Bs);
	kv_destroy(rec);

	if (kmer <= 31) {
		for (int i = 0; i < (1 << bsize); ++i) {
			kv_destroy(B[i]);
			kv_destroy(Bs[i]);
		}
	} else {
		for (int i = 0; i < (1 << bsize); ++i) {
			kv_destroy(BL[i]);
			kv_destroy(Bs[i]);
		}
	}

	// exit(0);
	
	// size_t *prid2;
	// int16_t *shift2;
	// bool *isrc2;
	// int *min_dif2;
	shift2 = new int16_t[max_rid];
	isrc2 = new bool[max_rid];
	min_dif2 = new int[max_rid];

	iscombined = new bool[max_rid];

	int pre_trees_cnt = 0;
	while (true) {
		countTree();
		if (abs(pre_trees_cnt - tree_num) < 10000) break;
		// if (abs(pre_trees_cnt - tree_num) < 100000) break;
		// if (abs(pre_trees_cnt - tree_num) < 1000000) break;
		// if (tree_num < 1000000) break;
		// if (abs(pre_trees_cnt - tree_num) < 1000000) break;
		// if (abs(pre_trees_cnt - tree_num) < 10000 || tree_num <= 10000000) break;
		// if (abs(pre_trees_cnt - tree_num) < 10000) break;
		// if (tree_num <= 1000000) break;
		// if (abs(pre_trees_cnt - tree_num) < 100) break;
		pre_trees_cnt = tree_num;
		stopwatch.resume();
	
		// kmer = max_kmer;
		// fprintf(stderr, "kmer: %lu\n", kmer);
		// while (kmer >= min_kmer) {
		for (int kid = 0; kid < kmervecsize; ++kid) {
			kmer = kmervec[kid];
			cout << "kmer: " << kmer << endl;
			// stopwatch.resume();
			// B = (mm128_v*)calloc(1 << bsize, sizeof(mm128_v));
			for (int i = 0; i < (1 << bsize); ++i) {
				kv_init(Bs[i]);
				// kv_resize(mm128_t, B[i], 1024);
			}
			// cout << "before mm_bucket_load()\n";
			// calcMinimizers2();
			// cout << "after calcMinimizers2()\n";
			// return ;
			// sortBuckets();
			// cout << "after sortBuckets()\n";
			if (bucket_load(fnv[kmer], Bs) == false) {
				cout << "bucket_load() fail...\n";
				exit(1);
			}

			// if (kmer == max_kmer) {
				
			// }
			// cout << kmer << " before processBuckets2()\n";
			processBuckets2();
			// cout << kmer << " after processBuckets2()\n";

			// processLargeBuckets2();
			// cout << kmer << " after processLargeBuckets2()\n";

			// if (kmer == max_kmer) countTree();
			for (int i = 0; i < (1 << bsize); ++i) {
				kv_destroy(Bs[i]);
			}
			// cout << "Time of findPrid2 kmer = " << kmer << " = " << stopwatch.stop() << std::endl;
			// stopwatch.resume();
			-- kmer;
		}
		// cout << "111\n";
		tstopwatch.resume();

		connectWithPrid2();

		cout << "Time of connectWithPrid2() = " << tstopwatch.stop() << std::endl;
		tstopwatch.resume();
		// cout << "after connectWithPrid2()\n";
		// cutCircle();
		// cout << "after cutCircle()\n";
		// countTree();
		labelRoot();

		cout << "Time of findPrid2 = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
	}
	free(B);

	delete[] co;
	delete[] prid2;
	delete[] shift2;
	delete[] isrc2;
	delete[] min_dif2;

	delete[] iscombined;

	countTree();

	return;
	//for test
	bool *visited = (bool*)alloca(max_rid * sizeof(bool));
	// (char*)alloca(20 * sizeof(char));
	// memset(visited, false, sizeof(bool)*max_rid);
	for (size_t rid = 0; rid < max_rid; ++rid) {
		visited[rid] = false;
	}
	char *temp_str = (char*)alloca((L + 1) * sizeof(char));
	FILE *fout = fopen("out.txt", "w");
	char *enstr = (char*)alloca((L<<1) * sizeof(char));
	for (size_t rid = 0; rid < max_rid; ++rid) {
		if (!visited[rid] && reads[rid].prid == rid && reads[rid].crid.n > 0) { // a root node
			queue<size_t> q;
			q.push(rid);
			while (!q.empty()) {
				rid = q.front();
				q.pop();

				visited[rid] = true;
				if (reads[rid].prid != rid) {
					strcpy(temp_str, seq[rid].seq);
					if (reads[rid].isrc) {
						reverseComplement(temp_str);
					}
					encode(seq[reads[rid].prid].seq, temp_str, reads[rid].shift, enstr);
					
					fprintf(fout, "--encode--\n%s\n%s\n%s\n", seq[reads[rid].prid].seq, temp_str, enstr);
					// fprintf(fout, "%s\n", enstr);
				}

				// fprintf(stderr, "%lu: ", rid);
				for (size_t i = 0; i < reads[rid].crid.n; ++i) {
					// fprintf(stderr, "%lu, ", reads[rid].crid.a[i]);
					q.push(reads[rid].crid.a[i]);
				}
				// fprintf(stderr, "\n");
			}
		}
	}
	fclose(fout);
}

struct ROOTNODE_t {
	uint32_t rid, nodecnt;
	ROOTNODE_t() {}
	ROOTNODE_t(uint32_t _rid, uint32_t _nodecnt): rid(_rid), nodecnt(_nodecnt) {}
};

bool cmp(const ROOTNODE_t &a, const ROOTNODE_t &b) {
	return a.nodecnt < b.nodecnt;
}

bool cmp2(const ROOTNODE_t &a, const ROOTNODE_t &b) {
	return a.nodecnt > b.nodecnt;
}

int rootnodevecsize;

inline void reRootNode(const size_t &_rid, const size_t &root);
inline void reRootNode(const size_t &_rid);


inline void removeChild(size_t prid, size_t rid) {
	size_t cridsize = reads[prid].crid.n;
	// if (debug) fprintf(stderr, "cridsize: %lu\n", reads[prid].crid.size());
	for (size_t i = 0; i < cridsize; ++i) {
		if (reads[prid].crid.a[i] == rid) {
			reads[prid].crid.a[i] = reads[prid].crid.a[cridsize - 1];
			break;
		}
	} 
	// reads[prid].crid.resize(cridsize - 1);
	-- reads[prid].crid.n; //resize(cridsize - 1);
}

inline int overlap(char *parent, char *child, const int16_t &b) { //str1 is the short one string
	int eq_char_num = 0, dif = 0;
	
	if (b >= 0) {
		for (int i = b, j = 0; i < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				if (eq_char_num > 1) {
					dif += intnum[eq_char_num];
					// dif += 2;
					eq_char_num = 0;
				} else {
	                dif += eq_char_num;
	                eq_char_num = 0;
				}
				++ dif;
			} else ++ eq_char_num;
		}
	} else {
		for (int i = -b, j = 0; i < L; ++i, ++j) {
			if (child[i] != parent[j]) {
				if (eq_char_num > 1) {
					dif += intnum[eq_char_num];
					// dif += 2;
					eq_char_num = 0;
				} else {
	                dif += eq_char_num;
	                eq_char_num = 0;
				}
				++ dif;
			} else ++ eq_char_num;
		}
	}
	return dif;
}

inline void reRootNode(const uint32_t &_rid, const uint32_t &root) { // set rid as root node
	if (reads[_rid].prid == _rid) return; // rid is root node; do nothing
	// fprintf(stderr, "reRootNode()...\n");
	// fprintf(stderr, "%lu %lu\n", _rid, reads[_rid].prid);

	bool debug = false;
	// if (_rid == 49755) debug = true;
	if (debug) fprintf(stderr, "------\n");

	uint32_t rid = _rid;
	uint32_t oriprid = reads[rid].prid, tprid;
	bool isrc = reads[rid].isrc, tisrc;
	int16_t shift = reads[rid].shift, tshift;
	while (true) {
		if (debug) fprintf(stderr, "oriprid: %lu; rid: %lu\n", oriprid, rid);
		removeChild(oriprid, rid);
		tprid = reads[oriprid].prid;
		tisrc = reads[oriprid].isrc;
		tshift = reads[oriprid].shift;
		
		reads[oriprid].prid = rid; // reverse the edge
		if (isrc) {
			reads[oriprid].shift = shift;
		} else {
			reads[oriprid].shift = - shift;
		}
		reads[oriprid].isrc = isrc;
		// reads[rid].crid.push_back(oriprid);
		kv_push(uint32_t, reads[rid].crid, oriprid);

		rid = oriprid;
		oriprid = tprid;
		
		if (oriprid == rid) break;

		isrc = tisrc;
		shift = tshift;
	}
	uint32_t noderid;
	reads[_rid].prid = _rid;
	reads[_rid].isrc = 0;
	reads[_rid].shift = 0;
}

inline void connectWithPrid2Ori() {
	uint32_t noderid, frid;
	int fmin_dif;

	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid && reads[rid].crid.n > 0) { //a root node
			bool debug = false;
			// if (rid == 94) debug = true;
			// debug = true;
			queue<uint32_t> q;
			q.push(rid);
			frid = rid;

			fmin_dif = 1024;
			while (!q.empty()) {
				noderid = q.front();
				q.pop();
				if (debug) {
					cout << "noderid: " << noderid << " parent: " << reads[noderid].prid << endl;
				}
				// if (noderid == 76 || noderid == 129 || noderid == 95) debug = true;
				if (isnextrnd[noderid] && prid2[noderid] != noderid && min_dif2[noderid] < fmin_dif) { // no itself
				// if (reads[noderid].prid2 != noderid && reads[noderid].min_dif2 < fmin_dif) { // no itself
					fmin_dif = min_dif2[noderid];
					frid = noderid;
				}
				for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
					// if (reads[reads[noderid].crid.a[i]].isnextrnd) {
					q.push(reads[noderid].crid.a[i]);
					// }
				}
			}
			if (debug) printf("rid in connectWithPrid2: %lu\n", rid);

			if (fmin_dif <= max_dif_thr) {
				if (debug) printf("frid: %lu\n", frid);
				uint32_t fprid = prid2[frid];
				if (debug) printf("fprid: %lu\n", fprid);
				reRootNode(frid, reads[fprid].root); //set frid as root
				// reads[frid].prid = frid;
				// findParent2Ex(frid);
			}

		}
	}
	// cout << "xxxx\n";
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid && prid2[rid] != rid && reads[rid].crid.n > 0 && min_dif2[rid] <= max_dif_thr) { // root node
			reads[rid].prid = prid2[rid];
			reads[rid].shift = shift2[rid];
			reads[rid].isrc = isrc2[rid];
			// reads[reads[rid].prid].crid.push_back(rid);
			kv_push(uint32_t, reads[reads[rid].prid].crid, rid);
		}
	}
}

bool edgecmp(Edge_t a, Edge_t b) {
	return a.dif < b.dif;
}

inline void reRootNode(const uint32_t &_rid) { // set rid as root node
	if (reads[_rid].prid == _rid) return; // rid is root node; do nothing

	uint32_t rid = _rid;
	uint32_t oriprid = reads[rid].prid, tprid;
	bool isrc = reads[rid].isrc, tisrc;
	int16_t shift = reads[rid].shift, tshift;
	while (true) {
		removeChild(oriprid, rid);
		tprid = reads[oriprid].prid;
		tisrc = reads[oriprid].isrc;
		tshift = reads[oriprid].shift;
		
		reads[oriprid].prid = rid; // reverse the edge
		if (isrc) {
			reads[oriprid].shift = shift;
		} else {
			reads[oriprid].shift = - shift;
		}
		reads[oriprid].isrc = isrc;
		// reads[rid].crid.push_back(oriprid);
		kv_push(uint32_t, reads[rid].crid, oriprid);

		rid = oriprid;
		oriprid = tprid;
		
		if (oriprid == rid) break;

		isrc = tisrc;
		shift = tshift;
	}
	uint32_t noderid;
	reads[_rid].prid = _rid;
	reads[_rid].isrc = 0;
	reads[_rid].shift = 0;
}

inline void reLabelRoot(const uint32_t &rid, const uint32_t &root) {
	// cout << "begin reLabelRoot()\n";
	uint32_t noderid = rid;
	queue<uint32_t> q;
	q.push(noderid);
	while (!q.empty()) {
		noderid = q.front();
		// reads[noderid].root = rid;
		reads[noderid].root = root;
		q.pop();
		for (uint32_t i = 0; i < reads[noderid].crid.n; ++i) {
			q.push(reads[noderid].crid.a[i]);
		}
	}
	// cout << "end reLabelRoot()\n";
}

/*inline uint32_t getRoot(const uint32_t &rid) {
	uint32_t p = rid, root = reads[rid].root;
	while (reads[root].root != root) {
		root = reads[root].root;
	}
	return root;
}
*/

// inline uint32_t getRoot(const uint32_t &rid) {
// 	if (reads[rid].root != rid) {
// 		reads[rid].root = getRoot(reads[rid].root);
// 	}
// 	return reads[rid].root;
// }

/*inline uint32_t getRoot(const uint32_t &rid) {
	uint32_t trid = rid, tt, root = reads[rid].root;
	while (root != reads[root].root) {
		root = reads[root].root;
	}
	while (reads[trid].root != root) {
		tt = reads[trid].root;
		reads[trid].root = root;
		trid = tt;
	}
	return root;
}
*/

inline uint32_t getRoot(const uint32_t &rid) {
	uint32_t trid = rid, tt, root = reads[rid].root;
	while (root != reads[root].root) {
		root = reads[root].root;
	}
	while (reads[trid].root != root) {
		tt = reads[trid].root;
		reads[trid].root = root;
		trid = tt;
	}
	return root;
}

#define getRoot(x, r) do { \
		uint32_t trid = (x), tt; \
		(r) = reads[x].root;  \
		while ((r) != reads[(r)].root) { \
			(r) = reads[(r)].root; \
		} \
		while (reads[trid].root != (r)) { \
			tt = reads[trid].root; \
			reads[trid].root = (r); \
			trid = tt; \
		} \
	} while(0)

inline void connectWithPrid2() {
	uint32_t prid, rid;

	CStopWatch tstopwatch;
	tstopwatch.start();

	Edge_v edges;
	kv_init(edges);
	kv_resize(Edge_t, edges, max_rid);

	Edge_t teds;

	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[prid2[rid]].root != reads[rid].root) {
			teds.rid = rid, teds.dif = min_dif2[rid];
			kv_push(Edge_t, edges, teds);
		}
	}

	radix_sort_edge(edges.a, edges.a + edges.n);

	// for (int i = 0; i < 5; ++i) {
	// 	cout << edges.a[i].dif << endl;
	// }
	// sort(edges.a, edges.a + edges.n, edgecmp);

	cout << "Time of radix_sort_edge = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();

	cout << "edges.n: " << edges.n << endl;
	// edges.n: 2283538
	uint32_t root1, root2;

	memset(iscombined, false, sizeof(bool)*max_rid);

	for (size_t i = 0; i < edges.n; ++i) {
		// cout << "i: " << i << endl;

		rid = edges.a[i].rid;
		prid = prid2[rid];
	
		if (!iscombined[reads[rid].root]) { //|| !iscombined[reads[prid].root]) {
		// if (getRoot(rid) != getRoot(prid)) { // not in the same tree
			// cout << "before reRootNode()..\n";
			reRootNode(rid); //
			// cout << "after reRootNode()..\n";
			// 
			reads[rid].prid = prid;
			reads[rid].shift = shift2[rid];
			reads[rid].isrc = isrc2[rid];

			kv_push(uint32_t, reads[prid].crid, rid);
			//
			iscombined[reads[rid].root] = true;
			iscombined[reads[prid].root] = true;
		}
	}

	/*for (size_t i = 0; i < edges.n; ++i) {
		// cout << "i: " << i << endl;

		rid = edges.a[i].rid;
		prid = prid2[rid];
		// if (reads[rid].root != reads[prid].root) { // not in the same tree
		// root1 = getRoot(rid);
		// root2 = getRoot(prid);
		getRoot(rid, root1);
		getRoot(prid, root2);

		if (root1 != root2) {
		// if (getRoot(rid) != getRoot(prid)) { // not in the same tree
			// cout << "before reRootNode()..\n";
			reRootNode(rid); //
			// cout << "after reRootNode()..\n";
			// 
			reads[rid].prid = prid;
			reads[rid].shift = shift2[rid];
			reads[rid].isrc = isrc2[rid];

			kv_push(uint32_t, reads[prid].crid, rid);
			//
			// reLabelRoot(rid, reads[prid].root);
			reads[root1].root = root2;
		}
	}*/

	kv_destroy(edges);

	cout << "Time of connectWithPrid2 = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();
}

inline void cutCircle() {
	bool *visited = new bool[max_rid];
	// memset(visited, false, sizeof(bool)*max_rid);
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		visited[rid] = false;
	}
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid != rid && !visited[rid] && reads[rid].crid.n == 0) { //a node
			uint32_t noderid = rid;
			// visited[noderid] = true;
			while (true) {
				if (visited[reads[noderid].prid]) {
					uint32_t prid = reads[noderid].prid;
					size_t cridsize = reads[prid].crid.n;

					// cout << "noderid: " << noderid << " cridsize: " << cridsize << endl;
					for (size_t i = 0; i < cridsize; ++i) {
						if (reads[prid].crid.a[i] == noderid) {
							reads[prid].crid.a[i] = reads[prid].crid.a[cridsize - 1];
							break;
						}
					}
					// reads[prid].crid.resize(cridsize - 1);
					-- reads[prid].crid.n;

					reads[noderid].prid = noderid;
					reads[noderid].shift = 0;
					reads[noderid].isrc = false;
					reads[noderid].root = noderid;

					break;
				}
				visited[noderid] = true;
				if (noderid == reads[noderid].prid) break; // obtain a root node
				noderid = reads[noderid].prid;
			}
			queue<uint32_t> q;
			q.push(noderid);
			while (!q.empty()) {
				noderid = q.front();
				q.pop();
				visited[noderid] = true;

				for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}
		}
	}

	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid != rid && !visited[rid]) { //a node
			uint32_t noderid = rid;
			uint32_t prid = reads[noderid].prid;
			size_t cridsize = reads[prid].crid.n;

			// cout << "noderid: " << noderid << " cridsize: " << cridsize << endl;
			for (size_t i = 0; i < cridsize; ++i) {
				if (reads[prid].crid.a[i] == noderid) {
					reads[prid].crid.a[i] = reads[prid].crid.a[cridsize - 1];
					break;
				}
			}
			--reads[prid].crid.n;

			reads[noderid].prid = noderid;
			reads[noderid].shift = 0;
			reads[noderid].isrc = false;
			reads[noderid].root = noderid;

			queue<uint32_t> q;
			q.push(noderid);
			while (!q.empty()) {
				noderid = q.front();
				q.pop();
				visited[noderid] = true;

				for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}
		}
	}

	delete[] visited;
}

typedef struct { size_t n; uint8_t a; } uint8bit_v;
#define bit_push(v, fp, x) do {									\
		(v).a += ((x) << ((v).n));										\
		++ (v).n;										\
		if ((v).n == 8) {										\
			fp.write((char*)&(v).a, sizeof(uint8_t));							\
			(v).a = (v).n = 0;							\
		}															\
	} while (0)

#define DNA_push(v, fp, x) do {									\
		(v).a += ((x) << (2*(v).n));										\
		++ (v).n;										\
		if ((v).n == 4) {										\
			fp.write((char*)&(v).a, sizeof(uint8_t));							\
			(v).a = (v).n = 0;							\
		}															\
	} while (0)

void outputORI() {
	stopwatch.resume();
	// char *str = (char*)alloca((L + 1) * sizeof(char));
	char *rcstr = (char*)alloca((L + 1) * sizeof(char));

	FILE *fpcnt = fopen("count.txt", "w");
	FILE *fptrees = fopen("trees.txt", "w");
	FILE *fpenstr = fopen("encodestr.txt", "w");
	FILE *fpenstrcopy = fopen("encodestr.txt.copy", "w");
	FILE *fpenstrcopycopy = fopen("encodestr.txt.copy.copy", "w");
	FILE *fproot = fopen("rootstr.txt", "w");
	char name[100]; 
	// uint8bit_v rootbin;
	// rootbin.a = rootbin.n = 0;
	// sprintf(name, "rootstr.bin");
	// std::ofstream rootOfs(name, std::ios::binary);

	FILE *fporder = NULL;
	if (isorder) {
		string orderfn = folder + string("order.bin");
		// string orderfn = string("order.bin");  // for test
		fporder = fopen(orderfn.c_str(), "wb");
	}

	vector<ROOTNODE_t> rootnodevec;
	for (size_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid && reads[rid].crid.n > 0) { //is the root node && not a leaf node == not a singleton reads
			queue<size_t> q;
			q.push(rid);
			size_t cnt = 0;
			while (!q.empty()) {
				size_t noderid = q.front();		
				q.pop();
				++cnt;
				for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}
			rootnodevec.push_back(ROOTNODE_t(rid, cnt));
		}
	}
	sort(rootnodevec.begin(), rootnodevec.end(), cmp);

	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	// size_t outputnum = 0;

	uint8bit_v dirbin;
	// uint8bit_v rootbin, dirbin;
	// rootbin.n = rootbin.a = 0;
	dirbin.n = dirbin.a = 0;
	sprintf(name, "dir.bin");
	std::ofstream fpdir(name, std::ios::binary);

	size_t ss = sizeof(size_t);
	char *en_str = (char*)alloca((L + 1) * sizeof(char));

	// int max_shift = 0;

	for (size_t i = 0; i < rootnodevec.size(); ++i) {
			size_t cnt = rootnodevec[i].nodecnt;
			fprintf(fpcnt, "%lu\n", cnt);

			size_t rootid = rootnodevec[i].rid;
			size_t num = 0;
			while (rootid < max_rid) {
				// fprintf(stderr, "%lu\n", rootid);
				++ num;
				visited[rootid] = true;

				if (fporder) fwrite(&rootid, ss, 1, fporder);

				if (num == 1) {
					// fprintf(fproot, "%s\n", reads[rootid].str.c_str());
					fprintf(fproot, "%s\n", seq[rootid].seq);
					// for (int j = 0; j < L; ++j) {
					// 	DNA_push(rootbin, rootOfs, seq_nt4_table[(uint8_t)seq[rootid].seq[j]]);
					// }
					rootid = reads[rootid].getChildren();
				} else {
					int dir = 0;
					fprintf(fpenstrcopy, "%s\n", seq[rootid].seq);
					// if (abs(reads[rootid].shift) > max_shift) {
					// 	max_shift = abs(reads[rootid].shift);
					// }
					if (reads[rootid].isrc) {
						dir = 1;
						strcpy(rcstr, seq[rootid].seq);
						reverseComplement(rcstr);
						encode(seq[reads[rootid].prid].seq, rcstr, reads[rootid].shift, en_str);
						fprintf(fpenstrcopycopy, "%s\n%s\n%d\n", seq[reads[rootid].prid].seq, rcstr, reads[rootid].shift);
						fprintf(fpenstrcopycopy, "%s\n", en_str);
						fprintf(fpenstr, "%s\n", en_str);
					} else {
						encode(seq[reads[rootid].prid].seq, seq[rootid].seq, reads[rootid].shift, en_str);
						fprintf(fpenstrcopycopy, "%s\n%s\n%d\n", seq[reads[rootid].prid].seq, seq[rootid].seq, reads[rootid].shift);
						fprintf(fpenstrcopycopy, "%s\n", en_str);
						fprintf(fpenstr, "%s\n", en_str);
					}

					bit_push(dirbin, fpdir, dir);

					if (num >= cnt) {
						rootid = max_rid + 1;
						break;
					}
					size_t child = reads[rootid].getChildren();
					// fprintf(stderr, "child: %lu\n", child);
					int back_step = 0;
					while (child > max_rid) { // no leaf node or all children are visited
						// back
						rootid = reads[rootid].prid;
						// fprintf(stderr, "rootid: %lu\n", rootid);
						child = reads[rootid].getChildren();
						++back_step;
						// ++num;
						// if (num > 20) exit(0);
					}
					// fprintf(stderr, "child: %lu\n", child);
					if (back_step > 0) {
						// fprintf(stderr, "--%d\n", back_step);
						fprintf(fpenstr, "%d\n", back_step);
					}
					rootid = child;
					// if (num > 20) exit(0);
				}
			}
			// exit(0);
			fprintf(fpenstr, "-\n");
		// } 
	}

	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	fclose(fproot);
	// if (rootbin.n > 0) {
	// 	rootOfs.write((char*)&rootbin.a, sizeof(uint8_t));
	// }
	// rootOfs.close();
	
	fclose(fpcnt);
	fclose(fptrees);
	fclose(fpenstr);

	// fprintf(stderr, "max_shift: %d\n", max_shift);
	fclose(fpenstrcopy);
	fclose(fpenstrcopycopy);

	// size_v nrid;
	// if (isquality) {
	// 	kv_init(nrid);
	// }

	uint8bit_v singlebin;
	singlebin.a = singlebin.n = 0;

	sprintf(name, "singleton.bin");
	std::ofstream singleOfs(name, std::ios::binary);

	FILE *fpN = fopen("readsN.txt", "w");
	int ncnt;

	for (size_t rid = 0; rid < max_rid; ++rid) {
		if (!visited[rid]) {
			ncnt = 0;
			for (int i = 0; i < L; ++i) {
				if (seq[rid].seq[i] == 'N') ++ncnt;
			}
			if (ncnt > 0) {
				fprintf(fpN, "%s\n", seq[rid].seq);
				// if (isquality) kv_push(size_t, nrid, rid);
				if (fporder) fwrite(&rid, ss, 1, fporder);
			} else {
				// fprintf(fpsingle, "%s\n", seq[rid].seq);
				for (int i = 0; i < L; ++i) {
					DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
				}
			}
		}
	}
	fclose(fpN);
	if (singlebin.n > 0) {
		singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
	}
	singleOfs.close();

	if (isorder) {
		fclose(fporder);
	}

	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();
}

void outputSingle() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	string rootstrfn = folder + "rootstr.txt";
	FILE *fproot = fopen(rootstrfn.c_str(), "w");

	std::ofstream fporder;
	string orderfn;
	if (isorder) {
		orderfn = folder + string("order.bin");
		// string orderfn = string("order.bin");  // for test
		// fporder = fopen(orderfn.c_str(), "wb");
		fporder.open(orderfn, std::ios::binary);
	}

	vector<ROOTNODE_t> rootnodevec;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid && reads[rid].crid.n > 0) { //is the root node && not a leaf node == not a singleton reads
			queue<uint32_t> q;
			q.push(rid);
			uint32_t cnt = 0;
			while (!q.empty()) {
				uint32_t noderid = q.front();		
				q.pop();
				++cnt;
				for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}
			rootnodevec.push_back(ROOTNODE_t(rid, cnt));
		}
	}
	sort(rootnodevec.begin(), rootnodevec.end(), cmp);

	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;

	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	size_t ss = sizeof(size_t);
	char *en_str = (char*)alloca((L + 1) * sizeof(char));

	for (uint32_t i = 0; i < rootnodevec.size(); ++i) {
			size_t cnt = rootnodevec[i].nodecnt;

			uint32_t rootid = rootnodevec[i].rid;
			uint32_t num = 0;
			while (rootid < max_rid) {
				++ num;
				visited[rootid] = true;

				if (isorder) fporder.write((char*)&rootid, sizeof(uint32_t));

				if (num == 1) {
					fprintf(fproot, "%s\n", seq[rootid].seq);
					rootid = reads[rootid].getChildren();
				} else {
					int dir = 0;
					if (reads[rootid].isrc) {
						dir = 1;
						strcpy(rcstr, seq[rootid].seq);
						reverseComplement(rcstr);
						encode(seq[reads[rootid].prid].seq, rcstr, reads[rootid].shift, en_str);
						fprintf(fpenstr, "%s\n", en_str);
					} else {
						encode(seq[reads[rootid].prid].seq, seq[rootid].seq, reads[rootid].shift, en_str);
						fprintf(fpenstr, "%s\n", en_str);
					}

					bit_push(dirbin, fpdir, dir);

					if (num >= cnt) {
						rootid = max_rid + 1;
						break;
					}
					uint32_t child = reads[rootid].getChildren();

					int back_step = 0;
					while (child > max_rid) { // no leaf node or all children are visited
						// back
						rootid = reads[rootid].prid;
						child = reads[rootid].getChildren();
						++back_step;
					}
					if (back_step > 0) {
						fprintf(fpenstr, "%d\n", back_step);
					}
					rootid = child;
				}
			}
			fprintf(fpenstr, "-\n");
	}

	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	// delete
	rootnodevec.clear();
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		kv_destroy(reads[rid].crid);
	}
	delete[] reads;

	fclose(fproot);
	fclose(fpenstr);

	uint8bit_v singlebin;
	singlebin.a = singlebin.n = 0;

	// sprintf(name, "singleton.bin");
	string singletonfn = folder + "singleton.bin"; 
	std::ofstream singleOfs(singletonfn.c_str(), std::ios::binary);

	string readsnfn = folder + "readsN.txt";
	FILE *fpN = fopen(readsnfn.c_str(), "w");
	int ncnt;

	uint32_t sgcnt = 0;
	for (size_t rid = 0; rid < max_rid; ++rid) {
		if (!visited[rid]) {
			ncnt = 0;
			for (int i = 0; i < L; ++i) {
				if (seq[rid].seq[i] == 'N') ++ncnt;
			}
			if (ncnt > 0) {
				fprintf(fpN, "%s\n", seq[rid].seq);
				if (fporder) fporder.write((char*)&rid, sizeof(uint32_t));
			} else {
				++ sgcnt;
				for (int i = 0; i < L; ++i) {
					DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
				}
			}
		}
	}
	fclose(fpN);
	if (singlebin.n > 0) {
		singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
	}
	singleOfs.close();

	delete[] visited;
	
	cout << "sgcnt: " << sgcnt << endl;
	if (isorder) {
		// fclose(fporder);
		fporder.close();
	}

	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	mstcom::bsc::BSC_compress(encodestrfn.c_str(), (encodestrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(rootstrfn.c_str(), (rootstrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(singletonfn.c_str(), (singletonfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(dirbinfn.c_str(), (dirbinfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(readsnfn.c_str(), (readsnfn+".bsc").c_str(), 64);

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc rootstr.txt.bsc singleton.bin.bsc dir.bin.bsc readsN.txt.bsc";

	if (isorder) {
		mstcom::bsc::BSC_compress(orderfn.c_str(), (orderfn+".bsc").c_str(), 64);
		tarcmd += " order.bin.bsc";
	}

	cout << tarcmd << endl;

	int status = system(tarcmd.c_str());
	if (status < 0) {
		fprintf(stderr, "cmd: %s\t error: %s", tarcmd, strerror(errno));
	}
	if(WIFEXITED(status)) {
		    fprintf(stderr, "normal termination, exit status = %d\n", WEXITSTATUS(status)); //取得cmdstring执行结果
	}
	else if(WIFSIGNALED(status)) {
		fprintf(stderr, "abnormal termination,signal number =%d\n", WTERMSIG(status)); //如果cmdstring被信号中断，取得信号值
	} 
	else if(WIFSTOPPED(status)) {
		fprintf(stderr, "process stopped, signal number =%d\n", WSTOPSIG(status)); //如果cmdstring被信号暂停执行，取得信号值
	}
	cout << "Time of bsc files = " << stopwatch.stop() << std::endl;

}

void outputPE() {
	stopwatch.resume();
	// char *str = (char*)alloca((L + 1) * sizeof(char));
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	// FILE *fpcnt = fopen("count.txt", "w");
	// FILE *fptrees = fopen("trees.txt", "w");
	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	// FILE *fpenstrcopy = fopen("encodestr.txt.copy", "w");
	// FILE *fpenstrcopycopy = fopen("encodestr.txt.copy.copy", "w");
	string rootstrfn = folder + "rootstr.txt";
	FILE *fproot = fopen(rootstrfn.c_str(), "w");
	// char name[100]; 
	// uint8bit_v rootbin;
	// rootbin.a = rootbin.n = 0;
	// sprintf(name, "rootstr.bin");
	// std::ofstream rootOfs(name, std::ios::binary);

	// FILE *fporder = NULL;
	std::ofstream fporder;
	string orderfn = folder + string("order.bin");
		// string orderfn = string("order.bin");  // for test
		// fporder = fopen(orderfn.c_str(), "wb");
	fporder.open(orderfn, std::ios::binary);

	vector<ROOTNODE_t> rootnodevec;
	for (size_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid && reads[rid].crid.n > 0) { //is the root node && not a leaf node == not a singleton reads
			queue<size_t> q;
			q.push(rid);
			size_t cnt = 0;
			while (!q.empty()) {
				size_t noderid = q.front();		
				q.pop();
				++cnt;
				for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}
			rootnodevec.push_back(ROOTNODE_t(rid, cnt));
		}
	}
	sort(rootnodevec.begin(), rootnodevec.end(), cmp);

	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	// size_t outputnum = 0;

	uint8bit_v dirbin;
	// uint8bit_v rootbin, dirbin;
	// rootbin.n = rootbin.a = 0;
	dirbin.n = dirbin.a = 0;
	// sprintf(name, "dir.bin");
	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	size_t ss = sizeof(size_t);
	char *en_str = (char*)alloca((L + 1) * sizeof(char));

	// int max_shift = 0;
	size_t frid;

	size_t *leftmap, *rightmap, leftmapid, rightmapid;
	if (!isorder) {
		leftmap = new size_t[max_rid >> 1];
		rightmap = new size_t[max_rid >> 1];
		leftmapid = rightmapid = 0;
	}

	for (size_t i = 0; i < rootnodevec.size(); ++i) {
			size_t cnt = rootnodevec[i].nodecnt;
			// fprintf(fpcnt, "%lu\n", cnt);

			size_t rootid = rootnodevec[i].rid;
			size_t num = 0;
			while (rootid < max_rid) {
				// fprintf(stderr, "%lu\n", rootid);
				++ num;
				visited[rootid] = true;

				// if (isorder) fwrite(&rootid, ss, 1, fporder);
				int file = 0;
				if (rootid >= (max_rid>>1)) file = 1; // right ends
				bit_push(dirbin, fpdir, file);

				if (isorder) {
					if (file) {
						frid = rootid - (max_rid >> 1);
						fporder.write((char*)&frid, sizeof(uint32_t));
					} else {
						fporder.write((char*)&rootid, sizeof(uint32_t));
					}
				} else { //
					if (file) { // right ends
						rightmap[rightmapid++] = rootid - (max_rid >> 1);
					} else {
						leftmap[rootid] = leftmapid++;
					}
				}


				if (num == 1) {
					// fprintf(fproot, "%s\n", reads[rootid].str.c_str());
					fprintf(fproot, "%s\n", seq[rootid].seq);

					// for (int j = 0; j < L; ++j) {
					// 	DNA_push(rootbin, rootOfs, seq_nt4_table[(uint8_t)seq[rootid].seq[j]]);
					// }
					rootid = reads[rootid].getChildren();
				} else {
					int dir = 0;
					// fprintf(fpenstrcopy, "%s\n", seq[rootid].seq);
					// if (abs(reads[rootid].shift) > max_shift) {
					// 	max_shift = abs(reads[rootid].shift);
					// }
					if (reads[rootid].isrc) {
						dir = 1;
						strcpy(rcstr, seq[rootid].seq);
						reverseComplement(rcstr);
						encode(seq[reads[rootid].prid].seq, rcstr, reads[rootid].shift, en_str);
						// fprintf(fpenstrcopycopy, "%s\n%s\n%d\n", seq[reads[rootid].prid].seq, rcstr, reads[rootid].shift);
						// fprintf(fpenstrcopycopy, "%s\n", en_str);
						fprintf(fpenstr, "%s\n", en_str);
					} else {
						encode(seq[reads[rootid].prid].seq, seq[rootid].seq, reads[rootid].shift, en_str);
						// fprintf(fpenstrcopycopy, "%s\n%s\n%d\n", seq[reads[rootid].prid].seq, seq[rootid].seq, reads[rootid].shift);
						// fprintf(fpenstrcopycopy, "%s\n", en_str);
						fprintf(fpenstr, "%s\n", en_str);
					}

					bit_push(dirbin, fpdir, dir);
					// dir = 0;
					// if (rootid >= (max_rid>>1)) dir = 1;
					// bit_push(dirbin, fpdir, dir);

					if (num >= cnt) {
						rootid = max_rid + 1;
						break;
					}
					size_t child = reads[rootid].getChildren();
					// fprintf(stderr, "child: %lu\n", child);
					int back_step = 0;
					while (child > max_rid) { // no leaf node or all children are visited
						// back
						rootid = reads[rootid].prid;
						// fprintf(stderr, "rootid: %lu\n", rootid);
						child = reads[rootid].getChildren();
						++back_step;
						// ++num;
						// if (num > 20) exit(0);
					}
					// fprintf(stderr, "child: %lu\n", child);
					if (back_step > 0) {
						// fprintf(stderr, "--%d\n", back_step);
						fprintf(fpenstr, "%d\n", back_step);
					}
					rootid = child;
					// if (num > 20) exit(0);
				}
			}
			// exit(0);
			fprintf(fpenstr, "-\n");
		// } 
	}

	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	fclose(fproot);
	// if (rootbin.n > 0) {
	// 	rootOfs.write((char*)&rootbin.a, sizeof(uint8_t));
	// }
	// rootOfs.close();
	
	// fclose(fpcnt);
	// fclose(fptrees);
	fclose(fpenstr);

	// fprintf(stderr, "max_shift: %d\n", max_shift);
	// fclose(fpenstrcopy);
	// fclose(fpenstrcopycopy);

	// size_v nrid;
	// if (isquality) {
	// 	kv_init(nrid);
	// }

	uint8bit_v singlebin;
	singlebin.a = singlebin.n = 0;

	// sprintf(name, "singleton.bin");
	string singletonfn = folder + "singleton.bin"; 
	std::ofstream singleOfs(singletonfn.c_str(), std::ios::binary);

	string leftreadsnfn = folder + "readsN.txt.left";
	FILE *leftfpN = fopen(leftreadsnfn.c_str(), "w");
	int ncnt;

	for (size_t rid = 0; rid < (max_rid>>1); ++rid) {
		if (!visited[rid]) {
			ncnt = 0;
			for (int i = 0; i < L; ++i) {
				if (seq[rid].seq[i] == 'N') ++ncnt;
			}
			if (ncnt > 0) {
				fprintf(leftfpN, "%s\n", seq[rid].seq);
				// if (isquality) kv_push(size_t, nrid, rid);
				if (isorder) {
					fporder.write((char*)&rid, sizeof(uint32_t));
				} else {
					leftmap[rid] = leftmapid++;
					// fwrite(&rid, ss, 1, fporder);
				}
				visited[rid] = true;
			} 
		}
	}
	fclose(leftfpN);

	uint32_t sgcnt = 0;
	for (size_t rid = 0; rid < (max_rid>>1); ++rid) {
		if (!visited[rid]) {
			++ sgcnt;
			// leftmap[rid] = leftmapid++;
			// fprintf(fpsingle, "%s\n", seq[rid].seq);
			for (int i = 0; i < L; ++i) {
				DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
			}
			if (!isorder) {
				leftmap[rid] = leftmapid++;
			}
		}
	}
	cout << "sgcnt: " << sgcnt << endl;

	/////////////
	string rightreadsnfn = folder + "readsN.txt.right";
	FILE *rightfpN = fopen(rightreadsnfn.c_str(), "w");

	for (size_t rid = (max_rid>>1); rid < max_rid; ++rid) {
		if (!visited[rid]) {
			ncnt = 0;
			for (int i = 0; i < L; ++i) {
				if (seq[rid].seq[i] == 'N') ++ncnt;
			}
			if (ncnt > 0) {
				fprintf(rightfpN, "%s\n", seq[rid].seq);
				// if (isquality) kv_push(size_t, nrid, rid);
				if (isorder) {	
					frid = rid - (max_rid >> 1);
					fporder.write((char*)&frid, sizeof(uint32_t));
				} else {
					rightmap[rightmapid++] = rid - (max_rid >> 1);
				}
				visited[rid] = true;
				// fwrite(&rid, ss, 1, fporder);
			} 
		}
	}
	fclose(rightfpN);

	for (size_t rid = (max_rid>>1); rid < max_rid; ++rid) {
		if (!visited[rid]) {
			if (!isorder) {
				rightmap[rightmapid++] = rid - (max_rid >> 1);
			}
			// fprintf(fpsingle, "%s\n", seq[rid].seq);
			for (int i = 0; i < L; ++i) {
				DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
			}
		}
	}
	if (singlebin.n > 0) {
		singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
	}
	singleOfs.close();

	if (!isorder) {
		for (size_t rid = 0; rid < (max_rid>>1); ++rid) {
			fporder.write((char*)&leftmap[rightmap[rid]], sizeof(uint32_t));
		}
		delete[] leftmap;
		delete[] rightmap;
	}
	fporder.close();

	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	mstcom::bsc::BSC_compress(encodestrfn.c_str(), (encodestrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(rootstrfn.c_str(), (rootstrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(dirbinfn.c_str(), (dirbinfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(leftreadsnfn.c_str(), (leftreadsnfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(rightreadsnfn.c_str(), (rightreadsnfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(singletonfn.c_str(), (singletonfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(orderfn.c_str(), (orderfn+".bsc").c_str(), 64);

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc rootstr.txt.bsc dir.bin.bsc singleton.bin.bsc readsN.txt.left.bsc readsN.txt.right.bsc  order.bin.bsc";

	// if (isorder) {
		// tarcmd += " order.bin.bsc";
	// }

	cout << tarcmd << endl;

	int status = system(tarcmd.c_str());
	if (status < 0) {
		fprintf(stderr, "cmd: %s\t error: %s", tarcmd, strerror(errno));
	}
	if(WIFEXITED(status)) {
		    fprintf(stderr, "normal termination, exit status = %d\n", WEXITSTATUS(status)); //取得cmdstring执行结果
	}
	else if(WIFSIGNALED(status)) {
		fprintf(stderr, "abnormal termination,signal number =%d\n", WTERMSIG(status)); //如果cmdstring被信号中断，取得信号值
	} 
	else if(WIFSTOPPED(status)) {
		fprintf(stderr, "process stopped, signal number =%d\n", WSTOPSIG(status)); //如果cmdstring被信号暂停执行，取得信号值
	}
	cout << "Time of bsc files = " << stopwatch.stop() << std::endl;

}

inline void build(uint32_t *&tr, uint32_t &M, uint32_t n) {
    for (M=1; M<=n+1; M<<=1); 
    tr = new uint32_t[M + n + 1];
	memset(tr, 0, sizeof(uint32_t)*(M + n + 1));
}

inline void up(uint32_t *tr, uint32_t x) {
    tr[x] = tr[x<<1] + tr[x<<1|1];
}

inline void update(uint32_t *tr, uint32_t M, uint32_t x, int y) {
    for(tr[x+=M]+=y,x>>=1; x; x>>=1)
        up(tr, x);
}

inline uint32_t query(uint32_t *tr, uint32_t M, uint32_t s, uint32_t t) {
    uint32_t ans = 0;
    s = s+M-1;
    t = t+M+1;
    for(; s^t^1; s>>=1,t>>=1) {
        if(~s&1)ans += tr[s^1];
        if(t&1) ans += tr[t^1];
    }
    return ans;
}

#if false
void outputPEX() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	string rootstrfn = folder + "rootstr.txt";
	FILE *fproot = fopen(rootstrfn.c_str(), "w");

	std::ofstream fporder;
	string orderfn = folder + string("order.bin");
	if (isorder) {
		fporder.open(orderfn, std::ios::binary);
	}

	std::ofstream fpdist;
	string distfn = folder + string("dist.bin");
	fpdist.open(distfn, std::ios::binary);

	vector<ROOTNODE_t> rootnodevec;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid && reads[rid].crid.n > 0) { //is the root node && not a leaf node == not a singleton reads
			queue<uint32_t> q;
			q.push(rid);
			uint32_t cnt = 0;
			while (!q.empty()) {
				uint32_t noderid = q.front();		
				q.pop();
				++cnt;
				for (uint32_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}
			rootnodevec.push_back(ROOTNODE_t(rid, cnt));
		}
	}
	sort(rootnodevec.begin(), rootnodevec.end(), cmp);

	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;
	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	char *en_str = (char*)alloca((L + 1) * sizeof(char));

	uint32_t *ids = new uint32_t[max_rid + 1];
	uint32_t idsid = 0;
	uint32_t *od = new uint32_t[max_rid + 3];

	// FILE *fpid = fopen("idx.txt", "w");
	for (uint32_t i = 0; i < rootnodevec.size(); ++i) {
			uint32_t cnt = rootnodevec[i].nodecnt;

			uint32_t rootid = rootnodevec[i].rid;
			uint32_t num = 0;
			while (rootid < max_rid) {
				ids[idsid] = rootid;
				od[rootid] = idsid ++;

				// fprintf(fpid, "%lu\n", rootid);
				++ num;
				visited[rootid] = true;

				if (num == 1) {
					fprintf(fproot, "%s\n", seq[rootid].seq);
					rootid = reads[rootid].getChildren();
				} else {
					int dir = 0;
					if (reads[rootid].isrc) {
						dir = 1;
						strcpy(rcstr, seq[rootid].seq);
						reverseComplement(rcstr);
						encode(seq[reads[rootid].prid].seq, rcstr, reads[rootid].shift, en_str);
						fprintf(fpenstr, "%s\n", en_str);
					} else {
						encode(seq[reads[rootid].prid].seq, seq[rootid].seq, reads[rootid].shift, en_str);
						fprintf(fpenstr, "%s\n", en_str);
					}

					bit_push(dirbin, fpdir, dir);

					if (num >= cnt) {
						rootid = max_rid + 1;
						break;
					}
					size_t child = reads[rootid].getChildren();
					int back_step = 0;
					while (child > max_rid) { // no leaf node or all children are visited
						// back
						rootid = reads[rootid].prid;
						child = reads[rootid].getChildren();
						++back_step;
					}
					if (back_step > 0) {
						fprintf(fpenstr, "%d\n", back_step);
					}
					rootid = child;
				}
			}
			fprintf(fpenstr, "-\n");
	}

	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	// rootnodevec.clear();

	fclose(fproot);
	
	fclose(fpenstr);

	uint8bit_v singlebin;
	singlebin.a = singlebin.n = 0;

	string singletonfn = folder + "singleton.bin"; 
	std::ofstream singleOfs(singletonfn.c_str(), std::ios::binary);

	string readsnfn = folder + "readsN.txt";
	FILE *fpN = fopen(readsnfn.c_str(), "w");
	int ncnt;

	for (uint32_t rid = 0; rid < (max_rid); ++rid) {
		if (!visited[rid]) {
			ncnt = 0;
			for (int i = 0; i < L; ++i) {
				if (seq[rid].seq[i] == 'N') ++ncnt;
			}
			if (ncnt > 0) {
				ids[idsid] = rid;
				od[rid] = idsid ++;
				// fprintf(fpid, "%lu\n", rid);
				fprintf(fpN, "%s\n", seq[rid].seq);
				visited[rid] = true;
			} 
		}
	}
	fclose(fpN);

	for (uint32_t rid = 0; rid < (max_rid); ++rid) {
		if (!visited[rid]) {
			ids[idsid] = rid;
			od[rid] = idsid ++;
			// fprintf(fpid, "%lu\n", rid);
			for (int i = 0; i < L; ++i) {
				DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
			}
		}
	}
	// fclose(fpid);

	if (singlebin.n > 0) {
		singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
	}
	singleOfs.close();

	// calc dist
	FILE *fp = fopen("dist.txt", "w");
	uint32_t *tr, M;
	build(tr, M, max_rid);
	memset(visited, false, sizeof(bool)*max_rid);

	uint32_t half = max_rid >> 1;
	uint32_t id = 0, num = 0, p, r, dist, orione;
	int32_t res;
	// 在线段数上 使用的是id+1, p+1
	while (id < max_rid && num < half) {
		if (!visited[id]) {
			uint32_t dis = 0;
			if (ids[id] < half) {
				p = od[half + ids[id]];

				r = query(tr, M, id + 1 + 1, p + 1);
				dist = p + 1 - (id + 1 + 1) + 1 - r;
				res = (int32_t)dist;
				if (isorder) {
					fporder.write((char *)&ids[id], sizeof(uint32_t));
				}
			} else {
				p = od[ids[id] - half];
				r = query(tr, M, id + 1 + 1, p + 1);
				dist = p + 1 - (id + 1 + 1) + 1 - r;
				res = 0 - (int32_t)dist;
				if (isorder) {
					orione = ids[id] - half;
					fporder.write((char *)&orione, sizeof(uint32_t));
				}
			}
			fprintf(fp, "%d\n", res);

			fpdist.write((char *)&res, sizeof(int32_t));

			update(tr, M, p + 1, 1);
			visited[p] = true;

			++ num;
		}
		++ id;
	}

	fpdist.close();
	fporder.close();

	delete[] visited;
	delete[] od;
	delete[] ids;

	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	mstcom::bsc::BSC_compress(encodestrfn.c_str(), (encodestrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(rootstrfn.c_str(), (rootstrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(dirbinfn.c_str(), (dirbinfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(readsnfn.c_str(), (readsnfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(distfn.c_str(), (distfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(singletonfn.c_str(), (singletonfn+".bsc").c_str(), 64);
	if (isorder) {
		mstcom::bsc::BSC_compress(orderfn.c_str(), (orderfn+".bsc").c_str(), 64);
	}

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc rootstr.txt.bsc dir.bin.bsc singleton.bin.bsc readsN.txt.bsc dist.bin.bsc";

	if (isorder) {
		tarcmd += " order.bin.bsc";
	}

	cout << tarcmd << endl;

	system(tarcmd.c_str());
	cout << "Time of bsc files = " << stopwatch.stop() << std::endl;

}
#endif

void outputPEX() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	string rootstrfn = folder + "rootstr.txt";
	FILE *fproot = fopen(rootstrfn.c_str(), "w");

	std::ofstream fporder;
	string orderfn = folder + string("order.bin");
	if (isorder) {
		fporder.open(orderfn, std::ios::binary);
	}

	std::ofstream fpdist;
	string distfn = folder + string("dist.bin");
	fpdist.open(distfn, std::ios::binary);

	vector<ROOTNODE_t> rootnodevec;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid && reads[rid].crid.n > 0) { //is the root node && not a leaf node == not a singleton reads
			queue<uint32_t> q;
			q.push(rid);
			uint32_t cnt = 0;
			while (!q.empty()) {
				uint32_t noderid = q.front();		
				q.pop();
				++cnt;
				for (uint32_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}
			rootnodevec.push_back(ROOTNODE_t(rid, cnt));
		}
	}
	sort(rootnodevec.begin(), rootnodevec.end(), cmp);

	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;
	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	char *en_str = (char*)alloca((L + 1) * sizeof(char));

	uint32_t *ids = new uint32_t[max_rid + 1];
	uint32_t idsid = 0;
	uint32_t *od = new uint32_t[max_rid + 3];

	// FILE *fpid = fopen("idx.txt", "w");
	for (uint32_t i = 0; i < rootnodevec.size(); ++i) {
			uint32_t cnt = rootnodevec[i].nodecnt;

			uint32_t rootid = rootnodevec[i].rid;
			uint32_t num = 0;
			while (rootid < max_rid) {
				ids[idsid] = rootid;
				od[rootid] = idsid ++;

				// fprintf(fpid, "%lu\n", rootid);
				++ num;
				visited[rootid] = true;

				if (num == 1) {
					fprintf(fproot, "%s\n", seq[rootid].seq);
					rootid = reads[rootid].getChildren();
				} else {
					int dir = 0;
					if (reads[rootid].isrc) {
						dir = 1;
						strcpy(rcstr, seq[rootid].seq);
						reverseComplement(rcstr);
						encode(seq[reads[rootid].prid].seq, rcstr, reads[rootid].shift, en_str);
						fprintf(fpenstr, "%s\n", en_str);
					} else {
						encode(seq[reads[rootid].prid].seq, seq[rootid].seq, reads[rootid].shift, en_str);
						fprintf(fpenstr, "%s\n", en_str);
					}

					bit_push(dirbin, fpdir, dir);

					if (num >= cnt) {
						rootid = max_rid + 1;
						break;
					}
					size_t child = reads[rootid].getChildren();
					int back_step = 0;
					while (child > max_rid) { // no leaf node or all children are visited
						// back
						rootid = reads[rootid].prid;
						child = reads[rootid].getChildren();
						++back_step;
					}
					if (back_step > 0) {
						fprintf(fpenstr, "%d\n", back_step);
					}
					rootid = child;
				}
			}
			fprintf(fpenstr, "-\n");
	}

	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	// rootnodevec.clear();

	fclose(fproot);
	
	fclose(fpenstr);

	uint8bit_v singlebin;
	singlebin.a = singlebin.n = 0;

	string singletonfn = folder + "singleton.bin"; 
	std::ofstream singleOfs(singletonfn.c_str(), std::ios::binary);

	string readsnfn = folder + "readsN.txt";
	FILE *fpN = fopen(readsnfn.c_str(), "w");
	int ncnt;

	for (uint32_t rid = 0; rid < (max_rid); ++rid) {
		if (!visited[rid]) {
			ncnt = 0;
			for (int i = 0; i < L; ++i) {
				if (seq[rid].seq[i] == 'N') ++ncnt;
			}
			if (ncnt > 0) {
				ids[idsid] = rid;
				od[rid] = idsid ++;
				// fprintf(fpid, "%lu\n", rid);
				fprintf(fpN, "%s\n", seq[rid].seq);
				visited[rid] = true;
			} 
		}
	}
	fclose(fpN);

	for (uint32_t rid = 0; rid < (max_rid); ++rid) {
		if (!visited[rid]) {
			ids[idsid] = rid;
			od[rid] = idsid ++;
			// fprintf(fpid, "%lu\n", rid);
			for (int i = 0; i < L; ++i) {
				DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
			}
		}
	}
	// fclose(fpid);

	if (singlebin.n > 0) {
		singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
	}
	singleOfs.close();

	// calc dist
	FILE *fp = fopen("dist.txt", "w");
	uint32_t *tr, M;
	build(tr, M, max_rid);
	memset(visited, false, sizeof(bool)*max_rid);

	uint32_t half = max_rid >> 1;
	uint32_t id = 0, num = 0, p, r, dist, orione;
	int32_t res;

	uint8bit_v filebin, rangebin;
	filebin.a = filebin.n = 0;
	rangebin.a = rangebin.n = 0;

	string filefn = folder + "file.bin"; 
	std::ofstream fileOfs(filefn.c_str(), std::ios::binary);
	string rangefn = folder + "range.bin"; 
	std::ofstream rangeOfs(rangefn.c_str(), std::ios::binary);

	// 在线段数上 使用的是id+1, p+1
	bool debug = false;
	int file;
	uint16_t disttemp;
	uint32_t small = 0;
	while (id < max_rid && num < half) {
		if (!visited[id]) {
			size_t dis = 0;
			if (ids[id] < half) { // from first file
				file = 0;
				p = od[half + ids[id]];

				r = query(tr, M, id + 1 + 1, p + 1);
				if (debug) cout << "r: " << r << endl;
				dist = p + 1 - (id + 1 + 1) + 1 - r;
				res = (int32_t)dist;

				if (isorder) {
					fporder.write((char *)&ids[id], sizeof(uint32_t));
				}
			} else { // from second file
				file = 1;
				p = od[ids[id] - half];

				r = query(tr, M, id + 1 + 1, p + 1);
				dist = p + 1 - (id + 1 + 1) + 1 - r;
				res = 0 - (int32_t)dist;

				if (isorder) {
					orione = ids[id] - half;
					fporder.write((char *)&orione, sizeof(uint32_t));
				}
			}

			bit_push(filebin, fileOfs, file);

			fprintf(fp, "%d\n", res);

			if (dist < 65536) {
				++ small;
				file = 0;
				disttemp = dist;
				fpdist.write((char *)&disttemp, sizeof(uint16_t));
			} else {
				file = 1;
				fpdist.write((char *)&dist, sizeof(uint32_t));
			}
			bit_push(rangebin, rangeOfs, file);

			update(tr, M, p + 1, 1);
			visited[p] = true;

			++ num;
		}
		++ id;
	}

	fprintf(stderr, "%u\n", small);
	if (filebin.n > 0) { 
		fileOfs.write((char*)&filebin.a, sizeof(uint8_t));
	}
	fileOfs.close();

	if (rangebin.n > 0) { 
		rangeOfs.write((char*)&rangebin.a, sizeof(uint8_t));
	}
	rangeOfs.close();

	fpdist.close();
	fporder.close();

	delete[] visited;
	delete[] od;
	delete[] ids;

	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	mstcom::bsc::BSC_compress(encodestrfn.c_str(), (encodestrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(rootstrfn.c_str(), (rootstrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(dirbinfn.c_str(), (dirbinfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(readsnfn.c_str(), (readsnfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(distfn.c_str(), (distfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(singletonfn.c_str(), (singletonfn+".bsc").c_str(), 64);
	if (isorder) {
		mstcom::bsc::BSC_compress(orderfn.c_str(), (orderfn+".bsc").c_str(), 64);
	}

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc rootstr.txt.bsc dir.bin.bsc singleton.bin.bsc readsN.txt.bsc dist.bin.bsc";

	if (isorder) {
		tarcmd += " order.bin.bsc";
	}

	cout << tarcmd << endl;

	system(tarcmd.c_str());
	cout << "Time of bsc files = " << stopwatch.stop() << std::endl;

}

int mainC() {
	mm192_v a;
	kv_init(a);
	kv_push(mm192_t, a, mm192_t(1, 2, 3));
	kv_push(mm192_t, a, mm192_t(3, 1, 3));
	kv_push(mm192_t, a, mm192_t(1, 2, 3));
	kv_push(mm192_t, a, mm192_t(4, 2, 3));
	kv_push(mm192_t, a, mm192_t(4, 3, 3));

	sort(a.a, a.a + a.n, mm192cmp);

	for (int i = 0; i < a.n; ++i) {
		cout << a.a[i].x.x << ", " << a.a[i].x.y << ", " << a.a[i].y << endl;
	}

	kv_destroy(a);

	return 0;
}

int mainccc() { // test bsc compress file
	string infile = "mstcom_3K49nHxPuo/11.out"; 
	char outfile[] = "11.out.bsc"; 
	char out[] = "11.out.dec"; 
	mstcom::bsc::BSC_compress(infile.c_str(), outfile, 64);
	mstcom::bsc::BSC_decompress(outfile, out);
}

int main(int argc, char *argv[]) {
	cmd = new char[512];
	stopwatch.start();

	init();
	getPars(argc, argv);

	L = getReadsLength(infile.c_str());

	max_dif_thr = L * 0.7;
	// max_dif_thr = L * 0.35;
	// max_dif_thr = L * 0.5;
	cout << "L: " << L << endl;

	getReads(infile.c_str());

	fprintf(stderr, "ispe: %d\nisorder: %d\n", ispe, isorder);

	if (ispe) {
		cout << "max_rid: " << max_rid << endl;
		getRightReads(infile1.c_str());
	}

	// RN 
	if (RN > ((1UL << 32)-1)) {
		cout << "reads numer is to large\n" <<endl;
		exit(1);
	}
	max_rid = (uint32_t)RN;
	// outputfn = argv[3];

	string parfn = folder + "par.txt";
	FILE *fp = fopen(parfn.c_str(), "w");
	fprintf(fp, "%d %d %d\n", L, ispe, isorder);
	// fprintf(fp, "%d\n", ispe);
	// fprintf(fp, "%d\n", isorder);
	if (isorder || ispe) {
		fprintf(fp, "%lu\n", max_rid);
	}
	fclose(fp);

	cout << "max_rid: " << max_rid << endl;
	cout << "Time of read file = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	// fprintf(stderr, "%s\n", seq[max_rid - 1].seq);

/*	kmervecsize = 15;
	kmervec = new int[kmervecsize];

	kmervec[0] = 37;
	for (int i = 1; i < 5; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}
	kmervec[5] = 25;
	for (int i = 6; i < kmervecsize; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}*/

	// kmervecsize = 20;
	// int largekmernum = 5;
	/*kmervecsize = 25;
	int largekmernum = 10;
	kmervec = new int[kmervecsize];

	kmervec[0] = 39;
	for (int i = 1; i < largekmernum; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}
	kmervec[largekmernum] = 29;
	for (int i = largekmernum + 1; i < kmervecsize; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}*/

	/*kmervecsize = 30;
	int largekmernum = 15;
	kmervec = new int[kmervecsize];

	kmervec[0] = 49;
	for (int i = 1; i < largekmernum; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}
	kmervec[largekmernum] = 29;
	for (int i = largekmernum + 1; i < kmervecsize; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}*/

/*	kmervecsize = 30;
	int largekmernum = 15;
	kmervec = new int[kmervecsize];

	kmervec[0] = 55;
	for (int i = 1; i < largekmernum; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}
	kmervec[largekmernum] = 29;
	for (int i = largekmernum + 1; i < kmervecsize; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}
*/

/*	kmervecsize = 30;
	int largekmernum = 15;
	kmervec = new int[kmervecsize];

	kmervec[0] = 60;
	for (int i = 1; i < largekmernum; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}
	kmervec[largekmernum] = 29;
	for (int i = largekmernum + 1; i < kmervecsize; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}
*/

	// kmervecsize = 15;
	kmervecsize = 5;
	int largekmernum = 0;
	kmervec = new int[kmervecsize];

	kmervec[0] = 29;
	for (int i = 1; i < kmervecsize; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}

	
	// -rw-r--r-- 1 yuansliu Research 1398589440 Jul  3 19:38 ERR174310_cb_mstcom10
	/*kmervecsize = 20;
	int largekmernum = 10;
	kmervec = new int[kmervecsize];

	kmervec[0] = 60;
	for (int i = 1; i < largekmernum; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}
	kmervec[largekmernum] = 29;
	for (int i = largekmernum + 1; i < kmervecsize; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}*/

	/*kmervecsize = 30;
	kmervec = new int[kmervecsize];

	kmervec[0] = 45;
	for (int i = 1; i < kmervecsize; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}*/

	/*kmervecsize = 27;
	kmervec = new int[kmervecsize];

	kmervec[0] = 55;
	kmervec[1] = 54;
	kmervec[2] = 53;
	kmervec[3] = 52;

	kmervec[4] = 35;
	kmervec[5] = 34;
	kmervec[6] = 33;
	kmervec[7] = 32;

	kmervec[8] = 29;
	for (int i = 9; i < kmervecsize; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}*/

// 45-15
	// test
	/*kmervecsize = 2;
	int largekmernum = 1;
	kmervec = new int[kmervecsize];

	kmervec[0] = 61;
	kmervec[1] = 15;*/
	
	max_kmer = kmervec[0];
	min_kmer = kmervec[kmervecsize - 1];

	cout << "max_kmer: " << max_kmer << endl;
	cout << "min_kmer: " << min_kmer << endl;

	fprintf(stderr, "kmervecsize: %d\n", kmervecsize);
	for (int i = 0; i < kmervecsize; ++i) {
		fprintf(stderr, "%d, ", kmervec[i]);
	}
	fprintf(stderr, "\n");


	startp[0] = L >= 100 ? L/2 - 32: L/2 - L*32/100;
	endp[0] = L/2 - 1;
	sublen[0] = endp[0] - startp[0] + 1;

	startp[1] = L/2;
	endp[1] = L >= 100 ? L/2 - 1 + 32 : L/2 - 1 + L*32/100;
	sublen[1] = endp[1] - startp[1] + 1;

	cout << "[" << startp[0] << ":" << endp[0] << "]" << endl;
	cout << "[" << startp[1] << ":" << endp[1] << "]" << endl;
	// exit(0);
	// return 0;
	// char str[1<<10];
	// cout << "pls input a string: ";
	// scanf("%s", str);
	reads = new READS_t[max_rid];
	rmtx = new mutex[mask + 1];
	isnextrnd = new bool[max_rid];
	
	treeConstuction();

	delete[] isnextrnd;
	delete[] rmtx;
	
	labelRoot();

	if (ispe) {
		// outputPE();
		outputPEX();
	} else {
		outputSingle();
	}

	sprintf(cmd, "rm -rf %s", folder.c_str());
	system(cmd);
	delete[] cmd;

	return 0;
}

/*************
 * bucket I/O *
 *************/

#define MM_BUCKET_MAGIC "MBI\1"

void bucket_dump(std::string fn, uint64_v *BB) {
	// cout << "begin dump - fn.c_str(): " << fn.c_str() << endl;
	FILE *fp = fopen(fn.c_str(), "wb");
	size_t i, k;
	fwrite(MM_BUCKET_MAGIC, 1, 4, fp);
	size_t tt[3], ss = sizeof(size_t);
	// fprintf(stderr, "sizeof(size_t): %lu\n", sizeof(size_t));
	// res
	// fprintf(stderr, "rec.n: %lu\n", rec.n);
	fwrite(&rec.n, ss, 1, fp);
	for (i = 0; i < rec.n; ++i) {
		// if (i < 5) {
		// 	fprintf(stderr, "%lu %lu %lu\n", rec.a[i].x, rec.a[i].y, rec.a[i].z);
		// }
		tt[0] = rec.a[i].x, tt[1] = rec.a[i].y, tt[2] = rec.a[i].z; 
		// if (i > rec.n - 5) {
		// if (i < 5) {
		// 	fprintf(stderr, "%lu %lu %lu\n", tt[0], tt[1], tt[2]);
		// }
		fwrite(tt, ss, 3, fp);
	}

	for (i = 0; i < (1<<bsize); ++i) {
		uint64_v *b = &BB[i];
		fwrite(&b->n, ss, 1, fp);
		// if (i > (1<<bsize) - 5) {
		// if (i < 5) {
		// 	fprintf(stderr, "\n%lu: ", b->n);
		// }
		for (k = 0; k < b->n; ++k) {
			fwrite(&b->a[k], 8, 1, fp);
			// if (i > (1<<bsize) - 5 && k < 5) {
			// if (i < 5 && k < 5) {
			// 	fprintf(stderr, "%llu, ", b->a[k]);
			// }
		}
	}
	fclose(fp);
	// fprintf(stderr, "end dump\n");
}

bool bucket_load(std::string fn, uint64_v *BB) {
	// cout << "fn.c_str(): " << fn.c_str() << endl;
	FILE *fp = fopen(fn.c_str(), "rb");
	char magic[4];
	if (fread(magic, 1, 4, fp) != 4) return false;
	if (strncmp(magic, MM_BUCKET_MAGIC, 4) != 0) return false;
	
	size_t i;	
	kv_init(rec);
	kv_resize(mmrec_t, rec, 1<<16);
	mmrec_t trec;
	// size_t tt[3] = {0, 0, 0}, recn;
	size_t tt[3], recn, ss = sizeof(size_t);
	fread(&recn, ss, 1, fp);
	// fprintf(stderr, "recn: %lu\n", recn);

	for (i = 0; i < recn; ++i) {
		fread(tt, ss, 3, fp);
		trec.x = tt[0], trec.y = tt[1], trec.z = tt[2];
		kv_push(mmrec_t, rec, trec);
		// if (i > recn - 5) {
		// // if (i < 5) {
		// 	fprintf(stderr, "%lu %lu %lu --\n", tt[0], tt[1], tt[2]);
		// }
	}

	mm128_t tv;
	uint64_t x;
	size_t j;
	size_t bn;

	for (i = 0; i < 1<<bsize; ++i) {
		uint64_v *b = &BB[i];
		fread(&bn, ss, 1, fp);
		// if (i > (1<<bsize) - 5) {
		// if (i < 5) {
		// 	// fprintf(stderr, "%lu\n", bn);
		// 	fprintf(stderr, "\n%lu: ", bn);
		// }
		for (j = 0; j < bn; ++j) {
			fread(&x, 8, 1, fp);
			kv_push(uint64_t, *b, x);
			// if (i > (1<<bsize) - 5 && j < 5) {
			// if (i < 5 && j < 5) {
			// 	fprintf(stderr, "%llu, ", x);
			// }
		}
	}
	fclose(fp);
	// fprintf(stderr, "\n");
	return true;
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
