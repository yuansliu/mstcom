#ifndef AC_MSTCOM_H
#define AC_MSTCOM_H

# include <string>
# include <cstring>
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
# include <sys/stat.h>
# include "bseq.h"
# include "kvec.h"
# include "libbsc/bsc.h"
# include "liblzma/lzma.h"
using namespace std;

typedef struct mm128_t {
 	uint64_t x, y;

 	mm128_t() {}
 	mm128_t(uint64_t _x, uint64_t _y): x(_x), y(_y) {}

 	bool operator < (const mm128_t a) const {
 		// return x < a.x && y < a.y;
 		return y < a.y || (y == a.y && x < a.x);
 	}

 	bool operator == (const mm128_t a) const {
 		return x == a.x && y == a.y;
 	}

 	bool operator != (const mm128_t a) const {
 		return x != a.x || y != a.y;
 	}
} mm128_t;

typedef struct mm192_t {
	mm128_t x;
 	uint64_t y;

 	mm192_t() {}
 	mm192_t(mm128_t _x, uint64_t _y): x(_x), y(_y) {}
 	mm192_t(uint64_t _x1, uint64_t _x2, uint64_t _y): x(_x1, _x2), y(_y) {}
} mm192_t;

typedef struct {
	uint32_t x, y, z;
} mmrec_t;

typedef struct {
	size_t a, b;
} segmemt_t;

typedef struct Edge_t {
	uint32_t rid;
	int dif;

	Edge_t(uint32_t _rid, int _dif): rid(_rid), dif(_dif) {}
	Edge_t() {};
} Edge_t;

typedef struct { size_t n, m; bool *a; } bool_v;
typedef struct { size_t n, m; Edge_t *a; } Edge_v;
typedef struct { size_t n, m; uint32_t *a; } uint32_v;
typedef struct { size_t n, m; uint64_t *a; } uint64_v;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; mm192_t *a; } mm192_v;
typedef struct { size_t n, m; mmrec_t *a; } mmrec_v;
typedef struct { size_t n, m; segmemt_t *a; } segmemt_v;

struct ROOTNODE_t {
	uint32_t rid, nodecnt;
	ROOTNODE_t() {}
	ROOTNODE_t(uint32_t _rid, uint32_t _nodecnt): rid(_rid), nodecnt(_nodecnt) {}
};

extern unsigned char seq_nt4_table[256];
extern char complement[256];

extern int nthreads;
extern int max_dif_thr;
extern int L;
extern size_t RN;
extern uint32_t max_rid;
extern bseq1_t *seq;
 
// input parameters
extern bool isorder;
extern bool ispe;

extern std::string folder;
// string outputfn;

// --- for minimizer
extern int bsize;
// #define bsize 20
extern mm128_v *B;
extern uint64_v *Bs;
extern mm192_v *BL;

extern int kmer, max_kmer, min_kmer;
extern int *kmervec, kmervecsize;

extern mutex *bmtx;

extern uint32_t rid_pthread;
extern int mask;

extern mm128_t *min128vec;
extern mm192_t *min192vec;
extern uint64_t *mini;
extern uint32_t minisize;

extern int intnum[256];

extern mmrec_v rec, lrec;
extern mutex recmtx;

extern string fnv[65];

typedef struct dup_t {
	uint32_t id;
	bool isrc;
	dup_t() {}
	dup_t(uint32_t _id, bool _isrc): id(_id), isrc(_isrc) {}
} dup_t;

struct READS_t {
	uint32_t prid; //reads id of parent
	uint32_t root, nextid, prechildrid; //root node of current tree 
	uint32_v crid;
	uint32_t dn;
	dup_t *dup;
	// int dif;
	int16_t shift;
	bool isrc;
	READS_t(): nextid(0), prechildrid(0), shift(0), isrc(false) {}
	READS_t(int16_t _shift, bool _isrc): nextid(0), prechildrid(0), shift(_shift), isrc(_isrc) {}

	uint32_t getChildren() {
		if (nextid < crid.n) {
			nextid++;
			return crid.a[nextid - 1];
		}
		return max_rid + 1;
	}
};

extern READS_t *reads;
extern bool *isnextrnd;
extern bool *isupdate;
extern uint32_t *prid2;
extern int16_t *shift2;
extern bool *isrc2;
extern int *min_dif2;
extern bool *iscombined;
extern string infile;
extern string infile1;
extern string outfile;
extern char *cmd;
extern bool debug;

void countTree();
// char complement[256];
// void mm_sketch(const char *str, int len, int k, uint32_t rid, mm128_t *p);
// void mm_large_sketch(const char *str, int len, int k, uint32_t rid, mm192_t *p);
void mm_sketch_wk(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p);
void radix_sort_128x(mm128_t *beg, mm128_t *end);
void radix_sort_edge(Edge_t *beg, Edge_t *end);

// util
void reverseComplement(char* start);
void getPars(int argc, char* argv[]);
void show_usage(const char* prog);

// for reads input 
bool getReads(char const *fn);
bool getRightReads(char const *fn);
int getReadsLength(char const *fn);

// for reads string encode
int16_t diffstrlen(char *parent, char *child, const int16_t &b);
int16_t diffstrlen(char *parent, char *_child, const int16_t &b, const bool &isrc);
void encode(char *parent, char *child, const int16_t &_shift, char *en_str);
int overlap(char *parent, char *child, const int16_t &b);

// for minimizers && bucket sort
void calcMinimizers();
void calcMaximizers();
void calcMinimizersDup();
void sortBuckets();

void calcMinimizers2();

void obtainMiniIdx();

void removeDuplicate();
void collectNext();

// for bucket;
void processBuckets();
void processLargeBuckets();

void processBuckets2();
void processLargeBuckets2();

void removeChild(size_t prid, size_t rid);

// for output
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

void outputSingle();
void outputSingleOrder();
void outputPEX();

// for index dump/load
#define MM_INDEX_MAGIC "MIM\1"
void idxDump(std::string fn);
bool idxLoad(std::string fn);

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

extern CStopWatch stopwatch;

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

inline void min128sketch(const char *str, int len, int k, uint32_t rid, mm128_t *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer = 0;
	int i, l;
	mm128_t min = { UINT64_MAX, UINT64_MAX };

	assert(len > 0 && k > 0);
	for (i = 0, l = 0; i < len; i++) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };

		kmer = (kmer << 2 | c) & mask;           // forward k-mer

		if (++l >= k) {
			// info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
			info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<2 | 0 << 1 | 0;
		}
		if (info.x < min.x) {
			min = info;
		}
	}

	char *rcstr = (char*)alloca((len + 3) * sizeof(char));
	strcpy(rcstr, str);
	reverseComplement(rcstr);

	for (i = 0, l = 0; i < len; i++) {
		int c = seq_nt4_table[(uint8_t)rcstr[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };

		kmer = (kmer << 2 | c) & mask;           // forward k-mer

		if (++l >= k) {
			// info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 1;
			info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<2 | 1 << 1 | 0;
		}
		if (info.x < min.x) {
			min = info;
		}
	}

	*p = min;
}

inline void min192sketch(const char *str, int len, int k, uint32_t rid, mm192_t *p)
{ // k > 31
	uint64_t maskX = (1ULL<<(2*31)) - 1, maskY = (1ULL<<(2*(k-31))) - 1;
	mm128_t kmer = {0, 0};
	int i, l;
	mm192_t min(UINT64_MAX, UINT64_MAX, UINT64_MAX);

	// assert(len > 0 && k > 0);
	for (i = 0, l = 0; i < len; i++) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm192_t info(UINT64_MAX, UINT64_MAX, UINT64_MAX);

		kmer.y = ((kmer.y << 2) | (kmer.x >> 60)) & maskY;
		kmer.x = (kmer.x << 2 | c) & maskX;           // forward k-mer

		if (++l >= k) {
			// info.x.x = hash64(kmer.x, maskX), info.x.y = hash64(kmer.y, maskY), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
			info.x.x = hash64(kmer.x, maskX), info.x.y = hash64(kmer.y, maskY), info.y = (uint64_t)rid<<32 | (uint32_t)i<<2 | 0 << 1 | 0;
		}
		if (info.x < min.x) {
			min = info;
		}
	}

	char *rcstr = (char*)alloca((len + 3) * sizeof(char));
	strcpy(rcstr, str);
	reverseComplement(rcstr);

	for (i = 0, l = 0; i < len; i++) {
		int c = seq_nt4_table[(uint8_t)rcstr[i]];
		mm192_t info(UINT64_MAX, UINT64_MAX, UINT64_MAX);

		kmer.y = ((kmer.y << 2) | (kmer.x >> 60)) & maskY;
		kmer.x = (kmer.x << 2 | c) & maskX;           

		if (++l >= k) {
			// info.x.x = hash64(kmer.x, maskX), info.x.y = hash64(kmer.y, maskY), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 1;
			info.x.x = hash64(kmer.x, maskX), info.x.y = hash64(kmer.y, maskY), info.y = (uint64_t)rid<<32 | (uint32_t)i<<2 | 1 << 1 | 0;
		}
		if (info.x < min.x) {
			min = info;
		}
	}

	*p = min;
}

inline void max128sketch(const char *str, int len, int k, uint32_t rid, mm128_t *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer = 0;
	int i, l;
	mm128_t max = { 0, 0 };

	assert(len > 0 && k > 0);
	for (i = 0, l = 0; i < len; i++) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { 0, 0 };

		kmer = (kmer << 2 | c) & mask;           // forward k-mer

		if (++l >= k) {
			// info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
			info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<2 | 0 << 1 | 0;
		}
		if (max.x < info.x) {
			max = info;
		}
	}

	char *rcstr = (char*)alloca((len + 3) * sizeof(char));
	strcpy(rcstr, str);
	reverseComplement(rcstr);

	for (i = 0, l = 0; i < len; i++) {
		int c = seq_nt4_table[(uint8_t)rcstr[i]];
		mm128_t info = { 0, 0 };

		kmer = (kmer << 2 | c) & mask;           // forward k-mer

		if (++l >= k) {
			// info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 1;
			info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<2 | 1 << 1 | 0;
		}
		if (max.x < info.x) {
			max = info;
		}
	}

	*p = max;
}

inline void max192sketch(const char *str, int len, int k, uint32_t rid, mm192_t *p)
{ // k > 31
	uint64_t maskX = (1ULL<<(2*31)) - 1, maskY = (1ULL<<(2*(k-31))) - 1;
	mm128_t kmer = {0, 0};
	int i, l;
	mm192_t max(0, 0, 0);

	// assert(len > 0 && k > 0);
	for (i = 0, l = 0; i < len; i++) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm192_t info(0, 0, 0);

		kmer.y = ((kmer.y << 2) | (kmer.x >> 60)) & maskY;
		kmer.x = (kmer.x << 2 | c) & maskX;           // forward k-mer

		if (++l >= k) {
			// info.x.x = hash64(kmer.x, maskX), info.x.y = hash64(kmer.y, maskY), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
			info.x.x = hash64(kmer.x, maskX), info.x.y = hash64(kmer.y, maskY), info.y = (uint64_t)rid<<32 | (uint32_t)i<<2 | 0 << 1 | 0;
		}
		if (max.x < info.x) {
			max = info;
		}
	}

	char *rcstr = (char*)alloca((len + 3) * sizeof(char));
	strcpy(rcstr, str);
	reverseComplement(rcstr);

	for (i = 0, l = 0; i < len; i++) {
		int c = seq_nt4_table[(uint8_t)rcstr[i]];
		mm192_t info(0, 0, 0);

		kmer.y = ((kmer.y << 2) | (kmer.x >> 60)) & maskY;
		kmer.x = (kmer.x << 2 | c) & maskX;           

		if (++l >= k) {
			// info.x.x = hash64(kmer.x, maskX), info.x.y = hash64(kmer.y, maskY), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 1;
			info.x.x = hash64(kmer.x, maskX), info.x.y = hash64(kmer.y, maskY), info.y = (uint64_t)rid<<32 | (uint32_t)i<<2 | 1 << 1 | 0;
		}
		if (max.x < info.x) {
			max = info;
		}
	}

	*p = max;
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

int compress_main(int argc, char *argv[]);
int decompress_main(int argc, char *argv[]);

inline int checkIsOriReads(char *seq, int L) {
	int res = 1, len = strlen(seq);

	if (len < L) res = 0;
	else 
	if (len == L) {
		while (*seq != '\0') {
			if (!(*seq >= 'A' && *seq <= 'Z')) {
				res = 0;
				break;
			}
			++ seq;
		}
	} else
	if (len > L) {
		if (seq[L] == '|') res = 2; // reads|consesus reads
		else res = 0;
	}
	// cout << seq << endl;
	// cout << "len: " << len << endl;
	// cout << "res: " << res << endl;
	return res;
}

#endif