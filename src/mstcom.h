#ifndef AC_MSTCOM_H
#define AC_MSTCOM_H

# include <string>
# include <cstring>
# include <cstdio>
# include <iostream>
# include <fstream>
# include <vector>
# include <queue>
# include <tuple>
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
# include <sys/types.h>
# include "bseq.h"
# include "kvec.h"
# include "khash.h"
# include "libbsc/bsc.h"
# include "liblzma/lzma.h"
using namespace std; 

#define _ENCODEORI
// #define _ENCODECONTIG

typedef struct mm128_t {
 	uint64_t x, y;

 	mm128_t(): x(0), y(0) {}
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

// typedef struct Edge_t {
// 	uint32_t rid;
// 	int dif;

// 	Edge_t(uint32_t _rid, int _dif): rid(_rid), dif(_dif) {}
// 	Edge_t() {};
// } Edge_t;

typedef struct Edge_t {
	uint32_t rid, prid;
	int16_t shift;
	int dif;
	bool isrc;

	Edge_t(uint32_t _rid, uint32_t _prid, int16_t _shift, int _dif, bool _isrc): rid(_rid), prid(_prid), shift(_shift), dif(_dif), isrc(_isrc) {}
	Edge_t() {};
} Edge_t;

typedef struct MstEdge_t {
	uint32_t rid, prid;
	int16_t shift;
	bool isrc;

	MstEdge_t(uint32_t _rid, uint32_t _prid, int16_t _shift, bool _isrc): rid(_rid), prid(_prid), shift(_shift), isrc(_isrc) {}
	MstEdge_t() {};
} MstEdge_t;

typedef struct { size_t n, m; char *a; } char_v;
typedef struct { size_t n, m; bool *a; } bool_v;
typedef struct { size_t n, m; Edge_t *a; } Edge_v;
typedef struct { size_t n, m; MstEdge_t *a; } MstEdge_v;
typedef struct { size_t n, m; uint8_t *a; } uint8_v;
typedef struct { size_t n, m; uint32_t *a; } uint32_v;
typedef struct { size_t n, m; uint64_t *a; } uint64_v;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; mm192_t *a; } mm192_v;
typedef struct { size_t n, m; mmrec_t *a; } mmrec_v;
typedef struct { size_t n, m; segmemt_t *a; } segmemt_v;

struct ROOTNODE_t {
	uint32_t rid, nodecnt;
	double weight;
	ROOTNODE_t() {}
	ROOTNODE_t(uint32_t _rid, uint32_t _nodecnt): rid(_rid), nodecnt(_nodecnt) {}
	ROOTNODE_t(uint32_t _rid, uint32_t _nodecnt, double _w): rid(_rid), nodecnt(_nodecnt), weight(_w) {}
};

typedef struct {
	mm128_v a;   // (hash value, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for hash value appearing >1 times
	void *h;     // hash table indexing _p_ and hash value appearing once
} mm_idx_bucket_t;

typedef struct {
	uint32_t n, b;  // number of reference sequences
	mm_idx_bucket_t *B;
} mm_idx_t;

const int HASH_BITS = 20;
const int HASH_SIZE = 1 << HASH_BITS;
const int HASH_SIZE_MINUS_ONE = HASH_SIZE - 1;

/*********************************
 * Sort and generate hash tables *
 *********************************/

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

mm_idx_t *mm_idx_init(int b);
const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n);
void mm_idx_destroy(mm_idx_t *mi);

// inline uint64_t maRushPrime1HashSimplified(const char *str, int K) {
inline uint64_t hashFunc(const char *str, int K) {
	std::uint64_t hash = K;
	for (std::uint32_t j = 0; j < K/4; ) {
		std::uint32_t k;
		memcpy(&k, str, 4);
		k += j++;
		hash ^= k;
		hash *= 171717;
		str += 4;
	}
	return hash;
}

void selectNoMismatchEdges();

void constructCotigIndex(char *ctgstr, uint32_t len, mm_idx_t *ht, int K, int k1);
tuple<uint32_t, uint32_t> reconstructtree();

extern uint32_v gridvec; // store reads IDs in a group (grid)

extern unsigned char seq_nt4_table[256];
extern char complement[256];

extern Edge_v edges;
extern int nthreads;
extern int max_dif_thr;
extern int L, L1, L2;
extern int nk;
extern size_t RN;
extern uint32_t max_rid;
extern bseq1_t *seq;

// input parameters
extern bool isorder;
extern bool ispe;
extern FILE *fppar;

extern std::string folder;
extern std::string hammingedgfolder;
extern std::string shiftedgfolder;
// string outputfn;

// --- for minimizer
extern int bsize;
// #define bsize 20
extern mm128_v *B;
extern uint64_v *Bs;
extern mm192_v *BL;

extern int kmer, max_kmer, min_kmer, gkmer;
extern int *kmervec, kmervecsize;

extern mutex *bmtx;

extern uint32_t rid_pthread;
extern uint32_t begrid;
extern uint32_t endrid;
extern int mask;

extern mm128_t *min128vec;
extern mm192_t *min192vec;
extern uint64_t *mini;
extern uint32_t *snum; // the number of minimizers at the same position
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

bool cmpdupid0(const dup_t &a, const dup_t &b);
bool cmpdupid1(const dup_t &a, const dup_t &b);

struct READS_t {
	uint32_t prid; //reads id of parent
	uint32_t nextid; //root node of current tree 
	uint32_v crid;
	uint32_t dn;
	dup_t *dup;
	int16_t shift;
	bool isrc;
	READS_t(): nextid(0), shift(0), isrc(false) {}
	READS_t(int16_t _shift, bool _isrc): nextid(0), shift(_shift), isrc(_isrc) {}

	uint32_t getChildren() {
		// if (prid == 73367) {
		// 	cout << "in getChildren()...\n";
		// 	cout << "nextid: " << nextid << endl;
		// 	cout << "crid.n: " << crid.n << endl;
		// }
		if (nextid < crid.n) {
			nextid++;
			return crid.a[nextid - 1];
		}
		// if (prid == 73367) {
		// 	cout << "----------in getChildren()...\n";
		// 	cout << "return: " << max_rid + 1 << endl;
		// }
		return max_rid + 1;
	}
};

// #define nDFSENCODING
#define DFSENCODING

extern READS_t *reads;
extern bool *isnextrnd;
extern string infile;
extern string infile1;
extern string outfile;
extern char *cmd;
extern bool debug;

extern uint32_t maxedges;

typedef struct shiftEdgesFile_t {
	FILE *fp, *fpshift;
	string fn, shiftfn;
	int16_t dis, *shift; // shift is composed of shift and isrc
	uint64_t *a;
	uint32_t n, m, maxreadnum; // m is total number of reads in this file

	void init(string lfd, int16_t _dis) {
		dis = _dis;
		fn = lfd + to_string(dis) + ".dis";
		shiftfn = lfd + to_string(dis) + ".shift";
		n = m = 0;
		a = new uint64_t[maxedges];
		shift = new int16_t[maxedges];
		fp = NULL;
	}

	void push(uint32_t _a, uint32_t _b, int16_t _shift, bool isrc) { // an edge a-->b  _b is parent
		// cout << "dis: " << dis << "; n: " << n << endl;
		a[n] = (uint64_t)_a<<32 | _b;
		shift[n] = _shift << 1 | isrc;
		++ n;
		++ m;

		if (n >= maxedges) {
			n = 0;
			fp = fopen(fn.c_str(), "ab");
			fpshift = fopen(shiftfn.c_str(), "ab");
			fwrite(a, sizeof(uint64_t), maxedges, fp);
			fwrite(shift, sizeof(int16_t), maxedges, fpshift);
			fclose(fp);
			fclose(fpshift);
		}
	}

	void reopen(uint32_t _maxreadnum = (1<<30)) {
		// fclose(fp);
		// cout << "fn: " << fn << endl;
		// cout << "begin reopen()..." << endl;
		if (n > 0) {
			fp = fopen(fn.c_str(), "ab");
			fpshift = fopen(shiftfn.c_str(), "ab");
			fwrite(a, sizeof(uint64_t), n, fp);
			fwrite(shift, sizeof(int16_t), n, fpshift);
			fclose(fp);
			fclose(fpshift);
			n = 0;
		}
		n = 0;
		delete[] a; a = NULL;
		delete[] shift; shift = NULL;

		maxreadnum = _maxreadnum;
		fp = fopen(fn.c_str(), "rb");
		fpshift = fopen(shiftfn.c_str(), "rb");
		// cout << "end reopen()..." << endl;
	}

	bool get() {
		if (m == 0) return false;
		uint32_t curreadnum = maxreadnum;
		// cout << "m: " << m << endl;
		if (curreadnum > m) {
			curreadnum = m;
		}
		// cout << "curreadnum: " << curreadnum << endl;

		if (NULL != a) {
			delete[] a;
		}
		if (NULL != shift) {
			delete[] shift;
		}
		a = new uint64_t[curreadnum + 2];
		shift = new int16_t[curreadnum + 2];
		// cout << "111\n";
		if (NULL == fp || fread(a, sizeof(uint64_t), curreadnum, fp) == 0) {
		// cout << "111222\n";
			n = 0;
			return false;
		}
		// cout << "after read a" << endl;
		m -= curreadnum;
		fread(shift, sizeof(int16_t), curreadnum, fpshift);
		// cout << "after read shift" << endl;
		n = curreadnum;
		return true;
	}

	// ~shiftEdgesFile_t() {
	void clear() {
		if (NULL != fp) fclose(fp);
		if (NULL != fpshift) fclose(fpshift);
		if (NULL != a) delete[] a;
		if (NULL != shift) delete[] shift;
	}
} shiftEdgesFile_t;

typedef struct threadShiftEdges_t {
	shiftEdgesFile_t *a;
	string lfd;

	void init(uint32_t tid) {
		lfd = shiftedgfolder + "tid" + to_string(tid);
		mode_t mode = 0777;
		int nError = mkdir(lfd.c_str(), mode);
		if (nError != 0) {
			fprintf(stderr, "Failed to creat the folder '%s'\n", lfd.c_str());
			exit(EXIT_FAILURE);
		}
		lfd += "/";
		a = new shiftEdgesFile_t[max_dif_thr];
		for (int i = 0; i < max_dif_thr; ++i) {
			a[i].init(lfd, i);
		}
	}

	void removefd() {
		string cmd = "rm -rf " + lfd;
		system(cmd.c_str());
	}

	~threadShiftEdges_t() {
		delete[] a;
	}
} threadShiftEdges_t;

extern threadShiftEdges_t *threadshifteds;

void reRootNode(const uint32_t &_rid);
void reRootNodeL1L2(const uint32_t &_rid);
void countTree();
// char complement[256];
// void mm_sketch(const char *str, int len, int k, uint32_t rid, mm128_t *p);
// void mm_large_sketch(const char *str, int len, int k, uint32_t rid, mm192_t *p);
void mm_sketch_wk(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p);
void radix_sort_128x(mm128_t *beg, mm128_t *end);
void radix_sort_edge(Edge_t *beg, Edge_t *end);
void radix_sort_mstedge(MstEdge_t *beg, MstEdge_t *end);

// util
void reverseComplement(char* start, int L);
void reverseReads(char* start, int L);
void getPars(int argc, char* argv[]);
void show_usage(const char* prog);

// for reads input 
bool getReads(char const *fn);
bool getRightReads(char const *fn);
bool getRightReads(char const *fn, int LL);
int getReadsLength(char const *fn);

// for reads string encode
// int16_t diffstrlen(char *parent, char *child, const int16_t &b);
// int16_t diffstrlen(char *parent, char *_child, const int16_t &b, const bool &isrc);
void encode(char *parent, char *child, const int16_t &_shift, char *en_str);
uint32_t encode(uint32_t rid, char *en_str);
int encode_v4(uint32_t rid, char *en_str);
int encode_v3(uint32_t rid, char *en_str);
int overlap(char *parent, char *child, const int16_t &b);

// for minimizers && bucket sort
void calcMinimizers();
void calcMinimizersIdx(int);
void calcGroupMinimizers();
void calcMaximizers();
void calcGroupMaximizers();
void calcMinimizersDup();
void sortBuckets();

void calcMinimizers2();

void obtainMiniIdx();

void removeDuplicate();
void collectNext();
void collectShiftEdges();
void collectHammingEdges();
void collectEdges(string kinds);
void constructSubstrIndexEdges();
void collectGroupNext();
void reRootNode(const uint32_t &_rid);
void labelRoot();

void printTree(uint32_t root); // for test

// for bucket;
void processBuckets();
void processLargeBuckets();

void processBuckets2();
void processLargeBuckets2();

void removeChild(size_t prid, size_t rid);

void constructContigs();

mm_idx_t *mm_idx_init(int b);
const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n);
void mm_idx_destroy(mm_idx_t *mi);
void mm_idx_generation(int n_threads, mm_idx_t *mi);

// for output
typedef struct { size_t n, m; uint8_t a; } uint8bit_v;
typedef struct { size_t n, m; uint16_t a; } uint16bit_v;

#define bit_push(v, fp, x) do {									\
		(v).a += ((x) << ((v).n));										\
		++ (v).m; \
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

#define base_push(v, fp, x) do {									\
		(v).a += ((x) << (3*(v).n));										\
		++ (v).n;										\
		if ((v).n == 5) {										\
			fp.write((char*)&(v).a, sizeof(uint16_t));							\
			(v).a = (v).n = 0;							\
		}															\
	} while (0)

bool constructContigs(uint32_t &rid, vector<uint32_t> &gridvec, FILE *fpenstr, std::ofstream &fpdig, std::ofstream &fpdir, uint8bit_v &dirbin);

void outputSingle();
void outputSingleDFS();
void outputSingleDFS_v1();
void outputSingleDFS_v2();
void outputSingleV1();
void outputSingleOrder();
void outputPEX();
void outputPEOrder();

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
	reverseComplement(rcstr, len);

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

inline void min128sketch(const char *str, int len, int k, uint32_t rid, mm128_t *p, int idx)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer = 0;
	int i, l;
	mm128_t min = { UINT64_MAX, UINT64_MAX };
	mm128_v kmerarr;
	kv_init(kmerarr);

	assert(len > 0 && k > 0);
	for (i = 0, l = 0; i < len; i++) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };

		kmer = (kmer << 2 | c) & mask;           // forward k-mer

		if (++l >= k) {
			// info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
			info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<2 | 0 << 1 | 0;
			kv_push(mm128_t, kmerarr, info);
		}
	}

	char *rcstr = (char*)alloca((len + 3) * sizeof(char));
	strcpy(rcstr, str);
	reverseComplement(rcstr, len);

	for (i = 0, l = 0; i < len; i++) {
		int c = seq_nt4_table[(uint8_t)rcstr[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };

		kmer = (kmer << 2 | c) & mask;           // forward k-mer

		if (++l >= k) {
			// info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 1;
			info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<2 | 1 << 1 | 0;
			kv_push(mm128_t, kmerarr, info);
		}
	}

	radix_sort_128x(kmerarr.a, kmerarr.a + kmerarr.n);

	*p = kmerarr.a[idx];
	kv_destroy(kmerarr);
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
	reverseComplement(rcstr, len);

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
	reverseComplement(rcstr, len);

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
	reverseComplement(rcstr, len);

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