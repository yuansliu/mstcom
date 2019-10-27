#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "kvec.h"
#include "mstcom.h"

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	 //a  b  c   d  e  f  g   h  i  j  k   l  m  n  o
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

/*static inline uint64_t hash64(uint64_t key, uint64_t mask)
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
*/
/*static inline uint64_t hash64_(uint64_t key)
{
	key = (~key + (key << 21)); // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)); // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)); // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31));
	return key;
}
*/
/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param p      minimizers; p->a[i].x is the 2k-bit hash value;
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */

void mm_sketch_wk(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, buf_pos, min_pos;
	mm128_t *buf, min = { UINT64_MAX, UINT64_MAX };

	assert(len > 0 && w > 0 && k > 0);
	buf = (mm128_t*)alloca(w * 16);
	memset(buf, 0xff, w * 16);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			if (++l >= k)
				info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		} else l = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k) kv_push(mm128_t, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1) kv_push(mm128_t, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX)
		kv_push(mm128_t, *p, min);
}

void mm_sketch_x(const char *str, int len, int k, uint32_t rid, mm128_t *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l;
	mm128_t min = { UINT64_MAX, UINT64_MAX };

	// int debug = 1;
	// if (debug) fprintf(stderr, "in mm_sketch(): %s\n", str);
	// if (debug) fprintf(stderr, "mask: %llu\n", mask);
	// len -= 3;
	assert(len > 0 && k > 0);
	// w += f;
	// debug = 0;
	// for (i = 0, l = 0; i < len; i++) {
	for (i = 0, l = 0; i < len; i++) {
		// if (i >= k - 1) debug = 1;
		// if (rid == 13148457)  fprintf(stderr, "%2d, \n", i);
		// if (debug) fprintf(stderr, "%2d:\n", i);

		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer

		// if (debug) fprintf(stderr, "%llu %llu\n", kmer[0], kmer[1]);

		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand

		// if (debug) fprintf(stderr, "%llu %llu %llu\n", kmer[0], kmer[1], kmer[z]);

		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		// if (debug) fprintf(stderr, "%llu %llu\n", info.x, info.y);
		// if (debug) fprintf(stderr, "%lu %lu\n", info.x, info.y);
		if (info.x < min.x) {
			min = info;
		}
		// if (debug) fprintf(stderr, "%lu, %lu\n", min.x, min.y);
	}
	// if (debug) fprintf(stderr, "%lu\n", UINT64_MAX);
	// if (min.x != UINT64_MAX) {
		// kv_push(mm128_t, *p, min);
	// }
	*p = min;
	// fprintf(stderr, "%llu, %llu\n", min.x, min.y);
	// exit(0);
	// if (debug) fprintf(stderr, "\n");

}

// void reverseComplement(char* start) {
// 	char* left = start; // sequence starts
// 	char* right = start + L - 1;
// 	while (right > left) {
// 		char tmp = complement[(uint8_t)*left];
// 		*left = complement[(uint8_t)*right];
// 		*right = tmp;
// 		++left;
// 		--right;
// 	}
// 	if (left == right)
// 		*left = complement[(uint8_t)*left];
// }

// void mm_sketch(const char *str, int &len, int &k, uint32_t &rid, mm128_t *p)
/*void mm_sketch(const char *str, int len, int k, uint32_t rid, mm128_t *p)
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
			info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
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
			info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 1;
		}
		if (info.x < min.x) {
			min = info;
		}
	}

	*p = min;
}

void mm_large_sketch(const char *str, int len, int k, uint32_t rid, mm192_t *p)
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
			info.x.x = hash64(kmer.x, maskX), info.x.y = hash64(kmer.y, maskY), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
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
			info.x.x = hash64(kmer.x, maskX), info.x.y = hash64(kmer.y, maskY), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 1;
		}
		if (info.x < min.x) {
			min = info;
		}
	}

	*p = min;
}
*/
#include "ksort.h"
#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8) 
KSORT_INIT_GENERIC(uint32_t)

#define sort_key_edge(a) ((a).dif)
KRADIX_SORT_INIT(edge, Edge_t, sort_key_edge, 4) 
