# include <cstdio>
# include <vector>
# include <assert.h>
# include <thread>
# include <cstdlib>

# include "kvec.h"
# include "khash.h"
# include "mstcom.h"
using namespace std;


#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

mm_idx_t *mm_idx_init(int b) {
	mm_idx_t *mi;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	return mi;
}

const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n) {
	int mask = (1<<mmb) - 1;
	khint_t k;
	mm_idx_bucket_t *b = &mi->B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) return 0;
	k = kh_get(idx, h, minier>>mmb<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) {
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}

void mm_idx_destroy(mm_idx_t *mi) {
	int i;
	if (mi == 0) return;
	for (i = 0; i < 1<<mmb; ++i) {
		free(mi->B[i].p);
		free(mi->B[i].a.a);
		kh_destroy(idx, (idxhash_t*)mi->B[i].h);
	}
	free(mi->B);
	free(mi);
}

size_t mm_idx_gen_idx;

inline void mm_idx_genEx(size_t i) {
	int j, start_a, start_p, n, n_keys;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t*)mmi;
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0) return;

	// sort by minimizer
	radix_sort_128x(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	// b->n is the number of minimizers appearing >1 times
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x != b->a.a[j-1].x) {
			++n_keys;
			if (n > 1) b->n += n;
			n = 1;
		} else ++n;
	}
	h = kh_init(idx); 
	kh_resize(idx, h, n_keys);
	b->p = (uint64_t*)calloc(b->n, 8);//uint64_t *p; // position array for minimizers appearing >1 times

	// create the hash table
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x != b->a.a[j-1].x) {
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			itr = kh_put(idx, h, p->x>>mmb<<1, &absent);
			assert(absent && j - start_a == n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;
				start_p += n;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;
	assert(b->n == start_p);

	// deallocate and clear b->a
	free(b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}

inline void mm_idx_genFun(size_t n) {
	while (1) {
		size_t i = __sync_fetch_and_add(&mm_idx_gen_idx, 1);
		if (i >= n) break;
		// fprintf(stderr, "%lu\n", i);
		mm_idx_genEx(i);
	}
}
 
void mm_idx_gen(int nthreads, mm_idx_t *mi, size_t n) {
	// fprintf(stderr, "mm_idx_gen()\n");
	mm_idx_gen_idx = 0;
	// nthreads = 1;
	std::vector<thread> threadVec;
	for (int i = 0; i < nthreads; ++i) {
		threadVec.push_back(std::thread(mm_idx_genFun, n));
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
}

