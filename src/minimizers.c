#include "mstcom.h"

void calc128MinimizersDupFun() { // before calling this method, MUST set rid_pthread = 0
	mm128_t minimizer;
	int mask = (1<<bsize) - 1, bucketidx, ncnt;
	READS_t *tr;
	while (1) {
		uint32_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;
		// cout << "rid: " << rid << endl;
		tr = &reads[rid];
		tr->prid = rid;
		kv_init(tr->crid);

		reverseReads(seq[rid].seq);

		min128sketch(seq[rid].seq, L, kmer, rid, &minimizer);

		// if (rid == 1634593 || rid == 4) {
		// 		cout << rid << " - " << minimizer.x << endl;
		// 		cout << seq[rid].seq << endl;
		// 	}
			
		bucketidx = minimizer.x & mask;
		mm128_v *p = &B[bucketidx];
		bmtx[bucketidx].lock(); 
		kv_push(mm128_t, *p, minimizer);
		bmtx[bucketidx].unlock();
	}
}

void calc192MinimizersDupFun() { // before calling this method, MUST set rid_pthread = 0
	mm192_t minimizer;
	int mask = (1<<bsize) - 1, bucketidx, ncnt;
	READS_t *tr;
	while (1) {
		uint32_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;
		// cout << "rid: " << rid << endl;
		// if (kmer == max_kmer) {
		// isnextrnd[rid] = true;
		tr = &reads[rid];
		tr->prid = rid;		
		kv_init(tr->crid);

		reverseReads(seq[rid].seq);

		min192sketch(seq[rid].seq, L, kmer, rid, &minimizer);
		bucketidx = ((minimizer.x.x & mask) + (minimizer.x.y & mask)) & mask;
		mm192_v *p = &BL[bucketidx];

		bmtx[bucketidx].lock(); // the lock can be remove by first counting the number of this bucketidx
		kv_push(mm192_t, *p, minimizer);
		bmtx[bucketidx].unlock();
	}
}

void calcMinimizersDup() {
	// bmtx = new mutex[1 << bsize];
	rid_pthread = 0;
	std::vector<thread> threadVec;
	if (kmer <= 31) {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(calc128MinimizersDupFun));
		} 
	} else {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(calc192MinimizersDupFun));
		}
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
	// delete[] bmtx;
}

///
void calc128MinimizersFun() { // before calling this method, MUST set rid_pthread = 0
	mm128_t minimizer;
	int mask = (1<<bsize) - 1, bucketidx;
	while (1) {
		uint32_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;
		
		if (!isnextrnd[rid]) {
			min128sketch(seq[rid].seq, L, kmer, rid, &minimizer);
			bucketidx = minimizer.x & mask;
			mm128_v *p = &B[bucketidx];
			bmtx[bucketidx].lock();
			kv_push(mm128_t, *p, minimizer);
			bmtx[bucketidx].unlock();
		}
	}
}

void calc192MinimizersFun() { // before calling this method, MUST set rid_pthread = 0
	mm192_t minimizer;
	int mask = (1<<bsize) - 1, bucketidx, ncnt;
	while (1) {
		uint32_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;

		if (!isnextrnd[rid]) {
			// __sync_fetch_and_add(&minicount, 1);
			min192sketch(seq[rid].seq, L, kmer, rid, &minimizer);
			bucketidx = ((minimizer.x.x & mask) + (minimizer.x.y & mask)) & mask;
			mm192_v *p = &BL[bucketidx];
			bmtx[bucketidx].lock(); // the lock can be remove by first counting the number of this bucketidx
			kv_push(mm192_t, *p, minimizer);
			bmtx[bucketidx].unlock();
		}
	}
}

void calcMinimizers() {
	// minicount = 0;
	rid_pthread = 0;
	std::vector<thread> threadVec;
	if (kmer <= 31) {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(calc128MinimizersFun));
		}
	} else {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(calc192MinimizersFun));
		}
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
	// cout << "minicount: " << minicount << endl;
}

///
void calc128MaximizersFun() { // before calling this method, MUST set rid_pthread = 0
	mm128_t maximizer;
	int mask = (1<<bsize) - 1, bucketidx, ncnt;
	while (1) {
		uint32_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;
		// cout << "rid: " << rid << endl;
		
		if (!isnextrnd[rid]) {
			// __sync_fetch_and_add(&minicount, 1);
			max128sketch(seq[rid].seq, L, kmer, rid, &maximizer);
			bucketidx = maximizer.x & mask;
			mm128_v *p = &B[bucketidx];
			bmtx[bucketidx].lock();
			kv_push(mm128_t, *p, maximizer);
			bmtx[bucketidx].unlock();
		}
	}
}

void calc192MaximizersFun() { // before calling this method, MUST set rid_pthread = 0
	mm192_t maximizer;
	int mask = (1<<bsize) - 1, bucketidx, ncnt;
	while (1) {
		uint32_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;

		if (!isnextrnd[rid]) {
			// __sync_fetch_and_add(&minicount, 1);
			max192sketch(seq[rid].seq, L, kmer, rid, &maximizer);
			bucketidx = ((maximizer.x.x & mask) + (maximizer.x.y & mask)) & mask;
			mm192_v *p = &BL[bucketidx];
			bmtx[bucketidx].lock(); // the lock can be remove by first counting the number of this bucketidx
			kv_push(mm192_t, *p, maximizer);
			bmtx[bucketidx].unlock();
		}
	}
}

void calcMaximizers() {
	// minicount = 0;
	// bmtx = new mutex[1 << bsize];
	rid_pthread = 0;
	std::vector<thread> threadVec;
	if (kmer <= 31) {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(calc128MaximizersFun));
		}
	} else {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(calc192MaximizersFun));
		}
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
	// delete[] bmtx;
	// cout << "minicount: " << minicount << endl;
}

///
void sort128BucketsFun() {
	uint32_t j, start_a, n;
	mmrec_t ttmmrec;

	while (1) {
		uint32_t bid = __sync_fetch_and_add(&rid_pthread, 1);
		if (bid >= (1ul << bsize)) break;

		mm128_v *b = &B[bid];
		if (b->n > 0) {
			// cout << "b->n: " << b->n << endl;
			radix_sort_128x(b->a, b->a + b->n);

			for (j = 1, n = 1, start_a = 0; j <= b->n; ++j) {
				if (j == b->n || b->a[j].x != b->a[j-1].x) {
					// fprintf(stderr, "j: %lu;\nstart_a: %lu\nn: %lu", j, start_a, n);
					assert(j - start_a == n);

					if (n >= 2) { // the size of a cluster > 1
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

void sort192BucketsFun() {
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

// inline 
void sortBuckets() {
	kv_init(rec);
	kv_resize(mmrec_t, rec, 1<<16);

	rid_pthread = 0;
	std::vector<thread> threadVec;
	if (kmer <= 31) {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(sort128BucketsFun));
		}
	} else {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(sort192BucketsFun));
		}
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();

	// sort rec
	sort (rec.a, rec.a + rec.n, reccmp);
	
	// fprintf(stderr, "rec.n: %lu\n", rec.n);
	// for (int i = 0; i < 10 && i < rec.n; ++i) {
	// 	fprintf(stderr, "%lu, ", rec.a[i].z);
	// }
	// fprintf(stderr, "\n");
}
