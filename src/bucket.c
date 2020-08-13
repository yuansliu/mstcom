#include "mstcom.h"

void process128BucketsFun() {
	uint32_t bid, start_a, n;
	mm128_t aa;

	// char *temp_str = (char*)alloca((L + 1) * sizeof(char));
	// cout << "text xxx\n";

	mm128_v *iv = new mm128_v[L];
	uint32_t ridvecid;
	
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

		mm128_v *b = &B[bid];
		for (uint32_t k = 0; k < n; ++k) {
			aa = b->a[start_a + k];
			kv_push(mm128_t, iv[(uint32_t)aa.y >> 2], aa);
		}

		ridvecid = __sync_fetch_and_add(&minisize, n);

		for (uint32_t k = L - 1; k >= 0; --k) {
			if (iv[k].n > 0) {
				for (int i1 = 0; i1 < iv[k].n; ++i1) {
					min128vec[ridvecid ++] = iv[k].a[i1];
				}
				kv_destroy(iv[k]);
			}
			if (k == 0) break;
		}
	}
	// cout << "1222\n";
	delete[] iv;
}

void process192BucketsFun() {
	uint32_t bid, start_a, n;
	mm192_t aa;

	// char *temp_str = (char*)alloca((L + 1) * sizeof(char));
	// cout << "text xxx\n";

	mm192_v *iv = new mm192_v[L];
	uint32_t ridvecid;
	
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

		mm192_v *b = &BL[bid];
		for (uint32_t k = 0; k < n; ++k) {
			aa = b->a[start_a + k];
			kv_push(mm192_t, iv[(uint32_t)aa.y >> 2], aa);
		}

		ridvecid = __sync_fetch_and_add(&minisize, n);

		// for (uint32_t k = 0; k < L; ++k) {
		for (uint32_t k = L - 1; k >= 0; --k) {
			if (iv[k].n > 0) {
				for (int i1 = 0; i1 < iv[k].n; ++i1) {
					min192vec[ridvecid ++] = iv[k].a[i1];
				}
			}
			if (k == 0) break;
		}

		for (int k = 0; k < L; ++k) {
			kv_destroy(iv[k]);
		}

	}
	// cout << "1222\n";
	delete[] iv;
}

// inline 
void processBuckets() {
	// if (kmer == max_kmer) nthreads = 1; // for test
	rid_pthread = 0;
	// if (kmer < max_kmer) nthreads = 1;

	std::vector<thread> threadVec;
	minisize = 0;
	// nthreads = 1;

	if (kmer <= 31) {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(process128BucketsFun));
		}
	} else {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(process192BucketsFun));
		}
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();

	kv_destroy(rec);
} 

void obtainMiniIdx() {
	if (minisize < 2) return;
	// 2^24=16777216
	// 2^25=33554432
	uint64_t y;
	// memset(minidx, 0xff,sizeof(uint32_t)*max_rid);

	// uint32_v bucketsizevec;
	// kv_init(bucketsizevec);
 	
 	memset(snum, 0, sizeof(uint32_t)*max_rid);
	int d = 0;
	// if (max_rid <= (1U << 24)) {
		if (kmer <= 31) {
			// uint32_t bucketsize = 1;
			mini[0] = min128vec[0].y | d;
			for (uint32_t i = 1; i < minisize; ++i) {
				if (min128vec[i].x != min128vec[i-1].x) {
					d ^= 1;
					// kv_push(uint32_t, bucketsizevec, bucketsize);
					// bucketsize = 1;
				}
				mini[i] = min128vec[i].y | d;

				// the same minimizer and the two minimizer having the same position
				// pos = (uint32_t)y >> 2;
				if (min128vec[i].x == min128vec[i-1].x && ((uint32_t)min128vec[i].y >> 2) == ((uint32_t)min128vec[i-1].y >> 2) ) {
					snum[i] = snum[i-1] + 1;
					// ++ bucketsize;
				}

				// if (min128vec[i].x == min128vec[i-1].x && ((uint32_t)min128vec[i].y >> 2) != ((uint32_t)min128vec[i-1].y >> 2) ) {
					// kv_push(uint32_t, bucketsizevec, bucketsize);
					// bucketsize = 1;
				// }
				// if ((mini[i] >> 32) == 4) cout << "kmer: " << kmer << "; rid 4 exist\n";
			}
			// kv_push(uint32_t, bucketsizevec, bucketsize);
		} else {
			// uint32_t bucketsize = 1;
			mini[0] = min192vec[0].y | d;
			for (uint32_t i = 1; i < minisize; ++i) {
				if (min192vec[i].x != min192vec[i-1].x) {
					d ^= 1;
					// kv_push(uint32_t, bucketsizevec, bucketsize);
					// bucketsize = 1;
				}
				mini[i] = min192vec[i].y | d;

				if (min192vec[i].x == min192vec[i-1].x && ((uint32_t)min192vec[i].y >> 2) == ((uint32_t)min192vec[i-1].y >> 2) ) {
					snum[i] = snum[i-1] + 1;
					// ++ bucketsize;
				}

				// if (min192vec[i].x == min192vec[i-1].x && ((uint32_t)min192vec[i].y >> 2) == ((uint32_t)min192vec[i-1].y >> 2) ) {
					// kv_push(uint32_t, bucketsizevec, bucketsize);
					// bucketsize = 1;
				// }
			}
		}
	// } else {
		
	// }

	/*FILE *fp = fopen("rebucket.txt.28.minimizer", "w");
	fprintf(fp, "kmer-%d\n", kmer);
	for (uint32_t i = 0; i < minisize; ++i) {
		char *str = (char*)alloca((L + 1) * sizeof(char));
		uint32_t rid = mini[i] >> 32;
		strcpy(str, seq[rid].seq);
		if ((mini[i]>>1) & 1) {
			reverseComplement(str);
		}
		fprintf(fp, "%s %u %d\n", str, ((uint32_t)mini[i] >> 2), (mini[i]&1));
	}
	fclose(fp);*/
	/*for (uint32_t i = 0; i < minisize; ++i) {
		if ((mini[i] >> 32) == 4) {
			cout << "kmer = " << kmer << "; rid 4 is found\n";
			cout << "i: " << i << endl;
			// for (int k = i + 1; k >= 0; --k) {
			// 	if ((mini[i]&1) != (mini[k]&1)) break;
			// 	cout << (mini[i] >> 32) << "; " << seq[mini[i] >> 32].seq << "; " << ((uint32_t)mini[i] >> 2) << endl;
			// }
		}
	}*/

	// sort(bucketsizevec.a, bucketsizevec.a + bucketsizevec.n);
	/*uint32_t count = 0;
	for (uint32_t i = 0; i < bucketsizevec.n; ++i) {
		count += bucketsizevec.a[i];
	}
	cout << "count: " << count << endl;
	cout << "minisize: " << minisize << endl;
	for (uint32_t i = bucketsizevec.n - 1, j = 0; j < 30; ++j, --i) {
		fprintf(stderr, "%u, ", bucketsizevec.a[i]);
	}
	fprintf(stderr, "\n");*/
	// exit (0);
}

void idxDump(std::string fn) {
	// cout << "---idxDump---\n";
	// cout << "file name: " << fn << endl;
	FILE *fp = fopen(fn.c_str(), "wb");
	fwrite(MM_INDEX_MAGIC, 1, 4, fp);
	fwrite(&minisize, sizeof(uint32_t), 1, fp);
	// if (kmer == 21) {
		// cout << "minisize: " << minisize << endl;
		// for (int i = 3279250; i < 3279260; ++i) {
		// for (int i = minisize-1; i > minisize-10; --i) {
			// cout << mini[i] << "; " << (mini[i]>>32) << endl;
		// }

	// }
	fwrite(mini, sizeof(uint64_t), minisize, fp);
	fwrite(snum, sizeof(uint32_t), minisize, fp);

	// for (int i = minisize-1; i > minisize-10; --i) {
	// // for (int i = 0; i < 10; ++i) {
	// 	cout << snum[i] << endl;
	// }
	// fclose(fp);
}

bool idxLoad(std::string fn) {
	// cout << "---idxLoad---\n";
	// cout << "file name: " << fn << endl;
	// cout << "rid_pthread: " << rid_pthread << endl;
	FILE *fp = fopen(fn.c_str(), "rb");
	// cout << "111\n";
	char *magic = new char[4];
	// cout << "111\n";
	if (fread(magic, 1, 4, fp) != 4) return false;
	// cout << "111\n";
	if (strncmp(magic, MM_INDEX_MAGIC, 4) != 0) return false;

	// cout << "rid_pthread: " << rid_pthread << endl;
	fread(&minisize, sizeof(uint32_t), 1, fp);
	// cout << "rid_pthread: " << rid_pthread << endl;
	// cout << "minisize: " << minisize << endl;
	fread(mini, sizeof(uint64_t), minisize, fp);
	fread(snum, sizeof(uint32_t), minisize, fp);
	// if (kmer == 21) {
		// cout << "minisize: " << minisize << endl;
		// for (int i = 3279250; i < 3279260; ++i) {
		// for (int i = 0; i < 10; ++i) {
		// for (int i = minisize-1; i > minisize-10; --i) {
			// cout << mini[i] << "; " << (mini[i]>>32) << endl;
		// }

	// }
	// for (int i = 0; i < 5; ++i) cout << mini[i] << endl;
	// for (int i = minisize-1; i > minisize-10; --i) {
	// // for (int i = 0; i < 10; ++i) {
	// 	cout << snum[i] << endl;
	// }
	fclose(fp);
	delete[] magic;
	// cout << "rid_pthread: " << rid_pthread << endl;
	return true;
}

