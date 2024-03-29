#include "mstcom.h"

bool rgcmp(const mm128_t &a, const mm128_t &b) {
	if (a.x == b.x)
		return ((uint32_t)a.y >> 2) < ((uint32_t)b.y >> 2);
	return a.x < b.x;
}

bool *fmm; //flag for maximizer or minimizer

void process128DupBucketsFun() {
	uint32_t bid, start_a, n, j;
	mm128_t aa;

	int rL = max(L1, L2);
	mm128_v *iv = new mm128_v[rL];
	uint32_t ridvecid, rid;
	mm128_t maximizer;
	
	while (1) {
		uint32_t i = __sync_fetch_and_add(&rid_pthread, 1);
		if (i >= rec.n) break;

		bid = rec.a[i].x;
		start_a = rec.a[i].y;
		n = rec.a[i].z;

		// fprintf(stderr, "%lu %lu %lu\n", rec.a[i].x, rec.a[i].y, rec.a[i].z);
		// if (i > 10) exit(0);

		for (int k = 0; k < rL; ++k) {
			kv_init(iv[k]);
		}
		mm128_v *b = &B[bid];
		for (uint32_t k = 0; k < n; ++k) {
			aa = b->a[start_a + k];
			// if (4 == (aa.y >> 32)) {
			// 	for (int w = 0; w < n; ++w) {
			// 		mm128_t bb = b->a[start_a + w];
			// 		cout << seq[bb.y>>32].seq << "; " << ((uint32_t)bb.y >> 2) << endl;
			// 	}
			// }
			kv_push(mm128_t, iv[(uint32_t)aa.y >> 2], aa);
		}
		ridvecid = __sync_fetch_and_add(&minisize, n);

		for (uint32_t k = rL - 1; k > 0; --k) {
			b = &iv[k];
			if (b->n > 0 ) {
				if (b->n < 5000) {
					for (int i1 = 0; i1 < b->n; ++i1) {
						min128vec[ridvecid ++] = b->a[i1];
					}
				} else {
					for (int i1 = 0; i1 < b->n; ++i1) {
						rid = b->a[i1].y >> 32;

						int trL = (rid < (max_rid >> 1))? L1: L2;
						// if (rid == 20387193) cout << "in maximizer\n";
						max128sketch(seq[rid].seq, trL, kmer, rid, &maximizer);
						b->a[i1] = maximizer;
					}
					//
					sort (b->a, b->a + b->n, rgcmp);

					for (int i1 = 0; i1 < b->n; ++i1) {
						fmm[ridvecid] = true; // label as maximizer
						min128vec[ridvecid ++] = b->a[i1];
					}
				}
				kv_destroy(iv[k]);
			}
		}
	}
	// cout << "1222\n";
	delete[] iv;
}

// inline 
void processDupBuckets() {
	// if (kmer == max_kmer) nthreads = 1; // for test
	rid_pthread = 0;
	// if (kmer < max_kmer) nthreads = 1;

	std::vector<thread> threadVec;
	minisize = 0;
	// nthreads = 1;

	for (int i = 0; i < nthreads; ++i) {
		threadVec.push_back(std::thread(process128DupBucketsFun));
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();

	kv_destroy(rec);
} 

void obtainMiniDupIdx() {
	if (minisize < 2) return;
	// 2^24=16777216
	// 2^25=33554432
	uint64_t y;
	// memset(minidx, 0xff,sizeof(uint32_t)*max_rid);

	int d = 0;
	// if (max_rid <= (1U << 24)) {
		mini[0] = min128vec[0].y | d;
		for (uint32_t i = 1; i < minisize; ++i) {
			if (min128vec[i].x != min128vec[i-1].x || (uint32_t)min128vec[i].y >> 2 != (uint32_t)min128vec[i-1].y >> 2 || fmm[i] != fmm[i-1]) {
				d ^= 1;
			}
			mini[i] = min128vec[i].y | d;
		}
	// } else {
		
	// }
}

// int tt;

void findDupFun() {
	uint64_t y, ty;
	uint32_t rid, trid;
	int dir, tdir;
	int maxL = max(L1, L2);
	char *stra = (char*)alloca((maxL + 1) * sizeof(char));
	char *strb = (char*)alloca((maxL + 1) * sizeof(char));
	// char *stra = new char[L+1];
	// char *strb = new char[L+1];
	while (1) {
		uint32_t idx = __sync_fetch_and_add(&rid_pthread, 1);
		if (idx >= minisize) break;

		// cout << "rid_pthread: " << rid_pthread << endl;
		y = mini[idx];
		rid = y >> 32;
		dir = (y>>1) & 1;
		strcpy(stra, seq[rid].seq);

		int rL = (rid < (max_rid>>1))?L1: L2;
		if (dir) {
			reverseComplement(stra, rL);
		}
		// cout << "idx: " << idx << "; rid: " << rid << endl;
		// bool debug = false;
		// if (rid == 20387193 || rid == 10536560) debug = true;
		// if (debug) {
		// 	cout << "rid: " << rid << endl;
		// 	cout << stra << endl;
		// }
		// for (uint32_t i = idx + 1; i < minisize; ++i) {
		for (uint32_t i = idx - 1; i >= 0; --i) {
			// cout << (mini[idx]&1) << "; " << (mini[i]&1) << endl;
			if ((mini[idx]&1) != (mini[i]&1)) break;

			ty = mini[i];
			trid = ty >> 32;

			// if (debug) cout << trid << endl;
			// cout << "i: " << i << "; trid: " << trid << endl;
			tdir = (ty>>1) & 1;
			strcpy(strb, seq[trid].seq);

			int trL = (trid < (max_rid>>1))?L1: L2;
			if (tdir) {
				reverseComplement(strb, trL);
			}
			// if (debug) cout << strb << endl;
	
			if (rL == trL && strcmp(stra, strb) == 0) {
				isnextrnd[rid] = true;
				kv_push(uint32_t, reads[trid].crid, rid);
				reads[rid].prid = trid;
				reads[rid].isrc = tdir ^ dir;
				reads[rid].shift = 0; 
				break;
			}
			if (i == 0) break;
		}
		// cout << "---\n";
	}
	// delete[] stra;
	// delete[] strb;
}

// inline 
void findDuplicate() {
	// tt = 0;
	rid_pthread = 1;
	std::vector<thread> threadVec;
	// nthreads = 1;
	// cout << "nthreads: " << nthreads << endl; 
	for (int i = 0; i < nthreads; ++i) {
		threadVec.push_back(std::thread(findDupFun));
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
	// cout << "tt: " << tt << endl;
}

bool cmpdupid0(const dup_t &a, const dup_t &b) {
	return a.id < b.id; // original order; small ID first
}

bool cmpdupid1(const dup_t &a, const dup_t &b) {
	return a.id > b.id; // new order; 
}

void processDup() {
	uint32_t trid, dn, prid;
	// cout << "reads[4].crid.n: " << reads[4].crid.n << endl;
	// cout << "reads[4].prid: " << reads[4].prid << endl;
	READS_t *tr, *ntr;
	bool isrc;
	// if (L1 == L2) {
		if (isorder) {
			for (uint32_t rid = 0; rid < max_rid; ++rid) {
				tr = &reads[rid];
				tr->dn = dn = 0;
				// if (reads[rid].prid == rid && reads[rid].crid.n > 0) {
				if (tr->prid == rid && tr->crid.n > 0) {
					trid = rid;
					while (reads[trid].crid.n > 0) {
						++ dn;
						trid = reads[trid].crid.a[0];
						// if (trid == 4) {
						// 	cout << "rid: " << rid << endl;
						// 	cout << seq[trid].seq << endl;
						// 	cout << seq[rid].seq << endl;
						// 	exit(0);
						// }
					}

					tr->dn = dn;
					tr->dup = new dup_t[dn + 1];
					trid = rid;
					dn = 0;
					while (reads[trid].crid.n > 0) {
						prid = trid;
						trid = reads[trid].crid.a[0];
						
						// tr->dup[dn ++] = dup_t(trid, reads[trid].isrc ^ reads[prid].isrc);
						tr->dup[dn ++] = dup_t(trid, false);

						kv_destroy(reads[prid].crid);
						kv_init(reads[prid].crid);
					}

					tr->dup[dn] = dup_t(rid, false);

					sort(tr->dup, tr->dup + dn + 1, cmpdupid0);

					/*for (uint32_t i = 0; i < dn + 1; ++i) {
						// if (tr->dup[i].id == 58344) cout << "HERERERERE!!!!!!2222222" << endl;
						if (tr->dup[i].id == 464744) {
							cout << "HERERERERE!!!!!!22222" << endl;
							cout << "rid: " << rid << endl;
							cout << seq[rid].seq << endl;
							cout << "tr->isrc: " << tr->isrc << endl;
							for (uint32_t w = 0; w < dn + 1; w++) {
								cout << "w: " << w << endl;
								cout << "tr->dup[w].id: " << tr->dup[w].id << endl;
								cout << "tr->dup[w].isrc: " << tr->dup[w].isrc << endl;
								cout << seq[tr->dup[w].id].seq << endl;
							}
							// exit(0);
						}
					}*/

					if (rid != tr->dup[0].id) {
						// move all to reads[rid].dup[0].id
						ntr = &reads[tr->dup[0].id];
						tr->dn = 0;
						ntr->dn = dn;
						ntr->dup = new dup_t[dn + 1];

						ntr->dup[0] = dup_t(tr->dup[0].id, false);
						// if (reads[rid].dup[0].isrc != reads[rid].isrc) {
						// if (strcmp(seq[tr->dup[0].id].seq, seq[rid].seq) != 0) {
						// if (tr->dup[0].id == 20371) cout << "HERERERERE!!!!!!2222222" << endl;
							for (uint32_t i = 1; i < dn + 1; ++i) {
								// if (tr->dup[i].id == 20371) cout << "HERERERERE!!!!!!2222222" << endl;
								isrc = true;
								if (strcmp(seq[tr->dup[0].id].seq, seq[tr->dup[i].id].seq) == 0) isrc = false;

								ntr->dup[i] = dup_t(tr->dup[i].id, isrc);

								/*if (tr->dup[i].id == 77186) {
									cout << "HERERERERE!!!!!!22222" << endl;
									cout << "rid: " << rid << endl;
									cout << seq[rid].seq << endl;
									cout << "tr->isrc: " << tr->isrc << endl;
									for (uint32_t w = 0; w < dn + 1; w++) {
										cout << "tr->dup[w].id: " << tr->dup[w].id << endl;
										cout << "tr->dup[w].isrc: " << tr->dup[w].isrc << endl;
										cout << seq[tr->dup[w].id].seq << endl;
									}
									exit(0);
								}*/
							}
						/*} else {
							for (uint32_t i = 0; i < dn + 1; ++i) {
								if (tr->dup[i].id == 57185) {
									cout << "HERERERERE!!!!!!33333333" << endl;
									cout << "rid: " << rid << endl;
									cout << seq[rid].seq << endl;
									cout << "tr->isrc: " << tr->isrc << endl;
									cout << "tr->dup[0].id: " << tr->dup[0].id << endl;
									cout << "tr->dup[0].isrc: " << tr->dup[0].isrc << endl;
									cout << seq[tr->dup[0].id].seq << endl;
									exit(0);
								}
								ntr->dup[i] = tr->dup[i];
							}
						}*/
						isnextrnd[tr->dup[0].id] = false;
						ntr->prid = tr->dup[0].id;
						isnextrnd[rid] = true;

						kv_destroy(tr->crid);
						delete[] tr->dup;

						/*for (uint32_t i = 0; i < dn + 1; ++i) {
						// if (tr->dup[i].id == 58344) cout << "HERERERERE!!!!!!2222222" << endl;
						if (ntr->dup[i].id == 20371) {
							cout << "HERERERERE!!!!!!22222------" << endl;
							// cout << "rid: " << rid << endl;
							// cout << seq[rid].seq << endl;
							// cout << "tr->isrc: " << tr->isrc << endl;
							for (uint32_t w = 0; w < dn + 1; w++) {
								cout << "w: " << w << endl; 
								cout << "tr->dup[w].id: " << ntr->dup[w].id << endl;
								cout << "tr->dup[w].isrc: " << ntr->dup[w].isrc << endl;
								cout << seq[ntr->dup[w].id].seq << endl;
							}
							exit(0);
						}
						}*/
					} else {
						for (uint32_t i = 1; i < dn + 1; ++i) {
							isrc = true;
							if (strcmp(seq[tr->dup[0].id].seq, seq[tr->dup[i].id].seq) == 0) isrc = false;

							tr->dup[i].isrc = isrc;
						}
					}
					// change to the smallest idx

				}
			}
		} else {
			for (uint32_t rid = 0; rid < max_rid; ++rid) {
				dn = 0;
				tr = &reads[rid];
				// if (reads[rid].prid == rid && reads[rid].crid.n > 0) {
				if (tr->prid == rid && tr->crid.n > 0) {
					trid = rid;
					while (reads[trid].crid.n > 0) {
						++ dn;
						trid = reads[trid].crid.a[0];
						// if (trid == 4) {
						// 	cout << "rid: " << rid << endl;
						// 	cout << seq[trid].seq << endl;
						// 	cout << seq[rid].seq << endl;
						// 	exit(0);
						// }
					}

					tr->dup = new dup_t[dn];
					trid = rid;
					dn = 0;
					while (reads[trid].crid.n > 0) {
						prid = trid;
						trid = reads[trid].crid.a[0];
						reads[trid].prid = trid;
						// tr->dup[dn ++] = dup_t(trid, reads[trid].isrc);
						tr->dup[dn ++] = dup_t(trid, false);
						kv_destroy(reads[prid].crid);
						kv_init(reads[prid].crid);
					}
					for (uint32_t i = 0; i < dn; ++i) {
						isrc = true;
						if (strcmp(seq[rid].seq, seq[tr->dup[i].id].seq) == 0) isrc = false;

						tr->dup[i].isrc = isrc;
					}
				}
				tr->dn = dn;
			}
		}
	// } 
	#ifdef false
	else {
		for (uint32_t rid = 0; rid < max_rid; ++rid) {
			tr = &reads[rid];
			tr->dn = dn = 0;
			// if (reads[rid].prid == rid && reads[rid].crid.n > 0) {
			if (tr->prid == rid && tr->crid.n > 0) {
				trid = rid;
				while (reads[trid].crid.n > 0) {
					++ dn;
					trid = reads[trid].crid.a[0];
					// if (trid == 4) {
					// 	cout << "rid: " << rid << endl;
					// 	cout << seq[trid].seq << endl;
					// 	cout << seq[rid].seq << endl;
					// 	exit(0);
					// }
				}

				tr->dn = dn;
				tr->dup = new dup_t[dn + 1];
				trid = rid;
				dn = 0;
				while (reads[trid].crid.n > 0) {
					prid = trid;
					trid = reads[trid].crid.a[0];
					
					// tr->dup[dn ++] = dup_t(trid, reads[trid].isrc ^ reads[prid].isrc);
					tr->dup[dn ++] = dup_t(trid, false);

					kv_destroy(reads[prid].crid);
					kv_init(reads[prid].crid);
				}

				tr->dup[dn] = dup_t(rid, false);

				if (L1 > L2) {
					sort(tr->dup, tr->dup + dn + 1, cmpdupid0);
				} else
					sort(tr->dup, tr->dup + dn + 1, cmpdupid1);

				if (rid != tr->dup[0].id) {
					// move all to reads[rid].dup[0].id
					ntr = &reads[tr->dup[0].id];
					tr->dn = 0;
					ntr->dn = dn;
					ntr->dup = new dup_t[dn + 1];

					ntr->dup[0] = dup_t(tr->dup[0].id, false);
					// if (reads[rid].dup[0].isrc != reads[rid].isrc) {
					// if (strcmp(seq[tr->dup[0].id].seq, seq[rid].seq) != 0) {
					// if (tr->dup[0].id == 20371) cout << "HERERERERE!!!!!!2222222" << endl;
					for (uint32_t i = 1; i < dn + 1; ++i) {
						// if (tr->dup[i].id == 20371) cout << "HERERERERE!!!!!!2222222" << endl;
						isrc = true;
						int l1 = (tr->dup[0].id < (max_rid>>1))?L1: L2;
						int l2 = (tr->dup[i].id < (max_rid>>1))?L1: L2;
						if (strncmp(seq[tr->dup[0].id].seq, seq[tr->dup[i].id].seq, min(l1, l2)) == 0) isrc = false;

						ntr->dup[i] = dup_t(tr->dup[i].id, isrc);
					}
					isnextrnd[tr->dup[0].id] = false;
					ntr->prid = tr->dup[0].id;
					isnextrnd[rid] = true;

					kv_destroy(tr->crid);
					delete[] tr->dup;
				} else {
					for (uint32_t i = 1; i < dn + 1; ++i) {
						isrc = true;
						int l1 = (tr->dup[0].id < (max_rid>>1))?L1: L2;
						int l2 = (tr->dup[i].id < (max_rid>>1))?L1: L2;
						if (strncmp(seq[tr->dup[0].id].seq, seq[tr->dup[i].id].seq, min(l1, l2)) == 0) isrc = false;

						tr->dup[i].isrc = isrc;
					}
				}
				// change to the smallest idx

			}
		}
	}
	#endif
}

void removeDuplicate() {
	CStopWatch tstopwatch;
	tstopwatch.start();
	stopwatch.start();

	kmer = 31; 

	B = (mm128_v*)calloc(1 << bsize, sizeof(mm128_v));
	min128vec = new mm128_t[max_rid];
	fmm = new bool[max_rid];
	memset(fmm, false, sizeof(bool)*max_rid);
	for (int i = 0; i < (1 << bsize); ++i) {
		kv_init(B[i]);
	}

	// cout << "kmer: " << kmer << endl;
	calcMinimizersDup();

	cout << "Time of calcMinimizersDup() = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();
	// cout << "after calcMinimizers()\n";
	// cout << "kmer: " << kmer << endl;
	sortBuckets(); 

	cout << "Time of sortBuckets() = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();
	// cout << "after sortBuckets()\n";
	// cout << "kmer: " << kmer << endl;
	// processBuckets();
	processDupBuckets();

	cout << "Time of processBuckets() = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();

	obtainMiniDupIdx();
	cout << "after obtainMiniDupIdx()\n";
	// cout << "kmer: " << kmer << endl;

	for (int i = 0; i < (1 << bsize); ++i) {
		kv_destroy(B[i]);
	}
	free(B); B = NULL;
	delete[] min128vec; min128vec = NULL;
	delete[] fmm; fmm = NULL;

	cout << "before findDuplicate()" << endl;
	findDuplicate();
	cout << "after findDuplicate()" << endl;

	processDup();
	// cout << a << endl;
	cout << "Time of removeDuplicate(): " << stopwatch.stop() << std::endl;
	stopwatch.resume();
}
