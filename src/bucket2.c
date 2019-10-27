#include "mstcom.h"

// inline 
void findPrid2(uint64_v *p) {
	
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

// inline 
void findPrid2WithIndex(uint64_v *p) {
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

	int maxnumsearched = 600;

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
		// if (reads[prid2[rid]].root == reads[rid].root) { // now in the same tree
		if (isupdate[rid]) {	
			// if (reads[rid].prid2 == rid) {
			// mindiff = L;
			mindiff = min_dif2[rid];
			
			// if (prid2[rid] == rid) {
			// } else {
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
				if (numsearched >= maxnumsearched) break;
			}

			// if (false && k > 300) { // search in hash table
			if (k > maxnumsearched) { // search in hash table

				for (int shift = 0; shift <= 50; ++shift) { // large than 15 ???
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

			// if (mindiff <= max_dif_thr && (reads[prid2[rid]].root == reads[rid].root || (reads[prid2[rid]].root != reads[rid].root && mindiff < min_dif2[rid]))) {
			if (mindiff < max_dif_thr && mindiff < min_dif2[rid]) {
				 // && reads[prid].prid != rid
				prid2[rid] = prid;
				isrc2[rid] = pdir ^ dir;
				shift2[rid] = ppos - pos;
				min_dif2[rid] = mindiff;

				if (pdir) {
					shift2[rid] = pos - ppos;
				}
			} 
			// else {
			// 	min_dif2[rid] = L;
			// }
		}

		idxstr = curstr.substr(startp[0], sublen[0]);
		addStrMap(&smap[0], idxstr, y);
		idxstr = curstr.substr(startp[1], sublen[1]);
		addStrMap(&smap[1], idxstr, y);
	}

	delStrMap(&smap[0]);
	delStrMap(&smap[1]);
	// smap[0].clear();
	// smap[1].clear();
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
		// if (rec.a[i].z > 10000) {
		// if (rec.a[i].z > 20000) {
		/*if (rec.a[i].z > 15000) {
			recmtx.lock();
			kv_push(mmrec_t, lrec, rec.a[i]);
			recmtx.unlock();
			continue;
		}*/

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

		// findPrid2WithIndex(p);
		findPrid2(p);
		
		kv_destroy(*p);
		free(p);
		// fprintf(stderr, "%lu\n", i);
	}
}

// inline 
void processBuckets2() {
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

//
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

	int maxnumsearched = 1000;

	while (1) {
		uint32_t k = __sync_fetch_and_add(&rid_pthread, 1);
		if (k >= p->n) break;
		// k
		// y = p->a[k];
		// rid = y >> 32;
		// pos = (uint32_t)y>>1;
		// dir = y & 1;
		rid = p->a[k] >> 32;
		pos = (uint32_t)p->a[k]>>1;
		dir = p->a[k] & 1;
		
		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra);
		}
		curstr = (string)(stra);

		if (isupdate[rid]) {
			mindiff = min_dif2[rid];
			
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

				if (numsearched >= maxnumsearched) break;
			}

			// if (false && k > 300) {
			if (k > maxnumsearched) {
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

			if (mindiff < max_dif_thr && mindiff < min_dif2[rid]) {
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

		if (isupdate[rid]) { ///????? why different?
			mindiff = min_dif2[rid];
			
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
									// if (co[trid] != k && abs(tpos - pos) <= mindiff && shift == abs(tpos - pos)) {
									if (abs(tpos - pos) <= mindiff && shift == abs(tpos - pos)) {

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

			if (mindiff < max_dif_thr && mindiff < min_dif2[rid]) {
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
}

void processLargeBuckets2() {
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

		if (n <= 50000) {
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
			
			rid_pthread = 0;
			for (int i = 0; i < nthreads; ++i) {
				threadVec.push_back(std::thread(processLargeBuckets2Fun, smap, &p));
			}

			std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
				thr.join();
			});
			threadVec.clear();

			kv_destroy(p);
			delStrMap(&smap[0]);
			delStrMap(&smap[1]);
			delete[] smap;
		} else {
			/*rid_pthread = rand()% (n - 20001);
			p.n = rid_pthread + 20000;

			StrMap *smap = new StrMap[2];

			kv_init(p);
			kv_resize(uint64_t, p, n);
			for (uint32_t k = 0; k < rid_pthread; ++k) {
				y = b->a[start_a + k];
				strcpy(stra, seq[y >> 32].seq);
				if (y & 1) reverseComplement(stra);

				// co[y >> 32] = k;

				curstr = (string)(stra);
				idxstr = curstr.substr(startp[0], sublen[0]);
				addStrMap(&smap[0], idxstr, y);
				idxstr = curstr.substr(startp[1], sublen[1]);
				addStrMap(&smap[1], idxstr, y);

				kv_push(uint64_t, p, y);
			}
				
			for (int i = 0; i < nthreads; ++i) {
				threadVec.push_back(std::thread(processLargeBuckets2RandomFun, smap, &p));
			}

			std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
				thr.join();
			});
			threadVec.clear();

			kv_destroy(p);
			smap[0].clear();
			smap[1].clear();
			delete[] smap;*/
		}
	}
	kv_destroy(lrec);

	cout << "Time of processLargeBuckets2() = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();
}
