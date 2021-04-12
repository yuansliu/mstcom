#include "mstcom.h"

int coutWeight(char *parent, char *child, int8_t _shift) {
	int weight = 0, i, j;
	int8_t shift = _shift;
	if (shift >= 0) {
		for (i = shift, j = 0; i < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				++ weight;
			}
		}
	} else {
		shift = 0 - shift;
		for (i = 0, j = shift; j < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				++ weight;
			}
		}
	}
	// weight = weight * 5 + shift; // SRR870667_2.mstcom.v1.94 // not good; weight*4 is best; 
	// maybe increase the maximum weight when setting weight*5
	// weight = weight * 3 + shift;
	weight = weight * 3 + shift;
	// weight = weight * 2 + shift; // not good
	/*if (_shift != 0) {
		weight = weight * 4 + shift;
	} else {
		weight = weight * 6;
	}*/
	// weight*(8 + 8 + 1) + shift*8 + 8
	// if (shift != 0) ++ weight;
	/*if (_shift < 0) {
		weight += 4;
	}*/

	return weight;
}

int coutmismatchno(char *parent, char *child, int8_t _shift) {
	int weight = 0, i, j;
	int8_t shift = _shift;
	if (shift >= 0) {
		for (i = shift, j = 0; i < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				++ weight;
			}
		}
	} else {
		shift = 0 - shift;
		for (i = 0, j = shift; j < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				++ weight;
			}
		}
	}	
	return weight;
}

// combine two types of edges
void collectEdgesFun(int threadid) {
	uint64_t y, ty;
	uint32_t rid, trid, prid;
	int pos, tpos, ppos;
	int dir, tdir, pdir;
	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));
	int16_t mindiff, oridiff, difflen, shift;
	int searchedNum, maxSearchedNum; 
	char *encstr = (char*)alloca((L + 1) * sizeof(char));

	Edge_v tedgs;

	while (1) {
		uint32_t idx = __sync_fetch_and_add(&rid_pthread, 1);
		if (idx >= minisize) break;

		bool debug = false;
		// maxSearchedNum = 200; // final version
		maxSearchedNum = 400;
		// int numberofedgestored = 50;
		// int numberofedgestored = 20;
		int numberofedgestored = 5;
		// maxSearchedNum = 1000;
		y = mini[idx];
		rid = y >> 32;

		if (debug) cout << idx << "; " << rid << "; " << seq[rid].seq << "; main\n";

		pos = (uint32_t)y >> 2;
		dir = (y >> 1) & 1;
		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra, L);
		}
		// cout << stra << endl;

		kv_init(tedgs);
		// forward search for hamming edges
		searchedNum = 0;
		for (uint32_t i = idx - 1; i >= 0; --i) {
			// cout << (mini[idx]&1) << "; " << (mini[i]&1) << endl;
			if ((mini[idx]&1) != (mini[i]&1)) break;

			// if (i - idx >= 400) break;
			ty = mini[i];
			trid = ty >> 32;
			tpos = (uint32_t)ty >> 2;
			tdir = (ty >> 1) & 1;

			if (tpos != pos) break;

			strcpy(strb, seq[trid].seq);
			if (tdir) {
				reverseComplement(strb, L);
			}
			if (debug) cout << i << "; " << trid << "; " << strb << "\n";

			difflen = coutWeight(strb, stra, tpos - pos);
			// cout << "difflen: " << difflen << endl;

			if (difflen < max_dif_thr) {
				prid = trid;
				ppos = tpos;
				pdir = tdir;

				// always no shift
				// cout << "rid: " << rid << endl;
				// cout << "trid: " << trid << endl;
				// cout << strb << endl;
				// cout << stra << endl;
				// cout << encstr << endl;
				// cout << "difflen - 1: " << difflen - 1 << endl;

				if (difflen > 0) {
					// ++ searchedNum;
					// threadeds[threadid].a[difflen - 1].push(rid, prid);
					if (rid > prid) {
						Edge_t teds(rid, prid, 0, difflen, pdir ^ dir);
						kv_push(Edge_t, tedgs, teds);
					} else {
						Edge_t teds(prid, rid, 0, difflen, pdir ^ dir);
						kv_push(Edge_t, tedgs, teds);
					}
				} else {
					if (rid > prid) {
						recmtx.lock();
						Edge_t teds(rid, prid, 0, 0, pdir ^ dir);
						kv_push(Edge_t, edges, teds);
						recmtx.unlock();	
					} else {
						recmtx.lock();
						Edge_t teds(prid, rid, 0, 0, pdir ^ dir);
						kv_push(Edge_t, edges, teds);
						recmtx.unlock();
					}
				}
			}

			++ searchedNum;
			if (searchedNum > maxSearchedNum) break;
			if (i == 0) break;
		}

		// maxSearchedNum = 400; // final version
		// numberofedgestored = 10;
		// maxSearchedNum = 10; // for test
		// shift edges; search in previous buckets
		searchedNum = 0;
		int fpos = pos; // fpos record the postion of minimizer in the neighbouring bucket
		int bucketsnum = 0;
		uint32_t finalsearchedi = idx;
		if (snum[idx] + 1 < idx)
		for (uint32_t i = idx - snum[idx] - 1; i >= 0; --i) {
			// cout << (mini[idx]&1) << "; " << (mini[i]&1) << endl;
			// cout << "i: " << i << endl;
			finalsearchedi = i;
			if ((mini[idx]&1) != (mini[i]&1)) break;

			// if (i - idx >= 400) break;
			ty = mini[i];
			trid = ty >> 32;

			bool tisleft = true;
			if (trid >= (max_rid >> 1)) tisleft = false;

			tpos = (uint32_t)ty >> 2;
			tdir = (ty >> 1) & 1;

			if (tpos == pos) continue;

			if (tpos != fpos) {
				fpos = tpos;
				++ bucketsnum;
			}
			// if (tpos != fpos) break;
			// if (fpos != pos)

			strcpy(strb, seq[trid].seq);
			if (tdir) {
				reverseComplement(strb, L);
			}
			if (debug) cout << i << "; " << trid << "; " << strb << "\n";

			difflen = coutWeight(strb, stra, tpos - pos);
			// cout << "difflen: " << difflen << endl;
			if (difflen < max_dif_thr) {
				prid = trid;
				ppos = tpos;
				pdir = tdir;
				int16_t shift = (pdir == 1) ? (pos - ppos):(ppos - pos);
				
				// ++ searchedNum;
				// if (mismatchno == 0) {
				// 	threadshifteds[threadid].a[difflen - 1].push(rid, prid, shift, pdir ^ dir);
				// } else {
				if (rid > prid) {
					Edge_t teds(rid, prid, shift, difflen, pdir ^ dir);
					kv_push(Edge_t, tedgs, teds);	
				} else {
					if (!(pdir ^ dir)) {
						shift = 0 - shift;
					}
					Edge_t teds(prid, rid, shift, difflen, pdir ^ dir);
					kv_push(Edge_t, tedgs, teds);
				}
				// }
			}
			// if (abs(tpos - pos) > 0) 
			++ searchedNum;
			
			if (searchedNum > maxSearchedNum) break;
			if (i == 0) break;
		}

		radix_sort_edge(tedgs.a, tedgs.a + tedgs.n);
		
		// for (int i = 0; i < tedgs.n && i < 500; ++i) {
		// for (int i = 0; i < tedgs.n && i < numberofedgestored; ++i) {
		for (int i = 0; i < tedgs.n && (i < numberofedgestored || tedgs.a[i].dif == tedgs.a[numberofedgestored - 1].dif); ++i) {
		// for (int i = 0; i < tedgs.n; ++i) {
			// threadshifteds[threadid].a[difflen - 1].push(rid, prid, shift, pdir ^ dir);
			threadshifteds[threadid].a[tedgs.a[i].dif - 1].push(tedgs.a[i].rid, tedgs.a[i].prid, tedgs.a[i].shift, tedgs.a[i].isrc);
		}

		kv_destroy(tedgs);

	}
}

int coutWeight(char *parent, int pl, char *child, int cl, int8_t _shift) {
	int weight = 0, i, j;
	int8_t shift = _shift;

	// bool debug = false;
	// if (strcmp(child, "CGGTAGACAGATACATAGATCGGTAGACAGACAGATAGATAGATCGGTAGATAGATAAATCGGTAGATAGATAGATAGATAGATAGGTAGGTAGATACAT") == 0) {
	// 	cout << parent << endl;
	// 	cout << child << endl;
	// 	cout << 0 + _shift << endl;
	// 	// exit(0);
	// 	// debug = true;
	// }
		
	if (shift == 0) {
		for (i = 0, j = 0; i < pl && j < cl; ++i, ++j) {
			if (parent[i] != child[j]) {
				++ weight;
			}
		}
	} else 
	if (shift > 0) {
		// case shift >= 0
		//    AATTGCATGC parent
		//      TTGCATGCGA
		for (i = shift, j = 0; i < pl && j < cl; ++i, ++j) {
			if (parent[i] != child[j]) {
				++ weight;
			}
		}
		// if (debug) {
		// 	cout << "pl: " << pl << endl;
		// 	cout << "cl: " << cl << endl;
		// 	cout << "shift: " << 0 + shift << endl;
		// }
		if (pl != cl) {
			if (cl > (pl - shift)) {
				if (cl - (pl - shift) > shift) shift = cl - (pl - shift);
			} else {
				shift = max(L1, L2);
			}
		}
	} else {
		// cast: shift < 0
		//       AATTGCATGC parent
		//    TCGAATTGCA
		shift = 0 - shift;
		for (i = 0, j = shift; i < pl && j < cl; ++i, ++j) {
			if (parent[i] != child[j]) {
				++ weight;
			}
		}
		if (pl != cl) {
			if (pl > (cl - shift)) {
				if (pl - (cl - shift) > shift) shift = pl - (cl - shift);
			} else {
				shift = max(L1, L2);
			}
		}
	} 
	// weight = weight * 5 + shift; // SRR870667_2.mstcom.v1.94 // not good; weight*4 is best; 
	// maybe increase the maximum weight when setting weight*5
	// weight = weight * 3 + shift;
	weight = weight * 3 + shift;
	return weight;
}

void collectEdgesFunL1L2(int threadid) {
	uint64_t y, ty;
	uint32_t rid, trid, prid;
	int pos, tpos, ppos;
	int dir, tdir, pdir, rL, trL;
	int maxL = max(L1, L2);
	char *stra = (char*)alloca((maxL + 1) * sizeof(char));
	char *strb = (char*)alloca((maxL + 1) * sizeof(char));
	int16_t mindiff, oridiff, difflen, shift;
	int searchedNum, maxSearchedNum; 
	char *encstr = (char*)alloca((maxL + 1) * sizeof(char));

	Edge_v tedgs;

	while (1) {
		uint32_t idx = __sync_fetch_and_add(&rid_pthread, 1);
		if (idx >= minisize) break;

		// bool debug = false;
		// maxSearchedNum = 200; // final version
		maxSearchedNum = 400;
		// int numberofedgestored = 50;
		// int numberofedgestored = 20;
		int numberofedgestored = 5;
		// maxSearchedNum = 1000;
		y = mini[idx];
		rid = y >> 32;
		rL = (rid < (max_rid>>1))?L1: L2;

		// if (debug) cout << idx << "; " << rid << "; " << seq[rid].seq << "; main\n";

		pos = (uint32_t)y >> 2;
		dir = (y >> 1) & 1;
		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra, rL);
		}
		// cout << stra << endl;

		kv_init(tedgs);
		// forward search for hamming edges
		// #ifdef false
		searchedNum = 0;
		for (uint32_t i = idx - 1; i >= 0; --i) {
			// cout << (mini[idx]&1) << "; " << (mini[i]&1) << endl;
			if ((mini[idx]&1) != (mini[i]&1)) break;

			// if (i - idx >= 400) break;
			ty = mini[i];
			trid = ty >> 32;
			tpos = (uint32_t)ty >> 2;
			tdir = (ty >> 1) & 1;

			if (tpos != pos) break;

			trL = (trid < (max_rid>>1))?L1: L2;

			strcpy(strb, seq[trid].seq);
			if (tdir) {
				reverseComplement(strb, trL);
			}
			// if (debug) cout << i << "; " << trid << "; " << strb << "\n";

			if (trL == rL) {
				difflen = coutWeight(strb, trL, stra, rL, tpos - pos);
				// cout << "difflen: " << difflen << endl;

				if (difflen < max_dif_thr) {
					prid = trid;
					ppos = tpos;
					pdir = tdir;

					// always no shift
					// cout << "rid: " << rid << endl;
					// cout << "trid: " << trid << endl;
					// cout << strb << endl;
					// cout << stra << endl;
					// cout << encstr << endl;
					// cout << "difflen - 1: " << difflen - 1 << endl;

					if (difflen > 0) {
						// ++ searchedNum;
						// threadeds[threadid].a[difflen - 1].push(rid, prid);
						if (rid > prid) {
							Edge_t teds(rid, prid, 0, difflen, pdir ^ dir);
							kv_push(Edge_t, tedgs, teds);
						} else {
							Edge_t teds(prid, rid, 0, difflen, pdir ^ dir);
							kv_push(Edge_t, tedgs, teds);
						}
					} else {
						if (rid > prid) {
							recmtx.lock();
							Edge_t teds(rid, prid, 0, 0, pdir ^ dir);
							kv_push(Edge_t, edges, teds);
							recmtx.unlock();	
						} else {
							recmtx.lock();
							Edge_t teds(prid, rid, 0, 0, pdir ^ dir);
							kv_push(Edge_t, edges, teds);
							recmtx.unlock();
						}
					}
				}
			}

			++ searchedNum;
			if (searchedNum > maxSearchedNum) break;
			if (i == 0) break;
		}
		// #endif
		// maxSearchedNum = 400; // final version
		// numberofedgestored = 10;
		// maxSearchedNum = 10; // for test
		// shift edges; search in previous buckets
		// #ifdef false
		// maxSearchedNum = 400;
		searchedNum = 0;
		int fpos = pos; // fpos record the postion of minimizer in the neighbouring bucket
		int bucketsnum = 0;
		uint32_t finalsearchedi = idx;
		if (snum[idx] + 1 < idx)
		for (uint32_t i = idx - snum[idx] - 1; i >= 0; --i) {
			// cout << (mini[idx]&1) << "; " << (mini[i]&1) << endl;
			// cout << "i: " << i << endl;
			finalsearchedi = i;
			if ((mini[idx]&1) != (mini[i]&1)) break;

			// if (i - idx >= 400) break;
			ty = mini[i];
			trid = ty >> 32;
			trL = (trid < (max_rid>>1))?L1: L2;

			bool tisleft = true;
			if (trid >= (max_rid >> 1)) tisleft = false;

			tpos = (uint32_t)ty >> 2;
			tdir = (ty >> 1) & 1;

			if (tpos == pos) continue;

			if (tpos != fpos) {
				fpos = tpos;
				++ bucketsnum;
			}
			// if (tpos != fpos) break;
			// if (fpos != pos)

			strcpy(strb, seq[trid].seq);
			if (tdir) {
				reverseComplement(strb, trL);
			}
			// if (debug) cout << i << "; " << trid << "; " << strb << "\n";

			if (tpos - pos != 0) { // always true ??
			// if (rid == 500002 && trL != rL) { // 
			// if (trL != rL && tpos - pos != 0) { // 
				difflen = coutWeight(strb, trL, stra, rL, tpos - pos);
				// cout << "difflen: " << difflen << endl;
				if (difflen < max_dif_thr) {
					prid = trid;
					ppos = tpos;
					pdir = tdir;
					int16_t shift = (pdir == 1) ? (pos - ppos):(ppos - pos);
					
					if (rid > prid) {
						// if (pdir) {
						// 	shift = ppos - pos; // > 0
						// 	if (pdir == 1) {
						// 		shift = 0 - (rL - (trL - shift));
						// 	} else {
						// 		shift = rL - (trL - shift);
						// 	}
						// }
						Edge_t teds(rid, prid, shift, difflen, pdir ^ dir);
						kv_push(Edge_t, tedgs, teds);	
						// if (debug) {
						// 	cout << "111..." << endl;
						// 	cout << seq[prid].seq << endl;
						// 	cout << seq[rid].seq << endl;
						// 	cout << shift << endl;
						// 	cout << (pdir ^ dir) << endl;
						// 	exit(0);
						// }
					} else {
						if (!(pdir ^ dir)) {
							shift = 0 - shift;
						}
						Edge_t teds(prid, rid, shift, difflen, pdir ^ dir);
						kv_push(Edge_t, tedgs, teds);
						// if (debug) {
						// 	cout << "222..." << endl;
						// 	cout << seq[rid].seq << endl;
						// 	cout << seq[prid].seq << endl;
						// 	cout << "shift: " << shift << endl;
						// 	cout << "isrc: " << (pdir ^ dir) << endl;
						// 	exit(0);
						// }
					}
					// }
				}
			}
			// if (abs(tpos - pos) > 0) 
			++ searchedNum;
			
			if (searchedNum > maxSearchedNum) break;
			if (i == 0) break;
		}
		// #endif

		radix_sort_edge(tedgs.a, tedgs.a + tedgs.n);
		
		// for (int i = 0; i < tedgs.n && i < 500; ++i) {
		// for (int i = 0; i < tedgs.n && i < numberofedgestored; ++i) {
		for (int i = 0; i < tedgs.n && (i < numberofedgestored || tedgs.a[i].dif == tedgs.a[numberofedgestored - 1].dif); ++i) {
		// for (int i = 0; i < tedgs.n; ++i) {
			// threadshifteds[threadid].a[difflen - 1].push(rid, prid, shift, pdir ^ dir);
			threadshifteds[threadid].a[tedgs.a[i].dif - 1].push(tedgs.a[i].rid, tedgs.a[i].prid, tedgs.a[i].shift, tedgs.a[i].isrc);
		}

		kv_destroy(tedgs);

	}
}

void collectEdges(string kinds) {
	CStopWatch tstopwatch;
	tstopwatch.start();
	tstopwatch.resume();

	rid_pthread = 1;
	// nthreads = 1;
	std::vector<thread> threadVec;
	// fprintf(stderr, "rec.n: %lu\n", rec.n);
	if (L1 == L2) {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(collectEdgesFun, i));
		}
	} else {
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(collectEdgesFunL1L2, i));
		}
	}
	// cout << "111xxxx22222\n";
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	// cout << "111xxxx222224444\n";
	threadVec.clear();
	// fclose(fptmp);

	cout << "Time of collectEdges() = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();
}
