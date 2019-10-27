#include "mstcom.h"

// inline 
/*void linkReads2TreeWithIndex(uint64_v *p) {
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

	// int maxnumsearched = 1000;
	int maxnumsearched = 300;

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

		bool tdebug = false;
		// if (rid == 7053579) tdebug = true;

		// if (reads[rid].prid == rid) { // || (reads[rid].prid != rid && abs(reads[rid].shift) > 1 && reads[rid].dif > 5)) {

			// if (reads[rid].prid == rid) {
			// 	mindiff = L;
			// } else {
			// 	mindiff = reads[rid].dif;
			// }
		// if (min_dif2[rid] > 0) {	
			mindiff = min_dif2[rid];

			int numsearched = 0;
			for (int q = k - 1; q >= 0; --q) {
				tpos = (uint32_t) p->a[q] >> 1;
				tdir = p->a[q] & 1;
				// if (k - q >= 100 && tpos - pos >= 5) {
				if (k - q >= maxnumsearched && tpos - pos >= 25) { // tpos > pos ?
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
				if (numsearched >= maxnumsearched) break;
			}

			// if (false && mindiff != 0 && k > 100) {
			if (mindiff != 0 && k > maxnumsearched) {
				for (int shift = 0; shift <= 35 && mindiff != 0; ++shift) {
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
			// if (mindiff <= max_dif_thr && (reads[rid].prid == rid || (reads[rid].prid != rid && mindiff < reads[rid].dif))) {
			if (mindiff < max_dif_thr && mindiff < min_dif2[rid]) {

				prid2[rid] = prid;
				isrc2[rid] = pdir ^ dir;
				shift2[rid] = ppos - pos;
				min_dif2[rid] = mindiff;

				if (pdir) {
					shift2[rid] = pos - ppos;
				}

				if (tdebug) {

				// if (rid == 1539203 || rid == 3443365) {
					cout << rid << "-->" << prid << endl;
					cout << seq[rid].seq << endl;
					cout << seq[prid].seq << endl;
					cout << isrc2[rid] << endl;
					cout << shift2[rid] << endl;
					cout << mindiff << endl;
					cout << "---\n";
				}
				
				if (mindiff == 0) {
					isnextrnd[rid] = false;

					kv_push(uint32_t, reads[prid].crid, rid);
					reads[rid].prid = prid;
					reads[rid].isrc = pdir ^ dir;
					reads[rid].shift = ppos - pos;
					// reads[rid].shift = (pdir == 1)? (pos - ppos):(ppos - pos);

					if (pdir) {
						reads[rid].shift = pos - ppos;
					}
				}
			}
		// }

		idxstr = curstr.substr(startp[0], sublen[0]);
		addStrMap(&smap[0], idxstr, y);
		idxstr = curstr.substr(startp[1], sublen[1]);
		addStrMap(&smap[1], idxstr, y);

	}

	delStrMap(&smap[0]);
	delStrMap(&smap[1]);
	delete[] smap;
}
*/

bool ridcmp(const uint64_t &a, const uint64_t &b) {
	return (a>>32) > (b>>32);
}

void find128NextFun() {
	uint64_t y, ty;
	uint32_t rid, trid, tval, idx, p;
	int pos, tpos;
	int dir, tdir;
	int16_t difflen, fshift;

	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));

	// cout << "rid_pthread in: " << rid_pthread << endl;
	while (1) {
		p = __sync_fetch_and_add(&rid_pthread, 1);
		// cout << "p: " << p << endl;
		// cout << "isnextrnd[p]: " << isnextrnd[p] << endl;
		if (p >= pend) break;


		if (!isnextrnd[p]) continue;

		if (kmer == max_kmer) {
			// cout << "new reads[" << p << "]\n";
			reads[p].next = new uint64_v[max_dif_thr];
			for (int i = 0; i < max_dif_thr; ++i) {
				kv_init(reads[p].next[i]);
			}
		}

		idx = minidx[p];
		if (idx >= ridvecsize) continue;

		y = min128vec[idx].y;
		rid = y >> 32;
		pos = (uint32_t)y >> 1;
		dir = y & 1;
		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra);
		}

		for (uint32_t i = idx + 1; i < ridvecsize; ++i) {
			if (min128vec[idx].x != min128vec[i].x) break;

			ty = min128vec[i].y;
			tpos = (uint32_t) ty >> 1;

			if (tpos - pos >= max_dif_thr) break;

			if (tpos - pos <= 30 && i - idx > (1U<<16)) break;

			tdir = ty & 1;
			trid = ty >> 32;
			strcpy(strb, seq[trid].seq);
			if (tdir) {
				reverseComplement(strb);
			}

			//
			difflen = diffstrlen(strb, stra, tpos - pos);
			if (difflen < max_dif_thr) {
				fshift = tpos - pos;
				
				tval = ((uint32_t)fshift<<2) | ((tdir ^ dir) << 1) | tdir;
				ty = (uint64_t)trid<<32 | tval;					

				// ty = (uint64_t)trid<<32 | tval << 16 | (uint32_t) difflen;
				// if (tdir) {
				// 	fshift = pos - tpos;
				// }
				// ty = (uint64_t)trid<<32 | (uint32_t)fshift<<1 | (tdir ^ dir);

				kv_push(uint64_t, reads[rid].next[difflen], ty);
				// if (tdir) {
				// 	fshift = pos - tpos;
				// }
				// ty = (uint64_t)trid<<32 | (uint32_t)fshift<<1 | (tdir ^ dir);
				// kv_push(uint64_t, reads[rid].next[difflen], ty);
			}

		}

		if (max_kmer - kmer >= 8 || kmer == min_kmer) {
			for (uint32_t i = 0; i < max_dif_thr; ++i) {
				if (reads[rid].next[i].n > 1) {
					sort(reads[rid].next[i].a, reads[rid].next[i].a + reads[rid].next[i].n, ridcmp);
					uint32_t s = 0;
					for (uint32_t j = 1; j < reads[rid].next[i].n; ++j) {
						if ((reads[rid].next[i].a[s]>>32) != (reads[rid].next[i].a[j]>>32)) {
							reads[rid].next[i].a[++s] = reads[rid].next[i].a[j];
						}
					}
					reads[rid].next[i].n = s + 1;
				}
			}
		}

	}
}

/*void find192NextFun() {
	uint64_t y, ty;
	uint32_t rid, trid, tval;
	int pos, tpos;
	int dir, tdir;
	int16_t difflen, fshift;

	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));

	while (1) {
		uint32_t idx = __sync_fetch_and_add(&rid_pthread, 1);
		if (idx >= ridvecsize) break;

		y = min192vec[idx].y;
		rid = y >> 32;
		pos = (uint32_t)y >> 1;
		dir = y & 1;
		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra);
		}

		for (uint32_t i = idx + 1; i < ridvecsize; ++i) {
			if (min192vec[idx].x != min192vec[i].x) break;

			ty = min192vec[i].y;
			tpos = (uint32_t) ty >> 1;

			if (tpos - pos >= max_dif_thr) break;

			tdir = ty & 1;
			trid = ty >> 32;
			strcpy(strb, seq[trid].seq);
			if (tdir) {
				reverseComplement(strb);
			}

			if (strcmp(stra, strb) == 0) {
				isnextrnd[rid] = false;

				kv_push(uint32_t, reads[trid].crid, rid);
				reads[rid].prid = trid;
				reads[rid].isrc = tdir ^ dir;
				reads[rid].shift = tpos - pos;
				if (tdir) {
					reads[rid].shift = pos - tpos;
				}

				if (reads[rid].next.n > 0) {
					kv_destroy(reads[rid].next);
				}

				break;
			} else {
				//
				difflen = diffstrlen(strb, stra, tpos - pos);
				if (difflen < max_dif_thr) {
					fshift = tpos - pos;
					
					tval = ((uint32_t)fshift<<2) | ((tdir ^ dir) << 1) | tdir;
					ty = (uint64_t)trid<<32 | tval << 16 | (uint32_t) difflen;					
					kv_push(uint64_t, reads[rid].next, ty);
					// if (tdir) {
					// 	fshift = pos - tpos;
					// }
					// ty = (uint64_t)trid<<32 | (uint32_t)fshift<<1 | (tdir ^ dir);
					// kv_push(uint64_t, reads[rid].next[difflen], ty);
				}
			}
		}
	}
}*/

void find192NextFun() {
	uint64_t y, ty;
	uint32_t rid, trid, tval, idx, p;
	int pos, tpos;
	int dir, tdir;
	int16_t difflen, fshift;

	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));

	while (1) {
		p = __sync_fetch_and_add(&rid_pthread, 1);
		// cout << "p: " << p << endl;
		// cout << "isnextrnd[p]: " << isnextrnd[p] << endl;
		if (p >= pend) break;

		if (!isnextrnd[p]) continue;

		if (kmer == max_kmer) {
			// cout << "new reads[" << p << "]\n";
			reads[p].next = new uint64_v[max_dif_thr];
			for (int i = 0; i < max_dif_thr; ++i) {
				kv_init(reads[p].next[i]);
			}
		}

		idx = minidx[p];
		if (idx >= ridvecsize) continue;

		y = min192vec[idx].y;
		rid = y >> 32;
		pos = (uint32_t)y >> 1;
		dir = y & 1;
		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra);
		}

		for (uint32_t i = idx + 1; i < ridvecsize; ++i) {
			if (min192vec[idx].x != min192vec[i].x) break;

			ty = min192vec[i].y;
			tpos = (uint32_t) ty >> 1;

			if (tpos - pos >= max_dif_thr) break;

			if (tpos - pos <= 30 && i - idx > (1U<<16)) break;

			tdir = ty & 1;
			trid = ty >> 32;
			strcpy(strb, seq[trid].seq);
			if (tdir) {
				reverseComplement(strb);
			}

			//
			difflen = diffstrlen(strb, stra, tpos - pos);
			if (difflen < max_dif_thr) {
				fshift = tpos - pos;
				
				tval = ((uint32_t)fshift<<2) | ((tdir ^ dir) << 1) | tdir;
				ty = (uint64_t)trid<<32 | tval;					

				kv_push(uint64_t, reads[rid].next[difflen], ty);
			}
		}

		if (max_kmer - kmer >= 8 || kmer == min_kmer) {
			for (uint32_t i = 0; i < max_dif_thr; ++i) {
				if (reads[rid].next[i].n > 1) {
					sort(reads[rid].next[i].a, reads[rid].next[i].a + reads[rid].next[i].n, ridcmp);
					uint32_t s = 0;
					for (uint32_t j = 1; j < reads[rid].next[i].n; ++j) {
						if ((reads[rid].next[i].a[s]>>32) != (reads[rid].next[i].a[j]>>32)) {
							reads[rid].next[i].a[++s] = reads[rid].next[i].a[j];
						}
					}
					reads[rid].next[i].n = s + 1;
				}
			}
		}
	}
}

// inline 
void collectNext() {
	std::vector<thread> threadVec;
	uint32_t avgnum = 1U<<15, ttid;

	avgnum *= 24;
	// avgnum = 1;

	if (min_kmer <= 31) {
		min128vec = new mm128_t[max_rid];
	}
	if (max_kmer >= 32) {
		min192vec = new mm192_t[max_rid];
	}

	nnsfn = folder + "index.bin";
	FILE *fp = fopen(nnsfn.c_str(), "wb");

	CStopWatch tstopwatch;
	tstopwatch.start();

	for (uint32_t tid = 0; tid < max_rid; tid += avgnum) {
		// nthreads = 1;
		ttid = tid;
		// cout << "ttid: " << ttid << endl;
		for (int kid = 0; kid < kmervecsize; ++kid) {
			kmer = kmervec[kid];
			// cout << "index.bin: " << kmer << endl;

			rid_pthread = ttid;
			pend = rid_pthread + avgnum;
			if (pend > max_rid) pend = max_rid;

			string idxfn = folder + to_string(kmer);
			// cout << "[rid_pthread, pend): " << rid_pthread << ": " << pend << endl;
			idxLoad(idxfn);
			// cout << "after idxLoad()...\n";

			// cout << "[rid_pthread, pend): " << rid_pthread << ": " << pend << endl;
			if (kmer <= 31) {
				// cout << "find128NextFun()\n";
				// cout << "[rid_pthread, pend): " << rid_pthread << ": " << pend << endl;
				for (int i = 0; i < nthreads; ++i) {
					threadVec.push_back(std::thread(find128NextFun));
				}
			} else {
				// cout << "rrrr find192NextFun()\n";
				for (int i = 0; i < nthreads; ++i) {
					threadVec.push_back(std::thread(find192NextFun));
				}
			}
			std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
				thr.join();
			});
			threadVec.clear();

			// cout << "Time of find1XXNextFun() = " << tstopwatch.stop() << std::endl;
			// tstopwatch.resume();
			// exit (0);
		}

		rid_pthread = ttid;
		pend = rid_pthread + avgnum;
		if (pend > max_rid) pend = max_rid;

		// cout << "after...\n";
		// exit(0);
		// write index
		uint32_t tn;
		for (uint32_t rid = rid_pthread; rid < pend; ++ rid) {
			if (isnextrnd[rid]) {
				// cout << "rid: " << rid << endl;
				tn = 0;
				for (int i = 0; i < max_dif_thr; ++i) {
					tn += reads[rid].next[i].n;
				}
				// cout << "tn: " << tn << endl;
				fwrite(&tn, sizeof(uint32_t), 1, fp);
				for (int i = 0; i < max_dif_thr; ++i) {
					if (reads[rid].next[i].n > 0) {
						fwrite(reads[rid].next[i].a, sizeof(uint32_t), reads[rid].next[i].n, fp);
						kv_destroy(reads[rid].next[i]);
					}
				}
				delete[] reads[rid].next;
			}
		}

		// cout << "Time of [ttid, pend) = " << tstopwatch.stop() << std::endl;
		cout << "Time of " << ttid << ", " << pend << ") = " << tstopwatch.stop() << std::endl;
		tstopwatch.resume();
		// cout << "write over..\n";
	}

	if (NULL != min128vec) {
		delete[] min128vec;
		min128vec = NULL;
	}
	if (NULL != min192vec) {
		delete[] min192vec;
		min192vec = NULL;
	}

	fclose(fp);
}

#ifdef false
void processLargeBucketsFun(StrMap *smap, uint64_v *p) {
	uint64_t y, ty;
	uint32_t prid, rid, trid;
	int ppos, pos, tpos;
	int pdir, dir, tdir;
	int16_t shift;
	int16_t mindiff, oridiff, difflen;

	// char *stra = (char*)alloca((L + 1) * sizeof(char));
	// char *strb = (char*)alloca((L + 1) * sizeof(char));
	char *stra = new char[L+1];
	char *strb = new char[L+1];

	string idxstr, curstr;
	int cstartp;

	// int maxnumsearched = 1000;
	int maxnumsearched = 300;

	while (1) {
		uint32_t k = __sync_fetch_and_add(&rid_pthread, 1);
		if (k >= p->n) break;

		// k
		y = p->a[k];
		rid = y >> 32;
		pos = (uint32_t)y>>1;
		dir = y & 1;

		// rid = p->a[k] >> 32;
		// pos = (uint32_t)p->a[k]>>1;
		// dir = p->a[k] & 1;

		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra);
		}
		curstr = (string)(stra);

		mindiff = min_dif2[rid];
		
		int numsearched = 0;
		for (int q = k - 1; q >= 0; --q) {
			tpos = (uint32_t) p->a[q] >> 1;
			tdir = p->a[q] & 1;

			if (k - q >= maxnumsearched && tpos - pos >= 25) { // tpos > pos ?
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
			if (numsearched >= maxnumsearched) break;
		}

		// cout << "111\n";

		// if (false && mindiff != 0 && k > 100) {
		if (mindiff != 0 && k > maxnumsearched) {
			for (int shift = 0; shift <= 35 && mindiff != 0; ++shift) {
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
		// cout << "222\n";

		if (mindiff < max_dif_thr && mindiff < min_dif2[rid]) {
			prid2[rid] = prid;
			isrc2[rid] = pdir ^ dir;
			shift2[rid] = ppos - pos;
			min_dif2[rid] = mindiff;

			if (pdir) {
				shift2[rid] = pos - ppos;
			}

			if (mindiff == 0) {
				isnextrnd[rid] = false;

				rmtx[prid & mask].lock();
				kv_push(uint32_t, reads[prid].crid, rid);
				rmtx[prid & mask].unlock();

				reads[rid].prid = prid;
				reads[rid].isrc = pdir ^ dir;
				reads[rid].shift = ppos - pos;
				// reads[rid].shift = (pdir == 1)? (pos - ppos):(ppos - pos);

				if (pdir) {
					reads[rid].shift = pos - ppos;
				}
			}
		}

	}

	delete[] stra;
	delete[] strb;
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

			mindiff = min_dif2[rid];
			
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
									// if (co[trid] != k && abs(tpos - pos) <= mindiff && shift == abs(tpos - pos)) {
									if (abs(tpos - pos) <= mindiff && shift == abs(tpos - pos)) {

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

			if (mindiff < max_dif_thr && mindiff < min_dif2[rid]) {
				
				prid2[rid] = prid;
				isrc2[rid] = pdir ^ dir;
				shift2[rid] = ppos - pos;
				min_dif2[rid] = mindiff;

				if (pdir) {
					shift2[rid] = pos - ppos;
				}

				if (mindiff == 0) {
					isnextrnd[rid] = false;

					kv_push(uint32_t, reads[prid].crid, rid);
					reads[rid].prid = prid;
					reads[rid].isrc = pdir ^ dir;
					reads[rid].shift = ppos - pos;
					// reads[rid].shift = (pdir == 1)? (pos - ppos):(ppos - pos);

					if (pdir) {
						reads[rid].shift = pos - ppos;
					}
				}
			}

		idxstr = curstr.substr(startp[0], sublen[0]);
		addStrMap(&smap[0], idxstr, y);
		idxstr = curstr.substr(startp[1], sublen[1]);
		addStrMap(&smap[1], idxstr, y);
	}
}

void processLargeBuckets() {
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

		// cout << "i: " << i << endl;
		// cout << "n: " << n << endl;

		if (n <= 50000) {

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

			// if (kmer == 26)	nthreads = 1;
			rid_pthread = 1;
			for (int i = 0; i < nthreads; ++i) {
				threadVec.push_back(std::thread(processLargeBucketsFun, smap, &p));
			}
			std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
				thr.join();
			});
			threadVec.clear();

			kv_destroy(p);
			delStrMap(&smap[0]);
			delStrMap(&smap[1]);
			// smap[0].clear();
			// smap[1].clear();
			delete[] smap;
		} else {
			/*rid_pthread = rand()% (n - 20001);
			p.n = rid_pthread + 20000;

			StrMap *smap = new StrMap[2];

			kv_init(p);
			kv_resize(uint64_t, p, n);
			for (uint32_t k = 0; k < rid_pthread; ++k) {
				// creat hash table here?
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
			// cout << "before processLargeBucketsFun() ...\n";

			// nthreads = 1;

			for (int i = 0; i < nthreads; ++i) {
				threadVec.push_back(std::thread(processLargeBucketsRandomFun, smap, &p));
			}
			std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
				thr.join();
			});
			threadVec.clear();
			kv_destroy(p);
			smap[0].clear();
			smap[1].clear();
			delete[] smap;
			*/
		}

		// cout << "after ...\n";
	}
	kv_destroy(lrec);

	cout << "Time of processLargeBuckets() = " << stopwatch.stop() << std::endl;
	stopwatch.resume();
}
#endif