#include "mstcom.h"

void collectNextFunOri() {
	uint64_t y, ty;
	uint32_t rid, trid, prid;
	int pos, tpos, ppos;
	int dir, tdir, pdir;
	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));
	int16_t mindiff, oridiff, difflen, shift;
	int searchedNum;

	while (1) {
		uint32_t idx = __sync_fetch_and_add(&rid_pthread, 1);
		if (idx >= minisize) break;

		y = mini[idx];
		rid = y >> 32;

		// cout << idx << "; " << rid << "; " << seq[rid].seq << "; main\n";

		if (isnextrnd[rid]) continue; // not for nexr round
		if (reads[prid2[rid]].root != reads[rid].root) continue; //prid2[rid] and rid are not in the same tree

		if (kmer == max_kmer) {
			min_dif2[rid] = max_dif_thr;
		}
		mindiff = min_dif2[rid];

		pos = (uint32_t)y >> 2;
		dir = (y >> 1) & 1;
		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra);
		}
		// cout << stra << endl;

		searchedNum = 0;
		for (uint32_t i = idx + 1; i < minisize; ++i) {
			// cout << (mini[idx]&1) << "; " << (mini[i]&1) << endl;
			if ((mini[idx]&1) != (mini[i]&1)) break;

			// if (i - idx >= 400) break;
			ty = mini[i];
			trid = ty >> 32;

			if (reads[trid].root == reads[rid].root) continue; //already in the same tree

			tpos = (uint32_t)ty >> 2;
			tdir = (ty>>1) & 1;

			if (abs(tpos - pos) >= mindiff) break;

			strcpy(strb, seq[trid].seq);
			if (tdir) {
				reverseComplement(strb);
			}
			// cout << i << "; " << trid << "; " << strb << "\n";

			difflen = diffstrlen(strb, stra, tpos - pos);

			// cout << "difflen: " << difflen << endl;

			if (difflen < mindiff) {
				prid = trid;
				ppos = tpos;
				pdir = tdir;

				mindiff = difflen;
			}
			++ searchedNum;

			// if (searchedNum > 400) break;
			if (searchedNum > 1000) break;

		}
		if (mindiff < max_dif_thr && mindiff < min_dif2[rid]) {
			prid2[rid] = prid;
			isrc2[rid] = pdir ^ dir;
			shift2[rid] = (pdir == 1) ? (pos - ppos):(ppos - pos);
			min_dif2[rid] = mindiff; 
		} 
		// exit(0);
		// fprintf(stderr, "%lu\n", i);
	}
}

void collectNextFun() {
	uint64_t y, ty;
	uint32_t rid, trid, prid;
	int pos, tpos, ppos;
	int dir, tdir, pdir;
	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));
	int16_t mindiff, oridiff, difflen, shift;
	int searchedNum, maxSearchedNum; 

	while (1) {
		uint32_t idx = __sync_fetch_and_add(&rid_pthread, 1);
		if (idx >= minisize) break;

		bool debug = false;

		maxSearchedNum = 400;
		// maxSearchedNum = 2000;
		y = mini[idx];
		rid = y >> 32;

		if (debug) cout << idx << "; " << rid << "; " << seq[rid].seq << "; main\n";

		if (isnextrnd[rid]) continue; // not for next round
		// if (reads[prid2[rid]].root != reads[rid].root) continue; //prid2[rid] and rid are not in the same tree

		if (reads[prid2[rid]].root == reads[rid].root) {
			min_dif2[rid] = max_dif_thr;
		}
		mindiff = min_dif2[rid];

		pos = (uint32_t)y >> 2;
		dir = (y >> 1) & 1;
		strcpy(stra, seq[rid].seq);
		if (dir) {
			reverseComplement(stra);
		}
		// cout << stra << endl;

		searchedNum = 0;
		for (uint32_t i = idx - 1; i >= 0; --i) {
			// cout << (mini[idx]&1) << "; " << (mini[i]&1) << endl;
			if ((mini[idx]&1) != (mini[i]&1)) break;

			// if (i - idx >= 400) break;
			ty = mini[i];
			trid = ty >> 32;

			if (reads[trid].root == reads[rid].root) {
				if (i == 0) break;
				continue; //already in the same tree
			}

			tpos = (uint32_t)ty >> 2;
			tdir = (ty >> 1) & 1;

			if (abs(tpos - pos) > mindiff) break;

			strcpy(strb, seq[trid].seq);
			if (tdir) {
				reverseComplement(strb);
			}
			if (debug) cout << i << "; " << trid << "; " << strb << "\n";

			// difflen = diffstrlen(strb, stra, tpos - pos);
			difflen = (diffstrlen(strb, stra, tpos - pos) + diffstrlen(stra, strb, pos - tpos))/2;

			// cout << "difflen: " << difflen << endl;

			if (difflen < mindiff) {
				prid = trid;
				ppos = tpos;
				pdir = tdir;

				mindiff = difflen;
			}
			++ searchedNum;

			if (searchedNum > maxSearchedNum) break;
			// if (searchedNum > 800) break;
			// if (searchedNum > 1000) break;
			if (i == 0) break;

		}
		if (mindiff < max_dif_thr && mindiff < min_dif2[rid]) {
			prid2[rid] = prid;
			isrc2[rid] = pdir ^ dir;
			shift2[rid] = (pdir == 1) ? (pos - ppos):(ppos - pos);
			min_dif2[rid] = mindiff; 
		} 
	}
}

// inline 
void collectNext() {
	/*FILE *fp = fopen("seq.txt", "w");
	char *str = (char*)alloca((L + 1) * sizeof(char));
	uint64_t y;
	for (uint32_t i = 0; i < minisize; ++i) {
		y = mini[i];
		strcpy(str, seq[y >> 32].seq);
		if ((y>>1) & 1) {
			reverseComplement(str);
		}
		fprintf(fp, "%s %u %u %d\n", str, y>>32, (uint32_t)y >> 2, y & 1);
	}
	fclose(fp);
	exit(0);*/

	CStopWatch tstopwatch;
	tstopwatch.start();
	tstopwatch.resume();

	rid_pthread = 1;
	// nthreads = 1;
	std::vector<thread> threadVec;
	// fprintf(stderr, "rec.n: %lu\n", rec.n);
	for (int i = 0; i < nthreads; ++i) {
		threadVec.push_back(std::thread(collectNextFun));
	}
	// cout << "111xxxx22222\n";
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	// cout << "111xxxx222224444\n";
	threadVec.clear();
	// fclose(fptmp);

	cout << "Time of collectNext() = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();
}
