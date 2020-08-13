# include "mstcom.h"
using namespace std;

int nthreads;
int max_dif_thr;
int L;
int nk;
size_t RN;
uint32_t max_rid;
bseq1_t *seq;
FILE *fppar;

// input parameters
bool isorder;
bool ispe;
// string outputfn;

// --- for minimizer
int bsize = 20;
mm128_v *B;
uint64_v *Bs;
mm192_v *BL;
// --- for minimizer

int kmer, max_kmer, min_kmer, gkmer;
int *kmervec, kmervecsize;

Edge_v edges;
mm128_t *min128vec;
mm192_t *min192vec;
uint64_t *mini;
uint32_t minisize;
uint32_t *snum;

mutex *bmtx;

std::string folder;
std::string hammingedgfolder;
std::string shiftedgfolder;

READS_t *reads;
bool *isnextrnd;

// uint32_t maxedges = 1<<20;
uint32_t maxedges;
uint32_t *rootarr;

char complement[256];

CStopWatch stopwatch;

uint32_t encode(uint32_t rid, char *en_str);
uint32_t getRoot(uint32_t rid);
inline void reRootNode(const size_t &_rid, const size_t &root);
// inline void reRootNode(const size_t &_rid);

inline void init() {
	memset(complement, 125, 256);
	complement['A'] = 'T';
	complement['a'] = 'T';
	complement['C'] = 'G';
	complement['c'] = 'G';
	complement['G'] = 'C';
	complement['g'] = 'C';
	complement['T'] = 'A';
	complement['t'] = 'A';
	complement['N'] = 'N';
	complement['n'] = 'N';
}

string infile, infile1, outfile;

uint32_t rid_pthread;
int mask = 1<<25 - 1;

bool debug = false;

void iscircle() {
	bool *v = new bool[max_rid], flag = true;
	memset(v, 0, sizeof(bool) * max_rid);
	for (uint32_t rid = 0; flag && rid < max_rid; ++rid) {
		if (reads[rid].prid == rid) {
			uint32_t noderid = rid;
			queue<uint32_t> q;
			q.push(noderid);
			while (!q.empty()) {
				noderid = q.front();
				q.pop();
				if (v[noderid]) {
					flag = false;
					break;
				}
				v[noderid] = true;
				for (uint32_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}
		}
	}
	delete[] v;
	if (flag) {
		cout << "No Circle\n";
	} else {
		cout << "Is Circle\n";
	}
}

inline void labelRoot() {
	size_t num_of_nodes = 0, num_of_root = 0;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (!isnextrnd[rid] && reads[rid].prid == rid) { //a root node 
			++ num_of_root;
			uint32_t noderid = rid;
			queue<uint32_t> q;
			q.push(noderid);
			while (!q.empty()) {
				noderid = q.front();
				rootarr[noderid] = rid;

				if (noderid == 1539203) {
					// cout << "root rid: " << rid << endl;
				}

				q.pop();
				++ num_of_nodes;
				for (uint32_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}
		}
	}
	cout << "num_of_root: " << num_of_root << endl;
	cout << "num_of_nodes: " << num_of_nodes << endl;
	// cout << "Time of labelRoot() = " << stopwatch.stop() << std::endl;
	// stopwatch.resume();
}

mmrec_v rec, lrec;
mutex recmtx;

int tree_num;

void countTreeFun() {
	while (1) {
		uint32_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;

		if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) {
			__sync_fetch_and_add(&tree_num, 1);
		}
	}
}

inline void countTree() {
	tree_num = 0;
	// rid_pthread = 0;
	// std::vector<thread> threadVec;
	// for (int i = 0; i < nthreads; ++i) {
	// 	threadVec.push_back(std::thread(countTreeFun));
	// }
	// std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
	// 	thr.join();
	// });
	// threadVec.clear();
	cout << "begin countTree()..\n";
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (!isnextrnd[rid]) {
			if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) {
				++ tree_num;
			}
		}
	}
	fprintf(stderr, "trees number: %d\n", tree_num);
}

uint32_t coutIsNext() {
	uint32_t cnt = 0;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (!isnextrnd[rid]) {
			++ cnt;
		}
	}
	return cnt;
	// cout << "isnextrnd: " << cnt << endl;
}

void mstConstruction(Edge_v &edges);

void reRootTree() {

}

threadShiftEdges_t *threadshifteds;

void indexConstruction(string kinds) {
	CStopWatch tstopwatch;
	tstopwatch.start();
	stopwatch.start();

	if (kinds == "minimizer")
		calcMinimizers();
	else 
	if (kinds == "maximizer") 
		calcMaximizers();
	cout << "Time of calc" + kinds + "() = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();

	sortBuckets();
	cout << "Time of sortBuckets() = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();

	processBuckets();
	obtainMiniIdx();

	if (kmer <= 31) {
		for (int i = 0; i < (1 << bsize); ++i) {
			kv_destroy(B[i]);
			kv_init(B[i]);
		}
	} else
	if (kmer > 31) {
		for (int i = 0; i < (1 << bsize); ++i) {
			kv_destroy(BL[i]);
			kv_init(BL[i]);
		}
	}

	cout << "Time of indexConstruction() = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();
}

void edgesConstruction() {
	kv_init(edges); // store some duplicate can not be detected by minimizer

	// hammingedgfolder = folder + "edges";
	mode_t mode = 0777;
	// int nError = mkdir(hammingedgfolder.c_str(), mode);
	// if (nError != 0) {
	// 	fprintf(stderr, "Failed to creat the folder '%s'\n", hammingedgfolder.c_str());
	// 	exit(EXIT_FAILURE);
	// }
	// hammingedgfolder += "/";
	//
	shiftedgfolder = folder + "shiftedges";
	int nError = mkdir(shiftedgfolder.c_str(), mode);
	if (nError != 0) {
		fprintf(stderr, "Failed to creat the folder '%s'\n", shiftedgfolder.c_str());
		exit(EXIT_FAILURE);
	}
	shiftedgfolder += "/";

	maxedges = 1<<20;
	// threadeds = new threadEdges_t[nthreads];
	threadshifteds = new threadShiftEdges_t[nthreads];

	for (int threadid = 0; threadid < nthreads; ++threadid) {
		// threadeds[threadid].init(threadid);
		threadshifteds[threadid].init(threadid);
	}
	CStopWatch tstopwatch;
	tstopwatch.start();
	// cout << "maxedges: " << maxedges << endl;

	min128vec = NULL;
	min192vec = NULL;
	if (min_kmer <= 31) {
		B = (mm128_v*)calloc(1 << bsize, sizeof(mm128_v));
		min128vec = new mm128_t[max_rid];
		for (int i = 0; i < (1 << bsize); ++i) {
			kv_init(B[i]);
		}
	} 
	if (max_kmer > 31) {
		BL = (mm192_v*)calloc(1 << bsize, sizeof(mm192_v));
		min192vec = new mm192_t[max_rid];
		for (int i = 0; i < (1 << bsize); ++i) {
			kv_init(BL[i]);
		}
	}

	for (int kid = 0; kid < kmervecsize; ++kid) {
		kmer = kmervec[kid];
		// cout << "kmer: " << kmer << endl;
		tstopwatch.resume();
		
		indexConstruction("minimizer");
		collectEdges("min");

		// if (kid >= 10) {
			indexConstruction("maximizer");
			collectEdges("max");
		// }

		cout << "All time of kmer = " << kmer << ": " << tstopwatch.stop() << std::endl;
		tstopwatch.resume();
	}

	if (min_kmer <= 31) {
		free(B);
	}
	if (max_kmer > 31) {
		free(BL);
	}
	if (NULL != min128vec) {
		delete[] min128vec;
	}
	if (NULL != min192vec) {
		delete[] min192vec;
	}

	cout << "after edgesConstruction().." << endl;
}

void removeChild(size_t prid, size_t rid);

typedef std::pair<uint32_t, uint32_t> DuplicateEdgePair;
std::vector<DuplicateEdgePair> duplicateEdges;

void removeDuplicateEdges_() {
	// cout << "duplicateEdges.size(): " << duplicateEdges.size() << endl;

	for (size_t k = 0; k < duplicateEdges.size(); ++k) {
		uint32_t rid = duplicateEdges[k].first;
		uint32_t prid = duplicateEdges[k].second;
		bool isfind = false;
		for (uint32_t i = 0; i < reads[prid].crid.n; ++i) {
			if (reads[prid].crid.a[i] == rid) { // find the edge
				isfind = true;
				break;
			}
		}

		if (!isfind) {
			rid = duplicateEdges[k].second;
			prid = duplicateEdges[k].first;
			for (uint32_t i = 0; i < reads[prid].crid.n; ++i) {
				if (reads[prid].crid.a[i] == rid) { // find the edge
					isfind = true;
					break;
				}
			}
		}

		if (isfind) {
			// cout << "prid: " << prid << ", rid: " << rid << endl;
			removeChild(prid, rid);

			for (uint32_t i = 0; i < reads[rid].crid.n; ++i) {
				uint32_t crid = reads[rid].crid.a[i];
				kv_push(uint32_t, reads[prid].crid, crid);
				reads[crid].prid = prid;

				if (reads[rid].isrc) {
					reads[crid].isrc ^= 1;
					reads[crid].shift = 0 - reads[crid].shift;
				}
			}

			READS_t *tr = &reads[prid];

			dup_t *tempdupvec = new dup_t[tr->dn + 2];
			if (isorder) {
				uint32_t dn = tr->dn;
				for (uint32_t i = 1; i < dn + 1; ++i) {
					tempdupvec[i] = tr->dup[i];
				}
				reads[prid].dn ++;
				tempdupvec[dn + 1] = dup_t(rid, reads[rid].isrc);
			} else {
				uint32_t dn = tr->dn;
				for (uint32_t i = 0; i < dn; ++i) {
					tempdupvec[i] = tr->dup[i];
				}
				reads[prid].dn ++;
				tempdupvec[dn ++] = dup_t(rid, reads[rid].isrc);
			}

			delete[] tr->dup;
			// tr->dup = new dup_t[dn];
			tr->dup = tempdupvec;

			// isnextrnd[rid] = true;
			reads[rid].prid = prid;
			kv_destroy(reads[rid].crid);
			reads[rid].crid.n = 0;
		}
	}
}

void removeDuplicateEdges() {
	// cout << "duplicateEdges.size(): " << duplicateEdges.size() << endl;

	for (size_t k = 0; k < duplicateEdges.size(); ++k) {
		uint32_t rid = duplicateEdges[k].first;
		uint32_t prid = duplicateEdges[k].second;
		bool isfind = false;
		for (uint32_t i = 0; i < reads[prid].crid.n; ++i) {
			if (reads[prid].crid.a[i] == rid) { // find the edge
				isfind = true;
				break;
			}
		}

		if (!isfind) {
			rid = duplicateEdges[k].second;
			prid = duplicateEdges[k].first;
			for (uint32_t i = 0; i < reads[prid].crid.n; ++i) {
				if (reads[prid].crid.a[i] == rid) { // find the edge
					isfind = true;
					break;
				}
			}
		}

		if (isfind) {
			// cout << "prid: " << prid << ", rid: " << rid << endl;
			removeChild(prid, rid);

			for (uint32_t i = 0; i < reads[rid].crid.n; ++i) {
				uint32_t crid = reads[rid].crid.a[i];
				kv_push(uint32_t, reads[prid].crid, crid);
				reads[crid].prid = prid;

				if (reads[rid].isrc) {
					reads[crid].isrc ^= 1;
					reads[crid].shift = 0 - reads[crid].shift;
				}
			}

			READS_t *tr = &reads[prid];

			dup_t *tempdupvec = new dup_t[tr->dn + 2];
			if (isorder) {
				uint32_t dn = tr->dn;
				for (uint32_t i = 1; i < dn + 1; ++i) {
					tempdupvec[i] = tr->dup[i];
				}
				reads[prid].dn ++;
				tempdupvec[dn + 1] = dup_t(rid, reads[rid].isrc);
			} else {
				uint32_t dn = tr->dn;
				for (uint32_t i = 0; i < dn; ++i) {
					tempdupvec[i] = tr->dup[i];
				}
				reads[prid].dn ++;
				tempdupvec[dn ++] = dup_t(rid, reads[rid].isrc);
			}

			delete[] tr->dup;
			// tr->dup = new dup_t[dn];
			tr->dup = tempdupvec;

			// isnextrnd[rid] = true;
			reads[rid].prid = prid;
			kv_destroy(reads[rid].crid);
			kv_init(reads[rid].crid);
			reads[rid].crid.n = 0;
		}
	}
}

void mstConstruction() {
	CStopWatch tstopwatch;
	tstopwatch.start();

	uint32_t rid, prid, aroot, broot;
	size_t edgesidx = 0;

	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *encstr = (char*)alloca(((L<<2) + 1) * sizeof(char));

	rootarr = new uint32_t[max_rid];

	for (uint32_t i = 0; i < max_rid; ++i) {
		rootarr[i] = i;
	}
	cout << "begin mstConstruction() ..." << endl;

	MstEdge_v mstedges;
	kv_init(mstedges);

	// process dif == 0;
	while (edgesidx < edges.n && edges.a[edgesidx].dif == 0) {
		rid = edges.a[edgesidx].rid;
		prid = edges.a[edgesidx].prid;
		if (getRoot(rid) != getRoot(prid)) { //these two reads not in a tree
			aroot = getRoot(rid);
			broot = getRoot(prid);
			rootarr[aroot] = broot;

			MstEdge_t teds(rid, prid, 0, edges.a[edgesidx].isrc);
			kv_push(MstEdge_t, mstedges, teds);

			duplicateEdges.push_back(make_pair(rid, prid));
		}
		++ edgesidx;
	}
	kv_destroy(edges);

	for (int dif = 1; dif <= max_dif_thr; ++ dif) {
		// if (false) // for debug
		// for shift edges
		for (int threadid = 0; threadid < nthreads; ++threadid) {
			shiftEdgesFile_t *edf = &threadshifteds[threadid].a[dif - 1];
			// edf->reopen(maxedges);
			edf->reopen();

			while (edf->get()) {
			// bool getres;
			// while (getres = edf->get()) {
				// cout << "edf->fn: " << edf->fn << endl;
				// cout << "getres: " << getres << endl;
				// cout << "edf->n: " << edf->n << endl;
				// cout << "edf->m: " << edf->m << endl;

				for (uint32_t i = 0; i < edf->n; ++i) {
					uint64_t v = edf->a[i];
					rid = v >> 32;
					prid = (uint32_t) v;
					int16_t shiftisrc = edf->shift[i];

					if (getRoot(rid) != getRoot(prid)) { //these two reads not in a tree
						MstEdge_t teds(rid, prid, shiftisrc >> 1, shiftisrc & 1);
						kv_push(MstEdge_t, mstedges, teds);
						//
						aroot = getRoot(rid);
						broot = getRoot(prid);
						rootarr[aroot] = broot;
					}
				}
			}
			// cout << "before clear\n";
			edf->clear();
		}
	}
	cout << endl;
	for (int threadid = 0; threadid < nthreads; ++threadid) {
		threadshifteds[threadid].removefd();
		// threadeds[threadid].removefd();
	}
	// edgfolder = folder + "edges";
	// string cmd = "rm -rf " + folder + "edges " + folder + "shiftedges " + folder + "shiftedgesbysubstridx";
	// string cmd = "rm -rf " + hammingedgfolder + " " + shiftedgfolder ;
	string cmd = "rm -rf " + shiftedgfolder ;
	system(cmd.c_str());

	cout << "Time of edges of mstConstruction() = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();

	// radix_sort_mstedge(mstedges.a, mstedges.a + mstedges.n);

	bool *v = new bool[max_rid], isrc;
	memset(v, false, sizeof(bool) * max_rid);
	int16_t shift;
	// FILE *fp = fopen("1.txt", "w");
	cout << "Number of edges in MST: " << mstedges.n << endl;
	for (size_t i = 0; i < mstedges.n; ++i) {
		rid = mstedges.a[i].rid;
		prid = mstedges.a[i].prid;
		isrc = mstedges.a[i].isrc;
		shift = mstedges.a[i].shift;
		// fprintf(fp, "%u %u\n", prid, rid);

		if (v[rid] && !v[prid]) {
			rid = prid;
			prid = mstedges.a[i].rid;
			if (!isrc) {
				shift = 0 - shift; 
			}
		} 
		else 
		if (v[rid] && v[prid]) {
			// cout << "prid: " << prid << "; " << "rid: " << rid << endl;
			reRootNode(rid);
		}
		reads[rid].prid = prid;
		reads[rid].shift = shift;
		reads[rid].isrc = isrc;
		kv_push(uint32_t, reads[prid].crid, rid);
		v[rid] = v[prid] = true;
	}
	delete[] v;
	delete[] rootarr;
	cout << "Time of mstConstruction() = " << tstopwatch.stop() << std::endl;
	// fclose(fp);

	removeDuplicateEdges();

	cout << "end mstConstruction() ...\n";
}

bool cmp2(const ROOTNODE_t &a, const ROOTNODE_t &b) {
	return a.nodecnt > b.nodecnt;
}

// inline 
void removeChild(size_t prid, size_t rid) {
	size_t cridsize = reads[prid].crid.n;
	// if (debug) fprintf(stderr, "cridsize: %lu\n", reads[prid].crid.size());
	for (size_t i = 0; i < cridsize; ++i) {
		if (reads[prid].crid.a[i] == rid) {
			reads[prid].crid.a[i] = reads[prid].crid.a[cridsize - 1];
			break;
		}
	} 
	// reads[prid].crid.resize(cridsize - 1);
	-- reads[prid].crid.n; //resize(cridsize - 1);
}
 
inline void reRootNode(const uint32_t &_rid, const uint32_t &root) { // set rid as root node
	if (reads[_rid].prid == _rid) return; // rid is root node; do nothing
	// fprintf(stderr, "reRootNode()...\n");
	// fprintf(stderr, "%lu %lu\n", _rid, reads[_rid].prid);

	bool debug = false;
	// if (_rid == 49755) debug = true;
	if (debug) fprintf(stderr, "------\n");

	uint32_t rid = _rid;
	uint32_t oriprid = reads[rid].prid, tprid;
	bool isrc = reads[rid].isrc, tisrc;
	int16_t shift = reads[rid].shift, tshift;
	while (true) {
		if (debug) fprintf(stderr, "oriprid: %lu; rid: %lu\n", oriprid, rid);
		removeChild(oriprid, rid);
		tprid = reads[oriprid].prid;
		tisrc = reads[oriprid].isrc;
		tshift = reads[oriprid].shift;
		
		reads[oriprid].prid = rid; // reverse the edge
		if (isrc) {
			reads[oriprid].shift = shift;
		} else {
			reads[oriprid].shift = - shift;
		}
		reads[oriprid].isrc = isrc;
		// reads[rid].crid.push_back(oriprid);
		kv_push(uint32_t, reads[rid].crid, oriprid);

		rid = oriprid;
		oriprid = tprid;
		
		if (oriprid == rid) break;

		isrc = tisrc;
		shift = tshift;
	}
	uint32_t noderid;
	reads[_rid].prid = _rid;
	reads[_rid].isrc = 0;
	reads[_rid].shift = 0;
}

bool edgecmp(Edge_t a, Edge_t b) {
	return a.dif < b.dif;
}

void reRootNode(const uint32_t &_rid) { // set rid as root node
	if (reads[_rid].prid == _rid) return; // rid is root node; do nothing

	uint32_t rid = _rid;
	uint32_t oriprid = reads[rid].prid, tprid;
	bool isrc = reads[rid].isrc, tisrc;
	int16_t shift = reads[rid].shift, tshift;
	while (true) {
		removeChild(oriprid, rid);
		tprid = reads[oriprid].prid;
		tisrc = reads[oriprid].isrc;
		tshift = reads[oriprid].shift;
		
		reads[oriprid].prid = rid; // reverse the edge
		if (isrc) {
			reads[oriprid].shift = shift;
		} else {
			reads[oriprid].shift = - shift;
		}
		reads[oriprid].isrc = isrc;
		// reads[rid].crid.push_back(oriprid);
		kv_push(uint32_t, reads[rid].crid, oriprid);

		rid = oriprid;
		oriprid = tprid;
		
		if (oriprid == rid) break;

		isrc = tisrc;
		shift = tshift;
	}
	uint32_t noderid;
	reads[_rid].prid = _rid;
	reads[_rid].isrc = 0;
	reads[_rid].shift = 0;
}

uint32_t getRoot(uint32_t rid) {
	uint32_t root = rid, k, j;
	while (rootarr[root] != root) {
		root = rootarr[root];
	}
	k = rid;
	while (k != root)  {
		j = rootarr[k];
		rootarr[k] = root;
		k = j;
	}
	return root;
}

int mainccc() { // test bsc compress file
	string infile = "mstcom_3K49nHxPuo/11.out"; 
	char outfile[] = "11.out.bsc"; 
	char out[] = "11.out.dec"; 
	mstcom::bsc::BSC_compress(infile.c_str(), outfile, 64);
	mstcom::bsc::BSC_decompress(outfile, out);
}

int mainXX() { // test lzma compress file
	string infile = "idx.txt"; 
	char outfile[] = "idx.txt.lzma0"; 
	char out[] = "idx.txt.dec"; 
	mstcom::lzma::lzma_compress(infile.c_str(), outfile);
	mstcom::lzma::lzma_decompress(outfile, out);
}

// inline 
uint32_t encode_(uint32_t rid, char *en_str) {
	if (reads[rid].prid == rid) return L;
	uint32_t enclen = 0;
	// char *en_str = (char*)alloca((L + 1) * sizeof(char));

	if (reads[rid].isrc) {
		char *rcstr = (char*)alloca((L + 3) * sizeof(char));
		strncpy(rcstr, seq[rid].seq, L);
		rcstr[L] = '\0';
		reverseComplement(rcstr);
		encode(seq[reads[rid].prid].seq, rcstr, reads[rid].shift, en_str);
	} else {
		encode(seq[reads[rid].prid].seq, seq[rid].seq, reads[rid].shift, en_str);
	}
	enclen = strlen(en_str);
	return enclen;
}

uint32_t encode(uint32_t rid, char *en_str) {
	if (reads[rid].prid == rid) return L;
	uint32_t enclen = 0;
	// char *en_str = (char*)alloca((L + 1) * sizeof(char));

	if (reads[rid].isrc) {
		char *rcstr = (char*)alloca((L + 3) * sizeof(char));
		strncpy(rcstr, seq[rid].seq, L);
		rcstr[L] = '\0';
		reverseComplement(rcstr);
		encode(seq[reads[rid].prid].seq, rcstr, reads[rid].shift, en_str);
	} else {
		encode(seq[reads[rid].prid].seq, seq[rid].seq, reads[rid].shift, en_str);
	}
	enclen = strlen(en_str);
	return enclen;
}

int encode_v3(char *parent, char *child, const int16_t &_shift, char *en_str) { //str1 is the short one string
	char *int_str = (char*)alloca(20 * sizeof(char));
	int16_t shift = _shift;
	int mismatchno = 0;
	// bool debug = false;

	// if (debugcnt >= mindebugcnt && debugcnt < mindebugcnt + 10) {
	// 	debug = true;
	// }
	// ++debugcnt;

	// if (strcmp(child, "GGTCTCAAACTGCTGACTTCAAGTGATCTGCCCGCCTTGGCCTCCCAAAGTGCTGAGATTACGGATGTGAGCCACTGTGCCCAAATTTTTTTTTTTTTTTT") == 0) {
		// debug = true;
	// }
	int en_str_len = 0;
	int eq_char_num = 0;

	if (shift >= 0) {
		// case shift >= 0
		//    AATTGCATGC parent
		//      TTGCATGCGA
		if (shift > 0) {
			// for (int i = 1; i <= shift; ++i) {
			// 	en_str[en_str_len++] = child[L - i];
			// }
			for (int i = shift; i >= 1; --i) {
				en_str[en_str_len++] = child[L - i];
			}
		}
		int i, j;
		for (i = shift, j = 0; i < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				++ mismatchno;

				sprintf(int_str, "%d", eq_char_num);
				for (char *tk = int_str; *tk != '\0'; ++tk) {
					en_str[en_str_len++] = *tk;
				}
				eq_char_num = 0;
				en_str[en_str_len++] = child[j]; 
			} else ++ eq_char_num;
		}
		en_str[en_str_len] = '\0';
	} else {
		// cast: shift < 0
		//       AATTGCATGC parent
		//    TCGAATTGCA
		int i, j;
		shift = 0 - shift;
		for (j = 0; j < shift; ++j) {
			en_str[en_str_len++] = child[j];
		}
		for (i = 0, j = shift; j < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				++ mismatchno;

				sprintf(int_str, "%d", eq_char_num);
				for (char *tk = int_str; *tk != '\0'; ++tk) {
					en_str[en_str_len++] = *tk;
				}
				eq_char_num = 0;
				en_str[en_str_len++] = child[j]; 
			} else ++ eq_char_num;	
		}
		en_str[en_str_len] = '\0';
 	}

 	// if (debug) {
 	// 	cout << "parent: " << parent << endl;
 	// 	cout << "child: " << child << endl;
 	// 	cout << "shift: " << + _shift << endl;
 	// 	cout << en_str << endl;
 	// }

 	return mismatchno;
}

int encode_v3(uint32_t rid, char *en_str) {
	if (reads[rid].prid == rid) return L;
	uint32_t enclen = 0;
	int mismatchno;

	if (reads[rid].isrc) {
		char *rcstr = (char*)alloca((L + 3) * sizeof(char));
		strncpy(rcstr, seq[rid].seq, L);
		rcstr[L] = '\0';
		reverseComplement(rcstr);
		mismatchno = encode_v3(seq[reads[rid].prid].seq, rcstr, reads[rid].shift, en_str);
	} else {
		mismatchno = encode_v3(seq[reads[rid].prid].seq, seq[rid].seq, reads[rid].shift, en_str);
	}
	return mismatchno;
}

void encodeReadsFun() {
	char *encstr = (char*)alloca((L + 1) * sizeof(char));
	uint32_t noderid;

	while (1) {
		uint32_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;

		// if (reads[rid].prid == rid) {
		if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) {
			uint32_v v;
			kv_init(v);
			kv_push(uint32_t, v, rid);
			uint32_t i = 0;
			while (i < v.n) {
				noderid = v.a[i];
				for (uint32_t j = 0; j < reads[noderid].crid.n; ++j) {
					kv_push(uint32_t, v, reads[noderid].crid.a[j]);
				}
				++ i;
			}
			for (i = v.n - 1; i >= 0; --i) {
				noderid = v.a[i];
				int encstrlen =	encode_v3(noderid, encstr);
				if (encstrlen < L) {
					free(seq[noderid].seq);
					seq[noderid].seq = strdup(encstr);
				}

				if (isorder)
				// for duplicate reads
				for (uint32_t w = 1; w < reads[noderid].dn + 1; ++w) {
					uint32_t trid = reads[noderid].dup[w].id;
					reads[trid].prid = noderid;		
					reads[trid].isrc = reads[noderid].dup[w].isrc;
					seq[trid].seq[0] = '\0'; 
				}

				if (i == 0) break;
			}
			kv_destroy(v);
		}
	}
}

void encodeReads() {
	cout << "begin encodeReads()" << endl;
	rid_pthread = 0;
	std::vector<thread> threadVec;
	for (int i = 0; i < nthreads; ++i) {
		threadVec.push_back(std::thread(encodeReadsFun));
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
	cout << "after encodeReads()" << endl;
}

void encodeReads_() {
	cout << "begin encodeReads()" << endl;
	char *encstr = (char*)alloca((L + 1) * sizeof(char));
	uint32_t noderid;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		// if (reads[rid].prid == rid) {
		if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) {
			// cout << "rid: " << rid << endl;
			uint32_v v;
			kv_init(v);
			kv_push(uint32_t, v, rid);
			uint32_t i = 0;
			while (i < v.n) {
				noderid = v.a[i];
				for (uint32_t j = 0; j < reads[noderid].crid.n; ++j) {
					kv_push(uint32_t, v, reads[noderid].crid.a[j]);
				}
				++ i;
			}
			// cout << "v.n: " << v.n << endl;
			for (i = v.n - 1; i >= 0; --i) {
				noderid = v.a[i];
				// cout << "noderid: " << noderid << endl;
				int encstrlen =	encode_v3(noderid, encstr);
				// cout << "111" << endl;
				if (encstrlen < L) {
					free(seq[noderid].seq);
					seq[noderid].seq = strdup(encstr);
				}
				// cout << "2222222" << endl;

				if (isorder)
				// for duplicate reads
				for (uint32_t w = 1; w < reads[noderid].dn + 1; ++w) {
					uint32_t trid = reads[noderid].dup[w].id;
					reads[trid].prid = noderid;		
					reads[trid].isrc = reads[noderid].dup[w].isrc;
					seq[trid].seq[0] = '\0'; 
				}

				if (i == 0) break;
			}
			kv_destroy(v);
		}
	}
	cout << "after encodeReads()" << endl;
}

const char invert_code_rule[5] = {'A', 'C', 'G', 'T', 'N'};

void storeTrees(string fn);
void loadTrees(string fn);

int compress_main(int argc, char *argv[]) {
	
	stopwatch.start();

	max_dif_thr = 0;
	init();
	getPars(argc, argv);

	L = getReadsLength(infile.c_str());

	if (max_dif_thr == 0) 
		max_dif_thr = L - 10;
	// cout << "max_dif_thr: " << max_dif_thr << endl;
	cout << "reads length: " << L << endl;

	getReads(infile.c_str());

	fprintf(stderr, "ispe: %d\nisorder: %d\n", ispe, isorder);

	if (ispe) {
		// cout << "max_rid: " << max_rid << endl;
		getRightReads(infile1.c_str());
	}

	// RN 
	if (RN > ((1UL << 32)-20)) {
		cout << "Large number of reads; the program cannot compression them.\n" <<endl;
		exit(1);
	}
	max_rid = (uint32_t)RN;
	// outputfn = argv[3];

	string parfn = folder + "par.txt";
	fppar = fopen(parfn.c_str(), "w");
	fprintf(fppar, "%d %d %d\n", L, ispe, isorder);
	// fprintf(fp, "%d\n", ispe);
	// fprintf(fp, "%d\n", isorder);
	if (isorder || ispe) {
		fprintf(fppar, "%lu\n", max_rid);
	}
	fclose(fppar);

	if (isorder && !ispe) {
	}
	// reorderreads();

	cout << "Number of reads: " << max_rid << endl;
	cout << "Time of read file = " << stopwatch.stop() << std::endl;
	stopwatch.resume();
	// nthreads = 1;
// #define _LOADTREE
#ifndef _LOADTREE
	// ----------
	int largekmernum = 0;
	// int maxkmer = 33;
	// kmervecsize = 20;
	// int maxkmer = 31;

	int maxkmer = 29;
	kmervecsize = 20;
	// kmervecsize = 15;

	// int maxkmer = 49;
	// kmervecsize = 40;

	// kmervecsize = 10;
	// kmervecsize = 20;
	// kmervecsize = 3;

	if (L < 75) {
		// maxkmer = 24;
		// maxkmer = 27;
		// kmervecsize = 20;
	}
	if (L > 120) {
		// maxkmer = 50;
		// kmervecsize = 20;
	}
	kmervec = new int[kmervecsize];

	kmervec[0] = maxkmer;
	for (int i = 1; i < kmervecsize; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}
	
	max_kmer = kmervec[0];
	min_kmer = kmervec[kmervecsize - 1];

	cout << "max_kmer: " << max_kmer << endl;
	cout << "min_kmer: " << min_kmer << endl;

	fprintf(stderr, "kmervecsize: %d\n", kmervecsize);
	for (int i = 0; i < kmervecsize; ++i) {
		fprintf(stderr, "%d, ", kmervec[i]);
	}
	fprintf(stderr, "\n");

	// char str[1<<10];
	// cout << "pls input a string: ";
	// scanf("%s", str);
	reads = new READS_t[max_rid + 1];
	
	isnextrnd = new bool[max_rid];
	memset(isnextrnd, 0, sizeof(bool) * max_rid);

	mini = new uint64_t[max_rid];
	snum = new uint32_t[max_rid];
	
	bmtx = new mutex[1 << bsize];
	removeDuplicate(); // kmer = 31

	cout << "isnextrnd: " << coutIsNext() << endl;

	edgesConstruction();

	cout << "after edgesConstruction()..." << endl;

	delete[] bmtx;

	delete[] mini; mini = NULL;
	delete[] snum; snum = NULL;
	delete[] isnextrnd; isnextrnd = NULL;

	mstConstruction();

	// labelRoot();

	// countTree();
	// store trees
	// storeTrees("test.tree");
	/*if (ispe) {
		if (isorder)
			storeTrees(infile + ".tree.pe.order");
		else storeTrees(infile + ".tree.pe");
	} else {
		if (isorder)
			storeTrees(infile + ".tree.order");
		else storeTrees(infile + ".tree");

	}*/
	// storeTrees("SRR870667_2.noN.tree");
	// storeTrees("SRR870667_2.noN.noDup.tree");

	// return 0;

#endif

#ifdef _LOADTREE
	// loadTrees("SRR870667_2.noN.tree");
	if (ispe) {
		if (isorder)
			loadTrees(infile + ".tree.pe.order");
		else	
			loadTrees(infile + ".tree.pe");
	} else {
		if (isorder)
			loadTrees(infile + ".tree.order");
		else	
			loadTrees(infile + ".tree");
	}

	cout << "end _LOADTREE" << endl;
	// labelRoot();
	// countTree();
#endif

	// changeRootNode_v3();

	// selectNoMismatchEdges();
#ifdef DFSENCODING
	// encodeAndSortChild();
	// if (isorder) swappIDforOrder();
	encodeReads();
#endif

#ifndef _LOADTREE
	// cout << "after labelRoot()..." << endl;
#endif	
	// cutSomeEdges();
	// if (isorder) {
	// 	encodeNode();
	// } else 

	if (ispe) {
		// outputPE();
		if (isorder) {
			outputPEOrder();
		} else
			outputPEX();
	} else {
		if (isorder) {
			outputSingleOrder();
		} else
			// outputSingle();
#ifdef DFSENCODING
			outputSingleDFS();
			// outputSingleDFS_v2();
#endif

#ifndef DFSENCODING			
			// outputSingleDFS_v1();
			outputSingleV1();
#endif
	}

	string cmd = "rm -rf " + folder;
	system(cmd.c_str());

	return 0;
}

int main(int argc, char *argv[]){
	if (argc < 6) {
		show_usage(argv[0]);
		exit(1);
	}
	
	if (strcmp(argv[1], "-h") == 0) {
		show_usage(argv[0]);
		exit(0);
	}


	if (strcmp(argv[1], "e") == 0) {
		compress_main(argc, argv);
	} else 
	if (strcmp(argv[1], "d") == 0) {
		decompress_main(argc, argv);
	} else {
		show_usage(argv[0]);
		exit(1);
	}

	return 0;
}


void storeTrees(string fn) {
	cout << "store fn: " << fn << endl;
	FILE *fp = fopen(fn.c_str(), "w");
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		fprintf(fp, "%u %d %d %d", reads[rid].prid, reads[rid].isrc, reads[rid].shift, reads[rid].dn);
		if (reads[rid].dn > 0) {
			if (isorder) {
				for (uint32_t w = 1; w < reads[rid].dn + 1; ++w) {
					fprintf(fp, " %u %d", reads[rid].dup[w].id, reads[rid].dup[w].isrc);
				}
			} else {
				for (uint32_t w = 0; w < reads[rid].dn; ++w) {
					fprintf(fp, " %u %d", reads[rid].dup[w].id, reads[rid].dup[w].isrc);
				}
			}
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

int coutMismatch(char *parent, char *child, int8_t shift) {
	int weight = 0, i, j;
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
	// weight = weight * 2 + shift;
	// weight*(8 + 8 + 1) + shift*8 + 8
	// if (shift != 0) ++ weight;
	return weight;
}

int courMismatch(uint32_t rid) {
	char *str = (char*)alloca((L + 3) * sizeof(char));
	strcpy(str, seq[rid].seq);
	if (reads[rid].isrc) {
		reverseComplement(str);
	}
	return coutMismatch(seq[reads[rid].prid].seq, str, reads[rid].shift);
}

void loadTrees(string fn) {
	cout << "load fn: " << fn << endl;
	int isrc, shift;
	FILE *fp = fopen(fn.c_str(), "r");
	reads = new READS_t[max_rid + 1];	
	cout << "max_rid: " << max_rid << endl;
	uint32_t did;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		reverseReads(seq[rid].seq);
		reads[rid].prid = rootarr[rid] = rid;
		kv_init(reads[rid].crid);
	}
	isnextrnd = new bool[max_rid];
	memset(isnextrnd, 0, sizeof(bool) * max_rid);
	char *str = new char[L + 1];

	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		fscanf(fp, "%u%d%d%d", &reads[rid].prid, &isrc, &shift, &reads[rid].dn);
		reads[rid].isrc = isrc;
		reads[rid].shift = shift;
		if (reads[rid].dn > 0) {
			reads[rid].dup = new dup_t[reads[rid].dn + 1];
			if (isorder) {
				for (uint32_t w = 1; w < reads[rid].dn + 1; ++w) {
					fscanf(fp, "%u%d", &did, &isrc);
					reads[rid].dup[w] = dup_t(did, isrc);
					isnextrnd[did] = true;
				}
			} else {
				for (uint32_t w = 0; w < reads[rid].dn; ++w) {
					fscanf(fp, "%u%d", &did, &isrc);
					reads[rid].dup[w] = dup_t(did, isrc);
					isnextrnd[did] = true;
				}
			}
		}

		strcpy(str, seq[rid].seq);
		bool flag = true;
		if (strcmp(str, seq[reads[rid].prid].seq) == 0) {
			flag = false;
		}

		reverseComplement(str);
		if (strcmp(str, seq[reads[rid].prid].seq) == 0) {
			flag = false;
		}

		if (reads[rid].prid != rid && flag) {
			kv_push(uint32_t, reads[reads[rid].prid].crid, rid);
		}
	}
	fclose(fp);
	cout << "loadTrees() over" << endl;

#ifdef false
	FILE *fpout = fopen("a.txt", "w");
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		// if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads	
		if (reads[rid].prid == rid && reads[rid].crid.n > 0) { //is the root node && not a leaf node == not a singleton reads	
			queue<uint32_t> q;
			q.push(rid);
			uint32_t cnt = 0;
			while (!q.empty()) {
				uint32_t noderid = q.front();		
				q.pop();

				++ cnt;
				for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}

			if (cnt > 500) {
				uint32_t noderid = rid;
				q.push(rid);
				while (!q.empty()) {
					noderid = q.front();
					q.pop();

					for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
						uint32_t crid = reads[noderid].crid.a[i];
						q.push(crid);

						int misno = courMismatch(crid);

						fprintf(fpout, "%d %d\n", + reads[crid].shift, misno);
					}
					fprintf(fpout, "----\n");
				}

				exit(0);
			}
		}
	}
	fclose(fpout);
#endif
}

