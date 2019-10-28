# include "mstcom.h"
using namespace std;

int nthreads;
int max_dif_thr;
int L;
size_t RN;
uint32_t max_rid;
bseq1_t *seq;
 
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

int kmer, max_kmer, min_kmer;
int *kmervec, kmervecsize;

mm128_t *min128vec;
mm192_t *min192vec;
uint64_t *mini;
uint32_t minisize;

mutex *bmtx;

std::string folder;

READS_t *reads;
bool *isnextrnd;
uint32_t *prid2;
int16_t *shift2;
bool *isrc2;
int *min_dif2;
bool *iscombined;
bool *isupdate;

char complement[256];

CStopWatch stopwatch;

inline void connectWithPrid2();

char *cmd;

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
	size_t num_of_nodes = 0;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid) { //a root node 
			uint32_t noderid = rid;
			queue<uint32_t> q;
			q.push(noderid);
			while (!q.empty()) {
				noderid = q.front();
				reads[noderid].root = rid;

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
	rid_pthread = 0;
	std::vector<thread> threadVec;
	for (int i = 0; i < nthreads; ++i) {
		threadVec.push_back(std::thread(countTreeFun));
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
	fprintf(stderr, "trees number: %d\n", tree_num);
}

string fnv[65];

inline void checkNextRnd() {
	memset(isupdate, 0, sizeof(bool)*max_rid);
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[prid2[rid]].root == reads[rid].root) {
			isupdate[rid] = true;
			min_dif2[rid] = L;
		}
	}
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

inline void indexConstruction() {
	CStopWatch tstopwatch;
	tstopwatch.start();
	stopwatch.start();

	if (min_kmer <= 31) {
		B = (mm128_v*)calloc(1 << bsize, sizeof(mm128_v));
	}
	if (max_kmer > 31) {
		BL = (mm192_v*)calloc(1 << bsize, sizeof(mm192_v));
	}

	min128vec = NULL;
	min192vec = NULL;
	string idxfn;

	for (int kid = 0; kid < kmervecsize; ++kid) {
		kmer = kmervec[kid];

		cout << "kmervecsize: " << kmervecsize << endl;
		cout << "kid: " << kid << endl;
		cout << "kmer: " << kmer << endl;
		tstopwatch.resume();

		if (kmer <= 31) {
			if (NULL == min128vec) {
				min128vec = new mm128_t[max_rid];
			}
			for (int i = 0; i < (1 << bsize); ++i) {
				kv_init(B[i]);
			}
		} else {
			if (NULL == min192vec) {
				min192vec = new mm192_t[max_rid];
			}
			for (int i = 0; i < (1 << bsize); ++i) {
				kv_init(BL[i]);
			}
		}
		// cout << "kmer: " << kmer << endl;
		calcMinimizers();

		cout << "Time of calcMinimizers() = " << tstopwatch.stop() << std::endl;
		tstopwatch.resume();
		// cout << "after calcMinimizers()\n";
		sortBuckets();

		cout << "Time of sortBuckets() = " << tstopwatch.stop() << std::endl;
		tstopwatch.resume();
		// cout << "after sortBuckets()\n";
		processBuckets();

		cout << "Time of processBuckets() = " << tstopwatch.stop() << std::endl;
		tstopwatch.resume();
		// cout << "after processBuckets()\n";

		// if (kmer <= 31) { 
		// 	for (int i = 0; i < (1 << bsize); ++i) {
		// 		kv_destroy(B[i]);
		// 	}
		// } else {
		// 	for (int i = 0; i < (1 << bsize); ++i) {
		// 		kv_destroy(BL[i]);
		// 	}	
		// }

		obtainMiniIdx();

		// storeIdx();
		idxfn = folder + to_string(kmer) + ".min";
		idxDump(idxfn);

		// idxLoad(idxfn);
		// exit(0);
		if (kmer <= 31) { 
			for (int i = 0; i < (1 << bsize); ++i) {
				kv_destroy(B[i]);
				kv_init(B[i]);
			}
		} else {
			for (int i = 0; i < (1 << bsize); ++i) {
				kv_destroy(BL[i]);
				kv_init(BL[i]);
			}	
		}

		calcMaximizers();
		sortBuckets();
		processBuckets();
		obtainMiniIdx();

		idxfn = folder + to_string(kmer) + ".max";
		idxDump(idxfn);

		if (kmer <= 31) { 
			for (int i = 0; i < (1 << bsize); ++i) {
				kv_destroy(B[i]);
			}
		} else {
			for (int i = 0; i < (1 << bsize); ++i) {
				kv_destroy(BL[i]);
			}	
		}

		if (kmer <= 31) {
			if (NULL != min192vec) {
				delete[] min192vec;
				min192vec = NULL;
			}
		} else {
			if (NULL != min128vec) {
				delete[] min128vec;
				min128vec = NULL;
			}
		}

		cout << "kmer: " << kmer << "; isnextrnd: " << coutIsNext() << endl;

		//
		// clearDuplicate();

		// cout << "Time of initialTreeConstuction kmer =" << kmer << " = " << stopwatch.stop() << std::endl;
		// stopwatch.resume();
		// -- kmer;
	}
	tstopwatch.resume();

	if (min_kmer <= 31) {
		free(B);
	}
	if (max_kmer > 31) {
		free(BL);
	}

	if (NULL != min128vec) {
		delete[] min128vec;
		min128vec = NULL;
	}
	if (NULL != min192vec) {
		delete[] min192vec;
		min192vec = NULL;
	}

	cout << "Time of indexConstruction(): " << stopwatch.stop() << std::endl;
	stopwatch.resume();
}

// #ifdef false
inline void treeConstuction() {

	CStopWatch tstopwatch;
	tstopwatch.start();
	stopwatch.start();

	// prid2 = new uint32_t[max_rid];
	shift2 = new int16_t[max_rid];
	isrc2 = new bool[max_rid];
	iscombined = new bool[max_rid];
	isupdate = new bool[max_rid];
	string idxfn;
	int pre_trees_cnt = 0;
	while (true) {
		stopwatch.resume();
	
		checkNextRnd();
		// cout << "before prid2[1539203]: " << prid2[1539203] << endl;	
		// kmer = max_kmer;
		// fprintf(stderr, "kmer: %lu\n", kmer);
		// while (kmer >= min_kmer) {
		for (int kid = 0; kid < kmervecsize; ++kid) {
			kmer = kmervec[kid];
			cout << "kmer: " << kmer << endl;
			// stopwatch.resume();
			
			idxfn = folder + to_string(kmer) + ".min";
			// cout << "before mm_bucket_load()\n";
			if (idxLoad(idxfn) == false) {
				cout << "the function idxLoad() fail...\n";
				string cmd = "rm -rf " + folder;
				system(cmd.c_str());
				exit(EXIT_FAILURE);
			}

			// cout << kmer << " before processBuckets2()\n";
			collectNext();
			// cout << kmer << " after processBuckets2()\n";
			
			idxfn = folder + to_string(kmer) + ".max";
			// cout << "before mm_bucket_load()\n";
			if (idxLoad(idxfn) == false) {
				cout << "the function idxLoad() fail...\n";
				string cmd = "rm -rf " + folder;
				system(cmd.c_str());
				exit(EXIT_FAILURE);
			}

			// cout << kmer << " before processBuckets2()\n";
			collectNext();
			// processLargeBuckets2();
			// cout << kmer << " after processLargeBuckets2()\n";

			// cout << "Time of findPrid2 kmer = " << kmer << " = " << stopwatch.stop() << std::endl;
			// stopwatch.resume();
		}
		/*FILE *fp = fopen("dis.txt", "w");
		for (uint32_t rid = 0; rid < max_rid; ++rid) {
			fprintf(fp, "%d\n", min_dif2[rid]);
		}
		fclose(fp);*/
		// exit(0);

		// cout << "111\n";
		tstopwatch.resume();

		connectWithPrid2();
		// connectWithPrid2(1);

		cout << "Time of connectWithPrid2() = " << tstopwatch.stop() << std::endl;
		tstopwatch.resume();
		// cout << "after connectWithPrid2()\n";
		// cutCircle();
		// cout << "after cutCircle()\n";
		// countTree();
		labelRoot();

		cout << "Time of findPrid2 = " << stopwatch.stop() << std::endl;
		stopwatch.resume();

		countTree();

		cout << "coutIsNext(): " << coutIsNext() << endl;
		if ((max_rid > 33554432 && abs(pre_trees_cnt - tree_num) < 100000) || // 2^25
			(max_rid <= 33554432 && max_rid > 8388608 && abs(pre_trees_cnt - tree_num) < 10000) || // 2^23
			(max_rid <= 8388608 && abs(pre_trees_cnt - tree_num) < 1000)) break;

		// if (tree_num < 20000) break;
		// if (tree_num < 200000) break;
		// if (tree_num < 2000000) break;
		// if (abs(pre_trees_cnt - tree_num) < 10000) break;
		// if (abs(pre_trees_cnt - tree_num) < 1000000) break;
		// if (tree_num < 1000000) break;
		// if (abs(pre_trees_cnt - tree_num) < 1000000) break;
		// if (abs(pre_trees_cnt - tree_num) < 10000 || tree_num <= 10000000) break;
		// if (abs(pre_trees_cnt - tree_num) < 10000) break;
		// if (tree_num <= 1000000) break;
		// if (abs(pre_trees_cnt - tree_num) < 100) break;
		pre_trees_cnt = tree_num;
	}

	delete[] shift2;
	delete[] isrc2;
	delete[] iscombined;
	delete[] isupdate;

	countTree();

	return;
}
// #endif

bool cmp2(const ROOTNODE_t &a, const ROOTNODE_t &b) {
	return a.nodecnt > b.nodecnt;
}

inline void reRootNode(const size_t &_rid, const size_t &root);
inline void reRootNode(const size_t &_rid);

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

inline void reRootNode(const uint32_t &_rid) { // set rid as root node
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

inline void connectWithPrid2() {
	uint32_t prid, rid;
	CStopWatch tstopwatch;
	tstopwatch.start();

	Edge_v edges;
	kv_init(edges);
	kv_resize(Edge_t, edges, max_rid);
	Edge_t teds;

	// FILE *fp = fopen("dif.txt", "w");

	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		// if (reads[prid2[rid]].root != reads[rid].root && min_dif2[rid] > 0 && min_dif2[rid] < max_dif_thr) {
		if (reads[prid2[rid]].root != reads[rid].root && !isnextrnd[rid]) {
			if (min_dif2[rid] < max_dif_thr) {
				teds.rid = rid, teds.dif = min_dif2[rid];
				kv_push(Edge_t, edges, teds);
			} else {
				isnextrnd[rid] = true; // hard reads; do not consider never
			}
			// fprintf(fp, "%u: %u; %d\n", rid, prid2[rid], min_dif2[rid]);
		} 
		else {
			isnextrnd[rid] = true;
		}
	}
	// fclose(fp);
	// exit(0);

	radix_sort_edge(edges.a, edges.a + edges.n);

	// for (int i = 0; i < 5; ++i) {
	// 	cout << edges.a[i].dif << endl;
	// }
	// sort(edges.a, edges.a + edges.n, edgecmp);

	cout << "Time of radix_sort_edge = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();

	cout << "edges.n: " << edges.n << endl;
	// edges.n: 2283538
	uint32_t root1, root2;

	memset(iscombined, false, sizeof(bool)*max_rid);

	for (size_t i = 0; i < edges.n; ++i) {
		// cout << "i: " << i << endl;
		rid = edges.a[i].rid;
		prid = prid2[rid];
	
		if (!iscombined[reads[rid].root]) { //|| !iscombined[reads[prid].root]) {
			reRootNode(rid); //
			if (debug) {
				cout << "after reRootNode()..\n";
			}
			// 
			reads[rid].prid = prid;
			reads[rid].shift = shift2[rid];
			reads[rid].isrc = isrc2[rid];

			kv_push(uint32_t, reads[prid].crid, rid);
			//
			iscombined[reads[rid].root] = true;
			iscombined[reads[prid].root] = true;
		}
		// if (edges.a[i].dif > 0.6*L) break;
	}

	kv_destroy(edges);

	cout << "Time of connectWithPrid2 = " << tstopwatch.stop() << std::endl;
	tstopwatch.resume();
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

inline uint32_t encode(uint32_t rid, char *en_str) {
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

const char invert_code_rule[5] = {'A', 'C', 'G', 'T', 'N'};

void cutSomeEdges() {
	char *encstr = (char*)alloca(((L<<2) + 2) * sizeof(char));
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (rid != reads[rid].prid) {
			encode(rid, encstr);
			if (strlen(encstr) > 0.6*L) {
				removeChild(reads[rid].prid, rid);
				reads[rid].prid = rid;
			}
		}
	}
}

void changeNodeXX() {
	// return;
	uint32_t root, noderid, i, trid;
	uint32_t orienclen, aftenclen;
	int **cnt = new int*[5], max_count;
	for (i = 0; i < 5; ++i) {
		cnt[i] = new int[L];
	}
	char *str, tc;
	bool ischanged;

	char *newstr = (char*)calloc((L + 2), sizeof(char));
	char *oristr = (char*)calloc(((L<<2) + 2), sizeof(char));
	char *encstr = (char*)calloc(((L<<2) + 2), sizeof(char));
	cout << "begin changeNode()\n";

	FILE *fpencoded = fopen("encode.txt", "w");

	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid) { // root node 
			// vector<uint32_t> v;
			// v.push_back(rid);	
			uint32_v v;
			// cout << "xxx\n";
			kv_init(v);
			// cout << "rid: " << rid << endl;
			kv_push(uint32_t, v, rid);
			i = 0;
			while (i < v.n) {
				noderid = v.a[i];
				// cout << "noderid: " << noderid << endl;
				for (uint32_t j = 0; j < reads[noderid].crid.n; ++j) {
					// v.push_back(reads[noderid].crid.a[j]);
					kv_push(uint32_t, v, reads[noderid].crid.a[j]);
				}
				++ i;
				// cout << "2 noderid: " << noderid << endl;
			}

			// cout << "v.size(): " << v.n << endl;
			// for (i = 0; i < v.n; ++i) {
			// 	cout << v.a[i] << ", ";
			// }
			// cout << endl;

			// i = v.n - 1;
			// cout << "i: " << i << endl;

			// while (i >= 0) {
			for (i = v.n - 1; i >= 0; --i) {
				// cout << v.n << "-i: " << i << endl;
				root = v.a[i];
				bool debug = false;
				if(strcmp(seq[root].seq, "ATTACATTGGATTCCATTCTATGATTCCATTCAATTCCATTCGTTGATGATTCGATTCCATTCAATGATGATTCCATTTGAGTTCATTCGATGATTCTATT") == 0)
						// debug = true;
						debug = false;
				// cout << "root: " << root << endl;
				if (reads[root].crid.n > 1) { // at least 2 child
					// - - do change
					for (int k = 0; k < 5; ++k) {
						memset(cnt[k], 0, sizeof(int) * L);
					}
					// cout << v.n << "-i: " << i << " - -" << endl;

					for (uint32_t j = 0; j < reads[root].crid.n; ++j) {
						trid = reads[root].crid.a[j];
						// cout << "j:" << j << "; " << v.n << "| before -i: " << i << endl;
						// str = seq[trid].seq;

						strncpy(newstr, seq[trid].seq, L);
						newstr[L] = '\0';
						// cout << "j:" << j << "; " << v.n << "| before -i: " << i << endl;
						if (reads[trid].isrc) {
							reverseComplement(newstr);
						}
						str = newstr;

						// cout << "j:" << j << "; " << v.n << "| before -i: " << i << endl;

						// cout << "shift: " << reads[trid].shift << endl;
						if (reads[trid].shift >= 0) {
							for (int k = 0; k < L - reads[trid].shift; ++k) {
								// cout << str[k] << "-" << (uint8_t)str[k] << " - [" << (int)seq_nt4_table[(uint8_t)str[k]] << ", " << reads[trid].shift + k << "]\n";
								++ cnt[(int)seq_nt4_table[(uint8_t)str[k]]][reads[trid].shift + k];
							}
						} else {
							for (int k = 0-reads[trid].shift; k < L; ++k) {
								// cout << str[k] << "-" << (uint8_t)str[k] << " - [" << (int)seq_nt4_table[(uint8_t)str[k]] << ", " << k << "]\n";
								++ cnt[(int)seq_nt4_table[(uint8_t)str[k]]][k];
							}
						}
						// cout << v.n << "| after -i: " << i << endl;
					}

					// cout << v.n << "-i: " << i << endl;

					if (debug && false) {
						for (int k = 0; k < 5; ++k) {
							cout << invert_code_rule[k] << ": ";
							for (int s = 0; s < 7; ++s) {
								cout << cnt[k][s] << " ";
							}
							cout << endl;
						}
					}

					// cout << v.n << "-i: " << i << endl;

					str = seq[root].seq;
					strcpy(newstr, str);
					ischanged = false;
					for (int j = 0; j < L; ++j) {
						max_count = cnt[0][j];
						tc = 'A';
						for (int k = 1; k < 5; ++k) {
							if (cnt[k][j] > max_count) {
								max_count = cnt[k][j];
								tc = invert_code_rule[k];
							}
						}
						if (tc != newstr[j]) {
							newstr[j] = tc;
							ischanged = true;
						}
					}
					// cout << v.n << "-i: " << i << endl;
					// cout << "zzz\n";
					if (ischanged) {
						// original encoding
						orienclen = 0;
						if (reads[root].prid != root) {
							orienclen += encode(root, encstr);
						} else {
							orienclen += L;
						}
						// cout << "zzz0000\n";
						if (debug) cout << "oristr: " << seq[root].seq << endl;
						for (uint32_t j = 0; j < reads[root].crid.n; ++j) {
							orienclen += encode(reads[root].crid.a[j], encstr);
							if (debug) cout << seq[reads[root].crid.a[j]].seq << " | " << encstr << endl;
						}
						// cout << "zzz000011111\n";
						//
						aftenclen = 0;
						strcpy(oristr, str);
						strcpy(seq[root].seq, newstr);
						if (reads[root].prid != root) {
							aftenclen += encode(root, encstr);
						} else {
							aftenclen += L;
						}
						if (debug) {
							cout << "-----\n";
							cout << "newstr: " << seq[root].seq << endl;
						}
						for (uint32_t j = 0; j < reads[root].crid.n; ++j) {
							aftenclen += encode(reads[root].crid.a[j], encstr);
							if (debug) cout << seq[reads[root].crid.a[j]].seq << " | " << encstr << endl;
						}

						// cout << "zzz111\n";

						if (aftenclen < orienclen) {
							// encode children 
							for (uint32_t j = 0; j < reads[root].crid.n; ++j) {
								trid = reads[root].crid.a[j];
								str = seq[trid].seq;
								encode(trid, encstr);
								if (str[L] == '|') {
									strcat(encstr, str+L);
								}
								free(seq[trid].seq);
								seq[trid].seq = strdup(encstr);
								// fprintf(fpencoded, "%s\n", encstr);
								// strcpy(seq[trid].seq, encstr);
							}
							//
							free(seq[root].seq);
							// cout << "zzz11100000000\n";

							encode(oristr, newstr, 0, encstr);
							oristr[L + 1] = '\0'; 
							oristr[L] = '|'; 
							strcat(oristr, encstr);

							seq[root].seq = strdup(oristr);
							fprintf(fpencoded, "%s\n", seq[root].seq);
							// exit(0);
						}

						// cout << "zzz222\n";
					} 
					// cout << "zzz333\n";
					// else do nothing
				}
				// cout << v.n << "-i: " << i << endl;
				// cout << "zzz444\n";
				if (i == 0) break;
			}
			// cout << "zzz555\n";
			// v.clear();
			// cout << "v.n: " << v.n << endl;
			kv_destroy(v);
			// cout << "zzz556666\n";
		}
	}

	free(newstr);
	free(oristr);
	free(encstr);
	fclose(fpencoded);
	cout << "after changeNode()\n";
}

// no consuses
void encodeNode() {
	// return;
	uint32_t root, noderid, i, trid;
	uint32_t orienclen, aftenclen;
	int **cnt = new int*[5], max_count;
	for (i = 0; i < 5; ++i) {
		cnt[i] = new int[L];
	}
	char *str, tc;
	bool ischanged;

	// char *newstr = (char*)calloc((L + 2), sizeof(char));
	// char *oristr = (char*)calloc(((L<<2) + 2), sizeof(char));
	// char *encstr = (char*)calloc(((L<<2) + 2), sizeof(char));	
	char *newstr = (char*)alloca((L + 2) * sizeof(char));
	char *oristr = (char*)alloca(((L<<2) + 2) * sizeof(char));
	char *encstr = (char*)alloca(((L<<2) + 2) * sizeof(char));

	cout << "begin encodeNode()\n";

	FILE *fpencoded = fopen("encode.txt", "w");

	// #pragma omp parallel for
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid) { // root node 
			// vector<uint32_t> v;
			// v.push_back(rid);	
			uint32_v v;
			// cout << "xxx\n";
			kv_init(v);
			// cout << "rid: " << rid << endl;
			kv_push(uint32_t, v, rid);
			i = 0;
			while (i < v.n) {
				noderid = v.a[i];
				// cout << "noderid: " << noderid << endl;
				for (uint32_t j = 0; j < reads[noderid].crid.n; ++j) {
					// v.push_back(reads[noderid].crid.a[j]);
					kv_push(uint32_t, v, reads[noderid].crid.a[j]);
				}
				++ i;
				// cout << "2 noderid: " << noderid << endl;
			}

			// cout << "v.size(): " << v.n << endl;
			// for (i = 0; i < v.n; ++i) {
			// 	cout << v.a[i] << ", ";
			// }
			// cout << endl;

			// i = v.n - 1;
			// cout << "i: " << i << endl;

			// while (i >= 0) {
			for (i = v.n - 1; i >= 0; --i) {
				// cout << v.n << "-i: " << i << endl;
				root = v.a[i];
				bool debug = false;
				// if(strcmp(seq[root].seq, "ATTACATTGGATTCCATTCTATGATTCCATTCAATTCCATTCGTTGATGATTCGATTCCATTCAATGATGATTCCATTTGAGTTCATTCGATGATTCTATT") == 0)
						// debug = true;
				// if (v.n == 2060)

				for (uint32_t j = 0; j < reads[root].crid.n; ++j) {
					trid = reads[root].crid.a[j];
					str = seq[trid].seq;
					encode(trid, encstr);
					if (str[L] == '|') {
						strcat(encstr, str + L);
					}
					free(seq[trid].seq);
					seq[trid].seq = strdup(encstr);
				}
				// cout << v.n << "-i: " << i << endl;
				if (debug) cout << "zzz444\n";

				// for duplicate reads
				for (uint32_t w = 1; w < reads[root].dn + 1; ++w) {
					trid = reads[root].dup[w].id;
					reads[trid].prid = root;		
					reads[trid].isrc = reads[root].dup[w].isrc;
					seq[trid].seq[0] = '\0'; 
				}

				if (i == 0) break;
			}
			// cout << "zzz555\n";
			// v.clear();
			// cout << "v.n: " << v.n << endl;
			kv_destroy(v);
			// cout << "zzz556666\n";
		}
	}
 	// cout << "999999\n";
	// free(newstr);
	// free(oristr);
	// free(encstr);
	fclose(fpencoded);
	cout << "after encodeNode()\n";
}

// change some node to consuses read
void changeNode() {
	// return;
	uint32_t root, noderid, i, trid;
	uint32_t orienclen, aftenclen;
	int **cnt = new int*[5], max_count;
	for (i = 0; i < 5; ++i) {
		cnt[i] = new int[L];
	}
	char *str, tc;
	bool ischanged;

	// char *newstr = (char*)calloc((L + 2), sizeof(char));
	// char *oristr = (char*)calloc(((L<<2) + 2), sizeof(char));
	// char *encstr = (char*)calloc(((L<<2) + 2), sizeof(char));	
	char *newstr = (char*)alloca((L + 2) * sizeof(char));
	char *oristr = (char*)alloca(((L<<2) + 2) * sizeof(char));
	char *encstr = (char*)alloca(((L<<2) + 2) * sizeof(char));

	cout << "begin changeNode()\n";

	FILE *fpencoded = fopen("encode.txt", "w");

	// #pragma omp parallel for
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid) { // root node 
			// vector<uint32_t> v;
			// v.push_back(rid);	
			uint32_v v;
			// cout << "xxx\n";
			kv_init(v);
			// cout << "rid: " << rid << endl;
			kv_push(uint32_t, v, rid);
			i = 0;
			while (i < v.n) {
				noderid = v.a[i];
				// cout << "noderid: " << noderid << endl;
				for (uint32_t j = 0; j < reads[noderid].crid.n; ++j) {
					// v.push_back(reads[noderid].crid.a[j]);
					kv_push(uint32_t, v, reads[noderid].crid.a[j]);
				}
				++ i;
				// cout << "2 noderid: " << noderid << endl;
			}

			// cout << "v.size(): " << v.n << endl;
			// for (i = 0; i < v.n; ++i) {
			// 	cout << v.a[i] << ", ";
			// }
			// cout << endl;

			// i = v.n - 1;
			// cout << "i: " << i << endl;

			// while (i >= 0) {
			for (i = v.n - 1; i >= 0; --i) {
				// cout << v.n << "-i: " << i << endl;
				root = v.a[i];
				bool debug = false;
				// if(strcmp(seq[root].seq, "ATTACATTGGATTCCATTCTATGATTCCATTCAATTCCATTCGTTGATGATTCGATTCCATTCAATGATGATTCCATTTGAGTTCATTCGATGATTCTATT") == 0)
						// debug = true;
				// if (v.n == 2060)
						debug = false;

				if (debug) cout << "root: " << root << endl;
				if (reads[root].crid.n > 1) { // at least 2 child
					// - - do change
					for (int k = 0; k < 5; ++k) {
						memset(cnt[k], 0, sizeof(int) * L);
					}
					// cout << v.n << "-i: " << i << " - -" << endl;

					for (uint32_t j = 0; j < reads[root].crid.n; ++j) {
						trid = reads[root].crid.a[j];
						// cout << "j:" << j << "; " << v.n << "| before -i: " << i << endl;
						// str = seq[trid].seq;

						strncpy(newstr, seq[trid].seq, L);
						newstr[L] = '\0';
						// cout << "j:" << j << "; " << v.n << "| before -i: " << i << endl;
						if (reads[trid].isrc) {
							reverseComplement(newstr);
						}
						str = newstr;

						// cout << "j:" << j << "; " << v.n << "| before -i: " << i << endl;

						// cout << "shift: " << reads[trid].shift << endl;
						if (reads[trid].shift >= 0) {
							for (int k = 0; k < L - reads[trid].shift; ++k) {
								// cout << str[k] << "-" << (uint8_t)str[k] << " - [" << (int)seq_nt4_table[(uint8_t)str[k]] << ", " << reads[trid].shift + k << "]\n";
								++ cnt[(int)seq_nt4_table[(uint8_t)str[k]]][reads[trid].shift + k];
							}
						} else {
							for (int k = 0-reads[trid].shift; k < L; ++k) {
								// cout << str[k] << "-" << (uint8_t)str[k] << " - [" << (int)seq_nt4_table[(uint8_t)str[k]] << ", " << k << "]\n";
								++ cnt[(int)seq_nt4_table[(uint8_t)str[k]]][k];
							}
						}
						// cout << v.n << "| after -i: " << i << endl;
					}

					// cout << v.n << "-i: " << i << endl;

					if (debug && false) {
						for (int k = 0; k < 5; ++k) {
							cout << invert_code_rule[k] << ": ";
							for (int s = 0; s < 7; ++s) {
								cout << cnt[k][s] << " ";
							}
							cout << endl;
						}
					}

					// cout << v.n << "-i: " << i << endl;

					str = seq[root].seq;
					strcpy(oristr, str);
					strcpy(newstr, str);

					ischanged = false;
					for (int j = 0; j < L; ++j) {
						max_count = cnt[0][j];
						tc = 'A';
						for (int k = 1; k < 5; ++k) {
							if (cnt[k][j] > max_count) {
								max_count = cnt[k][j];
								tc = invert_code_rule[k];
							}
						}
						if (tc != newstr[j]) {
							newstr[j] = tc;
							ischanged = true;
						}
					}
					// if (strcmp(oristr, "GAATGTAGTCGAATGGAATCATCATCGAATGGAGTCGAATGGAAACATCACTGAATGGAATCAAATGAAATCACCGAATTGAATCAAATGGAATGATCATC") == 0) {
					// 	cout << oristr << endl;
					// 	cout << newstr << endl;
					// }
					// cout << v.n << "-i: " << i << endl;
					// cout << "zzz\n";
					if (ischanged) {
						// original encoding
						orienclen = 0;
						/*if (reads[root].prid != root) {
							orienclen += encode(root, encstr);
						} else {
							orienclen += L;
						}*/
						// cout << "zzz0000\n";
						if (debug) cout << "oristr: " << seq[root].seq << endl;
						for (uint32_t j = 0; j < reads[root].crid.n; ++j) {
							trid = reads[root].crid.a[j];
							orienclen += encode(trid, encstr);
							// if (seq[trid].seq[L] == '|') {
							// 	orienclen += strlen(seq[trid].seq + L + 1);
							// }
							if (debug) cout << seq[reads[root].crid.a[j]].seq << " | " << encstr << endl;
						}
						// cout << "zzz000011111\n";
						//
						aftenclen = 0;

						strcpy(seq[root].seq, newstr);
						/*if (reads[root].prid != root) {
							aftenclen += encode(root, encstr);
						} else {
							aftenclen += L;
						}*/
						if (debug) {
							cout << "-----\n";
							cout << "newstr: " << seq[root].seq << endl;
						}
						for (uint32_t j = 0; j < reads[root].crid.n; ++j) {
							aftenclen += encode(reads[root].crid.a[j], encstr);

							if (debug) cout << seq[reads[root].crid.a[j]].seq << " | " << encstr << endl;
						}

						encode(oristr, newstr, 0, encstr);
						aftenclen += strlen(encstr);
						// #ifdef false					
						// cout << "zzz111\n";
						if (aftenclen + 1 < orienclen) {
							// encode children 
							for (uint32_t j = 0; j < reads[root].crid.n; ++j) {
								trid = reads[root].crid.a[j];
								str = seq[trid].seq;
								encode(trid, encstr);
								if (str[L] == '|') {
									strcat(encstr, str + L);
								}
								free(seq[trid].seq);
								seq[trid].seq = strdup(encstr);
								// fprintf(fpencoded, "%s\n", encstr);
								// strcpy(seq[trid].seq, encstr);
							}
							//
							// cout << "zzz11100000000\n";

							encode(oristr, newstr, 0, encstr);
							oristr[L] = '|'; 
							oristr[L + 1] = '\0'; 
							strcat(oristr, encstr);

							free(seq[root].seq);
							seq[root].seq = strdup(oristr);
							fprintf(fpencoded, "%s\n", seq[root].seq);

							// exit(0);
						} /*else {
							strcpy(seq[root].seq, oristr);
						}*/
						else 
						// #endif
						{
							strcpy(seq[root].seq, oristr);
							// encode children 
							for (uint32_t j = 0; j < reads[root].crid.n; ++j) {
								trid = reads[root].crid.a[j];
								str = seq[trid].seq;
								encode(trid, encstr);
								if (str[L] == '|') {
									strcat(encstr, str + L);
								}
								free(seq[trid].seq);
								seq[trid].seq = strdup(encstr);
							}
						}
						
						// cout << "zzz222\n";
					} else {
						for (uint32_t j = 0; j < reads[root].crid.n; ++j) {
							trid = reads[root].crid.a[j];
							str = seq[trid].seq;
							encode(trid, encstr);
							if (str[L] == '|') {
								strcat(encstr, str + L);
							}
							free(seq[trid].seq);
							seq[trid].seq = strdup(encstr);
						}
					}
					// cout << "zzz333\n";
					// else do something
				} else {
					for (uint32_t j = 0; j < reads[root].crid.n; ++j) {
						trid = reads[root].crid.a[j];
						str = seq[trid].seq;
						encode(trid, encstr);
						if (str[L] == '|') {
							strcat(encstr, str + L);
						}
						free(seq[trid].seq);
						seq[trid].seq = strdup(encstr);
					}
				}
				// cout << v.n << "-i: " << i << endl;
				if (debug) cout << "zzz444\n";

				// for duplicate reads
				if (isorder) {
					for (uint32_t w = 1; w < reads[root].dn + 1; ++w) {
						trid = reads[root].dup[w].id;
						reads[trid].prid = root;		
						reads[trid].isrc = reads[root].dup[w].isrc;
						// seq[trid].seq[0] = '\0'; 
						// seq[trid].seq[0] = '-';
						strcpy(seq[trid].seq, "-");
					}
				}

				if (i == 0) break;
			}
			// cout << "zzz555\n";
			// v.clear();
			// cout << "v.n: " << v.n << endl;
			kv_destroy(v);
			// cout << "zzz556666\n";
		}
	}
 	// cout << "999999\n";
	// free(newstr);
	// free(oristr);
	// free(encstr);
	fclose(fpencoded);
	cout << "after changeNode()\n";
}


int compress_main(int argc, char *argv[]) {
	cmd = new char[1024];
	stopwatch.start();

	init();
	getPars(argc, argv);

	L = getReadsLength(infile.c_str());

	max_dif_thr = L * 0.7;
	// max_dif_thr = L * 0.35;
	// max_dif_thr = L * 0.5;
	cout << "max_dif_thr: " << max_dif_thr << endl;
	cout << "L: " << L << endl;

	getReads(infile.c_str());

	fprintf(stderr, "ispe: %d\nisorder: %d\n", ispe, isorder);

	if (ispe) {
		cout << "max_rid: " << max_rid << endl;
		getRightReads(infile1.c_str());
	}

	// RN 
	if (RN > ((1UL << 32)-1)) {
		cout << "reads numer is too large\n" <<endl;
		exit(1);
	}
	max_rid = (uint32_t)RN;
	// outputfn = argv[3];

	string parfn = folder + "par.txt";
	FILE *fp = fopen(parfn.c_str(), "w");
	fprintf(fp, "%d %d %d\n", L, ispe, isorder);
	// fprintf(fp, "%d\n", ispe);
	// fprintf(fp, "%d\n", isorder);
	if (isorder || ispe) {
		fprintf(fp, "%lu\n", max_rid);
	}
	fclose(fp);

	cout << "max_rid: " << max_rid << endl;
	cout << "Time of read file = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	// fprintf(stderr, "%s\n", seq[max_rid - 1].seq);

	// kmervecsize = 1;
	// kmervec = new int[kmervecsize];
	// kmervec[0] = 21;
	// kmervecsize = 7;
	// kmervecsize = 5;
	// int maxkmer = 62;
	// kmervecsize = 47;

	/*kmervecsize = 22;
	kmervec = new int[kmervecsize];

	int maxkmer= 55;
	
	kmervec[0] = 55;
	kmervec[1] = 54;
	kmervec[2] = 53;
	kmervec[3] = 52;
	kmervec[4] = 45;
	kmervec[5] = 44;
	kmervec[6] = 43;
	kmervec[7] = 42;

	kmervec[8] = 29;
	for (int i = 9; i < kmervecsize; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}*/

	/*int largekmernum = 0;
	// int maxkmer = 29;
	// // int maxkmer = 54;
	// kmervecsize = 15;

	int maxkmer = 62;
	// int maxkmer = 54;
	kmervecsize = 42;

	if (L < 75) {
		maxkmer = 62;
		// kmervecsize = 20;
	}
	kmervec = new int[kmervecsize];

	kmervec[0] = maxkmer;
	for (int i = 1; i < kmervecsize; ++i) {
		kmervec[i] = kmervec[i-1] - 1;
	}*/
	// ----------
	int largekmernum = 0;
	int maxkmer = 29;
	kmervecsize = 15;

	if (L < 75) {
		maxkmer = 24;
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
	reads = new READS_t[max_rid];
	prid2 = new uint32_t[max_rid];
	min_dif2 = new int[max_rid];
	
	isnextrnd = new bool[max_rid];
	memset(isnextrnd, 0, sizeof(bool) * max_rid);

	mini = new uint64_t[max_rid];
	
	bmtx = new mutex[1 << bsize];
	removeDuplicate(); // kmer = 31

	cout << "isnextrnd: " << coutIsNext() << endl;
	indexConstruction();
	delete[] bmtx;

	// sprintf(cmd, "rm -rf %s", folder.c_str());
	// system(cmd);
	// exit(0);

	treeConstuction();

	delete[] isnextrnd;
	delete[] prid2;
	delete[] min_dif2;
	
	labelRoot();

	// cutSomeEdges();
	// if (isorder) {
	// 	encodeNode();
	// } else 
		changeNode();

	if (ispe) {
		// outputPE();
		outputPEX();
	} else {
		if (isorder) {
			outputSingleOrder();
		} else
			outputSingle();
	}

	sprintf(cmd, "rm -rf %s", folder.c_str());
	system(cmd);
	delete[] cmd;

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
