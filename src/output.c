#include "mstcom.h"
using namespace mstcom;

void compressFiles(vector<string> &f1s, vector<string> &f2s);
string getExt(string file);

uint32_t MAXNO = 1<<20;

// void outputSingleDFS_best_version_v0() {
void outputSingleDFS() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));
	char *encstr = (char*)alloca((L + 1) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	
	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;

	uint8bit_v isleftbin;
	isleftbin.n = isleftbin.a = 0;
	
	uint8bit_v isdupbin;
	isdupbin.n = isdupbin.a = 0;

	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string isleftbinfn = folder + "isleft.bin";
	std::ofstream fpisleft(isleftbinfn.c_str(), std::ios::binary);

	string isdupfn = folder + "isdup.bin";
	std::ofstream fpisdup(isdupfn.c_str(), std::ios::binary);

	string dupfn = folder + "dup.bin";
	std::ofstream fpdup(dupfn, std::ios::binary);

	// FILE *fpw = fopen("eenc.seq", "w");

	// uint32_t numberofbackstep = 0;
	uint32_t dupno = 0, rcdupno = 0;
	// uint32_t dupnum = 0, printnum = 0;
	// cout << "111" << endl;
	for (uint32_t rid = 0; rid < max_rid; ++rid) 
	if (!visited[rid] && reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
		// cout << "134904's parent " << reads[134904].prid << endl;
		// cout << "reads no. in this tree: " << rootnodevec[i].nodecnt << endl;
		queue<uint32_t> q;
		q.push(rid);
		uint32_t cnt = 0;
		while (!q.empty()) {
			uint32_t noderid = q.front();		
			q.pop();
			++cnt;
			visited[noderid] = true;
			// gridvec.push_back(noderid);
			for (size_t i = 0; i < reads[noderid].dn; ++i) {
				visited[reads[noderid].dup[i].id] = true;
			}
			for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
				q.push(reads[noderid].crid.a[i]);
			}
		}

		uint32_t rootid = rid;
		int dir;
		uint32_t dn;
		int preshiftdirection = 0; // 0 for left shift; 1 for right shift;

		uint32_t num = 0;
		while (rootid < max_rid) {
			++ num;
			dn = reads[rootid].dn;

			if (num == 1) {
				fprintf(fpenstr, "%s\n", seq[rootid].seq);
				// ++ printnum;

				// fprintf(fpw, "%s\n", seq[rootid].seq);
				{
					dupno = 0; rcdupno = 0;
					for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
						// fprintf(fpw, "%s\n", seq[rootid].seq);
						if (reads[rootid].dup[w].isrc) ++ rcdupno;
						else ++ dupno;
					}
					// for (uint32_t w = 0; w < dupno; ++w) {
						// fprintf(fpw, "%s\n", seq[rootid].seq);
					// 	++ dupnum;
					// 	++ printnum;
					// }
					// for (uint32_t w = 0; w < rcdupno; ++w) {
						// fprintf(fpw, "%s\n", seq[rootid].seq);
					// 	++ dupnum;
					// 	++ printnum;
					// }

					if (0 == dupno && 0 == rcdupno) {
						bit_push(isdupbin, fpisdup, 0);
					} else {
						bit_push(isdupbin, fpisdup, 1);
						fpdup.write((char*)&dupno, sizeof(uint32_t));
						fpdup.write((char*)&rcdupno, sizeof(uint32_t));
					}
				}

				rootid = reads[rootid].getChildren();
			} else {
				dir = reads[rootid].isrc ? 1 : 0;
				bit_push(dirbin, fpdir, dir);

				// encode_v3(rootid, encstr);
				strcpy(encstr, seq[rootid].seq);
				
				// fprintf(fpw, "%s\n", seq[rootid].seq);
				// ++ printnum;

				fprintf(fpenstr, "%s\n", encstr);

				if (reads[rootid].shift != 0) {
					// bit for shift offset
					if (reads[rootid].shift > 0) {
						bit_push(isleftbin, fpisleft, 1);
					} else {
						bit_push(isleftbin, fpisleft, 0);
					}
				}

				{
					dupno = 0; rcdupno = 0;
					for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
						// fprintf(fpw, "%s\n", seq[rootid].seq);
						if (reads[rootid].dup[w].isrc) ++ rcdupno;
						else ++ dupno;
					}
					// for (uint32_t w = 0; w < dupno; ++w) {
					// 	fprintf(fpw, "%s\n", seq[rootid].seq);
					// 	++ dupnum;
					// 	++ printnum;
					// }
					// for (uint32_t w = 0; w < rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", seq[rootid].seq);
					// 	++ dupnum;
					// 	++ printnum;
					// }

					if (0 == dupno && 0 == rcdupno) {
						bit_push(isdupbin, fpisdup, 0);
					} else {
						bit_push(isdupbin, fpisdup, 1);
						fpdup.write((char*)&dupno, sizeof(uint32_t));
						fpdup.write((char*)&rcdupno, sizeof(uint32_t));
					}
				} 
			
				if (num >= cnt) {
					rootid = max_rid + 1;
					break;
				}
				uint32_t child = reads[rootid].getChildren();

				int back_step = 0;
				while (child > max_rid) { // no leaf node or all children are visited
					// back
					rootid = reads[rootid].prid;
					child = reads[rootid].getChildren();
					++back_step;
				}
				if (back_step > 0) {
					fprintf(fpenstr, "%d\n", back_step);
					// ++ numberofbackstep;
				}
				rootid = child;
			}
		}
	}
	
	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	if (isdupbin.n > 0) { 
		fpisdup.write((char*)&isdupbin.a, sizeof(uint8_t));
	}
	fpisdup.close();

	if (isleftbin.n > 0) { 
		fpisleft.write((char*)&isleftbin.a, sizeof(uint8_t));
	}
	fpisleft.close();

	// cout << "debugcnt: " << debugcnt << endl;
	// cout << "printnum: " << printnum << endl;
	// cout << "dupnum: " << dupnum << endl;
	// fclose(fpw);
	// fpdig.close();
	fpdup.close();

	// cout << "numberofbackstep: " << numberofbackstep << endl;

	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].crid.n > 0)
			kv_destroy(reads[rid].crid);
	}
	delete[] reads;
	fprintf(fpenstr, "-\n");

	int ncnt;

	// FILE *fpencsg = fopen("encsg.txt", "w");
	uint32_t sgcnt = 0, prerid = 0, nrid;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (!visited[rid]) {
			// fprintf(fpencsg, "%s\n", seq[rid].seq);
			fprintf(fpenstr, "%s\n", seq[rid].seq);
			// ++ sgcnt;
		}
	}

	// fclose(fpencsg);
	fclose(fpenstr);
	delete[] visited;
	
	// cout << "sgcnt: " << sgcnt << endl;
	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	vector<string> f1s, f2s;
	f1s.push_back(encodestrfn);
	f1s.push_back(dirbinfn);
	// f1s.push_back(parentbinfn);
	f1s.push_back(dupfn);
	f1s.push_back(isdupfn);
	f1s.push_back(isleftbinfn);

	compressFiles(f1s, f2s);

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc dir.bin.bsc "; 
	tarcmd += " dup.bin.bsc";
	tarcmd += " isdup.bin.bsc";
	tarcmd += " isleft.bin.bsc";

	int status = system(tarcmd.c_str());
	if (status < 0) {
		fprintf(stderr, "cmd: %s\t error: %s", tarcmd, strerror(errno));
	}
	if(WIFEXITED(status)) {
		    fprintf(stderr, "normal termination, exit status = %d\n", WEXITSTATUS(status)); //取得cmdstring执行结果
	}
	else if(WIFSIGNALED(status)) {
		fprintf(stderr, "abnormal termination,signal number =%d\n", WTERMSIG(status)); //如果cmdstring被信号中断，取得信号值
	} 
	else if(WIFSTOPPED(status)) {
		fprintf(stderr, "process stopped, signal number =%d\n", WSTOPSIG(status)); //如果cmdstring被信号暂停执行，取得信号值
	}
	cout << "Time of bsc files = " << stopwatch.stop() << std::endl;
}

// store singleton in a seperate file
void outputSingleDFS_best_version_v0() {
// void outputSingleDFS() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));
	char *encstr = (char*)alloca((L + 1) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	
	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;

	uint8bit_v isleftbin;
	isleftbin.n = isleftbin.a = 0;
	
	uint8bit_v isdupbin;
	isdupbin.n = isdupbin.a = 0;

	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string isleftbinfn = folder + "isleft.bin";
	std::ofstream fpisleft(isleftbinfn.c_str(), std::ios::binary);

	string isdupfn = folder + "isdup.bin";
	std::ofstream fpisdup(isdupfn.c_str(), std::ios::binary);

	string dupfn = folder + "dup.bin";
	std::ofstream fpdup(dupfn, std::ios::binary);

	// FILE *fpw = fopen("eenc.seq", "w");

	// uint32_t numberofbackstep = 0;
	uint32_t dupno = 0, rcdupno = 0;
	// uint32_t dupnum = 0, printnum = 0;
	// cout << "111" << endl;
	for (uint32_t rid = 0; rid < max_rid; ++rid) 
	if (!visited[rid] && reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
		// cout << "134904's parent " << reads[134904].prid << endl;
		// cout << "reads no. in this tree: " << rootnodevec[i].nodecnt << endl;
		queue<uint32_t> q;
		q.push(rid);
		uint32_t cnt = 0;
		while (!q.empty()) {
			uint32_t noderid = q.front();		
			q.pop();
			++cnt;
			visited[noderid] = true;
			// gridvec.push_back(noderid);
			for (size_t i = 0; i < reads[noderid].dn; ++i) {
				visited[reads[noderid].dup[i].id] = true;
			}
			for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
				q.push(reads[noderid].crid.a[i]);
			}
		}

		uint32_t rootid = rid;
		int dir;
		uint32_t dn;
		int preshiftdirection = 0; // 0 for left shift; 1 for right shift;

		uint32_t num = 0;
		while (rootid < max_rid) {
			++ num;
			dn = reads[rootid].dn;

			if (num == 1) {
				fprintf(fpenstr, "%s\n", seq[rootid].seq);
				// ++ printnum;

				// fprintf(fpw, "%s\n", seq[rootid].seq);
				{
					dupno = 0; rcdupno = 0;
					for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
						// fprintf(fpw, "%s\n", seq[rootid].seq);
						if (reads[rootid].dup[w].isrc) ++ rcdupno;
						else ++ dupno;
					}
					// for (uint32_t w = 0; w < dupno; ++w) {
						// fprintf(fpw, "%s\n", seq[rootid].seq);
					// 	++ dupnum;
					// 	++ printnum;
					// }
					// for (uint32_t w = 0; w < rcdupno; ++w) {
						// fprintf(fpw, "%s\n", seq[rootid].seq);
					// 	++ dupnum;
					// 	++ printnum;
					// }

					if (0 == dupno && 0 == rcdupno) {
						bit_push(isdupbin, fpisdup, 0);
					} else {
						bit_push(isdupbin, fpisdup, 1);
						fpdup.write((char*)&dupno, sizeof(uint32_t));
						fpdup.write((char*)&rcdupno, sizeof(uint32_t));
					}
				}

				rootid = reads[rootid].getChildren();
			} else {
				dir = reads[rootid].isrc ? 1 : 0;
				bit_push(dirbin, fpdir, dir);

				// encode_v3(rootid, encstr);
				strcpy(encstr, seq[rootid].seq);
				
				// fprintf(fpw, "%s\n", seq[rootid].seq);
				// ++ printnum;

				fprintf(fpenstr, "%s\n", encstr);

				if (reads[rootid].shift != 0) {
					// bit for shift offset
					if (reads[rootid].shift > 0) {
						bit_push(isleftbin, fpisleft, 1);
					} else {
						bit_push(isleftbin, fpisleft, 0);
					}
				}

				{
					dupno = 0; rcdupno = 0;
					for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
						// fprintf(fpw, "%s\n", seq[rootid].seq);
						if (reads[rootid].dup[w].isrc) ++ rcdupno;
						else ++ dupno;
					}
					// for (uint32_t w = 0; w < dupno; ++w) {
					// 	fprintf(fpw, "%s\n", seq[rootid].seq);
					// 	++ dupnum;
					// 	++ printnum;
					// }
					// for (uint32_t w = 0; w < rcdupno; ++w) {
					// 	fprintf(fpw, "%s\n", seq[rootid].seq);
					// 	++ dupnum;
					// 	++ printnum;
					// }

					if (0 == dupno && 0 == rcdupno) {
						bit_push(isdupbin, fpisdup, 0);
					} else {
						bit_push(isdupbin, fpisdup, 1);
						fpdup.write((char*)&dupno, sizeof(uint32_t));
						fpdup.write((char*)&rcdupno, sizeof(uint32_t));
					}
				} 
			
				if (num >= cnt) {
					rootid = max_rid + 1;
					break;
				}
				uint32_t child = reads[rootid].getChildren();

				int back_step = 0;
				while (child > max_rid) { // no leaf node or all children are visited
					// back
					rootid = reads[rootid].prid;
					child = reads[rootid].getChildren();
					++back_step;
				}
				if (back_step > 0) {
					fprintf(fpenstr, "%d\n", back_step);
					// ++ numberofbackstep;
				}
				rootid = child;
			}
		}
	}
	
	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	if (isdupbin.n > 0) { 
		fpisdup.write((char*)&isdupbin.a, sizeof(uint8_t));
	}
	fpisdup.close();

	if (isleftbin.n > 0) { 
		fpisleft.write((char*)&isleftbin.a, sizeof(uint8_t));
	}
	fpisleft.close();

	// cout << "debugcnt: " << debugcnt << endl;
	// cout << "printnum: " << printnum << endl;
	// cout << "dupnum: " << dupnum << endl;
	// fclose(fpw);
	// fpdig.close();
	fpdup.close();

	// cout << "numberofbackstep: " << numberofbackstep << endl;

	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].crid.n > 0)
			kv_destroy(reads[rid].crid);
	}
	delete[] reads;
	fprintf(fpenstr, "-\n");

	uint8bit_v singlebin;
	singlebin.a = singlebin.n = 0;

	string singletonfn = folder + "singleton.bin"; 
	std::ofstream singleOfs(singletonfn.c_str(), std::ios::binary);

	uint32_t sgcnt = 0, ncnt;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (!visited[rid]) {
			ncnt = 0;
			for (int i = 0; i < L; ++i) {
				if (seq[rid].seq[i] == 'N') ++ncnt;
			}
			if (ncnt > 0) {
				fprintf(fpenstr, "%s\n", seq[rid].seq);
			} else {
				for (int i = 0; i < L; ++i) {
					DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
				}
				++ sgcnt;
			}
		}
	}
	if (singlebin.n > 0) {
		singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
	}
	singleOfs.close();

	// fclose(fpencsg);
	fclose(fpenstr);
	delete[] visited;
	
	cout << "sgcnt: " << sgcnt << endl;
	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	vector<string> f1s, f2s;
	f1s.push_back(encodestrfn);
	f1s.push_back(dirbinfn);
	// f1s.push_back(parentbinfn);
	f1s.push_back(dupfn);
	f1s.push_back(isdupfn);
	f1s.push_back(isleftbinfn);
	f1s.push_back(singletonfn);

	compressFiles(f1s, f2s);

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc dir.bin.bsc "; 
	tarcmd += " dup.bin.bsc";
	tarcmd += " isdup.bin.bsc";
	tarcmd += " isleft.bin.bsc";
	tarcmd += " singleton.bin.bsc";

	int status = system(tarcmd.c_str());
	if (status < 0) {
		fprintf(stderr, "cmd: %s\t error: %s", tarcmd, strerror(errno));
	}
	if(WIFEXITED(status)) {
		    fprintf(stderr, "normal termination, exit status = %d\n", WEXITSTATUS(status)); //取得cmdstring执行结果
	}
	else if(WIFSIGNALED(status)) {
		fprintf(stderr, "abnormal termination,signal number =%d\n", WTERMSIG(status)); //如果cmdstring被信号中断，取得信号值
	} 
	else if(WIFSTOPPED(status)) {
		fprintf(stderr, "process stopped, signal number =%d\n", WSTOPSIG(status)); //如果cmdstring被信号暂停执行，取得信号值
	}
	cout << "Time of bsc files = " << stopwatch.stop() << std::endl;
}


vector<string> files;

void compressFile() {
	while(1) {
		uint32_t i = __sync_fetch_and_add(&rid_pthread, 1);
		if (i >= files.size()) break;

		size_t e = files[i].find(".lzma");
		if (e != std::string::npos) { // compressed by lzma
			mstcom::lzma::lzma_compress(files[i].substr(0, e).c_str(), files[i].c_str());
		} else {
			e = files[i].find(".bsc");
			if (e != std::string::npos) { // compressed by bsc
				int bbs = 64; // bsc block size
				if (files[i].find("parent.bin") != std::string::npos) bbs = 128;
				mstcom::bsc::BSC_compress(files[i].substr(0, e).c_str(), files[i].c_str(), bbs);
			}
		}
	}
}

string getExt(string file) {
	std::ifstream f0(file + ".lzma", std::ios::ate | std::ios::binary);
	if (f0.fail()) {
		std::cerr << "\nFile '" << file + ".lzma" << "' does not exist. Quit.";
		exit(1);
	}
	size_t N0 = f0.tellg();
	// cout << "size of " << file + ".lzma" << ": " << N0 << endl;

	std::ifstream f1(file + ".bsc", std::ios::ate | std::ios::binary);
	if (f1.fail()) {
		std::cerr << "\nFile '" << file + ".bsc" << "' does not exist. Quit.";
		exit(1);
	}
	size_t N1 = f1.tellg();
	// cout << "size of " << file + ".bsc" << ": " << N1 << endl;
	if (N0 < N1) return "lzma";
	return "bsc";
}

//f1s -- all files are compressed by only bsc; f2s -- all files are compressed by both bsc and lzma
void compressFiles(vector<string> &f1s, vector<string> &f2s) { 
	for (size_t i = 0; i < f2s.size(); ++i) {
		files.push_back(f2s[i] + ".lzma");
		files.push_back(f2s[i] + ".bsc");
	}
	for (size_t i = 0; i < f1s.size(); ++i) {
		files.push_back(f1s[i] + ".bsc");
	}
	for (size_t i = 0; i < files.size(); ++i) {
		// cout << files[i] << endl;
	}
	rid_pthread = 0;
	vector<thread> threadVec;
	for (int i = 0; i < nthreads; ++i) {
		threadVec.push_back(std::thread(compressFile));
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
}

// remove singleton reads by adding bit stream
// void outputSingleOrder_v1_work_well_adding_bit() {
void outputSingleOrder() {
	stopwatch.resume();

	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	uint32_t noderid;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (!visited[rid] && reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
			queue<uint32_t> q;
			q.push(rid);
			uint32_t cnt = 0;
			while (!q.empty()) {
				noderid = q.front();		
				q.pop();
				++cnt;
				visited[noderid] = true;
				for (size_t i = 1; i < reads[noderid].dn + 1; ++i) {
					visited[reads[noderid].dup[i].id] = true;
				}
				for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}
		}
	}

	uint8bit_v isrootbin;
	isrootbin.n = isrootbin.a = 0;

	string isrootbinfn = folder + "isroot.bin";
	std::ofstream fpisroot(isrootbinfn.c_str(), std::ios::binary);

	uint32_t *newid = new uint32_t[max_rid], newidnum = 0;
	// FILE *fprid2newid = fopen("rid2newid.comp.txt", "w");
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		int len = strlen(seq[rid].seq);
		if (len == L) { // root node or singleton
			int isroot = visited[rid] ? 1 : 0;
			bit_push(isrootbin, fpisroot, isroot);

			if (isroot == 1) { // root node
				newid[rid] = newidnum ++;
				// fprintf(fprid2newid, "%u %u\n", rid, newid[rid]);
			}
			// singleton is ignored
		} else 
		if (len > 0) {
		// if (seq[rid].seq[0] != '\0') {
			newid[rid] = newidnum ++;
			// fprintf(fprid2newid, "%u %u\n", rid, newid[rid]);
		}
	}
	// fclose(fprid2newid);
	if (isrootbin.n > 0) { 
		fpisroot.write((char*)&isrootbin.a, sizeof(uint8_t));
	}
	fpisroot.close();
	// cout << "newidnum: " << newidnum << endl;
	// exit(0);

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	
	std::ofstream fporder;
	string orderfn = folder + string("parent.bin");
	fporder.open(orderfn, std::ios::binary);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;

	uint8bit_v isleftbin;
	isleftbin.n = isleftbin.a = 0;
	
	uint8bit_v isdupbin;
	isdupbin.n = isdupbin.a = 0;

	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string isleftbinfn = folder + "isleft.bin";
	std::ofstream fpisleft(isleftbinfn.c_str(), std::ios::binary);

	// FILE *fporiid = fopen("oriid.txt", "w");
	// FILE *fppid = fopen("pid.comp.txt", "w");

	uint32_t dupno = 0, rcdupno = 0;

	int dir, isori;
	// FILE *fpencoded = fopen("encoded.txt", "w");
	/// current is OK.........................................
	// cout << "1111 \n";
	uint32_t outputcnt = 0;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		fprintf(fpenstr, "%s\n", seq[rid].seq);
		// fprintf(fpencoded, "%s", seq[rid].seq);

		if (strlen(seq[rid].seq) < L) {
			if (outputcnt < 10) {
				cout << newid[reads[rid].prid] << endl;
			}
			++ outputcnt;
			fporder.write((char*)&newid[reads[rid].prid], sizeof(uint32_t));
			// fprintf(fppid, "%u %u\n", reads[rid].prid, newid[reads[rid].prid]);

			dir = reads[rid].isrc ? 1 : 0;
			bit_push(dirbin, fpdir, dir);

			if (reads[rid].shift != 0) {
				// bit for shift offset
				if (reads[rid].shift > 0) {
					bit_push(isleftbin, fpisleft, 1);
				} else {
					bit_push(isleftbin, fpisleft, 0);
				}
			}
		} /*else {
			fprintf(fppid, "\n");
		}*/
		// fprintf(fpencoded, " %u", reads[rid].prid);
		// fprintf(fpencoded, "\n");
	}
	// fclose(fppid);
	
	// fclose(fpencoded);
	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	if (isleftbin.n > 0) { 
		fpisleft.write((char*)&isleftbin.a, sizeof(uint8_t));
	}
	fpisleft.close();

	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].crid.n > 0)
			kv_destroy(reads[rid].crid);
	}
	delete[] reads;
	fclose(fpenstr);

	fporder.close();

	delete[] visited;
	delete[] newid;

	vector<string> f1s, f2s;
	f1s.push_back(encodestrfn);
	f1s.push_back(dirbinfn);

	f1s.push_back(orderfn);
	f1s.push_back(isleftbinfn);
	f1s.push_back(isrootbinfn);

	compressFiles(f1s, f2s);

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc dir.bin.bsc ";

	tarcmd += " parent.bin.bsc";
	tarcmd += " isleft.bin.bsc";
	tarcmd += " isroot.bin.bsc";

	cout << tarcmd << endl;

	int status = system(tarcmd.c_str());
	if (status < 0) {
		fprintf(stderr, "cmd: %s\t error: %s", tarcmd, strerror(errno));
	}
	if(WIFEXITED(status)) {
		    fprintf(stderr, "normal termination, exit status = %d\n", WEXITSTATUS(status)); //取得cmdstring执行结果
	}
	else if(WIFSIGNALED(status)) {
		fprintf(stderr, "abnormal termination,signal number =%d\n", WTERMSIG(status)); //如果cmdstring被信号中断，取得信号值
	} 
	else if(WIFSTOPPED(status)) {
		fprintf(stderr, "process stopped, signal number =%d\n", WSTOPSIG(status)); //如果cmdstring被信号暂停执行，取得信号值
	}
	cout << "Time of bsc files = " << stopwatch.stop() << std::endl;
}

inline void build(uint32_t *&tr, uint32_t &M, uint32_t n) {
    for (M=1; M<=n+1; M<<=1);
    tr = new uint32_t[M + n + 1];
	memset(tr, 0, sizeof(uint32_t)*(M + n + 1));
}

inline void up(uint32_t *tr, uint32_t x) {
    tr[x] = tr[x<<1] + tr[x<<1|1];
}

inline void update(uint32_t *tr, uint32_t M, uint32_t x, int y) {
    for(tr[x+=M]+=y,x>>=1; x; x>>=1)
        up(tr, x);
}

inline uint32_t query(uint32_t *tr, uint32_t M, uint32_t s, uint32_t t) {
    uint32_t ans = 0;
    s = s+M-1;
    t = t+M+1;
    for(; s^t^1; s>>=1,t>>=1) {
        if(~s&1)ans += tr[s^1];
        if(t&1) ans += tr[t^1];
    }
    return ans;
}

void outputPEX() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");

	std::ofstream fpdist;
	string distfn = folder + string("dist.bin");
	fpdist.open(distfn, std::ios::binary);

	std::ofstream smfpdist;
	string smdistfn = folder + string("smdist.bin");
	smfpdist.open(smdistfn, std::ios::binary);

	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;

	uint8bit_v isleftbin;
	isleftbin.n = isleftbin.a = 0;
	
	uint8bit_v isdupbin;
	isdupbin.n = isdupbin.a = 0;

	uint8bit_v filebin;
	filebin.n = filebin.a = 0;

	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string isleftbinfn = folder + "isleft.bin";
	std::ofstream fpisleft(isleftbinfn.c_str(), std::ios::binary);

	string isdupfn = folder + "isdup.bin";
	std::ofstream fpisdup(isdupfn.c_str(), std::ios::binary);

	string dupfn = folder + "dup.bin";
	std::ofstream fpdup(dupfn, std::ios::binary);

	string filefn = folder + "file.bin";
	std::ofstream fpfile(filefn, std::ios::binary);
	
	char *en_str = (char*)alloca((L + 1) * sizeof(char));

	uint32_t *ids = new uint32_t[max_rid + 1];
	uint32_t idsid = 0, dn, curid;
	uint32_t *od = new uint32_t[max_rid + 3];
	int dir;

	uint32_t dupno = 0, rcdupno = 0;
	uint32_t dupnum = 0, printnum = 0;

	// FILE *fpid = fopen("idx.txt", "w");
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (!visited[rid] && reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
			queue<uint32_t> q;
			q.push(rid);
			uint32_t cnt = 0;
			while (!q.empty()) {
				uint32_t noderid = q.front();		
				q.pop();
				++cnt;
				visited[noderid] = true;
				for (size_t i = 0; i < reads[noderid].dn; ++i) {
					visited[reads[noderid].dup[i].id] = true;
				}
				for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}

			uint32_t rootid = rid;
			uint32_t num = 0;
			while (rootid < max_rid) {
				ids[idsid] = rootid;
				od[rootid] = idsid ++;

				// fprintf(fpid, "%lu\n", rootid);
				++ num;
				visited[rootid] = true;

				dn = reads[rootid].dn;

				if (num == 1) {
					fprintf(fpenstr, "%s\n", seq[rootid].seq);

					{
						dupno = 0; rcdupno = 0;
						for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
							// fprintf(fpw, "%s\n", seq[rootid].seq);
							if (reads[rootid].dup[w].isrc) ++ rcdupno;
							else ++ dupno;
						}

						for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
							if (!reads[rootid].dup[w].isrc) {
								curid = reads[rootid].dup[w].id;
								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						}

						for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
							if (reads[rootid].dup[w].isrc) {
								curid = reads[rootid].dup[w].id;
								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						}

						if (0 == dupno && 0 == rcdupno) {
							bit_push(isdupbin, fpisdup, 0);
						} else {
							bit_push(isdupbin, fpisdup, 1);
							fpdup.write((char*)&dupno, sizeof(uint32_t));
							fpdup.write((char*)&rcdupno, sizeof(uint32_t));
						}
					}

					rootid = reads[rootid].getChildren();
				} else {
					dir = reads[rootid].isrc ? 1 : 0;
					bit_push(dirbin, fpdir, dir);

					fprintf(fpenstr, "%s\n", seq[rootid].seq);

					if (reads[rootid].shift != 0) {
						// bit for shift offset
						if (reads[rootid].shift > 0) {
							bit_push(isleftbin, fpisleft, 1);
						} else {
							bit_push(isleftbin, fpisleft, 0);
						}
					}

					{
						dupno = 0; rcdupno = 0;
						for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
							// fprintf(fpw, "%s\n", seq[rootid].seq);
							if (reads[rootid].dup[w].isrc) ++ rcdupno;
							else ++ dupno;
						}

						for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
							if (!reads[rootid].dup[w].isrc) {
								curid = reads[rootid].dup[w].id;
								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						}

						for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
							if (reads[rootid].dup[w].isrc) {
								curid = reads[rootid].dup[w].id;
								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						}

						if (0 == dupno && 0 == rcdupno) {
							bit_push(isdupbin, fpisdup, 0);
						} else {
							bit_push(isdupbin, fpisdup, 1);
							fpdup.write((char*)&dupno, sizeof(uint32_t));
							fpdup.write((char*)&rcdupno, sizeof(uint32_t));
						}
					}

					if (num >= cnt) {
						rootid = max_rid + 1;
						break;
					}
					size_t child = reads[rootid].getChildren();
					int back_step = 0;
					while (child > max_rid) { // no leaf node or all children are visited
						// back
						rootid = reads[rootid].prid;
						child = reads[rootid].getChildren();
						++back_step;
					}
					if (back_step > 0) {
						fprintf(fpenstr, "%d\n", back_step);
					}
					rootid = child;
				}
			}
		}
	}

	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	if (isdupbin.n > 0) { 
		fpisdup.write((char*)&isdupbin.a, sizeof(uint8_t));
	}
	fpisdup.close();
	fpdup.close();

	if (isleftbin.n > 0) { 
		fpisleft.write((char*)&isleftbin.a, sizeof(uint8_t));
	}
	fpisleft.close();

	fprintf(fpenstr, "-\n");
	for (uint32_t rid = 0; rid < (max_rid); ++rid) {
		if (!visited[rid]) {
			ids[idsid] = rid;
			od[rid] = idsid ++;
			// fprintf(fpid, "%lu\n", rid);
			fprintf(fpenstr, "%s\n", seq[rid].seq);
		}
	}
	fclose(fpenstr);

	// calc dist
	// FILE *fp = fopen("dist.txt", "w");
	uint32_t *tr, M;
	build(tr, M, max_rid);
	memset(visited, false, sizeof(bool)*max_rid);

	uint32_t half = max_rid >> 1;
	uint32_t id = 0, num = 0, p, r, dist, orione;
	int32_t res;

	// cout << "before segment array\n";
	// in the segment array, (id+1, p+1)
	bool debug = false;
	int file;
	uint16_t disttemp;
	uint32_t small = 0;
	while (id < max_rid && num < half) {
		if (!visited[id]) {
			size_t dis = 0;
			if (ids[id] < half) { // from first file
				file = 0;
				p = od[half + ids[id]];

				r = query(tr, M, id + 1 + 1, p + 1);
				if (debug) cout << "r: " << r << endl;
				dist = p + 1 - (id + 1 + 1) + 1 - r;
				res = (int32_t)dist;
			} else { // from second file
				file = 1;
				p = od[ids[id] - half];

				r = query(tr, M, id + 1 + 1, p + 1);
				dist = p + 1 - (id + 1 + 1) + 1 - r;
				res = 0 - (int32_t)dist;
			}
			bit_push(filebin, fpfile, file);
			// bit_push(filebin, fileOfs, file);

			// fprintf(fp, "%d\n", res);

			if (dist < 65536) {
				++ small;
				file = 0;
				disttemp = dist;
				smfpdist.write((char *)&disttemp, sizeof(uint16_t));
			} else {
				file = 1;
				fpdist.write((char *)&dist, sizeof(uint32_t));
			}
			bit_push(filebin, fpfile, file);
			// fpdist.write((char *)&res, sizeof(int32_t));
			// fpdist.write((char *)&dist, sizeof(uint32_t));
			// bit_push(rangebin, rangeOfs, file);

			update(tr, M, p + 1, 1);
			visited[p] = true;

			++ num;
		}
		++ id;
	}
	
	if (filebin.n > 0) { 
		fpfile.write((char*)&filebin.a, sizeof(uint8_t));
	}
	fpfile.close();

	// cout << "end segment array\n";
	fpdist.close();
	smfpdist.close();

	delete[] visited;
	delete[] od;
	delete[] ids;

	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	vector<string> f1s, f2s;
	f1s.push_back(encodestrfn);
	f1s.push_back(dirbinfn);

	f1s.push_back(dupfn);
	f1s.push_back(isdupfn);
	f1s.push_back(isleftbinfn);
	f1s.push_back(filefn);
	f2s.push_back(distfn);
	f2s.push_back(smdistfn);

	compressFiles(f1s, f2s);

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc dir.bin.bsc ";

	tarcmd += " dup.bin.bsc";
	tarcmd += " isdup.bin.bsc";
	tarcmd += " isleft.bin.bsc";
	tarcmd += " file.bin.bsc";
	tarcmd += " dist.bin." + getExt(distfn);
	tarcmd += " smdist.bin." + getExt(smdistfn);

	cout << tarcmd << endl;

	system(tarcmd.c_str());
	cout << "Time of bsc files = " << stopwatch.stop() << std::endl;

}

// v0 work well
// void outputPEOrder_work_well_v0() {
void outputPEOrder() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");

	std::ofstream fporder;
	string orderfn = folder + string("order.bin");
	fporder.open(orderfn, std::ios::binary);

	std::ofstream fpdist;
	string distfn = folder + string("dist.bin");
	fpdist.open(distfn, std::ios::binary);

	std::ofstream smfpdist;
	string smdistfn = folder + string("smdist.bin");
	smfpdist.open(smdistfn, std::ios::binary);

	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;

	uint8bit_v isleftbin;
	isleftbin.n = isleftbin.a = 0;
	
	uint8bit_v isdupbin;
	isdupbin.n = isdupbin.a = 0;

	uint8bit_v filebin;
	filebin.n = filebin.a = 0;

	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	string isleftbinfn = folder + "isleft.bin";
	std::ofstream fpisleft(isleftbinfn.c_str(), std::ios::binary);

	string isdupfn = folder + "isdup.bin";
	std::ofstream fpisdup(isdupfn.c_str(), std::ios::binary);

	string dupfn = folder + "dup.bin";
	std::ofstream fpdup(dupfn, std::ios::binary);

	string filefn = folder + "file.bin";
	std::ofstream fpfile(filefn, std::ios::binary);
	
	char *en_str = (char*)alloca((L + 1) * sizeof(char));

	uint32_t *ids = new uint32_t[max_rid + 1];
	uint32_t idsid = 0, dn, curid;
	uint32_t *od = new uint32_t[max_rid + 3];
	int dir;

	uint32_t dupno = 0, rcdupno = 0;
	uint32_t dupnum = 0, printnum = 0;
	// FILE *fprootdegree = fopen("rootdegree.txt", "w");

	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (!visited[rid] && reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
			queue<uint32_t> q;
			q.push(rid);
			uint32_t cnt = 0;
			while (!q.empty()) {
				uint32_t noderid = q.front();		
				q.pop();
				++cnt;
				visited[noderid] = true;
				for (size_t i = 1; i < reads[noderid].dn + 1; ++i) {
					visited[reads[noderid].dup[i].id] = true;
				}
				for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}

			// fprintf(fprootdegree, "%u\n", reads[rid].crid.n);

			uint32_t rootid = rid;
			uint32_t num = 0;
			while (rootid < max_rid) {
				ids[idsid] = rootid;
				od[rootid] = idsid ++;

				// fprintf(fpid, "%lu\n", rootid);
				++ num;
				visited[rootid] = true;

				dn = reads[rootid].dn;

				if (num == 1) {
					fprintf(fpenstr, "%s\n", seq[rootid].seq);

					{
						dupno = 0; rcdupno = 0;
						for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
							// fprintf(fpw, "%s\n", seq[rootid].seq);
							if (reads[rootid].dup[w].isrc) ++ rcdupno;
							else ++ dupno;
						}

						for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
							if (!reads[rootid].dup[w].isrc) {
								curid = reads[rootid].dup[w].id;
								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						}

						for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
							if (reads[rootid].dup[w].isrc) {
								curid = reads[rootid].dup[w].id;
								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						}

						if (0 == dupno && 0 == rcdupno) {
							bit_push(isdupbin, fpisdup, 0);
						} else {
							bit_push(isdupbin, fpisdup, 1);
							fpdup.write((char*)&dupno, sizeof(uint32_t));
							fpdup.write((char*)&rcdupno, sizeof(uint32_t));
						}
					}

					rootid = reads[rootid].getChildren();
				} else {
					dir = reads[rootid].isrc ? 1 : 0;
					bit_push(dirbin, fpdir, dir);

					fprintf(fpenstr, "%s\n", seq[rootid].seq);

					if (reads[rootid].shift != 0) {
						// bit for shift offset
						if (reads[rootid].shift > 0) {
							bit_push(isleftbin, fpisleft, 1);
						} else {
							bit_push(isleftbin, fpisleft, 0);
						}
					}

					{
						dupno = 0; rcdupno = 0;
						for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
							// fprintf(fpw, "%s\n", seq[rootid].seq);
							if (reads[rootid].dup[w].isrc) ++ rcdupno;
							else ++ dupno;
						}

						for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
							if (!reads[rootid].dup[w].isrc) {
								curid = reads[rootid].dup[w].id;
								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						}

						for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
							if (reads[rootid].dup[w].isrc) {
								curid = reads[rootid].dup[w].id;
								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						}

						if (0 == dupno && 0 == rcdupno) {
							bit_push(isdupbin, fpisdup, 0);
						} else {
							bit_push(isdupbin, fpisdup, 1);
							fpdup.write((char*)&dupno, sizeof(uint32_t));
							fpdup.write((char*)&rcdupno, sizeof(uint32_t));
						}
					}

					if (num >= cnt) {
						rootid = max_rid + 1;
						break;
					}
					size_t child = reads[rootid].getChildren();
					int back_step = 0;
					while (child > max_rid) { // no leaf node or all children are visited
						// back
						rootid = reads[rootid].prid;
						child = reads[rootid].getChildren();
						++back_step;
					}
					if (back_step > 0) {
						fprintf(fpenstr, "%d\n", back_step);
					}
					rootid = child;
				}
			}
		}
	}
	// fclose(fprootdegree);

	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	if (isdupbin.n > 0) { 
		fpisdup.write((char*)&isdupbin.a, sizeof(uint8_t));
	}
	fpisdup.close();
	fpdup.close();

	if (isleftbin.n > 0) { 
		fpisleft.write((char*)&isleftbin.a, sizeof(uint8_t));
	}
	fpisleft.close();

	fprintf(fpenstr, "-\n");
	for (uint32_t rid = 0; rid < (max_rid); ++rid) {
		if (!visited[rid]) {
			ids[idsid] = rid;
			od[rid] = idsid ++;
			// fprintf(fpid, "%lu\n", rid);
			fprintf(fpenstr, "%s\n", seq[rid].seq);
		}
	}
	fclose(fpenstr);

	// calc dist
	// FILE *fp = fopen("dist.txt", "w");
	uint32_t *tr, M;
	build(tr, M, max_rid);
	memset(visited, false, sizeof(bool)*max_rid);

	uint32_t half = max_rid >> 1;
	uint32_t id = 0, num = 0, p, r, dist, orione;
	int32_t res;

	// FILE *fpdisttxt = fopen("dist.txt", "w");
	// cout << "before segment array\n";
	// in the segment array, (id+1, p+1)
	bool debug = false;
	int file;
	uint16_t disttemp;
	uint32_t small = 0;
	while (id < max_rid && num < half) {
		if (!visited[id]) {
			if (ids[id] < half) { // from first file
				file = 0;
				p = od[half + ids[id]];

				r = query(tr, M, id + 1 + 1, p + 1);
				if (debug) cout << "r: " << r << endl;
				dist = p + 1 - (id + 1 + 1) + 1 - r;
				res = (int32_t)dist;

				fporder.write((char *)&ids[id], sizeof(uint32_t));
			} else { // from second file
				file = 1;
				p = od[ids[id] - half];

				r = query(tr, M, id + 1 + 1, p + 1);
				dist = p + 1 - (id + 1 + 1) + 1 - r;
				res = 0 - (int32_t)dist;

				orione = ids[id] - half;
				fporder.write((char *)&orione, sizeof(uint32_t));
			}
			bit_push(filebin, fpfile, file);
			// bit_push(filebin, fileOfs, file);

			// fprintf(fp, "%d\n", res);
			-- dist;
			if (dist < 65536) {
				++ small;
				file = 0;
				disttemp = dist;
				smfpdist.write((char *)&disttemp, sizeof(uint16_t));
			} else {
				file = 1;
				fpdist.write((char *)&dist, sizeof(uint32_t));
			}
			bit_push(filebin, fpfile, file);
			// fprintf(fpdisttxt, "%u\n", dist);
			// fpdist.write((char *)&res, sizeof(int32_t));
			// bit_push(rangebin, rangeOfs, file);

			update(tr, M, p + 1, 1);
			visited[p] = true;

			++ num;
		}
		++ id;
	}
	// fclose(fpdisttxt);

	if (filebin.n > 0) { 
		fpfile.write((char*)&filebin.a, sizeof(uint8_t));
	}
	fpfile.close();

	cout << "end segment array\n";

	fpdist.close();
	smfpdist.close();
	fporder.close();

	delete[] visited;
	delete[] od;
	delete[] ids;

	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	vector<string> f1s, f2s;
	f1s.push_back(encodestrfn);
	f1s.push_back(dirbinfn);

	f1s.push_back(dupfn);
	f1s.push_back(isdupfn);
	f1s.push_back(isleftbinfn);
	f1s.push_back(filefn);
	f1s.push_back(orderfn);
	f2s.push_back(distfn);
	f2s.push_back(smdistfn);

	compressFiles(f1s, f2s);

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc dir.bin.bsc ";

	tarcmd += " dup.bin.bsc";
	tarcmd += " isdup.bin.bsc";
	tarcmd += " isleft.bin.bsc";
	tarcmd += " file.bin.bsc";
	tarcmd += " order.bin.bsc";
	tarcmd += " dist.bin." + getExt(distfn);
	tarcmd += " smdist.bin." + getExt(smdistfn);

	cout << tarcmd << endl;

	system(tarcmd.c_str());
	cout << "Time of bsc files = " << stopwatch.stop() << std::endl;

}
