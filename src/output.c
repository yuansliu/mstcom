#include "mstcom.h"

bool cmp(const ROOTNODE_t &a, const ROOTNODE_t &b) {
	return a.nodecnt < b.nodecnt;
}

void outputORI() {
	stopwatch.resume();
	// char *str = (char*)alloca((L + 1) * sizeof(char));
	char *rcstr = (char*)alloca((L + 1) * sizeof(char));

	FILE *fpcnt = fopen("count.txt", "w");
	FILE *fptrees = fopen("trees.txt", "w");
	FILE *fpenstr = fopen("encodestr.txt", "w");
	FILE *fpenstrcopy = fopen("encodestr.txt.copy", "w");
	FILE *fpenstrcopycopy = fopen("encodestr.txt.copy.copy", "w");
	FILE *fproot = fopen("rootstr.txt", "w");
	char name[100]; 
	// uint8bit_v rootbin;
	// rootbin.a = rootbin.n = 0;
	// sprintf(name, "rootstr.bin");
	// std::ofstream rootOfs(name, std::ios::binary);

	FILE *fporder = NULL;
	if (isorder) {
		string orderfn = folder + string("order.bin");
		// string orderfn = string("order.bin");  // for test
		fporder = fopen(orderfn.c_str(), "wb");
	}

	vector<ROOTNODE_t> rootnodevec;
	for (size_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid && reads[rid].crid.n > 0) { //is the root node && not a leaf node == not a singleton reads
			queue<size_t> q;
			q.push(rid);
			size_t cnt = 0;
			while (!q.empty()) {
				size_t noderid = q.front();		
				q.pop();
				++cnt;
				for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}
			rootnodevec.push_back(ROOTNODE_t(rid, cnt));
		}
	}
	sort(rootnodevec.begin(), rootnodevec.end(), cmp);

	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	// size_t outputnum = 0;

	uint8bit_v dirbin;
	// uint8bit_v rootbin, dirbin;
	// rootbin.n = rootbin.a = 0;
	dirbin.n = dirbin.a = 0;
	sprintf(name, "dir.bin");
	std::ofstream fpdir(name, std::ios::binary);

	size_t ss = sizeof(size_t);
	char *en_str = (char*)alloca((L + 1) * sizeof(char));

	// int max_shift = 0;

	for (size_t i = 0; i < rootnodevec.size(); ++i) {
			size_t cnt = rootnodevec[i].nodecnt;
			fprintf(fpcnt, "%lu\n", cnt);

			size_t rootid = rootnodevec[i].rid;
			size_t num = 0;
			while (rootid < max_rid) {
				// fprintf(stderr, "%lu\n", rootid);
				++ num;
				visited[rootid] = true;

				if (fporder) fwrite(&rootid, ss, 1, fporder);

				if (num == 1) {
					// fprintf(fproot, "%s\n", reads[rootid].str.c_str());
					fprintf(fproot, "%s\n", seq[rootid].seq);
					// for (int j = 0; j < L; ++j) {
					// 	DNA_push(rootbin, rootOfs, seq_nt4_table[(uint8_t)seq[rootid].seq[j]]);
					// }
					rootid = reads[rootid].getChildren();
				} else {
					int dir = 0;
					fprintf(fpenstrcopy, "%s\n", seq[rootid].seq);
					// if (abs(reads[rootid].shift) > max_shift) {
					// 	max_shift = abs(reads[rootid].shift);
					// }
					if (reads[rootid].isrc) {
						dir = 1;
						strcpy(rcstr, seq[rootid].seq);
						reverseComplement(rcstr);
						encode(seq[reads[rootid].prid].seq, rcstr, reads[rootid].shift, en_str);
						fprintf(fpenstrcopycopy, "%s\n%s\n%d\n", seq[reads[rootid].prid].seq, rcstr, reads[rootid].shift);
						fprintf(fpenstrcopycopy, "%s\n", en_str);
						fprintf(fpenstr, "%s\n", en_str);
					} else {
						encode(seq[reads[rootid].prid].seq, seq[rootid].seq, reads[rootid].shift, en_str);
						fprintf(fpenstrcopycopy, "%s\n%s\n%d\n", seq[reads[rootid].prid].seq, seq[rootid].seq, reads[rootid].shift);
						fprintf(fpenstrcopycopy, "%s\n", en_str);
						fprintf(fpenstr, "%s\n", en_str);
					}

					bit_push(dirbin, fpdir, dir);

					if (num >= cnt) {
						rootid = max_rid + 1;
						break;
					}
					size_t child = reads[rootid].getChildren();
					// fprintf(stderr, "child: %lu\n", child);
					int back_step = 0;
					while (child > max_rid) { // no leaf node or all children are visited
						// back
						rootid = reads[rootid].prid;
						// fprintf(stderr, "rootid: %lu\n", rootid);
						child = reads[rootid].getChildren();
						++back_step;
						// ++num;
						// if (num > 20) exit(0);
					}
					// fprintf(stderr, "child: %lu\n", child);
					if (back_step > 0) {
						// fprintf(stderr, "--%d\n", back_step);
						fprintf(fpenstr, "%d\n", back_step);
					}
					rootid = child;
					// if (num > 20) exit(0);
				}
			}
			// exit(0);
			fprintf(fpenstr, "-\n");
		// } 
	}

	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	fclose(fproot);
	// if (rootbin.n > 0) {
	// 	rootOfs.write((char*)&rootbin.a, sizeof(uint8_t));
	// }
	// rootOfs.close();
	
	fclose(fpcnt);
	fclose(fptrees);
	fclose(fpenstr);

	// fprintf(stderr, "max_shift: %d\n", max_shift);
	fclose(fpenstrcopy);
	fclose(fpenstrcopycopy);

	// size_v nrid;
	// if (isquality) {
	// 	kv_init(nrid);
	// }

	uint8bit_v singlebin;
	singlebin.a = singlebin.n = 0;

	sprintf(name, "singleton.bin");
	std::ofstream singleOfs(name, std::ios::binary);

	FILE *fpN = fopen("readsN.txt", "w");
	int ncnt;

	for (size_t rid = 0; rid < max_rid; ++rid) {
		if (!visited[rid]) {
			ncnt = 0;
			for (int i = 0; i < L; ++i) {
				if (seq[rid].seq[i] == 'N') ++ncnt;
			}
			if (ncnt > 0) {
				fprintf(fpN, "%s\n", seq[rid].seq);
				// if (isquality) kv_push(size_t, nrid, rid);
				if (fporder) fwrite(&rid, ss, 1, fporder);
			} else {
				// fprintf(fpsingle, "%s\n", seq[rid].seq);
				for (int i = 0; i < L; ++i) {
					DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
				}
			}
		}
	}
	fclose(fpN);
	if (singlebin.n > 0) {
		singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
	}
	singleOfs.close();

	if (isorder) {
		fclose(fporder);
	}

	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();
}

bool cmpduprc0(const dup_t &a, const dup_t &b) {
	return a.isrc < b.isrc;
}

bool cmpduprc1(const dup_t &a, const dup_t &b) {
	return a.isrc > b.isrc;
}

void outputSingleOri() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	string rootstrfn = folder + "rootstr.txt";
	FILE *fproot = fopen(rootstrfn.c_str(), "w");

	std::ofstream fporder;
	string orderfn;
	if (isorder) {
		orderfn = folder + string("order.bin");
		// string orderfn = string("order.bin");  // for test
		// fporder = fopen(orderfn.c_str(), "wb");
		fporder.open(orderfn, std::ios::binary);
	}

	// cout << "begin outputSingle() \n";

	vector<ROOTNODE_t> rootnodevec;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
			queue<uint32_t> q;
			q.push(rid);
			uint32_t cnt = 0;
			while (!q.empty()) {
				uint32_t noderid = q.front();		
				q.pop();
				++cnt;
				for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}
			rootnodevec.push_back(ROOTNODE_t(rid, cnt));
		}
	}
	sort(rootnodevec.begin(), rootnodevec.end(), cmp);

	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;

	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	size_t ss = sizeof(size_t);
	char *en_str = (char*)alloca((L + 1) * sizeof(char));
	int dir;
	uint32_t dn, difid;
	// FILE *fpencoded = fopen("encoded.txt", "w");

	// cout << "1111 \n";
	if (isorder) {
		for (uint32_t i = 0; i < rootnodevec.size(); ++i) {

			/*if (isorder) {
				// sort the reads 
				uint32_t rootid = rootnodevec[i].rid, crid;				
				queue<uint32_t> q;
				q.push(rootid);
				while (!q.empty()) {
					rootid = q.front(); 
					q.pop();
					for (uint32_t i = 0; i < reads[rootid].crid.n; ++i) {
						crid = reads[rootid].crid.a[i];
						if (reads[crid].dn > 0) {

						}
					}

					dn = reads[rootid].dn;
					for (uint32_t i = 0; i < reads[noderid].crid.n; ++i) {
						q.push(reads[noderid].crid.a[i]);
					}
				}
			}*/

			bool tdebug = false;
			// tdebug = true;
			// cout << "i: " << i << endl;
			// if (i == 57185) tdebug = true;

			size_t cnt = rootnodevec[i].nodecnt;

			uint32_t rootid = rootnodevec[i].rid;

			if (rootid == 57185) tdebug = true;
			if (tdebug) {
				cout << "cnt: " << cnt << endl;
				cout << "rootid: " << rootid << endl;
			}
			uint32_t num = 0;
			while (rootid < max_rid) {
				++ num;
				visited[rootid] = true;

				tdebug = false;
				if (rootid == 464744) tdebug = true;

				if (tdebug) cout << "rootid: " << rootid << endl;

				dn = reads[rootid].dn;

				if (dn > 0) {
					// reads[rootid].dup[dn] = dup_t(rootid, reads[rootid].isrc);
					reads[rootid].dup[0].isrc = reads[rootid].isrc;
					// if (dn > 5) {
					// 	for (uint32_t w = 0; w < reads[rootid].dn + 1; ++w) {
					// 		cout << reads[rootid].dup[w].id << endl;
					// 	}
					// 	cout << "rootid: " << rootid << endl;
					// 	exit(0);
					// }
					// sort(reads[rootid].dup, reads[rootid].dup + dn + 1, cmpdupid);
				}

				// if (rootid == 0) {
				if (tdebug) {
					cout << seq[rootid].seq << endl;
					cout << reads[rootid].dn << endl;
					cout << "dn: " << dn << endl;
					cout << "isrc: " << reads[rootid].isrc << endl;
					// cout << "num: " << num << endl;
					// cout << rootnodevec[i].rid << endl;
				}

				if (num == 1) {
					if (dn > 0) {
						fprintf(fproot, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
						// fprintf(stdout, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
						// dir = 0; // 0 not changed
						// if (reads[rootid].isrc) dir = 1; // changed
						// bit_push(dirbin, fpdir, dir);

						for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
							dir = 0;
							if (reads[rootid].dup[w].isrc) dir = 1;
							bit_push(dirbin, fpdir, dir);
							
							visited[reads[rootid].dup[w].id] = true;
						}

						fporder.write((char*)&reads[rootid].dup[0].id, sizeof(uint32_t));
						for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
							difid = reads[rootid].dup[w].id - reads[rootid].dup[w-1].id;
							fporder.write((char*)&difid, sizeof(uint32_t));
							// if (reads[rootid].dup[w].id == 464744) cout << "HERERERERE1111" << endl;
						}
					} else {
						fprintf(fproot, "%s\n", seq[rootid].seq);
						fporder.write((char*)&rootid, sizeof(uint32_t));
					}

					if (tdebug) {
						fprintf(stderr, "---root---\n");
						fprintf(stderr, "%s\n", seq[rootid].seq);
						fprintf(stderr, "------\n");
					}

					rootid = reads[rootid].getChildren();
				} else {					
					if (dn > 0) {
						fprintf(fpenstr, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
						// fprintf(stdout, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);

						fporder.write((char*)&reads[rootid].dup[0].id, sizeof(uint32_t));
						for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
							difid = reads[rootid].dup[w].id - reads[rootid].dup[w-1].id;
							fporder.write((char*)&difid, sizeof(uint32_t));
						}

						// dir = 0; // 
						// if (reads[rootid].isrc) dir = 1; // 
						// bit_push(dirbin, fpdir, dir);
						for (uint32_t w = 0; w < reads[rootid].dn + 1; ++w) {
							dir = 0;
							if (reads[rootid].dup[w].isrc) dir = 1;
							bit_push(dirbin, fpdir, dir);
							
							// if (reads[rootid].dup[w].id == 20371) cout << "HERERERERE00000" << endl;
							visited[reads[rootid].dup[w].id] = true;
						}
					} else {
						dir = 0;
						if (reads[rootid].isrc) dir = 1;
						bit_push(dirbin, fpdir, dir);

						fprintf(fpenstr, "%s\n", seq[rootid].seq);
						fporder.write((char*)&rootid, sizeof(uint32_t));
					}
					// fprintf(fpencoded, "%s\n", seq[rootid].seq);

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
					}
					rootid = child;
				}

				// if (rootid == 0) exit(0);
			}
			fprintf(fpenstr, "-\n");
		}
	} else {
		for (uint32_t i = 0; i < rootnodevec.size(); ++i) {

			bool tdebug = false;
			// tdebug = true;
			// cout << "i: " << i << endl;
			// if (i == 419719) tdebug = true;

			size_t cnt = rootnodevec[i].nodecnt;

			uint32_t rootid = rootnodevec[i].rid;

			// if (rootid == 211517) tdebug = true;
			if (tdebug) {
				cout << "cnt: " << cnt << endl;
				cout << "rootid: " << rootid << endl;
			}
			uint32_t num = 0;
			while (rootid < max_rid) {
				++ num;
				visited[rootid] = true;

				if (tdebug) cout << "rootid: " << rootid << endl;

				dn = reads[rootid].dn;

				if (dn > 0) {
					// reads[rootid].dup[dn] = dup_t(rootid, reads[rootid].isrc);
					if (reads[rootid].isrc) { // sort 1 to 0
						sort(reads[rootid].dup, reads[rootid].dup + dn, cmpduprc1);
					} else { //small value 0 first
						sort(reads[rootid].dup, reads[rootid].dup + dn, cmpduprc0);
					}
				}

				// if (rootid == 0) {
				if (tdebug) {
					cout << seq[rootid].seq << endl;
					cout << reads[rootid].dn << endl;
					// cout << "num: " << num << endl;
					// cout << rootnodevec[i].rid << endl;
				}

				if (num == 1) {
					if (dn > 0) {
						fprintf(fproot, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
						// fprintf(stdout, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
						for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
							dir = 0;
							if (reads[rootid].dup[w].isrc) dir = 1;
							bit_push(dirbin, fpdir, dir);
							visited[reads[rootid].dup[w].id] = true;
						}
					} else {
						fprintf(fproot, "%s\n", seq[rootid].seq);
						if (isorder) {
							fporder.write((char*)&rootid, sizeof(uint32_t));
						}
					}

					if (tdebug) {
						fprintf(stderr, "---root---\n");
						fprintf(stderr, "%s\n", seq[rootid].seq);
						fprintf(stderr, "------\n");
					}

					rootid = reads[rootid].getChildren();
				} else {
					dir = 0;
					if (reads[rootid].isrc) dir = 1; // 
					bit_push(dirbin, fpdir, dir);
					
					if (dn > 0) {
						fprintf(fpenstr, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
						// fprintf(stdout, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);

						for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
							dir = 0;
							if (reads[rootid].dup[w].isrc) dir = 1;
							bit_push(dirbin, fpdir, dir);
							
							visited[reads[rootid].dup[w].id] = true;
						}
					} else {
						// dir = 0;
						// if (reads[rootid].isrc) dir = 1;
						// bit_push(dirbin, fpdir, dir);

						fprintf(fpenstr, "%s\n", seq[rootid].seq);
					}
					// fprintf(fpencoded, "%s\n", seq[rootid].seq);

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
					}
					rootid = child;
				}

				// if (rootid == 0) exit(0);
			}
			fprintf(fpenstr, "-\n");
		}
	}
	
	// cout << "2222111 \n";

	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	// cout << "2222111222 \n";
	// delete
	rootnodevec.clear();
	// cout << "2222111222333 \n";
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		kv_destroy(reads[rid].crid);
	}
	// cout << "2222111222333444 \n";
	delete[] reads;
	// cout << "2222111222333444555 \n";

	fclose(fproot);
	// cout << "2222111222333444555666 \n";
	fclose(fpenstr);
	// cout << "2222111222333444555666777 \n";

	uint8bit_v singlebin;
	singlebin.a = singlebin.n = 0;

	// sprintf(name, "singleton.bin");
	string singletonfn = folder + "singleton.bin"; 
	std::ofstream singleOfs(singletonfn.c_str(), std::ios::binary);

	string readsnfn = folder + "readsN.txt";
	FILE *fpN = fopen(readsnfn.c_str(), "w");
	int ncnt;

	// cout << "3333 \n";

	// FILE *fpsg = fopen("single.txt", "w");

	uint32_t sgcnt = 0, prerid = 0, nrid;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (!visited[rid]) {
			ncnt = 0;
			for (int i = 0; i < L; ++i) {
				if (seq[rid].seq[i] == 'N') ++ncnt;
			}
			if (ncnt > 0) {
				fprintf(fpN, "%s\n", seq[rid].seq);
				if (isorder) {
					nrid = rid - prerid;
					fporder.write((char*)&nrid, sizeof(uint32_t));
					prerid = rid;
				}
			} else {
				++ sgcnt;
				// fprintf(fpsg, "%s\n", seq[rid].seq);
				for (int i = 0; i < L; ++i) {
					DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
				}
			}
		}
	}
	// fclose(fpsg);

	fclose(fpN);
	if (singlebin.n > 0) {
		singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
	}
	singleOfs.close();

	delete[] visited;
	
	cout << "sgcnt: " << sgcnt << endl;
	if (isorder) {
		// fclose(fporder);
		fporder.close();
	}

	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	mstcom::bsc::BSC_compress(encodestrfn.c_str(), (encodestrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(rootstrfn.c_str(), (rootstrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(singletonfn.c_str(), (singletonfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(dirbinfn.c_str(), (dirbinfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(readsnfn.c_str(), (readsnfn+".bsc").c_str(), 64);

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc rootstr.txt.bsc singleton.bin.bsc dir.bin.bsc readsN.txt.bsc";

	if (isorder) {
		// mstcom::bsc::BSC_compress(orderfn.c_str(), (orderfn+".bsc").c_str(), 64);
		mstcom::lzma::lzma_compress(orderfn.c_str(), (orderfn+".lzma").c_str());
		tarcmd += " order.bin.lzma";
	}

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

void outputSingle_() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	string rootstrfn = folder + "rootstr.txt";
	FILE *fproot = fopen(rootstrfn.c_str(), "w");

	std::ofstream fporder;
	string orderfn;
	if (isorder) {
		orderfn = folder + string("order.bin");
		// string orderfn = string("order.bin");  // for test
		// fporder = fopen(orderfn.c_str(), "wb");
		fporder.open(orderfn, std::ios::binary);
	}

	// cout << "begin outputSingle() \n";

	// vector<ROOTNODE_t> rootnodevec;
	// for (uint32_t rid = 0; rid < max_rid; ++rid) {
	// 	if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
	// 		queue<uint32_t> q;
	// 		q.push(rid);
	// 		uint32_t cnt = 0;
	// 		while (!q.empty()) {
	// 			uint32_t noderid = q.front();		
	// 			q.pop();
	// 			++cnt;
	// 			for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
	// 				q.push(reads[noderid].crid.a[i]);
	// 			}
	// 		}
	// 		rootnodevec.push_back(ROOTNODE_t(rid, cnt));
	// 	}
	// }
	// sort(rootnodevec.begin(), rootnodevec.end(), cmp);

	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;

	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	size_t ss = sizeof(size_t);
	char *en_str = (char*)alloca((L + 1) * sizeof(char));
	int dir;
	uint32_t dn, difid;
	// FILE *fpencoded = fopen("encoded.txt", "w");
	/// current is OK.........................................
	// cout << "1111 \n";
	if (isorder) {
		uint32_t prerootrid = 0;
		for (uint32_t rid = 0; rid < max_rid; ++rid) 
		if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
		// for (uint32_t i = 0; i < rootnodevec.size(); ++i) {
			queue<uint32_t> q;
			q.push(rid);
			uint32_t cnt = 0;
			while (!q.empty()) {
				uint32_t noderid = q.front();		
				q.pop();
				++cnt;
				// sort reads[rnoderid].crid
				sort(reads[noderid].crid.a, reads[noderid].crid.a + reads[noderid].crid.n);
				for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}

			bool tdebug = false;
			// tdebug = true;
			// cout << "i: " << i << endl;
			// if (i == 57185) tdebug = true;
			// size_t cnt = rootnodevec[i].nodecnt;

			uint32_t rootid = rid;

			// if (rootid == 57185) tdebug = true;
			if (tdebug) {
				cout << "cnt: " << cnt << endl;
				cout << "rootid: " << rootid << endl;
			}
			uint32_t num = 0;
			while (rootid < max_rid) {
				++ num;
				visited[rootid] = true;

				tdebug = false;
				// if (rootid == 464744) tdebug = true;

				if (tdebug) cout << "rootid: " << rootid << endl;

				dn = reads[rootid].dn;

				if (dn > 0) {
					// reads[rootid].dup[dn] = dup_t(rootid, reads[rootid].isrc);
					reads[rootid].dup[0].isrc = reads[rootid].isrc;
				}

				// if (rootid == 0) {
				if (tdebug) {
					cout << seq[rootid].seq << endl;
					cout << reads[rootid].dn << endl;
					cout << "dn: " << dn << endl;
					cout << "isrc: " << reads[rootid].isrc << endl;
					// cout << "num: " << num << endl;
					// cout << rootnodevec[i].rid << endl;
				}

				if (num == 1) {
					difid = rootid - prerootrid;
					// difid = rootid;
					fporder.write((char*)&difid, sizeof(uint32_t));

					prerootrid = rootid;

					if (dn > 0) {
						if (dn == 1)
							fprintf(fproot, "%s$\n", seq[rootid].seq);
						else 
							fprintf(fproot, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);

						for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
							dir = 0;
							if (reads[rootid].dup[w].isrc) dir = 1;
							bit_push(dirbin, fpdir, dir);
							
							difid = reads[rootid].dup[w].id - reads[rootid].dup[w-1].id;
							fporder.write((char*)&difid, sizeof(uint32_t));

							visited[reads[rootid].dup[w].id] = true;
						}
						// fporder.write((char*)&reads[rootid].dup[0].id, sizeof(uint32_t));
						// for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
							// if (reads[rootid].dup[w].id == 464744) cout << "HERERERERE1111" << endl;
						// }
					} else {
						fprintf(fproot, "%s\n", seq[rootid].seq);
					}

					if (tdebug) {
						fprintf(stderr, "---root---\n");
						fprintf(stderr, "%s\n", seq[rootid].seq);
						fprintf(stderr, "------\n");
					}

					rootid = reads[rootid].getChildren();
				} else {					
					dir = 0;
					if (reads[rootid].isrc) dir = 1;
					bit_push(dirbin, fpdir, dir);
					difid = rootid - reads[reads[rootid].prid].prechildrid;
					fporder.write((char*)&difid, sizeof(uint32_t));

					reads[reads[rootid].prid].prechildrid = rootid;

					if (dn > 0) {
						if (dn == 1)
							fprintf(fpenstr, "%s$\n", seq[rootid].seq);
						else 
							fprintf(fpenstr, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
						// fprintf(stdout, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
						for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
							dir = 0;
							if (reads[rootid].dup[w].isrc) dir = 1;
							bit_push(dirbin, fpdir, dir);
							// if (reads[rootid].dup[w].id == 20371) cout << "HERERERERE00000" << endl;
							difid = reads[rootid].dup[w].id - reads[rootid].dup[w-1].id;
							fporder.write((char*)&difid, sizeof(uint32_t));

							visited[reads[rootid].dup[w].id] = true;
						}
					} else {
						fprintf(fpenstr, "%s\n", seq[rootid].seq);
					}
					// fprintf(fpencoded, "%s\n", seq[rootid].seq);

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
					}
					rootid = child;
				}

				// if (rootid == 0) exit(0);
			}
			fprintf(fpenstr, "-\n");
		}
	} else {
		// for (uint32_t i = 0; i < rootnodevec.size(); ++i) {
		for (uint32_t rid = 0; rid < max_rid; ++rid) 
		if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
		// for (uint32_t i = 0; i < rootnodevec.size(); ++i) {
			queue<uint32_t> q;
			q.push(rid);
			uint32_t cnt = 0;
			while (!q.empty()) {
				uint32_t noderid = q.front();		
				q.pop();
				++cnt;
				for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}

			bool tdebug = false;
			// tdebug = true;
			// cout << "i: " << i << endl;
			// if (i == 419719) tdebug = true;

			// size_t cnt = rootnodevec[i].nodecnt;

			uint32_t rootid = rid;

			// if (rootid == 211517) tdebug = true;
			if (tdebug) {
				cout << "cnt: " << cnt << endl;
				cout << "rootid: " << rootid << endl;
			}
			uint32_t num = 0;
			while (rootid < max_rid) {
				++ num;
				visited[rootid] = true;

				if (tdebug) cout << "rootid: " << rootid << endl;

				dn = reads[rootid].dn;

				if (dn > 0) {
					// reads[rootid].dup[dn] = dup_t(rootid, reads[rootid].isrc);
					if (reads[rootid].isrc) { // sort 1 to 0
						sort(reads[rootid].dup, reads[rootid].dup + dn, cmpduprc1);
					} else { //small value 0 first
						sort(reads[rootid].dup, reads[rootid].dup + dn, cmpduprc0);
					}
				}

				// if (rootid == 0) {
				if (tdebug) {
					cout << seq[rootid].seq << endl;
					cout << reads[rootid].dn << endl;
					// cout << "num: " << num << endl;
					// cout << rootnodevec[i].rid << endl;
				}

				if (num == 1) {
					if (dn > 0) {
						if (dn == 1)
							fprintf(fproot, "%s$\n", seq[rootid].seq);
						else 
							fprintf(fproot, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
						// fprintf(stdout, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
						for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
							dir = 0;
							if (reads[rootid].dup[w].isrc) dir = 1;
							bit_push(dirbin, fpdir, dir);
							visited[reads[rootid].dup[w].id] = true;
						}
					} else {
						fprintf(fproot, "%s\n", seq[rootid].seq);
						if (isorder) {
							fporder.write((char*)&rootid, sizeof(uint32_t));
						}
					}

					if (tdebug) {
						fprintf(stderr, "---root---\n");
						fprintf(stderr, "%s\n", seq[rootid].seq);
						fprintf(stderr, "------\n");
					}

					rootid = reads[rootid].getChildren();
				} else {
					dir = 0;
					if (reads[rootid].isrc) dir = 1; // 
					bit_push(dirbin, fpdir, dir);
					
					if (dn > 0) {
						if (dn == 1)
							fprintf(fpenstr, "%s$\n", seq[rootid].seq);
						else 
							fprintf(fpenstr, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
						// fprintf(stdout, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);

						for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
							dir = 0;
							if (reads[rootid].dup[w].isrc) dir = 1;
							bit_push(dirbin, fpdir, dir);
							
							visited[reads[rootid].dup[w].id] = true;
						}
					} else {
						// dir = 0;
						// if (reads[rootid].isrc) dir = 1;
						// bit_push(dirbin, fpdir, dir);

						fprintf(fpenstr, "%s\n", seq[rootid].seq);
					}
					// fprintf(fpencoded, "%s\n", seq[rootid].seq);

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
					}
					rootid = child;
				}

				// if (rootid == 0) exit(0);
			}
			fprintf(fpenstr, "-\n");
		}
	}
	
	// cout << "2222111 \n";

	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	// cout << "2222111222 \n";
	// delete
	// rootnodevec.clear();
	// cout << "2222111222333 \n";
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		kv_destroy(reads[rid].crid);
	}
	// cout << "2222111222333444 \n";
	delete[] reads;
	// cout << "2222111222333444555 \n";

	// cout << "2222111222333444555666 \n";
	fclose(fpenstr);
	// cout << "2222111222333444555666777 \n";

	uint8bit_v singlebin;
	singlebin.a = singlebin.n = 0;

	// sprintf(name, "singleton.bin");
	string singletonfn = folder + "singleton.bin"; 
	// FILE *fpsg = fopen(singletonfn.c_str(), "w");
	std::ofstream singleOfs(singletonfn.c_str(), std::ios::binary);

	// string readsnfn = folder + "readsN.txt";
	// FILE *fpN = fopen(readsnfn.c_str(), "w");
	int ncnt;

	// cout << "3333 \n";

	// string sgfn = folder + "single.txt"; 
	// FILE *fpsg = fopen(sgfn.c_str(), "w");

	fprintf(fproot, "-\n");

	uint32_t sgcnt = 0, prerid = 0, nrid;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (!visited[rid]) {
			ncnt = 0;
			for (int i = 0; i < L; ++i) {
				if (seq[rid].seq[i] == 'N') ++ncnt;
			}
			if (ncnt > 0) {
				fprintf(fproot, "%s\n", seq[rid].seq);
				if (isorder) {
					nrid = rid - prerid;
					fporder.write((char*)&nrid, sizeof(uint32_t));
					prerid = rid;
				}
			} else {
				++ sgcnt;
				for (int i = 0; i < L; ++i) {
					DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
				}
			}
		}
	}

	// fclose(fpN);
	fclose(fproot);
	if (singlebin.n > 0) {
		singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
	}
	singleOfs.close();

	delete[] visited;
	
	cout << "sgcnt: " << sgcnt << endl;
	if (isorder) {
		// fclose(fporder);
		fporder.close();
	}

	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	mstcom::bsc::BSC_compress(encodestrfn.c_str(), (encodestrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(rootstrfn.c_str(), (rootstrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(singletonfn.c_str(), (singletonfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(dirbinfn.c_str(), (dirbinfn+".bsc").c_str(), 64);
	// mstcom::bsc::BSC_compress(readsnfn.c_str(), (readsnfn+".bsc").c_str(), 64);

	// mstcom::bsc::BSC_compress(sgfn.c_str(), (sgfn+".bsc").c_str(), 64);

	// string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc rootstr.txt.bsc singleton.bin.bsc dir.bin.bsc readsN.txt.bsc";
	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc rootstr.txt.bsc singleton.bin.bsc dir.bin.bsc";

	if (isorder) {
		// mstcom::bsc::BSC_compress(orderfn.c_str(), (orderfn+".bsc").c_str(), 64);
		mstcom::lzma::lzma_compress(orderfn.c_str(), (orderfn+".lzma").c_str());
		tarcmd += " order.bin.lzma";
	}

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

void outputSingle() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	
	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;

	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	size_t ss = sizeof(size_t);
	char *en_str = (char*)alloca((L + 1) * sizeof(char));
	int dir;
	uint32_t dn, difid;
	for (uint32_t rid = 0; rid < max_rid; ++rid) 
	if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
	// for (uint32_t i = 0; i < rootnodevec.size(); ++i) {
		queue<uint32_t> q;
		q.push(rid);
		uint32_t cnt = 0;
		while (!q.empty()) {
			uint32_t noderid = q.front();		
			q.pop();
			++cnt;
			for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
				q.push(reads[noderid].crid.a[i]);
			}
		}

		uint32_t rootid = rid;

		uint32_t num = 0;
		while (rootid < max_rid) {
			++ num;
			visited[rootid] = true;

			dn = reads[rootid].dn;
			if (dn > 0) {
				if (reads[rootid].isrc) { // sort 1 to 0
					sort(reads[rootid].dup, reads[rootid].dup + dn, cmpduprc1);
				} else { //small value 0 first
					sort(reads[rootid].dup, reads[rootid].dup + dn, cmpduprc0);
				}
			}

			if (num == 1) {
				if (dn > 0) {
					if (dn == 1)
						fprintf(fpenstr, "%s$\n", seq[rootid].seq);
					else 
						fprintf(fpenstr, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
					// fprintf(stdout, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
					for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
						dir = 0;
						if (reads[rootid].dup[w].isrc) dir = 1;
						bit_push(dirbin, fpdir, dir);
						visited[reads[rootid].dup[w].id] = true;
					}
				} else {
					fprintf(fpenstr, "%s\n", seq[rootid].seq);
				}
				rootid = reads[rootid].getChildren();
			} else {
				dir = 0;
				if (reads[rootid].isrc) dir = 1; // 
				bit_push(dirbin, fpdir, dir);
				
				if (dn > 0) {
					if (dn == 1)
						fprintf(fpenstr, "%s$\n", seq[rootid].seq);
					else 
						fprintf(fpenstr, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);

					for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
						dir = 0;
						if (reads[rootid].dup[w].isrc) dir = 1;
						bit_push(dirbin, fpdir, dir);
						
						visited[reads[rootid].dup[w].id] = true;
					}
				} else {
					fprintf(fpenstr, "%s\n", seq[rootid].seq);
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
				}
				rootid = child;
			}

		}
		// fprintf(fpenstr, "-\n");
	}
	
	// cout << "2222111 \n";

	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		kv_destroy(reads[rid].crid);
	}
	delete[] reads;
	fprintf(fpenstr, "-\n");

	uint8bit_v singlebin;
	singlebin.a = singlebin.n = 0;

	// sprintf(name, "singleton.bin");
	string singletonfn = folder + "singleton.bin"; 
	// FILE *fpsg = fopen(singletonfn.c_str(), "w");
	std::ofstream singleOfs(singletonfn.c_str(), std::ios::binary);

	int ncnt;

	uint32_t sgcnt = 0, prerid = 0, nrid;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (!visited[rid]) {
			ncnt = 0;
			for (int i = 0; i < L; ++i) {
				if (seq[rid].seq[i] == 'N') ++ncnt;
			}
			if (ncnt > 0) {
				fprintf(fpenstr, "%s\n", seq[rid].seq);
			} else {
				++ sgcnt;
				for (int i = 0; i < L; ++i) {
					DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
				}
			}
		}
	}

	if (singlebin.n > 0) {
		singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
	}
	singleOfs.close();

	fclose(fpenstr);
	delete[] visited;
	
	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	mstcom::bsc::BSC_compress(encodestrfn.c_str(), (encodestrfn+".bsc").c_str(), 64);
	// mstcom::bsc::BSC_compress(rootstrfn.c_str(), (rootstrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(singletonfn.c_str(), (singletonfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(dirbinfn.c_str(), (dirbinfn+".bsc").c_str(), 64);
	// mstcom::bsc::BSC_compress(readsnfn.c_str(), (readsnfn+".bsc").c_str(), 64);

	// mstcom::bsc::BSC_compress(sgfn.c_str(), (sgfn+".bsc").c_str(), 64);

	// string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc rootstr.txt.bsc singleton.bin.bsc dir.bin.bsc readsN.txt.bsc";
	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc dir.bin.bsc singleton.bin.bsc";

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

void outputSingleOrder_() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	
	std::ofstream fporder;
	string orderfn = folder + string("parent.bin");
	fporder.open(orderfn, std::ios::binary);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;

	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	int dir, isori;
	FILE *fpencoded = fopen("encoded.txt", "w");
	/// current is OK.........................................
	// cout << "1111 \n";
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (strlen(seq[rid].seq) == 0) { // duplicate reads caused by consuses reads
			fprintf(fpenstr, "-\n");
			fprintf(fpencoded, "-");
		}
		else 
		// if (strcmp(seq[rid].seq, "-") == 0) { // original duplicate reads
		if (seq[rid].seq[0] == '-') { // original duplicate reads
			fprintf(fpenstr, "\n");
			fprintf(fpencoded, "\n");
		} else { 
			fprintf(fpenstr, "%s\n", seq[rid].seq);
			fprintf(fpencoded, "%s", seq[rid].seq);
		}

		isori = checkIsOriReads(seq[rid].seq, L);
		// cout << isori << ",";
		if (isori == 0) {
			fporder.write((char*)&reads[rid].prid, sizeof(uint32_t));
			dir = 0;
			if (reads[rid].isrc) dir = 1;
			bit_push(dirbin, fpdir, dir);
			fprintf(fpencoded, " %u", reads[rid].prid);
		}
		fprintf(fpencoded, "\n");
	}
	fclose(fpencoded);
	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		kv_destroy(reads[rid].crid);
	}
	delete[] reads;
	fclose(fpenstr);

	fporder.close();

	mstcom::bsc::BSC_compress(encodestrfn.c_str(), (encodestrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(dirbinfn.c_str(), (dirbinfn+".bsc").c_str(), 64);
	// string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc rootstr.txt.bsc singleton.bin.bsc dir.bin.bsc readsN.txt.bsc";
	mstcom::lzma::lzma_compress(orderfn.c_str(), (orderfn+".lzma").c_str());
	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc dir.bin.bsc parent.bin.lzma";

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
	cout << "size of " << file + ".lzma" << ": " << N0 << endl;

	std::ifstream f1(file + ".bsc", std::ios::ate | std::ios::binary);
	if (f1.fail()) {
		std::cerr << "\nFile '" << file + ".bsc" << "' does not exist. Quit.";
		exit(1);
	}
	size_t N1 = f1.tellg();
	cout << "size of " << file + ".bsc" << ": " << N1 << endl;
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
		cout << files[i] << endl;
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

void outputSingleOrder() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	
	std::ofstream fporder;
	string orderfn = folder + string("parent.bin");
	fporder.open(orderfn, std::ios::binary);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;

	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	int dir, isori;
	FILE *fpencoded = fopen("encoded.txt", "w");
	/// current is OK.........................................
	// cout << "1111 \n";
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (strlen(seq[rid].seq) == 0) { // duplicate reads caused by consuses reads
			fprintf(fpenstr, "-\n");
			fprintf(fpencoded, "-");
		}
		else 
		// if (strcmp(seq[rid].seq, "-") == 0) { // original duplicate reads
		if (seq[rid].seq[0] == '-') { // original duplicate reads
			fprintf(fpenstr, "\n");
			fprintf(fpencoded, "\n");
		} else { 
			fprintf(fpenstr, "%s\n", seq[rid].seq);
			fprintf(fpencoded, "%s", seq[rid].seq);
		}

		isori = checkIsOriReads(seq[rid].seq, L);
		// cout << isori << ",";
		if (isori == 0) {
			fporder.write((char*)&reads[rid].prid, sizeof(uint32_t));
			dir = 0;
			if (reads[rid].isrc) dir = 1;
			bit_push(dirbin, fpdir, dir);
			fprintf(fpencoded, " %u", reads[rid].prid);
		}
		fprintf(fpencoded, "\n");
	}
	fclose(fpencoded);
	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		kv_destroy(reads[rid].crid);
	}
	delete[] reads;
	fclose(fpenstr);

	fporder.close();

	vector<string> f1s, f2s;
	f1s.push_back(encodestrfn);
	f1s.push_back(dirbinfn);

	f2s.push_back(orderfn);

	compressFiles(f1s, f2s);

	// mstcom::bsc::BSC_compress(encodestrfn.c_str(), (encodestrfn+".bsc").c_str(), 64);
	// mstcom::bsc::BSC_compress(dirbinfn.c_str(), (dirbinfn+".bsc").c_str(), 64);
	// // string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc rootstr.txt.bsc singleton.bin.bsc dir.bin.bsc readsN.txt.bsc";
	// mstcom::lzma::lzma_compress(orderfn.c_str(), (orderfn+".lzma").c_str());

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc dir.bin.bsc ";

	tarcmd += "parent.bin." + getExt(orderfn);

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

void outputPE() {
	stopwatch.resume();
	// char *str = (char*)alloca((L + 1) * sizeof(char));
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	// FILE *fpenstrcopy = fopen("encodestr.txt.copy", "w");
	// FILE *fpenstrcopycopy = fopen("encodestr.txt.copy.copy", "w");
	string rootstrfn = folder + "rootstr.txt";
	FILE *fproot = fopen(rootstrfn.c_str(), "w");

	std::ofstream fporder;
	string orderfn = folder + string("order.bin");
		// string orderfn = string("order.bin");  // for test
		// fporder = fopen(orderfn.c_str(), "wb");
	fporder.open(orderfn, std::ios::binary);

	string singletonfn = folder + "singleton.bin"; 
	std::ofstream singleOfs(singletonfn.c_str(), std::ios::binary);

	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	// size_t outputnum = 0;

	uint8bit_v dirbin;
	// uint8bit_v rootbin, dirbin;
	// rootbin.n = rootbin.a = 0;
	dirbin.n = dirbin.a = 0;
	// sprintf(name, "dir.bin");
	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	size_t ss = sizeof(size_t);
	char *en_str = (char*)alloca((L + 1) * sizeof(char));

	// int max_shift = 0;
	size_t frid;

	if (isorder) {
		for (uint32_t rid = 0; rid < max_rid; ++rid) {
			if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
				queue<uint32_t> q;
				q.push(rid);
				uint32_t cnt = 0;
				while (!q.empty()) {
					uint32_t noderid = q.front();		
					q.pop();
					++cnt;
					for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
						q.push(reads[noderid].crid.a[i]);
					}
				}

				uint32_t rootid = rid;
					uint32_t num = 0;
				while (rootid < max_rid) {
					// fprintf(stderr, "%lu\n", rootid);
					++ num;
					visited[rootid] = true;

					// if (isorder) fwrite(&rootid, ss, 1, fporder);
					int file = 0;
					if (rootid >= (max_rid>>1)) file = 1; // right ends
					bit_push(dirbin, fpdir, file);

					if (isorder) {
						if (file) {
							frid = rootid - (max_rid >> 1);
							fporder.write((char*)&frid, sizeof(uint32_t));
						} else {
							fporder.write((char*)&rootid, sizeof(uint32_t));
						}
					} else { //
						if (file) { // right ends
							// rightmap[rightmapid++] = rootid - (max_rid >> 1);
						} else {
							// leftmap[rootid] = leftmapid++;
						}
					}


					if (num == 1) {
						// fprintf(fproot, "%s\n", reads[rootid].str.c_str());
						fprintf(fproot, "%s\n", seq[rootid].seq);

						// for (int j = 0; j < L; ++j) {
						// 	DNA_push(rootbin, rootOfs, seq_nt4_table[(uint8_t)seq[rootid].seq[j]]);
						// }
						rootid = reads[rootid].getChildren();
					} else {
						int dir = 0;
						// fprintf(fpenstrcopy, "%s\n", seq[rootid].seq);
						// if (abs(reads[rootid].shift) > max_shift) {
						// 	max_shift = abs(reads[rootid].shift);
						// }

						fprintf(fpenstr, "%s\n", seq[rootid].seq);
						// fprintf(fpencoded, "%s\n", seq[rootid].seq);
						if (reads[rootid].isrc) {
							dir = 1;
						} 

						bit_push(dirbin, fpdir, dir);
						// dir = 0;
						// if (rootid >= (max_rid>>1)) dir = 1;
						// bit_push(dirbin, fpdir, dir);

						if (num >= cnt) {
							rootid = max_rid + 1;
							break;
						}
						size_t child = reads[rootid].getChildren();
						// fprintf(stderr, "child: %lu\n", child);
						int back_step = 0;
						while (child > max_rid) { // no leaf node or all children are visited
							// back
							rootid = reads[rootid].prid;
							// fprintf(stderr, "rootid: %lu\n", rootid);
							child = reads[rootid].getChildren();
							++back_step;
							// ++num;
							// if (num > 20) exit(0);
						}
						// fprintf(stderr, "child: %lu\n", child);
						if (back_step > 0) {
							// fprintf(stderr, "--%d\n", back_step);
							fprintf(fpenstr, "%d\n", back_step);
						}
						rootid = child;
						// if (num > 20) exit(0);
					}
				}
				// exit(0);
				fprintf(fpenstr, "-\n");
			} 
		}
		if (dirbin.n > 0) { 
			fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
		}
		fpdir.close();

		fclose(fpenstr);

		// fprintf(stderr, "max_shift: %d\n", max_shift);
		// fclose(fpenstrcopy);
		// fclose(fpenstrcopycopy);

		uint8bit_v singlebin;
		singlebin.a = singlebin.n = 0;

		// sprintf(name, "singleton.bin");
		string singletonfn = folder + "singleton.bin"; 
		std::ofstream singleOfs(singletonfn.c_str(), std::ios::binary);

		int ncnt;

		fprintf(fproot, "-\n");
		for (size_t rid = 0; rid < (max_rid>>1); ++rid) {
			if (!visited[rid]) {
				ncnt = 0;
				for (int i = 0; i < L; ++i) {
					if (seq[rid].seq[i] == 'N') ++ncnt;
				}
				if (ncnt > 0) {
					fprintf(fproot, "%s\n", seq[rid].seq);
					// if (isquality) kv_push(size_t, nrid, rid);
					if (isorder) {
						fporder.write((char*)&rid, sizeof(uint32_t));
					} else {
						// leftmap[rid] = leftmapid++;
						// fwrite(&rid, ss, 1, fporder);
					}
					visited[rid] = true;
				} 
			}
		}
		fclose(fproot);

		uint32_t sgcnt = 0;
		for (size_t rid = 0; rid < (max_rid>>1); ++rid) {
			if (!visited[rid]) {
				++ sgcnt;
				// leftmap[rid] = leftmapid++;
				// fprintf(fpsingle, "%s\n", seq[rid].seq);
				for (int i = 0; i < L; ++i) {
					DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
				}
				if (!isorder) {
					// leftmap[rid] = leftmapid++;
				}
			}
		}
		cout << "sgcnt: " << sgcnt << endl;

		/////////////
		string rightreadsnfn = folder + "readsN.txt.right";
		FILE *rightfpN = fopen(rightreadsnfn.c_str(), "w");

		for (size_t rid = (max_rid>>1); rid < max_rid; ++rid) {
			if (!visited[rid]) {
				ncnt = 0;
				for (int i = 0; i < L; ++i) {
					if (seq[rid].seq[i] == 'N') ++ncnt;
				}
				if (ncnt > 0) {
					fprintf(rightfpN, "%s\n", seq[rid].seq);
					// if (isquality) kv_push(size_t, nrid, rid);
					if (isorder) {	
						frid = rid - (max_rid >> 1);
						fporder.write((char*)&frid, sizeof(uint32_t));
					} else {
						// rightmap[rightmapid++] = rid - (max_rid >> 1);
					}
					visited[rid] = true;
					// fwrite(&rid, ss, 1, fporder);
				} 
			}
		}
		fclose(rightfpN);
		fclose(fproot);

		for (size_t rid = (max_rid>>1); rid < max_rid; ++rid) {
			if (!visited[rid]) {
				if (!isorder) {
					// rightmap[rightmapid++] = rid - (max_rid >> 1);
				}
				// fprintf(fpsingle, "%s\n", seq[rid].seq);
				for (int i = 0; i < L; ++i) {
					DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
				}
			}
		}
		if (singlebin.n > 0) {
			singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
		}
		singleOfs.close();
		fporder.close();
	} else {
		size_t *leftmap, *rightmap, leftmapid, rightmapid;
		leftmap = new size_t[max_rid >> 1];
		rightmap = new size_t[max_rid >> 1];
		leftmapid = rightmapid = 0;
		int file, dir;
		uint32_t curid, dn;

		for (uint32_t rid = 0; rid < max_rid; ++rid) {
			if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
				queue<uint32_t> q;
				q.push(rid);
				uint32_t cnt = 0;
				while (!q.empty()) {
					uint32_t noderid = q.front();		
					q.pop();
					++cnt;
					for (size_t i = 0; i < reads[noderid].crid.n; ++i) {
						q.push(reads[noderid].crid.a[i]);
					}
				}

				uint32_t rootid = rid;
				uint32_t num = 0;
				while (rootid < max_rid) {
					// fprintf(stderr, "%lu\n", rootid);
					++ num;
					visited[rootid] = true;

					// if (isorder) fwrite(&rootid, ss, 1, fporder);
					file = 0;
					if (rootid >= (max_rid>>1)) file = 1; // right ends
					bit_push(dirbin, fpdir, file);
					if (file) { // right ends
						rightmap[rightmapid++] = rootid - (max_rid >> 1);
					} else {
						leftmap[rootid] = leftmapid++;
					}

					dn = reads[rootid].dn;
					if (dn > 0) {
						if (reads[rootid].isrc) { // sort 1 to 0
							sort(reads[rootid].dup, reads[rootid].dup + dn, cmpduprc1);
						} else { //small value 0 first
							sort(reads[rootid].dup, reads[rootid].dup + dn, cmpduprc0);
						}
					}

					if (num == 1) {
						// fprintf(fproot, "%s\n", reads[rootid].str.c_str());
						if (dn > 0) {
							fprintf(fproot, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
							for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
								curid = reads[rootid].dup[w].id;
								visited[curid] = true;

								file = 0;
								if (curid >= (max_rid>>1)) file = 1; // right ends
								bit_push(dirbin, fpdir, file);
								if (file) { // right ends
									rightmap[rightmapid++] = curid - (max_rid >> 1);
								} else {
									leftmap[curid] = leftmapid++;
								}

								dir = 0;
								if (reads[rootid].dup[w].isrc) dir = 1;
								bit_push(dirbin, fpdir, dir);
							}
						} else {
							fprintf(fproot, "%s\n", seq[rootid].seq);
						}
						rootid = reads[rootid].getChildren();
					} else {
						dir = 0;
						if (reads[rootid].isrc) dir = 1; // 
						bit_push(dirbin, fpdir, dir);

						if (dn > 0) {
							fprintf(fpenstr, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
							for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
								curid = reads[rootid].dup[w].id;
								visited[curid] = true;

								file = 0;
								if (curid >= (max_rid>>1)) file = 1; // right ends
								bit_push(dirbin, fpdir, file);
								if (file) { // right ends
									rightmap[rightmapid++] = curid - (max_rid >> 1);
								} else {
									leftmap[curid] = leftmapid++;
								}

								dir = 0;
								if (reads[rootid].dup[w].isrc) dir = 1;
								bit_push(dirbin, fpdir, dir);								
							}
						} else {
							fprintf(fpenstr, "%s\n", seq[rootid].seq);
						}

						if (num >= cnt) {
							rootid = max_rid + 1;
							break;
						}
						size_t child = reads[rootid].getChildren();
						// fprintf(stderr, "child: %lu\n", child);
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
						// if (num > 20) exit(0);
					}
				}
				// exit(0);
				fprintf(fpenstr, "-\n");
			} 
		}

		if (dirbin.n > 0) { 
			fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
		}
		fpdir.close();

		fclose(fpenstr);

		uint8bit_v singlebin;
		singlebin.a = singlebin.n = 0;

		int ncnt;

		fprintf(fproot, "-\n");
		for (size_t rid = 0; rid < (max_rid>>1); ++rid) {
			if (!visited[rid]) {
				ncnt = 0;
				for (int i = 0; i < L; ++i) {
					if (seq[rid].seq[i] == 'N') ++ncnt;
				}
				if (ncnt > 0) {
					fprintf(fproot, "%s\n", seq[rid].seq);
					leftmap[rid] = leftmapid++;
					visited[rid] = true;
				} 
			}
		}

		uint32_t sgcnt = 0;
		for (size_t rid = 0; rid < (max_rid>>1); ++rid) {
			if (!visited[rid]) {
				++ sgcnt;
				// leftmap[rid] = leftmapid++;
				// fprintf(fpsingle, "%s\n", seq[rid].seq);
				for (int i = 0; i < L; ++i) {
					DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
				}
				leftmap[rid] = leftmapid++;
			}
		}
		cout << "sgcnt: " << sgcnt << endl;

		/////////////
		fprintf(fproot, "-\n");
		for (size_t rid = (max_rid>>1); rid < max_rid; ++rid) {
			if (!visited[rid]) {
				ncnt = 0;
				for (int i = 0; i < L; ++i) {
					if (seq[rid].seq[i] == 'N') ++ncnt;
				}
				if (ncnt > 0) {
					fprintf(fproot, "%s\n", seq[rid].seq);
					rightmap[rightmapid++] = rid - (max_rid >> 1);
					visited[rid] = true;
				} 
			}
		}
		fclose(fproot);

		for (size_t rid = (max_rid>>1); rid < max_rid; ++rid) {
			if (!visited[rid]) {
				rightmap[rightmapid++] = rid - (max_rid >> 1);
				for (int i = 0; i < L; ++i) {
					DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
				}
			}
		}
		if (singlebin.n > 0) {
			singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
		}
		singleOfs.close();

		for (size_t rid = 0; rid < (max_rid>>1); ++rid) {
			fporder.write((char*)&leftmap[rightmap[rid]], sizeof(uint32_t));
		}
		delete[] leftmap;
		delete[] rightmap;
		fporder.close();
	}


	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	mstcom::bsc::BSC_compress(encodestrfn.c_str(), (encodestrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(rootstrfn.c_str(), (rootstrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(dirbinfn.c_str(), (dirbinfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(singletonfn.c_str(), (singletonfn+".bsc").c_str(), 64);
	mstcom::lzma::lzma_compress(orderfn.c_str(), (orderfn+".lzma").c_str());

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc rootstr.txt.bsc dir.bin.bsc singleton.bin.bsc order.bin.lzma";

	// cout << tarcmd << endl;

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

#if false
void outputPEX() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	string rootstrfn = folder + "rootstr.txt";
	FILE *fproot = fopen(rootstrfn.c_str(), "w");

	std::ofstream fporder;
	string orderfn = folder + string("order.bin");
	if (isorder) {
		fporder.open(orderfn, std::ios::binary);
	}

	std::ofstream fpdist;
	string distfn = folder + string("dist.bin");
	fpdist.open(distfn, std::ios::binary);

	vector<ROOTNODE_t> rootnodevec;
	for (uint32_t rid = 0; rid < max_rid; ++rid) {
		if (reads[rid].prid == rid && reads[rid].crid.n > 0) { //is the root node && not a leaf node == not a singleton reads
			queue<uint32_t> q;
			q.push(rid);
			uint32_t cnt = 0;
			while (!q.empty()) {
				uint32_t noderid = q.front();		
				q.pop();
				++cnt;
				for (uint32_t i = 0; i < reads[noderid].crid.n; ++i) {
					q.push(reads[noderid].crid.a[i]);
				}
			}
			rootnodevec.push_back(ROOTNODE_t(rid, cnt));
		}
	}
	sort(rootnodevec.begin(), rootnodevec.end(), cmp);

	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;
	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	char *en_str = (char*)alloca((L + 1) * sizeof(char));

	uint32_t *ids = new uint32_t[max_rid + 1];
	uint32_t idsid = 0;
	uint32_t *od = new uint32_t[max_rid + 3];

	// FILE *fpid = fopen("idx.txt", "w");
	for (uint32_t i = 0; i < rootnodevec.size(); ++i) {
			uint32_t cnt = rootnodevec[i].nodecnt;

			uint32_t rootid = rootnodevec[i].rid;
			uint32_t num = 0;
			while (rootid < max_rid) {
				ids[idsid] = rootid;
				od[rootid] = idsid ++;

				// fprintf(fpid, "%lu\n", rootid);
				++ num;
				visited[rootid] = true;

				if (num == 1) {
					fprintf(fproot, "%s\n", seq[rootid].seq);
					rootid = reads[rootid].getChildren();
				} else {
					int dir = 0;
					if (reads[rootid].isrc) {
						dir = 1;
						strcpy(rcstr, seq[rootid].seq);
						reverseComplement(rcstr);
						encode(seq[reads[rootid].prid].seq, rcstr, reads[rootid].shift, en_str);
						fprintf(fpenstr, "%s\n", en_str);
					} else {
						encode(seq[reads[rootid].prid].seq, seq[rootid].seq, reads[rootid].shift, en_str);
						fprintf(fpenstr, "%s\n", en_str);
					}

					bit_push(dirbin, fpdir, dir);

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
			fprintf(fpenstr, "-\n");
	}

	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	// rootnodevec.clear();

	fclose(fproot);
	
	fclose(fpenstr);

	uint8bit_v singlebin;
	singlebin.a = singlebin.n = 0;

	string singletonfn = folder + "singleton.bin"; 
	std::ofstream singleOfs(singletonfn.c_str(), std::ios::binary);

	string readsnfn = folder + "readsN.txt";
	FILE *fpN = fopen(readsnfn.c_str(), "w");
	int ncnt;

	for (uint32_t rid = 0; rid < (max_rid); ++rid) {
		if (!visited[rid]) {
			ncnt = 0;
			for (int i = 0; i < L; ++i) {
				if (seq[rid].seq[i] == 'N') ++ncnt;
			}
			if (ncnt > 0) {
				ids[idsid] = rid;
				od[rid] = idsid ++;
				// fprintf(fpid, "%lu\n", rid);
				fprintf(fpN, "%s\n", seq[rid].seq);
				visited[rid] = true;
			} 
		}
	}
	fclose(fpN);

	for (uint32_t rid = 0; rid < (max_rid); ++rid) {
		if (!visited[rid]) {
			ids[idsid] = rid;
			od[rid] = idsid ++;
			// fprintf(fpid, "%lu\n", rid);
			for (int i = 0; i < L; ++i) {
				DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
			}
		}
	}
	// fclose(fpid);

	if (singlebin.n > 0) {
		singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
	}
	singleOfs.close();

	// calc dist
	FILE *fp = fopen("dist.txt", "w");
	uint32_t *tr, M;
	build(tr, M, max_rid);
	memset(visited, false, sizeof(bool)*max_rid);

	uint32_t half = max_rid >> 1;
	uint32_t id = 0, num = 0, p, r, dist, orione;
	int32_t res;
	// 在线段数上 使用的是id+1, p+1
	while (id < max_rid && num < half) {
		if (!visited[id]) {
			uint32_t dis = 0;
			if (ids[id] < half) {
				p = od[half + ids[id]];

				r = query(tr, M, id + 1 + 1, p + 1);
				dist = p + 1 - (id + 1 + 1) + 1 - r;
				res = (int32_t)dist;
				if (isorder) {
					fporder.write((char *)&ids[id], sizeof(uint32_t));
				}
			} else {
				p = od[ids[id] - half];
				r = query(tr, M, id + 1 + 1, p + 1);
				dist = p + 1 - (id + 1 + 1) + 1 - r;
				res = 0 - (int32_t)dist;
				if (isorder) {
					orione = ids[id] - half;
					fporder.write((char *)&orione, sizeof(uint32_t));
				}
			}
			fprintf(fp, "%d\n", res);

			fpdist.write((char *)&res, sizeof(int32_t));

			update(tr, M, p + 1, 1);
			visited[p] = true;

			++ num;
		}
		++ id;
	}

	fpdist.close();
	fporder.close();

	delete[] visited;
	delete[] od;
	delete[] ids;

	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	mstcom::bsc::BSC_compress(encodestrfn.c_str(), (encodestrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(rootstrfn.c_str(), (rootstrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(dirbinfn.c_str(), (dirbinfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(readsnfn.c_str(), (readsnfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(distfn.c_str(), (distfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(singletonfn.c_str(), (singletonfn+".bsc").c_str(), 64);
	if (isorder) {
		mstcom::bsc::BSC_compress(orderfn.c_str(), (orderfn+".bsc").c_str(), 64);
	}

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc rootstr.txt.bsc dir.bin.bsc singleton.bin.bsc readsN.txt.bsc dist.bin.bsc";

	if (isorder) {
		tarcmd += " order.bin.bsc";
	}

	cout << tarcmd << endl;

	system(tarcmd.c_str());
	cout << "Time of bsc files = " << stopwatch.stop() << std::endl;

}

void outputPEX_backup() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	string rootstrfn = folder + "rootstr.txt";
	FILE *fproot = fopen(rootstrfn.c_str(), "w");

	std::ofstream fporder;
	string orderfn = folder + string("order.bin");
	if (isorder) {
		fporder.open(orderfn, std::ios::binary);
	}

	std::ofstream fpdist;
	string distfn = folder + string("dist.bin");
	fpdist.open(distfn, std::ios::binary);

	bool *visited = new bool[max_rid];
	memset(visited, false, sizeof(bool)*max_rid);

	int trees_cnt = 0;
	uint8bit_v dirbin;
	dirbin.n = dirbin.a = 0;
	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	char *en_str = (char*)alloca((L + 1) * sizeof(char));

	uint32_t *ids = new uint32_t[max_rid + 1];
	uint32_t idsid = 0, dn, curid;
	uint32_t *od = new uint32_t[max_rid + 3];
	int dir;

	// FILE *fpid = fopen("idx.txt", "w");
	if (isorder) {
		for (uint32_t rid = 0; rid < max_rid; ++rid) {
			if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
			// for (uint32_t i = 0; i < rootnodevec.size(); ++i) {
				queue<uint32_t> q;
				q.push(rid);
				uint32_t cnt = 0;
				while (!q.empty()) {
					uint32_t noderid = q.front();		
					q.pop();
					++cnt;
					sort(reads[noderid].crid.a, reads[noderid].crid.a + reads[noderid].crid.n);
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

					if (dn > 0) {
						reads[rootid].dup[0].isrc = reads[rootid].isrc;
					}

					if (num == 1) {
						if (dn > 0) {
							fprintf(fproot, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
							// fprintf(stdout, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
							for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
								dir = 0;
								if (reads[rootid].dup[w].isrc) dir = 1;
								bit_push(dirbin, fpdir, dir);

								curid = reads[rootid].dup[w].id;
								visited[curid] = true;

								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						} else {
							fprintf(fproot, "%s\n", seq[rootid].seq);
						}

						rootid = reads[rootid].getChildren();
					} else {
						dir = 0;
						if (reads[rootid].isrc) dir = 1;
						bit_push(dirbin, fpdir, dir);

						if (dn > 0) {
							fprintf(fpenstr, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);

							for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
								dir = 0;
								if (reads[rootid].dup[w].isrc) dir = 1;
								bit_push(dirbin, fpdir, dir);

								curid = reads[rootid].dup[w].id;
								visited[curid] = true;

								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						} else {
							fprintf(fpenstr, "%s\n", seq[rootid].seq);
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
				fprintf(fpenstr, "-\n");
			}
		}
	} else {
		for (uint32_t rid = 0; rid < max_rid; ++rid) {
			if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
			// for (uint32_t i = 0; i < rootnodevec.size(); ++i) {
				queue<uint32_t> q;
				q.push(rid);
				uint32_t cnt = 0;
				while (!q.empty()) {
					uint32_t noderid = q.front();		
					q.pop();
					++cnt;
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

					if (dn > 0) {
						// reads[rootid].dup[dn] = dup_t(rootid, reads[rootid].isrc);
						if (reads[rootid].isrc) { // sort 1 to 0
							sort(reads[rootid].dup, reads[rootid].dup + dn, cmpduprc1);
						} else { //small value 0 first
							sort(reads[rootid].dup, reads[rootid].dup + dn, cmpduprc0);
						}
					}

					if (num == 1) {
						if (dn > 0) {
							fprintf(fproot, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
							// fprintf(stdout, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
							for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
								dir = 0;
								if (reads[rootid].dup[w].isrc) dir = 1;
								bit_push(dirbin, fpdir, dir);

								curid = reads[rootid].dup[w].id;
								visited[curid] = true;

								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						} else {
							fprintf(fproot, "%s\n", seq[rootid].seq);
						}

						rootid = reads[rootid].getChildren();
					} else {
						dir = 0;
						if (reads[rootid].isrc) dir = 1;
						bit_push(dirbin, fpdir, dir);

						if (dn > 0) {
							fprintf(fpenstr, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);

							for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
								dir = 0;
								if (reads[rootid].dup[w].isrc) dir = 1;
								bit_push(dirbin, fpdir, dir);

								curid = reads[rootid].dup[w].id;
								visited[curid] = true;

								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						} else {
							fprintf(fpenstr, "%s\n", seq[rootid].seq);
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
				fprintf(fpenstr, "-\n");
			}
		}
	}

	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();
	// rootnodevec.clear();
	fclose(fpenstr);
	// string readsnfn = folder + "readsN.txt";
	// FILE *fpN = fopen(readsnfn.c_str(), "w");
	int ncnt;
	fprintf(fproot, "-\n");
	for (uint32_t rid = 0; rid < (max_rid); ++rid) {
		if (!visited[rid]) {
			ncnt = 0;
			for (int i = 0; i < L; ++i) {
				if (seq[rid].seq[i] == 'N') ++ncnt;
			}
			if (ncnt > 0) {
				ids[idsid] = rid;
				od[rid] = idsid ++;
				// fprintf(fpid, "%lu\n", rid);
				fprintf(fproot, "%s\n", seq[rid].seq);
				visited[rid] = true;
			} 
		}
	}
	fclose(fproot);

	uint8bit_v singlebin;
	singlebin.a = singlebin.n = 0;

	string singletonfn = folder + "singleton.bin"; 
	std::ofstream singleOfs(singletonfn.c_str(), std::ios::binary);
	for (uint32_t rid = 0; rid < (max_rid); ++rid) {
		if (!visited[rid]) {
			ids[idsid] = rid;
			od[rid] = idsid ++;
			// fprintf(fpid, "%lu\n", rid);
			for (int i = 0; i < L; ++i) {
				DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
			}
		}
	}
	// fclose(fpid);

	if (singlebin.n > 0) {
		singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
	}
	singleOfs.close();

	// calc dist
	FILE *fp = fopen("dist.txt", "w");
	uint32_t *tr, M;
	build(tr, M, max_rid);
	memset(visited, false, sizeof(bool)*max_rid);

	uint32_t half = max_rid >> 1;
	uint32_t id = 0, num = 0, p, r, dist, orione;
	int32_t res;

	// uint8bit_v filebin, rangebin;
	// filebin.a = filebin.n = 0;
	// rangebin.a = rangebin.n = 0;

	// string filefn = folder + "file.bin"; 
	// std::ofstream fileOfs(filefn.c_str(), std::ios::binary);
	// string rangefn = folder + "range.bin"; 
	// std::ofstream rangeOfs(rangefn.c_str(), std::ios::binary);

	// cout << "before segment array\n";
	// in the segment array, (id+1, p+1)
	bool debug = false;
	int file;
	uint16_t disttemp;
	uint32_t small = 0;
	if (isorder) {
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

				// bit_push(filebin, fileOfs, file);

				fprintf(fp, "%d\n", res);

				/*if (dist < 65536) {
					++ small;
					file = 0;
					disttemp = dist;
					fpdist.write((char *)&disttemp, sizeof(uint16_t));
				} else {
					file = 1;
					fpdist.write((char *)&dist, sizeof(uint32_t));
				}*/
				fpdist.write((char *)&res, sizeof(int32_t));
				// bit_push(rangebin, rangeOfs, file);

				update(tr, M, p + 1, 1);
				visited[p] = true;

				++ num;
			}
			++ id;
		}
	} else {
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

				// bit_push(filebin, fileOfs, file);

				fprintf(fp, "%d\n", res);

				/*if (dist < 65536) {
					++ small;
					file = 0;
					disttemp = dist;
					fpdist.write((char *)&disttemp, sizeof(uint16_t));
				} else {
					file = 1;
					fpdist.write((char *)&dist, sizeof(uint32_t));
				}*/
				fpdist.write((char *)&res, sizeof(int32_t));
				// bit_push(rangebin, rangeOfs, file);

				update(tr, M, p + 1, 1);
				visited[p] = true;

				++ num;
			}
			++ id;
		}
	}
	cout << "end segment array\n";
	// fprintf(stderr, "%u\n", small);
	// if (filebin.n > 0) { 
	// 	fileOfs.write((char*)&filebin.a, sizeof(uint8_t));
	// }
	// fileOfs.close();

	// if (rangebin.n > 0) { 
	// 	rangeOfs.write((char*)&rangebin.a, sizeof(uint8_t));
	// }
	// rangeOfs.close();

	fpdist.close();
	if (isorder) fporder.close();

	delete[] visited;
	delete[] od;
	delete[] ids;

	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	mstcom::bsc::BSC_compress(encodestrfn.c_str(), (encodestrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(rootstrfn.c_str(), (rootstrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(dirbinfn.c_str(), (dirbinfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(singletonfn.c_str(), (singletonfn+".bsc").c_str(), 64);
	mstcom::lzma::lzma_compress(distfn.c_str(), (distfn+".lzma").c_str());
	if (isorder) {
		mstcom::lzma::lzma_compress(orderfn.c_str(), (orderfn+".lzma").c_str());
	}

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc rootstr.txt.bsc dir.bin.bsc singleton.bin.bsc dist.bin.lzma";

	if (isorder) {
		tarcmd += " order.bin.lzma";
	}

	cout << tarcmd << endl;

	system(tarcmd.c_str());
	cout << "Time of bsc files = " << stopwatch.stop() << std::endl;

}
#endif

void outputPEX() {
	stopwatch.resume();
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));

	string encodestrfn = folder + "encodestr.txt";
	FILE *fpenstr = fopen(encodestrfn.c_str(), "w");
	string rootstrfn = folder + "rootstr.txt";
	FILE *fproot = fopen(rootstrfn.c_str(), "w");

	std::ofstream fporder;
	string orderfn = folder + string("order.bin");
	if (isorder) {
		fporder.open(orderfn, std::ios::binary);
	}

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
	string dirbinfn = folder + "dir.bin";
	std::ofstream fpdir(dirbinfn.c_str(), std::ios::binary);

	char *en_str = (char*)alloca((L + 1) * sizeof(char));

	uint32_t *ids = new uint32_t[max_rid + 1];
	uint32_t idsid = 0, dn, curid;
	uint32_t *od = new uint32_t[max_rid + 3];
	int dir;

	// FILE *fpid = fopen("idx.txt", "w");
	if (isorder) {
		for (uint32_t rid = 0; rid < max_rid; ++rid) {
			if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
			// for (uint32_t i = 0; i < rootnodevec.size(); ++i) {
				queue<uint32_t> q;
				q.push(rid);
				uint32_t cnt = 0;
				while (!q.empty()) {
					uint32_t noderid = q.front();		
					q.pop();
					++cnt;
					sort(reads[noderid].crid.a, reads[noderid].crid.a + reads[noderid].crid.n);
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

					if (dn > 0) {
						reads[rootid].dup[0].isrc = reads[rootid].isrc;
					}

					if (num == 1) {
						if (dn > 0) {
							if (dn == 1)
								fprintf(fproot, "%s$\n", seq[rootid].seq);
							else 
								fprintf(fproot, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
							// fprintf(stdout, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
							for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
								dir = 0;
								if (reads[rootid].dup[w].isrc) dir = 1;
								bit_push(dirbin, fpdir, dir);

								curid = reads[rootid].dup[w].id;
								visited[curid] = true;

								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						} else {
							fprintf(fproot, "%s\n", seq[rootid].seq);
						}

						rootid = reads[rootid].getChildren();
					} else {
						dir = 0;
						if (reads[rootid].isrc) dir = 1;
						bit_push(dirbin, fpdir, dir);

						if (dn > 0) {
							if (dn == 1)
								fprintf(fpenstr, "%s$\n", seq[rootid].seq);
							else 
								fprintf(fpenstr, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);

							for (uint32_t w = 1; w < reads[rootid].dn + 1; ++w) {
								dir = 0;
								if (reads[rootid].dup[w].isrc) dir = 1;
								bit_push(dirbin, fpdir, dir);

								curid = reads[rootid].dup[w].id;
								visited[curid] = true;

								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						} else {
							fprintf(fpenstr, "%s\n", seq[rootid].seq);
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
				fprintf(fpenstr, "-\n");
			}
		}
	} else {
		for (uint32_t rid = 0; rid < max_rid; ++rid) {
			if (reads[rid].prid == rid && (reads[rid].crid.n > 0 || reads[rid].dn > 0)) { //is the root node && not a leaf node == not a singleton reads
			// for (uint32_t i = 0; i < rootnodevec.size(); ++i) {
				queue<uint32_t> q;
				q.push(rid);
				uint32_t cnt = 0;
				while (!q.empty()) {
					uint32_t noderid = q.front();		
					q.pop();
					++cnt;
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

					if (dn > 0) {
						// reads[rootid].dup[dn] = dup_t(rootid, reads[rootid].isrc);
						if (reads[rootid].isrc) { // sort 1 to 0
							sort(reads[rootid].dup, reads[rootid].dup + dn, cmpduprc1);
						} else { //small value 0 first
							sort(reads[rootid].dup, reads[rootid].dup + dn, cmpduprc0);
						}
					}

					if (num == 1) {
						if (dn > 0) {
							if (dn == 1)
								fprintf(fproot, "%s$\n", seq[rootid].seq);
							else 
								fprintf(fproot, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
							// fprintf(stdout, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);
							for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
								dir = 0;
								if (reads[rootid].dup[w].isrc) dir = 1;
								bit_push(dirbin, fpdir, dir);

								curid = reads[rootid].dup[w].id;
								visited[curid] = true;

								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						} else {
							fprintf(fproot, "%s\n", seq[rootid].seq);
						}

						rootid = reads[rootid].getChildren();
					} else {
						dir = 0;
						if (reads[rootid].isrc) dir = 1;
						bit_push(dirbin, fpdir, dir);

						if (dn > 0) {
							if (dn == 1)
								fprintf(fpenstr, "%s$\n", seq[rootid].seq);
							else 
								fprintf(fpenstr, "%s$%u\n", seq[rootid].seq, reads[rootid].dn);

							for (uint32_t w = 0; w < reads[rootid].dn; ++w) {
								dir = 0;
								if (reads[rootid].dup[w].isrc) dir = 1;
								bit_push(dirbin, fpdir, dir);

								curid = reads[rootid].dup[w].id;
								visited[curid] = true;

								ids[idsid] = curid;
								od[curid] = idsid ++;
							}
						} else {
							fprintf(fpenstr, "%s\n", seq[rootid].seq);
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
				fprintf(fpenstr, "-\n");
			}
		}
	}

	// rootnodevec.clear();
	fclose(fpenstr);
	// string readsnfn = folder + "readsN.txt";
	// FILE *fpN = fopen(readsnfn.c_str(), "w");
	int ncnt;
	fprintf(fproot, "-\n");
	for (uint32_t rid = 0; rid < (max_rid); ++rid) {
		if (!visited[rid]) {
			ncnt = 0;
			for (int i = 0; i < L; ++i) {
				if (seq[rid].seq[i] == 'N') ++ncnt;
			}
			if (ncnt > 0) {
				ids[idsid] = rid;
				od[rid] = idsid ++;
				// fprintf(fpid, "%lu\n", rid);
				fprintf(fproot, "%s\n", seq[rid].seq);
				visited[rid] = true;
			} 
		}
	}
	fclose(fproot);

	uint8bit_v singlebin;
	singlebin.a = singlebin.n = 0;

	string singletonfn = folder + "singleton.bin"; 
	std::ofstream singleOfs(singletonfn.c_str(), std::ios::binary);
	for (uint32_t rid = 0; rid < (max_rid); ++rid) {
		if (!visited[rid]) {
			ids[idsid] = rid;
			od[rid] = idsid ++;
			// fprintf(fpid, "%lu\n", rid);
			for (int i = 0; i < L; ++i) {
				DNA_push(singlebin, singleOfs, seq_nt4_table[(uint8_t)seq[rid].seq[i]]);
			}
		}
	}
	// fclose(fpid);

	if (singlebin.n > 0) {
		singleOfs.write((char*)&singlebin.a, sizeof(uint8_t));
	}
	singleOfs.close();

	// calc dist
	// FILE *fp = fopen("dist.txt", "w");
	uint32_t *tr, M;
	build(tr, M, max_rid);
	memset(visited, false, sizeof(bool)*max_rid);

	uint32_t half = max_rid >> 1;
	uint32_t id = 0, num = 0, p, r, dist, orione;
	int32_t res;

	// uint8bit_v filebin, rangebin;
	// filebin.a = filebin.n = 0;
	// rangebin.a = rangebin.n = 0;

	// string filefn = folder + "file.bin"; 
	// std::ofstream fileOfs(filefn.c_str(), std::ios::binary);
	// string rangefn = folder + "range.bin"; 
	// std::ofstream rangeOfs(rangefn.c_str(), std::ios::binary);

	// cout << "before segment array\n";
	// in the segment array, (id+1, p+1)
	bool debug = false;
	int file;
	uint16_t disttemp;
	uint32_t small = 0;
	if (isorder) {
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
				bit_push(dirbin, fpdir, file);
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
				bit_push(dirbin, fpdir, file);
				// fpdist.write((char *)&res, sizeof(int32_t));
				// bit_push(rangebin, rangeOfs, file);

				update(tr, M, p + 1, 1);
				visited[p] = true;

				++ num;
			}
			++ id;
		}
	} else {
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
				bit_push(dirbin, fpdir, file);
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
				bit_push(dirbin, fpdir, file);
				// fpdist.write((char *)&res, sizeof(int32_t));
				// fpdist.write((char *)&dist, sizeof(uint32_t));
				// bit_push(rangebin, rangeOfs, file);

				update(tr, M, p + 1, 1);
				visited[p] = true;

				++ num;
			}
			++ id;
		}
	}

	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	fpdir.close();

	cout << "end segment array\n";
	// fprintf(stderr, "%u\n", small);
	// if (filebin.n > 0) { 
	// 	fileOfs.write((char*)&filebin.a, sizeof(uint8_t));
	// }
	// fileOfs.close();

	// if (rangebin.n > 0) { 
	// 	rangeOfs.write((char*)&rangebin.a, sizeof(uint8_t));
	// }
	// rangeOfs.close();

	fpdist.close();
	smfpdist.close();
	if (isorder) fporder.close();

	delete[] visited;
	delete[] od;
	delete[] ids;

	// cout << "trees_cnt: " << trees_cnt << "\n";
	cout << "Time of exploring trees = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	mstcom::bsc::BSC_compress(encodestrfn.c_str(), (encodestrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(rootstrfn.c_str(), (rootstrfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(dirbinfn.c_str(), (dirbinfn+".bsc").c_str(), 64);
	mstcom::bsc::BSC_compress(singletonfn.c_str(), (singletonfn+".bsc").c_str(), 64);
	mstcom::lzma::lzma_compress(distfn.c_str(), (distfn+".lzma").c_str());
	mstcom::lzma::lzma_compress(smdistfn.c_str(), (smdistfn+".lzma").c_str(), 1);
	if (isorder) {
		mstcom::lzma::lzma_compress(orderfn.c_str(), (orderfn+".lzma").c_str());
	}

	string tarcmd = "tar -cvf " + outfile + " -C " + folder + " par.txt encodestr.txt.bsc rootstr.txt.bsc dir.bin.bsc singleton.bin.bsc dist.bin.lzma smdist.bin.lzma";

	if (isorder) {
		tarcmd += " order.bin.lzma";
	}

	cout << tarcmd << endl;

	system(tarcmd.c_str());
	cout << "Time of bsc files = " << stopwatch.stop() << std::endl;

}
