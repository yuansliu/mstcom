# include <cstdio>
# include <iostream>
# include <fstream>
# include <vector>
# include <queue>
# include <assert.h>
# include <chrono>
# include <thread>
# include <string.h>
# include <mutex>
# include <cstdlib>
# include <ctime>
# include <algorithm>
# include <map>
# include <cmath>
# include <unistd.h>
// # include "mstcom.h"

using namespace std;

inline void build(size_t *&tr, size_t &M, size_t n) {
    for (M=1; M<=n+1; M<<=1); 
    tr = new size_t[M + n + 1];
	memset(tr, 0, sizeof(size_t)*(M + n + 1));
}

inline void up(size_t *tr, size_t x) {
	// cout << "up(x): " << x << endl;
    tr[x]=tr[x<<1]+tr[x<<1|1];
}

inline void update(size_t *tr, size_t M, size_t x, int y) {
    for(tr[x+=M]+=y,x>>=1;x;x>>=1)
        up(tr, x);
}

inline size_t query(size_t *tr, size_t M, size_t s, size_t t) {
	bool debug = false;
	if (s == 96677 && t == 96785) {
		// debug = true;
	}
    size_t ans=0;
    s=s+M-1;t=t+M+1;
    if (debug) cout << "s: " << s << "; t: " << t << endl;
    for(;s^t^1;s>>=1,t>>=1)
    {
        if(~s&1) {
        	if (debug) cout << "s^1: " << (s^1) << "; " << tr[s^1] << endl;
        	ans+=tr[s^1];
        }
        if(t&1) {
        	if (debug) cout << "t^1: " << (t^1) << "; " << tr[t^1] << endl;
        	ans+=tr[t^1];
        }
    }
    if (debug) cout << "ans: " << ans << endl;
    return ans;
}

// size_t binarySearch(size_t *tr, size_t M, size_t s, size_t e, bool *flag, size_t dist) { // find a 
inline size_t binarySearch(size_t *tr, size_t M, size_t s, size_t max_rid, bool *flag, size_t dist, bool debug) { // find a 
	size_t m, num, e = max_rid, oe = max_rid, tempe, rdist;
	int rnd = 0;
	while (e > s) {
		num = query(tr, M, s, e);
		if (debug) cout << s << ", " << e << ": " << num << " oe: " << oe << endl;
		// if (dist == e - s + 1 - num) break;
		tempe = e;
		if (e - s + 1 >= num) {
			rdist = e - s + 1 - num;
			if (debug) cout << "rdist: " << rdist << endl;
			if (dist < rdist) {
				// oe = e;
				e = s + ((e - s) >> 1) - 1;
				if (e == s) {
					++ e;
				} else 
				if (e < s) {
					++ e;
					num = query(tr, M, s, e);
					if (dist == e - s + 1 - num) {
						break;
					}
				}
			} else 
			if (dist > rdist) {
				e += ((oe - e) >> 1);
			} else {
				if (debug)	cout << "e: " << e << "; rdist: " << rdist << endl;
				while (e >= s && flag[e - 1]) {
					-- e;
				} 
				if (debug) cout << "xx e: " << endl;
				num = query(tr, M, s, e);
				if (!(e - s + 1 >= num && e - s + 1 - num == dist)) {
					++ e;
				}
				break;
			}
			if (tempe > e) oe = tempe;
		} else {
			e += ((oe - e) >> 1);
		}

		rnd++;
		if (debug) {
			// if (rnd > 30) exit(0);
		}
	}
	if (debug) cout << "s: " << s << ", e: " << e << endl;

	while (e <= max_rid && flag[e - 1]) {
		++ e;
	}
	if (debug) cout << "e: " << e << endl;
	return e;
}

int main(int argc, char const *argv[]) {
	// size_t max_rid = 500000;
	size_t max_rid = 33111272;
	size_t *ids = new size_t[max_rid + 3];
	size_t *od = new size_t[max_rid + 3];

	size_t *tr, M;
	build(tr, M, max_rid);

	bool *flag = new bool[max_rid + 3];
	memset(flag, false, sizeof(bool)*max_rid); 

	size_t half = max_rid >> 1, pid;
	ifstream fporder("order.bin", std::ios::binary);
	FILE *fpdist = fopen("dist.txt", "r");
	uint32_t tempid;
	int32_t dist;

	bool debug = false;
	for (size_t num = 0; num < max_rid; ++num) {
		if (!flag[num]) {
			fporder.read((char*)&tempid, sizeof(uint32_t));
			fscanf(fpdist, "%d", &dist);

			if (tempid == 246295 && dist == -3) {
				// debug = true;
			} else {
				// debug = false;
			}
			// if (num == 12867) debug = true;
			if (debug) 
				cout << num << " " << tempid << " " << dist << endl;

			ids[num] = tempid;
			if (dist < 0) {
				dist = 0 - dist;
				ids[num] += half;
			}
			if (debug) cout << "dist: " << dist << endl;

			pid = binarySearch(tr, M, num + 1 + 1, max_rid, flag, dist, debug);
			// pid = binarySearch(tr, M, num, max_rid, flag, dist, false);

			if (debug) cout << "pid: " << pid << endl;

			if (pid - 1 == 96) {
				cout << "num: " << num << endl;
				// exit(0);
			}

			if (ids[num] < half) {
				ids[pid - 1] = half + ids[num];
			} else {
				ids[pid - 1] = ids[num] - half;
			}

			if (debug) cout << "ids[pid - 1]: " << ids[pid-1] << endl;

			update(tr, M, pid, 1);
			flag[pid - 1] = true;
		}
		// if (num >= 5)  break;
		if (debug)	break;
	}

	FILE *fpres = fopen("decorder.txt", "w");
	for (size_t i = 0; i < max_rid; ++i) {
		fprintf(fpres, "%lu\n", ids[i]);
	}
	fclose(fpres);

	return 0;
}

int mainX(int argc, char const *argv[]) {
	size_t max_rid = 33111272;
	size_t *ids = new size_t[max_rid + 3];
	size_t *od = new size_t[max_rid + 3];

	size_t *tr, M;

	FILE *fp = fopen("idx.txt", "r");
	for (size_t i = 0; i < max_rid; ++i) {
		fscanf(fp, "%lu", &ids[i]);
		od[ids[i]] = i;
		// cout << ids[i] << endl;
		if (ids[i] == 1661232) cout << "i: " << i << endl;
	}
	fclose(fp);

	build(tr, M, max_rid);
	cout << "M: " << M << endl;

	fp = fopen("order1.txt", "w");

	bool *flag = new bool[max_rid + 3];
	memset(flag, false, sizeof(bool)*max_rid); 

	size_t half = max_rid >> 1;
	cout << "half: " << half << endl;

	size_t id = 0, num = 0, p, r;
	while (id < max_rid && num < half) {
		// cout << "id: " << id << " ids[id]: " << ids[id] << endl;
		if (!flag[id]) {
			size_t dis = 0;
			if (ids[id] < half) {
				p = od[half + ids[id]];
				// cout << "yyy " << ids[id] + half << endl;
				r = query(tr, M, id + 1, p - 1);
				fprintf(fp, "%lu\n", p - 1 - (id + 1) + 1 - r);
				// fprintf(stderr, "%lu\n", p - id - r);
				// size_t i = id + 1;
				// while (i < max_rid) {
				// 	if (!flag[i]) {
				// 		// cout << "-ids[i]: " << ids[i] << endl;
				// 		if (ids[i] == half + ids[id]) break;
				// 		++ dis;
				// 	}
				// 	++i;
				// }
				// // ids[i] is the pair of ids[id]
				// // fprintf(fp, "%lu\n", dis);
				// fprintf(stderr, "dis: %lu\n", dis);
			} else {
				p = od[ids[id] - half];
				// cout << "xxx " << ids[id] - half << endl;
				// cout << "p: " << p << endl;

				r = query(tr, M, id + 1, p - 1);
				// cout << "r: " << r << endl;
				// cout << "xxx 111\n";
				// fprintf(fp, " %lu\n", p - id - r);
				fprintf(fp, "-%lu\n", p - 1 - (id + 1) + 1 - r);
				// fprintf(stderr, " %lu\n", p - id - r);
				// cout << "xxx 1112222\n";
				// cout << "xxx 1112222333\n";

				// size_t i = id + 1;
				// while (i < max_rid) {
				// 	if (!flag[i]) {
				// 		++ dis;
				// 		// cout << "ids[i]: " << ids[i] << endl;
				// 		if (ids[i] == ids[id] - half) break;
				// 	}
				// 	++i;
				// }
				// // ids[i] is the pair of ids[id]
				// // fprintf(fp, " %lu\n", dis);
				// fprintf(stderr, " %lu\n", dis);
			}
			update(tr, M, p, 1);
			flag[p] = true;
			++ num;
			// if (num > 300) break;
		}
		++ id;
	} 
	fclose(fp);

	return 0;
}

int mainO(int argc, char const *argv[]) {
	size_t max_rid = 33111272;
	size_t *ids = new size_t[max_rid + 3];
	FILE *fp = fopen("idx.txt", "r");
	for (size_t i = 0; i < max_rid; ++i) {
		fscanf(fp, "%lu", &ids[i]);
		// cout << ids[i] << endl;
	}
	fclose(fp);

	fp = fopen("order.txt", "w");

	bool *flag = new bool[max_rid + 3];
	memset(flag, false, sizeof(bool)*max_rid); 

	size_t half = max_rid >> 1;
	cout << "half: " << half << endl;

	size_t id = 0, num = 0;
	while (id < max_rid && num < half) {
		// cout << "id: " << id << " ids[id]: " << ids[id] << endl;
		if (!flag[id]) {
			size_t dis = 0;
			if (ids[id] < half) {
				size_t i = id + 1;
				while (i < max_rid) {
					if (!flag[i]) {
						// cout << "-ids[i]: " << ids[i] << endl;
						if (ids[i] == half + ids[id]) break;
						++ dis;
					}
					++i;
				}
				// ids[i] is the pair of ids[id]
				fprintf(fp, "%lu\n", dis);
				// fprintf(stderr, "%lu\n", dis);
				flag[i] = true;
			} else {
				size_t i = id + 1;
				while (i < max_rid) {
					if (!flag[i]) {
						++ dis;
						// cout << "ids[i]: " << ids[i] << endl;
						if (ids[i] == ids[id] - half) break;
					}
					++i;
				}
				// ids[i] is the pair of ids[id]
				fprintf(fp, " %lu\n", dis);
				// fprintf(stderr, " %lu\n", dis);
				flag[i] = true;
			}
			++ num;
			// if (num > 10) break;
		}
		++ id;
	} 
	fclose(fp);

	return 0;
}