# include <cstdio>
# include <iostream>
# include <fstream>
# include <vector>
# include <queue>
# include <assert.h>
# include <chrono>
# include <thread>
# include <mutex>
# include <cstdlib>
# include <ctime>
# include <algorithm>
# include <map>
# include <unistd.h>
using namespace std;

uint32_t dignum[64];

int preCalcx(uint32_t max_rid) {
	int reskbits = 8;
	dignum[8] = 256;
	uint32_t total = 256;
	max_rid -= 256;
	for (int i = 9; i <= 32; ++i) {
		dignum[i] = (1U << i) - (1U << (i-1));

		// cout << dignum[i] << "; " << total << endl;
		if (max_rid < dignum[i]) {
			dignum[i] = max_rid;
			reskbits = i;
			break;
		}
		max_rid -= dignum[i];
	}
	for (int i = 8; i <= reskbits; ++i) {
		cout << dignum[i] << endl;
	}
	return reskbits;
}

int preCalc(uint32_t max_rid) {
	int reskbits = 8;
	dignum[8] = 256;
	uint32_t total = 256;
	for (int i = 9; i <= 32; ++i) {
		dignum[i] = (1U << i) - (1U << (i-1));

		// cout << dignum[i] << "; " << total << endl;
		if (max_rid <= (1U << i)) {
			cout << (1U << (i-1)) << "; " << max_rid << endl;
			dignum[i] = max_rid - (1U << (i-1));
			reskbits = i;
			break;
		}
		// max_rid -= dignum[i];
	}
	for (int i = 8; i <= reskbits; ++i) {
		cout << i << "; " << dignum[i] << endl;
	}
	return reskbits;
}

typedef struct { int n; uint32_t a; } uint32bit_v;

// 256+256+512+1024+2048+4096+8192+8616

struct BIT32POOL {
	int n, a[32];
};

inline int getOneBit(std::ifstream& fpdist, BIT32POOL& curpool) {
	if (curpool.n >= 32) {
		uint32_t bin;
		fpdist.read((char*)&bin, sizeof(uint32_t));
		// cout << "bin: " << bin << endl;
		for (int i = 0; i < 32; ++i) {
			curpool.a[i] = bin&1;
			bin >>= 1;
		}
		curpool.n = 0;
	}
	int res = curpool.a[curpool.n];
	++curpool.n;
	// cout << res << endl;
	return res;
}

inline uint32_t getDist(std::ifstream& fpdist, BIT32POOL& curpool, int kbits) {
	uint32_t res = 0;
	int bits;
	for (int i = 0; i < kbits; ++i) {
		bits = getOneBit(fpdist, curpool);
		// res <<= 1;
		// cout << bits << endl;
		if (bits) {
			res += (1U << i);
			// res += 1;
		}
	}
	// cout << "--\n";
	return res;
}

int main(int argc, char const *argv[]) {
	uint32_t max_rid = 500000;
	int kbits = preCalc(max_rid);
	std::ifstream fpdist("dist.bin", std::ios::binary);

	FILE *fp = fopen("recdist.txt", "w");
	BIT32POOL curpool;
	curpool.n = 32;

	// cout << "kbits: " << kbits << endl;

	int num = 0;
	for (uint32_t i = 0; i < (max_rid >> 1); ++i) {
		int bits = getOneBit(fpdist, curpool);
		uint32_t dist = getDist(fpdist, curpool, kbits);
		if (bits == 0) fprintf(fp, "-");
		fprintf(fp, "%u\n", dist);
		// cout << "dist: " << dist << endl;

		++ num;
		// if (num >= 10)
		// return 0;
		-- dignum[kbits];
		if (dignum[kbits] <= 0) {
			-- kbits;
		}
	}
	fpdist.close();
	fclose(fp);
	return 0;
}