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

int preCalc(uint32_t max_rid) {
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

typedef struct { int n; uint32_t a; } uint32bit_v;
#define bit1_push32(v, fp) do {									\
		(v).a += (1U << ((v).n));										\
		++ (v).n;										\
		if ((v).n == 32) {										\
			fp.write((char*)&(v).a, sizeof(uint32_t));							\
			(v).a = (v).n = 0;							\
		}															\
	} while (0)

#define bit0_push32(v, fp) do {									\
		++ (v).n;										\
		if ((v).n == 32) {										\
			fp.write((char*)&(v).a, sizeof(uint32_t));							\
			(v).a = (v).n = 0;							\
		}															\
	} while (0)

void outKbits(uint32_t a, int kbits, uint32bit_v &v, ofstream &fp) {
	for (int i = 0; i < kbits; ++i) {
		if (a&1) {
			cout << 1 << endl;
			bit1_push32(v, fp);
		} else {
			cout << 0 << endl;
			bit0_push32(v, fp);
		}
		a >>= 1;
	}
}

int main(int argc, char const *argv[]) {
	// uint32_t max_rid = 207579467;//25000;
	uint32_t max_rid = 250000;
	int kbits = preCalc(max_rid);

	FILE *fp = fopen("dist.txt", "r");
	int tempdist;
	uint32bit_v distbin;
	distbin.a = distbin.n = 0;

	std::ofstream fpdist("dist.bin", std::ios::binary);

	while (fscanf(fp, "%d", &tempdist) != EOF) {
		//
		if (tempdist >= 0) {
			cout << 1 << endl;
			bit1_push32(distbin, fpdist);
		} else {
			cout << 0 << endl;
			bit0_push32(distbin, fpdist);
			tempdist = 0 - tempdist;
		}

		outKbits(tempdist, kbits, distbin, fpdist);

		return 0;
		-- dignum[kbits];
		if (dignum[kbits] <= 0) {
			-- kbits;
		}
	}
	fclose(fp);
	if (distbin.n > 0) {
		fpdist.write((char*)&distbin.a, sizeof(uint32_t));
	}
	fpdist.close();

	return 0;
}

// 256+256+512+1024+2048+4096+8192+8616
