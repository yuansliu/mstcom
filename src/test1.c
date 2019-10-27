# include <cstdio>
# include <iostream>
# include <fstream>
# include <cstring>
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
			// cout << 1 << endl;
			bit1_push32(v, fp);
		} else {
			// cout << 0 << endl;
			bit0_push32(v, fp);
		}
		a >>= 1;
	}
}

int mainx(int argc, char const *argv[]) {
	uint32_t max_rid = 33111272;//25000;
	// uint32_t max_rid = 500000;
	int kbits = preCalc(max_rid);

	FILE *fp = fopen("dist.txt", "r");
	int tempdist;
	uint32bit_v distbin;
	distbin.a = distbin.n = 0;

	std::ofstream fpdist("dist.bin", std::ios::binary);

	int num = 0;
	cout << "kbits: " << kbits << endl;

	while (fscanf(fp, "%d", &tempdist) != EOF) {
		//
		if (tempdist >= 0) {
			// cout << 1 << endl;
			bit1_push32(distbin, fpdist);
		} else {
			// cout << 0 << endl;
			bit0_push32(distbin, fpdist);
			tempdist = 0 - tempdist;
		}

		outKbits(tempdist, kbits, distbin, fpdist);

		++ num;
		// if (num >= 10) break;

		-- dignum[kbits];
		if (dignum[kbits] <= 0) {
			-- kbits;
		}
	}
	fclose(fp);
	// cout << ": " << distbin.a << endl;
	if (distbin.n > 0) {
		fpdist.write((char*)&distbin.a, sizeof(uint32_t));
	}
	fpdist.close();

	return 0;
}

int mainX() {
	int s = 3, t = 1;
	// fprintf(stderr, "%d\n", s^t^1);
	fprintf(stderr, "%d\n", RAND_MAX);
}

// 256+256+512+1024+2048+4096+8192+8616

void split2TwoEncodedStr(char *str, char *ori, char *aft) {
	char *s = str;
	bool res = false;
	*aft = '\0';
	for (; (*ori = *s) != '\0'; ++ori, ++s) {
		if (*s == '|') {
			*ori = '\0';
			++ s;
			for (; (*aft = *s) != '\0'; ++aft, ++s);
			break;
		}
	}
}

int mainccc() {
	char a[] = "ab12N|abN";
	char ori[100], aft[100];
	split2TwoEncodedStr(a, ori, aft);
	cout << ori << "; " << strlen(ori) << endl;
	cout << aft << "; " << strlen(aft) <<  endl;
}

int main() {
	int a[] = {3, 4, 10, 9, 7, 8, 6, 5, 4, 3, 2, 1};
	for (int i = 0; i < 12; ++i) {
		cout << a[i] << " ";
	}
	cout << endl;
	sort (a + 2, a + 2 + 5);
	for (int i = 0; i < 12; ++i) {
		cout << a[i] << " ";
	}
	cout << endl;
}

