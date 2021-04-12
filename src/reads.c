#include "mstcom.h"

int intnum[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

// inline 
bool getReads(char const *fn) {
	bseq_file_t *fp;
	fp = bseq_open(fn);
	if (fp == 0) return false;
	seq = bseq_read(fp, &RN, L);
	bseq_close(fp);
	return true;
}

// inline 
bool getRightReads(char const *fn) {
	bseq_file_t *fp;
	fp = bseq_open(fn);
	if (fp == 0) return false;
	fprintf(stdout, "Loading the second file %s ...\n", fn);
	size_t leftmaxrid = RN;
	bseq_read_second(seq, fp, &RN, L);
	bseq_close(fp); 
	if (leftmaxrid != (RN>>1)) {
		fprintf(stderr, "Two files contains different number of reads.\n");
		return false;
	}
	return true;
}

bool getRightReads(char const *fn, int LL) {
	bseq_file_t *fp;
	fp = bseq_open(fn);
	if (fp == 0) return false;
	fprintf(stdout, "Loading the second file %s ...\n", fn);
	size_t leftmaxrid = RN;
	bseq_read_second(seq, fp, &RN, LL);
	bseq_close(fp); 
	if (leftmaxrid != (RN>>1)) {
		fprintf(stderr, "Two files contains different number of reads.\n");
		return false;
	}
	return true;
}

// inline 
int getReadsLength(char const *fn) {
	bseq_file_t *fp;
	fp = bseq_open(fn);
	if (fp == 0) return false;
	int readslen = bseq_read_len(fp);
	bseq_close(fp);
	return readslen;
}

#ifdef false
// inline 
int16_t diffstrlen(char *parent, char *child, const int16_t &b) {
	int eq_char_num = 0;
	int16_t dif = 0;
	// if (debug) {
	// 	cout << "b: " << b << endl;
	// 	cout << parent << endl;
	// 	cout << child << endl;
	// }
	// AAATGC --> parent
	//   ATGCAA --> child
	if (b >= 0) {
		dif = b;
		for (int i = b, j = 0; i < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				if (eq_char_num > 1) {
					dif += intnum[eq_char_num];
					// dif += 2;
					eq_char_num = 0;
				} else {
	                dif += eq_char_num;
	                eq_char_num = 0;
				}
				++ dif;
			} else ++ eq_char_num;
		}
	} else {
		dif = -b;
		for (int i = -b, j = 0; i < L; ++i, ++j) {
			if (child[i] != parent[j]) {
				if (eq_char_num > 1) {
					dif += intnum[eq_char_num];
					// dif += 2;
					eq_char_num = 0;
				} else {
	                dif += eq_char_num;
	                eq_char_num = 0;
				}
				++ dif;
			} else ++ eq_char_num;
		}
	}
	return dif;
}

// inline 
int16_t diffstrlen(char *parent, char *_child, const int16_t &b, const bool &isrc) {
	char *child = (char*)alloca((L + 1) * sizeof(char));
	strcpy(child, _child);
	if (isrc) {
		reverseComplement(child);
	}
	int eq_char_num = 0;
	int16_t dif = 0;
	// if (debug) {
	// 	cout << "b: " << b << endl;
	// 	cout << parent << endl;
	// 	cout << child << endl;
	// }
	if (b >= 0) {
		dif = b;
		for (int i = b, j = 0; i < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				if (eq_char_num > 1) {
					dif += intnum[eq_char_num];
					// dif += 2;
					eq_char_num = 0;
				} else {
	                dif += eq_char_num;
	                eq_char_num = 0;
				}
				++ dif;
			} else ++ eq_char_num;
		}
	} else {
		dif = -b;
		for (int i = -b, j = 0; i < L; ++i, ++j) {
			if (child[i] != parent[j]) {
				if (eq_char_num > 1) {
					dif += intnum[eq_char_num];
					// dif += 2;
					eq_char_num = 0;
				} else {
	                dif += eq_char_num;
	                eq_char_num = 0;
				}
				++ dif;
			} else ++ eq_char_num;
		}
	}
	return dif;
}
#endif
// inline 
void encode(char *parent, char *child, const int16_t &_shift, char *en_str) { //str1 is the short one string
	// char *en_str = (char*)alloca(((L<<2) + 1) * sizeof(char));
	// en_str = (char*)calloc(L<<1, sizeof(char));
	char *int_str = (char*)alloca(20 * sizeof(char));
	int16_t shift = _shift;
	bool debug = false;
	if (strcmp(child, "GGTCTCAAACTGCTGACTTCAAGTGATCTGCCCGCCTTGGCCTCCCAAAGTGCTGAGATTACGGATGTGAGCCACTGTGCCCAAATTTTTTTTTTTTTTTT") == 0) {
		// debug = true;
	}
	int en_str_len = 0;
	int eq_char_num = 0;

	if (shift >= 0) {
		// case shift >= 0
		//    AATTGCATGC parent
		//      TTGCATGCGA
		int i, j;
		for (i = shift, j = 0; i < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				if (eq_char_num > 1) {
					sprintf(int_str, "%d", eq_char_num);
					for (char *tk = int_str; *tk != '\0'; ++tk) {
						en_str[en_str_len++] = *tk;
					}
					eq_char_num = 0;
				} else {
					for (int i0 = j - eq_char_num; i0 < j; ++i0) {
	                    en_str[en_str_len++] = child[i0];
	                }
	                eq_char_num = 0;
				}
				en_str[en_str_len++] = child[j]; 
			} else ++ eq_char_num;
		}
		if (shift > 0)
		en_str[en_str_len++] = ' ';
		while (j < L) {
			en_str[en_str_len++] = child[j++];
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
		en_str[en_str_len++] = ' ';
		int spidx = en_str_len - 1;
		for (i = 0, j = shift; j < L; ++i, ++j) {
			if (parent[i] != child[j]) {
				if (eq_char_num > 1) {
					sprintf(int_str, "%d", eq_char_num);
					for (char *tk = int_str; *tk != '\0'; ++tk) {
						en_str[en_str_len++] = *tk;
					}
					eq_char_num = 0;
				} else {
					for (int i0 = j - eq_char_num; i0 < j; ++i0) {
	                    en_str[en_str_len++] = child[i0];
	                }
	                eq_char_num = 0;
				}
				en_str[en_str_len++] = child[j]; 
			} else ++ eq_char_num;	
		}
		en_str[en_str_len] = '\0';

		// ???
		if (spidx > 0 && en_str_len - 1 != spidx) {
			bool isdig = false;
			for (int i = 0; i < en_str_len; ++i) {
				if (en_str[i] >= '0' && en_str[i] <= '9') {
					isdig = true;
					break;
				}
			}
			if (!isdig) {
				en_str[en_str_len++] = '0';
				en_str[en_str_len] = '\0';
				// cout << "parent: " << parent << endl;
 			// 	cout << "child: " << child << endl;
 			// 	cout << en_str << endl;
			}
		}
 	}
 	// if (strcmp(parent, "ACTACCAGACTTCCTGTGAGTTCTTGAGCCATAGCTCCAGAACATTCAGGAAAGATCCATGTTTGCTTTCTCTTCTTTCTTTCTTTATTTTGTTTTTGAGA") == 0
 	// if (strcmp(child, "GGTCTCAAACTGCTGACTTCAAGTGATCTGCCCGCCTTGGCCTCCCAAAGTGCTGAGATTACGGATGTGAGCCACTGTGCCCAAATTTTTTTTTTTTTTTT") == 0) {
 	if (debug) {
 		cout << "shift: " << shift << endl;
 		cout << "parent: " << parent << endl;
 		cout << "child: " << child << endl;
 		cout << en_str << endl;
 	}
	// en_str[en_str_len] = '\0';
	// cout << en_str << endl;
	// cout << "dif: " << dif << "\n";
	// return dif;
	// return string(en_str);
	// return en_str_len;
}

