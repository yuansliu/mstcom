# include <cstdio>
# include <iostream>
# include <fstream>
# include <string>
# include <cstring>
# include <cstdlib>
using namespace std;

bool checkStrOrDig(char *str) {
	int len = strlen(str);
	// if (str[0] == ' ' && str[1] == '\n') { // ' \n'
	if (str[0] == '\n') { // ' \n'
		return true;
	}
	for (int i = 0; i < len - 1; ++i) {
		if (str[i] >= 'A' && str[i] <= 'Z') {
			return true;
		}
	}
	return false;
}


int main(int argc, char const *argv[]) {
	
	char *enstr = new char[1024];
	FILE *fp = fopen(argv[1], "r");

	FILE *p0fp = fopen("part0.txt", "w");
	FILE *p1fp = fopen("part1.txt", "w");
	
	while (fgets(enstr, 1024, fp) != NULL) {
		if (enstr[0] == '-') {
			fprintf(p0fp, "-\n");
			continue;
		}

		if (checkStrOrDig(enstr)) { // is encode string
			fprintf(p0fp, "%s", enstr);
		} else { // back
			fprintf(p1fp, "%s", enstr);
		}
	}
	fclose(p0fp);
	fclose(p1fp);
	fclose(fp);
	return 0;
}

int mainX(int argc, char const *argv[]) {
	
	char *enstr = new char[1024];
	FILE *fp = fopen(argv[1], "r");

	FILE *p0fp = fopen("part3.txt", "w");
	FILE *p1fp = fopen("part2.txt", "w");
	
	while (fgets(enstr, 1024, fp) != NULL) {
		if (enstr[0] == '-') {
			fprintf(p0fp, "-\n");
			continue;
		} else {
			fprintf(p1fp, "%s", enstr);
		}
	}
	fclose(p0fp);
	fclose(p1fp);
	fclose(fp);
	return 0;
}
