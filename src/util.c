#include "mstcom.h"

CStopWatch::CStopWatch(){
	running = 0;
}

CStopWatch::~CStopWatch(){
}

void CStopWatch::start(){
	running = 1;
	elapsed.clear();
	t1 = std::chrono::steady_clock::now();
}

double CStopWatch::stop(){
	t2  = std::chrono::steady_clock::now();
	running = 0;
	double el = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0;
	elapsed.push_back(el);
	return el;
}

void CStopWatch::resume(){
	running = 1;
	t1 = std::chrono::steady_clock::now();
}

double CStopWatch::totalTime(){
	double sum = 0;
	for (std::vector<double>::iterator it = elapsed.begin(); it != elapsed.end(); ++it) {
		sum += *it;
	}
	return sum;
}

// inline 
void show_usage(const char* prog) {
	printf("This is mstcom, Reads Compressor. Version v0.1, October 2019.\n");
	printf("Author: Yuansheng Liu <yyuanshengliu@gmail.com>\n");
	printf("\nUsage: %s <e|d> -i <inputfile '*.fastq'> -o <outputfile> <other options>\n", prog);
	printf("\nMstcom options:\n");
	// printf("-----------\n");
	printf("\t -i is a input file, fastq/fasta format file contains reads\n");
	printf("\t -f is the second input file, fastq/fasta format file; only for paired-end reads\n");
	printf("\t -o is the output file\n");
	printf("\t -t is the number of threads, default: 24\n");
	printf("\t -h print help message\n");

	printf("Example:\n");
	printf("\t./mstcom e -i SRR445718.fastq -o SRR445718.mstcom\n");
	printf("\t./mstcom e -i SRR1265495_1.fastq -f SRR1265495_2.fastq -o SRR1265495_pe.mstcom\n");
	printf("\t./mstcom d -i SRR445718.mstcom -o dec.seq\n\n");
}

// static const char alphanum[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
// inline std::string generateString(const std::string &chr, int length = 5) {
// 	srand(time(0));
// 	string res = chr + "_";
// 	for (int i = 0; i < length; ++i) {
// 		res += alphanum[rand() % 62];
// 	}
// 	return res;
// }

// inline 
void getPars(int argc, char* argv[]) {
	nthreads = 24;
	isorder = false;
	ispe = false;
	bool isinfile = false, isoutfile = false; 
	int oc;
	while ((oc = getopt(argc, argv, "i:f:o:t:ph")) >= 0) {
		switch (oc) {
			case 'i':
				infile = optarg;
				isinfile = true;
				break;
			case 'f':
				infile1 = optarg;
				ispe = true;
				// cout << "f: ispe\n";
				break;
			case 'o':
				outfile = optarg;
				isoutfile = true;
				break;
			case 'p':
				isorder = true;
				break;
			case 't':
				nthreads = atoi(optarg);
				break;
			case 'h':
				show_usage(argv[0]);
				exit(0);
			case '?':
				std::cerr << "Error parameters.\n Please run 'mstcom -h'\n";
				exit(1);
				break;
		}
	}

	if (!isinfile || !isoutfile) {
		fprintf(stderr, "Required parameters (input or output file) are not provided!!\n\n");
		exit(1);
	}
	
	std::ifstream f;

	f.open(infile);
	if (f.fail()) {
		fprintf(stderr, "Input file '%s' does not exist.\n", infile.c_str());
		exit(1);
	}
	f.close();

	if (ispe) {
	// fprintf(stderr, "ispe: %d\n", ispe);
		// fprintf(stderr, "%s\n", );
		f.open(infile1);
		if (f.fail()) {
			fprintf(stderr, "Input file '%s' does not exist.\n", infile1.c_str());
			exit(1);
		}
		f.close();
	}

	f.open(outfile);
	if (f.is_open()) {
		fprintf(stderr, "The output file '%s' exist.\n", outfile.c_str());
		exit(1);
	}
	f.close();

	folder = generateString("mstcom", 10); //creat a temp folder for current input

	// sprintf(cmd, "mkdir -p %s", folder.c_str());
	// fprintf(stderr, "%s\n", folder.c_str());
	// system(cmd);
	mode_t mode = 0777;
	int nError = mkdir(folder.c_str(), mode);
	if (nError != 0) {
		fprintf(stderr, "Failed to creat the folder '%s'\n", folder.c_str());
		exit(EXIT_FAILURE);
	}
	folder += "/";
}

// inline 
void reverseComplement(char* start) {
	char* left = start; // sequence starts
	char* right = start + L - 1;
	while (right > left) {
		char tmp = complement[(uint8_t)*left];
		*left = complement[(uint8_t)*right];
		*right = tmp;
		++left;
		--right;
	}
	if (left == right)
		*left = complement[(uint8_t)*left];
}
