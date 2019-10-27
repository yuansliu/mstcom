
#ifndef MSTCOM_LIBLZMA_LZMA_H_
#define MSTCOM_LIBLZMA_LZMA_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace mstcom {
namespace lzma {
 
const int BSC_BLOCK_SIZE = 64;  // 64 MB

void lzma_compress(const char *infile, const char *outfile, int _pb = 2);

void lzma_decompress(const char *infile, const char *outfile);

}  // namespace lzma
}  // namespace mstcom

#endif  // 
