/*
 * DO NOT MODIFY THE CONTENTS OF THIS FILE.
 * IT WILL BE REPLACED DURING GRADING
 */
#ifndef CONST_H
#define CONST_H

#include <stdio.h>

#include "huff.h"

#define USAGE(program_name, retcode) do{ \
fprintf(stderr, "USAGE: %s %s\n", program_name, \
"[-h] [-c|-d] [-b BLOCKSIZE]\n" \
"    -h       Help: displays this help menu.\n" \
"    -c       Compress: read raw data, output compressed data\n" \
"    -d       Decompress: read compressed data, output raw data\n" \
"    -b       For compression, specify blocksize in bytes (range [1024, 65536])\n"); \
exit(retcode); \
} while(0)

/* Options info, set by validargs. */
int global_options;

/*
 * The following are functions that you are required to implement.
 * Refer to the assignment handout and the documentation associated with the
 * stubs in huff.c for detailed specification.
 */

void emit_huffman_tree();
int read_huffman_tree();
int compress_block();
int decompress_block();
int compress();
int decompress();

int validargs(int argc, char **argv);

#endif
