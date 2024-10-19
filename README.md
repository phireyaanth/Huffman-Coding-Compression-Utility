This project implements a command-line utility for data compression and decompression using the Huffman coding algorithm. Huffman coding minimizes storage by generating prefix-free codes based on the frequency of symbols in the data. This utility reads raw binary data, compresses it block-by-block, and writes compressed output, with the ability to decompress it back to the original form.

Features
Huffman coding for efficient compression and decompression using prefix-free codes.
Block-based processing: Reads and writes data in configurable block sizes.
Memory-efficient operations: Uses pointer arithmetic and bitwise operations without dynamic memory allocation.
Cross-platform support: Works with standard input and output streams for flexible usage.
Unit-tested with Criterion to ensure correctness and stability.
Optimized build system using Makefiles for seamless debugging and compilation.
