/*
 * DO NOT MODIFY THE CONTENTS OF THIS FILE.
 * IT WILL BE REPLACED DURING GRADING
 */
#ifndef HUFF_H
#define HUFF_H

/*
 * Maximum number of "symbols" represented in the Huffman tree used for
 * compression and decompression.  Each possible byte value (0-255)
 * represents a symbol, and in addition there is a special symbol used
 * to signal the end of a data block.
 */
#define MAX_SYMBOLS (256+1)

/*
 * Data to be compressed will be read in "blocks" that consist of a fixed
 * number bytes, called the "block size", which may be specified by the
 * user as a parameter.  All the blocks read by the program will be of the
 * specified block size, except for the last, which might be shorter due
 * to end-of-file being encountered before full block can be read.
 */

/*
 * The smallest block size that may be specified by the user.
 */
#define MIN_BLOCK_SIZE (1024)

/*
 * The largest block size that may be specified by the user.
 */
#define MAX_BLOCK_SIZE (65536)

/*
 * Buffer to hold the block currently being compressed.
 */
unsigned char current_block[MAX_BLOCK_SIZE];

/*
 * Structure of the nodes used for the Huffman tree and for the priority
 * queue used during its construction.
 */
typedef struct node {
    struct node *left;      // Pointer to left child
    struct node *right;     // Pointer to right child
    struct node *parent;    // Pointer to parent
    int weight;             // Weight (priority) of node
    short symbol;           // Symbol at leaf node
} NODE;

/*
 * Array of nodes used for the Huffman tree and priority queue.
 * Since each symbol that occurs in the input will correspond to a leaf
 * of the Huffman tree, and since a binary tree with N leaves has exactly
 * 2*N-1 nodes, we will never need more than 2*MAX_SYMBOLS-1 nodes
 * for the tree.
 */
NODE nodes[2*MAX_SYMBOLS - 1];

/*
 * Number of nodes that are actually in the current Huffman tree.
 * Note that this may be less than all of them, since leaf nodes are only
 * created for symbols that actually occur in an input block.
 * You *must* make sure that your Huffman tree is built in such a way that
 * it occupies a contiguous sequence of nodes in the nodes array,
 * with the root node at index 0, and that you have set the num_nodes
 * variable to correctly indicate the number of nodes in the tree.
 */
int num_nodes;

/*
 * Auxiliary array used to represent a mapping from symbols to nodes.
 * This is used during the reading of an uncompressed data block, in order
 * to determine if a symbol has previously been encountered in the current
 * block, and if so, to find the node that has been assigned to that symbol.
 * It is also used when a Huffman tree is read in during decompression,
 * to keep track of the node that has been assigned to each input symbol so
 * that the actual symbol values can be filled in once the tree structure
 * has been built.
 */
NODE *node_for_symbol[MAX_SYMBOLS];

#endif
