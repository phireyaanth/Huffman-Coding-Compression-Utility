#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "debug.h"

int main(int argc, char **argv)
{
    int ret;
    if(validargs(argc, argv))
        USAGE(*argv, EXIT_FAILURE);
    debug("Options: 0x%x", global_options);
    if(global_options & 1)
        USAGE(*argv, EXIT_SUCCESS);

    //call compress or decompress depending on global_options
    int error = 0;
    if(global_options & 2)
        error = compress();
    else if(global_options & 4)
        error = decompress();

    //if compress or decompress has an error, EXIT_FAILURE
    if(error == -1)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

/*
 * Just a reminder: All non-main functions should
 * be in another file not named main.c
 */
