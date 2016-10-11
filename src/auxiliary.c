/**
 *
 * @file auxiliary.c
 *
 *  @brief BBLAS auxiliary routines, should be unused.
 *
 *  BBLAS is a software package provided by Univ. of Tennessee,
 *  Univ. of Manchester Univ. of California Berkeley
 *  and Univ. of Colorado Denver
 *
 * @version 1.0.0
 * @date 2016-02-20
 *
 **/
#include <stdio.h>
#include <stdlib.h>
#include "auxiliary.h"

/**
 *
 *  Indicates a recoverable problem.
 *  User's erroneous action without severe consequences.
 *  Problems occuring while BBLAS is being used correctly.
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 **/
void bblas_warning(const char *func_name, char* msg_text)
{

    fprintf(stderr, "BBLAS WARNING: %s(): %s\n", func_name, msg_text);
}

/**
 *
 *  Indicates a recoverable problem.
 *  User's erroneous action with potentially severe consequences.
 *  Problems occuring due to incorrect use of BBLAS.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 **/
void bblas_error(const char *func_name, char* msg_text)
{

    fprintf(stderr, "BBLAS ERROR: %s(): %s\n", func_name, msg_text);
}

/**
 *
 *  Unexpected behavior within the library.
 *  Unrecoverable user errors.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 **/
void bblas_fatal_error(const char *func_name, char* msg_text)
{
    fprintf(stderr, "BBLAS FATAL ERROR: %s(): %s\n", func_name, msg_text);
    exit(0);
}

/**
 * Check whether a malloc has been performed correctly.
 **/

void bblas_malloc_check(void *ptr, char* msg_text)
{
    if( ptr == NULL)
    {
	fprintf(stderr, "BBLAS MALLOC ERROR: %s\n", msg_text);
	exit(0);
    }
}
