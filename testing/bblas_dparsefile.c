/**
 * @file bblas_dparsefile.c
 *
 * @brief Parses input file to setup the current test.
 *
 * BBLAS is a software package provided by Univ. of Manchester,
 * Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author  Samuel  D. Relton
 * @author  Pedro   V. Lara
 * @author  Mawussi Zounon
 * @date    2016-02-20
 *
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @generated from bblas_zparsefile.c normal z -> d, Mon Jun  6 09:44:13 2016
 **/
#endif

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "bblas_dtesting.h"

/** Max length of parameter name **/
#define MAXLEN 80

/**
 * Gets rid of trailing and leading whitespace including
 * the annoying "\n" from fgets().
 **/
char * trim (char * s)
{
  /* Initialize start, end pointers */
  char *s1 = s, *s2 = &s[strlen (s) - 1];

  /* Trim and delimit right side */
  while ( (isspace (*s2)) && (s2 >= s1) )
    s2--;
  *(s2+1) = '\0';

  /* Trim left side */
  while ( (isspace (*s1)) && (s1 < s2) )
    s1++;

  /* Copy finished string */
  memmove(s, s1, strlen(s));
  return s;
}


/*Prototype of bblas_checkconfig */
void bblas_checkconfig(bblas_dtest_t *test);

/**
 * parse_dconfigfile reads the parameter file
 * regardless of the order of the parameters
 * and fills the parameter structure.
 * The parameter file must be opened
 * and closed outside of this function.
 **/

void parse_dconfigfile (FILE *file, bblas_dtest_t *test)
{
  /* Local variables */
  char *s, buff[256];
  int num_value;

  /* Check if file pointer is valid */
  if (file == NULL)
    {
      return;
    }


  /* Read next line */
  while ((s = fgets (buff, sizeof buff, file)) != NULL)
  {
      /* Skip blank lines and comments */
      if (buff[0] == '\n' || buff[0] == '#')
	  continue;

      /* Parse name/value pair from line */
      char name[MAXLEN], value[MAXLEN];
      s = strtok (buff, "=");
      if (s==NULL)
	  continue;
      else
	  strncpy (name, s, MAXLEN);
      s = strtok (NULL, "=");
      if (s==NULL)
	  continue;
      else
	  strncpy (value, s, MAXLEN);
      trim (value);
      trim(name);
      num_value = atoi(value);

      /* Copy into correct entry in parameters struct */
      if (strcmp(name, "gen_uplo")==0)
	  test->gen_uplo = num_value;
      else if (strcmp(name, "gen_transA")==0)
	  test->gen_transA = num_value;
      else if (strcmp(name, "gen_transB")==0)
	  test->gen_transB = num_value;
      else if (strcmp(name, "gen_trans")==0)
	  test->gen_trans = num_value;
      else if (strcmp(name, "gen_side")==0)
	  test->gen_side = num_value;
      else if (strcmp(name, "gen_diag")==0)
	  test->gen_diag = num_value;
      else if (strcmp(name, "minM")==0)
	  test->minM = num_value;
      else if (strcmp(name, "maxM")==0)
	  test->maxM = num_value;
      else if (strcmp(name, "minN")==0)
	  test->minN = num_value;
      else if (strcmp(name, "maxN")==0)
	  test->maxN = num_value;
      else if (strcmp(name, "minK")==0)
	  test->minK = num_value;
      else if (strcmp(name, "maxK")==0)
	  test->maxK = num_value;
      else if (strcmp(name, "minbatch_count")==0)
	  test->minbatch_count = num_value;
      else if (strcmp(name, "maxbatch_count")==0)
	  test->maxbatch_count = num_value;
      else if (strcmp(name, "batch_opts")==0)
	  test->batch_opts = num_value;
      else if (strcmp(name, "routine")==0)
	  test->routine = num_value;
      else if (strcmp(name, "nb_test")==0)
	  test->nb_test = num_value;
      else if (strcmp(name, "target")==0)
	  test->target = num_value;
      else if (strcmp(name, "tolerance")==0)
	  test->tolerance = num_value;
      else if (strcmp(name, "nb_threads")==0)
	  test->nb_threads = num_value;
      else if (strcmp(name, "set_error")==0)
	  test->set_error = num_value;
      else if (strcmp(name, "global_error")==0)
	  test->global_error = num_value;
      else if (strcmp(name, "faulty_iter")==0)
	  test->faulty_iter = num_value;

      else if (strcmp(name, "mkl_sequential")==0)
	  test->mkl_sequential = num_value;
      else if (strcmp(name, "new_accuracy")==0)
	  test->new_accuracy = num_value;
      else
	  printf ("WARNING: %s/%s: Unknown parameter/value pair!\n",
		  name, value);
  }

  /*check configuration */
  bblas_checkconfig(test);
}


/**
 * Checks if parameter values are correct else, exit(0).
 **/
void bblas_checkconfig(bblas_dtest_t *test)
{
    /*Local variables */
    char *filename="bblas_dparsefile.c";

	if ( (test->target == 2) && !(test->routine == 1 || test->routine == 9) )
	{
		bblas_error(filename, "CuBLAS does not implement this function\n");
		exit(0);
	}
	if ( (test->target == 1) && !(test->routine == 1) )
	{
		bblas_error(filename, "MKL does not implement this function\n");
		exit(0);
	}
    if ( test->minM < 0)
    {
		bblas_error(filename, "wrong minM value\n");
		exit(0);
    }
    if ( test->maxM < 0)
    {
		bblas_error(filename, "wrong maxM value\n");
		exit(0);
    }
    if ( test->minN < 0)
    {
		bblas_error(filename, "wrong minN value\n");
		exit(0);
    }
    if ( test->maxN < 0)
    {
		bblas_error(filename, "wrong maxN value\n");
		exit(0);
    }
    if ( test->minK < 0)
    {
		bblas_error(filename, "wrong minK value\n");
		exit(0);
    }
    if ( test->maxK < 0)
    {
		bblas_error(filename, "wrong maxK value\n");
		exit(0);
    }
    if ( test->minbatch_count < 0)
    {
		bblas_error(filename, "wrong minbatch_count value\n");
		exit(0);
    }
    if ( test->maxbatch_count < 0)
    {
		bblas_error(filename, "wrong maxbatch_count value\n");
		exit(0);
    }
    if (!(( test->batch_opts == 0) || ( test->batch_opts == 1)))
    {
		bblas_error(filename, "wrong batch_opts value\n");
		exit(0);
    }
    if ((test->routine < 1) || (test->routine > 9))
    {
		bblas_error(filename, "wrong routine value\n");
		exit(0);
    }
    if (test->nb_test < 0)
    {
		bblas_error(filename, "wrong  nb_test value\n");
		exit(0);
    }

    if ((test->target < 1) || (test->target > 4))
    {
		bblas_error(filename, "wrong  target value\n");
		exit(0);
    }

    /* Variable batch computation is not currently supported in CuBLAS. */
    if ((test->target == 2) && (test->batch_opts == 1))
    {
		bblas_error(filename, "Variable batch option is not implemented by CuBLAS\n");
		exit(0);
    }

    if ((test->set_error !=0) && (test->set_error!= 1))
    {
		bblas_error(filename, "set_error should be 0 or 1 \n");
		exit(0);
    }

    if ((test->set_error) && ((test->faulty_iter < 0) || (test->faulty_iter >= test->nb_test)))

    {
		bblas_error(filename, "Wrong value of faulty_iter, correct value in [0, nb_test] \n");
		exit(0);
    }
}


/**
 * Prints out all the parameters used during the current test.
 **/

void bblas_dprintconfig(bblas_dtest_t *test)
{

  printf(" minM           = %d\n", test->minM);
  printf(" maxM           = %d\n", test->maxM);
  printf(" minN           = %d\n", test->minN);
  printf(" maxN           = %d\n", test->maxN);
  printf(" minK           = %d\n", test->minK);
  printf(" maxK           = %d\n", test->maxK);
  printf(" maxbatch_count = %d\n", test->maxbatch_count);
  printf(" minbatch_count = %d\n", test->minbatch_count);
  printf(" batch_opts     = %d\n", test->batch_opts);
  printf(" bblas routine  = %d\n", test->routine);
  printf(" number of test = %d\n", test->nb_test);

}
