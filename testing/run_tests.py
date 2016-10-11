#!/usr/bin/python
##
# @file run_tests.py
# @author Samuel D. Relton
# @author Pedro V. Lara
# @author Mawussi Zounon
# @date 2016-04-08
#
# @brief Python script to automatically run tests.
#
# BBLAS is a software package provided by Univ. of Manchester,
# Univ. of Tennessee.
#
# Purpose
# -------
# This script is designed to facilitate the running of multiple tests
# to compare different implementations of the BBLAS standard.
# It provides a command-line interface to run tests using
# multiple routines and parameters,
# obtain the output, and save the results for later investigation.
#
# As a simple example we could run the command
# @code{.sh} ./run_tests.py --precisions s --gemm --alltargets --batch_type f @endcode
# which would run SGEMM with all implementations for a fixed batch using the
# default matrix sizes and batch sizes.
#
# Parameters
# ----------
# @param[in] -p,\--precisions
#            Run tests with the given precisions
#            (initals e.g. -p sd for single and double),
#            Defaults to "sdcz".
#
# @param[in] -q,\--quiet
#            Run tests where all output destined for STDOUT is piped to /dev/null.
#            This is useful for checking that all tests run without errors.
#            Turned off by default.
#
# @param[in] -f,\--folder
#            Folder to store all input (and any output) files.
#            Defaults to "testinfo".
#
# @param[in] -o,\--output
#            Create output files to store the results of the tests.
#            These are stored in the same location as the input files
#            with the extension *.dat.
#            To change the storage location use the -f parameter.
#            Turned off by default.
#
# @param[in] \--allroutines
#            Test all level-3 BLAS routines with the given precisions and parameters.
#            Turned on by default unless specific routines are passed as parameters.
#
# @param[in] \--gemm
#            Run all GEMM related tests with the given precisions and parameters.
#
# @param[in] \--hemm
#            Run all HEMM related tests with the given precisions and parameters.
#
# @param[in] \--her2k
#            Run all HER2K related tests with the given precisions and parameters.
#
# @param[in] \--herk
#            Run all HERK related tests with the given precisions and parameters.
#
# @param[in] \--symm
#            Run all SYMM related tests with the given precisions and parameters.
#
# @param[in] \--syr2k
#            Run all SYR2K related tests with the given precisions and parameters.
#
# @param[in] \--syrk
#            Run all SYRK related tests with the given precisions and parameters.
#
# @param[in] \--trmm
#            Run all TRMM related tests with the given precisions and parameters.
#
# @param[in] \--trsm
#            Run all TRSM related tests with the given precisions and parameters.
#
# @param[in] \--allroutines
#            Test using all available BBLAS implementations.
#            Turned on by default unless specific implementations are passed as parameters.
#
# @param[in] \--other
#            Run tests comparing against other BBLAS implementations.
#            Defaults to an OpenMP version of the reference implementation.
#
# @param[in] \--mkl
#            Run tests comparing against the Intel MKL BBLAS implementation.
#
# @param[in] \--mklsequential
#            Force Intel MKL to run sequentially. 0=False, 1=True.
#            Defaults to 0.
#
# @param[in] \--cublas
#            Run tests comparing against the NVIDIA CuBLAS BBLAS implementation.
#
# @param[in] \--magma
#            Run tests comparing against the MAGMA BBLAS implementation.
#
# @param[in] \--allparams
#            Run tests using all possible combinations of the level-3 BLAS parameters.
#            Turned off by default.
#
# @param[in] \--uplo
#            Value of uplo to be used: 1="U", 2="L", 3=random.
#            Defaults to 1.
#
# @param[in] \--transA
#            Value of transA to be used: 1="N", 2="T", 3="C", 4=random.
#            Defaults to 1.
#
# @param[in] \--transB
#            Value of transB to be used: 1="N", 2="T", 3="C", 4=random.
#            Defaults to 1.
#
# @param[in] \--trans
#            Value of trans to be used: 1="N", 2="T", 3="C", 4=random.
#            Defaults to 1.
#
# @param[in] \--side
#            Value of side to be used: 1="L", 2="R", 3=random.
#            Defaults to 1.
#
# @param[in] \--diag
#            Value of diag to be used: 1="N", 2="U", 3=random.
#            Defaults to 1.
#
# @param[in] \--minM
#            Specifies the minimum value of M to be used.
#            Each subproblem in a generated batch has matrix sizes
#            chosen randomly from the range of M, N, and K given.
#            Defaults to 32. See also the parameter -M.
#
# @param[in] \--maxM
#            Specifies the maximum value of M to be used.
#            Each subproblem in a generated batch has matrix sizes
#            chosen randomly from the range of M, N, and K given.
#            Defaults to 64. See also the parameter -M.
#
# @param[in] \--minN
#            Specifies the minimum value of N to be used.
#            Each subproblem in a generated batch has matrix sizes
#            chosen randomly from the range of M, N, and K given.
#            Defaults to 32. See also the parameter -N.
#
# @param[in] \--maxN
#            Specifies the maximum value of N to be used.
#            Each subproblem in a generated batch has matrix sizes
#            chosen randomly from the range of M, N, and K given.
#            Defaults to 64. See also the parameter -N.
#
# @param[in] \--minK
#            Specifies the minimum value of K to be used.
#            Each subproblem in a generated batch has matrix sizes
#            chosen randomly from the range of M, N, and K given.
#            Defaults to 32. See also the parameter -K.
#
# @param[in] \--maxK
#            Specifies the maximum value of K to be used.
#            Each subproblem in a generated batch has matrix sizes
#            chosen randomly from the range of M, N, and K given.
#            Defaults to 64. See also the parameter -K.
#
# @param[in] -M
#            Value of M to use for all tests.
#            Overwrites the value of `--minM` and `--maxM`.
#
# @param[in] -N
#            Value of N to use for all tests.
#            Overwrites the value of `--minN` and `--maxN`.
#
# @param[in] -K
#            Value of K to use for all tests.
#            Overwrites the value of `--minK` and `--maxK`.
#
# @param[in] -n,\--nbtest
#            Specifies the number of times that each test should be run.
#            Each run increases the size of the batch by `--minbatchcount`.
#            Defaults to 10.
#
# @param[in] \--minbatchcount
#            Specifies the smallest batch size used in the tests.
#            For each test this is increased until it reaches
#            `--nbtest` * `--minbatchcount`.
#            Defaults to 100.
#
# @param[in] \--batch_opts
#            Specifies if tests are performed using fixed or variables batches.
#            Give the initial "f" or "v" to specify fixed or variable, respectively.
#            Defaults to "f".
#
# @param[in] \--tolerance
#            Specifies the tolerance which must be exceeded before an accuracy test
#            counts as a failure.
#            Defaults to 10.
#
# @param[in] \--relativeerror
#            Specifies whether to use accuracy testing based upon the relative error of
#            the computation (as compared with the reference implementation),
#            instead of basing the accuracy tests on rigorous error bounds.
#            Turned off by default.
#
# @param[in] \--seterror
#            Specifies whether to introduce errors into the level-3 BLAS calls,
#            in order to check the error handling of the different implementations.
#            If this is not set then `--globalerror` and `--faultyiter` have no effect.
#            Turned off by default.
#
# @param[in] \--globalerror
#            Specifies whether to introduce an error in the batch_count or batch_opts
#            parameters.
#            If this is not selected then an error will be introduced into one of
#            the other parameters such that only one subproblem out of an entire batch
#            should be problematic.
#            Turned off by default.
#
# @param[in] \--faultyiter
#            Specified which iteration of the test should have a fault injected into it.
#            Must be a value between 0 and nb_test.
#            Set this to -1 to introduce a fault in a random iteration.
#            Defaults to -1.
#
# Note
# ----
# This is based upon a similar script for the MAGMA project by Mark Gates.
#
##

import os
import re
import sys
import time
import subprocess
from subprocess import call
from optparse import OptionParser


##
# Parser to take all user input and determine which tests should be run.
##
parser = OptionParser()
# Generic options to describe how to code should run
parser.add_option('-p', '--precisions', action='store', dest='precisions',
                  help='Run given precisions (initials, e.g., "sd" for single and double)',
                  default='sdcz')
parser.add_option('-q', '--quiet', action='store_true', dest='quiet',
                  help='Send all test output to /dev/null instead of STDOUT. ' +
                  'Useful to check that all tests run without errors (default False).',
                  default=False)
parser.add_option('-f', '--folder', action='store', dest='outputfolder',
                  help='Folder to store all input (and any output) files (default testinfo).',
                  default='testinfo')
parser.add_option('-o', '--output', action='store_true', dest='create_output',
                  help='Create output files to store the results of each test with, ' +
                  'ignores --quiet (default False).',
                  default=False)

# Options to specify which routines should be tested.
parser.add_option('--allroutines', action='store_true', dest='all_routines',
                  help='Run ALL tests with given precisions.',
                  default=False)
parser.add_option('--gemm', action='store_true', dest='gemm_routines',
                  help='Run GEMM related tests with given precisions.',
                  default=False)
parser.add_option('--hemm', action='store_true', dest='hemm_routines',
                  help='Run HEMM related tests with given precisions.',
                  default=False)
parser.add_option('--her2k', action='store_true', dest='her2k_routines',
                  help='Run HER2K related tests with given precisions.',
                  default=False)
parser.add_option('--herk', action='store_true', dest='herk_routines',
                  help='Run HERK related tests with given precisions.',
                  default=False)
parser.add_option('--symm', action='store_true', dest='symm_routines',
                  help='Run SYMM related tests with given precisions.',
                  default=False)
parser.add_option('--syr2k', action='store_true', dest='syr2k_routines',
                  help='Run SYR2K related tests with given precisions.',
                  default=False)
parser.add_option('--syrk', action='store_true', dest='syrk_routines',
                  help='Run SYRK related tests with given precisions.',
                  default=False)
parser.add_option('--trmm', action='store_true', dest='trmm_routines',
                  help='Run TRMM related tests with given precisions.',
                  default=False)
parser.add_option('--trsm', action='store_true', dest='trsm_routines',
                  help='Run TRSM related tests with given precisions.',
                  default=False)

# Options to specify which targets to use.
parser.add_option('--alltargets', action='store_true', dest='all_targets',
                  help='Run tests comparing against all other BBLAS implementations.',
                  default=False)
parser.add_option('--other', action='store_true', dest='other_target',
                  help='Run tests comparing against other BBLAS implementation (default OpenMP).',
                  default=False)
parser.add_option('--mkl', action='store_true', dest='mkl_target',
                  help='Run tests comparing against MKL BBLAS implementation.',
                  default=False)
parser.add_option('--mklsequential', action='store', dest='mkl_sequential',
                  help='Force MKL to run sequentially 0=False 1=True (default 0).',
                  default='0')
parser.add_option('--cublas', action='store_true', dest='cublas_target',
                  help='Run tests comparing against CuBLAS BBLAS implementation.',
                  default=False)
parser.add_option('--magma', action='store_true', dest='magma_target',
                  help='Run tests comparing against MAGMA BBLAS implementation.',
                  default=False)

# Options to specify how the function arguments should be set.
parser.add_option('--allparams', action='store_true', dest='all_params',
                  help='Run tests using all possible combinations of parameters.',
                  default=False)
parser.add_option('--uplo', action='store', dest='uplo',
                  help='Value of UPLO to be used: 1="U" 2="L" 3=random (default 1).',
                  default='1')
parser.add_option('--transA', action='store', dest='transA',
                  help='Value of TRANSA to be used: 1="N" 2="T" 3="C" 4=random (default 1).',
                  default='1')
parser.add_option('--transB', action='store', dest='transB',
                  help='Value of TRANSB to be used: 1="N" 2="T" 3="C" 4=random (default 1).',
                  default='1')
parser.add_option('--trans', action='store', dest='trans',
                  help='Value of TRANS to be used: 1="N" 2="T" 3="C" 4=random (default 1).',
                  default='1')
parser.add_option('--side', action='store', dest='side',
                  help='Value of SIDE to be used: 1="L" 2="R" 3=random (default 1).',
                  default='1')
parser.add_option('--diag', action='store', dest='diag',
                  help='Value of DIAG to be used: 1="N" 2="U" 3=random (default 1).',
                  default='1')

# Options to specify the size of the matrices and batch
parser.add_option('--minM', action='store', dest='minM',
                  help='Minimum value of M to be used (default 32).',
                  default='32')
parser.add_option('--maxM', action='store', dest='maxM',
                  help='Maximum value of M to be used (default 64).',
                  default='64')
parser.add_option('--minN', action='store', dest='minN',
                  help='Minimum value of N to be used (default 32).',
                  default='32')
parser.add_option('--maxN', action='store', dest='maxN',
                  help='Maximum value of N to be used (default 64).',
                  default='64')
parser.add_option('--minK', action='store', dest='minK',
                  help='Minimum value of K to be used (default 32).',
                  default='32')
parser.add_option('--maxK', action='store', dest='maxK',
                  help='Maximum value of K to be used (default 64).',
                  default='64')
parser.add_option('-M', action='store', dest='M',
                  help='Value of M to overwrite both --minM and --maxM (default -1 for unused).',
                  default='-1')
parser.add_option('-N', action='store', dest='N',
                  help='Value of N to overwrite both --minN and --maxN (default -1 for unused).',
                  default='-1')
parser.add_option('-K', action='store', dest='K',
                  help='Value of K to overwrite both --minK and --maxK (default -1 for unused).',
                  default='-1')
parser.add_option('-n', '--nbtest', action='store', dest='nbtest',
                  help='Number of times which each test should run, ' +
                  'increases size of batchcount in each run (default 10).',
                  default='10')
parser.add_option('--batch_opts', action='store', dest='batch_opts',
                  help='Type of batch operation f=FIXED v=VARIABLE (default f).',
                  default='f')
parser.add_option('--minbatchcount', action='store', dest='minbatchcount',
                  help='Mininum batch count to use, increases up to ' +
                  'nbtest * minbatchcount (default 100).',
                  default='100')

# Options to specify the accuracy testing parameters
parser.add_option('--tolerance', action='store', dest='tolerance',
                  help='Tolerance used to decide whether a test passes or fails (default 10).',
                  default='10')
parser.add_option('--relativeerror', action='store_true', dest='relativeerror',
                  help='Use accuracy testing based upon the relative error instead '+
                  'of rigorous error bounds (default False).',
                  default=False)


# Options to specify whether or not to insert errors into the tests.
parser.add_option('--seterror', action='store_true', dest='seterror',
                  help='Introduce an error into the tests to see how they ' +
                  'are handled (default False).',
                  default=False)
parser.add_option('--globalerror', action='store_true', dest='globalerror',
                  help='Introduce an error in batch_count or batch_opts (default False).',
                  default=False)
parser.add_option('--faultyiter', action='store', dest='faultyiter',
                  help='Generate a fault in the given iteration number, ' +
                  'or -1 for random (default -1).',
                  default='-1')

# Parse input and perform error checks.
(opts, args) = parser.parse_args()

# Is faultiter used?
if int(opts.faultyiter) == -1:
    opts.faultyiter = ''

# Parse batch_opts parameter
if opts.batch_opts is 'f':
    opts.batch_opts = 0
elif opts.batch_opts is 'v':
    opts.batch_opts = 1
else:
    print("Invalid argument in parameter --batch_opts.\n")
    exit(1)

# Overwrite M, N, K?
if int(opts.M) > 0:
    opts.minM = opts.M
    opts.maxM = opts.M

if int(opts.N) > 0:
    opts.minN = opts.N
    opts.maxN = opts.N

if int(opts.K) > 0:
    opts.minK = opts.K
    opts.maxK = opts.K

# Create lists of all parameters needed for the tests.

## List of all the precisions to be used
precisionlist = list(opts.precisions)
## List of all routines to be tested.
routinelist = list()

if opts.all_routines:
    routinelist = [1, 2, 3, 4, 5, 6, 7, 8, 9]
else:
    if opts.gemm_routines:
            routinelist.append(1)
    if opts.hemm_routines:
            routinelist.append(2)
    if opts.her2k_routines:
            routinelist.append(3)
    if opts.herk_routines:
            routinelist.append(4)
    if opts.symm_routines:
            routinelist.append(5)
    if opts.syr2k_routines:
            routinelist.append(6)
    if opts.syrk_routines:
            routinelist.append(7)
    if opts.trmm_routines:
            routinelist.append(8)
    if opts.trsm_routines:
            routinelist.append(9)
if len(routinelist) == 0:
    print("No routines chosen for testing. Using all as default.\n")
    routinelist = [1, 2, 3, 4, 5, 6, 7, 8, 9]

## List of all software targets to test against.
targetlist = list()
if opts.all_targets:
    targetlist = [1, 2, 3, 4]
else:
    if opts.mkl_target:
        targetlist.append(1)
    if opts.cublas_target:
        targetlist.append(2)
    if opts.magma_target:
        targetlist.append(3)
    if opts.other_target:
        targetlist.append(4)
if len(targetlist) == 0:
    print("No targets chosen for testing. Using all as default\n")
    targetlist = [1, 2, 3, 4]

## List of all uplo values to use
uplolist = list()
## List of all transA values to use
transAlist = list()
## List of all transB values to use
transBlist = list()
## List of all trans values to use
translist = list()
## List of all side values to use
sidelist = list()
## List of all diag values to use
diaglist = list()

if opts.all_params:
    uplolist = [1, 2]
    transAlist = [1, 2, 3]
    transBlist = [1, 2, 3]
    translist = [1, 2, 3]
    sidelist = [1, 2]
    diaglist = [1, 2]
else:
    uplolist.append(int(opts.uplo))
    transAlist.append(int(opts.transA))
    transBlist.append(int(opts.transB))
    translist.append(int(opts.trans))
    sidelist.append(int(opts.side))
    diaglist.append(int(opts.diag))


## Blank input file for the test program.
blankinputfile ="""
# Input file automatically generated by run_tests.py

gen_uplo = {uplo}
gen_transA = {transA}
gen_transB = {transB}
gen_trans = {trans}
gen_side = {side}
gen_diag = {diag}
minM = {minM}
maxM = {maxM}
minN = {minN}
maxN = {maxN}
minK = {minK}
maxK = {maxK}
tolerance = {tolerance}
nb_test = {nb_test}
minbatch_count = {minbatch_count}
batch_opts = {batch_opts}
routine = {routine}
target = {target}
new_accuracy = {new_accuracy}
mkl_sequential = {mkl_sequential}
set_error = {set_error}
global_error = {global_error}
faulty_iter = {faulty_iter}
"""

## Count how many tests have been performed so far.
testcounter = 0

## Used to create outputfolder if it doesn't yet exist.
mydir = opts.outputfolder
if not os.path.exists(mydir):
        os.makedirs(mydir)

# Run the required tests
for prec in precisionlist:
    for routine in routinelist:
    ## Skip HEMM/HERK/HER2K for single and double precision
        if (prec == 'd') or (prec == 's'):
            if (routine == 2) or (routine == 3) or (routine == 4):
                continue
        for target in targetlist:
            for uplo in uplolist:
                for transA in transAlist:
                    for transB in transBlist:
                        for trans in translist:
                            for side in sidelist:
                                for diag in diaglist:
                                    testcounter += 1
                                    ## Input file with blanks filled in.
                                    curtest = blankinputfile.format(
                                        uplo=uplo,
                                        transA=transA,
                                        transB=transB,
                                        trans=trans,
                                        side=side,
                                        diag=diag,
                                        minM=opts.minM,
                                        maxM=opts.maxM,
                                        minN=opts.minN,
                                        maxN=opts.maxN,
                                        minK=opts.minK,
                                        maxK=opts.maxK,
                                        tolerance=opts.tolerance,
                                        nb_test=opts.nbtest,
                                        minbatch_count=opts.minbatchcount,
                                        batch_opts=opts.batch_opts,
                                        routine=routine,
                                        target=target,
                                        new_accuracy=int(not(opts.relativeerror)),
                                        mkl_sequential=opts.mkl_sequential,
                                        set_error=opts.seterror,
                                        global_error=opts.globalerror,
                                        faulty_iter=opts.faultyiter)
                                    ## Name of the current test.
                                    curtestname = "input" + str(testcounter)
                                    ## File to write input to.
                                    f = open(opts.outputfolder + '/' + curtestname, 'w')
                                    f.write(curtest)
                                    f.close()
                                    ## Command to run the test from the input file.
                                    callcmd = './' + prec + 'test' + ' ' + opts.outputfolder + \
                                              '/' + curtestname
                                    if opts.create_output:
                                        callcmd += '>' + opts.outputfolder + '/' + curtestname + '.dat'
                                    elif opts.quiet:
                                        callcmd += '>/dev/null'
                                    print(callcmd)
                                    call(callcmd, shell=True)
