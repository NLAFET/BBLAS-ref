#!/usr/bin/python
##
# @file performance_tests.py
# @author Samuel D. Relton
# @author Pedro V. Lara
# @author Mawussi Zounon
# @date 2016-04-08
#
# @brief Python script to automatically run performance tests.
#
# BBLAS is a software package provided by Univ. of Manchester,
# Univ. of Tennessee.
#
# FILL IN REST LATER
#
##

import os
import re
import sys
import time
import subprocess
from subprocess import call, Popen, PIPE
from optparse import OptionParser

##
# Parser to take all user input and determine which tests should be run.
##
parser = OptionParser()
# Generic options to describe how to code should run
parser.add_option('-p', '--precisions', action='store', dest='precisions',
                  help='Run given precisions (initials, e.g., "sd" for single and double)',
                  default='sdcz')
parser.add_option('-f', '--folder', action='store', dest='outputfolder',
                  help='Folder to store all input (and any output) files (default testinfo).',
                  default='testinfo')
parser.add_option('-o', '--output', action='store', dest='create_output',
                  help='Create output file to store the result (default mytest.csv).',
                  default="mytest.csv")

# Options to specify which routines should be tested.
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

# Options to specify the size of the matrices and batch
parser.add_option('-M', action='store', dest='M',
                  help='Value of M to overwrite both --minM and --maxM (default -1 for unused).',
                  default='2,4,8,16,32,64,128,256,512')
parser.add_option('-N', action='store', dest='N',
                  help='Value of N to overwrite both --minN and --maxN (default -1 for unused).',
                  default='2,4,8,16,32,64,128,256,512')
parser.add_option('-K', action='store', dest='K',
                  help='Value of K to overwrite both --minK and --maxK (default -1 for unused).',
                  default='2,4,8,16,32,64,128,256,512')
parser.add_option('--batchcount', action='store', dest='batchcount',
                  help='Mininum batch count to use, increases up to ' +
                  'nbtest * minbatchcount (default 100).',
                  default='100')

# Parse input and perform error checks.
(opts, args) = parser.parse_args()

## List of all the precisions to be used
precisionlist = list(opts.precisions)

# Get sizes
N = opts.N.split(',')
N = [int(x) for x in N]
M = opts.M.split(',')
M = [int(x) for x in M]
K = opts.K.split(',')
K = [int(x) for x in K]

# See which routine to use
if opts.gemm_routines:
    routine = 1
elif opts.hemm_routines:
    routine = 2
elif opts.her2k_routines:
    routine = 3
elif opts.herk_routines:
    routine = 4
elif opts.symm_routines:
    routine = 5
elif opts.syr2k_routines:
    routine = 6
elif opts.syrk_routines:
    routine = 7
elif opts.trmm_routines:
    routine = 8
elif opts.trsm_routines:
    routine = 9
else:
    print("No routines chosen for testing.\n")
    exit(1)

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

## Blank input file for the test program.
blankinputfile ="""
# Input file automatically generated by run_tests.py

gen_uplo = 1
gen_transA = 1
gen_transB = 1
gen_trans = 1
gen_side = 1
gen_diag = 1
minM = {minM}
maxM = {maxM}
minN = {minN}
maxN = {maxN}
minK = {minK}
maxK = {maxK}
tolerance = 1
nb_test = 1
minbatch_count = {batch_count}
batch_opts = 0
routine = {routine}
target = {target}
new_accuracy = 1
mkl_sequential = {mkl_sequential}
set_error = 0
global_error = 0
faulty_iter = 5
"""

## Count how many tests have been performed so far.
testcounter = 0

## Regex for line of interest
myline = re.compile(r'(.*)PASSED ')

myoutput = "# batchcount = {batchcount}, target = {target}, routine = {routine}\n".format(
            batchcount=opts.batchcount,
            target=targetlist[0],
            routine=routine)
myoutput += "M,N,K,RefPerf,TargetPerf,MinErr,AvgErr,MaxErr,Std\n"

# Run the required tests
for prec in precisionlist:
    ## Skip HEMM/HERK/HER2K for single and double precision
    if (prec == 'd') or (prec == 's'):
        if (routine == 2) or (routine == 3) or (routine == 4):
            continue
        for target in targetlist:
            for sizes in range(len(N)):
                minM = M[sizes]
                maxM = M[sizes]
                minN = N[sizes]
                maxN = N[sizes]
                minK = K[sizes]
                maxK = K[sizes]
                testcounter += 1
                ## Input file with blanks filled in.
                curtest = blankinputfile.format(
                    minM=minM,
                    maxM=maxM,
                    minN=minN,
                    maxN=maxN,
                    minK=minK,
                    maxK=maxK,
                    batch_count=opts.batchcount,
                    routine=routine,
                    target=target,
                    mkl_sequential=opts.mkl_sequential
                )
                ## Name of the current test.
                curtestname = "input" + str(testcounter)
                ## File to write input to.
                f = open(opts.outputfolder + '/' + curtestname, 'w')
                f.write(curtest)
                f.close()
                ## Command to run the test from the input file.
                callcmd = './' + prec + 'test' + ' ' + opts.outputfolder + \
                          '/' + curtestname
                output,error = Popen(callcmd, shell=True,
                                     stdout = subprocess.PIPE,
                                     stderr= subprocess.PIPE).communicate()
                # Parse input to get necessary part and create CSV file
                output = output.split("\n")
                for line in output:
                    if myline.search(line):
                        line = myline.search(line).groups()[0]
                        break
                line = line.split()
                line = [x.replace(")", "") for x in line]
                line = [x.replace("(", "") for x in line]
                while "" in line:
                    line.remove("")
                myoutput += str(minM) + "," + str(minN) + "," + str(minK) + "," + line[2] + \
                            "," + line[4] + "," + line[6] + "," + line[7] + "," + line[8] + \
                            "," + line[9] + "\n"


f = open(opts.outputfolder + '/' + opts.create_output, 'w')
f.write(myoutput)
f.close()
