#! /usr/bin/env python
# -*- coding: utf-8 -*-

##
# @file netlib.py
#
# @brief Holds parameters required for compiling and linking BBLAS.
#
# BBLAS is a software package provided by Univ. of Manchester,
# Univ. of Tennessee.
#
# @version 1.0.0
# @author Julie Langou
# @author Mathieu Faverge
# @author Samuel D. Relton
# @date 2016-04-14
#
##

##
# This class holds parameters required to compile and link BBLAS with BLAS, CBLAS,
# LAPACK and LAPACKE.
##
class Config:
  compiler    = "GNU"
  blasname    = "Unknown"
  blaslib     = ""
  cblasdir    = ""                # the CBLAS library
  cblaslib    = ""
  lapackdir   = ""                # the Lapack library
  lapacklib   = ""
  lapcdir     = ""                # the Lapack C interface
  lapclib     = ""
  cudadir     = ""
  magmadir    = ""
  lapackinc   = ""                # the Lapack headers
  cc          = "gcc"             # the C compiler for BBLAS
  fc          = "gfortran"        # the Fortran compiler for core_lapack
  ranlib      = ""                # Ranlib
  arflags     = "rc"              # ar flags
  ldflags_c   = ""                # loader flags when main program is in C
  ldflags_fc  = ""                # loader flags when main program is in Fortran
  ld_fcmain   = ""                # the option to link C main with fortran linker
  withf90     = 1                 # Compile the f90 interface
  ccflags     = "-O2"
  fcflags     = "-O2"
  noopt       = "-O0"
  lapinstalled = 0

  def __init__(self, version):
    self.version = version
