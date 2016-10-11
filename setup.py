#! /usr/bin/env python
# -*- coding: utf-8 -*-

##
#
# @file setup.py
#
# @brief Installs linear algebra dependencies automatically.
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
# The script will try to figure out some of the features of your system
# and the location of few other tools required for
# the installation of the BBLAS package and the other software
# packages that it requires.
# It is strongly recommended that you help the script by providing info
# through the flags listed below.
#
# @param
# \--help,-h          : display this help and exit
#
# @param
# \--prefix=[DIR]      : install files in DIR [./install]
#
# @param
# \--build=[DIR]       : libraries are built in DIR [./build]
#                       Contains log, downloads and builds.
#
# @param
# \--cc=[CMD]          : the C compiler. [gcc]
#
# @param
# \--fc=[CMD]          : the Fortran compiler. [gfortran]
#
# @param
# \--cflags=[FLAGS]    : the flags for the C compiler [-02]
#
# @param
# \--fflags=[FLAGS]    : the flags for the Fortran compiler [-O2]
#
# @param
# \--ldflags_c=[flags] : loader flags when main program is in C. Some
#                       compilers (e.g. PGI) require different
#                       options when linking C main programs to
#                       Fortran subroutines and vice-versa
#
# @param
# \--ldflags_fc=[flags]: loader flags when main program is in
#                       Fortran. Some compilers (e.g. PGI) require
#                       different options when linking Fortran main
#                       programs to C subroutines and vice-versa.
#                       If not set, ldflags_fc = ldflags_c.
#
# @param
# \--make=[CMD]        : the make command [make]
#
# @param
# \--blaslib=[LIB]     : a BLAS library
#                       (path should be absolute if \--prefix is used)
#
# @param
# \--cblaslib=[LIB]    : a CBLAS library
#                       (path should be absolute if \--prefix is used)
#
# @param
# \--lapacklib=[LIB]   : a Lapack library
#                       (path should be absolute if \--prefix is used)
#
# @param
# \--lapclib=[LIB]     : path to a LAPACK C interface.
#                        (path should be absolute if \--prefix is used)
#
# @param
# \--cudadir=[DIR]     : path to CUDA installation folder e.g. "/usr/local/cuda".
#                        (If not used then software is compiled without CUDA).
#
# @param
# \--downblas          : Download and install reference BLAS.
#
# @param
# \--downcblas         : Download and install reference CBLAS.
#
# @param
# \--downlapack        : Download and install reference LAPACK.
#
# @param
# \--downlapc          : Download and install reference LAPACK C Interface.
#
# @param
# \--downall           : Download and install all missing external libraries.
#                       If you don't have access to wget or no network
#                       connection, you can provide the following packages
#                       in the directory build/download:
#                       http://netlib.org/blas/blas.tgz
#                       http://www.netlib.org/blas/blast-forum/cblas.tgz
#                       http://www.netlib.org/lapack/%s.tgz
#
# @param
# \--[no]testing       : enables/disables the testings. All externals
#                       libraries are required and tested if enabled.
#                       Enabled by default.
#
# @param
# \--[no]documentation : enables/disables the generation of documentation.
#                       Requires doxygen and graphviz to be installed.
#                       Enabled by default.
#
# @param
# \--disable-f90       : to disable the compilation of the f90 interface.
#
# @param
# \--clean             : cleans up the installer directory.
#
# @param
# \--src               : Generates a make.inc for BBLAS developers with
#                       options given. If some external libraries are
#                       not available, they are automatically downloaded
#                       and installed in the prefix directory.
#                       Testing step is also deactivated by default.
#
#
# University of Tennessee ICL License
#
# \-- Innovative Computing Laboratory
# \-- Electrical Engineering and Computer Science Department
# \-- University of Tennessee
# \-- (C) Copyright 2008-2010
#
# Redistribution  and  use  in  source and binary forms, with or without
# modification,  are  permitted  provided  that the following conditions
# are met:
#
# * Redistributions  of  source  code  must  retain  the above copyright
#   notice,  this  list  of  conditions  and  the  following  disclaimer.
# * Redistributions  in  binary  form must reproduce the above copyright
#   notice,  this list of conditions and the following disclaimer in the
#   documentation  and/or other materials provided with the distribution.
# * Neither  the  name of the University of Tennessee, Knoxville nor the
#   names of its contributors may be used to endorse or promote products
#   derived from this software without specific prior written permission.
#
# THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# ``AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
# LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
##

import sys
##
# Major version number
##
VERSION_MAJOR = 1
##
# Minor version number
##
VERSION_MINOR = 0
##
# Micro version number
##
VERSION_MICRO = 0

from script.blas        import Blas
from script.cblas       import CBlas
from script.lapack      import Lapack
from script.tmg         import Tmg
from script.lapcwrapper import Lapcwrapper
from script.bblas       import BBLAS

import script.netlib as netlib

##
# Start installation script
##
def main(argv):

  ### Store history of executed commands in config.log
  cmd = ""
  for arg in argv:
      cmd += arg+" "
  cmd += "\n"
  fp = open("history.log",'a')
  fp.write(cmd)
  fp.close()
  ### END

  config = netlib.Config((VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO))

  bblas = BBLAS(argv, config)

  if bblas.testing or bblas.src or bblas.downblas :
    Blas(config, bblas)

  if bblas.testing or bblas.src or bblas.downcblas :
    CBlas(config, bblas)

  if bblas.testing or bblas.src or bblas.downlapack :
    Lapack(config, bblas)

  # bblas.downtmg set to 1 by lapack if necessary
  if bblas.needtmg :
    Tmg(config, bblas)

  # Always required for the lapack.h
  Lapcwrapper(config, bblas)

  bblas.resume()

  return 0

if "__main__" == __name__:
  sys.exit(main(sys.argv))
