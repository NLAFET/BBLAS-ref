#! /usr/bin/env python
# -*- coding: utf-8 -*-

##
#
# @file blas.py
#
# @brief Installs Netlib BLAS if optimized version is not provided.
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
#
# University of Tennessee ICL License
#
# -- Innovative Computing Laboratory
# -- Electrical Engineering and Computer Science Department
# -- University of Tennessee
# -- (C) Copyright 2008-2010
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

from utils import writefile, runShellCommand, killfiles, downloader, getURLName
import sys
import os
import shutil
import re

##
# This class is responsible for testing, downloading, and verifying that
# the system has a suitable version of BLAS for compiling and running BBLAS.
##
class Blas:
    blasurl      = "http://netlib.org/blas/blas.tgz"

    """ This class takes care of the BLAS installation. """
    def __init__(self, config, bblas):
        print "\n","="*40
        print "  BLAS installation/verification"
        print "="*40

        self.config  = config
        self.downcmd = bblas.downcmd
        self.prefix  = bblas.prefix
        downblas = bblas.downblas;

        if config.blaslib:
            self.check_blas()
        else:
            if downblas:
                self.down_install_blas()
            else:
                print """
Please provide a working BLAS library. If a BLAS library
is not present on the system, the reference BLAS library it can be
automatically downloaded and installed by adding the --downblas flag.
Be aware that a reference BLAS library will be installed with the --downblas
flag so don't expect good performance.
You can specify the optimized BLAS library **installed** on your machine using the --blaslib option
For instance
    --blaslib="-lessl" for IBM BLAS,
    --blaslib="-framework veclib" for Mac OS/X,
    --blaslib="-lgoto" for GOTO BLAS
    --blaslib="-lf77blas -lcblas -latlas" for ATLAS
      For MKL we recommend using the link line generated by the MKL Link Line Adviser:
      https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor.
      For instance:
    --blaslib=" -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -lpthread -lm -ldl"
                """
                sys.exit()

    def check_blas(self):

        print "Checking if provided BLAS works...",
        # This function simply generates a FORTRAN program
        # that contains few calls to BLAS routine and then
        # checks if compilation, linking and execution are succesful

        # Try to detect which BLAS is used
        if re.search('mkl', self.config.blaslib, re.IGNORECASE):
            self.config.blasname = "mkl"
        elif re.search('acml', self.config.blaslib, re.IGNORECASE):
            self.config.blasname = "acml"

        if self.config.blasname != "Unknown":
            if self.config.compiler == "Intel":
                self.config.ccflags    += " -openmp"
                self.config.ldflags_c  += " -openmp"
                self.config.ldflags_fc += " -openmp"
            elif self.config.compiler == "GNU":
                self.config.ccflags    += " -fopenmp"
                self.config.ldflags_c  += " -fopenmp"
                self.config.ldflags_fc += " -fopenmp"

        sys.stdout.flush()
        writefile('tmpf.f',"""
      program ftest
      double precision da, dx(1)
      dx(1)=1
      da = 2
      call dscal(1,da,dx,1)
      stop
      end\n""")

        fcomm = self.config.fc+' -o tmpf '+'tmpf.f '+self.config.blaslib+' '+self.config.ldflags_fc+' -lm'
        (output, error, retz) = runShellCommand(fcomm)

        if(retz != 0):
            print '\n\nBLAS: provided BLAS cannot be used! aborting...'
            print 'error is:\n','*'*40,'\n',fcomm,'\n',error,'\n','*'*40
            sys.exit()

        comm = './tmpf'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print '\n\nBLAS: provided BLAS cannot be used! aborting...'
            print 'error is:\n','*'*40,'\n',comm,'\n',error,'\n','*'*40
            sys.exit()

        killfiles(['tmpf.f','tmpf'])
        print 'yes'

        return 0;


    def down_install_blas(self):

        print """
The reference BLAS library is being installed.
Don't expect high performance from this reference library!
If you want performance, you need to use an optimized BLAS library and,
to avoid unnecessary complications, if you need to compile this optimized BLAS
library, use the same compiler you're using here."""
        sys.stdout.flush()

        savecwd = os.getcwd()

        # creating the build,lib and log dirs if don't exist
        if not os.path.isdir(os.path.join(self.prefix,'lib')):
            os.mkdir(os.path.join(self.prefix,'lib'))

        if not os.path.isdir(os.path.join(os.getcwd(),'log')):
            os.mkdir(os.path.join(os.getcwd(),'log'))

        # Check if blas.tgz is already present in the working dir
        # otherwise download it
        if not os.path.isfile(os.path.join(os.getcwd(),getURLName(self.blasurl))):
            print "Downloading reference BLAS...",
            downloader(self.blasurl,self.downcmd)
            print "done"

        # unzip and untar
        print 'Unzip and untar reference BLAS...',
        comm = 'gunzip -f blas.tgz'
        (output, error, retz) = runShellCommand(comm)
        if retz:
            print '\n\nBLAS: cannot unzip blas.tgz'
            print 'stderr:\n','*'*40,'\n',comm,'\n',error,'\n','*'*40
            sys.exit()

        comm = 'mkdir BLAS && tar x --strip-components=1 -C BLAS -f blas.tar'
        (output, error, retz) = runShellCommand(comm)
        if retz:
            print '\n\nBLAS: cannot untar blas.tgz'
            print 'stderr:\n','*'*40,'\n',comm,'\n',error,'\n','*'*40
            sys.exit()
        os.remove('blas.tar')
        print 'done'

        # change to BLAS dir
        os.chdir(os.path.join(os.getcwd(),'BLAS'))

        # compile and generate library
        print 'Compile and generate reference BLAS...',
        sys.stdout.flush()
        comm = self.config.fc+' '+self.config.fcflags+" -c *.f"
        (output, error, retz) = runShellCommand(comm)
        if retz:
            print "\n\nBLAS: cannot compile blas"
            print "stderr:\n","*"*40,"\n",comm,'\n',error,"\n","*"*40
            sys.exit()

        log = output+error

        comm = "ar cr librefblas.a *.o"
        (output, error, retz) = runShellCommand(comm)
        if retz:
            print "\n\nBLAS: cannot create blas library"
            print "stderr:\n","*"*40,"\n",comm,'\n',error,"\n","*"*40
            sys.exit()
        print "done"

        log = log+output+error

        comm = self.config.ranlib+" librefblas.a"
        (output, error, retz) = runShellCommand(comm)
        if retz:
            print "\n\nBLAS: cannot create table of contents for blas library"
            print "stderr:\n","*"*40,"\n",comm,'\n',error,"\n","*"*40
            sys.exit()
        print "done"

        # write the log on a file
        log = log+output+error
        fulllog = os.path.join(savecwd,'log/blaslog')
        writefile(fulllog, log)
        print 'Installation of reference BLAS successful.'
        print '(log is in ',fulllog,')'

        # move librefblas.a to the lib directory
        shutil.copy('librefblas.a',os.path.join(self.prefix,'lib/libblas.a'))

        # set framework variables to point to the freshly installed BLAS library
        self.config.blaslib  = '-L'+os.path.join(self.prefix,'lib')+' -lblas '
        self.config.blasdir  = self.prefix
        os.chdir(savecwd)
