#! /usr/bin/env python
# -*- coding: utf-8 -*-

##
#
# @file utils.py
#
# @brief Utility functions to help installation routines.
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

import os
import re
import shutil
import subprocess
import select

##
# Writes file fname with content fill
##
def writefile(fname, fill):
    """ writes the file fname with content fill """
    fp = open(fname,'w')
    fp.write(fill)
    fp.close()

##
# Deletes a list of files
##
def killfiles(lst):
    """ deletes a list of files """
    for i in lst:
        if(os.path.isfile(i)):
            os.remove(i)

##
# Open pipe
##
def openPipe(command):
    pipe = None
    if hasattr(popen2, 'Popen3'):
        pipe   = popen2.Popen3(command, 1)
        input  = pipe.tochild
        output = pipe.fromchild
        err    = pipe.childerr
    else:
        (input, output, err) = os.popen3(command)

    return (input, output, err, pipe)

##
# This class is responsible for dealing with Null files that may be encountered.
##
class NullFile:
    def write(self, s): pass
    def flush(self): pass

##
# Runs a shell command and obtains the output and error code
##
def runShellCommand(command, outfile=None):
    """ runs a shell command """
    if not outfile:
      outfile = NullFile()

    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = proc.communicate()
    return (out, err, proc.returncode)

##
# Gets URL name
##
def getURLName(url):
    directory=os.curdir

    name="%s%s%s" % (
        directory,
        os.sep,
        url.split("/")[-1]
    )

    return name


##
# Downloads the content of a URL
##
def downloader(uri,cmd):
    """ downloads the content of an URL """

    savedir = os.getcwd()
    downdir = savedir+"/download"
    if(not os.path.isdir(downdir)):
        print"Creating directory", downdir
        os.mkdir(downdir)
    os.chdir(downdir)

    name = getURLName(uri)
    try:
        if(os.path.isfile(downdir+"/"+name)):
            print "Package "+name+" already downloaded"
        elif(cmd == 'urllib2'):
            import urllib2
            url = urllib2.urlopen(uri)
            f = open(name,'w')
            for line in url.readlines():
                f.write(line)
            url.close()
            f.close()
        elif(cmd == 'wget'):
            comm = 'wget '+uri
            (output, error, retz) = runShellCommand(comm)
        else:
            raise
    except:
        print " "
        print "================================================================================="
        print "ERROR: Cannot download "+name
        print "Make sure the network is reachable."
        print "Packages may be downloaded with wget though a proxy; in order to"
        print "enable this feature it is enough the set the http_proxy environment"
        print "variable (read the wget man pages for more details)."
        print "If you still have troubles, you can manually download "+name+" from this URL:"
        print uri
        print "into the current directory:"
        print os.getcwd()
        print ""

    os.chdir(savedir)
    shutil.copy('download/'+name, './')


##
# Fixes inputted path
##
def fixpaths(inpath):
    lst = inpath.split(" ")

    outpath = ""

    if ((len(lst) == 1) and (inpath[0] != "download") and (inpath[0] != '-')):
        outpath = os.path.abspath(inpath)
        return outpath
    else:
        for i in lst:
            if re.search("^-L",i):
                p = "-L"+os.path.abspath(i[2:])
            else:
                p = i

            outpath = outpath+p+" "

    return outpath
