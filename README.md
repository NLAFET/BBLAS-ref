# BBLAS-ref
Reference implementation of Batched BLAS routines

Please note that this is a draft specification which is currently open to suggestions from the community 
and may be subject to large changes in the future.

Compilation is achieved using the script setup.py. 
For instance to compile a basic version using only MKL we can run the following.
</br>
<code>
source /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64
</code>
</br>
<code>
./setup.py --blaslib="-Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl"
</code>

By running <code>./setup.py -h</code> you can see the options required to compile with CUDA and MAGMA.

This compilation will also create all the documentation in the "docs" folder. After compilation is complete 
(don't worry that "make install" fails) tests can be run using the script testing/run_tests.py.
By default running ./run_tests.py with no input arguments will run all available tests.
