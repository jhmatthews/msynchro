Installation
-------------------------------

To install msynchro, run -- using python3 only! -- 

.. code:: bash

    python setup.py install

Depending on your system, you may need to run ``sudo python setup.py install`` or ``python setup.py install --user``. You may wish to install in a virtual environment. 

Requirements
====================================

Strict requirements:

	* python (version 3.5 or later)
	* numpy (version 1.17.2 or later)
	* scipy (version 1.3.1 or later)
	* setuptools (version 40.8.0 or later)

These version numbers are approximate and based on the earliest versions I have run the code with. I expect the code will work with earlier versions, but I can't guarantee it. 

You will also need a C compiler, such as gcc, clang or icc to compile the C extension. I have tested with the following compilers:

	* clang Apple LLVM version 10.0.1 (clang-1001.0.46.4)
	* gcc-10 (Homebrew GCC 10.2.0) 10.2.0


.. todo:: give more detailed requirements.

