Installation           
============           
               
The AUM package depends upon the following             
               
    1. `swig <http://www.swig.org>`_           
    2. `GNU scientific library <http://www.gnu.org/software/gsl>`_             
    3. C++ compiler            
               
Make sure to put the path to swig, g++, in your PATH environment variable, the         
path to the gsl include files in your INCLUDE environment variable and the part            
to the gsl library in your LD_LIBRARY_PATH variable. Then,             
               
.. sourcecode:: bash           
               
    $ python setup.py install --prefix=`pwd`/install           
    $ python setup.py install --prefix=`pwd`/install           
               
should install the python package AUM correctly. Yes you have to run the same          
command twice, so that the swig generated python modules are also copied               
correctly.             
               
If all goes well, you should have a working python library in the subdirectory         
install/lib/python2.7/...              
               
To test the installation, run:         
               
.. code-block:: python         
               
    >>> import sys             
    >>> sys.path.append('PATH_TO_COSMOLOGY.PY')                
    >>> import cosmology as cc         
    >>> # This is the default constructor with some basic cosmological parameters              
    >>> a = cc.cosmology()             
    >>> # Prints out the comoving distance in the fiducial cosmology           
    >>> print a.Dcofz(2.0)             
    >>> # Prints the abundance of 1e9 Msun halos at z=0.0              
    >>> print a.nofm(1e9, 0.0)         
    >>> # Print all functions          
    >>> help(a)                
               
Step by step instructions              
=========================              
               
Download source code           
               
    - `Source code for AUM <http://github.com/surhudm/aum>`_           
               
.. sourcecode:: bash           
               
    $ git clone http://github.com/surhudm/aum          
    $ cd aum           
    $ python setup.py install --prefix=`pwd`/install           
    $ python setup.py install --prefix=`pwd`/install           
               
Add the directory where aum was installed to your PYTHONPATH variable.         
               
.. sourcecode:: bash           
               
    $ export PYTHONPATH=$PYTHONPATH:`pwd`/install/lib/python2.7/site-packages          
               
Please modify the python path in the code above depending upon your system             
install.
