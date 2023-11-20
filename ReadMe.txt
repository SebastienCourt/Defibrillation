Sebastien Court, University of Innsbruck, Copyright 2023.

0/ This C++ code was written for performing the numerical experiments that are presented in the following preprint: 
https://arxiv.org/abs/2309.12973

Requirements: Linux system, libraries "qhull" and "MUMPS" installed.

1/ For compiling and running the code under Linux, first install the Getfem++ library:
        1.1/ Download: http://download-mirror.savannah.gnu.org/releases/getfem/stable/getfem-5.4.1.tar.gz
        1.2/ Uncompress the file: "tar xvfz getfem-5.4.1.tar.gz"
        1.3/ Create the folder "defibrillation" in "./getfem-5.4.1/contrib/" 
        1.4/ Copy the "Makefile.am" in "./getfem-5.4.1/contrib/defibrillation/"
        1.5/ Modify the file Makefile.am located in ./getfem-5.4.1/contrib/ by adding the sub-folder "defibrillation"
        1.6/ Modify the file "./getfem-5.4.1/configure.ac" by adding around line 1200 the string 
                "contrib/defibrillation/Makefile                                      \"

2/ Install Getfem as follows: In "./getfem-5.4.1/"
        2.1/ run successively "aclocal", "autoconf", "autoheader" and "automake --add-missing"
        2.2/ run "./configure --enable-qhull"
        2.3/ run "make"
        2.4/ run "sudo make install"

3/ Compile the code "SFB_maxmax2":
        3.1/ Copy files located from the folder "maxmax_elasticity" to "/getfem-5.4.1/contrib/defibrillation/"
        3.2/ First create a "MATLAB_SFB" folder in "./getfem-5.4.1/contrib/defibrillation/"
        4.2/ Run "make SFB_maxmax2"

4/ Execute the code by running "./SFB_maxmax2 SFB_maxmax2.param". The results are contained in the "MATLAB_SFB" folder.

