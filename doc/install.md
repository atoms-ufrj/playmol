Installation      {#install}
============

Playmol is distributed as a git repository. To download it, just run:

    git clone https://github.com/atoms-ufrj/playmol

To compile the source code and install Playmol in your system, you can do:

    cd playmol
    make
    sudo make install

Local documentation can be generated (if [Doxygen](http://www.doxygen.org) is available) by:

    make doc

To update Playmol, enter the playmol directory and execute the following command before recompilation and reinstallation:

    git pull

