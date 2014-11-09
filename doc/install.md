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


Using Playmol
-------------

Once Playmol is installed, you can execute a series of input scripts by typing:

    playmol file-1 [file-2 ...]

This will execute the files in sequence as if they were a unique script. To execute the scripts one at a time, just run playmol multiple times.

Another way of runnig a playmol script is by starting it with the following line and then making it executable (e.g. via chmod +x):

    #!/usr/local/bin/playmol


