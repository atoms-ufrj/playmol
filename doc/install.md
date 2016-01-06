Installation      {#install}
============

Playmol is distributed as a git repository. To download it, just run:

    git clone https://github.com/atoms-ufrj/playmol

To compile the source code and install Playmol in your system, you can do:

    cd playmol
    make
    sudo make install

To update Playmol, enter the playmol directory and execute the following commands (including
recompilation and reinstallation):

    git pull
    make
    sudo make install

----------------------------------------------------------------------------------------------------
User's Manual
----------------------------------------------------------------------------------------------------

The Playmol User's Manual is available online [here](http://atoms.peq.coppe.ufrj.br/playmol). You
can also generate a local version if you have [Doxygen](http://www.doxygen.org) (version 1.8 or
later) installed in your system. If you do not have Doxygen, you can download and install it by:

    wget http://ftp.stack.nl/pub/users/dimitri/doxygen-1.8.10.src.tar.gz
    tar -zxvf doxygen-1.8.10.src.tar.gz
    cd doxygen-1.8.10/
    sudo apt-get install cmake flex bison
    mkdir build && cd build
    cmake -G "Unix Makefiles" ../
    make
    sudo make install

In order to generate the local User's Manual, please go to the playmol directory and execute:

    make doc

The manual will be available as a file _playmol/doc/html/index.html_, which you can open using your
favorite web browser.

----------------------------------------------------------------------------------------------------
Using Playmol
----------------------------------------------------------------------------------------------------

Once Playmol is installed, you can execute a series of input scripts by typing:

    playmol file-1 [file-2 ...]

This will execute the files in sequence as if they were a unique script. To execute the scripts one
at a time, just run playmol multiple times.

Another way of runnig a playmol script is by starting it with the following line and then making it
executable (e.g. via chmod +x):

    #!/usr/local/bin/playmol

