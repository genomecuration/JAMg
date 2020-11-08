# Installing Trinity

After [downloading](https://github.com/trinityrnaseq/trinityrnaseq/releases) the software to a Linux server, simply type 
   
    %  make 

in the base installation directory.  This should build Inchworm and Chrysalis, both written in C++. Butterfly should not require any special compilation, as its written in Java and already provided as portable precompiled software, but *Java-1.8* (or higher) is required.

>Note, starting with Trinity-v2.8, [cmake](https://cmake.org/) is required for building the software.

Afterwards, you may want to build the additional plugin components that provide support for downstream analyses in which case you would then type:

    %  make plugins


Additional tools required for running Trinity include:

- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [jellyfish](http://www.genome.umd.edu/jellyfish.html)
- [salmon](http://salmon.readthedocs.io/en/latest/salmon.html)

If you want to install Trinity in a central location (not required), you can

    %  make install

and it'll copy the software package to /usr/local/bin/trinityrnaseq-version

You can set the environmental variable TRINITY_HOME to point to this, which will make it easy to access both Trinity as well as supported downstream applications that come bundled with Trinity.  

    %   export TRINITY_HOME=/path/to/trinity/installation/dir

>You can put the above command in your ~/.bashrc file so it'll be available to you by default.

Additional installation requirements:

*  python 2.7 or 3.*  with numpy


Trinity has been tested and is supported on Linux.

To test your installation of Trinity, try assembling the small sample data set provided with Trinity like so:

    cd sample_data/test_Trinity_Assembly/
    
    ./runMe.sh


