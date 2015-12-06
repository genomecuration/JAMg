#!/bin/bash
echo What is the PATH for JAMg? e.g. `pwd`
read JAMG_PATH

echo What is the maximum intron length? e.g. 70000
read MAX_INTRON_LENGTH

echo What is the path to the genome?
read GENOME_PATH

echo What is the name you want to give to the genome - computer friendly, no spaces?
read GENOME_NAME

echo How many CPUs would like to use for some of the multi-threaded scripts? e.g. 5
read LOCAL_CPUS

echo How much memory would you like to be able to use? In Gigabytes e.g. 50
read MAX_MEMORY_G

#system variables
LD_LIBRARY_PATH=$JAMG_PATH/3rd_party/lib:$JAMG_PATH/3rd_party/lib64:$JAMG_PATH/3rd_party/mysql/lib:$JAMG_PATH/3rd_party/lib:$JAMG_PATH/3rd_party/lib64:$JAMG_PATH/3rd_party/mysql/lib:$LD_LIBRARY_PATH
LDFLAGS="$LDFLAGS -L$JAMG_PATH/3rd_party/lib -L$JAMG_PATH/3rd_party/lib64 -L$JAMG_PATH/3rd_party/mysql/lib -L$JAMG_PATH/3rd_party/lib -L$JAMG_PATH/3rd_party/lib64 -L$JAMG_PATH/3rd_party/mysql/lib"
CPPFLAGS="$CPPFLAGS -I$JAMG_PATH/3rd_party/include -I$JAMG_PATH/3rd_party/mysql/include -fPIC -I$JAMG_PATH/3rd_party/include -I$JAMG_PATH/3rd_party/mysql/include -fPIC"
CFLAGS="$CFLAGS -I$JAMG_PATH/3rd_party/include -I$JAMG_PATH/3rd_party/mysql/include -fPIC -I$JAMG_PATH/3rd_party/include -I$JAMG_PATH/3rd_party/mysql/include -fPIC"
PATH=$JAMG_PATH/bin:$JAMG_PATH/3rd_party/bin:$PATH

echo " 
#!/bin/bash
export JAMG_PATH=$JAMG_PATH
export MAX_INTRON_LENGTH=$MAX_INTRON_LENGTH
export GENOME_PATH=$GENOME_PATH
export GENOME_NAME=$GENOME_NAME
export LOCAL_CPUS=$LOCAL_CPUS
export MAX_MEMORY_G=$MAX_MEMORY_G

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH

export LDFLAGS=\"$LDFLAGS\"

export CPPFLAGS=\"$CPPFLAGS\"

export CFLAGS=\"$CFLAGS\"

export PATH=$JAMG_PATH/bin:$JAMG_PATH/3rd_party/bin:$PATH

" > env.source

echo Variables exported. In the future you can run this:
echo '
 $ source env.source
'

source env.source

