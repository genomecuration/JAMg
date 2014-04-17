rsync --exclude=sync.sh --exclude=cutting_a_release.sh --exclude="dbs/ecoli_pseudomonas*" --exclude="human_genome.fast*"  --times -rl --progress --cvs-exclude --partial --perms --copy-unsafe-links --hard-links preprocess_reads/ preprocess_reads_rel16Jan2014
tar --preserve-permissions --preserve-order --verify --exclude-vcs -cf preprocess_reads_rel16Jan2014.tar preprocess_reads_rel16Jan2014
pbzip2 preprocess_reads_rel16Jan2014.tar
rm -rf preprocess_reads_rel16Jan2014 preprocess_reads_rel17NOV2013 preprocess_reads_rel17NOV2013.tar.gz
sftp alpapan,justpreprocessmyreads@frs.sf.net

