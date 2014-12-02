results/$(GENOME_NAME): $(GENOME_PATH)
	ln -s $^ $@

results/%.translation: results/%.fa
	cd results && ~/src/RepeatModeler/BuildDatabase -engine ncbi -name $* $*.fa && cd -

results/RM_%: results/%.translation
	cd results && nohup ~/src/RepeatModeler/RepeatModeler -database $* -engine ncbi -pa $(LOCAL_CPU) &> $*.out &
	cd -
$(REPEAT_LIBRARY):  results/RM_*/consensi.fa.classified $(MIPS).classified 
	cat $^ > $@



results/$(GENOME_NAME).out.gff: $(REPEAT_LIBRARY)
	. src/env_vash.sh
	cd results && $(JAMG_PATH)/bin/prepare_domain_exon_annotation.pl -verbose -genome $(GENOME_PATH) -repthreads $(LOCAL_CPUS) -engine localmpi -mpi $(LOCAL_CPUS) -repeatoptions="-lib $(REPEAT_LIBRARY) -nolow" -uniprot_db $(JAMG_PATH)/databases/hhblits/refseq_plant.nr95_50kclust -trans $(JAMG_PATH)/databases/hhblits/transposons -only_repeat

