## Trinity process died due to 'std::bad_alloc'

This is an indicator that the process ran out of available RAM. If you have more RAM resources to make available to Trinity, then simply rerun your original Trinity command with the altered resources allocated and it should resume approximately where it left off.  

If you are resource limited, please consider running Trinity [here](Accessing-Trinity-on-Publicly-Available-Compute-Resources).  If you want to continue to try to run Trinity given your available resources, you can reduce the total RAM requirements by running Trinity with parameter '--min_kmer_cov 2'. Although the assembly should still be of high quality and require less RAM, lowly expressed transcripts may be more highly fragmented in the assembly.


<a name="ques_butterfly_GC_thread_fail"></a>
## Butterfly fails with java Error: Cannot create GC thread. Out of system resources.

There are a couple reasons why this error message might creep up.

1.  *all memory has been consumed on the machine*.  Each butterfly process wants to reserve 10G of maximum heap space.  If there's less than 10G of free memory on the machine per butterfly (--CPU setting), then java may not be able to initialize (depends on your OS configuration).  Try reducing the --CPU setting and rerunning your Trinity command. It should resume where it left off.

2.  Your server might be configured to allow only limited numbers of processes per user.  Check your hardware configuration like so:

```
    %  cat /etc/security/limits.conf 

     soft    memlock         unlimited
     hard    memlock         unlimited
     soft    stack           unlimited
     hard    stack           unlimited
        -    nofile          131072
     soft    nproc           unlimited
     hard    nproc           unlimited
```

There are various ways to deal with restricted settings (Google can help).


3.  *NUMA architecture*:  one of our users found that the java invocation required: -XX:ParallelGCThreads=<Numerical Thread Count>, otherwise it would try to use too many threads.

<a name="share_data_trin_devel"></a>
## Last resort - sharing your data privately with Trinity developers for debugging

A great way to share your data with the Trinity developers is through:
<https://mega.nz/>
as, a command-line interface can be used to easily download data from there.

Alternatively, [google drive](https://www.google.com/drive/) or [dropbox](http://www.dropbox.com) can be used.

Get the Trinity developers direct email address, rather than posting links to your data on the google forum for privacy concerns.  Trinity developers will not share your data, so no worries!

