#!/usr/bin/env python

import subprocess

def run_trinity_seqtk(inputfile, readdir):

    cmd = "../seqtk-trinity seq -A -R {} {} > /dev/null".format(readdir, inputfile)

    return subprocess.run(cmd, shell=True, check=False).returncode


## old format: @61DFRAAXX100204:1:100:10494:3070/1

def test_oldformat_1_success():
    ret = run_trinity_seqtk("oldformat_1.fq", 1)
    assert(ret==0)


def test_oldformat_1_failure():
    ret = run_trinity_seqtk("oldformat_1.fq", 2)
    assert(ret==2)


def test_oldformat_2_success():
    ret = run_trinity_seqtk("oldformat_2.fq", 2)
    assert(ret==0)


def test_oldformat_2_failure():
    ret = run_trinity_seqtk("oldformat_2.fq", 1)
    assert(ret==2)


## new format:  @M01581:927:000000000-ARTAL:1:1101:19874:2078 1:N:0:1

def test_newformat_1_success():
    ret = run_trinity_seqtk("newformat_1.fq", 1)
    assert(ret==0)


def test_newformat_1_failure():
    ret = run_trinity_seqtk("newformat_1.fq", 2)
    assert(ret==2)


    
def test_newformat_2_success():
    ret = run_trinity_seqtk("newformat_2.fq", 2)
    assert(ret==0)


def test_newformat_2_failure():
    ret = run_trinity_seqtk("newformat_2.fq", 1)
    assert(ret==2)

    
## corrupt format tests

def test_corrupt_fastq():
    ret = run_trinity_seqtk("corrupt_1.fq", 1)
    assert(ret==4)

def test_unequal_seq_n_quals_fastq():
    ret = run_trinity_seqtk("unequal_seq_n_quals_1.fq", 1)
    assert(ret==3)


