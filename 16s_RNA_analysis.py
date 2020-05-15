"""
code
====
"""
import time
from ruffus import *
from ruffus.combinatorics import product
import sys
import os
import shutil
import sqlite3
import subprocess
import glob
from cgatcore import experiment as E
import cgat.Sra as Sra
from cgatcore import pipeline as P
import cgatcore.iotools as iotools
import tempfile

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

###################################################

@follows(mkdir("greengene.dir"))
@transform(PARAMS["database_fasta"],
	   formatter(),
	   "greengene.dir/{basename[0]}.qza")
def qzadb(infile, outfile):
    job_mem = PARAMS["qzadb_mem"]
    statement = '''qiime tools import 
                   --input-path %(infile)s
                   --output-path %(outfile)s
                   --type FeatureData[Sequence]''' 
    P.run(statement, job_condaenv="qiime2-2017.8")

##################################################

@follows(mkdir("manifest.dir"))
@merge("input_files.dir/*",
         "manifest.dir/manifest")
def createmani():
    statement = '''Rscript /data/md1nsk/Softwares/manifest_create.R'''
    P.run(statement)

#################################################

@follows(createmani, mkdir("mani_qza.dir"))
@transform("manifest.dir/manifest",
           formatter(),
           "mani_qza.dir/community_raw_sequence.qza")
def mani_qza(infile, outfile):
    job_threads = PARAMS["mani_qza_threads"]
    job_memory = PARAMS["mani_qza_mem"]
    statement = '''qiime tools import  
                    --type 'SampleData[SequencesWithQuality]' 
                    --output-path %(outfile)s
                    --input-format SingleEndFastqManifestPhred33
                    --input-path %(infile)s
		    '''
    P.run(statement, job_condaenv="qiime2-2017.8")
    #while not os.path.exists(outfile):
     #   time.sleep(1)

#################################################               
@follows(mani_qza, mkdir ("dereplicate.dir"))
@transform(mani_qza,
           formatter(),
           ["dereplicate.dir/community_raw.seqs.qza", "dereplicate.dir/community_raw.table.qza"])
def dereplicate(infile,outfiles):
    job_threads = PARAMS["dereplicate_threads"]
    job_memory = PARAMS["dereplicate_mem"]
    raw_seq, raw_table = outfiles   
    statement = ''' qiime vsearch 
                    dereplicate-sequences 
                    --i-sequences %(infile)s 
                    --o-dereplicated-table %(raw_table)s 
                    --o-dereplicated-sequences %(raw_seq)s 
                    --verbose '''
    P.run(statement, job_condaenv="qiime2-2017.8")
    
#################################################
@follows(dereplicate, mkdir("Alignment.dir"))
@product(qzadb, formatter(),
         dereplicate, formatter(),
         ["Alignment.dir/community_raw_sequence_derep_table.qza","Alignment.dir/community_raw_sequence_derep_seq.qza",
          "Alignment.dir/community_raw_sequence_derep_unmatched.qza"])
def alignment(infiles, outfiles):
    job_threads = PARAMS["alignment_threads"]
    job_memory = PARAMS["alignment_mem"]
    reference = infiles[0]
    sequence, table = infiles[1]
    derep_table, derep_seq, unmatched = outfiles
    statement = ''' qiime vsearch
                    cluster-features-closed-reference
                    --p-threads 12 
                    --i-table %(table)s
                    --i-sequences %(sequence)s
                    --i-reference-sequences %(reference)s
                    --p-perc-identity 0.85 
                    --o-clustered-table %(derep_table)s 
                    --o-clustered-sequences %(derep_seq)s
                    --o-unmatched-sequences %(unmatched)s
                    --verbose
                    '''
    P.run(statement, job_condaenv="qiime2-2017.8")

#################################################

@follows(qzadb, createmani, mani_qza, dereplicate, alignment)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))


