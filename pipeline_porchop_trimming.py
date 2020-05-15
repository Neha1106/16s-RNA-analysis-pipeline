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
    statement = '''qiime tools import 
                   --input-path %(infile)s
                   --output-path %(outfile)s
                   --type FeatureData[Sequence]''' 
    P.run(statement, job_condaenv="qiime2-2017.8")

##################################################
@follows (mkdir("trimmed.dir"))
@transform ("input_files.dir/*.fastq",formatter(),
            r"trimmed.dir/{basename[0]}.fastq")
         
def trimmed(infile,outfile):
    job_mem = PARAMS["trimmed_mem"]
    job_threads = PARAMS["trimmed_threads"]
    statement = ''' porechop 
                    -i %(infile)s 
                    -o %(outfile)s '''
    P.run(statement)
 
####################################################
@follows (trimmed,mkdir("filter.dir"))
@transform ("trimmed.dir/*.fastq",formatter(),
            r"filter.dir/{basename[0]}.fastq")

def filter (infile,outfile):
    ref = PARAMS["database_fasta"]
    statement = ''' filtlong 
                   %(infile)s
                  --assembly %(ref)s
                  --min_length 1400
                  --trim 
                  --window_size 50
                  --keep_percent 90 > %(outfile)s '''
    
    P.run(statement,job_condaenv="qiime2-2020.2")
    


###################################################
     
@follows(filter,mkdir("manifest.dir"))
@merge("filter.dir/*",
         "manifest.dir/manifest")
def createmani(infile,outfile):
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
    job_threads = 3
    job_memory = "16G"
    reference = infiles[0]
    sequence, table = infiles[1]
    derep_table, derep_seq, unmatched = outfiles
    statement = ''' qiime vsearch
                    cluster-features-closed-reference
                    --p-threads 12 
                    --i-table %(table)s
                    --i-sequences %(sequence)s
                    --i-reference-sequences %(reference)s
                    --p-perc-identity 0.95 
                    --o-clustered-table %(derep_table)s 
                    --o-clustered-sequences %(derep_seq)s
                    --o-unmatched-sequences %(unmatched)s
                    --verbose
                    '''
    P.run(statement, job_condaenv="qiime2-2017.8")


#################################################
@follows (alignment,mkdir("Chimeras.dir"))
@transform (alignment, 
            formatter(),
            ["Chimeras.dir/nonchimeras.qza", "Chimeras.dir/chimeras.qza", "Chimeras.dir/stats.qza"])
def chimeras(infiles,outfiles):
    job_threads =  PARAMS["chimeras_threads"]
    job_memory = PARAMS["chimeras_mem"]
    table, seq, unmatch = infiles
    nonchimera,chimera,stats = outfiles
    statement = '''qiime vsearch uchime-denovo
                   --i-table %(table)s
                   --i-sequences %(seq)s
                   --o-chimeras %(chimera)s
                   --o-nonchimeras %(nonchimera)s
                   --o-stats %(stats)s 
                   --verbose '''
    P.run(statement, job_condaenv="qiime2-2017.8")

#################################################

@follows (chimeras,mkdir("Feature_table_filter.dir"))
@product(chimeras, formatter(),
         alignment, formatter(),
         "Feature_table_filter.dir/table-nonchimeric-wo-borderline.qza")

def feature_table(infiles,outfile):
    job_thread= PARAMS["feature_table_threads"]
    job_memory= PARAMS["feature_table_mem"]
    nonchimeras,chimeras,stats = infiles[0]
    table,seq,unmatch = infiles[1]
    statement = ''' qiime feature-table filter-features 
                    --i-table %(table)s
                    --m-metadata-file %(nonchimeras)s
                    --o-filtered-table %(outfile)s
                '''

    P.run(statement, job_condaenv="qiime2-2017.8")

##################################################

@follows (chimeras,mkdir("Feature_filter_sequence.dir"))
@product(chimeras, formatter(),
         alignment, formatter(),
         "Feature_filter_sequence.dir/table-nonchimeric-wo-borderline.qza")

def filterseq(infiles,outfile):

    job_thread=PARAMS["filterseq_threads"]
    job_memory=PARAMS["filterseq_mem"]
    nonchimeras,chimeras,stats = infiles[0]
    table,seq,unmatch = infiles[1]
    statement = ''' qiime feature-table filter-seqs
                    --i-data %(seq)s
                    --m-metadata-file %(nonchimeras)s
                    --o-filtered-data %(outfile)s
                 '''

    P.run(statement, job_condaenv="qiime2-2017.8")

##############################################

@follows (alignment,mkdir("Summary.dir"))
@transform(alignment,
           formatter(),
           "Summary.dir/table-nonchimeric-wo-borderline.qzv")

def summary(infiles,outfile):
    table,seq,unmatch = infiles
    statement = ''' qiime feature-table summarize
                     --i-table %(table)s
                     --o-visualization %(outfile)s '''
    

    P.run(statement, job_condaenv="qiime2-2017.8")


##################################################

@follows (feature_table, mkdir ("taxonomy.dir"))
@transform(PARAMS["database_taxonomy"],
           formatter(),
           "taxonomy.dir/{basename[0]}.qza")

def taxonomy (infile,outfile):
    statement = ''' qiime tools import 
                   --type 'FeatureData[Taxonomy]'
                   --input-format 'HeaderlessTSVTaxonomyFormat' 
                   --input-path %(infile)s 
                   --output-path %(outfile)s
                   '''

    P.run(statement, job_condaenv="qiime2-2017.8")


#################################################

@follows(taxonomy, mkdir("barplot.dir"))
@product(taxonomy, formatter(),
         feature_table, formatter(),
         "barplot.dir/table-nonchimeric-wo-borderline_barplot")

def barplot (infiles,outfile):
    reference = infiles[0]
    table = infiles[1]
    metadata = "metadata.dir/metadata.tsv"
    statement = '''qiime taxa barplot
                   --i-table %(table)s
                   --i-taxonomy %(reference)s
                   --m-metadata-file %(metadata)s
                   --o-visualization %(outfile)s
                  '''
    P.run(statement, job_condaenv="qiime2-2020.2")

######################################################

@follows(qzadb, trimmed,filter,createmani, mani_qza, dereplicate, alignment,
chimeras,feature_table,filterseq,summary,taxonomy, barplot)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))


