##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################
"""===========================
Pipeline APA
===========================

:Author: Ian Sudbery
:Release: 0.01
:Date: 21/11/16
:Tags: Python


Overview
========

This pipeline aims to identify cases of alternate exon useage from
RNAseq data.  It uses two different appraoched. The DaPars program
will be applied, which bulids models of read depth over final exons to
identify cases of APA. The second is to use DEXSeq to identify cases
of alternate exon usage where the exon in question is the last exon in
a transcript.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_apa.py config

By default the pipeline will try to guess the experimental design
but a design file can be provided, called :file:`design.tsv` to
contain a different design. The file has three columns, a column
with the comparison name, and two columns with regular expressions
that match file in condition1 and condition2 respectively. 
e.g.

    #name    pattern1           pattern2
    tissue   heart-control.+    brain-control.+
    kd       heart-kd.+         heart-control.+


If a design file is not present, files with control in the second
part of the file name will be matched as controls for those with
same first part, but different second part. 

e.g. 

if heart-control-r1 and heart-kd-r1 are present, the frist will be used
as the control for the second. 

Input files
-----------

The input files are indexed bam files, named with three part
names, seperated by a dash. Traditionally part 1 is the tissue
or cell type, or experiment name, part 2 is the condition, and
part 3 is the replicate. e.g.

heart-control-R1.bam

would be the heart control from replicate one.
Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

Requirements:

* samtools >= 1.1
* DaPars
* R
* DEXSeq
* ExperimentR
* bedtools
* bgzip & tabix

Pipeline output
===============

Most of the output is in the sqlite database associated with the
pipeline (csvdb by default). Also exported are the last exon chunks
found to be differentially used by DEXSeq in the export directory.

Diagram
=======

..image:: pipeline_diagram.png

 Glossary ========

.. glossary::


Code
====

"""
from ruffus import *

import sys
import os
import sqlite3
import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
from CGATPipelines import PipelineRnaseq
from CGAT import IOTools
from CGAT import GTF

# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

# add configuration values from associated pipelines
#
# 1. pipeline_annotations: any parameters will be added with the
#    prefix "annotations_". The interface will be updated with
#    "annotations_dir" to point to the absolute path names.
PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True))


# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
    '''

  
    dbh = sqlite3.connect(PARAMS["database_name"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        os.path.join(PARAMS["annotations_dir"],
                     PARAMS["annotations_database"]))
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

# -----------------------------------------------------------------
@active_if(PARAMS["build_geneset"] == 1)
@follows(mkdir("geneset.dir"))
@transform("*.bam", formatter(),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"])),
           "geneset.dir/{basename[0]}.gtf.gz")
def assembleWithStringTie(infiles, outfile):

    infile, reference = infiles

    job_threads = PARAMS["stringtie_threads"]
    job_memory = PARAMS["stringtie_memory"]

    statement = '''stringtie %(infile)s
                           -p %(stringtie_threads)s
                           -G <(zcat %(reference)s)
                           %(stringtie_options)s
                           2> %(outfile)s.log
                   | gzip > %(outfile)s '''

    P.run()


# -----------------------------------------------------------------
@active_if(PARAMS["build_geneset"] == 1)
@merge([assembleWithStringTie,
        os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"])],
       "geneset.dir/agg-agg-agg.gtf.gz")
def mergeAllAssemblies(infiles, outfile):

    infiles = ["<(zcat %s)" % infile for infile in infiles]
    infiles, reference = infiles[:-1], infiles[-1]

    job_threads = PARAMS["stringtie_merge_threads"]

    infiles = " ".join(infiles)

    statement = '''stringtie --merge
                             -G %(reference)s
                             -p %(stringtie_merge_threads)s
                             %(stringtie_merge_options)s
                             %(infiles)s
                            2> %(outfile)s.log
                   | python %(scriptsdir)s/gtf2gtf.py --method=sort
                           --sort-order=gene+transcript
                            -S %(outfile)s -L %(outfile)s.log'''

    P.run() 


# -----------------------------------------------------------------
@transform(mergeAllAssemblies if PARAMS["build_geneset"] == 1 else
           os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"]),
           regex(".+"),
           "geneset.bed")
def getGenesetBed12(infile, outfile):
    '''Convert geneset to BED12 format'''

    statement = '''python %(scriptsdir)s/gff2bed.py 
                           --bed12-from-transcripts
                            -I %(infile)s
                            -S %(outfile)s
                            -L %(outfile)s.log '''

    P.run()


# -----------------------------------------------------------------
@transform(mergeAllAssemblies if PARAMS["build_geneset"] == 1 else
           os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"]),
           formatter(),
           "transcripts_to_genes.txt")
def generateDaParsTranscriptsToGenes(infile, outfile):

    import CGAT.GTF as GTF

    outlines = []

    for transcript in GTF.transcript_iterator(
            GTF.iterator(IOTools.openFile(infile))):
        try:
            gene_id = transcript[0].ref_gene_id
        except AttributeError:
            gene_id = transcript[0].gene_id
            
        outlines.append((transcript[0].transcript_id,
                         gene_id))

    outlines = list(set(outlines))

    IOTools.writeLines(outfile, outlines,
                       header=["#transcript_id", "gene_id"])


# -----------------------------------------------------------------
@merge([getGenesetBed12, generateDaParsTranscriptsToGenes],
       "geneset_extracted.bed")
def getDaParsGeneset(infiles, outfile):
    ''' Process geneset to generate the input file for DaPars '''

    geneset, symbols = infiles

    statement=''' DaPars_Extract_Anno.py -b %(geneset)s
                                         -s %(symbols)s
                                         -o %(outfile)s '''

    P.run()

# -----------------------------------------------------------------
@transform("*.bam", suffix(".bam"), ".bedGraph")
def bam_to_bedGraph(infile, outfile):
    '''Convert alignments into depth bedGraphs '''

    genome_file = os.path.join(PARAMS['annotations_dir'],
                               PARAMS['annotations_interface_contigs_tsv'])

    statement = ''' genomeCoverageBed -split -bg -ibam %(infile)s 
                                      -g %(genome_file)s > %(outfile)s 2> %(outfile)s.log;
                    '''

    P.run()


# -----------------------------------------------------------------
def generate_config(infiles, utrs, outfile):
    '''Generate DaPars config from template. The following parameters are 
    required in the ini:

    :param:dapars_num_least_in_group int
    :param:dapars_coverage_cufoff int
    :param:dapars_fdr_cutoff float
    :param:dapars_pdui_cutoff float
    :param:dapars_logfc_cufoff float

    '''

    template_file = os.path.join(P.snip(__file__, ".py"),
                                 "dapars_config_template.txt")
    config_template = IOTools.openFile(template_file).read()

    condition1_files, condition2_files = infiles
    condition1_files = ",".join([os.path.abspath(f) for f in condition1_files])
    condition2_files = ",".join([os.path.abspath(f) for f in condition2_files])

    dapars_outfile = os.path.join(P.snip(outfile, ".dapars_config.txt"),
                                  "dapars_out.tsv")
    outdir = os.path.dirname(os.path.abspath(dapars_outfile))
    dapars_outfile = os.path.basename(dapars_outfile)

    local_params = PARAMS.copy()
    local_params.update(locals())

    config = config_template % local_params

    with IOTools.openFile(outfile, "w") as outf:
        outf.write(config)


if not os.path.exists("design.tsv"):

    @follows(mkdir("dapars_out.dir"))
    @collate(bam_to_bedGraph,
             regex("(.+)-((?!Control).+)-(.+).bedGraph"),
             add_inputs(r"\1-Control-\3.bedGraph", getDaParsGeneset),
             r"dapars_out.dir/\1-\2.dapars_config.txt")
    def generate_dapars_config(infiles, outfile):

        condition1_files, condition2_files, utrs = zip(*infiles)
        generate_config([condition1_files, condition2_files], utrs[0], outfile)

else:

    @follows(mkdir("dapars_out.dir"))
    @subdivide("design.tsv",
               formatter(),
               add_inputs(bam_to_bedGraph, getDaParsGeneset,),
               ["dapars_out.dir/%s.dapars_config.txt" % line.split()[0]
                for line in IOTools.openFile("design.tsv")
                if not line.startswith("#")])
    def generate_dapars_config(infiles, outfile):

        bedgraphs = infiles[1:-1]
        utrs = infiles[-1]

        comparisons = [line.split() for line in IOTools.openFile("design.tsv")
                       if not line.startswith("#")]

        for name, pat1, pat2 in comparisons:
            condition1_files = [f for f in bedgraphs if re.match(f, pat1)]
            condition2_files = [f for f in bedgraphs if re.match(f, pat2)]
            generate_config([condition1_files, condition2_files], utrs, outfile)


# -----------------------------------------------------------------
@transform(generate_dapars_config,
           regex(".+/(.+).dapars_config.txt"),
           r"dapars_out.dir/\1/dapars_out.tsv_All_Prediction_Results.txt")
def run_DaPars(infile, outfile):

    job_memory = "6G"
    statement = '''DaPars_main.py %(infile)s > %(infile)s.log'''
    P.run()   
    

# -----------------------------------------------------------------
@merge(run_DaPars, "dapars.load")
def loadDapars(infiles, outfile):
    '''Munge the DaPars output to seperate transcript and gene_ids,
    and load into database'''

    infiles = " ".join(infiles)

    statement='''python %(scriptsdir)s/combine_tables.py
                   --cat=track
                   --use-file-prefix
                   --regex-filename='dapars_out.dir/(.+)/dapars_out'
                   %(infiles)s -L %(outfile)s
            |   sed 's/[|]/\\t/g'
            |   sed '1!b;s/Gene/transcript_id\\tgene_id\\tchrom\\tstrand/'
            |   %(load_statement)s
            > %(outfile)s'''


    load_statement=P.build_load_statement(
        P.toTable(outfile),
        options = "-i track -i gene_id -i transcript_id")

    P.run()


# ----------------------------------------------------------------
@active_if(PARAMS["build_geneset"] == 1)
@follows(mkdir("export"))
@transform(mergeAllAssemblies,
           formatter(),
           [r"export/agg-agg-agg.gtf.gz",
            r"export/agg-agg-agg.gtf.tbi"])
def export_geneset(infile, outfiles):

    gtffile, index = outfiles

    statement = '''zcat %(infile)s
                 | sort -k1,1 -k4,4n
                 | bgzip > %(gtffile)s;
     
                 checkpoint;
 
                 tabix %(gtffile)s -p gff'''
    P.run()

    
# ----------------------------------------------------------------    
@follows(mkdir("export"))
@transform(bam_to_bedGraph,
           formatter(),
           add_inputs(PARAMS["annotations_interface_contigs_tsv"]),
           r"export/{basename[0]}.bw")
def export_bigwigs(infiles, outfile):

    bedgraph, contigs = infiles
    statement = '''bedGraphToBigWig %(bedgraph)s %(contigs)s 
                                    %(outfile)s'''
    P.run()


# ----------------------------------------------------------------    
@follows(loadDapars,
         export_geneset,
         export_bigwigs)
def dapars():
    pass


# -----------------------------------------------------------------            
# Alternate 3' exon analysis
# -----------------------------------------------------------------            
@follows(mkdir("alt_utr_analysis.dir"))
@transform(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"]),
           formatter(),
           "geneset.dir/transcript_chunks.gtf.gz")
def get_transcript_chunks(infile, outfile):
    '''Turn the transcript models into 'chunks'. That is each
    chunks non-overalapping part of an exon or intron that may or may not be
    included in a transcript'''

    statement = '''python %(scriptsdir)s/gtf2gtf.py 
                          --method=genes-to-unique-chunks
                           -I %(infile)s
                           -L %(outfile)s.log
                 | bedtools intersect -u -a stdin  -s
                          -b <(zcat %(infile)s | awk '$3=="exon"')
                 | gzip > %(outfile)s'''

    P.run()


# -----------------------------------------------------------------
@transform(get_transcript_chunks, formatter(),
           "geneset.dir/filtered_{basename[0]}.gz")
def filter_overlapping_genes(infile, outfile):
    '''Filter out exons that overlapp with exons from another 
    gene'''


    tmp1=P.getTempFilename()
    tmp2=P.getTempFilename(shared=True)

    # the first command in the statment trancates exons that overlap
    # on oppsite strands. The second exons that overlap on the same
    # strand.  the first part of the second command identifies exons
    # that overlap on the same strand, the second part removes them
    # from the geneset.
    statement = ''' bedtools subtract -a %(infile)s -b %(infile)s -S > %(tmp1)s;

                    checkpoint;

                    bedtools merge -i <( sort -k1,1 -k4,4n %(tmp1)s)
                                   -c 6 -o count -d -2
                  | awk '$4>1'
                  | bedtools subtract -a %(tmp1)s -b stdin
                  | python %(scriptsdir)s/gtf2gtf.py
                           --method=set-transcript-to-gene
                           -L %(outfile)s.log
                  | python %(scriptsdir)s/gtf2gtf.py
                           --method=sort -L %(outfile)s.log
                  | gzip > %(tmp2)s;

                    checkpoint;

                    rm %(tmp1)s'''

    P.run()

    # renumber exons as new exons have probably been created.
    with IOTools.openFile(outfile, "w") as outf:
        for transcript in GTF.transcript_iterator(
                GTF.iterator(IOTools.openFile(tmp2))):
            nexon = 0
            for exon in transcript:
                nexon += 1
                exon = GTF.Entry().fromGTF(exon)
                exon["exon_id"] = int(nexon)
                outf.write(str(exon) + "\n")

    os.unlink(tmp2)
# -----------------------------------------------------------------
@follows(mkdir("export"))
@transform(filter_overlapping_genes, formatter(),
           "export/{basename[0]}.gz")
def export_chunks(infile, outfile):

    statement = '''zcat %(infile)s
                 | python %(scriptsdir)s/gtf2gtf.py
                          --method=set-transcript-to-gene
                           -L %(outfile)s.log
                 | sort -k1,1 -k4,4n
                 | bgzip > %(outfile)s;

                 checkpoint;
                
                 tabix -p gff %(outfile)s'''

    P.run()


# -----------------------------------------------------------------
if os.path.exists("design.tsv"):
    @follows(mkdir("alt_utr_analysis.dir"))
    @split("design.tsv", "alt_utr_anlysis.dir/*.design.tsv")
    def generate_dexseq_design_files(infile, outfiles):
        '''take the design specification for the pipeline and convert 
        into dexseq design matricies'''

        bamfiles = glob.glob("*.bam")
        bamfiles = [P.snip(os.path.basename(f), ".bam") for f in bamfiles]
        comparisons = [line.split() for line in IOTools.openFile("design.tsv")
                       if not line.startswith("#")]

        
        for name, pat1, pat2 in comparisons:
            condition1_files = [(f, "test") for f in bamfiles
                                if re.match(f, pat1)]
            condition2_files = [(f, "control") for f in bamfiles
                                if re.match(f, pat2)]
            IOTools.writeLines("alt_utr_anlysis.dir/%s.design.tsv" % name,
                               condition1_files + condition2_files,
                                header=["track","condition"])

else:

    @follows(mkdir("alt_utr_analysis.dir"))
    @collate("*.bam",
             regex("(.+)-((?!Control).+)-(.+).bam"),
             add_inputs(r"\1-Control-\3.bam"),
             r"alt_utr_analysis.dir/\1-\2.design.txt")
    def generate_dexseq_design_files(infiles, outfile):

        track = os.path.basename(infiles[0][0]).split("-")[0]
        files = [P.snip(os.path.basename(f), ".bam")
                 for p in infiles for f in p]
        files = [(f, f.split("-")[1]) for f in files]
        IOTools.writeLines(outfile, files, header=["track","condition"])


# -----------------------------------------------------------------
@transform("*.bam", formatter(),
           add_inputs(filter_overlapping_genes),
           "alt_utr_analysis.dir/{basename[0]}.tsv.gz")
def count_chunks(infiles, outfile):

    gtffile = infiles[1]
    bamfile = infiles[0]

    PipelineRnaseq.runFeatureCounts(
        gtffile,
        bamfile,
        outfile,
        job_threads=2,
        strand=0,
        options=' -f -O -T 2 --primary -p -B -C')


# -----------------------------------------------------------------
@merge(count_chunks,
       "alt_utr_analysis.dir/chunk_counts.tsv.gz")
def merge_chunk_counts(infiles, outfile):

    infiles = " ".join(infiles)
    job_memory = "10G"
    statement=''' python %(scriptsdir)s/combine_tables.py
                         -c 1,2,3,4,5,6
                         -k 7
                         --regex-filename='(.+).tsv'
                         --use-file-prefix
                         --merge-overlapping
                         %(infiles)s
                         -L %(outfile)s.log
               | gzip > %(outfile)s '''

    P.run()


# -----------------------------------------------------------------
@transform(generate_dexseq_design_files,
           regex(".+/(.+).design.txt"),
           add_inputs(merge_chunk_counts, filter_overlapping_genes),
           r"alt_utr_analysis.dir/\1.dexseq.tsv")
def run_dexseq(infiles, outfile):
    '''run dexseq on the chunks'''

    design, counts, models = infiles

    infiles = ",".join([models, counts, design])
    outfile = P.snip(outfile, ".tsv")

    job_threads = 3
    job_memory="10G"

    pipeline_src = os.path.dirname(__file__)
    script = os.path.join(pipeline_src, "run_dexseq_all.R")
    statement = ''' Rscript %(script)s 
                            --infiles %(infiles)s
                            --outfiles %(outfile)s.tsv,%(outfile)s.gtf.gz,%(outfile)s.RData
                             -p 3
                    &> %(outfile)s.log '''

    P.run()


# -----------------------------------------------------------------
@merge(run_dexseq,
       "alt_utr_analysis.dir/dexseq_results.load")
def load_dexseq(infiles, outfile):

    statement = " checkpoint;".join(
        [" sed 's/log2fold_\S+/log2fold/' %s > %s.tmp;" % (f, f)
         for f in infiles])

    P.run()

    infiles = ["%s.tmp" % f for f in infiles]
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).dexseq.tsv.tmp",
                         options="-i groupID -i featureID -i track -i padj",
                         job_memory="6G")

    for f in infiles:
        os.unlink(f)

        
@transform(load_dexseq, suffix(".load"), ".index")
def joint_index_dexseq(infile, outfile):

    db = connect()
    db.executescript('''
             DROP INDEX IF EXISTS dexseq_results_joint;
             CREATE INDEX dexseq_results_joint
                    ON dexseq_results(groupID,featureID);''')
    P.touch(outfile)


# -----------------------------------------------------------------
@follows(joint_index_dexseq,
         load_dexseq)
def dexseq():
    pass


# -----------------------------------------------------------------
@transform(os.path.join(
    PARAMS["annotations_dir"],
    PARAMS["annotations_interface_geneset_all_gtf"]),
           formatter(),
           "geneset.dir/last_exons.gtf.gz")
def get_last_exons(infile, outfile):
    '''identify exons that are the last exon in a transcript'''

    from CGAT import GTF

    gtfs = GTF.iterator(IOTools.openFile(infile))
    with IOTools.openFile(outfile, "w") as outf:

        for transcript in GTF.transcript_iterator(gtfs):
            exons = [e for e in transcript
                     if e.feature == "exon"]

            if exons[0].strand == "+":
                exons.sort(key=lambda x: x.end)
            else:
                exons.sort(key=lambda x: x.start, reverse=True)

            last_exon = exons[-1]
            outf.write(str(last_exon) + "\n")


# -----------------------------------------------------------------
@transform(get_last_exons,
           suffix(".gtf.gz"),
           add_inputs(get_transcript_chunks),
           "_chunks.gtf.gz")
def get_last_exon_chunks(infiles, outfile):
    '''Overlap the last exons with the chunks to get
    chunks that are from a last exon'''

    last_exons, chunks = infiles

    statement = '''bedtools intersect -u -a %(chunks)s -b %(last_exons)s -s
                   | gzip > %(outfile)s'''

    P.run()


# -----------------------------------------------------------------
@transform(get_last_exon_chunks, suffix(".gtf.gz"), ".load")
def load_last_exon_chunks(infile, outfile):
    '''Load gene and exon_ids for last exons into database'''

    from CGAT import GTF
    
    with P.getTempFile(shared=True) as tmpfile:
        tmpfile.write("gene_id\tchunk_id\n")
        for exon in GTF.iterator(IOTools.openFile(infile)):
            tmpfile.write("\t".join([exon.gene_id, exon["exon_id"]])+"\n")
        tmpfn = tmpfile.name

    P.load(tmpfn, outfile, options="-i gene_id -i exon_id")
    os.unlink(tmpfn)


# -----------------------------------------------------------------
@transform(load_last_exon_chunks, suffix(".load"), ".index")
def joint_index_on_last_exon_chunks(infile, outfile):
    db = connect()
    db.executescript('''
             DROP INDEX IF EXISTS last_exon_chunks_joint;
             CREATE INDEX last_exon_chunks_joint
                    ON last_exons_chunks(gene_id,chunk_id);''')
    P.touch(outfile)


# -----------------------------------------------------------------
@follows(mkdir("export"))
@transform(run_dexseq,
           formatter(),
           inputs([r"alt_utr_analysis.dir/{basename[0]}.gtf.gz",
                   get_last_exon_chunks]),
           "export/{basename[0]}.last_exons.gtf.gz")
def export_diff_last_exons(infiles, outfile):
    ''' Overlap the differential chunks with the last exon annotations.
    also convert transcript name to gene name so that browsers arn't
    confused'''

    diff, chunks = infiles

    statement = '''bedtools intersect -u -a %(chunks)s -b %(diff)s -s
                 | python %(scriptsdir)s/gtf2gtf.py
                         --method=set-gene-to-transcript
                          -S %(outfile)s -L %(outfile)s.log'''

    P.run()


# -----------------------------------------------------------------
@follows(export_diff_last_exons,
         joint_index_on_last_exon_chunks)
def last_exons():
    pass


# -----------------------------------------------------------------
@follows(last_exons,
         dexseq)
def alt_utr_analysis():
    pass


# -----------------------------------------------------------------
@follows(alt_utr_analysis,
         dapars)
def full():
    pass


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.

    Any existing report will be overwritten.
    '''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.

    This will update a report with any changes inside the report
    document or code. Note that updates to the data will not cause
    relevant sections to be updated. Use the cgatreport-clean utility
    first.
    '''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report in the CGAT downloads directory.'''

    E.info("publishing report")
    P.publish_report()

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
