/*
 *
 * Author: Peeyush Sahu
 * This is a script for nextflow integration into NGS pipeline.
 *
 */


/*
 * Enable DSL 2 syntax
 */

nextflow.enable.dsl = 2


/*
 * Define parameters
 */
params.genomefasta = "data/genome/genome.fa"
params.genomegff = "data/genome/genes.gtf"
params.sequencefile = "data/input/*.fa"
params.tcga = "data/TCGA"
params.bwt2ind = "data/bwt2ind/"
params.buildinx = false
params.bowtie = "tools/bowtie2-2.5.1-linux-x86_64"
params.outpath = "data/output/"

log.info """\
QC NGS  -  N F    v 2.3
================================
GenomeFasta : $params.genomefasta
GenomeGFF   : $params.genomegff
InputReads  : $params.sequencefile
TcgaData    : $params.tcga
BowtiePath  : $params.bowtie
Bowtie2Index: $params.bwt2ind
BuildIndex  : $params.buildinx
OutputDir   : $params.outpath
"""

/*
 * Import modules
 */
 include {
    bt2index_build;
    bowtie_align;
    seq_processing;
 } from './sgRNA_alignment.nf'


/*
 * Pipeline logic
 */

workflow {
    sequencefile = Channel.fromPath(params.sequencefile) // incase of multiple files
    genomefasta = Channel.fromPath(params.genomefasta)
    genomegff = Channel.fromPath(params.genomegff)
    bowtie = Channel.fromPath(params.bowtie)
    tcga = Channel.fromPath(params.tcga)
    outpath = Channel.fromPath(params.outpath)
    project_path = Channel.fromPath('./')
    count_bwtidx = new File(params.bowtie).listFiles().size()

    if (params.buildinx == true){
        bt2index_build(bowtie, genomefasta, project_path)
        bowtie_align(bowtie, bt2index_build.out, sequencefile, outpath)
        seq_processing(project_path, genomegff, bowtie_align.out, outpath, tcga)
        }
    else {
        bwt2ind = Channel.fromPath(params.bwt2ind)
        bowtie_align(bowtie, bwt2ind, sequencefile, outpath)
        seq_processing(project_path, genomegff, bowtie_align.out, outpath, tcga)
    }

}
