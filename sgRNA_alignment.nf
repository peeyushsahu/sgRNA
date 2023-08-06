/*
 * This will have various process for sgRNA processing pipeline
 */

process bt2index_build {
    tag "Building bowtie2 index for given fasta"

    input:
    path bowtie
    path genomepath
    path project_path

    output:
    path bwt2ind_path

    script:
    bwt2ind_path = project_path.toString() + '/data/bwt2ind'
    index_path = bwt2ind_path + '/bwt2'
    """
    mkdir -p $bwt2ind_path
    $bowtie/bowtie2-build \
    -f \
    --threads $task.cpus \
    $genomepath \
    $index_path
    """
}

process bowtie_align {
    tag "Aligning $sequencefile with bowtie2"

    input:
    path bowtie
    path bwt2ind
    path sequencefile
    path outpath

    output:
    path samout

    script:
    filename = sequencefile.toString().split('/')[-1] + '.sam'
    samout = outpath.toString() +  '/' + filename
    """
    ${bowtie}/bowtie2 \
    -L 12 \
    -p $task.cpus \
    --un $outpath \
    -x $bwt2ind/bwt2 \
    -f $sequencefile \
    -S $samout 2>&1 | tee $outpath/$filename-bowtie2.log  #--no-unal
    """
}

process seq_processing {
    tag "Annotation mapping and expression retrieval for aligned sgRNA"

    input:
    path project_path
    path genomegff
    path samout
    path outpath
    path tcga

    script:
    filename = samout.toString().split('/')[-1]

    """
    python3 $project_path/sgRNA_processing.py $filename $genomegff $samout $outpath $tcga
    """

}
