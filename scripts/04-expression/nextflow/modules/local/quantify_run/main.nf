/*
 * QUANTIFY_RUN: one sequencing run, all references, in a single task so the FASTQ exists only once and only while
 * the task runs (on scratch). Settings are those of the pilot (scripts/03-rnaseq-pilot/align-quantify.sbatch).
 *   fasterq-dump -> read pairs must equal ENA's count (meta.read_pairs)
 *   per reference: STAR (GeneCounts + transcriptome BAM) -> salmon (alignment mode) -> idxstats
 *                  -> keep alignments in keep-regions.bed -> delete the genome BAM
 * Output folder <outdir>/<cohort>/<run>/: fastq-check.tsv, log.txt, <ref>/{Log.final.out, Log.out,
 * ReadsPerGene.out.tab, SJ.out.tab, idxstats.tsv, salmon/, regions.bam(.bai), time-star.txt, time-salmon.txt},
 * mapping-summary.tsv, and complete.txt (written last; runs that have it are skipped).
 */
process QUANTIFY_RUN {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.cohort}", mode: 'copy'

    input:
    tuple val(meta), path(sra)
    val refs

    output:
    path "${meta.id}",                         emit: results
    path "${meta.id}/mapping-summary.tsv",     emit: summary

    script:
    def run = meta.id
    def idx = "${params.star_prefix}${meta.read_length - 1}"
    """
    set -euo pipefail
    export TMPDIR=${params.tmp}/nf-\${SLURM_JOB_ID:-local}-${run}
    mkdir -p \$TMPDIR; trap 'rm -rf \$TMPDIR' EXIT
    mkdir -p ${run}
    echo "start \$(date -u +%FT%TZ) ${run} ${meta.sample} ${meta.cohort}" > ${run}/log.txt

    fasterq-dump --split-files --threads ${task.cpus} --temp \$TMPDIR -O \$TMPDIR/fq ${sra}
    n=\$(( \$(wc -l < \$TMPDIR/fq/${run}_1.fastq) / 4 ))
    printf 'run\\tsample\\tread_pairs\\texpected\\n${run}\\t${meta.sample}\\t%s\\t${meta.read_pairs}\\n' \$n > ${run}/fastq-check.tsv
    [[ \$n -eq ${meta.read_pairs} ]] || { echo "${run}: \$n read pairs, expected ${meta.read_pairs}" >&2; exit 1; }

    printf 'run\\tsample\\treference\\tinput_pairs\\tunique_pct\\tmulti_pct\\ttoo_many_pct\\tunmapped_short_pct\\tunmapped_other_pct\\n' > ${run}/mapping-summary.tsv
    for ref in ${refs.join(' ')}; do
        R=${params.ref_dir}/\$ref
        O=${run}/\$ref; mkdir -p \$O
        /usr/bin/time -v -o \$O/time-star.txt STAR --runThreadN ${task.cpus} --genomeDir \$R/${idx} \\
            --readFilesIn \$TMPDIR/fq/${run}_1.fastq \$TMPDIR/fq/${run}_2.fastq \\
            --outFileNamePrefix \$TMPDIR/\$ref. --outTmpDir \$TMPDIR/star-\$ref \\
            --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 \\
            --outFilterMismatchNoverReadLmax 0.04 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \\
            --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \\
            --outSAMattributes NH HI AS NM MD --outSAMunmapped Within \\
            --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 20000000000 \\
            --quantMode TranscriptomeSAM GeneCounts
        /usr/bin/time -v -o \$O/time-salmon.txt salmon quant -p ${task.cpus} -l A -t \$R/transcripts.fa \\
            -a \$TMPDIR/\$ref.Aligned.toTranscriptome.out.bam -o \$O/salmon
        samtools index -@ ${task.cpus} \$TMPDIR/\$ref.Aligned.sortedByCoord.out.bam
        samtools idxstats \$TMPDIR/\$ref.Aligned.sortedByCoord.out.bam > \$O/idxstats.tsv
        samtools view -@ ${task.cpus} -b -M -L \$R/keep-regions.bed -o \$O/regions.bam \$TMPDIR/\$ref.Aligned.sortedByCoord.out.bam
        samtools index \$O/regions.bam
        for f in Log.final.out Log.out ReadsPerGene.out.tab SJ.out.tab; do cp \$TMPDIR/\$ref.\$f \$O/\$f; done
        rm -f \$TMPDIR/\$ref.Aligned.*
        awk -F'|' -v r=${run} -v s=${meta.sample} -v ref=\$ref '
            { gsub(/^[ \\t%]+|[ \\t%]+\$/, "", \$1); gsub(/^[ \\t]+|[ \\t%]+\$/, "", \$2); v[\$1] = \$2 }
            END { OFS = "\\t"; print r, s, ref, v["Number of input reads"], v["Uniquely mapped reads"],
                  v["of reads mapped to multiple loci"], v["of reads mapped to too many loci"],
                  v["of reads unmapped: too short"], v["of reads unmapped: other"] }' \$O/Log.final.out >> ${run}/mapping-summary.tsv
    done
    echo "end \$(date -u +%FT%TZ)" >> ${run}/log.txt
    date -u +%FT%TZ > ${run}/complete.txt
    """
}
