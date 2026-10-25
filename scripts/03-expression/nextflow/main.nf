#!/usr/bin/env nextflow
/*
 * t2t-eqtl step 3: quantify LCL RNA-seq runs on several references from the same reads.
 *
 *   nextflow run scripts/03-expression/nextflow/main.nf -profile slurm \
 *       --input config/samplesheets/mage.csv --batch 1 -resume
 *   nextflow run scripts/03-expression/nextflow/main.nf -profile slurm,test -resume      (2-run test)
 *
 * Normally launched per batch by scripts/03-expression/run-batch.sbatch. Inputs are the validated .sra files
 * from scripts/03-expression/fetch-sra-batch.sbatch (samplesheet column `sra`). Layout:
 *   main.nf                          entry point (this file)
 *   nextflow.config                  parameters and profiles (slurm, test)
 *   conf/base.config                 per-process resources
 *   workflows/expression.nf          samplesheet checks, skip completed runs, run QUANTIFY_RUN, summaries
 *   modules/local/quantify_run/      one task per run: fasterq-dump, then STAR + salmon per reference
 * Documentation: docs/t2t-eqtl.md section 13A.
 */
nextflow.enable.dsl = 2

include { EXPRESSION } from './workflows/expression'

workflow {
    EXPRESSION()
}
