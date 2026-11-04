/*
 * Read and check the samplesheet, skip runs whose results are complete, quantify the rest on every reference,
 * and collect a per-batch mapping summary.
 */
include { QUANTIFY_RUN } from '../modules/local/quantify_run/main'


workflow EXPRESSION {
    main:
    if (!params.input) error "--input <samplesheet.csv> is required (config/samplesheets/<cohort>.csv)"
    def sheet = file(params.input, checkIfExists: true)
    def header = sheet.readLines()[0].split(',') as List
    def required = ['sample', 'run', 'cohort', 'batch', 'read_pairs', 'read_length', 'sra']
    def missing = required - header
    if (missing) error "samplesheet ${sheet} lacks columns ${missing}"

    def refs = params.references.tokenize(',')
    refs.each { ref ->
        if (!file("${params.ref_dir}/${ref}/genome.fa").exists()) error "reference ${ref} not found in ${params.ref_dir}"
    }

    runs = channel.fromPath(sheet)
        .splitCsv(header: true)
        .filter { row -> params.batch == 'all' || row.batch == params.batch.toString() }
        .filter { row -> !file("${params.outdir}/${row.cohort}/${row.run}/complete.txt").exists() }
        .filter { row ->
            def ok = file("${file(row.sra).parent}/${row.run}.validated").exists()
            if (!ok) log.warn "${row.run}: no validated .sra at ${row.sra}; skipped (rerun fetch-sra-batch.sbatch, then this batch)"
            ok
        }
        .map { row ->
            def sra = file(row.sra)
            def meta = [id: row.run, sample: row.sample, cohort: row.cohort,
                        read_pairs: row.read_pairs as long, read_length: row.read_length as int]
            tuple(meta, sra)
        }

    QUANTIFY_RUN(runs, refs)

    QUANTIFY_RUN.out.summary
        .collectFile(name: "mapping-summary-batch${params.batch}.tsv", keepHeader: true, skip: 1, sort: true,
                     storeDir: "${params.outdir}/summaries")
}
