# t2t-eqtl: regulatory variation that GRCh38 cannot see

**T2T-native eQTLs in lymphoblastoid cells, linked to GM12878 regulatory elements.** A computational study using public data only, started 13 October 2026. It extends the sibling project `~/research/t2t-atlas` (GM12878 ChIP-seq on T2T-CHM13 vs GRCh38) from regulatory elements to the genetic variants that act through them. Target: ASHG or ISMB 2027.

This is the project's only documentation file: aim, plan, rules, step log, step details, history and sources. Inputs that have been checked are labelled as observations; everything else is a proposal.

**Status (2026-10-27):**
- Steps 1–2b done; the user decided **go**, with MAGE (731 individuals) as discovery and Geuvadis (451) as replication.
- MAGE: expression on three references (T2T, GRCh38, FixItFelix-patched GRCh38), cis-eQTLs (chr1–22) and the comparison are done (13B.3). About 15,000 eGenes per reference. Positive control reproduced (8 of 9 patch-changed genes, all in FixItFelix regions). **Robust H0 set: 57 genes changed by T2T vs patched GRCh38, none in FixItFelix regions**, mostly SD paralog pairs where reads move between copies; 13 reversed eQTLs; about 15 well-genotyped H1 examples.
- Geuvadis: expression running (642 of 652 runs on 2026-10-27).
- Next: Geuvadis replication, ctyper copy number as a covariate, chrX, allele-specific expression at the reversal loci (15).
- How the scripts fit together: 11.4.

**Framing changed on 2026-10-17** after a literature search (2.1, 2.2):
- FixItFelix (Behera et al. 2023) already showed artefactual GRCh38 eQTLs at the chr21 false duplication, in Geuvadis.
- This study's main claim is therefore what a *complete* reference changes **beyond** those known, patched errors (H0, section 3). The chr21 result is a replication.

Contents: 1 Aim · 2 Background and gap · 3 Hypotheses · 4 Data and sources (4.2 where every file comes from, 4.3 download commands) · 5 Design · 6 Extensions · 7 Go/no-go points · 8 Risks and controls · 9 Compute and storage · 10 Outputs and venue · 11 Execution rules (11.4 how the code fits together) · 12 Step log · 13 Step 2: RNA-seq pilot · 13A Steps 3–4: expression and genotypes · 13B Steps 5–6: eQTLs and comparison · 14 History · 15 Open items · 16 Rule changes · 17 Sources

---

## 1. Aim

Find out whether a complete reference genome changes the map of genetic effects on gene expression, and tie the changes to regulatory elements in matched cells.

**In plain terms.** An eQTL is a DNA variant whose alleles go with different expression levels of a nearby gene. Every input to an eQTL study depends on the reference genome: which variants can be called, where RNA-seq reads are placed, and which genes are annotated. GRCh38 lacks some sequence, falsely duplicates some regions (e.g. part of chr21) and collapses others. So a reference change can add real eQTLs and remove artefactual ones, even for variants in sequence both references share.

**What is already known.** FixItFelix (Behera et al. 2023) patched about 20 known errors in GRCh38: false copies masked, missing paralogs added as decoys. In 449 Geuvadis LCLs this removed artefactual eQTLs, e.g. PWP2 and GATD3A, at the chr21 false duplication. That analysis is limited to a few Mb and kept GRCh38 genotypes and gene models.

**What this study adds.** The whole genome on a complete reference: genotypes, alignments and gene models all native to T2T-CHM13v2.0, compared with GRCh38 from the same raw reads.

**Central claim to test (H0):** most reference-driven expression and eQTL differences lie **outside** the regions FixItFelix patched. The pilots support this for expression: 93–96% of changed genes are outside them (2.2).

**Recommendation.**
1. Map cis-eQTLs in **MAGE** (731 LCLs, 26 populations; discovery) and **Geuvadis** (451; replication), each natively on T2T-CHM13v2.0 and on GRCh38.
2. Also on FixItFelix's patched GRCh38 (third reference, decided 2026-10-17). That separates what known patches fix from what completeness adds.
3. Classify every association as shared, T2T-only or GRCh38-only.
4. Label each difference by region (FixItFelix-patched, other SD, T2T-only sequence, rest) and by mechanism (H1–H4).
5. Test whether T2T-only eQTLs fall in GM12878 regulatory elements (from `t2t-atlas`) and show allele-specific binding in NA12878.

## 2. Background and gap

| Work | What it showed | What it did not do |
|---|---|---|
| **FixItFelix**, Behera et al. 2023, *Genome Biology* | Patched GRCh38 (1.2 Mbp false duplications masked, decoys for collapsed genes); in 449 Geuvadis LCLs, false paralogs compete for reads and some GRCh38 eQTLs (PWP2, GATD3A) are artefacts | **Closest prior work** (2.2): patched regions only (~20), GRCh38 genotypes and gene models, no new sequence, no regulatory link |
| Nurk et al. 2022, *Science*; Wagner et al. 2022, *Nat Biotechnol* | T2T adds or corrects 238 Mbp; GRCh38 false duplications (CBS, CRYAA, KCNE1) cause missed variants (2.1) | No expression |
| Aganezov et al. 2022, *Science* | 1000 Genomes (1KGP) recalled on CHM13; better mapping and variant calling; fewer false variants in medically relevant genes | No expression or regulatory analysis |
| Rhie et al. 2023, *Nature*; Lalli et al. (phased panel, Zenodo 7612953; bioRxiv 2025) | Phased CHM13v2.0 1KGP panel and recombination map; better phasing, gains mostly for rare variants in imputation | No QTLs |
| Gershman et al. 2022, *Science* | ENCODE ChIP-seq remapped to CHM13; new regulatory elements, e.g. at NBPF | No genetic variation |
| Swedish cohort on T2T, *Genome Research* 2025 | T2T removes artefactual heterozygous calls from collapsed segmental duplications (SDs) | No expression |
| T2T mapping bias at disease variants, bioRxiv Dec 2025 | ENCODE TF clusters **lifted over** to T2T near ClinVar variants; little change | No raw-read reanalysis, no QTLs |
| *Widespread Regulatory Turnover Across Human Segmental Duplications* (Morrissey, Brown, Yang, Mahony), bioRxiv 2026-01-22 (read in full 2026-10-23) | T2T-CHM13v2.0; SDs re-called with BISER (≥ 75% identity; 24.9% in non-syntenic regions); ENCODE DNase plus 6 histone marks in 23 cell types, 134 GM12878 TF ChIP-seq, RNA-seq, Hi-C (lifted hg38 → T2T). Multimapping reads allocated with Allo (bowtie `-k 25`; STAR `--outSAMmultNmax 25`). TF binding conserved across SD copies on average only 16.5% (0–39.2%); paralog expression correlated more than random; open chromatin conserved in GM12878 | **No eQTL, genotype or individual-level analysis.** Overlaps our step 7 (regulatory elements in SDs, GM12878); code repository not public (404) |
| HSat3 TF binding platforms, bioRxiv (v2) | 1,067 ChIP targets in 13 cell lines on T2T | High overlap with a satellite angle |
| TE "arms race" on T2T, bioRxiv 2026 / R Soc Open Sci | T2T ENCODE tracks, 7 modalities, 12 cell lines | A TE-only angle is largely taken |
| Cell Reports 2025, methylation with complete reference and pangenome | "T2T improves methylation" is established | — |
| NOTCH2NL Fiber-seq (Cell Genomics 2026), SD rewiring papers | Paralog-specific regulation; Fiber-seq is the standard | — |
| **ctyper**, Ma & Chaisson, *Nature Genetics* 2025 (doi 10.1038/s41588-025-02346-4; read in full 2026-10-23) | Paralog-specific copy numbers (paCN) from short reads by k-mer matching to pangenome alleles: 3,351 CNV genes plus 273 challenging medical genes, from 114 assemblies plus CHM13 and GRCh38. Expression: Geuvadis (n = 386 PEER-corrected; 450 in their matrix) and GTEx v8 (879 samples, 57 tissues), quantified on unique exons of 4,512 genes (MANE transcripts clustered at > 98% similarity); RNA alignment reference not stated, annotation GRCh38/GENCODE v39. ANOVA in Geuvadis: paCN explains 10.3% of expression variance vs 2.14% for known eQTL variants (lifted from hg19), hence the "4.81-fold", on **639 transcripts**. Divergent expression in 7.94% of paralogs | **Closest competitor for H2–H3**, but no cis-eQTL scan, no fine-mapping, no T2T-vs-GRCh38 comparison (CHM13 only as a database haplotype), small gene set. Their **1000 Genomes copy numbers are public** (3,202 samples; Zenodo 13683015), so we use them as a covariate rather than compete (section 6) |
| **MAGE**, Taylor et al., *Nature* 2024 | LCL RNA-seq for 731 1KGP individuals from all 26 populations, one lab; >15 k putatively causal eQTLs and >16 k sQTLs; raw reads SRA PRJNA851328 | Mapped on GRCh38. A larger, more diverse cohort than Geuvadis, same cell type and 1KGP genotypes |
| Schloissnig et al., *Nature* 2025 (1KGP ONT, 1,019 samples); HGSVC3 (65 samples, 130 haplotype assemblies) | >100 k sequence-resolved SVs and ~300 k VNTRs in 1KGP LCL samples | SV-eQTLs on T2T not reported in what was found; T2T-coordinate availability not yet checked |
| **FIRE / Fiber-seq in GM12878**, *A haplotype-resolved view of human gene regulation* (bioRxiv 2024, published 2025) | Long-read, single-molecule, haplotype-resolved accessibility; 204,965 peaks in GM12878; haplotype-selective elements | Long reads map uniquely where short reads fail: an orthogonal check on T2T-only short-read peaks and allele-specific effects in NA12878 |
| ENCODE4 cCRE registry, Moore et al., *Nature* 2026 | 2.35 M human cCREs; SCREEN offers a T2T set **lifted over** from GRCh38 | The lifted T2T set is the baseline that direct T2T peaks are compared against |
| PBRM1 pericentromere k-mer pipeline (*Nat Commun* 2025); excluderanges | Prior art for the `t2t-atlas` k-mer filter | — |

**Gap.** Literature searches on 2026-10-13, 2026-10-17 and 2026-10-23 (web and Europe PMC) found no eQTL study run natively on T2T, in LCLs, GTEx or elsewhere. Adjacent 2025–2026 work: the CHM13 1000 Genomes phasing panel, the Swedish mapping study, the CHM13 mapping-bias preprint, Shew et al. 2026 (MPRA of SD paralogs in GM12878) and a 1000 Chinese Pangenome pan-variant eQTL study (*Nature* 2026, doi 10.1038/s41586-026-10315-y; pangenome, not T2T). The closest work, FixItFelix (2.2), compared GRCh38 with a *patched* GRCh38 in Geuvadis. Recheck before claiming priority, and avoid the word "first".

### 2.1 Evidence that GRCh38 misplaces sequence, and what T2T fills in (searched 2026-10-17)

| Claim | Evidence | Source |
|---|---|---|
| GRCh38 is incomplete, and some of it is wrong | 151 Mbp of unknown sequence (pericentromeres, subtelomeres, recent SDs, ampliconic arrays, rDNA). The chr21 short arm "appears falsely duplicated and poorly assembled". | Nurk et al. 2022, *Science* |
| T2T fills the gaps | T2T-CHM13 adds or corrects 238 Mbp that is not co-linear with GRCh38 (76% centromeric satellites, 19% non-satellite SDs, 4% rDNA); 182 Mbp has no primary alignment to GRCh38. 219 complete rDNA copies (9.9 Mbp). | Nurk et al. 2022 |
| Genes | Draft T2T annotation: 63,494 genes, 19,969 protein-coding. Only 263 GENCODE genes are GRCh38-exclusive, and 23 of these correspond to known GRCh38 false duplications. | Nurk et al. 2022 |
| Better mapping and fewer false variants | Across 3,202 1000 Genomes samples: better mapping and calling, hundreds of thousands of variants per sample in previously unresolved regions, tens of thousands of spurious variants per sample removed, false positives in 269 medically relevant genes reduced up to 12-fold. | Aganezov et al. 2022, *Science* (abstract; full text not open access) |
| False duplications cause missed variants | Reference-specific missed variants in CBS, CRYAA and KCNE1 for short and long reads; masking the false copies raises variant recall from 8% to 100%. | Wagner et al. 2022, *Nat Biotechnol* (GIAB CMRG benchmark) |
| False duplications and collapses in GRCh38 | 1.2 Mbp falsely duplicated and 8.04 Mbp collapsed sequence, affecting 33 protein-coding genes (12 medically relevant). Their patched set: 12 falsely duplicated regions with 22 protein-coding genes (2.03 Mbp incl. the correct copies) and 9 collapsed regions with 9 genes (0.84 Mbp). In gnomAD, coverage is too low for falsely duplicated genes (CBS, KCNE1) and too high for collapsed ones (KCNJ18, GPRIN2). | Errors found by Wagner et al. 2022, *Nat Biotechnol* (GIAB CMRG; the chr21 false duplications) and Aganezov et al. 2022, *Science* (T2T vs GRCh38: further false duplications and the collapses); patched set and gnomAD coverage from Behera et al. 2023, *Genome Biol* (FixItFelix). See 2.2 |
| ... and false eQTLs | See 2.2. | Behera et al. 2023 |
| Short-read RNA-seq is affected | Immunoglobulin isotypes in scRNA-seq are misassigned on GRCh38 and resolved on T2T. In long-read whole-blood RNA-seq, GRCh38 gives more isoforms in intronic and intergenic regions, while T2T is likely more accurate in repeats. | bioRxiv 10.1101/2023.05.24.542206; PMC12469794 |
| SD paralogs differ in regulation | In an MPRA in GM12878 and SH-SY5Y, 14–24% of SD sequences differ in activity from the chimpanzee ortholog or between human paralogs (validated for CHRFAM7A, HYDIN2, SRGAP2C). SDs are usually excluded from standard pipelines because of low mapping quality. | Shew et al. 2026, *Genome Res* 36:1318 |
| Paralog copy number drives SD-gene expression | ctyper: divergent ASE in 7.94% of paralogs; better expression prediction than known eQTL SNPs. | Ma & Chaisson 2025, *Nat Genet* |
| Limits of T2T | Dark regions remain on CHM13 for short reads even with a complete reference. | bioRxiv 10.1101/2025.05.23.655776 (not yet read) |
| Quantification choice matters | The RNA-seq quantification method and the transcriptome reference substantially change eQTL detection, colocalisation and TWAS. | Head et al. 2025 (PubMed 40894657) |

### 2.2 Closest prior work: FixItFelix (Behera et al. 2023), and how this study differs

**What they did:**
- Mapped **449 Geuvadis LCL RNA-seq samples** (Lappalainen 2013, the same data as our replication cohort) with STAR to GRCh38 and to a *modified GRCh38*: the false copies masked, plus decoys for 5 collapsed genes (MAP2K3, KCNJ18, FANCD2, GPRIN2, DUSP22).
- Ran tensorQTL (±1 Mb, MAF ≥ 1%; sex, 4 genotype PCs, 30 expression PCs).
- Genotypes were GRCh38 in both arms.

**Where the errors and the patch come from** (full text read 2026-11-02, Europe PMC PMC9942314):
- FixItFelix did not discover the errors. The chr21 false duplications come from Wagner et al. 2022 (*Nat Biotechnol*, GIAB challenging medically relevant genes), "using a recently published study [13], we identified multiple genes that were accidentally duplicated". The further false duplications and all collapsed regions come from Aganezov et al. 2022 (*Science*, T2T consortium; several FixItFelix authors are on it): "from our previous work [11], we identified … 1.2 Mbp of falsely duplicated regions … and 8.04 Mbp of falsely collapsed regions".
- T2T is the template for the fix. A false copy is masked when it "is absent in T2T-CHM13". For collapsed genes, the missing paralogs were found with the T2T segmental duplication annotation (Vollger et al. 2022), BLAT and minimap2, and "T2T-CHM13 sequences of duplicate regions not represented in GRCh38 were included as decoys".
- FixItFelix's own contributions: the patched GRCh38, a tool that re-aligns only the reads in the affected regions, and the tests (variant calling, Geuvadis expression and eQTLs).

**Why patch GRCh38 instead of moving to T2T** (their stated reasons):
- Resources: "T2T-CHM13 currently lacks many of the resources and annotations that exist for GRCh38" (gnomAD, variant databases), so "the re-analysis by using T2T-CHM13 on these large consortia data may not be an efficient solution".
- Same coordinates: the patched reference keeps GRCh38 coordinates, so existing pipelines, annotations and databases still apply.
- Cost: re-aligning only the affected reads corrects a ~30× BAM in 4–5 CPU minutes, against ~24 CPU hours for a full remap.
- They also tried aligning these reads to T2T and lifting them back to GRCh38 (LevioSAM2). It was about as accurate but 50% slower.
- Their aim was a practical fix for known errors in existing GRCh38 data, not a T2T analysis. Our study asks the other question: what else changes when the whole analysis is on T2T. A whole-genome T2T run costs us a full remap (~6,200 CPU h for MAGE), T2T-native genotypes and the same RefSeq release on both assemblies.

**What they found:**
- 10,437 eGenes on GRCh38 vs 10,417 on the modified reference.
- In falsely duplicated regions, the true gene and its false paralog compete for reads, and their expression is often negatively correlated.
- Some GRCh38 eQTLs weaken strongly on the modified reference, most strikingly **PWP2 and GATD3A** (= GATD3 in RS_2025_08). For these genes the eVariant genotypes were identical between references, but the read loss depended on genotype, so the GRCh38 signal was likely artefactual. These are also top eVariants in GTEx.
- They urge caution with GRCh38 expression results in these regions.

**Consequence for us:** the pilot's chr21 result (U2AF1, GATD3, PWP2, ICOSLG, SIK1, CBS) **replicates FixItFelix**; it is not a new finding and must be cited as such.

**How this study differs:**
1. **Whole genome, not 21 patched regions.** Of the genes that change ≥ 2-fold in all five pilot samples, only 13 of 176 (Geuvadis, STAR) and 12 of 215 (MAGE, STAR) overlap FixItFelix's patched regions (Tables S1–S2, ~5.2 Mb in GRCh38 coordinates); **93–96% are outside**.
   - With the stricter classes (13A.4b), 44–64 of those outside lie in ≥ 98%-identity SDs and 10–13 in T2T-only sequence.
   - Many in `other` are STAR-only annotation-ambiguity artefacts.
   - So the robust count of new genes is those outside the patched regions where STAR and salmon agree; to be set in step 6.
2. **Native genotypes.** T2T genotypes, including variants in T2T-only sequence and without the false variants of collapsed regions (H1, H4). FixItFelix kept GRCh38 genotypes.
3. **New sequence and genes.** T2T-only genes (69–115 expressed in the pilots) and SD families not touched by the patch.
4. **Mechanism:** GM12878 regulatory elements and allele-specific binding (H5).
5. **Scale and diversity:** MAGE (731, 26 populations, 150 bp reads) as discovery.
6. **Disease:** GWAS colocalisation (H6).

FixItFelix's Table S17 (top eQTLs for genes in their regions, on both references) can serve as a direct benchmark for our chr21 results.

## 3. Question and hypotheses

How does a complete reference change the map of genetic effects on gene expression beyond the known, patchable errors of GRCh38, and can the changed associations be tied to regulatory elements in matched cells?

0. **H0, beyond patches (central claim, added 2026-10-17):** most genes whose expression or eQTLs differ between GRCh38 and T2T lie outside the GRCh38 regions FixItFelix patched. Pilot support (expression only): 93–96% of changed genes are outside them (2.2). To test on eQTLs in step 6, after the annotation and MHC filters.
1. **H1, new variants:** some T2T-only variants (in sequence absent from GRCh38) are eQTLs.
2. **H2, new or corrected genes:** genes absent from, falsely duplicated in or collapsed in GRCh38 gain measurable expression and eQTLs on T2T. *Already shown for the chr21 false duplication by FixItFelix (expression and eQTLs); new here: T2T-only genes (69–115 expressed in the pilots) and genes outside the patched regions.*
3. **H3, read placement:** for genes and variants present on both references, more accurate read assignment on T2T changes eQTL detection, mainly in SDs.
4. **H4, artefact removal:** a subset of GRCh38 eQTLs are artefacts that disappear on T2T. Two routes:
   - **(a) expression:** genotype-dependent read competition between a gene and its false or collapsed paralog. Shown by FixItFelix for PWP2 and GATD3A.
   - **(b) genotype:** paralogous sequence variants called as SNPs in collapsed GRCh38 regions. Not tested by FixItFelix, which kept GRCh38 genotypes; testable here with T2T-native genotypes.
5. **H5, mechanism:** T2T-only eQTLs are enriched in GM12878 T2T-only peaks and show allele-specific binding in NA12878 more often than matched background variants.
6. **H6, disease interpretation:** at some GWAS loci, the target gene implied by eQTL colocalisation changes between references. This tests GWAS without lifting GWAS data into new sequence: GWAS variants lie in shared sequence, but the expression they colocalise with is measured differently on T2T (H2–H4). Typical cases are loci next to SD paralog families or GRCh38 false duplications.

## 4. Data and sources

### 4.1 Overview

| Input | Source | Status |
|---|---|---|
| Genotypes | Phased 1KGP T2T-CHM13v2.0 panel, SHAPEIT5 v1.1, 3,202 samples (`Phased_SHAPEIT5_v1.1`, human-pangenomics S3) | **Downloaded 2026-10-13, 13.0 GB, sizes match S3.** 451/464 Geuvadis samples present; NA12878 present (NA12878 is not a Geuvadis sample) |
| Expression | Geuvadis LCL RNA-seq (Lappalainen et al. 2013), ENA PRJEB3366 | **Checked:** 667 paired-end runs (2 × 76 bp, unstranded), 464 samples (291 with one run, 168 with two, 5 with eight), ~2.5 TB FASTQ, mean 26.5 M pairs and 3.8 GB per run. Run report with md5s: `data/metadata/geuvadis/ena-PRJEB3366-fastq-md5-2026-10-17.tsv`. Sample sheet with sex, ancestry and lab: `data/metadata/geuvadis/E-GEUV-1.sdrf.txt` (BioStudies E-GEUV-1). Five pilot samples only downloaded |
| Gene models | NCBI RefSeq annotation release **RS_2025_08** on both references (decision 2026-10-17, section 13) | On disk: T2T 58,563 genes, GRCh38 59,257; 55,334 NCBI GeneIDs on both |
| GRCh38 genotypes | 1KGP high-coverage GRCh38 calls (NYGC) | To obtain. Must match the CHM13 calls' pipeline as closely as possible (section 8) |
| Regulatory elements | `t2t-atlas` GM12878 ChIP-seq peaks/domains on both references, with shared/T2T-only/new-sequence classes | Available (read-only link `data/external/t2t-atlas-analysis`) |
| Allele-specific binding | `t2t-atlas` GM12878 ChIP-seq alignments plus NA12878 phased T2T genotypes | Genotypes extracted; ASB not started |
| References | CHM13v2.0 (NCBI names) and GRCh38 no-alt analysis set, shared with `t2t-atlas` | Linked. Differences that matter for RNA-seq handled in step 2 (chrM, chrEBV, chrY PARs) |
| Replication expression cohort | MAGE LCL RNA-seq, 731 samples, SRA PRJNA851328 | Not checked. Check sample overlap with Geuvadis and the CHM13 panel; one known sample swap (SRR19762530 is NA11919, not HG00237, per the MAGE README 2026-03-27) |
| GWAS summary statistics | Immune-mediated traits relevant to LCLs (e.g. SLE, RA, MS, IBD, T1D), GRCh38, GWAS Catalog / Open Targets | Not obtained. Used only for colocalisation (H6) |
| Long-read GM12878 / NA12878 data | Fiber-seq FIRE peaks and reads (GM12878); GIAB HG001 ONT simplex reads and near-T2T assemblies (hifiasm-ONT, *Nature* 2026) | Availability and formats not yet checked |
| Paralog copy number | ctyper genotypes for 1KGP samples, or run ctyper on 1KGP short reads | Not checked whether precomputed genotypes are public |

### 4.2 Where every resource comes from

Every input the project uses, with its origin, where it is on disk, how it was obtained and how it was checked. Paths are relative to the project root (`data/…` is on `/storage/projects/yusuf.maulana/t2t-eqtl/`). "t2t-atlas" means the file was prepared in `~/research/t2t-atlas` and is read here through a link; its preparation is documented in `t2t-atlas/docs/execution.md` section 5.

**Samples and expression**

| Resource | Origin | Accession / URL | Retrieved | Local path | Obtained by | Checked by |
|---|---|---|---|---|---|---|
| Geuvadis RNA-seq reads (LCLs, 2 × 76 bp, unstranded) | Geuvadis consortium; Lappalainen et al. 2013, *Nature* | ENA study PRJEB3366 (667 runs, 464 individuals); same run IDs mirrored at NCBI SRA | Pilot: 5 runs, 2026-10-17, **from NCBI SRA** | `work/03-rnaseq-pilot/fastq/<run>/` (scratch; deleted after use in step 3) | `scripts/03-rnaseq-pilot/fetch-sra.sbatch` (`prefetch`, `fasterq-dump`) | `vdb-validate` (NCBI checksums) + read-pair count equal to ENA's |
| ENA run report (run, individual, read and byte counts) | ENA portal API | `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB3366&result=read_run` | 2026-10-13 | `data/metadata/geuvadis/ena-PRJEB3366-2026-10-13.tsv` | ENA API query (13 Oct session) | — |
| ENA run report with FASTQ URLs and md5s | ENA portal API (fields `fastq_ftp,fastq_md5,fastq_bytes,read_count,nominal_length`) | as above | 2026-10-17 | `data/metadata/geuvadis/ena-PRJEB3366-fastq-md5-2026-10-17.tsv` | `curl` | — |
| Geuvadis sample sheet (sex, ancestry, sequencing lab per run) | Geuvadis via EBI BioStudies / ArrayExpress | E-GEUV-1, `https://www.ebi.ac.uk/biostudies/files/E-GEUV-1/E-GEUV-1.sdrf.txt` | 2026-10-17 | `data/metadata/geuvadis/E-GEUV-1.sdrf.txt` | `curl` | — |
| MAGE run report (run, individual, population, sex, platform, read and byte counts) | ENA portal API | `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA851328&result=read_run` | 2026-10-17 | `data/metadata/mage/ena-PRJNA851328-2026-10-17.tsv` | `curl` | 779 runs, 731 individuals, 26 populations, NovaSeq 6000, 2 × 150 bp, median 28.4 M pairs, 2.51 TB; all 731 in the T2T panel; 93 also in Geuvadis |
| Expression individuals and whether they have T2T genotypes (MAGE 731 of 731; Geuvadis 451 of 464; 1,089 included) | Derived: MAGE and Geuvadis ENA run reports (sample_alias) and the panel's sample names | — | 2026-10-13 (Geuvadis lists, by hand); 2026-10-27 (one table, both cohorts) | `config/samples/expression-individuals.tsv` (sample, mage, geuvadis, t2t_panel, included) | `scripts/01-genotypes/make-sample-lists.sh` | identical to the earlier lists (451 / 13 / 1,089), checked 2026-10-27 |
| Pilot sample sheet (5 runs) | Derived: ENA report + E-GEUV-1; selection rule in 13.1 | — | 2026-10-17 | `config/samples/pilot-5.tsv` | `scripts/03-rnaseq-pilot/make-pilot-sheet.py` | reproduced exactly |

**Genotypes**

| Resource | Origin | Accession / URL | Retrieved | Local path | Obtained by | Checked by |
|---|---|---|---|---|---|---|
| Phased 1000 Genomes genotypes on T2T-CHM13v2.0 (3,202 individuals) | Lalli et al., SHAPEIT5 v1.1 (released 2024-06-19), phasing the CHM13 calls of Rhie et al. 2023 / Aganezov et al. 2022 | Zenodo doi:10.5281/zenodo.7612953; `s3://human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/` (file `1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf.gz`) | 2026-10-13 | `data/external/1kgp-chm13v2.0-phased-shapeit5-v1.1/` (with the authors' `README.txt`) | `scripts/01-genotypes/chm13v2.0/fetch-panel.sh` | file sizes equal S3 (no md5 published) |
| Expression individuals (MAGE + Geuvadis, 1,089) + NA12878 subset of the panel | Derived | — | 2026-10-13 (451 Geuvadis, job 356292); 2026-10-17 (1,089, job 360598) | `results/01-genotypes/chm13v2.0/` (`expression/`, `na12878/`, `summary.tsv`, `counts-by-chrom.tsv`) | `scripts/01-genotypes/chm13v2.0/subset-genotypes.sbatch` | — |

**Reference genomes and annotation**

| Resource | Origin | Accession / URL | Retrieved | Local path | Obtained by | Checked by |
|---|---|---|---|---|---|---|
| T2T-CHM13v2.0 genome (24 sequences, NCBI names, chrY from HG002, no chrM) | NCBI RefSeq | GCF_009914755.1, `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0` | cluster admin bundle `/storage/datasets/ref_genomes/t2t_chm13v2` (built by the admins' `build_t2t_full_reference_bundle.sh`); copied 2026-09-27 | `data/references/chm13v2.0-ncbi/genome/` (t2t-atlas) | t2t-atlas `scripts/01-references/prepare-references.sbatch` | t2t-atlas `CHECKSUMS.md5` |
| T2T RefSeq gene annotation, release RS_2025_08 (58,563 genes) | NCBI RefSeq | GCF_009914755.1-RS_2025_08 | same admin bundle | `data/references/chm13v2.0-ncbi/annotations/ncbi/T2T-CHM13v2.0.ncbiRefSeq.gtf` (t2t-atlas link to the bundle) | — | not checksummed (admin-managed) |
| GRCh38 no-alt analysis set (195 sequences: primary + unplaced, chrM = rCRS, chrEBV, chrY PARs masked; UCSC names) | NCBI; the set ENCODE uses | GCA_000001405.15, `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids` | 2026-10-04/05 | `data/references/grch38-noalt/genome/` (t2t-atlas) | t2t-atlas `scripts/01-references/prepare-noalt-reference.sbatch` | NCBI `md5checksums.txt`; t2t-atlas `CHECKSUMS.md5` |
| GRCh38 RefSeq gene annotation, release RS_2025_08 (59,257 genes), UCSC names | NCBI RefSeq (GCF_000001405.40-RS_2025_08, GRCh38.p14), from the admin bundle `/storage/datasets/ref_genomes/grch38`, renamed with UCSC `hg38.chromAlias.txt` | — | 2026-10-05 | `data/references/grch38-noalt/annotations/ncbi/GRCh38_noalt.ncbiRefSeq.gtf` (t2t-atlas) | same t2t-atlas script | — |
| chrM and chrEBV added to the T2T STAR reference | Copied from the GRCh38 no-alt FASTA above (chrM = rCRS); chrM genes from the GRCh38 GTF | — | 2026-10-17 | inside `data/star-references/chm13v2.0/` | `prepare-star-references.sbatch` | identical sequence on both references by construction |
| STAR references (genome, GTF, gene and transcript tables, transcript FASTA, STAR 2.7.11b index) | Derived from the four rows above (13.2) | — | 2026-10-17 | `data/star-references/<chm13v2.0\|grch38>/` | `scripts/02-references/prepare-star-references.sbatch` (job 360138) | `CHECKSUMS.md5` in each folder; chrY PAR1 all N checked |
| T2T sequence absent from GRCh38 (`chm13v2-unique_to_hg38.bed`, 616 regions, 251.3 Mb; the source of the panel's `SYNTENIC` flag) | T2T consortium | `https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed` | 2026-10-17 | `data/external/t2t-consortium/` (read-only) | `curl` with retry (6 tries) | no published checksum; local md5 63b48dfa9cee799c878c7e11e63b0af3 |
| T2T segmental duplications, satellites (CenSat v2.1), RepeatMasker 4.1.2p1 | T2T consortium | `https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/` | 2026-10-06/07 | `data/references/t2t-annotations/` (t2t-atlas, with README) | t2t-atlas | t2t-atlas `CHECKSUMS.md5` |
| Liftover chains (hg38 ↔ hs1) and T2T name table (`hs1.chromAlias.txt`) | UCSC | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/`, `…/hs1/liftOver/`, `…/hs1/bigZips/` | 2026-10-04 | `data/references/liftover/` (t2t-atlas, with README) | t2t-atlas | t2t-atlas `CHECKSUMS.md5` |
| T2T chromosome name map (`chr1` ↔ `NC_060925.1`, plus chrM, chrEBV) | Derived from `hs1.chromAlias.txt` | — | 2026-10-17 | `config/chrom-names-chm13v2.0.tsv` | `awk` on the alias file (13.4) | all 24 FASTA names mapped |

Not used: t2t-atlas also has `chm13v2.0-maskedY` (the T2T consortium analysis set with chrY PARs masked and rCRS chrM, no chrEBV). The STAR reference here is equivalent apart from adding chrEBV. The 1000 Genomes population table (`20130606_g1k_3202_samples_ped_population.txt`) could not be downloaded (connection resets, 2026-10-17); ancestry comes from E-GEUV-1 instead.

**Prior-work tables**

| Resource | Origin | URL | Retrieved | Local path | Checked by |
|---|---|---|---|---|---|
| ctyper genotypes for the 1000 Genomes Project (`PangenomeAlleles_1kg.tar.gz`: 3,202 samples × 84,346 pangenome-allele types; preprint-era v1.0 database) and their Geuvadis expression matrix (`geuvadis_matrix.tar.gz`: 450 samples × 4,449 transcript units) | Ma & Chaisson 2025, *Nat Genet* | Zenodo 10.5281/zenodo.13683015 | 2026-10-23 | `data/external/ctyper-ma2025/` (read-only) | md5 equal to Zenodo (`MD5SUMS`) |


| Resource | Origin | URL | Retrieved | Local path | Checked by |
|---|---|---|---|---|---|
| FixItFelix Additional file 2 (S1 falsely duplicated regions, S2 collapsed regions, S16 Geuvadis samples used, S17 eQTL results in affected genes, …) | Behera et al. 2023, *Genome Biol* | `https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-023-02863-7/MediaObjects/13059_2023_2863_MOESM2_ESM.xlsx` | 2026-10-17 | `data/external/fixitfelix-behera2023/13059_2023_2863_MOESM2_ESM.xlsx` (read-only) | — |
| FixItFelix patched GRCh38 (`GRCh38_masked_v2_decoy_gene.fasta.gz`, 873 MB; built on the GRCh38 no-alt analysis set) | Behera et al. 2023 | Zenodo 7532965 | 2026-10-17 | `data/external/fixitfelix-behera2023/GRCh38_masked_v2_decoy_gene.fasta.gz`; prepared as `data/star-references/grch38-ffx/` | md5 35d00c3c48d1849adb8e02086c3a1c70 (Zenodo), checked in `prepare-ffx-reference.sbatch` |

**Regulatory elements**

| Resource | Origin | Local path | Notes |
|---|---|---|---|
| GM12878 ChIP-seq peaks and domains on both references, with shared / T2T-only / new-sequence classes | ENCODE GM12878 ChIP-seq (20 experiments + 3 controls) reprocessed from raw reads in t2t-atlas | `data/external/t2t-atlas-analysis` → `t2t-atlas/results/analysis/` (folders `01-…` to `07-…`) | Read-only; ENCODE accessions are listed in t2t-atlas's manifests |

**Being obtained (2026-10-17):**

| Resource | Origin | URL | Local path | Obtained by | Checked by |
|---|---|---|---|---|---|
| GRCh38 phased 1000 Genomes panel (3,202 samples; SNV, indel, SV; chr1–22, chrX v2) | NYGC; Byrska-Bishop et al. 2022, *Cell* | `http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/` | `data/external/1kgp-grch38-nygc-phased-20220422/` | `scripts/01-genotypes/grch38/fetch-panel.sbatch` | md5 from the server's `20220804_manifest.txt` |
| MAGE RNA-seq reads (779 runs, 2 × 150 bp, unstranded) | Taylor et al. 2024, *Nature* | NCBI SRA PRJNA851328 | `work/04-expression/sra/mage/` (transient) | `scripts/04-expression/fetch-sra-batch.sbatch` | `vdb-validate` + read pairs equal to ENA |

**Planned, not yet obtained:** GWAS summary statistics, long-read GM12878 / HG001 data, ctyper copy numbers (see 4.1).

### 4.3 How each file was downloaded (commands)

Exact commands, grouped by reference. All downloads run on Karakoram. HTTPS to AWS S3, `ftp.1000genomes.ebi.ac.uk` and (at times) `ftp.sra.ebi.ac.uk` is reset during the TLS handshake (2026-10-13/17), so every download retries and resumes (`curl -C - --retry … --retry-all-errors`, or `prefetch`, which resumes). `$ROOT` = `/storage/home/yusuf.maulana/research/t2t-eqtl`.

**T2T-CHM13v2.0 side**

| File | How | Command (essentials) |
|---|---|---|
| Phased 1KGP genotypes on T2T (13.0 GB BCF + index + README) | `scripts/01-genotypes/chm13v2.0/fetch-panel.sh` | `curl -sSf -C - --connect-timeout 20 --retry 50 --retry-delay 1 --retry-all-errors -o <file> https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/<file>`, looped up to 2,000 times; size checked against S3 |
| T2T genome FASTA and RefSeq GTF (GCF_009914755.1, RS_2025_08) | Not downloaded by this project: the cluster admins' bundle `/storage/datasets/ref_genomes/t2t_chm13v2` (their `build_t2t_full_reference_bundle.sh` uses `wget -c` from `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0`); registered in t2t-atlas by `scripts/01-references/prepare-references.sbatch` | — |
| T2T SD, CenSat, RepeatMasker tracks | t2t-atlas (2026-10-06/07) | `curl` from `https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/<file>` |
| `chm13v2-unique_to_hg38.bed` | by hand, 2026-10-17 | `curl -sSf -C - --connect-timeout 20 --retry 20 --retry-delay 1 --retry-all-errors -o data/external/t2t-consortium/chm13v2-unique_to_hg38.bed https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-unique_to_hg38.bed`, looped |
| `hs1.chromAlias.txt`, `hs1ToHg38.over.chain.gz` | t2t-atlas (2026-10-04) | `curl` from `https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/` and `…/hs1/liftOver/` |

**GRCh38 side**

| File | How | Command (essentials) |
|---|---|---|
| GRCh38 no-alt analysis set FASTA (GCA_000001405.15) + NCBI BWA index | t2t-atlas `scripts/01-references/prepare-noalt-reference.sbatch` (2026-10-04/05) | `curl -sSfL --retry 5 -C - -o <file> https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/<file>`; md5 against the directory's `md5checksums.txt` |
| GRCh38 RefSeq GTF (GCF_000001405.40, RS_2025_08) | admin bundle `/storage/datasets/ref_genomes/grch38` (`wget -c` from `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14`), renamed to UCSC names in t2t-atlas with the bundle's `hg38.chromAlias.txt` | — |
| Phased 1KGP genotypes on GRCh38 (NYGC; chr1–22, chrX v2, `.tbi`, manifest, README, change log) | `scripts/01-genotypes/grch38/fetch-panel.sbatch` (2026-10-17) | `curl -sSf -C - --retry 5 --retry-all-errors --retry-delay 20 -o <file> http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/<file>` (plain HTTP; HTTPS and the AWS mirror reset), 4 in parallel; each file's md5 checked against `20220804_manifest.txt` |
| `hg38ToHs1.over.chain.gz` | t2t-atlas (2026-10-04) | `curl` from `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/` |

**Patched GRCh38 (FixItFelix, third reference)**

| File | How | Command (essentials) |
|---|---|---|
| `GRCh38_masked_v2_decoy_gene.fasta.gz` (873 MB) | `scripts/02-references/prepare-ffx-reference.sbatch` (2026-10-17) | `curl -sSfL -C - --retry 5 --retry-all-errors -o data/external/fixitfelix-behera2023/GRCh38_masked_v2_decoy_gene.fasta.gz "https://zenodo.org/records/7532965/files/GRCh38_masked_v2_decoy_gene.fasta.gz?download=1"`; md5 must be `35d00c3c48d1849adb8e02086c3a1c70` (Zenodo) |
| ctyper Zenodo files | by hand, 2026-10-23 (a research subagent downloaded them; checked in the project) | `curl -sSL -o <file> https://zenodo.org/records/13683015/files/<file>`; md5 compared with the Zenodo API (`https://zenodo.org/api/records/13683015`) |
| FixItFelix Additional file 2 (Tables S1–S19) | by hand, 2026-10-17 | `curl -sSfL -o data/external/fixitfelix-behera2023/13059_2023_2863_MOESM2_ESM.xlsx https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-023-02863-7/MediaObjects/13059_2023_2863_MOESM2_ESM.xlsx` |

**Both references (shared, used identically)**

| File | How | Command (essentials) |
|---|---|---|
| Geuvadis run report | 2026-10-13 session (command not saved; same API) | ENA portal API, `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB3366&result=read_run&format=tsv` |
| Geuvadis run report with URLs and md5s | by hand, 2026-10-17 | `curl -sSf -o data/metadata/geuvadis/ena-PRJEB3366-fastq-md5-2026-10-17.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB3366&result=read_run&fields=run_accession,sample_alias,fastq_ftp,fastq_md5,fastq_bytes,read_count,nominal_length&format=tsv"` |
| Geuvadis sample sheet E-GEUV-1 | by hand, 2026-10-17 | `curl -sSfL -o data/metadata/geuvadis/E-GEUV-1.sdrf.txt https://www.ebi.ac.uk/biostudies/files/E-GEUV-1/E-GEUV-1.sdrf.txt` |
| MAGE run report | by hand, 2026-10-17 | `curl -sSf -o data/metadata/mage/ena-PRJNA851328-2026-10-17.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA851328&result=read_run&fields=run_accession,sample_alias,sample_title,library_layout,library_strategy,library_source,library_selection,instrument_model,read_count,base_count,fastq_bytes,first_public&format=tsv"` |
| RNA-seq reads, pilots (Geuvadis, MAGE) | `scripts/03-rnaseq-pilot/fetch-sra.sbatch` | `prefetch --max-size 50G -O <dir> <run>` (retried up to 100 times) → `vdb-validate` → `fasterq-dump --split-files` → `pigz`; read pairs must equal ENA's `read_count` |
| RNA-seq reads, full cohorts (step 3) | `scripts/04-expression/fetch-sra-batch.sbatch <cohort> <batch>` | `prefetch --max-size 100G -O work/04-expression/sra/<cohort> <run>` (retried) → `vdb-validate` → `<run>.validated`; FASTQ made inside the pipeline and checked against ENA's `read_count` |
| RNA-seq reads via ENA (not used: too slow, ~0.75 MB/s, 2026-10-17) | `scripts/03-rnaseq-pilot/fetch-fastq.sbatch` | `curl -sSf -C - --retry 3 --retry-all-errors -o <file> https://ftp.sra.ebi.ac.uk/vol1/fastq/<…>`; md5 against ENA |

**Software not available through pixi's normal route**

| Package | How | Command (essentials) |
|---|---|---|
| tensorQTL 1.0.10, qtl 0.1.10 (PyPI wheels) | by hand, 2026-10-17 | `curl -sSfL -o envs/vendor/<wheel> https://files.pythonhosted.org/packages/…/<wheel>`; sha256 checked against the PyPI JSON (`93e0da92…` tensorqtl, `04b387c8…` qtl); unzipped to `envs/vendor/<name>-<version>/`; `PYTHONPATH` in `pixi.toml` `[activation.env]` |

Not obtained: the 1000 Genomes population table (`https://ftp.1000genomes.ebi.ac.uk/…/20130606_g1k_3202_samples_ped_population.txt`, reset on every attempt); ancestry comes from E-GEUV-1 and the MAGE run report instead.

### 4.4 Step 1 genotype counts

**Step 1 observations** (`results/01-genotypes/chm13v2.0/summary.tsv`, the 1,089 expression individuals (MAGE + Geuvadis), multiallelic sites split, MAF recomputed in this set):

| Variant class | MAF ≥ 0.05 | 0.01–0.05 | < 0.01 |
|---|---:|---:|---:|
| Sequence syntenic with GRCh38 | 7.76 M | 6.15 M | 49.3 M |
| T2T-only sequence (no `SYNTENIC` flag) | 18.2 k (9.5 k SNPs, 8.7 k indels) | 25.5 k | 157 k |

NA12878: 2.33 M heterozygous sites in syntenic sequence, 5.9 k in T2T-only sequence. T2T-only common variants cluster on chr9 (4.0 k) and chr15 (3.3 k).

**Interpretation.** Common variants (MAF ≥ 0.05) in new sequence are about 0.23% of the total, because most new sequence is satellite and short-read calls there are sparse. H1 can contribute examples but not the main result. The study's weight rests on H2–H4, which act on all 7.8 M shared variants as well. The share of rare variants among T2T-only calls is high, which suggests some calling error in difficult sequence; restrict association testing to common variants.

(The first run, 2026-10-13, used the 451 Geuvadis individuals only and gave the same picture: 17.5 k T2T-only vs 7.39 M syntenic common variants, 0.24%. Its files were retired on 2026-10-27, section 16.)

## 5. Design

The core is a **controlled comparison**: the same samples and raw reads, the same pipeline and versions, and matched gene models, so that the only intended difference is the reference.

| Step | What | Output |
|---|---|---|
| 1 | Genotypes: subset the CHM13 panel to Geuvadis and NA12878 | **done** (`results/01-…`) |
| 2 | Pilot: 5 Geuvadis samples aligned (STAR) to both references, quantified | Runtime/storage per sample; expression differences at candidate genes → **go/no-go 1** (section 13) |
| 3 | Full expression: MAGE (731 individuals, discovery) and Geuvadis (451, replication), both references; per run download → quantify → delete (section 13A) | Gene (and transcript) matrices on each reference |
| 3b | *Decided 2026-10-17:* the same runs on FixItFelix's patched GRCh38 (Zenodo 7532965, built on the same GRCh38 no-alt set) as a third reference, with GRCh38 genotypes | Three-way comparison: GRCh38 → patched GRCh38 → T2T. Separates what known patches fix from what completeness adds |
| 4 | GRCh38 genotypes for the same samples, harmonised filters | Matched genotype sets |
| 5 | cis-eQTL mapping (tensorQTL), ±1 Mb, covariates: genotype PCs, PEER/PCA expression factors, sex (choose settings comparable to FixItFelix: 4 genotype PCs, expression PCs at saturation, in-sample MAF ≥ 1%, so the chr21 results can be checked against theirs) | Per-gene permutation p-values, nominal statistics, fine-mapped credible sets (SuSiE) |
| 6 | Cross-reference comparison | Shared / T2T-only / GRCh38-only eGenes and variant–gene pairs, each assigned a mechanism (H1–H4) and a **region class**: FixItFelix-patched (Behera 2023 Tables S1–S2), other SD, T2T-only sequence, rest; annotation-difference and MHC-haplotype genes flagged (13.7). Tests H0. FixItFelix Table S17 serves as a benchmark for the patched regions |
| 7 | Regulatory link | Enrichment of T2T-only eQTL credible sets in GM12878 peak classes vs matched background; case loci |
| 8 | Allele-specific binding in NA12878 (WASP-filtered) at heterozygous sites in peaks, on both references | Independent support for H5 |
| 9 | Replication in Geuvadis (individuals not in MAGE; roles swapped 2026-10-17, section 13.8), same pipeline | Replication rate of T2T-only and GRCh38-only eQTLs |
| 10 | GWAS colocalisation (coloc / SuSiE-coloc) of each reference's eQTLs with immune GWAS loci | Loci where the implicated gene changes between references (H6) |

**Matching between references.** Genes are matched by NCBI GeneID (same RefSeq release on both). Associations are paired by converting coordinates (chain files) only for results computed independently from raw data; no input data are lifted over. Pairs that do not convert one-to-one are reported as their own category.

**Calling a difference.** Significant on one reference and not the other is not enough (threshold effects). Compare effect sizes and test both references jointly (π1 replication, mashr); call "T2T-only" only when the association is untestable or clearly absent on GRCh38, and record why (variant absent, gene absent, expression unmeasurable, effect differs).

## 6. Extensions, by priority

Searched 2026-10-13. **Core** items strengthen the main claim and should be in the first submission. **Optional** items are added if time allows, or kept for the journal version.

| Priority | Extension | What it adds | Cost |
|---|---|---|---|
| Core | **GWAS colocalisation (H6, step 10)** | The disease link, without liftover into new sequence: "the complete reference changes which gene a GWAS locus points to" | Low. Summary statistics plus coloc, after step 5 |
| Core | **Geuvadis replication (step 9)** (roles swapped 2026-10-17: MAGE is discovery) | Independent replication in the same cell type; reviewers will expect replication; also the dataset FixItFelix used, so the patched-region results are directly comparable | Medium. ~1,700 CPU h through the same pipeline |
| Core (decided 2026-10-17) | **Third reference: FixItFelix's patched GRCh38 (step 3b)** | Answers the obvious reviewer question "would patching GRCh38 be enough?" with data rather than argument; turns the comparison into GRCh38 → patched → complete | Medium. +50% alignment (~1,700 CPU h for MAGE), 0.9 GB download; no new genotypes (GRCh38 calls apply) |
| Core | **Long-read validation in GM12878** (Fiber-seq FIRE; HG001 ONT) | Long reads place uniquely in SDs and much new sequence, so they directly test whether T2T-only short-read peaks are real, and add haplotype-resolved accessibility at T2T-only eQTLs | Low to medium, depending on what is downloadable |
| Core (data in hand 2026-10-23) | **Paralog copy number as a covariate or genotype** (ctyper) | Separates "T2T measures the gene better" from "people differ in how many copies they have"; directly addresses the main competitor | Low. Precomputed ctyper copy numbers for all 3,202 1000 Genomes samples (84,346 allele types; `data/external/ctyper-ma2025/`). Still to do: map their gene groups (GRCh38 annotation) to our GeneIDs. Running ctyper ourselves would be ~1.3 CPU h and 20 GB per genome, on GRCh38 CRAMs |
| Optional | **sQTLs** (splicing QTLs), as in MAGE | Splicing in paralogs and corrected genes is also reference-dependent; a second molecular trait | Low once alignments exist (LeafCutter or similar) |
| Optional | **SV- and VNTR-eQTLs** (1KGP ONT 1,019; HGSVC3) | Many causal variants in SDs are structural; tests whether T2T-only SNP eQTLs tag SVs | Medium. Calls must be on, or carefully placed onto, T2T |
| Optional | **Sequence-model variant scoring** (e.g. ChromBPNet trained on GM12878 DNase/ATAC aligned to T2T) | Predicts the effect of T2T-only variants on accessibility; a computational component that suits ISMB | Medium to high. Needs a GPU; training on T2T alignments not yet checked for precedent |
| Optional | **Personal diploid reference for NA12878** (HG001 near-T2T assembly) for allele-specific binding | Removes reference-allele mapping bias in step 8 more thoroughly than WASP | Medium, depending on assembly availability |

**Not recommended:** heritability partitioning (S-LDSC) with T2T annotations. LD scores and GWAS are built on GRCh38 panels that exclude new sequence, so the result would mostly reflect GRCh38 again. Colocalisation answers the disease question more directly.

## 7. Go/no-go points

1. **After step 2 (5 samples):** do the expected loci show clear expression differences between references? Candidates: the chr21 genes falsely duplicated on GRCh38 (e.g. CBS, U2AF1, KCNE1), NBPF and NOTCH2NL on 1q21, and other SD paralog families. **No visible difference at these candidates → stop; return to the `t2t-atlas` methods-and-benchmark paper.** Criteria in section 13.6.
2. **After step 6:** are there enough T2T-only and GRCh38-only eGenes, with mechanisms assigned, to support a paper?
   - The count that matters is the eGenes changed **outside FixItFelix's patched regions**, after removing annotation and MHC artefacts. Changes inside the patched regions replicate published work.
   - A focused result ("tens to a few hundred genes beyond the known errors, mostly in SDs, explained by X") is acceptable.
   - "Only the known chr21 genes change" or "almost nothing changes" is not a paper on its own.

## 8. Risks and controls

| Risk | Control |
|---|---|
| Genotype pipeline differs between CHM13 and GRCh38 calls, so the comparison measures callers rather than references | Use the 1KGP GRCh38 calls whose pipeline the CHM13 recall was modelled on; harmonise filters; report the residual difference |
| Annotation confound: richer T2T gene models create apparent gains | Same annotation pipeline and release on both (RefSeq RS_2025_08); main comparison on genes with a GeneID on both references; genes on one reference only as a separate analysis |
| Reference assemblies differ in more than the intended way (chrM, EBV, chrY PAR masking, unplaced contigs) | Make chrM, chrEBV and chrY PAR masking identical (section 13.2); GRCh38 unplaced contigs stay, as part of GRCh38 as it is used |
| Short reads cannot resolve satellites or near-identical paralogs | Report unresolved genes as missing, not zero; quantify both unique-only (STAR) and multi-mapping-aware (salmon); k-mer uniqueness from `t2t-atlas` step 6 |
| CHM13 is one haploid genome; copy-number-variable SDs differ between people | State as a limitation; flag genes in CNV SDs; pangenome follow-up |
| Population structure (Geuvadis includes YRI) | Genotype PCs as covariates; sensitivity analysis on the European subset |
| Small effect overall (98% of GM12878 peaks were shared between references) | Go/no-go 1 before full compute; frame as targeted to SDs and corrected regions |
| Prior or concurrent work | FixItFelix (Behera 2023) already showed artefactual GRCh38 eQTLs at the chr21 false duplication in Geuvadis (2.2): cite it, present chr21 as replication, and lead with what lies outside its patched regions. Recheck the literature before submission; contrast explicitly with the Jan 2026 SD preprint and Shew et al. 2026 |
| ctyper (2025) already explains SD-gene expression with paralog copy number | Read it before step 3; position this study as the T2T-vs-GRCh38 comparison of SNP eQTLs plus the regulatory-element link, and use paralog copy number as a covariate rather than competing with it |

## 9. Compute and storage

- Each cohort's FASTQ is about 2.5 TB; never stored at once (batches of 60 runs, 13A).
- From the pilot ledger: about 7–9 min per sample per reference at 16 CPUs (~2 CPU h).
  - MAGE (779 runs × 3 references): ~4,700 CPU h (measured ~6 CPU h per run in the test), about 2 days at 96 cores.
  - Geuvadis (652 runs × 3 references, shorter reads): ~2,500 CPU h (estimate).
- Scratch had 2.8 TB free (shared) on 2026-10-17.
- **`/storage/projects` was 97% full on 2026-10-17 (5.4 TB free, shared by the group).** Keep only counts, quantifications and logs there; genome BAMs stay on scratch (4.0 TB free) or are not kept.
- STAR indexes: about 30 GB per reference, kept in `data/star-references/`.

## 10. Expected outputs and venue

- **Headline (to test):** "Patching known GRCh38 errors fixes the tip; a complete reference changes expression and eQTLs at many more loci." Supported by the three-way comparison (GRCh38 → patched GRCh38 → T2T) and by region classes.
- Catalogue of eQTLs and eGenes on T2T and GRCh38 with per-pair mechanism and region labels; T2T-native LCL eQTL summary statistics (MAGE, Geuvadis) as a resource.
- Replication of FixItFelix's chr21 artefactual eQTLs (PWP2, GATD3A) in an independent, larger cohort.
- Case loci linking a T2T-only eQTL, a GM12878 regulatory element and allele-specific binding.
- Venues: ASHG 2027 (abstract) or ISMB 2027 (Copenhagen, 18–22 Jul 2027; deadline not found, usually January). RECOMB 2027 (abstract 28 Oct, paper 4 Nov 2026) is too soon and expects a method contribution. Combining this with the `t2t-atlas` k-mer confidence method would strengthen an ISMB submission.

---

## 11. Execution rules

This project follows `~/research/t2t-atlas/docs/execution.md` sections 7–9 (scripts and software, running jobs on Karakoram, data handling) unless this section says otherwise. Cluster facts: `~/research/t2t-atlas/docs/karakoram-hpc-summary.md`. If a rule blocks necessary work, change it here first (section 16), then do the work. The inherited rules in short:

- Scripts are grouped by step (`scripts/NN-name/`), named `<verb>-<object>.<ext>`, safe to rerun, exit non-zero on failure.
- Software: **pixi only** since 2026-11-17 (`pixi.toml` + `pixi.lock`; environments `default`, `bio`, `longread`, `docs`, `lint`; see 11.5). Results before that date used cluster modules of the same versions.
- All compute through `sbatch`/`srun`, with `--mem` and `--time` (measured runtime + ~50%). Stay within **128 cores and 400 GB RAM in total**; cap arrays (`%K`). At most 6 parallel download streams.
- Slurm accounting is disabled: wrap the main command in `/usr/bin/time -v` and log `date -u` at start and end. Append every job to `data/metadata/resource-usage.tsv`; `alloc_cpu_hours` = CPUs × elapsed hours. Lower requests that use less than 50%.
- Slurm logs: `logs/<step>/<job-name>_<jobid>.{out,err}` (relative paths; submit from the project root with `scripts/submit.sh`, 11.5).
- Raw reads are immutable and read-only after the md5 check; scratch is disposable.

### 11.1 Layout

Code and documents in the home directory, bulk data on project storage, temporary files on scratch. Since 2026-11-17 the four links below are made by `pixi run setup` from `config/site.env` and are not tracked by git (11.5).

```
/storage/home/yusuf.maulana/research/t2t-eqtl/     project root (small files only)
├── docs/t2t-eqtl.md   this file (all documentation)
├── config/            sample lists (samples/), pipeline samplesheets (samplesheets/), chromosome-name maps
├── scripts/           one numbered folder per step (01-genotypes/ … 07-presentation/; see 11.4)
├── pixi.toml, pixi.lock, .pixi/   environment (python 3.11, pandas, gffread, liftOver, plink2, pytorch-cpu, R qvalue; see 13B.2)
├── envs/vendor/       tensorQTL 1.0.10 and qtl 0.1.10 from PyPI wheels (on PYTHONPATH via pixi.toml)
├── envs/              container records
├── data     -> /storage/projects/yusuf.maulana/t2t-eqtl/data
├── results  -> /storage/projects/yusuf.maulana/t2t-eqtl/results
├── logs     -> /storage/projects/yusuf.maulana/t2t-eqtl/logs     (Slurm logs by step)
└── work     -> /scratch/yusuf.maulana/t2t-eqtl/work               (FASTQs, BAMs; disposable)
```

```
config/
├── samples/expression-individuals.tsv   step 1: every MAGE/Geuvadis individual, cohort(s), in T2T panel, included (1,089)
├── samples/pilot-5.tsv            step 2: run, sample, ancestry, sex, lab, read pairs, ENA URLs and md5s
└── chrom-names-chm13v2.0.tsv      UCSC (chr1) <-> RefSeq (NC_060925.1) names, plus chrM, chrEBV
data/
├── external/
│   ├── 1kgp-chm13v2.0-phased-shapeit5-v1.1/   phased panel + index + README, read-only
│   └── t2t-atlas-analysis -> t2t-atlas results/analysis (GM12878 peaks), read-only
├── metadata/
│   ├── geuvadis/ena-PRJEB3366-2026-10-13.tsv            ENA run report (run, sample, layout, reads, bytes)
│   ├── geuvadis/ena-PRJEB3366-fastq-md5-2026-10-17.tsv  ENA run report with FASTQ URLs and md5s
│   ├── geuvadis/E-GEUV-1.sdrf.txt                       Geuvadis sample sheet (sex, ancestry, lab)
│   └── resource-usage.tsv                               job ledger, same columns as t2t-atlas
├── references/        links to t2t-atlas: chm13v2.0-ncbi, grch38-noalt, t2t-annotations, liftover
└── star-references/<chm13v2.0|grch38|grch38-ffx>/   stage 02: genome.fa, annotation.gtf, star-2.7.11b-oh75/ and -oh149/
                                          (STAR indexes for 76- and 150-bp reads), transcripts.fa, genes.tsv, tx2gene.tsv,
                                          keep-regions.bed, CHECKSUMS*.md5
results/
├── 01-genotypes/<chm13v2.0|grch38>/   genotypes of the 1,089 expression individuals on each reference: per-chromosome
│                                   BCFs (expression/, na12878/), counts-by-chrom.tsv, summary.tsv
└── 03-rnaseq-pilot/<ref>/<run>/   STAR logs and gene counts, salmon/, idxstats; comparison/
work/03-rnaseq-pilot/              fastq/<run>/, bam/<ref>/<run>.bam
```

### 11.2 Dependencies on `t2t-atlas`

- References, annotations, chain files and the GM12878 peak results are **read through links and never modified here**. Moving or renaming them in `t2t-atlas` breaks this project. Derived references (e.g. for STAR) are built here, in `data/star-references/`.
- **Chromosome names differ.** The 1KGP panel and the T2T annotation tracks (SD, censat, RepeatMasker) use UCSC names (`chr1`); `t2t-atlas`'s `chm13v2.0-ncbi` reference and its RefSeq GTF use `NC_…`. Convert explicitly with `config/chrom-names-chm13v2.0.tsv` (from `liftover/hs1.chromAlias.txt`); never assume they match. GRCh38 uses `chr1` throughout.
- The panel has no chrY or chrM. The NCBI CHM13v2.0 assembly has chrY (from HG002) but no chrM and no EBV.

### 11.3 Project rules

- **The comparison is controlled.** Both references are processed from the same raw reads with the same pipeline, tool versions and gene-model definitions. Any intended difference is written down in the step log.
- **Liftover is used only to pair results** computed independently on each reference, never to create input data.
- **Downloads from S3 and the 1000 Genomes FTP must retry and resume** (HTTPS resets on most TLS handshakes from Karakoram; plain HTTP is refused by the bucket). Use `curl -C - --retry … --retry-all-errors` as in `scripts/01-genotypes/chm13v2.0/fetch-panel.sh`. ENA (`ftp.sra.ebi.ac.uk` over HTTPS) and BioStudies work normally.
- **Geuvadis FASTQs are never all on disk at once:** per sample, download → verify (NCBI SRA: `vdb-validate` plus read-pair count equal to ENA's; ENA: md5) → align to both references → quantify → verify outputs → delete FASTQ. Record each deletion.
- Every job sets `TMPDIR` under `$T2T_TMP` with `t2t_tmpdir` (`scripts/lib/common.sh`; on Karakoram `/scratch/$USER/t2t-eqtl/tmp/`, set in `config/site.env`). `/tmp` is the node's small root disk and was filled on 2026-10-13. Every job is appended to `data/metadata/resource-usage.tsv`.
- Python scripts run with `pixi run --manifest-path <root>/pixi.toml python <script>`; add packages only with `pixi add`.

### 11.4 How the code fits together (written 2026-10-27)

A walk through the scripts in run order: what each step answers, what it reads and writes, and why it is built the way it is. Details and results are in sections 13–13B.

**The idea behind every step: change only the reference.** Same individuals, raw reads, tools and annotation release on every reference. Much of the code exists to remove differences that are *not* the reference: chrM/chrEBV/PAR fixes, genes matched by GeneID, shared genotype PCs, matched variants.

**Folders and step numbers.** Script, result and log folders share one number per stage. Sections 5, 12, 13 and 13A–13B were written earlier and keep the original step numbers:

| Folder (`scripts/`, `results/`, `logs/`) | Stage | Called in sections 5, 12, 13–13B |
|---|---|---|
| `01-genotypes/{chm13v2.0,grch38}` | genotypes on both references | steps 1 and 4 |
| `02-references` (results in `data/star-references/`) | STAR references, indexes, third reference, kept regions | parts of steps 2, 3 and 3b |
| `03-rnaseq-pilot` | five-sample pilots (go/no-go 1) | steps 2 and 2b |
| `04-expression` | expression for whole cohorts | steps 3 and 3c |
| `05-eqtl` | cis-eQTL mapping | step 5 |
| `06-compare` | comparison between references | step 6 |
| `07-presentation` | figures, workbook, deck, manuscript | step 7 |

**Stage 01, genotypes on both references (`scripts/01-genotypes/`; before 2026-10-27 the GRCh38 half was a separate step 4).** One sample list, then one subfolder per reference with the same two scripts, so both genotype sets are made the same way.
- `make-sample-lists.sh` (login node, seconds): MAGE and Geuvadis individuals (ENA run reports) checked against the T2T panel's sample names → `config/samples/expression-individuals.tsv`, one row per individual with its cohort(s), whether it is in the panel, and `included` (1,089 = 731 MAGE + 451 Geuvadis − 93 in both; 13 Geuvadis individuals lack T2T genotypes). Read by `make-cohort-sheets.py` and both subset scripts.
- `chm13v2.0/fetch-panel.sh` (login node): phased 1KGP panel on CHM13 (3,202 individuals), `curl -C - --retry` in a loop, because S3 resets most TLS handshakes from Karakoram; size checked against S3.
- `grch38/fetch-panel.sbatch`: the NYGC phased GRCh38 panel, per chromosome, md5 against the manifest.
- `<ref>/subset-genotypes.sbatch`: per chromosome, keep the `included` individuals, split multiallelic sites (eQTL tools test one ALT per row), recompute MAF in these individuals, drop monomorphic sites; also NA12878 (for allele-specific binding) and variant counts. The T2T script counts by `SYNTENIC` flag (absent = T2T-only sequence, the H1 candidates); the GRCh38 script also drops SVs (`-v snps,indels`), which the T2T panel lacks. Output: `results/01-genotypes/<ref>/expression/<chr>.bcf`. Both GRCh38 references use the `grch38` set.

**Stage 02, references (`scripts/02-references/`).** Everything reads are aligned to and counted against, built once and used by every later stage (pilot, expression, eQTL phenotypes, comparison, figures). Output: `data/star-references/<chm13v2.0|grch38|grch38-ffx>/`.
- `prepare-star-references.sbatch` (array: chm13v2.0, grch38): makes the two genomes comparable first (13.2: chrM and chrEBV added to T2T, T2T chrY PARs masked as on GRCh38), then writes `genome.fa`, `annotation.gtf`, `genes.tsv` (one row per gene with its NCBI GeneID, the key for matching genes across references), `tx2gene.tsv` (transcript → gene, for salmon), `transcripts.fa` (gffread, for salmon) and the STAR index for 76-bp reads, `star-2.7.11b-oh75/`.
- `build-star-index.sbatch <overhang>`: a STAR index for another read length, e.g. `star-2.7.11b-oh149/` for MAGE's 150-bp reads.
- `prepare-ffx-reference.sbatch`: the third reference, FixItFelix's patched GRCh38 (`grch38-ffx/`; md5 and contig checks, GRCh38 gene models, both indexes).
- `make-region-beds.py` → `keep-regions.bed`: candidate loci and the MHC, whose alignments each run keeps (genome BAMs are not kept).

**Stage 03, pilot (`scripts/03-rnaseq-pilot/`; go/no-go 1).** Does anything change, and do ordinary genes stay the same?
- `make-pilot-sheet.py`, `make-mage-pilot-sheet.py` → `config/samples/*pilot-5.tsv`.
- `fetch-sra.sbatch`: `prefetch` → `vdb-validate` → `fasterq-dump`; read pairs must equal ENA's count.
- `align-quantify.sbatch`: STAR GeneCounts (unique reads only) and salmon (ambiguous reads shared by EM). When both methods agree, the reference moved the reads; when only STAR changes, it is usually annotation ambiguity. This became the step-6 rule.
- `compare-expression.py`: join by GeneID → `results/03-rnaseq-pilot/comparison-geuvadis/`.

**Stage 04, expression for whole cohorts (`scripts/04-expression/`).** A cohort is ~2.5 TB of reads on a shared scratch, so reads are never all on disk: download a batch, process, delete.
- Once: `make-cohort-sheets.py` → `config/samplesheets/<cohort>.csv` (one row per run, batches of 60). The references, indexes and `keep-regions.bed` come from stage 02.
- `submit-cohort.sh <cohort>` chains, per batch, `fetch-sra-batch.sbatch` (download and validate, waits while scratch < 400 GB free, 30-min limit per `prefetch` attempt) and `run-batch.sbatch` (Nextflow head, then deletes the `.sra` of completed runs). fetch_k waits for run_{k-2}, so at most about two batches of `.sra` exist; all dependencies are `afterany`.
- `nextflow/modules/local/quantify_run/main.nf`, one task per run: `fasterq-dump` → read-pair check → for each reference STAR → salmon → idxstats → `regions.bam` → delete the genome BAM → `complete.txt` last. Both parts skip runs that have `complete.txt`, so rerunning a batch is always safe (that is the resume mechanism).
- `merge-expression.py --cohort <c>` → `results/04-expression/<c>/matrices/<ref>/{star-counts,salmon-counts,salmon-tpm}.tsv.gz` (genes × individuals; runs summed).

**Stage 06 preparation, needed by stage 05 (`scripts/06-compare/`).**
- `match-variants.sbatch` + `match-variants.py`: pair T2T and GRCh38 variants (lift each T2T variant, match position and alleles, allow swapped alleles, pair nearby indels, require MAF agreement). Liftover is used only here, only to pair independently called variants. The GRCh38 panel has ~1.2 M common SNPs the T2T panel lacks, mostly removed by the T2T panel's filters, so the main comparison uses matched variants only (13B.1). Output: `results/06-compare/variant-matching/<chr>.tsv.gz`.
- `make-gene-classes.py` → `gene-classes.tsv`: region class per GeneID (`ffx_patched` > `t2t_new_sequence` > `sd` (≥ 98% SDs) > `other`) and MHC/PAR/chrM/rRNA flags. This label tests H0.
- `compare-annotation-overlaps.sh` → `annotation-overlaps.tsv`: genes whose exon-overlap partners differ between references (STAR-only artefacts such as HCP5).

**Stage 05, eQTL mapping (`scripts/05-eqtl/`).** Per cohort: genotypes once, then per reference and method (star | salmon):
- `prepare-genotypes.sbatch <cohort> <t2t|grch38>` → pgen sets `.all` and `.matched`, IDs `chr:pos:ref:alt` (the variant-matching IDs).
- `prepare-phenotypes.py` → `phenotypes.bed.gz`: GTEx recipe (expressed filter, TMM CPM, inverse-normal, TSS); T2T `NC_…` names converted to `chrN`.
- `make-covariates.py`: 5 genotype PCs computed **once** from matched T2T SNPs and shared by all references (ancestry correction cannot differ); sex; expression PCs per reference. `prepare-phenotypes-covariates.sbatch` runs both.
- `run-tensorqtl.py` / `.sbatch`: per chromosome, `map_cis` (permutations: is the gene an eGene?) and `map_nominal` (every pair, so step 6 can look up any variant).
- `collect-eqtl.py`: q-values → `cis-eqtl.tsv.gz`.
- `make-simulated-cohort.py`: the `simtest` cohort with planted eQTLs, to test the pipeline before trusting it (13B.2).

**Stage 06, comparison (`scripts/06-compare/`): a funnel.**
1. `compare-eqtl.py` (via `compare-eqtl.sbatch`): eGene on A, on B; reference-specific only if the lead variant also fails on the other reference. Presence/absence misses the main artefact (PWP2 is an eGene everywhere, far too strong on GRCh38), hence step 2.
2. `compare-effects.py`: the same variant on both references; allele coding aligned by the sign of the genotype correlation (the matching table's "swapped" label was wrong for some indels, 13B.3); z_diff, Bonferroni → changed.
3. `combine-methods.py`: changed under STAR **and** salmon with the same reference stronger → `robust-matched/changed-<A>-vs-<B>.tsv` (the 57-gene H0 set).
3b. `paired-test.py`: the same pairs with a paired test that accounts for the shared individuals (z_paired), with bootstrap calibration, an optional ctyper copy-number covariate (`ctyper-copy-number.py`, `--copy-number`) and replication in another cohort (`--at`). The method is described in full in 13B.4.
4. `read-destination.py`: which gene's per-person change mirrors it (paralog swaps, reads taken by T2T-only genes).
5. `reversal-variant-location.py` (+ `.sbatch`): true reversals and where their variant lies.
6. `h1-variant-qc.py` (+ `.sbatch`, on `--variants all`): are the T2T-only lead variants well genotyped?

**Stage 07, presentation (`scripts/07-presentation/`).** `prepare-figure-data.py`, `pipeline.dot` (Fig 1), `make-figures.R` (Figs 2–13), `build-workbook.py`, `build-deck.py`; manuscript `.docx` by pandoc from `results/07-presentation/manuscript/t2t-eqtl-manuscript.md`.

**Data flow.**

```
config/samples ─► 01 chm13v2.0/subset-genotypes ─► results/01-genotypes/chm13v2.0 ─┐
                  01 grch38/subset-genotypes    ─► results/01-genotypes/grch38    ─┼─► 06 match-variants
02 references ─► data/star-references/<ref>/ (index, genes.tsv, transcripts.fa) │
config/samplesheets ─► 04 Nextflow ─► results/04-expression/<c>/<run>/   │
                       merge-expression ─► matrices/<ref>/               │
                                                                         ▼
     05 prepare-genotypes + prepare-phenotypes + make-covariates ─► run-tensorqtl ─► collect-eqtl
     06 make-gene-classes + compare-annotation-overlaps (labels) ─┐               │
                                                                  ▼               ▼
        compare-eqtl ─► compare-effects ─► combine-methods ─► read-destination / reversal / h1-qc
                                                                  ▼
                                                     07 figures, workbook, deck, manuscript
```

**Rules that recur in the code:**
1. Remove every difference except the reference.
2. A change counts only when STAR and salmon agree.
3. Compare effect sizes at the same variant, not eGene counts.
4. Keep a positive control (GRCh38 vs patched GRCh38 must reproduce FixItFelix).
5. Check labels against raw data (read counts against ENA, md5s, allele signs against genotype correlation).
6. Every step is safe to rerun (completion markers, skip finished work, resumable downloads).

### 11.5 Reproducibility set-up (branch `capstone`, 2026-11-17)

What makes the project run outside Karakoram, following the practices of the course "The Art of Reproducible Research" (CBIO7501). Nothing here changes a result: the default environment's packages are unchanged apart from removing python-pptx (the talk-deck builder is not part of this branch), and the bioinformatics tools are pinned to the module versions the results were made with.

- **No machine-specific paths in code.** Scripts find the project root themselves (`SLURM_SUBMIT_DIR` for jobs, git for shell scripts, `T2T_ROOT` to override) and source `scripts/lib/common.sh`, which reads `config/site.env` (not in git; template `config/site.env.example`): storage folders, `T2T_SCRATCH`, `T2T_TMP`, Slurm partitions (`T2T_PARTITION`, `T2T_PARTITION_LONG`) and account. Python and R scripts already used paths relative to the root. `scripts/00-setup/lint.sh` fails if `/storage/`, `/scratch/` or `/home/…` appears in a script.
- **Storage links.** `data/`, `results/`, `logs/`, `work/` are git-ignored; `pixi run setup` (`scripts/00-setup/link-storage.sh`) links them to the folders in `config/site.env` (or makes plain folders), plus the t2t-atlas inputs (`T2T_ATLAS_DIR`). After switching from `main` (which tracks the four links) to this branch, git removes the links; rerun `pixi run setup`.
- **Slurm.** `scripts/submit.sh <job.sbatch> [args]` (or `pixi run submit …`) submits from the root with the site's partition and account and creates the log folder. Batch scripts carry no partition; one that needs the long partition says `# t2t-partition: long`. `#SBATCH --output` paths are relative (`logs/<stage>/`). The cohort submitters (`submit-cohort.sh`, `submit-eqtl-cohort.sh`) go through `submit.sh`.
- **Software.** `module load` replaced by `t2t_env bio`, which puts the pixi environment `bio` first on PATH (STAR 2.7.11b, salmon 1.11.4, samtools 1.23.1, bcftools 1.23.1, bedtools 2.31.1, sra-tools 3.4.1, minimap2 2.31, Nextflow 25.10.4; a separate environment because these builds need older system libraries than the R/Python stack). Only PATH is changed, so `pixi run python …` in the same script still uses the default environment. Nextflow tasks get the same PATH through `beforeScript` (`conf/base.config`); `-profile modules` restores the cluster modules. Clair3 2.0.3, WhatsHap 2.8, samtools/bcftools 1.25 pinned in `longread`; Quarto 1.10.19 in `docs` (Quarto pins its own pandoc). The figure font (TeX Gyre Heros 2.004) is installed by `pixi run fonts` (sha256-checked; files identical to the ones used); `make-figures.R` falls back to sans if it is missing.
- **Samplesheets.** The `sra` column is relative to the root (`work/04-expression/sra/…`); Nextflow resolves it against `params.root`, which is now derived from the pipeline's location.
- **Data.** `config/data-sources.tsv`: every input with origin, URL or accession, version or date, local path, the script that fetched it and the published checksum. `pixi run checksums` writes sha256 sums of `data/external/` and `data/references/` to `config/checksums/inputs.sha256`; `pixi run verify-data` checks them.
- **One command per output.** pixi tasks: `setup`, `fonts`, `checksums`, `verify-data`, `submit`, `figure-data`, `figures`, `diagrams`, `workbook`, `reproduce` (figures, diagrams and workbook from the stage 05–06 results), and `-e docs data-journey`.
- **Checks.** `pixi run -e lint lint` and GitHub Actions (`.github/workflows/check.yml`): lock file current, Python and shell syntax, ShellCheck errors, no hard-coded paths, data-sources table well formed.
- **Sharing.** `LICENSE` (MIT, code), `LICENSE-docs.md` (CC BY 4.0, documents and figures; vendored tensorQTL keeps its own licence), `CITATION.cff`, README rewritten as the entry point. Rendered HTML is not tracked: `docs/data-journey/data-journey.html` is rebuilt with `pixi run -e docs data-journey` (to be published as a release asset or on GitHub Pages), and the old `docs/stage-04-expression.html` export (no source, not referenced) was removed from git.

## 12. Step log

| # | Step | Status | Where |
|---|---|---|---|
| 0 | Set up project; move eQTL work out of `t2t-atlas` | done 2026-10-13 | this file |
| 1 | Genotypes: download the phased 1KGP T2T-CHM13v2.0 panel; subset to Geuvadis and NA12878; count variants by class | done 2026-10-13. Panel (SHAPEIT5 v1.1, 3,202 samples, 12,949,615,208 bytes, size matches S3) downloaded with `fetch-1kgp-chm13-panel.sh` (now `scripts/01-genotypes/chm13v2.0/fetch-panel.sh`). 451/464 Geuvadis samples in panel (`config/samples/`; lists made by hand, reproduced exactly by `make-sample-lists.sh` on 2026-10-17). Job 356292 (3.4 min, 23 CPUs): MAF ≥ 0.05: 7.39 M syntenic, 17.5 k T2T-only variants; NA12878 5.9 k het sites in T2T-only sequence | `scripts/01-genotypes/`, `data/external/1kgp-chm13v2.0-phased-shapeit5-v1.1/`, `results/01-genotypes/chm13v2.0/` (renamed 2026-10-27, section 16) |
| 2 | Five-sample RNA-seq pilot on both references (go/no-go 1) | **done 2026-10-17**; result supports go (13.7). Decisions made, metadata fetched, scripts written, pixi environment created. Submitted 09:12 UTC: references 360122 (array 0–1), FASTQ download 360123, alignment 360124. 360123/360124 cancelled 09:14 (ENA refused most of 10 parallel HTTPS streams); download cut to 2 streams with longer backoff and resubmitted as 360127. 360122_0 failed after 1 min: gffread rejects RefSeq gene records (`transcript_id ""`); fixed by passing it only transcript/exon/CDS records, tested on both GTFs; 360122_1 and alignment 360128 cancelled. Resubmitted: references 360138, alignment 360139. 360138 done: GRCh38 09:52 (20 min, 52 GB peak; 59,257 genes, 195 contigs), T2T 09:54 (29 min, 48 GB; 58,600 genes incl. 37 chrM, 26 contigs, chrY PAR1 verified all N). Indexes 37 + 30 GB. ENA download 360127 got one connection at ~0.75 MB/s (~8 h for the pilot); NCBI SRA tested at ~12 MB/s, so at the user's choice 360127/360139 were cancelled, partial ENA files deleted, and reads fetched from NCBI SRA instead (`fetch-sra.sbatch`, job 360219); alignment 360220 (after 360138 and 360219). 360219 done 10:14 (18 min, 3 runs in parallel, no S3 resets): 5 runs, 21.2 GB FASTQ, read pairs equal to ENA for all five. 360220 done 10:52 (10 tasks, 5.3–9.9 min each; peak RSS STAR 32–33 GB, salmon ≤ 12 GB). `compare-expression.py` run 10:56. 360122 waits for memory: the `t2t-atlas` ChIP-seq tasks held 1,368 of 1,470 GB | section 13; `scripts/03-rnaseq-pilot/` |
| 2b | Five-sample MAGE pilot (+ HG00265 in both cohorts) | **done 2026-10-17**: effect replicates (13.8). Jobs 360465 (MAGE download, 22 min), 360466 (ERR188222), 360467 (MAGE alignments, 10 × 8–9 min), 360468 (ERR188222 alignments) | section 13.8 |
| 3 | Full expression: MAGE (discovery) then Geuvadis (replication), three references | 2026-10-17: pipeline restructured (Nextflow DSL2: entry, workflow, module, profiles) with CSV samplesheets; 150 bp STAR indexes 360593; 2-run test 360624 passed; **MAGE submitted 15:25; batches 2–13 resubmitted 17:16 after a download fix (13A.5)** | section 13A; `scripts/04-expression/`, `config/samplesheets/` |
| 3c | Full expression, Geuvadis (replication): 652 runs, 451 individuals, 2×76 bp, 11 batches, three references (oh75 indexes) | **submitted 2026-10-25 10:00 UTC** with `submit-cohort.sh geuvadis`: fetch/run pairs 365827/365828 … 365847/365848 (list in `logs/04-expression/submitted-geuvadis-20261025T100000Z.tsv`); scratch 5.2 TB free at launch. **2026-10-26: scratch full again** (t2t-atlas 4.7 TB): batches 5 and 6 lost tasks ('terminated by external system'; 63 runs without results), batch-7 download paused by the 400 GB guard at 297/652 runs. With the user's approval deleted our Nextflow work (MAGE 90 GB, Geuvadis 83 GB), the .sra of the 57 completed Geuvadis runs (152 GB) and tmp/vmtest + old-docs; scratch 125 → 449 GB free. The 63 unfinished runs keep their .sra. Download resumed 10:01 UTC; batches 5 and 6 resubmitted (user's approval) as 370298 → 370299 to process the 63 runs; both completed (15 + 48 runs). Batches 7–10 completed 2026-10-26/27 (60 runs each, 0 failed). **Batch 11 (365848) stopped at 629/652**: ERR204978 (NA20527) failed the read-count check twice, and Nextflow then aborted the batch, 4 tasks running. ERR204978 has 41.6 M spots, but 5.27 M have one empty mate. `fasterq-dump` writes 36,380,059 complete pairs to _1/_2 plus 5.27 M orphan reads (75 bp), confirmed by inspection job 370916. Batch 11 resubmitted without it (`config/samplesheets/geuvadis-without-ERR204978.csv`, job 370911). ERR204978 is quantified from its complete pairs only, as batch 12 (`config/samplesheets/geuvadis-ERR204978-paired-only.csv`, expected 36,380,059 pairs; job 370949) **2026-10-27, ERR204978 (NA20527, batch 11):** failed the read-pair check twice (Nextflow tasks at 05:33 and 05:58 UTC). An ad hoc job (370916, `fasterq-dump --split-3`) showed 36,380,059 pairs plus 5,266,444 unpaired reads (together ENA's 41,646,503), all 75 bp; the pipeline's `--split-files` check compares pairs with ENA's total and so fails. Fix: batch 11 rerun without it (`config/samplesheets/geuvadis-without-ERR204978.csv`, job 370911), and the run processed on its pairs only as batch 12 (`config/samplesheets/geuvadis-ERR204978-paired-only.csv`, expected 36,380,059 pairs; job 370949, done 08:54, 24 min, 5.9 CPU h); the 5.3 M unpaired reads are not used. | section 13A |
| 3b | Third reference: FixItFelix patched GRCh38 | decided 2026-10-17 (user); reference built (360623); in the MAGE run | section 13A.3 |
| 4 | Genotypes for 1,089 individuals: T2T subset (360598) and GRCh38 NYGC panel download (360594) | started 2026-10-17 | section 13A.2 |
| 5 | cis-eQTL mapping, MAGE (731 individuals), chr1–22, STAR counts, matched variants, three references | **done 2026-10-25**. Genotypes 365584/365585; tensorQTL arrays 365591 (T2T), 365592 (GRCh38), 365593 (patched GRCh38), 22 tasks each, 0 failed, 22:43–00:09 UTC, ~44 CPU h per reference. eGenes (FDR 0.05): 15,013 / 15,259 / 15,224 (MAGE paper: ~15,000). chrX and salmon still to run | section 13B.3; `results/05-eqtl/mage/` |
| 6 | Comparison between references (H0), MAGE | **first pass done 2026-10-25** (`compare-eqtl.py`, 15 min, 9 GB, run on the login node): positive control reproduced (PWP2, CBS, U2AF1); see 13B.3. Same-variant effect sizes (`compare-effects.py`, via `compare-eqtl.sbatch`, job 365661): patch changes 9 genes, all in its regions; T2T vs patched GRCh38 changes 127 unflagged genes, none in patched regions (56 SD). Salmon run: prep 365657 failed (`--wrap` runs `sh`, no pipefail), resubmitted as `prepare-phenotypes-covariates.sbatch` 365662, tensorQTL 365663–365665 (0 failed), comparison 365789. **Read destination done (job 366433; 33/57 with a clear partner). All-variants tensorQTL 366424–366426 done 2026-10-25 19:22 (H1: 59 genes led by T2T-only variants, 9 found only with them; job 366638); H1 lead-variant QC job 370504 on 2026-10-26: 22 solid / 25 borderline / 12 doubtful, 15 solid genes where the T2T-only variant carries signal). Robust set (STAR and salmon agree): 57 unflagged genes changed by T2T vs patched GRCh38, none in patched regions; patch control 8/9** (69 before the allele-sign fix, reruns 366132/366133). 2026-10-26: lncRNA breakdown (9 of 57; ad hoc) and reversed eQTLs: 13 true reversals, lead variant located with `reversal-variant-location.sbatch` (370515 failed, `--wrap` runs sh; 370516, 14 s) | section 13B.3; `results/06-compare/eqtl-comparison/mage/star-matched/` |
| 5–6g | eQTL mapping and comparison, Geuvadis (replication, 451 individuals), chr1–22, STAR and salmon, matched variants | **submitted 2026-10-27 09:45 UTC** after all 652 runs completed (09:24), with the new `scripts/05-eqtl/submit-eqtl-cohort.sh geuvadis`. It chains `merge-expression.sbatch` 371311 → `prepare-genotypes.sbatch` t2t/grch38 371312/371313 → phenotypes and covariates star 371314, then salmon 371322 (serial, so the genotype PCs are computed once) → tensorQTL arrays %2: star 371315/371317/371319, salmon 371323/371325/371327 → `collect-eqtl.sbatch` → `compare-eqtl.sbatch` 371321 (star), 371329 (salmon) → combine-methods 371330. All dependencies are afterok. Job list: `logs/05-eqtl/submitted-geuvadis-*.tsv`. **All done 13:17 UTC, 0 failed**; same-variant replication of the MAGE 57 genes: 50/54 testable genes nominal under both methods, 44 Bonferroni (13B.3). **H1 replication 2026-10-27/28:** STAR all-variants tensorQTL arrays 373575 / 373577 / 373579 (T2T / GRCh38 / patched), collects 373576 / 373578 / 373580, then `replicate-h1.py` (job 373588, after the T2T collect): MAGE's 59 H1 lead variants looked up at the same variant in Geuvadis, with the conditional test on MAGE's best shared variant. Done 2026-10-28 01:45, 0 failed: 13,481 / 13,676 / 13,712 eGenes with all variants; 8 of the 15 solid H1 examples replicate with their own signal (5 distinct variants; 13B.3). **Independent replication 2026-10-28** (job 377254, after a login-node check that the recomputation matches tensorQTL): the 93 individuals also in MAGE dropped, 358 left; results unchanged (54/54 same sign, 49 nominal, 42 Bonferroni; reversals 10/12; H1 9/15 and 8/15; 13B.3). **Paired test and copy number 2026-10-28** (`paired-test.py`, jobs 377313 MAGE, 377411 Geuvadis replication (377360 failed on a print, rerun); `ctyper-copy-number.py` 377320; `paired-test.py --copy-number` 377358): the paired test, calibrated by bootstrap, finds 255 unflagged changed genes T2T vs patched GRCh38 (all 57 included); in Geuvadis's 358 independent individuals 167 of the 193 new genes have the same sign and 122 replicate nominally. ctyper copy number as a covariate leaves 33 of the 34 covered genes of the 57 changed (13B.3). 2026-10-29: coloc with 6 GWAS on both references (job 380275; 8 extent genes colocalise on one reference only, e.g. PTGES3 and IFITM3-AS1); SuSiE-coloc at the 32 loci with a signal (job 380902; IFITM3-AS1, LOC124903630, SERF1B hold, PTGES3 does not; 380440, 380490 failed); ASE realignment (array 380491, 369 runs) and ONT SNP calling (array 380492, 260 individuals) running; GWAS Catalog lookup (job 380065; main 60/117 genes with a proxy hit); gene sets defined (main 117 / core 57 / extent 255; `define-gene-sets.py`, login node, 2 s); ASE feasibility (job 379963): 28 of 57 genes and 4–7 of 13 reversals testable with panel SNPs; the SD reversals have no exon SNPs (13B.3). Done by 10:27 UTC: matrices (59,257 GRCh38 genes × 451 individuals; job 371311, 24 min); genotypes 12,487,313 all / 10,834,722 matched variants (T2T) and 14,801,041 / 10,836,577 (GRCh38); covariates 5 genotype PCs, sex, 60 expression PCs. HG00237 and NA07000 have no sex in the Geuvadis metadata, so it was filled with the mean: harmless for chr1–22, but take it from the 1KGP pedigree before chrX. These jobs ran before the 2026-10-27 folder renames (section 16); the rest of the chain uses only `results/05-eqtl` and `results/06-compare`, which did not move. 2026-10-29: SuSiE-coloc 380902 (2 h 04 min, 8 GB; failed 380440, 380490); ASE realignment array 380491 (369/369 runs, done 2026-10-30 03:35 UTC, 91 GB) and ONT calling array 380492 (260/260, done 2026-10-29 17:58 UTC, 29 GB), 0 failed (13B.3). | section 13B.3 |
| 7 | Presentation materials (oral-talk overview) | **done 2026-10-25**: `scripts/07-presentation/` — `prepare-figure-data.py` (job 366441, 1 min, 9 GB), `pipeline.dot` (Graphviz, Fig 1), `make-figures.R` (ggplot2, Figs 2–10), `build-workbook.py` → `results/07-presentation/t2t-eqtl-project-overview.xlsx` (talk outline, data sources, samples, pipeline, one sheet per figure, key results, 57 genes, glossary, status). pixi: added graphviz, r-ggrepel, r-patchwork, r-data.table, r-scales. LibreOffice is not on the cluster, so formulas are not pre-calculated (set to recalculate on open; logic checked in Python). Talk deck: `build-deck.py` (python-pptx; node is not on the cluster) → `results/07-presentation/t2t-eqtl-talk.pptx`, 18 slides with speaker notes; pixi: added python-pptx, defusedxml. Journal-style write-up (title, abstract, methods with datasets, results with Figs 1-10, discussion, conclusion, next steps, references) as a Claude Doc: https://claude.ai/artifact/SKsKT2bh58uE62Sz41gVBb (references written from memory, DOIs to verify; not maintained after 2026-10-26, the Word file below is the current version). Word version: `results/07-presentation/manuscript/t2t-eqtl-manuscript.docx` (pandoc from `t2t-eqtl-manuscript.md`, 10 figures embedded; pixi: added pandoc). **Updated 2026-10-28** with the Geuvadis replication and H1: `prepare-figure-data.py` (job 375349 with the figures, 3.5 min, 8.5 GB) now writes `replication-h0.tsv`, `replication-h1.tsv` and `effect_class` in `robust-genes.tsv` (stronger on T2T / stronger on patched / reversed / opposite sign, one weak: 13 / 13 / 13 / 18); new Fig 11 (replication, H0 and H1); Fig 8b by effect pattern; workbook: Fig11 sheet, replication table in Key results, talk outline and status updated, 57-gene sheet 'Effect pattern' column; deck 19 slides (new slide 16, Replication in Geuvadis; slide 13, limits and conclusions updated); manuscript draft of 28 Oct with methods 9–10, H1 QC and replication in 3.5, new 3.6 and Fig 11 (11 figures). 'Reversed again' in Geuvadis = opposite-sign slopes with STAR and salmon, all \|slope\| ≥ 0.05 (10/12). **Updated 2026-10-29** with the paired test, gene sets (main 117 headline, core 57, extent 255), independent replication (358), ctyper, GWAS lookup, coloc and SuSiE-coloc, and the ASE status: `prepare-figure-data.py` (job 381190 with the figures, 3 min) adds `gene-sets.tsv`, `paired-effects.tsv`, `calibration.tsv`, `size-floor.tsv`, `replication-paired.tsv`, `copy-number.tsv`, `coloc.tsv`, `gwas-lookup.tsv`, adds `tier` to `effects.tsv`, and now reads the independent (358) replication. Figures renumbered: Fig 7 coloured by gene set; new Fig 8 paired test (z_diff vs z_paired, bootstrap calibration, replication by size floor); old 8–11 → 9–12; Fig 12 left = main set at the 358 with the paired test; new Fig 13 GWAS (coloc PP4 on both references, lookup vs background). The four old PNGs moved to `figures/archive-2026-10-28/`. Deck 23 slides (new: paired test, copy number, GWAS, ASE; numbers updated). Workbook: Fig8/Fig12/Fig13 sheets, new 'Gene sets' sheet (255 genes), Key results with gene sets, calibration, size floor, independent replication, copy number and coloc. Manuscript draft of 29 Oct (13 figures; methods 6 and 10–14 rewritten or new; results 3.3, 3.6, 3.8–3.10 new or rewritten; 25 references, statistics DOIs checked on Crossref). **Updated 2026-11-01**: talk v2 `results/07-presentation/t2t-eqtl-talk-v2.pptx` (`build-deck.py` now writes v2; the 29 Oct deck stays as `t2t-eqtl-talk.pptx`, its script version in git), 24 slides + 3 backup, story reordered: T2T vs GRCh38 first (121 main genes, 6 in FixItFelix regions), the patch alone as sanity check (7, all inside), then T2T vs patched GRCh38 (117, none inside; 114 of the 115 T2T-vs-GRCh38 genes outside the patch remain). New slides: PWP2 in HG00096/NA19717 (613/11/628 reads), steps 1-3 with native charts, backup 'how the tools work'; data sources moved to backup. Plain-language walkthrough with real file rows: `docs/data-journey/data-journey.qmd` (Quarto). Simple journal-style method overview: `scripts/07-presentation/method-overview.dot` → `results/07-presentation/figures/method-overview.{pdf,svg,png}` (Graphviz, left to right, flow only with tool versions, no numbers; a MAGE cohort box (731 individuals) splits into RNA-seq and genotypes (called on T2T and on GRCh38), side by side; expression lane align → eQTLs → compare, genotype lane prepare genotypes → match variants (patched GRCh38 uses the GRCh38 genotypes; left to the caption); each box = name + two short lines (tool and version, one key setting); reference genomes from above; output = one gene list for each of the three pairs (T2T vs GRCh38, patched vs GRCh38, T2T vs patched), marked with the reference colours; dashed replication box below the output (Geuvadis, 451 individuals); ~187 × 80 mm at natural size (placed at 183 mm, ~98%; 7–8 pt text); Okabe-Ito colours; Helvetica: TeX Gyre Heros (Helvetica metrics, GUST e-foundry) in `~/.local/share/fonts`, embedded in the PDF, SVG asks for Helvetica first; MAGE drawn as a borderless grey block with a gap before the split; transparent background; liftOver labelled 'UCSC liftOver v482' (UCSC kent source release 482)); the detailed pipeline stays `pipeline.dot` (Fig 1). Deck v2 slides 3–4 rewritten 2026-11-02: slide 3 = the gap and the research question (H0, plan); slide 4 = prior work (T2T used for variants, ChIP-seq, copy number, duplications, never for eQTLs) ending with FixItFelix, which used T2T only to patch ~20 known GRCh38 errors; the H1–H4 cards moved to a backup slide (28 slides in all); workbook talk-outline rows 3–4 match. Expression slides added 2026-11-02 (deck slides 9–10; 30 slides in all): `prepare-figure-data.py` `expression_tables()` → `expression-compare.tsv` (mean CPM over the 731 MAGE individuals per reference and method; genes annotated on both, one symbol per GeneID, no MHC/PAR/chrM/rRNA; expressed = mean CPM ≥ 1; changed = ≥ 2-fold under STAR and salmon, same direction) and `expression-eqtl-examples.tsv`; `make-figures.R` Figs 14–15 in journal style (183 mm, 7 pt, no titles, PNG + PDF). Results: the patch alone changes 5 of 13,449 genes, all in FixItFelix regions; T2T vs patched GRCh38 changes 97 of 13,563 (99.3% within 2-fold), 31% of them in segmental duplications vs 2.4% of all genes, 7% in T2T-only sequence. Expression change is not eQTL change: of 59 changed genes that were eQTL-tested only 15 change their eQTL, and 89 of the 104 main-set eQTL changes in this table come with < 2-fold expression change. Examples: PKD1P4-NPIPA8 (~3× fewer reads on T2T for everyone, slope 0.74 vs 0.70) and CCZ1 (same mean, slope +0.65 vs −0.25). Run interactively (~17 s, 1.6 GB). Fig 3 (mapping) redone 2026-11-02 in the same journal style (120 × 62 mm, PNG + PDF; theme `nat` and `save_nat()` now at the top of `make-figures.R`); slide 8 numbers are computed from `mapping-qc.tsv`: multi-mapped median 4.3% GRCh38, 4.2% patched, 2.1% T2T (the slide previously said 4.2% for GRCh38); uniquely mapped 84.5% vs 84.9%. Deck v2 reordered 2026-11-02 to follow `docs/data-journey` (32 slides: 27 + 5 backup): introduction (why the reference matters, gap and question, prior work, a new roadmap slide) → data and design (cohorts; `method-overview.png`) → 1 expression (mapping, expression, PWP2) → 2 genotypes (variant matching) → 3 eQTL mapping (eGenes; expression ≠ eQTL) → 4 comparing references (paired test, T2T vs GRCh38, patch alone, T2T vs patched, the 117, gene types, mechanism, examples) → 5 does it matter (Geuvadis, copy number, GWAS, allele-specific expression) → limits, conclusions; the detailed `pipeline.dot` figure moved to the backup. Every slide has a section label (top right) and speaker notes in three parts (what the slide shows / why it matters / where it fits in the story), via `notes3()`. Same day: white background on every slide (title and conclusions no longer navy; emphasis cards light grey), and the workbook talk outline rewritten to the new order (32 rows; column 'why it matters / link in the story'; about 19 minutes). Slides 7 and 9 annotated 2026-11-02: slide 7 says why STAR and salmon (they fail differently; a change needs both); slide 9 defines the duplication label (≥ 50% of the gene in ≥ 98%-identical copies, T2T SD map) and the label order (FixItFelix > T2T-only > duplication > elsewhere). All 7 changed T2T-only genes are also fully in ≥ 98% SDs, so 37 of the 97 changed genes (38%) lie in duplicated sequence. New slide 9 'Two ways of counting reads: STAR and salmon' (2026-11-02; deck now 33 slides) with Fig 16 (`make-figures.R`, journal style, from `expression-compare.tsv`): each gene's T2T/patched change under STAR vs salmon, and counts per comparison. T2T vs patched: 97 genes ≥ 2-fold with both, 114 STAR only, 49 salmon only; patch vs GRCh38: 5 both, 6 STAR only; T2T vs GRCh38: 100 both. Workbook outline row added. PWP2 'one gene, two people' slide moved to the backup (slide 13, expression vs eQTL, makes the same point with new genes); the expression section now ends with a take-home in the notes (27 main + 6 backup slides). Same day, all deck figures redone in journal style and simplified to one idea each (`make-figures.R`: Figs 2, 4–13 rewritten, new Fig 17 three comparisons by region, Fig 17a T2T vs GRCh38, Fig 18 copy number; PNG + PDF), and every figure slide rebuilt as one figure plus 4–5 short bullets (`fig_bullets()` in `build-deck.py`). Fig 8 keeps only the standard-vs-paired z panel (calibration and size floor in the notes); Fig 9 shows the 117 by direction and region; Fig 12 only the main-set replication (H1 in the bullets); Fig 13 only colocalisation. Fig 4 now MAGE-only (`variant_matching_mage()` → `variant-matching-mage.tsv`; chr1–22, MAF ≥ 1% in the 731): T2T 13.51 M (10.53 M SNPs, 2.98 M indels), GRCh38 16.02 M (12.67 M, 3.35 M); matched 11.79 M (9.96 M SNPs, 1.82 M indels); T2T variants not called on GRCh38 1.02 M, frequency differs 0.12 M, T2T-only sequence 0.08 M, never paired because < 1% across the 1,089 people 0.50 M (mostly SNPs; a limitation of pairing on the pooled cohorts), GRCh38 without a matched T2T partner 4.23 M. Run interactively (~100 s, 10 GB). 2026-11-03: deck reorganised into five sections (1 Expression, 2 Genotypes, 3 eQTL mapping, 4 Comparison incl. replication and copy number, 5 Link to GWAS), each opened by a methods slide (`method_slide()`: numbered steps, tools with versions, output); roadmap updated; 32 main + 6 backup slides; workbook outline matches (about 21 minutes). Same day: the comparison methods slide rewritten as a plain how-to of the paired test (one variant per gene → slope on each reference → difference ÷ paired error, which removes the noise shared by the same people → Bonferroni p < 0.05, ≥ 0.1 SD, both methods same direction); the paired-test justification slide (Fig 8: shared noise, standard test 9–15× too cautious, calibration) moved to the backup as 'Backup: why a paired test' (B7), and conclusion 4 now reads 'the changes are real: they replicate (100 of 110) and are not copy number'; 31 main + 7 backup slides; workbook outline renumbered. New backup slide 'how the paired test is calculated' (B8; 31 + 8 slides): Fig 19 (`make-figures.R`, from `paired-example.tsv`, written by `prepare-figure-data.py` `paired_example()` as `paired-test.py` computes it) shows IFITM2 (STAR, T2T vs patched GRCh38): per-person expression after covariates on both references (r = 0.98) and the per-person difference by genotype; the slide gives z = Δ / √(se_A² + se_B² − 2 cov) with the covariance formula, the IFITM2 numbers (z −3.1 independent → −25 paired) and the references (Zellner 1962; Weesie 1999; Lin & Sullivan 2009). Methods slide 17: the change criteria now stand in their own box, 'Call a change: all four must hold' (real: Bonferroni p < 0.05 over about 15,000 genes per method, |z| > 4.6; big enough: |Δ| ≥ 0.1 SD, against a typical eQTL of 0.30 SD, the median |lead slope| of T2T STAR eGenes; not one tool: STAR and salmon, same direction; clean: no annotation, MHC or PAR flag), numbers computed in `build-deck.py`. Same day, sections 4–5 simplified at the user's request (24 main + 12 backup slides, about 17 minutes): after the comparison methods slide come (1) how many eQTLs change, Fig 20 (`make-figures.R`, plain bars replacing the stacked Fig 17 on the slide: T2T vs GRCh38 121, 6 in patched regions; the patch alone 7, all inside; T2T vs patched 117, none inside) with the robustness claim (76 genes at a 0.2 SD floor; 100 of 110 replicate in Geuvadis); (2) the top 10 changed genes as a booktabs-style table (`gt_table()` in `build-deck.py`: gene, type, region, effect on each reference, change, replication; 9 of 10 in SDs, 8 reverse); (3) segmental duplications, Fig 21 (50% of changed genes vs 2.6% of the 14,924 compared); (4) why, Fig 22 (CCZ1/CCZ1B per-person change, r = −0.96, and the CCZ1 eQTL by genotype, +0.65 vs −0.25); (5) one GWAS slide with a table of the 6 main-set genes that colocalise on one reference only (the GWAS methods slide removed); then limits and conclusions (conclusion 2 adds the 0.2 SD count). The old step 1–3, 117-scatter, gene-type and mechanism slides were removed (figures kept on disk); replication, copy number, examples and allele-specific expression moved to the backup. Workbook outline rewritten to match. Same day: both slide tables now made in R with gt (pixi: added r-gt, r-gtsummary; rendered through `gt::as_gtable()` so no browser is needed): Fig 23 top 10 changed genes with the GRCh38 variant's dbSNP rsID (Ensembl REST overlap, GRCh38, allele-matched, indels matched by inserted/deleted sequence; cached in `data/metadata/rsid-cache.tsv`), positions on both references, effects, change, paired z and P (`prepare-figure-data.py` `top10_table()` → `top10-changed.tsv`); Fig 24 the 6 main-set genes that colocalise on one reference only. New slide 21 with Fig 25, a genome-browser-style RNA coverage plot (`coverage_example()` → `coverage-example.tsv`: `samtools depth -a -Q 255 -r` on the 369 ASE target BAMs, CCZ1 and CCZ1B exons on T2T and patched GRCh38, reads per million uniquely mapped, mean per genotype of the CCZ1 eQTL variant; 37 / 130 / 180 people with 0 / 1 / 2 copies): 2 vs 0 copies, patched GRCh38 CCZ1 0.88×, CCZ1B 1.79×; T2T CCZ1 1.61×, CCZ1B 0.87×. Run interactively (3.8 min, 4.3 GB; a first version with `-b`, which scans the whole BAM, was stopped). Deck 25 + 12 slides. Same day the coverage example was compared across four paralog pairs (2 vs 0 copies of the first gene's eQTL allele, total exon coverage, patched GRCh38 / T2T): CCZ1/CCZ1B 0.88, 1.79 / 1.61, 0.87; ASAH2/ASAH2B 0.47, 1.40 / 4.48, 0.87; PMS2CL/PMS2 1.41, 0.75 / 0.91, 0.94; SMN2/SMN1 0.80, 1.25 / 0.19, 1.11 (only 2 people with 0 copies). ASAH2/ASAH2B is the clearest and is the largest change in the top-10 table, so Fig 25 and Fig 22 now use it (`COVERAGE_PAIR`, `EXPR_EXAMPLES` += ASAH2); Fig 25 redrawn as 0 vs 2 copies (filled areas) with the 2 vs 0 ratio in each panel, x = exons joined 5′ to 3′; Fig 22a no longer draws a −1 line (ASAH2B's absolute gain is larger than ASAH2's loss). Run interactively (5.9 min, 4.4 GB). Top-10 table (Fig 23) extended the same day: GRCh38 REF/ALT; the effect allele (the allele the slopes refer to, i.e. the T2T ALT written on GRCh38 via the genotype correlation: GRCh38 ALT if geno_r > 0, else REF; 9 of 10 are the GRCh38 REF, because CHM13 carries the other base); its frequency from 1000 Genomes phase 3 ALL (Ensembl REST, cached in `data/metadata/af-cache.tsv`), with gnomAD genomes or MAGE where 1000 Genomes has none (GLIDR, RASA4DP); MAGE frequencies kept in `top10-changed.tsv` (they agree except STAG3L1, 0.34 MAGE vs 0.56 1000 Genomes phase 3, a duplicated region); region and T2T position dropped from the table (in the tsv and the slide text) so it fits; canvas sized to the table (190 × 50 mm). Figs 22 and 25 now label genotypes by allele with zygosity in the GRCh38/dbSNP convention (rs9794140: A/A hom ALT, A/C het, C/C hom REF; C is the effect allele). REF and ALT are swapped on T2T at 9 of the 10 top variants, so the labels always name the reference convention and keep the letters. Fig 23 marks genes with * (segmental duplication) and † (reversed: opposite sign of slope on the two references), explained in the first footnote line. Fig 21 redrawn as the rate of eQTL change inside vs outside segmental duplications (59 of 385 = 15% vs 58 of 14,539 = 0.4%, ~38×) instead of the share of changed genes in SDs; slide 20 retitled accordingly. Then redrawn again at the user's request as the 117 changed genes by likely reason (`changed_mechanism()` → `changed-mechanism.tsv`; first that applies: ≥ 98% SD 59; read partner annotated on one reference only 3; 90–98% SD ≥ 50% of the gene 14; read partner r ≤ −0.7 6; read destination run without a clear partner 10; not yet traced 25), plain bars with counts and examples; 82 of 117 have a near copy. Slide 20: 'Why they change: 82 of 117 changed genes have a near copy'. Replaced at the user's request by a segmental-duplication-only view (the mechanism mix used read destination for the core genes only, so it was not consistent): ≥ 98% SD 59 (of 385 compared genes), 90–98% SD covering ≥ 50% of the gene 15 (of 419), unique sequence 43 (of 14,120); written to `changed-by-sd.tsv` by `make-figures.R`. Slide 20: '74 of 117 changed genes lie in segmental duplications'. `changed-mechanism.tsv` is kept but no longer shown. Figs 22 and 25 back to genotype letters only (A/A, A/C, C/C, 'effect allele C'), without hom REF/ALT, because REF and ALT are swapped between the references; Fig 23 footnote says REF/ALT are as in GRCh38 and swapped on T2T at 9 of 10 variants. GWAS table (Fig 24) now gives the shared variant per gene (the top variant by combined eQTL + GWAS evidence on the reference where it colocalises, GRCh38 ID, rsID, GWAS P; `coloc_table()` → `coloc-table.tsv`); the table spans the slide with bullets below. Note SERF1B (GWAS P 7 × 10⁻⁵) and SMG1P6 (3 × 10⁻⁶) are below genome-wide significance at the shared variant. Supplementary Fig 25s (backup slide B13): the ASAH2/ASAH2B coverage with all three genotypes (57 A/A, 140 A/C, 150 C/C of the 347 re-aligned people); heterozygotes lie between the homozygotes (ASAH2 patched 1 : 0.81 : 0.47, T2T 1 : 2.35 : 4.48). Fig 21 simplified again to two bars, segmental duplication (≥ 90% identity covering ≥ half of the gene: 74 of 117) vs unique sequence (43), counts only. GWAS slide (2026-11-04): colocalisation defined as whether the eQTL and GWAS signals rise and fall together across the region's variants (one shared causal variant, PP4, vs two nearby ones, PP3), plus an 'In short' takeaway bullet. `docs/data-journey/data-journey.qmd` (2026-11-04): new subsection 'Swapped alleles, step by step (ASAH2)' in the comparison chapter (each reference mapped with its own dosages; 731-person genotype cross-table, geno_r −0.9987; slope per A = −slope per C; paired test flips the dosage column; the common allele is the T2T ALT); rendered with quarto. | `results/07-presentation/` |
| 8 | Reproducibility set-up (branch `capstone`) | **done 2026-11-17**: site config, shared shell library, submit wrapper, pixi environments `bio`/`docs`/`lint`, tasks, data-sources table and sha256 sums (job 395082, 149 files; verified by 395083, 0 failed; a planted bad and a missing file are caught), smoke test on a compute node (395113: tools from `bio`, Python and R from `default`, TMPDIR removed), every Slurm script dry-submitted through `submit.sh`, Nextflow `-preview` from another folder, lint + CI, licences, citation, README; talk-deck builder removed from the branch (11.5) | branch `capstone` |

## 13. Step 2: five-sample RNA-seq pilot (go/no-go 1)

**Question.** Before spending a week of compute on 451 samples: when the same RNA-seq reads are aligned to T2T and to GRCh38, do the genes we expect to differ actually differ, while ordinary genes stay the same?

### 13.0 Where the pilot data come from

Summary for the pilot; the full list of every resource, with URLs, dates and checks, is in section 4.2.

- **Expression (RNA-seq reads):** the Geuvadis project (Lappalainen et al. 2013, *Nature*): lymphoblastoid cell lines of 1000 Genomes individuals, paired-end 2 × 76 bp, unstranded. Raw reads are archived at ENA as study **PRJEB3366** (667 runs, 464 individuals) and mirrored at NCBI SRA under the same run IDs (`ERR…`). The pilot downloaded the 5 runs from **NCBI SRA** (ENA was too slow from Karakoram, 13.3) and checked that every run has exactly ENA's read-pair count.
- **Which individuals:** the Geuvadis individuals (464 in ENA's run report for PRJEB3366) that also appear in the T2T genotype panel's sample names (3,202): 451; 13 have no T2T genotypes. Made by hand on 2026-10-13; since 2026-10-27 recorded, with MAGE, in `config/samples/expression-individuals.tsv` by `scripts/01-genotypes/make-sample-lists.sh` (same 451).
- **Sample details** (sex, ancestry, sequencing lab): Geuvadis sample sheet E-GEUV-1 from EBI BioStudies (`data/metadata/geuvadis/E-GEUV-1.sdrf.txt`).
- **Pilot sheet:** `config/samples/pilot-5.tsv`, written by `scripts/03-rnaseq-pilot/make-pilot-sheet.py` from the ENA report (URLs, md5s, read pairs) and E-GEUV-1 (sex, ancestry, lab). Selection rule in 13.1.
- **Genotypes:** 1000 Genomes T2T-CHM13v2.0 phased panel (section 4), not used in the pilot itself.

### 13.1 Decisions (2026-10-17)

- **Gene models: NCBI RefSeq RS_2025_08 on both references**, replacing the JHU RefSeq/Liftoff v5.3 set named earlier. RefSeq annotates each assembly with the same pipeline in the same release; Liftoff projects GRCh38 genes onto T2T, which is close to the liftover this study avoids. Genes are matched by NCBI GeneID (55,334 shared; 3,875 GRCh38-only, 3,172 T2T-only). The GRCh38 false copies on chr21p are annotated as separate GRCh38-only genes (e.g. LOC102724560 "cystathionine beta-synthase like", LOC102724594 "U2AF1 like 5"), so read splitting between CBS and its false copy is directly visible.
  - T2T: `chm13v2.0-ncbi/annotations/ncbi/T2T-CHM13v2.0.ncbiRefSeq.gtf` (admin bundle, `NC_…` names)
  - GRCh38: `grch38-noalt/annotations/ncbi/GRCh38_noalt.ncbiRefSeq.gtf` (`chr…` names)
- **Aligner and quantification:** STAR 2.7.11b (module) with ENCODE long-RNA settings, quantified two ways:
  - STAR `GeneCounts`: uniquely mapped read pairs only. Shows where T2T rescues reads that GRCh38 discards or misplaces.
  - salmon 1.11.4 in alignment mode on STAR's transcriptome BAM: multi-mapping reads shared between transcripts by EM.
  
  Comparing both separates "the reference changes where reads go" from "the method handles ambiguous reads differently", the question reviewers will ask for every duplicated gene.
- **Samples:** five single-run samples near the median depth (29 M pairs), mostly from lab 1, three male and two female (to test chrY and PAR handling). Two Yoruba, because African genomes differ most from both references and vary most in paralog copy number.

  | Run | Sample | Ancestry | Sex | Lab | Read pairs |
  |---|---|---|---|---|---:|
  | ERR188298 | NA18867 | Yoruba | female | 1 | 28,283,821 |
  | ERR188139 | NA18917 | Yoruba | male | 1 | 26,971,208 |
  | ERR188107 | NA12342 | Utah (CEU) | male | 1 | 29,022,959 |
  | ERR188242 | HG00258 | British | female | 1 | 29,322,467 |
  | ERR188197 | HG00341 | Finnish | male | 2 | 29,464,398 |

### 13.2 Making the two references comparable

Checks on 2026-10-17 found three differences between the assemblies that are not about T2T completeness and would create false differences in LCL RNA-seq:

| Difference | Why it matters | Fix (in `prepare-star-references.sbatch`) |
|---|---|---|
| GRCh38_noalt has **chrM**; NCBI CHM13v2.0 has none | Mitochondrial transcripts are a large share of RNA-seq reads. On T2T they would land on nuclear mitochondrial insertions (NUMTs) and inflate nearby genes | Append the identical chrM (rCRS, from GRCh38_noalt.fa) to T2T, and its 37 RefSeq genes from the GRCh38 GTF |
| GRCh38_noalt has **chrEBV**; CHM13 has none | LCLs are EBV-transformed and express viral genes; without the EBV genome these reads may mis-map | Append the identical chrEBV to T2T (no genes annotated on either) |
| GRCh38_noalt hard-masks the **chrY PARs** (2.77 Mb + 0.33 Mb of N); CHM13 chrY PARs are unmasked | PAR genes would be unique on GRCh38 but split between X and Y on T2T | Hard-mask T2T chrY PARs (0-based `NC_060948.1` 0–2,458,320 and 62,122,809–62,460,029), as in the UCSC `chm13v2.0_maskedY` analysis set |

Not changed: GRCh38_noalt keeps its 171 unplaced and random contigs, since they are part of GRCh38 as it is used. Result: T2T 24 + 2 contigs, GRCh38 195 contigs.

### 13.3 Scripts (`scripts/03-rnaseq-pilot/`)

| Order | Script | What it does | Request |
|---|---|---|---|
| 1 | `prepare-star-references.sbatch` | Array 0–1 (chm13v2.0, grch38): builds `data/star-references/<ref>/` (genome with the fixes above, GTF, gene table with GeneIDs, transcript→gene map, transcript FASTA via gffread, STAR index with `--sjdbOverhang 75`, md5s) | 16 CPUs, 64 GB, 4 h each |
| 1 (parallel) | `fetch-sra.sbatch` | Gets the 5 runs from NCBI SRA: `prefetch` (retried, resumes after S3 resets) → `vdb-validate` (NCBI checksums) → `fasterq-dump` → pigz; read pairs must equal ENA's count; FASTQs made read-only, `.sra` deleted. Used since 2026-10-17 | 12 CPUs, 16 GB, 6 h |
| (unused) | `fetch-fastq.sbatch` | The ENA route: 10 FASTQs over HTTPS with resume, md5 checked against ENA. Too slow on 2026-10-17 (one connection at ~0.75 MB/s); kept in case ENA speeds up | 2 CPUs, 4 GB, 12 h |
| 2 | `align-quantify.sbatch` | Array 0–9 (%2) = 5 runs × 2 references: STAR alignment + GeneCounts, salmon, BAM index, idxstats | 16 CPUs, 48 GB, 3 h each |
| 3 | `compare-expression.py` | Joins both references by GeneID and writes the comparison (13.5) | `srun`, 2 CPUs, 16 GB, minutes |
| (once) | `make-pilot-sheet.py` | Writes `config/samples/pilot-5.tsv` from the ENA report and E-GEUV-1 (13.0) | seconds |

How to run (from the project root):

```
sbatch scripts/02-references/prepare-star-references.sbatch
sbatch scripts/03-rnaseq-pilot/fetch-sra.sbatch
sbatch --dependency=afterok:<starref job>:<sra job> scripts/03-rnaseq-pilot/align-quantify.sbatch
pixi run python scripts/03-rnaseq-pilot/compare-expression.py
```

Peak use is 2 × 16 cores, well within the 128-core limit, but it shares the node with the `t2t-atlas` ChIP-seq runs.

### 13.4 Supporting files

- `config/samples/pilot-5.tsv`: the five runs with ENA HTTPS URLs and md5s, built from the ENA report and the E-GEUV-1 sample sheet.
- `config/chrom-names-chm13v2.0.tsv`: 24 `chrN` ↔ `NC_…` rows plus `chrM` and `chrEBV`, checked against the T2T FASTA index.
- `pixi.toml` / `pixi.lock`: python 3.11, pandas, matplotlib, gffread 0.12.

### 13.5 Outputs (`results/03-rnaseq-pilot/comparison-geuvadis/`)

- `mapping-stats.tsv`: unique, multi-mapped and unmapped %, and chrM / chrEBV share of mapped records, per run and reference.
- `genes.tsv`: per GeneID, the loci on each reference, the share of the T2T gene covered by segmental duplications (`chm13v2.0_SD.bed`), and counts, CPM and mean log2(T2T/GRCh38) per method.
- `candidates.tsv`: the candidate families:
  - chr21 false duplication: CBS, U2AF1, CRYAA, SIK1, ICOSLG, DNMT3L, GATD3, PWP2, TRAPPC10, each with its GRCh38-only LOC copy, plus KCNE1
  - 1q21: NOTCH2, NOTCH2NLA/B/C/R and all NBPF genes
  - 5q13: SMN1/2, NAIP, SERF1A/B, GTF2H2/2C
  - 16q22: HYDIN/HYDIN2
- `summary.tsv`: expressed genes, genes changed up or down on T2T, and the SD share among changed vs all expressed genes, per method.
- `scatter-star.png`, `scatter-salmon.png`: mean expression on T2T vs GRCh38, with SD genes highlighted.

**SD caveat.** `chm13v2.0_SD.bed` also covers parents of processed pseudogenes (ACTB is 73% covered, GAPDH 49%), so 9,608 of 58,563 T2T genes are ≥ 50% "in SDs". Good enough for the pilot's enrichment check; the full analysis needs a stricter definition (e.g. high-identity SDs only).

**Observations from the first sample (NA18867, 2026-10-17), affecting how results are read:**

- **rRNA.** About 1.75% of reads are "unmapped: other" on T2T vs 0.02% on GRCh38. Of 20,000 sampled, almost all align on GRCh38 to its five rDNA copies (chr21 p-arm, `chrUn_GL000220v1`, `chr22_KI270733v1_random`). T2T has the full rDNA arrays (hundreds of copies), so STAR's seeds hit too many loci and the reads stay unmapped. Real biology, not an error. STAR counts are unaffected (on GRCh38 these reads multimap and get 0 unique counts), but salmon gives them to GRCh38's 45S genes (~10,000 CPM), while on T2T they stay unmapped (T2T annotates only 14 45S-type genes). So `compare-expression.py` uses genes on both references minus rRNA and chrM genes as the CPM library size.
- **Small RNAs missing from salmon.** STAR's transcriptome BAM leaves out 582 of the RefSeq `unassigned_transcript_*` records (NCBI's IDs for mature miRNAs and similar small RNAs, 8–2,866 bp). salmon warns and gives them zero. The count is identical on both references, so it does not bias the comparison.
- **Multi-mapping.** Reads mapped to several loci: 3.3% on T2T vs 6.8% on GRCh38.

### 13.6 Go/no-go 1 criteria

- **Go** if both hold:
  1. The candidates show clear, same-direction differences in all five samples. The expected direction: on GRCh38, reads for CBS, U2AF1 and the other chr21 genes are split with their false copies, so these genes rise on T2T in the unique-read counts. Some 1q21 / 5q13 paralogs should also change.
  2. Single-copy genes agree closely (median |log2 ratio| near 0), and changed genes are enriched in SDs. This shows the comparison is fair.
- **No-go** if the candidates do not differ or the differences are not consistent across samples: return to the `t2t-atlas` k-mer methods paper.
- Also record from the ledger: runtime and storage per sample per reference, to cost step 3.

### 13.7 Results (2026-10-17, five samples)

Outputs in `results/03-rnaseq-pilot/comparison-geuvadis/`. CPM = counts per million, library size as in 13.5. "STAR" = uniquely mapped reads only; "salmon" = ambiguous reads shared by EM.

**Mapping** (`mapping-stats.tsv`):

| | T2T | GRCh38 |
|---|---|---|
| Uniquely mapped | 85.3–92.8% | 84.5–91.8% (T2T higher in all 5, by 0.4–1.0 points) |
| Mapped to multiple loci | 3.3–4.8% | 6.4–10.0% (about half on T2T) |
| Unmapped "other" (mostly rRNA on T2T, 13.5 notes) | 0.4–3.9% | ≤ 0.02% |
| chrM share of mapped records | 4.5–8.3% | 3.6–7.8% |

**Genome-wide** (`summary.tsv`): about 13,100 expressed genes have a GeneID on both references.

| | STAR | salmon |
|---|---|---|
| Median \|log2 T2T/GRCh38\| of expressed shared genes | 0.009 | 0.002 |
| Changed ≥ 2-fold, same direction in all 5 (up / down on T2T) | 87 / 89 | 53 / 74 |
| Share of changed genes ≥ 50% in T2T SDs (vs 4.7–5.6% of all expressed) | 39% | 51% |

Ordinary genes agree almost exactly, so the comparison is fair; the genes that change are 8–9× enriched in segmental duplications.

**chr21 false duplication** (`candidates.tsv`; STAR CPM, mean of 5; per-sample log2 ratio range):

| Gene | T2T | GRCh38 | GRCh38 false copy (salmon) | log2 T2T/GRCh38 per sample (STAR) |
|---|---:|---:|---:|---|
| U2AF1 | 126.9 | 0.3 | 89.9 | +6.4 to +8.1 |
| GATD3 | 123.6 | 5.8 | 69.4 | +3.1 to +6.4 |
| PWP2 | 62.3 | 6.5 | 33.9 | +2.0 to +6.6 |
| ICOSLG | 70.4 | 6.0 | 37.5 | +2.6 to +5.7 |
| SIK1 | 35.4 | 3.3 | 15.2 | +2.6 to +4.3 |
| CBS | 3.6 | 0.4 | 0.3 | +1.5 to +2.9 |
| TRAPPC10 | 98.4 | 50.5 | 21.0 | +0.7 to +1.5 |
| CRYAA, DNMT3L, KCNE1 | not expressed in LCLs (< 0.2 CPM) | | | |

On GRCh38 reads of these genes match the real and the false copy equally: STAR discards them (unique counts collapse) and salmon splits them (much of the expression goes to the false copy). On T2T they are counted once, in the right gene. Same direction in all five samples.

**Other candidate families:** 5q13 shows consistent, moderate gains on T2T in unique counts (SMN1 +1.0 to +1.3, GTF2H2 +1.1 to +1.3, GTF2H2C +1.0 to +1.7 log2), but salmon largely evens these out (reads are shared, not lost). NOTCH2NL and most NBPF genes barely change (NBPF14 +0.7 to +1.5 and NBPF1 −0.4 to −1.7 in STAR are the exceptions). HYDIN/HYDIN2 are not expressed in LCLs.

**Further classes of difference seen** (from the top changed genes and genes annotated on one reference only); each needs its own category in step 6:

1. **Paralog / SD read placement** (expected): e.g. DDT, SMIM11, NUDT4P2, TBC1D3, NPIPB4, SULT1A4, POLR2J3, UPK3BL2 and the chr21 genes.
2. **Reference MHC haplotype:** HLA-DRB3 is annotated only on T2T (392 CPM), HLA-DRB5 only on GRCh38 (264 CPM), because the two references carry different DR haplotypes. Expression of these genes then depends on which haplotype the reference happens to carry, not on the sample. Important for LCL eQTLs (HLA is among the strongest eQTL regions) and needs care (e.g. HLA-aware quantification or excluding these genes).
3. **Annotation differences on the same sequence:** e.g. HCP5 (0 vs 162 CPM in STAR; a T2T-only lncRNA overlaps it, making reads ambiguous), NEAT1 (14 vs 288 CPM in STAR, not in salmon), and **P2RY8**, a B-cell gene in PAR1: T2T RefSeq annotates it only on chrY (masked here) and not on chrX, so it reads 0 on T2T vs 172 CPM on GRCh38. Annotation artefacts, not biology. STAR unique counting is more sensitive to them than salmon.
4. **rRNA** (13.5 notes): excluded from normalisation; not of interest for eQTLs.

**Cost for step 3** (from the ledger): download from NCBI ≈ 3.6 min per run (3 in parallel); alignment + salmon 5–10 min per sample per reference at 16 CPUs (≈ 1.9 CPU-h per sample per reference); peak memory 33 GB (STAR) per task. For 451 samples × 2 references ≈ 1,700 CPU-h, about 14 h at 128 cores. Per sample on disk: FASTQ ~4.2 GB (deleted after use), genome BAM ~4 GB per reference (scratch), results ~9 MB per reference (kept).

**In light of FixItFelix (added 2026-10-17, 2.2):**
- The chr21 result replicates Behera et al. 2023: same Geuvadis data, same mechanism; PWP2 and GATD3 are their eQTL examples.
- What is new in the pilot is the rest: only 13 of the 176 changed genes (STAR) overlap FixItFelix's patched regions.

**Assessment against 13.6:** both conditions are met. The chr21 candidates change strongly and in the same direction in all five samples, ordinary genes agree (median |log2| < 0.01), and changed genes are concentrated in SDs. The 1q21 NOTCH2NL/NBPF candidates do not change much, so H2/H3 effects are concentrated in specific loci, not every SD family. **Recommendation: go**, with three additions to the design before step 3: an annotation-difference category and filter (HCP5, NEAT1, P2RY8), an explicit decision on MHC genes, and a stricter SD definition.

---

### 13.8 MAGE pilot (step 2b, started 2026-10-17)

**Why.** MAGE (Taylor et al. 2024; SRA PRJNA851328) is newer and larger than Geuvadis:
- 731 individuals from all 26 1000 Genomes populations, all of them in the T2T genotype panel;
- one lab, NovaSeq 6000, **2 × 150 bp** reads (Geuvadis: 2 × 76 bp, 7 labs, 5 populations).

Only 93 individuals are in both cohorts (1,089 different individuals together). MAGE itself was analysed on GRCh38 only; "in the T2T panel" means its individuals have T2T genotypes, not that MAGE has T2T expression. Longer reads should place better in paralogs, so the T2T effect may be smaller for some families (e.g. SMN, NBPF) but not for the chr21 false duplication (identical sequence). If the effect holds, MAGE may become the discovery cohort and Geuvadis the replication (to decide after this pilot).

**Samples** (`config/samples/mage-pilot-5.tsv`, from `scripts/03-rnaseq-pilot/make-mage-pilot-sheet.py`): single-run individuals near MAGE's median depth (28.4 M pairs), excluding the known swap SRR19762530.

| Run | Individual | Population | Sex |
|---|---|---|---|
| SRR19762269 | HG03129 | ESN (Africa) | female |
| SRR19762285 | NA19982 | ASW (Africa, admixed) | male |
| SRR19762887 | HG00581 | CHS (East Asia) | female |
| SRR19762232 | HG04107 | STU (South Asia) | male |
| SRR19762408 | HG00265 | GBR (Europe) | male |

HG00265 is also in Geuvadis (ERR188222, lab 2; `config/samples/geuvadis-overlap-1.tsv`), so the same person is measured in both cohorts. None of the five Geuvadis pilot individuals is in MAGE. No AMR individual (five slots).

**Same pipeline, two differences:**
- The STAR indexes were built with `--sjdbOverhang 75` for 76 bp reads. For 150 bp reads the ideal is 149; 75 still works, slightly less sensitive at splice junctions near read ends, and it applies equally to both references. A MAGE-main study should build indexes with 149.
- Strandedness: MAGE is **unstranded** like Geuvadis (STAR forward and reverse gene counts 10.79 M vs 10.73 M in SRR19762269), so the default `--star-column unstranded` is correct.

**Jobs:** downloads 360465 (MAGE) and 360466 (ERR188222); alignments 360467 (array 0–9%2) and 360468 (array 0–1).

**Results (2026-10-17).** Outputs: `results/03-rnaseq-pilot/comparison-mage/` (5 MAGE samples) and `comparison-overlap-geuvadis/` (ERR188222). All runs: read pairs equal ENA's; 8–9 min per sample per reference.

*Mapping.* As in Geuvadis, multi-mapped reads halve on T2T (1.9–2.4% vs 3.8–8.3% on GRCh38), unique mapping is slightly higher on T2T (+0.2 to +0.5 points), and rRNA stays unmapped on T2T (0.8–2.8% "other"). MAGE has more "too short" unmapped reads than Geuvadis (8–16% vs 2–10%) on both references; HG00581 is the weakest sample (79% unique).

*Genome-wide.* About 13,600 expressed genes on both references; median |log2 T2T/GRCh38| 0.014 (STAR), 0.003 (salmon). Changed ≥ 2-fold in all five samples: 215 (STAR; 108 up, 107 down) and 164 (salmon); 42–46% of them in T2T SDs vs 5.3–5.6% of all expressed genes.

*Same genes as in Geuvadis.* Across the two cohorts, the per-gene T2T effect correlates at r = 0.89 (STAR) and 0.91 (salmon); every gene changed in both cohorts changes in the same direction (315 STAR, 299 salmon, counting genes on one reference only too).

*Candidates* (STAR, per-sample log2 T2T/GRCh38, MAGE vs Geuvadis):

| Gene | MAGE | Geuvadis | Note |
|---|---|---|---|
| U2AF1 | +3.7 to +6.3 | +6.4 to +8.1 | still 20× in MAGE; 150 bp reads rescue some reads on GRCh38 (7.6 vs 0.3 CPM) |
| GATD3, PWP2, ICOSLG, SIK1, CBS | +1.4 to +4.8 | +1.5 to +6.6 | same direction in all 10 samples |
| TRAPPC10 | +0.5 to +0.7 | +0.7 to +1.5 | smaller |
| SMN1 | +1.4 to +2.0 | +1.0 to +1.3 | larger in MAGE |
| GTF2H2 / GTF2H2C | +0.5 to +1.2 | +1.0 to +1.7 | |
| NBPF14 / NBPF1 | +0.7 to +1.6 / −1.0 to −2.3 | similar | |
| NOTCH2NLA/B, SMN2 | small, mixed | small, mixed | |
| P2RY8, HCP5, HLA-DRB3, NEAT1 | same as Geuvadis | | annotation / MHC-haplotype artefacts (13.7) recur, as expected |

*Same person in both cohorts* (HG00265; STAR CPM, T2T vs GRCh38): PWP2 50 vs 0.5 (Geuvadis) and 48 vs 2.2 (MAGE); U2AF1 115 vs 0.8 and 132 vs 10; GATD3 136 vs 1.3 and 121 vs 3.9. Whole-transcriptome agreement between the two cohorts for this person: r = 0.96 (log CPM) on either reference.

*Region check (2026-10-17):* 12 of 215 changed genes (STAR) and 6 of 164 (salmon) overlap FixItFelix's patched regions; 94–96% are outside (2.2).

*Assessment.* The T2T effect replicates in a newer, independent cohort with twice-as-long reads: same genes, same direction, similar SD enrichment. Longer reads shrink the chr21 effect (GRCh38 recovers some reads) but do not remove it. MAGE offers more individuals (731, all with T2T genotypes), 26 populations and one lab, so it is the stronger discovery cohort; Geuvadis then serves as replication (decision with the user).

## 13A. Steps 3–4: full expression (MAGE, Geuvadis) and genotypes

Decided 2026-10-17 after the two pilots: **MAGE is the discovery cohort, Geuvadis the replication cohort.** Both are quantified on **three references** from the same reads:
- `chm13v2.0` (T2T);
- `grch38` (no-alt analysis set);
- `grch38-ffx`: FixItFelix's patched GRCh38 (step 3b, decided 2026-10-17).

### 13A.1 Pipeline layout and samplesheets

```
scripts/04-expression/
├── make-cohort-sheets.py          writes the samplesheets (below)
├── fetch-sra-batch.sbatch         part 1: download one batch of .sra from NCBI SRA (prefetch, vdb-validate)
├── run-batch.sbatch               part 2: Nextflow head job for one batch, then delete .sra of completed runs
├── submit-cohort.sh               chains part 1 and 2 for every batch of a cohort
└── nextflow/
    ├── main.nf                    entry point
    ├── nextflow.config            params (input, batch, references, ref_dir, outdir, star_prefix); profiles slurm, test
    ├── conf/base.config           resources and modules for QUANTIFY_RUN
    ├── workflows/expression.nf    samplesheet checks, skip completed runs, QUANTIFY_RUN, per-batch mapping summary
    └── modules/local/quantify_run/main.nf   one task per run: fasterq-dump, then per reference STAR → salmon → idxstats → region BAM
```

The reference scripts used here (`make-region-beds.py`, `build-star-index.sbatch`, `prepare-ffx-reference.sbatch`) are in `scripts/02-references/` since 2026-10-27 (11.4).

**Samplesheets** (`config/samplesheets/`, written by `make-cohort-sheets.py`, one row per sequencing run):

| File | Runs | Individuals | Batches (60 runs) |
|---|---:|---:|---:|
| `mage.csv` (discovery) | 779 | 731 | 13 |
| `geuvadis.csv` (replication; individuals with T2T genotypes) | 652 | 451 | 11 |
| `test-mage-2.csv` (pipeline test, profile `test`) | 2 | 2 | 1 |

Columns:
- `sample`: 1000 Genomes individual;
- `run`: SRA/ENA run;
- `cohort`, `batch`;
- `read_pairs`: ENA's count, checked exactly after `fasterq-dump`;
- `read_length`: picks the STAR index `star-2.7.11b-oh<read_length − 1>`;
- `sex`, `ancestry`, `lab`;
- `sra`: where part 1 puts the validated `.sra`.

Notes on the samplesheet contents:
- Individuals with several runs (MAGE 24 with 3; Geuvadis 173 with 2 or 8) are quantified per run and merged afterwards.
- Two Geuvadis individuals have no sex in E-GEUV-1.
- MAGE's reported sample swap (SRR19762530) is already corrected in ENA: the run is withdrawn, and SRR37907958 = NA11919, SRR37907959 = HG00237.

### 13A.2 How a run works and how to resume

1. **Download** (`fetch-sra-batch.sbatch <cohort> <batch>`):
   - fetches one batch of 60 runs as `.sra` (~160 GB);
   - `prefetch` is retried and resumes after S3 resets; then `vdb-validate`, then `<run>.validated`.
   - Batches exist because scratch had 2.8 TB free (shared) and all MAGE `.sra` would be ~2 TB.
2. **Nextflow** (`run-batch.sbatch <cohort> <batch>` → `nextflow run scripts/04-expression/nextflow/main.nf -profile slurm --input config/samplesheets/<cohort>.csv --batch <k> -resume`):
   - one task per run; the FASTQ exists only inside the task, on scratch;
   - `fasterq-dump`, and the read pairs must equal ENA's count;
   - for each of the three references: STAR (pilot settings) → salmon → idxstats → keep alignments in `keep-regions.bed` → delete the genome BAM;
   - `complete.txt` last;
   - the head job then deletes the `.sra` of completed runs.
3. **Whole cohort:** `bash scripts/04-expression/submit-cohort.sh <cohort>` chains all batches. At most about two batches of `.sra` are on disk, and one Nextflow head runs at a time. The job IDs go to `logs/04-expression/submitted-<cohort>-<time>.tsv`.
4. **Resume:** rerun `sbatch scripts/04-expression/run-batch.sbatch <cohort> <k>`.
   - Each batch has its own launch folder (`work/04-expression/launch/<cohort>-b<k>/`) and shares the work folder (`work/04-expression/nf-work/<cohort>/`).
   - Runs with `complete.txt` are skipped by both parts, and validated `.sra` files are not fetched again.
5. **Test:** `sbatch scripts/04-expression/run-batch.sbatch mage 1 test` (profile `slurm,test`, the 2 runs of `test-mage-2.csv`).

**Resources:**
- QUANTIFY_RUN: 16 CPUs, 48 GB, ≤ 6 h; about 30 min per run for three references.
- `queueSize` 6 → 96 cores, 288 GB, within the 128-core / 400 GB limit.
- The head job runs on `long` (2 CPUs, 8 GB). Traces go to `logs/04-expression/trace-<cohort>-b<k>-<job>.tsv`.

**Outputs** (`results/04-expression/<cohort>/<run>/`):
- `fastq-check.tsv`, `log.txt`, `mapping-summary.tsv`, `complete.txt`;
- per reference: `Log.final.out`, `Log.out`, `ReadsPerGene.out.tab`, `SJ.out.tab`, `idxstats.tsv`, `salmon/`, `regions.bam(.bai)`, `time-star.txt`, `time-salmon.txt`;
- per batch: `results/04-expression/summaries/mapping-summary-batch<k>.tsv`.

Genome BAMs are not kept (~8 GB per sample per reference); only the region BAMs (~13 Mb of loci per reference).

### 13A.3 The third reference: patched GRCh38 (`grch38-ffx`)

- **Source:** FixItFelix's `GRCh38_masked_v2_decoy_gene.fasta.gz` (Zenodo 7532965; downloaded and md5-checked by `prepare-ffx-reference.sbatch`, stored in `data/external/fixitfelix-behera2023/`).
- **What it is:** the GRCh38 no-alt analysis set with the false copies masked to N and decoy contigs added for collapsed genes.
- **Gene models:** the same GRCh38 RefSeq GTF, so masked genes receive no reads and decoys are not counted. Genes and transcripts keep GRCh38 coordinates and GeneIDs.
- **Checks** (written to `compare-to-grch38.tsv`): every GRCh38 contig present with the same length; decoys listed; newly masked bases per contig.
- **STAR indexes:** `star-2.7.11b-oh75` and `-oh149`.
- **Genotypes:** GRCh38, the same as for `grch38`.

### 13A.4 Genotypes (steps 1b and 4)

- **T2T** (`subset-genotypes.sbatch`, now its default): the 1,089 expression individuals → `results/01-genotypes/chm13v2.0/expression/` (job 360598, done 2026-10-17; folder renamed 2026-10-27 from `01b-…/geuvadis/`, see 16).
- **GRCh38 download** (`scripts/01-genotypes/grch38/fetch-panel.sbatch`): NYGC phased panel → `data/external/1kgp-grch38-nygc-phased-20220422/` (30 GB).
  - Job 360594 verified 44 of 46 files. chr3 and chr7 had the right size but the wrong md5: bytes corrupted in transit, which plain HTTP cannot detect but the manifest md5 caught.
  - Both files were deleted and re-downloaded (job 360849); the 44 verified files are skipped.
  - All 1,089 expression individuals and NA12878 are in the panel (3,202 samples).
- **GRCh38 subset** (`scripts/01-genotypes/grch38/subset-genotypes.sbatch`, job 360850, after 360849) → `results/01-genotypes/grch38/`.
  - Same filters as the T2T subset: individuals, `norm -m-`, AC/AN/MAF recomputed, monomorphic sites dropped.
  - Plus SNVs and indels only, because the NYGC panel also has SVs and the T2T panel does not.
  - Next: match variants between references for step 6.

### 13A.4b Merging and gene classes (written 2026-10-17, tested on the 2 completed runs)

- **`scripts/04-expression/merge-expression.py --cohort <cohort>`** → `results/04-expression/<cohort>/matrices/`:
  - `runs.tsv`: every samplesheet run, whether it is included, and unique, multi and "unmapped other" % per reference;
  - per reference: `star-counts.tsv.gz` (STAR unstranded, uniquely mapped pairs), `salmon-counts.tsv.gz` (NumReads summed to genes) and `salmon-tpm.tsv.gz`;
  - genes × individuals: counts summed over an individual's runs, TPM averaged;
  - rows are that reference's `gene_id`; match references through `gene_ncbi` in its `genes.tsv`.
  - Rerun whenever more runs are complete.
- **`scripts/06-compare/make-gene-classes.py`** → `results/06-compare/gene-classes.tsv`: one row per GeneID.
  - Fields: `annotated_on` (both / t2t_only / grch38_only); `region_class` (first that applies: `ffx_patched` > `t2t_new_sequence` (≥ 50% in `chm13v2-unique_to_hg38.bed`) > `sd` (≥ 50% in T2T SDs) > `other`); and flags for MHC, PAR, chrM and rRNA.
  - Spot checks: U2AF1, PWP2 and LOC102724594 → ffx_patched; SMN1, NOTCH2NLA → sd; HLA-DRB1 → sd + MHC; CD99 → PAR.

  | Region class | All genes | Protein-coding |
  |---|---:|---:|
  | ffx_patched | 156 | 56 |
  | t2t_new_sequence | 2,798 | 351 |
  | sd | 7,320 | 1,222 |
  | other | 52,107 | 18,714 |

  - **Stricter SD definition (done 2026-10-17):** `sd` now uses only high-identity SDs (≥ 98%: UCSC colours yellow and orange in `chm13v2.0_SD.bed`; grey = 90–98%).
    - ACTB and GAPDH, which were covered only by grey SDs, move to `other`. SMN1, NOTCH2NLA, NBPF10 and HYDIN2 stay `sd`.
    - The any-identity fraction is kept as `t2t_sd_frac`.
    - Counts now: ffx_patched 156, t2t_new_sequence 2,798, sd 3,050 (636 protein-coding), other 56,377.

**Pilot changed genes by region class** (genes on both references, ≥ 2-fold in all 5 samples):

| Pilot / method | Changed | ffx_patched | t2t_new_sequence | sd (≥ 98%) | other | Expressed background: sd / new / ffx |
|---|---:|---:|---:|---:|---:|---|
| Geuvadis / STAR | 176 | 13 | 13 | 44 | 106 | 2.1% / 0.2% / 0.2% |
| Geuvadis / salmon | 127 | 5 | 11 | 46 | 65 | 3.0% / 0.3% / 0.2% |
| MAGE / STAR | 215 | 12 | 10 | 64 | 129 | 2.6% / 0.2% / 0.2% |
| MAGE / salmon | 164 | 6 | 11 | 54 | 93 | 3.0% / 0.3% / 0.2% |

Enrichment over background: high-identity SDs ~10×, T2T-only sequence ~30×, FixItFelix regions ~20–35×.

**The `other` class.** For MAGE / STAR, 129 changed genes are `other`:
- 15 lie in lower-identity SDs and 1 is in the MHC;
- **57 change in STAR but not in salmon** (|log2| < 0.5), e.g. S100A6 (27 vs 2.6 CPM in STAR, flat in salmon), HCP5, NDUFA2, INAFM1. That is the signature of annotation-overlap artefacts: on one reference a different or new overlapping gene makes STAR call the reads ambiguous, while salmon still assigns them;
- 59 change in both methods.

**Annotation-overlap check (`scripts/06-compare/compare-annotation-overlaps.sh` → `results/06-compare/annotation-overlaps.tsv`, 2026-10-17).** For each gene on both references, the script lists the genes whose exons overlap its exons, on each reference (bedtools on the two GTFs).
- 2,026 of 55,371 genes have different overlap partners (906 gain a partner on T2T, 1,206 lose one).
- This explains much of the STAR-only signal. Among pilot STAR-only changes, 41% (Geuvadis) and 56% (MAGE) have different partners, against 7.3–7.5% of expressed genes; among changes where STAR and salmon agree, 19–27%.
- Examples:
  - HCP5: a new overlapping lncRNA on T2T, LOC124906562;
  - S100A6: an overlapping gene only on GRCh38, LOC124904423;
  - INAFM1: overlaps CCDC9 only on GRCh38;
  - P2RY8: overlaps ASMTL on T2T chrY;
  - NEAT1: gains FRMD8 on T2T.
- PWP2 gains GATD3 as a partner on T2T, which is real biology at the chr21 locus.
- NDUFA2 changes without a partner difference, so other causes exist.

**Rule for step 6:**
- A gene is *changed* only when STAR and salmon agree in direction and size.
- STAR-only changes are reported as annotation ambiguity.
- The robust "beyond FixItFelix" set is: changed, outside the FixItFelix regions, and without an overlap difference (or shown separately).

### 13A.5 Status (2026-10-17)

| Job | What | State |
|---|---|---|
| 360593 | STAR indexes oh149 (T2T, GRCh38) | done (21–27 min, ≤ 51 GB) |
| 360623 | Patched GRCh38 reference + oh75/oh149 indexes | done 14:59 (44 min; 63 GB on disk) |
| 360596, 360624 | 2-run test (SRR19762165, SRR19762166) through the new pipeline on all three references | **passed** 15:24 |
| 360689, 360690 | **Full MAGE run, batch 1** (fetch, Nextflow) | done: 60/60 runs |
| 360821–360844 | MAGE batches 2–13, chain of 2026-10-17 17:16 | batch 2 38/60 done, then paused 2026-10-18 (360822 cancelled; 360824–360844 held, then cancelled at resume) |
| 362132–362145 | **MAGE batches 2–8, resumed 2026-10-22 06:02** (`logs/04-expression/submitted-mage-20261022T060200Z.tsv`) | batches 2–6 done; 7 running (27/60 at 2026-10-23 07:37); 8 downloaded |
| 365114, 365115 | MAGE batch 7 and 8 reruns after the disk-full failure (2026-10-23 09:20) | done |
| 362981–362990 | MAGE batches 9–13 (**done**; MAGE complete 2026-10-24), resubmitted 2026-10-23 07:39 with the free-space guard (the 2026-10-22 jobs 362146–362155 were cancelled); 362982 and 362983 re-pointed to wait for 365115 | queued |
| 360594, 360849 | GRCh38 genotype download | 44/46 verified; chr3 and chr7 re-downloading (md5 mismatch) |
| 360850 | GRCh38 genotype subset (1,089 individuals) | waits for 360849 |
| 360850, 360926 | GRCh38 subset; variant matching (MAF rule) | done 18:17, 18:35 |
| 360598 | T2T genotypes for 1,089 individuals | done |

**Patched-reference check** (`data/star-references/grch38-ffx/compare-to-grch38.tsv`):
- all 195 GRCh38 contigs are kept with identical lengths;
- 26 decoy contigs are added (MAP2K3, FANCD2, GPRIN2, DUSP22, KCNJ18 and 21 KMT2C paralogs);
- 1,237,765 bases are newly masked on chr9, chr13, chr16, chr21, chrX and two unlocalised contigs.

**Test results.**
- Both runs: read pairs equal ENA's; all files present for all three references.
- Region BAMs are 42–66 MB per reference.
- 22–24 min and 33.6 GB peak per run (3 references): 11.7 CPU h for 2 runs.
- First three-way look (STAR CPM, two MAGE individuals):

| Gene | T2T | GRCh38 | Patched GRCh38 |
|---|---:|---:|---:|
| U2AF1 | 154–185 | 7–9 | 155–187 |
| SIK1 | 102–146 | 4–10 | 103–147 |
| ICOSLG | 177–209 | 20–30 | 178–211 |
| GATD3 | 108–119 | 22–25 | 108–119 |
| ACTB (control) | 9,880–10,822 | 9,969–10,919 | 9,964–10,914 |

Multi-mapped reads: T2T 2.0–2.1%, GRCh38 3.4–3.9%, patched GRCh38 3.3–3.9%.

So the patch fixes the chr21 genes completely (it matches T2T there) but barely changes multi-mapping genome-wide, while T2T halves it. That is H0 in miniature, to be tested genome-wide.

**Paused 2026-10-18 05:24 UTC and resumed 2026-10-22 06:02 UTC (user requests).**
- *State at the pause:* batch 2's Nextflow job (360822) was cancelled; batches 3–13 (360824–360844) were held (`JobHeldUser`).
  - Completed: batch 1 60/60 runs (Nextflow 360690: 58 runs in 5 h 12 min, plus the 2 test runs), batch 2 38/60, **98 of 779 MAGE runs**.
  - Batch 3 was downloaded (fetch 360823) but not processed. Download times: batch 1 5 h 10 min (stuck-prefetch incident), batch 2 2 h, batch 3 2 h 20 min.
  - 120 validated `.sra` (338 GB) were kept on scratch.
- *Resume:*
  - The held chain was cancelled and a new chain submitted from batch 2: `bash scripts/04-expression/submit-cohort.sh mage 2` → jobs 362132–362155 (`logs/04-expression/submitted-mage-20261022T060200Z.tsv`).
  - Simply releasing the held jobs would have skipped batch 2's remaining 22 runs, whose Nextflow job was cancelled.
  - Batch 2's fetch (362132) found all 120 `.sra` already validated and finished in 30 s; its Nextflow job (362133) resumed the 22 remaining runs.
  - The node was idle at resume (the t2t-atlas ChIP-seq runs had finished); scratch was 85% full (1.1 TB free).

**Progress and scratch incident, 2026-10-23 07:37 UTC.**
- *Progress:* **387 of 779 MAGE runs complete** (batches 1–6 60/60, batch 7 27/60); batch 8 already downloaded.
- *Disk:* the shared scratch (`fastpool`) reached **99% (129 GB free)**.
  - Main user: t2t-atlas, 4.1 TB under `/scratch/yusuf.maulana/t2t-atlas`, with its ChIP-seq pipelines running again.
  - This project: 334 GB of `.sra` (batches 7–8), 121 GB pilot files, 86 GB Nextflow work, 81 GB task temp.
- *Actions (user-approved deletions):*
  - downloads for batches 9–13 put on hold;
  - deleted `work/03-rnaseq-pilot/` (121 GB: pilot FASTQs and full BAMs; all pilot results stay in `results/03-rnaseq-pilot/`);
  - deleted the Nextflow task folders of the 360 completed batch 1–6 runs (366 folders, 96.8 GB; outputs already published; completed runs are skipped on reruns).
  - Scratch then had 288 GB free (ZFS frees space with a delay).
- *Safeguard:* `fetch-sra-batch.sbatch` now waits before each download while scratch has less than `MIN_FREE_GB` (default 400 GB) free.
- *Requeue:* batches 9–13 were cancelled and resubmitted with it: `PREV_FETCH=362144 PREV_RUN=362145 bash scripts/04-expression/submit-cohort.sh mage 9` → jobs 362981–362990 (`logs/04-expression/submitted-mage-20261023T073956Z.tsv`).
- *Jobs:* batch 7 (Nextflow 362143) and batch 8 (fetch 362144 done, Nextflow 362145) continue from the 2026-10-22 chain.

**MAGE expression complete (by 2026-10-24 22:20 UTC): 779 of 779 runs, 731 individuals, three references.**
- Totals (from the Nextflow traces, in the ledger as one row per batch): 777 runs plus the 2 test runs, **6,157 CPU h** (7.9 CPU h per run for three references; above the 6 CPU h estimated from the test).
- Merged with `merge-expression.py --cohort mage`: 58,600 (T2T) / 59,257 (GRCh38, patched GRCh38) genes × 731 individuals.

**Disk full again, 2026-10-23 ~08:30–09:15 UTC; batches 7–8 rerun.**
- *What happened:* scratch filled to 100% (44 GB free). Our usage stayed flat; `t2t-atlas/work/gm12878-tf-pe100` grew from 4.1 to 4.4 TB.
- *Consequences:*
  - batch 7 (Nextflow 362143, then at 32/60) and batch 8 (362145) ended with tasks "terminated by the external system": a task cannot write its exit status on a full disk;
  - the t2t-atlas TF ChIP-seq runs on both references also ended with `exit_status 1` (head jobs 362434, 362435);
  - the new free-space guard held batch 9's download ("waiting: 44G free < 400G"), as intended.
- *Recovery:* by 09:20, 555 GB were free (cleaned outside this project, per the user).
  - Batch 7 rerun as 365114 and batch 8 as 365115 (after 365114); completed runs are skipped and their downloads were still on disk.
  - Batch 9's processing (362982) and batch 10's download (362983) were re-pointed to wait for 365115 (`scontrol update … Dependency=afterany:362981,afterany:365115`), so downloads and parallel processing don't stack up.
- **392 of 779 MAGE runs complete at 09:15.**

**Incident, 2026-10-17 15:51–17:15: stuck downloads.**
- *What happened:* in batch 1, 15 runs downloaded in 25 min. Then three `prefetch` processes (SRR19762180, SRR19762183, SRR19762184) kept restarting their own partial files at ~120 MB for 80 min. They never failed, so the retry loop never ran.
- *Network ruled out:* a fresh `prefetch` at the same time got ~11 MB/s.
- *Fixes:*
  - the three processes were killed and their partial files removed (they then restarted and progressed);
  - `fetch-sra-batch.sbatch` now limits each `prefetch` attempt to 30 min and, after a timeout, removes the partial files and starts that run fresh;
  - all chain dependencies are now `afterany`, so one failed download no longer blocks later batches;
  - Nextflow skips runs without a validated `.sra`, with a warning, instead of stopping the batch.
- *Requeue:* batches 2–13 were resubmitted with the fixed scripts, attached behind batch 1 (`PREV_FETCH`/`PREV_RUN` in `submit-cohort.sh`).
- *Batch 1 safeguard:* batch 1's fetch still runs the old script, so a watchdog kills any of its `prefetch` processes older than 40 min until it ends.

**MAGE cost** (from the test): ~6 CPU h per run for three references → ~4,700 CPU h for 779 runs. At 6 tasks at a time that is about 2 days if the node has room. Downloads (~4 min per run, 3 in parallel) keep ahead of processing.

The t2t-atlas ChIP-seq pipelines are left to finish (user, 2026-10-17); they share the node.

## 13B. Step 5–6 setup: variant matching and eQTL mapping (drafted 2026-10-17)

### 13B.1 Variant matching between references (`scripts/06-compare/match-variants.sbatch`, `match-variants.py`)

**Method.**
- Inputs: common variants (MAF ≥ 0.01 in the 1,089 individuals) of the T2T and GRCh38 genotype subsets.
- Each T2T variant is lifted to GRCh38 (UCSC `hs1ToHg38`, `liftOver -minMatch=0.95`) and paired on position and alleles: same REF/ALT, or swapped where the references carry different bases; on minus-strand chains, SNP alleles are reverse-complemented.
- Indels left unpaired are paired with the nearest unpaired GRCh38 indel within 20 bp that has the same length change and MAF within 0.01.
- Liftover is used only to pair variants called independently on each reference.

**Categories** (T2T side): `matched` (MAF within 0.02), `matched_maf_discordant`, `no_liftover` (T2T-only sequence), `lifted_not_called_grch38`; GRCh38 side: `not_called_t2t`.

**MAF-agreement rule (added 2026-10-17).** A true pair is the same variant in the same individuals from the same reads, so its MAFs agree; exact SNP matches agree to 6 digits. Genome-wide, 0.24% of SNP pairs and 5% of indel pairs differ by more than 0.02. Example: T2T chr22:9,664,292 G>A (MAF 0.056) lifted onto GRCh38 chr22:12,693,176 G>A (MAF 0.495), a paralogous position. Such pairs are `matched_maf_discordant` and are kept for the H4 analysis.

**Output:** `results/06-compare/variant-matching/<chr>.tsv.gz`, `summary.tsv`. Job 360926 (2026-10-17, replaces 360887, which lacked the MAF rule).

**Genome-wide result (MAF ≥ 0.01, 1,089 individuals):**

| | SNPs MAF 0.01–0.05 | SNPs MAF ≥ 0.05 | Indels 0.01–0.05 | Indels ≥ 0.05 |
|---|---:|---:|---:|---:|
| T2T matched | 4,706,353 | 5,970,193 | 910,796 | 1,051,377 |
| T2T matched, MAF discordant | 496 | 25,367 | 3,312 | 99,186 |
| T2T lifted, not called on GRCh38 | 63,734 | 57,910 | 415,603 | 562,478 |
| T2T no liftover (T2T-only sequence) | 17,351 | 18,159 | 28,795 | 25,896 |
| GRCh38 not called on T2T | 1,153,013 | 1,184,064 | 682,545 | 728,473 |

**chr22 test (2026-10-17, MAF ≥ 0.01):**

| | SNPs | Indels |
|---|---:|---:|
| T2T matched | 141,515 (132,571 same alleles, 28,686 swapped across both types) | 30,229 (incl. 10,487 nearby) |
| T2T lifted, not called on GRCh38 | 3,520 | 16,957 |
| T2T no liftover (T2T-only sequence) | 991 | 2,964 |
| GRCh38 not called on T2T | 47,575 | 20,950 |

- Matched SNPs agree in MAF at r = 0.998 (exact matches agree to 6 digits: same reads, same individuals); nearby-indel pairs agree at r = 0.9999.
- Swapped alleles (~20% of matched SNPs) are expected at common sites, because GRCh38 and CHM13 are different haplotypes.

**Finding: the GRCh38 panel calls many common SNPs the T2T panel lacks.**
- On chr22, 47,575 GRCh38 SNPs (median MAF 0.04) have no T2T partner.
- 46,424 of them lift to T2T, and the full T2T panel has a variant at only ~3,600 of those positions.
- They are spread along the whole chromosome, not clustered. 59% of their T2T positions are repeat-masked, against 32% for matched SNPs.
- Example: GRCh38 chr22:29,538,848 C>T (MAF 0.23) → T2T chr22:30,002,217. The sequence is identical, the neighbouring SNPs match perfectly (same lift offset, identical MAF), but there is no T2T call.
- **Excess-heterozygosity test (chr22, 2026-10-17).** Paralogous sequence variants from collapsed duplications make nearly everyone heterozygous, so they show extreme excess heterozygosity (`bcftools +fill-tags -t ExcHet,HWE`, 1,089 individuals):

  | chr22 SNPs, MAF ≥ 0.01 | n | ExcHet p < 1e-6 |
  |---|---:|---:|
  | matched | 141,515 | 0.01% |
  | all T2T | 147,932 | 0.02% |
  | GRCh38, not called on T2T | 47,575 | **4.7%** (MAF 0.1–0.2: 19.6%; 0.2–0.5: 11.0%) |

  So about 2,250 of them (on chr22 alone) look like **GRCh38 paralog artefacts (H4b)**. The other ~95% have normal heterozygosity.
- **Are they in the raw T2T calls?** The panel was built from the unphased Rhie et al. 2023 calls (`all_samples_3202/1KGP.CHM13v2.0.<chr>.recalibrated.snp_indel.pass.vcf.gz` on the human-pangenomics S3; chr22 alone is 33.7 GB).
  - Method: 200 unmatched chr22 GRCh38 SNPs (60 excess-het, 140 normal; stratified random) were probed at their lifted T2T positions, read remotely through the tabix index (HTTP range requests with retries). 157 resolved, 43 failed on S3 resets.
  - A position counts as "same variant" if the raw call has the same alleles and a MAF within 0.1 (raw calls cover 3,202 individuals).
  - Results in `results/06-compare/raw-t2t-probe/chr22-probe-200.tsv`.

  | Unmatched GRCh38 SNPs | Resolved | Same variant in raw T2T calls (PASS) | No raw T2T call | Different variant at the site |
  |---|---:|---:|---:|---:|
  | normal heterozygosity (~95%) | 110 | **56%** | 26% | 17% |
  | excess heterozygosity (~5%) | 47 | 13% | **74%** | 13% |

  (The first 25 positions resolved were all near the chr22 p-arm and mostly lacked raw calls; that early reading was not representative.)
- **Interpretation:**
  - Most unmatched GRCh38 SNPs are real variants that **were called on T2T but dropped when the phased panel was built** (Lalli et al.'s QC and phasing). That is a pipeline difference, not a reference effect.
  - The excess-heterozygosity subset mostly **has no T2T call at all**, as expected for GRCh38 paralog artefacts (**H4b**).
  - The rest have no raw call or a different variant at the site, and remain unexplained.
  - The main comparison stays on matched variants.
  - Using the raw T2T calls instead of the panel would recover the filtered variants but needs the per-chromosome 3,202-sample VCFs (~1.5–2 TB over S3); not planned unless needed.

**Consequence for step 6:**
- The **primary** eQTL comparison uses **matched variants only** on both references.
- Unmatched variants are analysed separately: T2T-only sequence, T2T calls without a GRCh38 call, and GRCh38 calls without a T2T call.
- The latter mostly reflect panel filtering and must not be counted as GRCh38 errors without further evidence (e.g. mapping quality, or collapsed-SD location for H4b).
- The T2T panel's filters (Lalli et al.; README: "filtered for variant call quality") are to be checked against this pattern.

### 13B.2 eQTL mapping (step 5) design

Per cohort (MAGE first) and per reference (`chm13v2.0`, `grch38`, `grch38-ffx`):

- **Expression:** `merge-expression.py` matrices.
  - One column per individual (runs merged).
  - Gene filter as in GTEx: TPM ≥ 0.1 and ≥ 6 reads in ≥ 20% of individuals, applied per reference; record which genes pass on which reference.
  - TMM-normalised CPM, then inverse-normal transform per gene (`qtl.norm`, as in the GTEx pipeline).
  - Phenotype BED at the gene TSS on that reference.
  - Main counts: STAR unique; salmon as the sensitivity analysis.
- **Genotypes:** T2T subset for `chm13v2.0`, GRCh38 subset for `grch38` and `grch38-ffx`. Converted to plink2 pgen; in-cohort MAF ≥ 0.01. Main analysis on matched variants (13B.1); secondary on all variants. Chromosome names must match the phenotype BED (T2T genes `NC_…` → `chr…` with `config/chrom-names-chm13v2.0.tsv`).
- **Covariates:**
  - genotype PCs computed **once** from LD-pruned, matched, common SNPs and used for all three references, so they cannot differ;
  - sex;
  - expression PCs per reference, with the number chosen where eGene discovery saturates (as in FixItFelix);
  - the same covariate set for all references otherwise.
- **Mapping:** tensorQTL `cis` with permutations (beta approximation), ±1 Mb, q-values (Storey) at FDR 5%; nominal statistics for all pairs; conditionally independent eQTLs; SuSiE fine-mapping (tensorQTL `susie`).
- **Comparison (step 6):**
  - genes matched by GeneID, variants by the matching table;
  - effect sizes compared jointly (π1, mashr) rather than by significance alone;
  - each difference labelled with region class (`gene-classes.tsv`), annotation-overlap status (`annotation-overlaps.tsv`), STAR/salmon agreement and variant category;
  - the three-way comparison separates what the FixItFelix patch fixes from what only T2T fixes.
- **Benchmark:** FixItFelix Table S17 (their Geuvadis eQTLs in patched regions) against our Geuvadis replication run.

**Scripts (`scripts/05-eqtl/`, written and tested 2026-10-17):**

| Script | Does |
|---|---|
| `prepare-phenotypes.py --cohort --ref [--method star\|salmon]` | Matrices → genes on chr1–22/X (T2T names mapped to chrN; rRNA dropped) → expressed (TPM ≥ 0.1 and ≥ 6 reads in ≥ 20%) → TMM CPM → inverse-normal → TSS BED. Output: `results/05-eqtl/<cohort>/<ref>/<method>/phenotypes.bed.gz`, `genes-tested.tsv` |
| `prepare-genotypes.sbatch <cohort> <t2t\|grch38> [chroms]` | Genotype subset → plink2 pgen for the cohort's individuals, MAF ≥ 0.01, IDs `chrN:pos:ref:alt`. Two sets: `<geno>.all` and `<geno>.matched` (variant-matching table). Output: `results/05-eqtl/<cohort>/genotypes/` |
| `make-covariates.py --cohort --ref [--method] [--n-expr-pcs auto]` | 5 genotype PCs, computed once per cohort from LD-pruned common matched T2T SNPs and shared by all references; sex; expression PCs per reference (GTEx sample-size rule by default) |
| `run-tensorqtl.py` / `run-tensorqtl.sbatch <cohort> <ref> [method] [variants] [nperm]` | One task per chromosome (array 0–22, 6 at a time): `map_cis` (permutations, beta approximation) and `map_nominal`. Output: `by-chrom/`, `nominal/` |
| `collect-eqtl.py --cohort --ref [--method] [--variants]` | Merges chromosomes and computes q-values (R `qvalue`) → `cis-eqtl.tsv.gz`, number of eGenes at FDR 5% |
| `make-simulated-cohort.py` | Test cohort `simtest` (below) |

**End-to-end test on a simulated cohort (`simtest`).**
- Setup: 300 real MAGE individuals with their real T2T chr22 genotypes, and 300 real chr22 genes. Expression was simulated (negative binomial counts with a hidden factor), with 60 planted eQTLs on real common SNPs within 200 kb of the TSS.
- Run: all scripts on chr22, 1,000 permutations. The run took under 1 min (permutations 14 s).

| Check | Result |
|---|---|
| Planted eQTLs found (FDR 5%) | 43 of 60 (missed ones have the weaker effects) |
| False eGenes among 231 null genes | 2 (observed FDR 2/45 = 4.4%) |
| Null genes with beta-approximated p < 0.05 | 5.6% (calibrated) |
| Effect direction | correct in 98% |
| Lead variant | the planted one in 19 of 43, otherwise a median 0.5 kb away (LD) |

Test fixtures: `config/samplesheets/simtest.csv`, `results/04-expression/simtest/`, `results/05-eqtl/simtest/` (incl. `planted.tsv`). Bugs found and fixed in the test:
- `zcat | head` under `pipefail`;
- plink2 writing `22` instead of `chr22` (`--output-chr chr26`).

**Running on MAGE (2026-10-23/24).**
- *Interim cohort:* a `mage420` cohort (401 individuals, 425 runs completed by 2026-10-23 12:30) was set up for an early look: links to the completed runs, its own samplesheet, `results/05-eqtl/mage420/`. Its genotype preparation failed twice:
  1. plink2 needs sex for chrX → `prepare-genotypes.sbatch` now writes `sex.txt` from the samplesheet and passes `--update-sex`;
  2. chrX pseudoautosomal variants need `--split-par` with reference-specific boundaries.
- *Switch to the full cohort:* MAGE then finished, so the analysis moved to the full cohort; `mage420` was removed on 2026-10-27 (results, links and samplesheet; section 16).
- **Full MAGE, autosomes first** (chrX to be added with the correct PAR boundaries for T2T and GRCh38):
  - phenotypes (STAR): 17,160 expressed genes on T2T, 17,580 on GRCh38, 17,573 on patched GRCh38 (731 individuals);
  - genotypes chr1–22: jobs 365584 (T2T) and 365585 (GRCh38).
  - covariates: 5 genotype PCs, sex, 60 expression PCs (66 in total);
  - tensorQTL 365591–365593 (array 0–21, two tasks at a time per reference), done 2026-10-25 00:09 UTC with no failures; results in 13B.3.

### 13B.3 First results: MAGE, autosomes, STAR, matched variants (2026-10-25)

**eGenes per reference** (FDR 0.05, `collect-eqtl.py`; 1 − π0 = 0.89–0.90):

| Reference | Genes tested | eGenes |
|---|---:|---:|
| T2T-CHM13 | 16,605 | 15,013 |
| GRCh38 | 17,026 | 15,259 |
| Patched GRCh38 | 17,021 | 15,224 |

The MAGE paper (Taylor et al. 2024) reports about 15,000 eGenes, so the pipeline behaves as expected.

**Genes tested on both T2T and GRCh38** (matched by GeneID): 16,236 genes.
- 14,338 are eGenes on both.
- 379 are eGenes on T2T only and 274 on GRCh38 only. Almost all are threshold effects: the lead variant still associates on the other reference.
- With the replication rule (the lead variant has p ≥ 0.05 on the other reference), 32 are T2T-only and 29 GRCh38-only. On the other reference these genes sit just above the cut-off (median q 0.06–0.07), so most are weak.
- The strongest cases, with q on the reference where they are eGenes and q on the other:
  - T2T-only: FAM27C (1e-25 vs 0.10; T2T-new sequence), TRIM74 (2e-21 vs 0.09; SD), NEURL2 (2e-14 vs 0.10).
  - GRCh38-only: LOC124904395 (8e-15; SD), U2AF1 (6e-06; FixItFelix region).

**Genes annotated on both references but tested (expressed) on only one.** These are a reference effect, because reads land elsewhere or nowhere.
- 138 eGenes are tested on T2T only: 102 without an artefact flag (71 other, 25 SD, 3 FixItFelix region, 3 T2T-new). Examples: PLGLB1/2, TBC1D3, NPIPA2, GUSBP3, APOBEC3A, GSTT2, OR2T5/OR2T35, SMIM11.
- 114 eGenes are tested on GRCh38 only: 80 without a flag (56 other, 21 SD, 3 T2T-new). Examples: NUTM2B, UGT2B15, PKD1L1, GOLGA2P7, AGAP10P, CTAGE6, TBC1D3D, NPIPB7, SNCA.

**Genes annotated on only one reference** account for most of the raw "reference-only" counts in `summary.tsv`: 533 GRCh38-only and 158 T2T-only, mostly lncRNAs. These reflect the RefSeq annotation, not the reference sequence, and are reported separately rather than counted for H0.

**Positive control: GRCh38 vs patched GRCh38.** Only 8 + 2 genes differ, almost all inside FixItFelix's patched regions: U2AF1 and five ffx_patched LOCs (GRCh38 only), and SMIM11 and DNMT3L-AS1 (patched only). The known artefacts are reproduced as inflated signals rather than present/absent eGenes:

| Gene | q, T2T | q, GRCh38 | q, patched GRCh38 |
|---|---:|---:|---:|
| PWP2 | 5e-34 | 2e-210 | 2e-32 |
| CBS | 7e-04 | 6e-19 | 1e-03 |
| U2AF1 | 0.06 | 6e-06 | 0.07 |

So T2T agrees with the patch where the patch applies. This is the FixItFelix result, reproduced in an independent and larger cohort.

**What this means for H0 (first pass).**
- Outside the patched regions, T2T and GRCh38 differ for about 60 genes tested on both references, plus about 180 clean eGenes that are expressed on only one reference.
- These genes concentrate in SDs and gene families (TBC1D3, NPIP, GOLGA, NUTM2, PLG-like, FAM27), not in FixItFelix's regions. That supports H0, but each gene still needs checking.
- The count of present/absent eGenes misses the main kind of artefact: PWP2 and CBS are eGenes on every reference but much too strong on GRCh38.

**Note on lead variants.** The "slope sign agreement" between the two references' own lead variants (0.71) is not meaningful: the leads are different variants, and their alleles can be coded the other way round. The same-variant comparison below replaces it.

**Same-variant effect sizes** (`scripts/06-compare/compare-effects.py`, job 365661, 2 min, 12 GB, 2026-10-25).

*Method.* For every gene tested on both references and an eGene on at least one:
- each reference's lead variant is looked up for the same gene on the other reference;
- the slope sign is flipped where the alleles are coded the other way round (`allele_match` = swapped);
- z_diff = (b_A − b_B) / √(se_A² + se_B²), keeping the lead with the larger |z|; this is conservative, because both estimates come from the same people and genotypes;
- a gene is "changed" if the Bonferroni-corrected p is below 0.05.

*Check of allele handling,* T2T vs patched GRCh38, all compared genes:

| Allele match | Pairs | Correlation | Same sign |
|---|---:|---:|---:|
| same | 9,825 | 0.992 | 99.6% |
| swapped | 4,009 | 0.952 | 98.4% |
| nearby indel | 1,097 | 0.943 | 97.5% |

*Results.*

| Comparison | Genes compared | Changed | Without artefact flag | Of these: other / SD / FixItFelix region |
|---|---:|---:|---:|---:|
| Patched GRCh38 vs GRCh38 (positive control) | 15,447 | 9 | 7 | 0 / 0 / 7 |
| T2T vs GRCh38 | 14,976 | 168 | 131 | 69 / 56 / 6 |
| T2T vs patched GRCh38 | 14,931 | 163 | 127 | 71 / 56 / 0 |

- **Positive control.** All 9 genes changed by the patch lie in FixItFelix's regions, and all are stronger on unpatched GRCh38: GATD3 (z −39), PWP2 (−28), LOC124905035, LINC01678, ICOSLG, CBS, MRPL23, U2AF1, LOC102725065. This is FixItFelix's result, reproduced.
- **H0.** Against the patched reference, T2T still changes 127 genes without artefact flags (the 36 flagged include 9 HLA genes), and none of them lies in a patched region:
  - 69 are stronger on T2T, e.g. SIRPB1, TCAF2/TCAF2P1, ASAH2, NBPF10, LRRC37A, NPIPB13, STAG3L3, CCZ1, UBC;
  - 58 are stronger on GRCh38, e.g. RPL9, RPL13, RPS2, STAG3L1, PMS2/PMS2CL, SMG1P6, RASA4B, TBC1D3E, GOLGA8N, NUTM2E, ASAH2B, CCZ1B.
- **Paralog pairs.** About 44% are in SDs. The changes often come in pairs of the same gene family going opposite ways: ASAH2 stronger on T2T while ASAH2B is stronger on GRCh38, and likewise CCZ1/CCZ1B and TCAF2/TCAF2P1. This fits reads moving between paralogs.
- **Opposite signs.** 70 of the 127 genes have opposite-sign effects at the same variant, which also fits reads moving between paralogs. Ribosomal-protein genes (RPL9, RPL13, RPS2) probably reflect processed pseudogene copies; their region class is "other" because the copies are not SDs.
- These are the first H0 candidates; the salmon check comes next.

**Salmon and the robust set** (2026-10-25).

*Jobs.*
- Phenotypes and covariates: job 365662.
- tensorQTL: jobs 365663–365665 (22 tasks each, 0 failed).
- Comparisons: `compare-eqtl.sbatch mage salmon`, job 365789 (14 + 2 min).
- Combination of STAR and salmon: `scripts/06-compare/combine-methods.py`.

*eGenes with salmon.* 14,574 on T2T (16,032 genes tested), 14,716 on GRCh38 (16,342) and 14,794 on patched GRCh38 (16,332).

*Changed genes.* A gene is robustly changed when STAR and salmon both call it changed with the same reference stronger. The counts are after the allele-sign correction below.

| Pair | STAR | salmon | Both | Both, no artefact flag (SD / other / FixItFelix) |
|---|---:|---:|---:|---:|
| Patched GRCh38 vs GRCh38 | 9 | 9 | 8 | 6 (0 / 0 / 6) |
| T2T vs GRCh38 | 145 | 563 | 79 | 63 (34 / 24 / 5) |
| T2T vs patched GRCh38 | 138 | 549 | 71 | **57 (33 / 24 / 0)** |

- Salmon calls many more genes changed (549), but most of its calls show no difference with STAR (before the correction: 440 salmon-only genes, median STAR |z| 0.2). They are mostly protein-coding genes, e.g. ATF4, MRPL18, CTDNEP1, HIF1A-AS2.
- Our reading: salmon splits reads between transcript models, and the RefSeq transcript sets differ between the two assemblies, so salmon also measures annotation differences. That is why a gene must be called changed by both methods.
- The positive control holds under both methods: 8 of 9 patched genes, including GATD3, PWP2, CBS, ICOSLG and U2AF1.

*Allele-sign correction (2026-10-25).* A first version used the variant-matching label `allele_match` = swapped to flip slope signs. Checking the 69 genes it reported against the genotypes themselves showed problems:
- Some indels are written differently by the two panels: T2T `C:CA` vs GRCh38 `CA:C`, labelled "swapped", although the genotypes correlate at r = +1.
- So 11 genes had identical effects turned into opposite-sign "changes": SLC25A25, MRPS21, HS2ST1, CHP1, OXR1, LOC124902367, ZFAND6, CARF, SLC50A1, SMIM29 and PIKFYVE. Their per-person expression correlated at r ≥ 0.998 between references.
- One gene (HTATIP2) paired two different indels, with genotype r = 0.03.
- Genome-wide, the label and the genotype correlation disagree in 133 of 90,760 lookups (0.15%), and 212 have |r| < 0.9. These few errors pile up among the "changed" genes, because a wrong sign looks like a large change.

`compare-effects.py` now takes the sign from the genotype correlation in the cohort and skips pairs with |r| < 0.9. Reruns: jobs 366132 (STAR) and 366133 (salmon).

**The robust H0 set: 57 genes** that T2T changes relative to patched GRCh38, without artefact flags. None lies in a FixItFelix region. File: `results/06-compare/eqtl-comparison/mage/robust-matched/changed-chm13v2.0-vs-grch38-ffx.tsv`.
- Composition: 35 protein-coding, 13 transcribed pseudogenes and 9 lncRNAs; 33 in SDs; 31 with opposite-sign effects at the same variant.
- Checks on these 57:
  - the genotypes of the compared variant agree between panels (|r| ≥ 0.96);
  - the gene's per-person expression differs between references (median r 0.70; 43 genes below 0.9; none above 0.99).
  - So each change comes from the expression measurement, i.e. from where reads are assigned, not from the genotypes.
- Strongest differences, with per-person expression correlation r between references: ASAH2 (r −0.26), STAG3L1 (0.15), CCZ1/CCZ1B (0.27/0.26), TCAF2P1 (0.27), TBC1D3E, RGPD4, RASA4DP/RASA4B, NUTM2D/NUTM2E, NPIPB13, PMS2/PMS2CL, SMN2.
- The remaining protein-coding genes outside SDs (UBC, ELOB, H3-3A, GPX1, RPL9, RPL13, RPS2, RPSA, EIF4B, PTGES3, PCCA, IFITM3) are highly expressed and keep r 0.86–0.99. They fit a small share of reads moving to processed-pseudogene copies; this needs checking per gene.
- **Paralog pairs that switch:** ASAH2/ASAH2B, CCZ1/CCZ1B, TCAF2/TCAF2P1, STAG3L1/STAG3L3, PMS2/PMS2CL, PDPR/PDPR2P, NUTM2A/D vs NUTM2E, RASA4B/RASA4DP, ANAPC1/ANAPC1P2.

**Read destination** (`scripts/06-compare/read-destination.py`, job 366433, 13 s, 2 GB; 2026-10-25).

*Method.* For each of the 57 genes:
- per individual, delta = CPM on T2T − CPM on patched GRCh38, with genes matched by GeneID; a gene missing on one reference counts as 0 there;
- the partner is the gene whose delta moves most strongly the other way across individuals (most negative correlation);
- partner_frac is the share of the gene's delta that the partner takes.

Output: `results/06-compare/eqtl-comparison/mage/robust-matched/read-destination-chm13v2.0-vs-grch38-ffx.tsv`.

*Results.* 33 of 57 genes have a partner with r < −0.5, and 16 have r < −0.8.
- **Paralog swaps** (r −0.7 to −0.96):
  - CCZ1↔CCZ1B (−0.96)
  - ASAH2↔ASAH2B (−0.91)
  - IFITM3↔IFITM2 (−0.92)
  - SDHAP1↔SDHAP2 (−0.88)
  - LOC642846↔DDX12P (−0.86)
  - WHAMMP1↔WHAMMP4 (−0.83)
  - BMP8A↔BMP8B (−0.81)
  - PMS2↔PMS2CL (−0.80)
  - DHRS4L2↔DHRS4 (−0.79)
  - GTF2IP4↔GTF2IP1 (−0.73)
  - STAG3L1↔STAG3L2/STAG3L3 (−0.70/−0.66)
- **Processed pseudogenes:** RPSA↔RPSA2 (−0.94; T2T gives RPSA 105 CPM fewer reads, 46% of them to RPSA2) and H3-3A↔H3P6 (−0.71). This is the mechanism suggested for the highly expressed genes.
- **Reads taken by genes that exist only on T2T:**
  - SMG1P6 → LOC124907836 in T2T-new sequence (−0.94);
  - PDPR and PDPR2P → LOC128966601 in T2T-new sequence (−0.82/−0.59). PDPR has 64.5 CPM on GRCh38 and 15.2 on T2T, but LOC128966601 takes only 7.5% of that difference (partner_frac 0.075; about 3.7 of 49 CPM), so most of PDPR's reads on T2T go elsewhere (not yet traced). For PDPR2P it takes 95%;
  - IFITM3-AS1 → LOC124906937, a T2T-only annotation (−0.96);
  - LGALS9B and GLIDR → T2T-new-sequence genes (−0.66).
  - In these cases GRCh38 lacks the sequence, so its reads pile onto the nearest paralog: a direct example of missing sequence creating an eQTL artefact. Note: "33 of 57 have a partner" means a single gene mirrors the change (r < −0.5), not that it takes all the reads; the share ranges from 7.5% (PDPR) to about 100% (CCZ1B 86%, PMS2CL 98%, PDPR2P 95%). Wording corrected in the manuscript and talk materials 2026-10-26.
- **No clear partner** (r > −0.4): RPL9, RPL13, RPS2, GPX1, PTGES3, UBC, EIF4B, ELOB, SMN2, NUTM2A/D/E, RASA4B/RASA4DP and TCAF2/TCAF2P1 (−0.40).
  - Their reads may go to several copies at once, or become multi-mapped and uncounted. Checking that needs the BAMs or per-contig counts.
  - SMN2↔SMN1 is only −0.25, as expected: the two copies are nearly identical and both references contain both.
- **Caution:** LOC124906110 (39 vs 3 CPM) correlates at about −0.84 with many unrelated genes (ANKRD17, SLAIN2, ZFC3H1). That looks like a sample-wide technical factor, not a paralog, so treat it with caution.

**All variants and H1** (2026-10-25). tensorQTL with all variants (`--variants all`, jobs 366424–366426, 0 failed): 14,954 / 15,213 / 15,119 eGenes on T2T / GRCh38 / patched GRCh38. Lead variants of the 14,954 T2T eGenes:

| Category of the lead variant | eGenes |
|---|---:|
| Matched on both panels | 13,111 |
| T2T variant not called in the GRCh38 panel (mostly indels removed by panel filters) | 1,319 |
| Not in the matching table (likely below MAF 0.01 in the 1,089 individuals, or truncated IDs of long alleles; not checked) | 270 |
| Paired, but the allele frequency differs (paralogous) | 195 |
| **In T2T-only sequence (H1)** | **59** |

- 59 genes have their lead variant in sequence GRCh38 lacks. 9 of them are not eGenes when only matched variants are used.
- Others become much stronger with the T2T-only variant, e.g. IVNS1ABP (p 1e-82 vs 1e-34), LOC128966600 in new chr9 sequence (q 1e-51 vs 0.007), LOC105373279 and OR2T34 (chr1 OR2T cluster), LOC101060632 and CBSLR.
- Most stay eGenes through linked shared variants. So H1 adds examples, not the main result, as predicted in 4.4.
- Script: ad hoc, job 366638. Table: `results/06-compare/eqtl-comparison/mage/h1-t2t-only-lead-variants.tsv`.

**H1 lead-variant quality check** (2026-10-26, `h1-variant-qc.py` via `h1-variant-qc.sbatch`, job 370504 after two failed attempts 370480/370482; 43 min).

The phased T2T panel carries no per-variant quality scores, only PASS sites and genotypes. Quality is therefore judged from three things:
- **Genotypes in the 731 MAGE individuals.** Heterozygosity within each population: an excess of heterozygotes means two paralogues were called as one site. Hardy–Weinberg p within populations, and whether the site is multiallelic in the 3202 panel.
- **Region.** censat satellite arrays, SD identity, RepeatMasker class, a length change in a tandem repeat, and whether the variant sits inside an exon of the gene.
- **Signal.** The best shared (matched) variant for the gene, and the lead's p-value after conditioning on it.

Calls:
- *doubtful*: heterozygote excess (F < −0.1), HWE p < 1e-6, or a length change inside a satellite array;
- *borderline*: any other flag;
- *solid*: no flag.

| Genotype call | Genes | T2T-only variant carries the signal (alone or with a shared variant) |
|---|---:|---:|
| Solid | 22 | 15 |
| Borderline | 25 | 14 |
| Doubtful | 12 | 5 |

- **The 15 solid genes with their own signal are the defensible H1 examples.** They are IVNS1ABP, CBSLR and LOC105371651 (one 5-bp deletion at chr1:184763743), the chr1 OR2T cluster (OR2T2, LOC105373279, LOC128966579; one SNV), SHC1, WASH8P, ADSS1, CAPN15, VPS35L, KIF1A, PDCD2, TRIP6 and MNX1-AS2.
- **Most borderline calls only reflect a multiallelic site in the panel.** In practice this means a nearby second allele, often a repeat.
- **Doubtful cases are mostly length changes in HSat3 (GGAAT)n arrays.** Examples: ZNF25, HM13-AS1, LOC105379426, FRYL, and three chr9 lncRNAs. The other is OR2T34, with heterozygote excess.
- **LOC128966600** (new chr9 sequence) is borderline: a (TG)n indel in a >99% SD. Its expression change is a reference effect regardless.
- **22 of 59 are explained by, or inseparable from, shared variants.** Examples: CIP2A, DLG1, LOC105374445, and EHMT1, a perfect proxy. In these genes the T2T-only variant is just the best tag.
- **The "9 eGenes only with T2T-only variants" are weak.**
  - Six have a matched-variant q of 0.051–0.070, just above the cutoff.
  - Three are genes in new sequence with no matched variant in cis (LOC124905534, LOC124908428, LOC124908450); two of these are doubtful.
  - The cleanest are VPS35L and KIF1A (solid SNVs), but both are marginal (q 0.014 and 0.005).
  - So H1 should be presented as **about 15 well-genotyped examples where the new variant carries signal**, not as new eGenes.
- Table: `results/06-compare/eqtl-comparison/mage/h1-variant-qc.tsv`. `--classify-only` re-applies the calls without recomputing.

**H1 replication in Geuvadis** (2026-10-28; `replicate-h1.py`, job 373588, 4 min, after the Geuvadis all-variants tensorQTL, STAR, jobs 373575–373580: 13,481 / 13,676 / 13,712 eGenes on T2T / GRCh38 / patched GRCh38).
- Design as for H0: each MAGE H1 lead variant is looked up at the same variant in Geuvadis (same panel, same IDs and alleles), and its effect is retested conditional on MAGE's best shared variant. *Replicated* = same sign and p < 0.05/53; *signal replicated* = also conditional p < 0.05/53.
- 53 of 59 testable. Untested (variant or gene not testable in Geuvadis): ADAMTSL4, OR2T5, SUMO1, TNNC2, LOC124908450, FRYL.

| MAGE genotype call | Genes | Tested | Same sign | Replicated | Signal replicated |
|---|---:|---:|---:|---:|---:|
| Solid | 22 | 22 | 18 | 14 | 10 |
| Borderline | 25 | 21 | 17 | 11 | 4 |
| Doubtful | 12 | 10 | 7 | 3 | 1 |
| *Solid, T2T-only variant carries the signal in MAGE* | *15* | *15* | *11* | *9* | *8* |

- **8 of the 15 solid examples replicate with their own signal, at 5 distinct variants.**
  - IVNS1ABP and CBSLR: the 5-bp deletion chr1:184763743.
  - The chr1 OR2T cluster: OR2T2, LOC105373279 and LOC128966579 share one SNV, chr1:247973939.
  - WASH8P (chr12:29314 deletion), KIF1A (chr2:241323426) and PDCD2 (chr6:171889515).
- TRIP6 replicates the association (p 8e-7) but not clearly beyond the shared variant (conditional p 0.025).
- **Not replicated:** SHC1, ADSS1, VPS35L, MNX1-AS2, CAPN15 and LOC105371651.
  - The first four were marginal in MAGE (gene-level permutation p 0.1–0.24).
  - CAPN15 (MAF 0.03) and LOC105371651 were strong in MAGE, but in Geuvadis their T2T-only variant has no effect, although both genes are eGenes there through other variants. LOC105371651 shares its variant with CBSLR and IVNS1ABP, which do replicate.
  - Allele frequencies are similar in the two cohorts (e.g. CAPN15 0.028 vs 0.029, VPS35L 0.19 vs 0.10), so lower power explains at most VPS35L.
- **H1 is therefore about 5 replicated loci (8 genes)**, all genotyped without quality flags, where a variant in sequence GRCh38 lacks carries eQTL signal in both cohorts. This is a set of examples, not a genome-wide effect, as expected in 4.4.
- Table: `results/06-compare/eqtl-comparison/replication-geuvadis/h1-genes.tsv`.

**Independent replication: without the 93 individuals in both cohorts** (2026-10-28; `replicate-changes.py --exclude-shared` and `replicate-h1.py --exclude-shared`, job 377254, 17 min).
- **Why.** 93 of Geuvadis's 451 individuals are also in MAGE. Their RNA was sequenced separately, but the genotypes are the same people, so the replication above is not fully independent.
- **How.** The 93 are dropped, and the nominal effect of each variant is recomputed on the other 358, as tensorQTL map_nominal computes it: genotype and expression residualised on the covariates, dof = n − 2 − covariates.
  - The step-5 phenotypes and covariates are reused. They were normalised, and the expression PCs computed, on all 451.
  - Check (`replicate-changes.py --check`, login node, 13 min, 72 GB): on all 451 individuals the recomputation matches the tensorQTL parquet values to 3×10⁻⁷ (slope and se) for all 12 reference × method sets.
- **Results: unchanged.** T2T vs patched GRCh38, 54 genes tested:

| | 451 individuals | 358 independent |
|---|---:|---:|
| Same sign (STAR and salmon) | 53 | 54 |
| Nominal (p < 0.05, both methods) | 50 | 49 |
| Bonferroni | 44 | 42 |
| Sign test p (STAR) | 3×10⁻¹⁵ | 6×10⁻¹⁷ |
| Reversals that reverse again (of 12) | 10 | 10 (TCAF2, ELOB one-sided again) |
| H1, solid with own signal (of 15): replicated / signal replicated | 9 / 8 | 9 / 8 (same genes) |

  - The genes that do not replicate are the same: UBC, PPP4R3B-DT, PTGES3, PDPR2P, plus IFITM3-AS1, which falls just short (STAR p 0.026, salmon p 0.073).
  - Outputs: `results/06-compare/eqtl-comparison/replication-geuvadis-independent/`.

**Paired test of effect changes** (2026-10-28; `scripts/06-compare/paired-test.py`; MAGE job 377313, 17 min; Geuvadis replication job 377411, 12 min).
- **Why.** z_diff divides by √(se_A² + se_B²), as if the two slopes came from different people. They come from the same individuals, and most of a gene's expression is measured identically on both references (median correlation of the two references' residual expression ρ = 0.99). The two estimates are therefore strongly correlated, and z_diff was far too conservative.
- **How.** At the variant pair compare-effects.py chose for each gene, both slopes are recomputed as tensorQTL does (they match to 5×10⁻⁴). The covariance of two OLS slopes whose errors are correlated within an individual is cov = [Σ e_A e_B / dof] · Σ g_A g_B / (Σ g_A² Σ g_B²). Then z_paired = (slope_A − slope_B) / √(se_A² + se_B² − 2 cov).
- **Calibration.**
  - Bootstrap over individuals (300 random genes × 200 resamples per pair and method): the bootstrap SD of slope_A − slope_B is 0.94–0.96 × the paired se, so the paired test is about 5% conservative.
  - The independent se is 9–15 × too large. Its genome-wide λ is 0.01–0.03, against 1.3–2.1 for the paired test.
  - Internal null: patched GRCh38 vs GRCh38 under salmon (they differ only in the patched regions) gives λ = 1.03.
  - T2T comparisons give λ 1.6–2.1. With the bootstrap agreeing, this reflects many small real shifts across the genome, not miscalibration.
- **Results, MAGE, T2T vs patched GRCh38:** 297 genes changed under both methods with the same sign (255 without artefact flag), against 71 (57) with z_diff.
  - All 57 are included. The 198 new genes have smaller changes (median |Δslope| 0.063, median 19% of the larger slope). 124 of the 255 have |Δslope| ≥ 0.1.
  - New genes are mostly protein-coding (156), and 35 are in SDs. Examples: SMN1, GOLGA8O, GOLGA6B, NPIPB6, LRRC37A2, SERF1B, ARMC10, GCSH.
- **Replication in Geuvadis** (`--at mage`: the 255 unflagged genes at MAGE's variant pairs, on the 358 individuals not in MAGE; 247 testable):

| Genes (T2T vs patched GRCh38) | Tested | Same sign (STAR and salmon) | Nominal | Bonferroni |
|---|---:|---:|---:|---:|
| The 57-gene set | 54 | 54 | 53 | 51 |
| New with the paired test | 193 | 167 (sign test p 1×10⁻²⁶) | 122 | 91 |
| new, MAGE \|Δslope\| ≥ 0.1 | 63 | 58 | 51 | 39 |
| new, MAGE \|Δslope\| < 0.1 | 130 | 109 | 71 | 52 |

  - With the paired test, the 57-gene set replicates more strongly: 51 Bonferroni, against 42 with z_diff.
  - The large new changes replicate about as well as the original set (51 of 63 nominal). The small ones replicate less often but still mostly (71 of 130).
- Outputs: `results/06-compare/eqtl-comparison/{mage/paired-matched/, geuvadis/paired-matched-at-mage/}` (effects per pair and method, robust-*.tsv, calibration.tsv, summary.tsv, replication-*.tsv).

**ctyper copy number as a covariate** (2026-10-28; `ctyper-copy-number.py`, job 377320, 20 s; `paired-test.py --copy-number`, job 377358, 17 min).
- **Copy numbers.** ctyper's 1KGP table (`data/external/ctyper-ma2025/PangenomeAlleles_1kg.tar.gz`) was extracted to `/scratch/$USER/t2t-eqtl/ctyper/` (1.1 GB, not kept in the project). From it, per gene symbol:
  - *gene* copy number: allele types whose main gene is the symbol;
  - *family* copy number: all allele types of every ctyper group that mentions the symbol.
- **Coverage.** 2,657 genes and 4,050 families, for all 1,089 cohort individuals (`results/06-compare/ctyper/copy-number.tsv.gz`).
  - 34 of the 57 genes are covered, including 10 of the 13 reversals.
  - Not covered: UBC, H3-3A, the RPL/RPS genes, STAG3L1, PPIP5K1P1-CATSPER2 and most LOC lncRNAs. ctyper's gene set does not include processed-pseudogene families or RefSeq-only genes.
- **Method.** Both copy numbers, where they vary, are added as covariates for that gene on both references, and the paired test is rerun (1,784 genes with copy number).
- **Result: copy number does not explain the changes.**
  - 32 of the 34 covered genes of the 57-gene set stay changed under both methods, with the same sign. The median Δslope with copy number / without is 1.00.
  - The lead variants do not tag copy number: copy number explains a median 0.3% of the genotype variance, and > 50% only for NUTM2D (69%), whose change shrinks from −0.54 to −0.36. It stays nominally significant (STAR z −3.5, salmon −4.8) but falls below Bonferroni with STAR, so it is lost. (Corrected 2026-10-29: this line first said NUTM2D "stays significant" and 33 of 34.)
  - Also lost: PDPR2P (z 9.1 → 3.1), which also failed replication in Geuvadis.
  - Main set (117, checked 2026-10-29 with `prepare-figure-data.py` → `results/07-presentation/data/copy-number.tsv`): 62 covered, 55 stay changed with the same sign, median ratio 1.00. Lost: NUTM2D, PDPR2P, SMN1, SERF1B, NPIPB6, GRAPL, GSTM2 (GSTM2, NPIPB6 and SMN1 also fail replication in Geuvadis).
  - All 10 covered reversals are unchanged, which argues against the idea (reversals above) that the reversal variants tag a copy-number difference. Caveat: ctyper counts gene copies and allele types, not every structural variant near the gene.
  - Output: `results/06-compare/eqtl-comparison/mage/paired-matched-cn/`.

**Allele-specific expression: feasibility** (2026-10-29; `scripts/06-compare/ase-feasibility.py` + `.sbatch`, job 379963, 4 min, 12 GB; no new data).
- **Idea.** In a person heterozygous at the lead variant, the haplotype that carries the alt allele should give more RNA reads (slope > 0) or fewer (slope < 0) than the other haplotype. The test is within a person, so covariates and trans effects cancel, and it can decide which reference has the right direction at a reversal. Reads are assigned to haplotypes at heterozygous exon SNPs, phased with the lead in the panel. The panel is phased (pgen hard-call phase present for all 731).
- **Check, per gene of the 57.**
  - Panel SNPs in the gene's T2T exons.
  - Whether a 150 bp window around each SNP, with either allele, maps back to its own position with the best other hit ≥ 2 mismatches worse (minimap2 `-x sr`, score gap ≥ 20). Checked on T2T, and on patched GRCh38 via the matching table.
  - Individuals heterozygous at the lead and at ≥ 1 such SNP ("informative").
  - Approximate read depth per exonic base (mean STAR count × 2 × 150 / exonic length).
- **Result.** Rerun on the main set (job 380040, 2026-10-29): 59 of the 117 main genes have ≥ 20 informative individuals on both references (core: 28 of 57). Of the 13 true reversals:
  - testable with good depth: ELOB (242 informative, ~600 reads per base), TCAF2 (344, 14–25), RPS2 (141, > 29,000) and PMS2CL (58, 37–51);
  - weak: H3-3A (7 individuals, deep), and ASAH2 (139) and PPIP5K1P1-CATSPER2 (149) at only ~0.5–2 reads per base, so usable only if pooled;
  - not testable: STAG3L1, TCAF2P1, ANAPC1P2, CCZ1, RASA4DP, NUTM2E. They have no panel SNP in their exons at all.
- **Why.** The 1KGP short-read panel calls almost no variants in segmental duplications: 25 of the 33 SD genes of the 57 have no exon SNP (3 of 24 outside SDs). For example, the CCZ1 gene body has 2 SNPs in 28 kb, and STAG3L1 has none. So the cohort ASE route covers mostly the non-SD genes. The SD reversals need genotypes from assemblies (HPRC), where paralogs are resolved.
- Outputs: `results/06-compare/ase/feasibility-mage.tsv`, `snps-mage.tsv.gz` (per SNP: score gaps on both references).
- **Long-read genotypes for the SD genes** (checked 2026-10-29).
  - The 1KGP ONT release (Schloissnig et al. 2025, `ftp.1000genomes.ebi.ac.uk/.../1KG_ONT_VIENNA/`) has SV calls only (phased Giggles VCF, 74 MB) and no small variants. It does have per-person alignments to T2T and GRCh38 (`t2t/<id>.t2t.cram`, about 8.6 GB each).
  - Coverage of our cohorts: 260 of 731 MAGE and 121 of 451 Geuvadis individuals.
  - Remote region streaming works with their T2T reference (`data/external/1kg-ont-vienna/1KG_ONT_VIENNA_t2t.fa`, chr names). For CCZ1 in HG00096 it took 8 s, with 101 primary reads (median 4.8 kb, 90th percentile 21 kb) and mean depth 16 at MAPQ ≥ 20, covering 99.8% of the gene at ≥ 8×.
  - So SNPs in the SD genes could be called from long reads in about 260 MAGE individuals. Read-backed phasing would link them to the lead where the two lie within a read length.
- **Data for the ASE test** (started 2026-10-29).
  - Targets (`ase-targets.py`, login node): the 117 main genes plus 38 read partners; RNA regions = gene ± 5 kb on chm13v2.0 and grch38-ffx; ONT regions = main genes ± 20 kb and T2T lead ± 20 kb (123 regions, 10.2 Mb). Individuals: all 260 MAGE individuals with an ONT T2T CRAM, plus 87 chosen greedily to bring testable genes to ≥ 30 informative people: 347 individuals, 369 runs (`results/06-compare/ase/targets/`).
  - RNA realignment (`ase-realign.sbatch`, array 380491, 369 tasks, 10 at a time, ~24 min each): SRA fetch, STAR with the stage-04 settings on chm13v2.0 and grch38-ffx, keep only reads in the target regions (`results/06-compare/ase/rna/<run>/<ref>/targets.bam`, ~150–210 MB per run). Tests 380422_1-2.
  - ONT SNP calls (`ont-call.sbatch`, array 380492, 260 tasks, 20 at a time, ~15 min each; pixi environment `longread`: Clair3 2.0.3 model r941_prom_sup_g5014, WhatsHap 2.8): target slices streamed from the remote T2T CRAMs, Clair3 with read-backed phasing (`results/06-compare/ase/ont/<id>/`; 21–25k PASS variants per person). Tests 380429_1-2.
  - Next: allele counts at phased heterozygous exon SNPs on both references, reference-bias correction, direction against the eQTL slope.

**GWAS lookup** (2026-10-29; `scripts/06-compare/gwas-lookup.py` + `.sbatch`, job 380065, 6.5 min, 9 GB).
- **What it does.** A first screen for H6 (colocalisation comes next). Query variants per gene, as GRCh38 IDs: the compared pair, and each reference's own lead (STAR, matched variants). For each, it looks for GWAS Catalog associations within 1 Mb in LD with it. LD is r² in MAGE's GRCh38 genotypes, either all 731 or the 142 EUR individuals; a proxy needs r² ≥ 0.8. Associations must be genome-wide (p < 5×10⁻⁸).
  - GWAS Catalog `latest` fetched 2026-10-29 to `data/external/gwas-catalog/`: 1,203,840 associations, 987,857 single-SNP rows with a GRCh38 position.
- **Result, T2T vs patched GRCh38.**
  - Genes with a proxy hit: core 37 of 57, main 60 of 117, extent 137 of 255.
  - Among the 60 main genes, hits by trait type (a gene can count more than once): disease or behaviour 37, molecular levels (protein, metabolite) 33, other quantitative traits 40, blood cell traits 13.
  - In 31 main genes the hit is reached through only one reference's own lead.
- **Background** (`--background 5`, job 380162, 7 min). 494 unchanged genes, drawn per main gene from the same region class and the same quintile of eQTL strength. There are fewer than 5 per gene because unchanged SD genes are few.
  - **No enrichment.** 229 of 494 background genes have a proxy hit (46%), against 60 of 117 main genes (51%): OR 1.2, Fisher p 0.36. Disease-like traits: OR 1.4, p 0.10. Most eQTLs tag some GWAS hit in a catalogue this large, so "changed genes sit at GWAS loci" is not a claim we can make.
  - **The link depends on the reference much more often.** The hit is reached through only one reference's own lead in 31 of 117 main genes (26%; 16 through T2T, 15 through GRCh38), against 21 of 494 background genes (4%): OR 8.1, p 1×10⁻¹¹. This partly follows from the definition, since changed genes often have different leads on the two references. It is the reason to run coloc on both references.
- **Caveat.** Proximity is not colocalisation.
- **Candidates for coloc** (disease traits; hits through the compared pair mean the eQTL effect changes, through one lead only that the locus maps differently):
  - autoimmune thyroid disease and hypothyroidism: IFITM3, IFITM2, IFITM3-AS1 (GRCh38 lead only), PPP4R3B-DT, PTGES3 (T2T lead only);
  - systemic lupus: LOC124906110 (GRCh38 lead);
  - ulcerative colitis, multiple sclerosis, Hodgkin's lymphoma: NPIPB13, SMG1P6;
  - type 2 diabetes: NBPF10, RPL13, PPIP5K1P1-CATSPER2 (reversal; also lung cancer);
  - atrial fibrillation: RPS2 (reversal), PTGES3, LOC124903630;
  - alcohol consumption: LRRC37A2 (17q21.31);
  - the PMS2/PMS2CL reversal: one weak hit only (Alzheimer's × APOE, p 3×10⁻⁸, r² 0.83).
- Outputs: `results/06-compare/gwas/{lookup,genes}-mage-chm13v2.0-vs-grch38-ffx.tsv`.

**Colocalisation with GWAS on each reference (H6)** (2026-10-29; `scripts/06-compare/coloc.py` + `.sbatch`; test job 380214 on hypothyroidism, full job 380275, 38 min, 9 GB).
- **Method.** coloc's approximate Bayes factor test (Giambartolomei et al. 2014), written out in the script and checked on simulated signals.
  - It assumes one causal variant per trait. Priors are coloc's defaults: p1 = p2 = 10⁻⁴, p12 = 10⁻⁵, prior effect SD 0.15 (eQTL) and 0.2 (log OR).
  - Each gene's whole cis window, STAR, matched variants. T2T IDs are converted to GRCh38, so both references are tested on the same variants.
  - PP4 = one shared causal variant; PP3 = two different causal variants. A call is "only" on a reference when PP4 ≥ 0.8 there and < 0.5 on the other.
- **GWAS** (European ancestry, harmonised GRCh38 summary statistics from the GWAS Catalog, `data/external/gwas-sumstats/`):

| Trait | Study | Accession |
|---|---|---|
| Hypothyroidism | White 2026 | GCST90627749 |
| Systemic lupus | Bentham 2015 | GCST003156 |
| Ulcerative colitis | de Lange 2017 | GCST004133 |
| Multiple sclerosis | Shigesi 2025 | GCST90558092 (1.6 M variants only) |
| Type 2 diabetes | Xue 2018 | GCST006867 |
| Atrial fibrillation | Yuan 2025 | GCST90624412 |

  - GCST007517 (Mahajan 2018) was downloaded first but is a targeted-array study (131 k variants). It is not used.
- **Result.**
  - Of the 255 extent genes, 15 colocalise (PP4 ≥ 0.8) with at least one trait on at least one reference, and 8 of these on one reference only. Of the 117 main genes, all 6 that colocalise do so on one reference only.
  - Where one reference loses the colocalisation, its PP3 is usually high: both signals are still there, but the eQTL now points to a different variant.
  - No colocalised gene changes the direction of the gene–trait link.

| Gene (tier) | Trait | Colocalises on | PP4 T2T / GRCh38-ffx | PP3 on the other | GWAS min p |
|---|---|---|---|---:|---:|
| PTGES3 (core) | Atrial fibrillation | T2T | 0.93 / 0.003 | 0.99 | 1×10⁻⁴⁰ |
| POM121C (core, SD) | Atrial fibrillation | T2T | 0.99 / 0.03 | 0.97 | 2×10⁻¹⁸ |
| LOC124903630 (core) | Atrial fibrillation | GRCh38-ffx | 0.04 / 0.999 | 0.65 | 1×10⁻¹⁴ |
| IFITM3-AS1 (core) | Hypothyroidism | GRCh38-ffx | 0.0002 / 0.999 | 1.00 | 8×10⁻¹⁸ |
| SMG1P6 (core, SD) | Ulcerative colitis | GRCh38-ffx | 0.0003 / 0.83 | 0.99 | 8×10⁻⁸ |
| SERF1B (main, SD) | Hypothyroidism | GRCh38-ffx | 0.01 / 0.81 | 0.05 | 4×10⁻⁶ |

- **The same locus pointing to different genes.** At 16p11.2, the ulcerative colitis variant chr16:30328432 is the shared variant for SMG1P6 on GRCh38 (PP4 0.83) and for its neighbour NPIPB13 on T2T (PP4 0.97, against 0.56 on GRCh38). The GWAS signal there is modest (min p 8×10⁻⁸).
- **Also of note (suggestive).**
  - RPS2 (reversal), atrial fibrillation: PP4 ≈ 0 on T2T and 0.53 on GRCh38; GWAS p 2×10⁻³¹.
  - GPX1, ulcerative colitis: 0.74 on T2T, ≈ 0 on GRCh38.
  - BMP8A, type 2 diabetes: ≈ 0 on T2T, 0.73 on GRCh38.
- **Caveats.**
  - One causal variant per trait is assumed. SuSiE-coloc (below) relaxes this at the 32 loci with a signal.
  - The GWAS are European and the eQTL cohort is multi-ancestry. This affects both references equally, so the comparison between them stays fair.
  - Which reference is right is not decided here; that needs allele-specific expression or own-genome evidence.
- Output: `results/06-compare/gwas/coloc-mage-chm13v2.0-vs-grch38-ffx.tsv` (gene × trait, both references); test run `...-test.tsv`.

**SuSiE-coloc at the colocalisation loci** (2026-10-29; `scripts/06-compare/coloc-susie.py` + `.sbatch`; job 380902, 2 h 04 min, 8 GB, after two failed attempts 380440 and 380490).
- **Why.** coloc.abf assumes one causal variant per trait. Where a region has several signals, it can call a colocalisation that is really two signals, or miss one. SuSiE (Wang et al. 2020) fine-maps each signal into credible sets, and coloc.susie (Wallace 2021) tests every pair of eQTL and GWAS credible sets.
- **How.** The 32 gene × trait pairs with PP4 ≥ 0.5 on either reference. Window: ±500 kb around the two references' top variants; SNPs only, MAF ≥ 0.01 in the cohort and in EUR. LD: the MAGE GRCh38 genotypes for the eQTL (in-sample), and 525 unrelated 1KGP EUR individuals for the GWAS (out-of-sample; the GWAS are European). R coloc 5.2.3, susieR 0.14 (`runsusie`, maxit 1000), coloc's default priors.
- **Fixes.** 380440: a pandas list & Series operation. 380490: R's `fread` cannot read `.gz` without R.utils, and read the numeric header of the LD matrices as data; LD is now written as plain TSV without a header (tested on a simulated locus first).
- **Result.** SuSiE-coloc can be run only where SuSiE finds a credible set for both the eQTL and the GWAS.

| Gene, trait | coloc.abf call | coloc.abf PP4 T2T / ffx | SuSiE-coloc PP4 T2T / ffx | Reading |
|---|---|---|---|---|
| IFITM3-AS1, hypothyroidism | ffx only | 0.0002 / 0.999 | 0.004 / 0.998 | **confirmed**; ffx eQTL credible set contains the GWAS lead chr11:324163 |
| LOC124903630, atrial fibrillation | ffx only | 0.04 / 0.999 | no T2T eQTL set / 0.91 | **confirmed** |
| SERF1B, hypothyroidism | ffx only | 0.01 / 0.81 | no T2T eQTL set / 0.86 | **confirmed** |
| PTGES3, atrial fibrillation | T2T only | 0.93 / 0.003 | 0.03 / 0.000 (PP3 0.87 / 1.00) | **not confirmed**: different variants on both |
| SCAF11, hypothyroidism | T2T only | 0.94 / 0.002 | 0.000 / 0.000 (PP3 0.99) | **not confirmed** |
| POM121C, AF; SMG1P6, UC; RMI1, hypothyroidism; SCAF11, AF | one reference | — | no GWAS credible set | not testable |
| IRF5, UC / hypothyroidism | both | 0.99 / 0.99; 0.998 / 0.998 | 0.98 / 0.98; 0.61 / 0.61 | consistent |
| TOP3B, hypothyroidism | both | 0.92 / 0.92 | 0.86 / 0.86 | consistent |
| SEC14L1, hypothyroidism | both | 0.997 / 0.95 | 0.02 / 0.03 (PP3 0.97) | not confirmed on either |
| PPP4R3B-DT, hypothyroidism | suggestive | 0.75 / 0.95 | 0.39 / 0.91 | leans to ffx |

- **Reading.** 3 of 5 testable one-reference colocalisations hold with several causal variants allowed, all on the patched GRCh38 (IFITM3-AS1, LOC124903630, SERF1B). PTGES3 (the strongest GWAS, p 1×10⁻⁴⁰) does not: its T2T-only call came from the one-variant assumption. The headline H6 example is therefore IFITM3-AS1, not PTGES3 (manuscript changed accordingly).
- **Caveat.** Out-of-sample LD for very large GWAS often stops SuSiE from finding credible sets (no GWAS set at 4 of the 9 one-reference pairs, RPS2 had 10). In-sample GWAS LD is not available.
- Output: `results/06-compare/gwas/coloc-susie-mage-chm13v2.0-vs-grch38-ffx.tsv` (per locus: abf and SuSiE PP3/PP4, credible-set counts and the top SNPs on each reference).

**SuSiE-coloc at the coloc loci** (2026-10-29; `scripts/06-compare/coloc-susie.py` + `.sbatch`; job 380902, 2 h 04 min, 8 GB; 380440 and 380490 failed, see the ledger).
- **Why.** coloc's ABF test assumes one causal variant per trait. When a region has two signals, it can call a colocalisation that is really two different variants, or miss one. SuSiE splits each trait's signal into credible sets (one per independent signal) using LD, and coloc then tests every eQTL set against every GWAS set (Wang et al. 2020; Zou et al. 2022; Wallace 2021).
- **How.**
  - Loci: the 32 gene × trait pairs with ABF PP4 ≥ 0.5 on either reference, ±500 kb around the top variants.
  - Variants: SNPs present in both references' eQTL results and in the GWAS, with MAF ≥ 1% in MAGE and in EUR.
  - LD: for the eQTL, from MAGE's own GRCh38 genotypes; for the GWAS (European), from the 525 unrelated 1KGP EUR individuals.
  - Run with R coloc 5.2.3 (`runsusie`, `coloc.susie`) and susieR 0.14, using coloc's default priors. The R step was checked on a simulated locus first (planted shared variant: PP4 0.996).
- **Controls agree on both references.** IRF5 and ulcerative colitis gives 0.98 / 0.98, TOP3B and hypothyroidism 0.86 / 0.86, and IRF5 and hypothyroidism 0.61 / 0.61. So SuSiE-coloc gives the same answer when the gene does not change.
- **The reference-dependent loci** (SuSiE PP4 T2T / patched GRCh38):

| Gene, trait | ABF call | SuSiE-coloc | Reading |
|---|---|---|---|
| IFITM3-AS1, hypothyroidism | GRCh38 only | 0.004 / 0.998 | **Holds, and is sharper.** On GRCh38 the eQTL credible set and the GWAS credible set share their lead, chr11:324163 C>T. On T2T the eQTL lead is chr11:324170 A>G, 7 bp away, and SuSiE calls the two signals distinct (PP3 0.996) |
| LOC124903630, atrial fibrillation | GRCh38 only | no eQTL set / 0.91 | Holds: on T2T there is no eQTL credible set in the region |
| SERF1B, hypothyroidism | GRCh38 only | no eQTL set / 0.86 | Holds, in the same way. The GWAS signal is weak |
| PPP4R3B-DT, hypothyroidism | suggestive (0.77 / 0.95) | 0.39 / 0.91 | New candidate, favouring GRCh38 |
| PTGES3, atrial fibrillation | T2T only (0.94 / 0.004) | 0.03 / < 10⁻⁵ | **Not supported.** With several signals allowed, the eQTL (T2T lead chr12:56639896) and the atrial fibrillation signal (chr12:56716008) are different variants on both references (PP3 0.87 / 1.00). The ABF call was a single-signal artefact; drop it as an H6 example |
| SCAF11, hypothyroidism | T2T only | ≈ 0 / ≈ 0 | Not supported: distinct signals on both references |
| POM121C, atrial fibrillation; SMG1P6, ulcerative colitis; RMI1, hypothyroidism; SCAF11, atrial fibrillation | one reference only | not testable | SuSiE found no GWAS credible set, or did not converge, with 1KGP EUR LD. That is common when reference-panel LD does not match the GWAS. The ABF result stands but is unconfirmed |

- **Takeaway.** The firmest H6 example is IFITM3-AS1 and hypothyroidism: a genome-wide GWAS signal (p 8×10⁻¹⁸) colocalises with the GRCh38 eQTL, while on T2T the eQTL moves to a variant 7 bp away that is not the GWAS variant. LOC124903630 and SERF1B are cases where the eQTL itself disappears on T2T. Allele-specific expression (in progress) should tell which reference's eQTL is right for these genes.
- Output: `results/06-compare/gwas/coloc-susie-mage-chm13v2.0-vs-grch38-ffx.tsv` (locus, ABF and SuSiE PP4/PP3, number of credible sets, lead variants of the best pair, per reference).

**lncRNAs** (2026-10-26; T2T vs patched GRCh38, STAR, matched variants; from `star-matched/genes.tsv`, `robust-matched/changed-chm13v2.0-vs-grch38-ffx.tsv` and the read-destination table; ad hoc on the login node, < 1 min).

*Share of lncRNAs in each gene list.* lncRNAs are 16% of eGenes but up to half of the reference-specific lists:

| Gene list | Genes (T2T / GRCh38) | lncRNAs |
|---|---:|---:|
| All eGenes | 15,013 / 15,224 | 2,382 / 2,736 (16% / 18%) |
| Robust changed set (unflagged) | 57 | 9 (16%) |
| eGene on one reference only, lead fails the replication rule (tested on both) | 33 / 28 | 11 / 12 (~37%) |
| Annotated on both, tested (expressed) on one only, unflagged | 101 / 76 | 49 / 41 (~50%) |
| Annotated on one reference only (eGenes) | 158 / 524 | 88 / 466 (56% / 89%) |

- One-reference-only, replication rule: T2T FAM27C, FAM27E3 (T2T-new sequence), ZBTB47-AS1, OR2A1-AS1, DDN-AS1; GRCh38 CCND2-AS1, BLK-AS1, LINC01410, LCMT1-AS1.
- Likely reasons: low expression puts lncRNAs near the expression filter and the eGene threshold (threshold effects); RefSeq predicts `LOC…` lncRNAs per assembly (annotation, not sequence; not counted for H0); a few sit in sequence T2T added or fixed.

*The 9 lncRNAs in the robust set* (slope T2T / GRCh38; CPM T2T / GRCh38; best read partner):

| lncRNA | Region | Stronger on | Slope | CPM | Partner (r) | Reading |
|---|---|---|---|---|---|---|
| IFITM3-AS1 | other | T2T | −0.45 / −0.12 | 1.5 / 4.7 | LOC124906937, T2T-only annotation (−0.96) | clear read move |
| LOC642846 | SD | T2T | −0.73 / −0.53 | 12.0 / 14.2 | DDX12P (−0.86) | paralog swap |
| GLIDR | SD | T2T | +0.59 / −0.34 | 0.6 / 2.0 | gene in T2T-new sequence (−0.66) | reads move to new sequence; sign flips |
| LOC124902471 | other | GRCh38 | −0.23 / −0.62 | 0.5 / 2.2 | TSPAN14-AS1 (−0.57) | weak partner |
| MIR570HG | other | GRCh38 | +0.04 / −0.45 | 0.2 / 1.5 | LOC124906327 (−0.53) | weak partner, very low expression |
| LOC105374445 | other | T2T | +1.37 / +0.72 | 1.0 / 2.6 | none (−0.50) | unclear |
| RAD51B-AS1 | other | T2T | +0.29 / −0.11 | 2.1 / 1.3 | none (−0.42) | unclear; sign flips |
| PPP4R3B-DT | other | T2T | −0.92 / −0.41 | 1.6 / 2.4 | none (−0.36) | unclear |
| LOC124906110 | other | GRCh38 | −0.03 / +0.29 | 42.2 / 2.9 | ANKRD17 (−0.88) | probably technical (correlates with many unrelated genes) |

- Usable examples: IFITM3-AS1, LOC642846, GLIDR. The broader lncRNA point ("the reference decides whether a lowly expressed lncRNA is detected") needs Geuvadis replication before it goes in the paper.

**Reversed eQTLs** (2026-10-26; robust set, STAR slopes at the compared variant; `scripts/06-compare/reversal-variant-location.py` via `reversal-variant-location.sbatch`, job 370516, 14 s, 0.7 GB; first attempt 370515 failed because `--wrap` runs `sh`).

*Which genes reverse.* 31 of the 57 genes have opposite-sign slopes, but only some are true reversals:

| Class (nominal p of each slope) | Genes |
|---|---:|
| Both strong (p < 1e-5 on both references): true reversal | 13 |
| Both nominal (p < 0.01 on both) | 6 |
| One side near zero: eQTL present on one reference only, not a reversal | 12 |

*The 13 true reversals* (slope T2T / patched GRCh38; z_diff STAR / salmon; where the compared variant lies, the same on both references):

| Gene | Type | Slope | z_diff | Read partner (r) | Variant location |
|---|---|---|---|---|---|
| STAG3L1 | pseudogene, SD | −0.58 / +0.63 | −23.7 / −18.6 | STAG3L2 (−0.70) | outside exons, 153 kb from the gene, in an SD |
| TCAF2P1 | pseudogene, SD | +0.78 / −0.45 | +22.4 / +49.8 | TCAF2 (−0.40) | **in an exon of its paralog TCAF2** |
| ASAH2 | protein-coding | +0.77 / −0.53 | +21.0 / +27.0 | ASAH2B (−0.91) | outside exons, 534 kb |
| ANAPC1P2 | pseudogene, SD | −0.47 / +0.36 | −20.4 / −21.9 | | outside exons, 134 kb, SD ≥ 98% |
| PMS2CL | pseudogene, SD | −0.24 / +0.63 | −17.9 / −19.2 | PMS2 (−0.80) | outside exons, 6 kb, in an SD |
| H3-3A | protein-coding | −0.14 / +0.32 | −17.8 / −23.9 | H3P6 (−0.71) | outside exons, 9 kb |
| RPS2 | protein-coding | −0.15 / +0.41 | −17.6 / −12.8 | | **in its own exon** |
| ELOB | protein-coding | −0.32 / +0.07 | −17.3 / −18.5 | | **in its own exon** (indel) |
| TCAF2 | protein-coding, SD | −0.33 / +0.27 | −16.9 / −18.1 | TCAF2P1 | **in its own exon** (same variant as TCAF2P1) |
| CCZ1 | protein-coding, SD | +0.65 / −0.25 | +14.7 / +14.6 | CCZ1B (−0.96) | outside exons, 894 kb, SD ≥ 98% |
| RASA4DP | pseudogene, SD | −0.33 / +0.74 | −13.5 / −9.5 | | in an intron, SD ≥ 98% |
| PPIP5K1P1-CATSPER2 | pseudogene, SD | +0.53 / −0.32 | +11.2 / +8.8 | PPIP5K1 (−0.53) | in an exon of STRC, 30 kb away, SD ≥ 98% |
| NUTM2E | protein-coding, SD | −0.38 / +0.43 | −8.9 / −6.1 | | in an intron, SD ≥ 98% |

*Reading.*
- **Variant in transcribed sequence (4 of 13: RPS2, ELOB, and the TCAF2/TCAF2P1 pair at one variant, chr7:145008621).** Candidate "mapping eQTLs": the allele itself can make a read match one copy better than another (a processed pseudogene for RPS2 and ELOB, the paralog for TCAF2). The sign then depends on which copies the reference contains. TCAF2/TCAF2P1 is the clearest case: one exonic variant, mirror-image effects on the two copies, flipping between references.
- **Variant outside exons (9 of 13), often far away (6–894 kb) and in SDs.** The allele cannot steer reads directly. More likely it tags a copy-number or structural difference between people (ctyper's model), which changes how many reads the copies receive; the references then split those reads differently. Supports adding ctyper copy number as a covariate.
- Why it matters: colocalisation, Mendelian randomisation and TWAS use the direction of an eQTL, so a reversal can flip whether a risk allele is read as raising or lowering a gene. PMS2 (Lynch syndrome) and its pseudogene PMS2CL are a known clinical testing problem.
- Not yet shown which reference is right. Allele-specific expression in heterozygotes would decide; it needs read alignments for these loci, which were not kept (only `keep-regions.bed` loci are).
- Output: `results/06-compare/eqtl-comparison/mage/robust-matched/reversal-variant-location-chm13v2.0-vs-grch38-ffx.tsv` (all 31 opposite-sign genes).

**Geuvadis replication** (2026-10-27; chain `submit-eqtl-cohort.sh geuvadis`, step log 5–6g; replication by `scripts/06-compare/replicate-changes.py`, login node, 34 s, 1.6 GB).

*Geuvadis on its own* (451 individuals, STAR / salmon, matched variants):
- eGenes, T2T / GRCh38 / patched GRCh38: 13,221 / 13,481 / 13,517 (STAR), 13,409 / 13,492 / 13,448 (salmon).
- The positive control repeats: patched vs GRCh38 changes 8 genes (STAR), all in FixItFelix regions (GATD3, PWP2, CBS, ICOSLG, U2AF1 and others).
- T2T vs patched GRCh38: 121 STAR, 295 salmon, 78 in both, 64 without artefact flag, none in patched regions.

*Same-variant replication of the MAGE set.* For each of the 57 MAGE genes, the MAGE STAR variant pair (var_a on T2T, var_b on patched GRCh38) is looked up in the Geuvadis nominal results under both methods. Alleles are aligned with the MAGE genotype correlation, and z is compared with the MAGE z. Variant IDs and alleles are the same in both cohorts, so the signs are directly comparable.

| T2T vs patched GRCh38 | Genes |
|---|---:|
| MAGE robust set (unflagged) | 57 |
| Testable in Geuvadis (both variants pass MAF 0.01, gene tested; untested: TCAF2P1, SPDYE2B, GLIDR) | 54 |
| Same sign, STAR | 53 (sign test p = 3e-15) |
| Same sign under STAR and salmon, p < 0.05 under both | **50** |
| Same, Bonferroni (p < 0.05/54) | **44** |

- Not replicated: UBC (salmon p = 0.15), PPP4R3B-DT and PTGES3 (STAR p ≈ 0.1; salmon strongly the same way), PDPR2P (no difference).
- Positive control (patched vs GRCh38, 6 unflagged MAGE genes): 5 replicate nominally and at Bonferroni; U2AF1 only under STAR (p = 0.012).
- T2T vs GRCh38: 55 of 60 testable replicate nominally, 49 at Bonferroni.
- **Reversed eQTLs replicate.**
  - 12 of the 13 true reversals are testable (TCAF2P1 is not), and 10 of those 12 again have opposite-sign slopes in Geuvadis: ASAH2, ANAPC1P2, H3-3A, STAG3L1, CCZ1, PMS2CL, RPS2, NUTM2E, PPIP5K1P1-CATSPER2 and RASA4DP.
  - In TCAF2 and ELOB the GRCh38 slope is about 0 in Geuvadis: still a difference, but one-sided.
- Table: `results/06-compare/eqtl-comparison/replication-geuvadis/{genes-<A>-vs-<B>.tsv, summary.tsv}`.

*Caveat found on the way: `stronger_on` and z signs from separately chosen leads are not comparable.*
- `compare-effects.py` signs z_diff by slope_A − slope_B in A's allele coding, while `stronger_on` compares |slope|. In 31 of the 57 genes the slopes have opposite signs, so the two disagree, and both depend on which lead variant was picked.
- Comparing `stronger_on` across cohorts gave a misleading 37/55. Hence the same-variant design above.
- Within a cohort, `combine-methods.py` requires the same z_diff sign under STAR and salmon even when their lead variants differ (18 of 57 in MAGE). Checked 2026-10-27: at the STAR lead pair, the MAGE salmon z has the same sign and Bonferroni significance for **all 57**, so the set stands.
- "Stronger on" in the tables and slides is only meaningful for same-sign genes. For opposite-sign genes, report the two slopes.

**Next steps.**
- For genes without a clear partner, look at multi-mapped reads in the BAMs and per-contig counts.
- Reversed eQTLs: allele-specific expression at the 13 loci on both references (needs realignment of a subset of samples with full BAMs).
- Add ctyper copy number as a covariate: 34 of the 57 genes (60%) are in ctyper's 1KGP copy-number table, against 13% of eGenes tested on both references (symbol match, so a minimum; checked 2026-10-26).
- Add chrX (first take sex of HG00237 and NA07000 from the 1KGP pedigree).
- Geuvadis replication done (above). H1 replicated 2026-10-28 (8 of 15 solid examples, 5 variants; above). Still to do: the lncRNA point and the read destination.
- ~~Make `combine-methods.py` and the figures use effects at the same variant pair under both methods, and label opposite-sign genes by their slopes rather than `stronger_on`.~~ Done 2026-10-27: `combine-methods.py` now adds `pattern` ("stronger on <ref>" when the STAR slopes share a sign, else "opposite sign") and salmon z at STAR's own variant pair. The selection rule and the gene sets are unchanged (checked against the previous tables). T2T vs patched GRCh38, unflagged: MAGE 13 stronger on T2T / 13 stronger on patched / 31 opposite sign, salmon agrees at STAR's pair for 57/57; Geuvadis 11 / 14 / 39, 64/64. Figure 8b and the workbook's 57-gene sheet use `pattern` (`make-figures.R`, `build-workbook.py`); the figures, workbook and deck are not yet regenerated.
- Run `compare-eqtl.py` as a Slurm job next time, not on the login node.

**Software** (pixi environment, 2026-10-17):
- tensorQTL 1.0.10 and its helper `qtl` 0.1.10 are vendored from their PyPI wheels (sha256 checked) in `envs/vendor/` and put on `PYTHONPATH` in `pixi.toml`. pixi's PyPI route fails here because its conda–PyPI mapping service (`conda-mapping.prefix.dev`) is reset from Karakoram.
- From conda-forge/bioconda: `pytorch-cpu` 2.14, `scipy`, `pgenlib`, `pandas-plink`, `pyarrow`, `seaborn`, `bx-python`, `pybigwig`, `plink2` 2.0a7.8, `bioconductor-qvalue`, `rpy2`.
- Tested: `import tensorqtl` works; q-values via R `qvalue` work.
- Compute: CPU only (the node's one GPU is used by the t2t-atlas Parabricks runs). Permutation mapping is split by chromosome into parallel jobs; cost to be measured on MAGE.

### 13B.4 Method: how a gene's eQTL effect is compared between two references (written 2026-10-29)

This section collects the comparison method in one place. 13B.3 gives the results and the history of each part.

**The question.** For one gene and two references A and B (for example T2T and patched GRCh38): does the same genetic variant have a different effect on the gene's expression when the reads are processed on A instead of B? The individuals, reads, genotype panel, software and settings are identical, so any real difference comes from the reference: where reads land, the gene models, and which variants exist.

**Step 1: map eQTLs separately on each reference** (stage 05, 13B.2).
- Expression is TMM-normalised CPM, then an inverse-normal transform per gene, as in the GTEx pipeline. Effect sizes (slopes) are therefore in standard deviations (SD) of expression per copy of the alternative allele.
- Covariates on each reference: 5 genotype PCs, sex and 60 expression PCs, the expression PCs computed on that reference's own expression.
- tensorQTL `map_cis` gives each gene's lead variant and q-value; `map_nominal` gives slope, se and p for every gene–variant pair within 1 Mb.

**Step 2: choose the variant pair to compare** (`compare-effects.py`).
- Only genes tested on both references and an eGene (q < 0.05) on at least one.
- Each reference's lead variant is translated to the other reference through the variant-matching table (13B.1). Between the two GRCh38 references the variant IDs are identical.
- The two variants' genotypes are correlated across the cohort (r). The pair is used only if |r| ≥ 0.9. B's slope is multiplied by sign(r), so that both slopes describe the same allele. Using the genotypes, rather than the matching table's "swapped" label, fixed wrong signs at some indels (13B.3).
- This gives two candidate pairs per gene (A's lead, B's lead). The one with the larger |z_diff| is kept. The selection is a mild multiple test; the internal null (below) includes it and shows no inflation.

**Step 3: the difference and its uncertainty.**
- Δ = slope_A − slope_B.
- *Old test, `compare-effects.py`:* z_diff = Δ / √(se_A² + se_B²). This is the standard test for two coefficients from **independent** samples (Clogg et al. 1995; Paternoster et al. 1998). Our two slopes are not independent, so it is far too conservative here. It produced the 57-gene set, which is therefore a set of very safe calls.
- *Paired test, `paired-test.py`:* z_paired = Δ / √(se_A² + se_B² − 2·cov(slope_A, slope_B)).
- **Why the covariance matters.** Both slopes come from the same individuals. When a person's expression is unusually high on A, it is almost always unusually high on B as well, because most of the reads are the same. That shared noise moves both slopes together, and it cancels in Δ. The median correlation of the two references' residual expression is ρ = 0.99. With equal se and correlation c between the slopes, the uncertainty of Δ is se·√(2(1 − c)): 1.41·se if independent, 0.14·se at c = 0.99.
- **How cov is estimated.** In each reference, expression (y) and genotype (g) are residualised on that reference's covariates plus an intercept, exactly as tensorQTL does; the recomputed slopes match tensorQTL to 5×10⁻⁴. Then slope = Σg·y / Σg². The errors e_A and e_B are correlated within a person and independent between people, so
  cov(slope_A, slope_B) = [Σ e_A e_B / dof] · Σ g_A g_B / (Σ g_A² · Σ g_B²), with dof = n − 2 − number of covariates.
  This is the covariance of two least-squares estimates whose errors are correlated, as in seemingly unrelated regressions (Zellner 1962) and Stata's `suest` (Weesie 1999). The same idea in genetics is the correlation of association statistics from overlapping samples (Lin & Sullivan 2009). The formula assumes the same error variance for every individual. The bootstrap below checks this assumption without making it.
- **When the gain is small.** The covariance is large only when both references measure the gene the same way. Where reads move between references (paralogs, false duplications), ρ is lower and the paired se approaches the independent one. Examples, T2T vs patched GRCh38, STAR:

| Gene | Slope T2T | Slope GRCh38-ffx | Δ | ρ | se of Δ, independent | se of Δ, paired | z_diff | z_paired |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ARMC10 | 0.004 | 0.440 | −0.436 | 0.93 | 0.062 | 0.017 | −7.0 | −25.1 |
| GCSH | 0.565 | 0.890 | −0.325 | 0.81 | 0.047 | 0.021 | −6.9 | −15.1 |
| SMN1 | 0.060 | 0.516 | −0.455 | 0.46 | 0.101 | 0.076 | −4.5 | −6.0 |
| NUTM2D | −0.559 | −0.022 | −0.537 | 0.40 | 0.074 | 0.057 | −7.3 | −9.4 |

**Step 4: calls.**
- *Changed in one quantification method:* two-sided p from z, Bonferroni over the genes compared (about 13,500; p < about 4×10⁻⁶).
- *Robust (the project rule):* changed under STAR **and** salmon, with the same sign of Δ. Requiring two independent quantifiers protects against artefacts of one aligner.
- *Excluded (artefact flag):* genes whose annotation differs between the two references (the gene models are not the same object), HLA/MHC genes, and PAR genes.
- *Size of the change (decided 2026-10-29):* |Δ| ≥ 0.1 SD under both STAR and salmon. The paired test detects changes as small as about 0.05 SD. Replication in the 358 independent Geuvadis individuals rises with the floor and then levels off, which is why the floor is 0.1:

| Definition (T2T vs patched GRCh38, unflagged) | Genes | Of the 57 | SD | Geuvadis tested | Same sign | Nominal | Bonferroni |
|---|---:|---:|---:|---:|---:|---:|---:|
| Paired, no floor (**extent**) | 255 | 57 | 68 | 247 | 221 | 175 (71%) | 142 |
| Paired, \|Δ\| ≥ 0.05 (STAR) | 171 | 57 | 66 | 163 | 153 | 137 (84%) | 116 |
| Paired, \|Δ\| ≥ 0.1 (STAR) | 124 | 57 | 62 | 117 | 112 | 104 (89%) | 90 |
| Paired, \|Δ\| ≥ 0.1 under both (**main**) | 117 | 57 | 59 | 110 | 106 | 100 (91%) | 86 |
| Paired, \|Δ\| ≥ 0.2 under both | 76 | 54 | 45 | 69 | 67 | 63 (91%) | 59 |
| Old test (**core**) | 57 | 57 | 33 | 54 | 54 | 53 (98%) | 51 |

- **Gene sets (decided 2026-10-29), `scripts/06-compare/define-gene-sets.py`:**
  - **main**: 117 genes, the headline result and the target list for allele-specific expression and GWAS;
  - **core**: the 57 genes of the old test, a high-confidence subset, all of them in main;
  - **extent**: 255 genes, the genome-wide extent.
  - Output: `results/06-compare/eqtl-comparison/mage/gene-sets/changed-<A>-vs-<B>.tsv`, with a `tier` column and the Geuvadis replication calls. T2T vs unpatched GRCh38: 63 / 121 / 246. Patch control: 6 / 7 / 8.

**Step 5: checks that the test is correct.**
- *Bootstrap over individuals* (Efron & Tibshirani 1993). For 300 random genes, the individuals are resampled 200 times, and the spread (SD) of Δ is measured directly, without any formula. It is 0.94–0.96 × the paired se, so the paired se is accurate and slightly conservative. The independent se is 9–15 × too large.
- *Genomic inflation λ* (Devlin & Roeder 1999): median z² / 0.455 over all genes; 1 means as expected under no change.
  - Independent test: λ = 0.01–0.03, so it misses almost everything.
  - Paired test, internal null (patched GRCh38 vs GRCh38 under salmon, which differ only in the patched regions): λ = 1.03.
  - Paired test, T2T comparisons: λ = 1.6–2.1. With the bootstrap agreeing, this reflects many small real shifts.
- *Copy number* (ctyper): adding paralog copy number as a covariate leaves 33 of 34 covered genes changed (13B.3).

**Step 6: replication** (`paired-test.py --at mage`). The discovery genes are tested in Geuvadis at exactly the MAGE variant pairs and allele coding, on the 358 individuals not in MAGE. A gene replicates if Δ has the MAGE sign under STAR and salmon, with p < 0.05 (nominal) or p < 0.05 / genes tested (Bonferroni).

**What is published and what is ours.**
- *Published, standard:* the eQTL mapping (tensorQTL, Taylor-Weiner et al. 2019; GTEx normalisation), the two-coefficient z-test, the covariance of correlated regressions (Zellner 1962; Weesie 1999), the bootstrap, λ and Bonferroni.
- *Our design:* applying these to the same individuals processed on two references. That means:
  - pairing variants across references through the matching table and the genotype correlation;
  - requiring STAR and salmon to agree;
  - flagging annotation and HLA/PAR artefacts;
  - the patched-vs-unpatched GRCh38 comparison as an internal null.
- FixItFelix (Behera et al. 2023) reported the eQTLs of its regions on both references side by side (Table S17) and described the ones that weaken, without a formal test of the difference. None of the reference comparisons read for this project (section 17) uses a paired test. This was not a systematic literature search; one should be done before claiming novelty.

## 14. History

### 14.1 How the project started (2026-10-13)

The project was created from a Claude Code session in `~/research/t2t-atlas`.

**Starting question.** The user (Yusuf Maulana, MBZUAI, CBCRG group) asked for recent studies on T2T and regulatory genomics, and what a study could achieve for a **conference submission**. Context: `t2t-atlas` is a pilot reprocessing ENCODE GM12878 ChIP-seq (20 experiments + 3 controls) on T2T-CHM13v2.0 and GRCh38. Pilot results (`t2t-atlas/docs/execution.md`, step 7):

- 98.1% of T2T narrow peaks are shared with GRCh38 no-alt.
- 74% of T2T-only peaks seen against GRCh38 *with alt contigs* were alt-contig artefacts (hotspots at MHC, 8p23.1, KIR).
- About 2,500 narrow peaks and 450 broad domains lie in new T2T sequence. Outside satellites: ×14 enriched, 50% reproduce in the other replicate. Inside satellite arrays: 18% reproduce, half from two libraries.
- Step 6 (`scripts/06-kmer-filter/`) is a unique-k-mer filter, still in progress in `t2t-atlas`.

**Paper options discussed for `t2t-atlas`:** (1) k-mer confidence score for peaks in new sequence + benchmark (best for RECOMB/ISMB); (2) "fair comparison" finding (alt-contig baseline); (3) resource/atlas (journal); (4) SD biology (needs distinction from the Jan 2026 preprint). A public GitHub repo `yusufdjuanda/encode-t2t-epigenome-reanalysis` (one CTCF experiment on T2T) appeared in the search; possibly the user's.

**How the direction was chosen.**

1. The user asked about **GWAS + T2T**. Finding: no study overlaps GWAS with raw-read T2T regulatory peaks; earlier work did variant calling (Aganezov 2022, Swedish cohort 2025) or liftover.
2. The user objected: **"it is just lifting over, so nothing will be important there."** Agreed: lifting GWAS SNPs adds no information, and a lifted SNP can never land in new sequence. With 98% of peaks shared, a GWAS-overlap study gives at most a few case loci.
3. Proposed instead: **do everything natively on T2T, no liftover.** 1KGP genotypes called on CHM13 + Geuvadis LCL RNA-seq remapped to T2T → T2T-native eQTLs → overlap with GM12878 T2T-only peaks → allele-specific binding in NA12878. The user agreed to try it.
4. The user asked for plain explanations ("I need to know and understand fully", "what does this mean actually"): what an eQTL is, why the reference affects genotypes, expression and gene models, the mechanisms, why both references are needed, and how ChIP peaks and allele-specific binding supply the mechanism. **Explain concepts plainly before acting.**
5. The study became its own project (moved from `t2t-atlas/scripts/07-eqtl`).
6. A second literature round strengthened the proposal: ctyper, MAGE, 1KGP ONT SVs and HGSVC3, Fiber-seq FIRE and HG001 ONT assemblies, the ENCODE4 cCRE registry, GWAS colocalisation (all now in sections 2 and 6). S-LDSC was judged not useful.

**Work done on 2026-10-13.** Novelty search (no T2T-native eQTL study found); Geuvadis metadata checked; genotype panel downloaded; subset job 356292 (step 1), whose result moved the weight from H1 to H2–H4.

**Incident (2026-10-13).** A `t2t-atlas` k-mer peak-calling array (job 355980) filled the node's root disk (`/tmp`, ~43 GB of MACS3 temp files). At the user's request it was cancelled, the user's temp files deleted, and `t2t-atlas/scripts/06-kmer-filter/call-peaks.sbatch` given a scratch `TMPDIR` and `--tempdir`. That array also asked for 576 GB in total, above the 400 GB limit; suggested `%24`. It still needs resubmitting in `t2t-atlas`.

### 14.2 2026-10-17

- All documentation merged into this file (previously `conversation-2026-10-13.md`, `study-proposal.md`, `execution.md`), at the user's request.
- Step 2 decisions, reference fixes and scripts (section 13). The 1000 Genomes FTP (`ftp.1000genomes.ebi.ac.uk`) also reset HTTPS connections; population labels came from the Geuvadis BioStudies sample sheet instead.
- Observed: the two `t2t-atlas` ChIP-seq pipelines (head jobs 359837, 359838) were using about 230 cores and 1,368 GB of allocated memory together, above the 128-core / 400 GB limit, which blocked the pilot's 64 GB jobs.
- ENA (`ftp.sra.ebi.ac.uk`) refused or reset most HTTPS and FTP connections when the pilot download opened 10 streams; it answered again a few minutes later. Downloads now use 2 streams. Later it allowed one connection at ~0.75 MB/s. NCBI SRA (`prefetch`, from AWS S3) gave ~12 MB/s with resets that `prefetch` resumes, so downloads moved to NCBI. Aspera would need ENA's key file, which the cluster module lacks.

- **Literature search (2026-10-17, sections 2.1–2.2):**
  - Confirmed the evidence for GRCh38 false duplications and collapses (Nurk 2022, Wagner 2022, Behera 2023) and what T2T fills in.
  - Found that FixItFelix (Behera et al. 2023) already compared GRCh38 with a patched GRCh38 for eQTLs in 449 Geuvadis LCLs, finding PWP2 and GATD3A artefacts.
  - Consequences:
    - framing changed to H0 ("beyond the known patches");
    - region classes added to step 6;
    - the go/no-go 2 criterion now counts only changes outside the patched regions;
    - a third reference (FixItFelix's patched GRCh38, step 3b) proposed;
    - candidate-gene sources corrected (Aganezov 2022 names no specific genes).
  - FixItFelix Additional file 2 (Tables S1–S19) saved to `data/external/fixitfelix-behera2023/` (4.2).

## 15. Open items

Rewritten 2026-10-29 (the earlier list was out of date; its done items are in the step log and 13B.3).

1. **Running:** ASE data, RNA realignment (380491, 369 runs) and ONT SNP calls (380492, 260 individuals); then write the ASE counting and test (13B.3, "Data for the ASE test"). Priority genes: IFITM3-AS1, ELOB, TCAF2, RPS2, PMS2CL, POM121C, SMG1P6, CCZ1, STAG3L1.
2. Read destination for the 60 further main genes (`read-destination.py` was run on the 57), and in Geuvadis.
3. chrX with reference-specific pseudo-autosomal boundaries (`--split-par`); take sex for HG00237 and NA07000 from the 1KGP pedigree first.
4. H5: link changed eQTLs to GM12878 regulatory elements from `t2t-atlas`.
5. The 24 core genes without a single read partner: multi-mapped reads in the alignments (the ASE realignment keeps target-region BAMs).
6. The higher STAR "unmapped: other" rate on T2T (0.93% vs 0.09%).
7. HPRC own-genome check for the SD reversals (assemblies resolve paralogs).
8. ctyper correction (doi 10.1038/s41588-026-02518-w, published 2026-02-16): content not yet checked.
9. Unused download `data/external/gwas-sumstats/GCST007517/` (targeted-array T2D study, replaced by GCST006867): delete only if the user agrees.
10. Before submission: a systematic literature search on paired tests of eQTL effects between references (13B.4), and checks of the remaining reference entries.
11. Commit only when the user asks. Since 2026-10-27 the repository is on GitHub (private): https://github.com/yusufdjuanda/t2t-eqtl; data, results, logs and work are symlinks there, so only code, config and docs are pushed.

## 16. Rule changes

| Date | Change | Reason |
|---|---|---|
| 2026-10-13 | Project created; step 1 moved from `t2t-atlas` (`scripts/07-eqtl` → `scripts/01-genotypes`, `config/eqtl` → `config/samples`, results `07-geuvadis-genotypes-chm13` → `01-…`, ledger row 356292) | Keep the eQTL study separate from the ChIP-seq atlas |
| 2026-10-17 | All documentation in one file, `docs/t2t-eqtl.md` | User request |
| 2026-10-17 | Gene models: NCBI RefSeq RS_2025_08 on both references instead of JHU Liftoff v5.3 | Same annotation pipeline on both; avoids projecting GRCh38 genes onto T2T |
| 2026-10-17 | Derived references allowed in `data/star-references/` (T2T with chrY PARs masked, chrM and chrEBV added) | Remove assembly differences unrelated to completeness (13.2) |
| 2026-10-17 | Project pixi environment (`pixi.toml` in the project root) | Python scripts and gffread |
| 2026-10-17 | Main eQTL comparison on matched variants only; tensorQTL and qtl vendored in `envs/vendor/` | The T2T and GRCh38 panels differ in filtering (13B.1); pixi cannot reach PyPI's mapping service |
| 2026-10-17 | Third reference (FixItFelix patched GRCh38) in step 3; Nextflow pipeline in `scripts/03-expression/nextflow/` with CSV samplesheets in `config/samplesheets/` | User decisions; tidy, resumable pipeline |
| 2026-10-17 | Central claim reframed as H0 ("beyond the regions FixItFelix patched"); region classes in step 6; go/no-go 2 counts changes outside the patched regions | FixItFelix (Behera et al. 2023) already published the chr21 part (2.2) |
| 2026-10-17 | MAGE is the discovery cohort, Geuvadis the replication | MAGE pilot (13.8): larger, 26 populations, one lab, 150 bp reads, effect replicates |
| 2026-10-17 | Reads come from NCBI SRA (`prefetch` + `fasterq-dump`), verified with `vdb-validate` and ENA's read-pair count, instead of ENA FASTQs checked by md5 | ENA gave ~0.75 MB/s from Karakoram (≈ 40 days for 2.5 TB); NCBI ~12 MB/s. Same reads; FASTQ headers differ, which does not affect alignment |
| 2026-10-27 | One T2T genotype folder: `results/01b-expression-individuals-genotypes-chm13/geuvadis/` renamed to `results/01-expression-individuals-genotypes-chm13/expression/` (set label in `summary.tsv` and `counts/` changed from `geuvadis` to `expression`); the superseded 451-person step-1 output `results/01-geuvadis-genotypes-chm13/` retired and 4.4 rewritten with the 1,089-person counts (final version only, user request); the per-chromosome `counts/` intermediate no longer kept, here and in `results/04-expression-individuals-genotypes-grch38/` (both subset scripts now write it to TMPDIR); `subset-genotypes.sbatch` defaults to the 1,089 list; `prepare-genotypes.sbatch`, `match-variants.sbatch`, `make-simulated-cohort.py` updated | The 1,089 set (731 MAGE + 451 Geuvadis − 93 in both) contains all 451; the old names said 'Geuvadis' for data that is both cohorts, and matched neither the GRCh38 folder (`04-…/expression/`) (user request) |
| 2026-10-27 | Interim cohort `mage420` removed (`results/05-eqtl/mage420/` 2.0 GB, `results/03-expression/mage420/` links and matrices, `config/samplesheets/mage420.csv`) | Superseded by the full MAGE cohort; nothing used it (user request) |
| 2026-10-27 | One sample table: `make-sample-lists.sh` now checks both cohorts (MAGE and Geuvadis ENA run reports) against the T2T panel and writes `config/samples/expression-individuals.tsv` (sample, mage, geuvadis, t2t_panel, included); it replaces `geuvadis-in-1kgp-chm13.txt`, `geuvadis-not-in-1kgp-chm13.txt` and `expression-individuals-in-1kgp-chm13.txt` (the last was written by `make-cohort-sheets.py`). `make-cohort-sheets.py` and both genotype subset scripts read the table; samplesheets and the default 1,089-person list are byte-identical | MAGE had been checked against the panel only by hand, and the 1,089 list was made in step 3 but used in step 1 (user request) |
| 2026-10-27 | Steps 1 and 4 merged into one genotype step: `scripts/01-genotypes/{make-sample-lists.sh, chm13v2.0/{fetch-panel.sh, subset-genotypes.sbatch}, grch38/{fetch-panel.sbatch, subset-genotypes.sbatch}}` (from `01-genotypes/fetch-1kgp-chm13-panel.sh`, `01-genotypes/subset-genotypes.sbatch`, `04-genotypes-grch38/fetch-1kgp-grch38-panel.sbatch`, `04-genotypes-grch38/subset-genotypes-grch38.sbatch`); results `results/01-genotypes/{chm13v2.0,grch38}/` (from `results/01-expression-individuals-genotypes-chm13/`, `results/04-expression-individuals-genotypes-grch38/`); logs `logs/01-genotypes/{chm13v2.0,grch38}/`. Files moved, nothing rerun; paths in scripts and this file updated (rows above keep the names of their time; the ledger keeps its step labels) | The two steps do the same thing on two references; mirrors the per-reference layout elsewhere (user request) |
| 2026-10-27 | References get their own stage; folders renumbered so scripts, results and logs share one number per stage: `scripts/02-references/` (`prepare-star-references.sbatch` from 02-rnaseq-pilot; `build-star-index.sbatch`, `prepare-ffx-reference.sbatch`, `make-region-beds.py` from 03-expression); `02-rnaseq-pilot` → `03-rnaseq-pilot` and `03-expression` → `04-expression` in `scripts/`, `results/`, `logs/` and `work/` (reference job logs to `logs/02-references/`); 05–07 unchanged. The 76-bp STAR index folder `star-2.7.11b/` (with alias `star-2.7.11b-oh75`) renamed `star-2.7.11b-oh75/` on chm13v2.0 and grch38, as on grch38-ffx; `CHECKSUMS.md5` updated. `results/03-rnaseq-pilot/comparison/` → `comparison-geuvadis/` (like `comparison-mage/`). Samplesheet `sra` paths now `work/04-expression/sra/`. Files moved, nothing rerun; Nextflow `-resume` caches in the old work folder no longer apply (all runs complete). Sections 5, 12, 13–13B keep their step numbers; 11.4 maps them | `data/star-references/` is used by every later stage, not only the pilot; one number per stage across folders (user request) |
| 2026-11-17 | (branch `capstone`) Software from pixi only (`bio` environment replaces cluster modules); no absolute paths in scripts (`config/site.env`, `scripts/lib/common.sh`, `scripts/submit.sh`); `data/`, `results/`, `logs/`, `work/` no longer tracked; samplesheet paths relative | Reproducible outside Karakoram; course CBIO7501 practices (11.5) |

## 17. Sources

- Aganezov et al. 2022, *A complete reference genome improves analysis of human genetic variation*, Science, doi:10.1126/science.abl3533
- Rhie et al. 2023, *The complete sequence of a human Y chromosome*, Nature, doi:10.1038/s41586-023-06457-y
- Lalli et al., Phased T2T 1KGP panel v1.1, Zenodo doi:10.5281/zenodo.7612953; bioRxiv doi:10.1101/2025.02.24.639687
- Gershman et al. 2022, *Epigenetic patterns in a complete human genome*, Science, doi:10.1126/science.abj5089
- T2T-CHM13 in the Swedish population, Genome Research 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC12581843/
- T2T-CHM13 reduces mapping bias at disease-associated variants, bioRxiv 2025, https://www.biorxiv.org/content/10.64898/2025.12.17.694618v2
- Widespread Regulatory Turnover Across Human Segmental Duplications, bioRxiv 2026, https://www.biorxiv.org/content/10.64898/2026.01.21.700947v1
- Nurk S et al. 2022, *The complete sequence of a human genome*, Science 376:44–53, doi:10.1126/science.abj6987 (PMC9186530)
- Wagner J et al. 2022, *Curated variation benchmarks for challenging medically relevant autosomal genes*, Nat Biotechnol 40:672–680, doi:10.1038/s41587-021-01158-1
- Vollger MR et al. 2022, *Segmental duplications and their variation in a complete human genome*, Science 376:eabj6965, doi:10.1126/science.abj6965 (T2T SD annotation; FixItFelix used it to build its decoys)
- Behera S et al. 2023, *FixItFelix: improving genomic analysis by fixing reference errors*, Genome Biol 24:31, doi:10.1186/s13059-023-02863-7 (PMC9942314; Additional file 2 Tables S1, S2, S16, S17 used)
- Shew CJ, Kaya G, McGinty SP, Dennis MY 2026, *Influence of cis-regulatory elements on expression divergence in human segmental duplications*, Genome Res 36:1318–1328, doi:10.1101/gr.281526.125
- T2T-CHM13 vs GRCh38 immunoglobulin isotypes in scRNA-seq, bioRxiv doi:10.1101/2023.05.24.542206
- Long-read RNA-seq in whole blood, GRCh38 vs T2T-CHM13, PMC12469794
- Dark genomic regions persist in CHM13, bioRxiv doi:10.1101/2025.05.23.655776
- Head et al. 2025, quantification method affects eQTL replicability, PubMed 40894657
- Lappalainen et al. 2013, *Transcriptome and genome sequencing uncovers functional variation in humans* (Geuvadis), Nature 501:506–511, doi:10.1038/nature12531. RNA-seq: ENA PRJEB3366; sample sheet BioStudies E-GEUV-1
- Byrska-Bishop et al. 2022, *High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios*, Cell, doi:10.1016/j.cell.2022.08.004 (the NYGC GRCh38 calls for step 4)
- Taylor-Weiner A et al. 2019, *Scaling computational genomics to millions of individuals with GPUs* (tensorQTL), Genome Biol 20:228, doi:10.1186/s13059-019-1836-7
- Zellner A 1962, *An efficient method of estimating seemingly unrelated regressions and tests for aggregation bias*, J Am Stat Assoc 57:348–368, doi:10.1080/01621459.1962.10480664
- Weesie J 1999, *Seemingly unrelated estimation and the cluster-adjusted sandwich estimator*, Stata Technical Bulletin 52:34–47 (sg121)
- Clogg CC, Petkova E, Haritou A 1995, *Statistical methods for comparing regression coefficients between models*, Am J Sociol 100:1261–1293, doi:10.1086/230638
- Paternoster R, Brame R, Mazerolle P, Piquero A 1998, *Using the correct statistical test for the equality of regression coefficients*, Criminology 36:859–866, doi:10.1111/j.1745-9125.1998.tb01268.x
- Lin DY, Sullivan PF 2009, *Meta-analysis of genome-wide association studies with overlapping subjects*, Am J Hum Genet 85:862–872, doi:10.1016/j.ajhg.2009.11.001
- Devlin B, Roeder K 1999, *Genomic control for association studies*, Biometrics 55:997–1004, doi:10.1111/j.0006-341x.1999.00997.x
- Efron B, Tibshirani RJ 1993, *An Introduction to the Bootstrap*, Chapman & Hall
- Wallace C 2021, *A more accurate method for colocalisation analysis allowing for multiple causal variants*, PLoS Genet 17:e1009440, doi:10.1371/journal.pgen.1009440
- Wang G et al. 2020, *A simple new approach to variable selection in regression, with application to genetic fine mapping*, J R Stat Soc B 82:1273–1300, doi:10.1111/rssb.12388
- Zou Y et al. 2022, *Fine-mapping from summary data with the "Sum of Single Effects" model*, PLoS Genet 18:e1010299, doi:10.1371/journal.pgen.1010299
- Wang G, Sarkar A, Carbonetto P, Stephens M 2020, *A simple new approach to variable selection in regression, with application to genetic fine mapping* (SuSiE), J R Stat Soc B 82:1273–1300, doi:10.1111/rssb.12388
- Wallace C 2021, *A more accurate method for colocalisation analysis allowing for multiple causal variants*, PLoS Genet 17:e1009440, doi:10.1371/journal.pgen.1009440
- Giambartolomei C et al. 2014, *Bayesian test for colocalisation between pairs of genetic association studies using summary statistics*, PLoS Genet 10:e1004383, doi:10.1371/journal.pgen.1004383
- Sollis E et al. 2023, *The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource*, Nucleic Acids Res 51:D977–D985, doi:10.1093/nar/gkac1010 (associations release 'latest' fetched 2026-10-29)
- Ma & Chaisson 2025, *Genotyping sequence-resolved copy number variation using pangenomes reveals paralog-specific global diversity and expression divergence of duplicated genes*, Nature Genetics; preprint doi:10.1101/2024.08.11.607269
- Taylor et al. 2024, *Sources of gene expression variation in a globally diverse human cohort* (MAGE), Nature, https://www.nature.com/articles/s41586-024-07708-2; https://github.com/mccoy-lab/MAGE
- Schloissnig et al. 2025, *Structural variation in 1,019 diverse humans based on long-read sequencing*, Nature, https://www.nature.com/articles/s41586-025-09290-7
- HGSVC3 and 1KGP long-read releases, https://www.internationalgenome.org/announcements/New-publications-HGSVC-and-IKG/
- *A haplotype-resolved view of human gene regulation* (Fiber-seq FIRE, GM12878), bioRxiv doi:10.1101/2024.06.14.599122; PubMed 40501892
- *Efficient near-telomere-to-telomere assembly of nanopore simplex reads* (hifiasm-ONT, incl. HG001), Nature 2026, https://www.nature.com/articles/s41586-026-10105-6
- Moore et al. 2026, *An expanded registry of candidate cis-regulatory elements*, Nature; preprint doi:10.1101/2024.12.26.629296
- ChromBPNet, RRID SCR_024806
- Population-level SV characterisation using pangenome graphs (Swave), Nature Genetics 2026, https://www.nature.com/articles/s41588-026-02538-6
- NCBI RefSeq annotation releases GCF_009914755.1-RS_2025_08 (T2T-CHM13v2.0) and GCF_000001405.40-RS_2025_08 (GRCh38.p14)
- UCSC CHM13v2.0 analysis set with masked chrY PARs (`chm13v2.0_maskedY`), https://github.com/marbl/CHM13
