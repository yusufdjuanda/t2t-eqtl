# t2t-eqtl

T2T-native eQTL mapping in lymphoblastoid cells: MAGE (731 individuals, discovery) and Geuvadis (451, replication), quantified on T2T-CHM13v2.0, GRCh38 and FixItFelix-patched GRCh38 from the same raw reads, and linked to GM12878 regulatory elements from the sibling project `~/research/t2t-atlas`.

All documentation (aim, plan, rules, how the code fits together, step log, history): [docs/t2t-eqtl.md](docs/t2t-eqtl.md)

Status (2026-10-27): MAGE expression, cis-eQTLs (chr1–22) and the reference comparison are done. Main claim (H0): what a complete reference changes beyond the GRCh38 errors FixItFelix (Behera et al. 2023) already patched; robust set of 57 genes changed by T2T vs patched GRCh38, none in the patched regions. Geuvadis expression is running.
