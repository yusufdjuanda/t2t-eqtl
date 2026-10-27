"""Step 2: write the pilot sample sheet config/samples/pilot-5.tsv.

    pixi run python scripts/02-rnaseq-pilot/make-pilot-sheet.py

Inputs (data/metadata/geuvadis/):
  ena-PRJEB3366-fastq-md5-2026-10-17.tsv  ENA run report: run, sample alias, FASTQ URLs and md5s, read pairs
  E-GEUV-1.sdrf.txt                       Geuvadis sample sheet (BioStudies E-GEUV-1): sex, ancestry, lab per run
The five runs were chosen on 2026-10-17 from Geuvadis individuals that (a) have one run only, (b) are in
config/samples/expression-individuals.tsv (Geuvadis, included), and (c) are near the median depth (29.3 M pairs): two Yoruba, one each
Utah (CEU), British and Finnish; three male, two female; four from lab 1. The run IDs are listed here so the sheet
is reproducible; this script only joins their metadata.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RUNS = ["ERR188298", "ERR188139", "ERR188107", "ERR188242", "ERR188197"]

s = pd.read_csv(ROOT / "data/metadata/geuvadis/E-GEUV-1.sdrf.txt", sep="\t")
s = s[["Comment[ENA_RUN]", "Characteristics[sex]", "Characteristics[ancestry category]",
       "Characteristics[laboratory]"]].drop_duplicates()
s.columns = ["run_accession", "sex", "ancestry", "lab"]
e = pd.read_csv(ROOT / "data/metadata/geuvadis/ena-PRJEB3366-fastq-md5-2026-10-17.tsv", sep="\t")
m = e[e.run_accession.isin(RUNS)].merge(s, on="run_accession")
url = m.fastq_ftp.str.split(";", expand=True)
md5 = m.fastq_md5.str.split(";", expand=True)
out = pd.DataFrame({"run": m.run_accession, "sample": m.sample_alias.str.replace("GEUV:", ""),
                    "ancestry": m.ancestry, "sex": m.sex, "lab": m.lab.astype(int), "read_pairs": m.read_count,
                    "url_1": "https://" + url[0], "url_2": "https://" + url[1], "md5_1": md5[0], "md5_2": md5[1]})
out = out.set_index("run").loc[RUNS].reset_index()
out.to_csv(ROOT / "config/samples/pilot-5.tsv", sep="\t", index=False)
print(out.drop(columns=["url_1", "url_2", "md5_1", "md5_2"]).to_string(index=False))
