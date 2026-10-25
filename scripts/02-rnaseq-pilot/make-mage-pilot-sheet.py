"""Step 2b: write the MAGE pilot sample sheets.

    pixi run python scripts/02-rnaseq-pilot/make-mage-pilot-sheet.py

Writes, in the same columns as config/samples/pilot-5.tsv (so fetch-sra.sbatch and align-quantify.sbatch work
unchanged; URL and md5 columns are NA because reads come from NCBI SRA):
  config/samples/mage-pilot-5.tsv          five MAGE runs
  config/samples/geuvadis-overlap-1.tsv    the Geuvadis run of HG00265, who is in both cohorts
Inputs: data/metadata/mage/ena-PRJNA851328-2026-10-17.tsv (ENA run report for MAGE; sample_title holds
"individual,population,super-population,sex") and the Geuvadis files used by make-pilot-sheet.py.
Selection (2026-10-17): single-run individuals near MAGE's median depth (28.4 M pairs), excluding the known sample
swap SRR19762530; two African (ESN, ASW), one East Asian (CHS), one South Asian (STU) and one British individual
who is also in Geuvadis (HG00265); three male, two female. None of the five Geuvadis pilot individuals is in MAGE.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
MAGE_RUNS = ["SRR19762269", "SRR19762285", "SRR19762887", "SRR19762232", "SRR19762408"]
GEUVADIS_OVERLAP_RUNS = ["ERR188222"]  # HG00265 in Geuvadis
COLS = ["run", "sample", "ancestry", "sex", "lab", "read_pairs", "url_1", "url_2", "md5_1", "md5_2"]

m = pd.read_csv(ROOT / "data/metadata/mage/ena-PRJNA851328-2026-10-17.tsv", sep="\t")
t = m.sample_title.str.split(",", expand=True)
mage = pd.DataFrame({"run": m.run_accession, "sample": m.sample_alias, "ancestry": t[1] + " (" + t[2] + ")",
                     "sex": t[3], "lab": "MAGE", "read_pairs": m.read_count,
                     "url_1": "NA", "url_2": "NA", "md5_1": "NA", "md5_2": "NA"})
mage = mage.set_index("run").loc[MAGE_RUNS].reset_index()[COLS]
mage.to_csv(ROOT / "config/samples/mage-pilot-5.tsv", sep="\t", index=False)

s = pd.read_csv(ROOT / "data/metadata/geuvadis/E-GEUV-1.sdrf.txt", sep="\t")
s = s[["Comment[ENA_RUN]", "Characteristics[sex]", "Characteristics[ancestry category]",
       "Characteristics[laboratory]"]].drop_duplicates()
s.columns = ["run", "sex", "ancestry", "lab"]
e = pd.read_csv(ROOT / "data/metadata/geuvadis/ena-PRJEB3366-fastq-md5-2026-10-17.tsv", sep="\t")
e = e.rename(columns={"run_accession": "run"}).merge(s, on="run")
url = e.fastq_ftp.str.split(";", expand=True)
md5 = e.fastq_md5.str.split(";", expand=True)
g = pd.DataFrame({"run": e.run, "sample": e.sample_alias.str.replace("GEUV:", ""), "ancestry": e.ancestry,
                  "sex": e.sex, "lab": e.lab.astype(int), "read_pairs": e.read_count,
                  "url_1": "https://" + url[0], "url_2": "https://" + url[1], "md5_1": md5[0], "md5_2": md5[1]})
g = g.set_index("run").loc[GEUVADIS_OVERLAP_RUNS].reset_index()[COLS]
g.to_csv(ROOT / "config/samples/geuvadis-overlap-1.tsv", sep="\t", index=False)

print(pd.concat([mage, g])[COLS[:6]].to_string(index=False))
