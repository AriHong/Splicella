# Splicella

Splicella is a Python framework providing standardized workflows for isoform diverisity analysis using single cell long read seuqnecing data. Two complementary analyses are availiable: (1) DCI (Distinct Combination of Isoforms) analysis identifying genes with significantly different isoform combinations between cell groups, and (2) PSI (Percent Spliced In) analysis for exon usage pattern detection across cell populations.

---

## Installation
Splicella can be easily installed from directly from the source.


```bash
git clone https://github.com/AriHong/Splicella.git
cd Splicella

conda create -n splicella python>=3.9.0
conda activate splicella
pip install -r requirements.txt
```

---

## Quick Start

We provide a demonstration notebook under `notebooks/Demo.ipynb` showing how to run MoFlow on a toy dataset and compute downstream scores.

---



