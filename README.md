# The reading frame in dentinogenesis imperfecta and the stop position in amelogenesis imperfecta decide which truncating variants cause disease

**N. Vindy Valentina<sup>1</sup>, Alper Karagöl<sup>2</sup>, Taner Karagöl<sup>3</sup>**

<sup>1</sup>Independent dental practice  
<sup>2</sup>Department of Mechanical Engineering, Massachusetts Institute of Technology, Cambridge, MA, USA  
<sup>3</sup>McGovern Institute for Brain Research, Massachusetts Institute of Technology, Cambridge, MA, USA

## Abstract

Truncating variants are interpreted under a loss-of-function model in which a premature stop removes the product and position matters only through escape from nonsense-mediated decay. We tested that model across the 31 genes that ClinVar links to an isolated dental condition. The pathogenic set is 905 truncating variants with the reported phenotype taken per submission, and the population control is 1,911 predicted loss-of-function variants from gnomAD v4 in the same genes. Position alone separates the phenotypes in one gene. In LAMB3 the variants reported with amelogenesis imperfecta fall in the escape window while those reported with junctional epidermolysis bullosa do not, P < 0.001, and no other gene separates under any of four window definitions. The reading frame separates them where position does not. In DSPP, 50 of 52 pathogenic truncating variants are frameshifts, a share of 0.962 against 0.519 among the gnomAD variants in the same gene (q < 0.001), so the excess is not what mutation supplies. Those frameshifts translate a median of 372.5 novel amino acids before the new stop against 9 for the population set, P < 0.001. The panel-wide form of that comparison is carried by DSPP alone and falls to P = 0.059 when it is removed. Across 1,719 genes with twenty or more pathogenic truncating variants in ClinVar, DSPP holds the second highest frameshift share and FAM83H the lowest, so both extremes of the genome-wide distribution fall on this panel. We give a per-gene table of which evidence applies where, for use when PVS1 is applied to a dental gene.

**Keywords:** Truncating variants, Nonsense-mediated decay, Frameshift, PVS1, Amelogenesis imperfecta, Dentinogenesis imperfecta

## Repository layout

```
dental-truncation-frame/
├── README.md
├── Supplementary_Information.pdf
├── code/
│   ├── 01_panel_v2.py                   Panel definition from the ClinVar gene-condition map
│   ├── 02_gnomad_pull.py                gnomAD v4 predicted loss-of-function variants and MANE exons
│   ├── 03_labels_submission_level.py    Reported phenotype per submission
│   ├── 04_review_filter.py              Per-gene counts by review status
│   ├── 05_frame_dependence.py           Frameshift against nonsense composition
│   ├── 06_fsd2.py                       Frameshift-to-stop distance on the full transcript
│   ├── 07_export_variant_table.py       One row per variant, the table the figures read
│   ├── 08_fsd_summary.py                Neo-tail summary statistics
│   ├── 09_genomewide_power.py           Frameshift share across ClinVar
│   └── 10_make_figures.py               Rebuilds every figure from the deposited tables
├── dataset/
│   ├── variant_table.csv                Per-variant table, both sets
│   ├── dental_panel_v2.csv              The panel as a flat table
│   └── dspp_features.json               UniProt domain boundaries for DSPP
├── analysis/
│   ├── derived_v2.json                  Panel and per-gene ClinVar counts
│   ├── derived_phenotypes.json          Genes carrying a second condition
│   ├── derived_labels_submission.json   Window counts by reported phenotype
│   ├── derived_review.json              Counts by review status
│   ├── derived_frame.json               Composition tests
│   ├── derived_fsd2.json                Neo-tail lengths
│   ├── derived_fsd_exclusion.json       Panel-wide neo-tail comparison
│   └── derived_genomewide_power.json    Genome-wide frameshift shares
└── figures/                             Figures as rebuilt by code/10_make_figures.py
```

## Dataset summary

| File | Rows | Description |
| --- | --- | --- |
| `dataset/variant_table.csv` | 2,822 | One row per variant: 905 pathogenic truncating variants from ClinVar and 1,917 predicted loss-of-function variants from gnomAD v4, of which 1,911 carry a codon position readable from the protein-level notation, with gene, consequence, codon position, window membership and neo-tail length. |
| `dataset/dental_panel_v2.csv` | 50 | The genes ClinVar associates with a dental condition, 31 of them carrying an isolated dental condition. |
| `analysis/derived_*.json` | - | Derived results each figure and table reads, one file per analysis step. |

## Reproducing the analysis

Every figure rebuilds from the deposited tables alone:

```
python3 code/10_make_figures.py
```

Rebuilding the derived files from source needs two downloads that are not redistributed here. Place the ClinVar `variant_summary.txt.gz` and `submission_summary.txt.gz` for the release you want in the repository root, then run, from the repository root:

```
python3 code/01_panel_v2.py
python3 code/02_gnomad_pull.py
python3 code/03_labels_submission_level.py
python3 code/04_review_filter.py
python3 code/05_frame_dependence.py
python3 code/06_fsd2.py
python3 code/07_export_variant_table.py
python3 code/08_fsd_summary.py
python3 code/09_genomewide_power.py
```

`code/02_gnomad_pull.py` queries the gnomAD GraphQL API and writes `analysis/gnomad_plof.json`. `code/06_fsd2.py` builds a transcript cache on first run. Results depend on the ClinVar release used, so the release the manuscript reports is stated in its Methods.

## Software dependencies

| Domain | Packages |
| --- | --- |
| Numerics | `numpy` |
| Statistics | `scipy` |
| Plotting | `matplotlib` |
| External APIs | `urllib` (standard library), gnomAD GraphQL, UniProt REST |

Tested on Python 3.10.12 with numpy 2.2.6, scipy 1.15.3 and matplotlib 3.10.9.

## Data availability

* ClinVar variant and submission summaries are available from NCBI: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/`. The classifications and submission records are cited and not redistributed here.
* gnomAD v4 variant records are queried through the gnomAD API: `https://gnomad.broadinstitute.org/`. They are not redistributed here.
* UniProt and RefSeq records are queried through their REST interfaces: `https://rest.uniprot.org/` and `https://eutils.ncbi.nlm.nih.gov/`.
* The per-variant table and the derived result files needed to rebuild every figure are included under `analysis/`.

## Citation

Valentina N.V., Karagöl A., Karagöl T. *The reading frame in dentinogenesis imperfecta and the stop position in amelogenesis imperfecta decide which truncating variants cause disease* (2026).

## Contact

Taner Karagöl, karagolt@mit.edu

## License

The code in this repository is released under the MIT License. The underlying ClinVar, gnomAD, UniProt and RefSeq data are subject to the terms of those providers.
