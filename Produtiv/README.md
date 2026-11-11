![PWS overview](https://github.com/user-attachments/assets/753aa507-2f9e-45c6-8dea-510f91ee2f2f)

# PWS: Parameter-Weighted Score

A Python tool to rank bacterial taxa by how their change in relative abundance between two conditions relates to an external numerical parameter measured in each experimental group (for example, crop yield, pH or temperature).

PWS was developed and first applied in the study below, where it ranked soil bacterial species by their yield-weighted change in abundance before planting and after harvest of garlic in four commercial farms in southern Brazil:

> Azevedo GPR, Fritsche Y, Giovanetti LK, da Silva D, Marques CEC, Mazzuco V, Welter LJ, Stefenon VM (2026). Changes in Soil Microbiome Linked to Management of *Allium sativum* Cultivation: Impacts on Soil Health and Agricultural Productivity. *Journal of Soil Science and Plant Nutrition*. https://doi.org/10.1007/s42729-026-03680-0

**On the name.** In the article the score is called *Productivity*-Weighted Score, because the parameter used was garlic yield. The tool accepts any numerical parameter measured on the same scale across groups, hence *Parameter*-Weighted Score. The calculation is identical.

**Scope.** PWS is a ranking and screening approach. It points to candidate taxa from field data so they can be evaluated in controlled, multi-site studies. It does not establish causal effects of any taxon on the parameter.

---

## Contents

- [How the score works](#how-the-score-works)
- [Installation](#installation)
- [Input data](#input-data)
- [Usage](#usage)
- [Reproducing the published analysis](#reproducing-the-published-analysis)
- [Output files](#output-files)
- [Statistical testing and filtering](#statistical-testing-and-filtering)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## How the score works

For each taxon *j* in group *i*, the relative change between the mean abundance after (*After_ij*) and before (*Before_ij*) is:

```
Δ_ij = (After_ij − Before_ij) / Before_ij
```

Each group's parameter value *P_i* is normalized by the mean across groups:

```
φ_i = P_i / mean(P)
```

The score is the mean of the group contributions (Δ_ij × φ_i), weighted by φ_i:

```
PWS_j = Σ_i ( Δ_ij × φ_i² ) / Σ_i φ_i
```

The sums run over the groups in which the taxon was detected in at least one condition. Groups with above-average parameter values therefore weigh more in the score.

| Situation in a group | Handling |
|---|---|
| Detected only after | Δ = 1 (the ratio is undefined) |
| Absent before and after | Group excluded from the sums |
| Missing values | Treated as 0 |

**Worked example (from the article).** *Flavobacterium buctense* had Δ = +4.089, +13.959, −0.944 and +0.163 in farms 3 to 6, with φ = 0.444, 1.206, 0.952 and 1.397:

```
PWS = (4.089×0.444² + 13.959×1.206² − 0.944×0.952² + 0.163×1.397²) / (0.444 + 1.206 + 0.952 + 1.397)
    = 20.583 / 4.000 = +5.15
```

**Interpretation.** A positive PWS means a net increase in relative abundance, with the increases in higher-parameter groups weighing more. A negative PWS means a net decrease of the same kind. Because Δ is a ratio, a large proportional change from a very low starting abundance can dominate a taxon's score. Check the per-group values (`{group}_Variation_Rel`, `{group}_Before_Mean`, `{group}_After_Mean` in `parameter_weighted_analysis.csv`) before interpreting a high score.

---

## Installation

Requirements: Python 3.8 or higher, and:

```
pandas >= 1.3.0
numpy >= 1.21.0
scipy >= 1.7.0
matplotlib >= 3.4.0
```

The tool is in the `Produtiv/` folder of the laboratory repository:

```bash
git clone https://github.com/LFDGV-UFSC/LFDGV.git
cd LFDGV/Produtiv
pip install -r requirements.txt
```

Using a virtual environment is recommended:

```bash
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

All scripts read their inputs from, and write their outputs to, the current directory. Run them from inside `Produtiv/` (or from the folder that holds your own input files).

---

## Input data

Two files are required. The files included in this folder are the data from the published study (see [Reproducing the published analysis](#reproducing-the-published-analysis)).

### 1. `abundance_data.csv`: abundance matrix

Six taxonomic columns followed by one column per sample.

```
Phylum,Class,Order,Family,Genus,Species,F3A1,F3A2,F3A3,F3B1,F3B2,F3B3,...
Bacteroidota,Flavobacteriia,Flavobacteriales,Flavobacteriaceae,Flavobacterium,Flavobacterium buctense,0.0104,0.0027,0,0.0398,0.0017,0.0254,...
```

- Taxonomic columns, with exactly these names: `Phylum, Class, Order, Family, Genus, Species`
- Every other column is treated as a sample and must appear in the metadata
- Values: **relative abundance as proportions (0–1)**. The abundance filter is fixed at 0.005 (0.5%), so a table in percentages (0–100) would make the filter ineffective. Divide by 100 before running.

### 2. `metadata.tsv`: sample design

Tab-separated, one row per sample. `metadata.txt` or `metadata.csv` (comma-separated) are also accepted.

```
sample-id	group	condition	parameter
F3A1	F3	A	3500
F3A2	F3	A	3500
F3B1	F3	B	3500
F4A1	F4	A	9500
```

| Column | Description |
|---|---|
| `sample-id` | Must match the sample column names in `abundance_data.csv` exactly |
| `group` | Experimental unit (farm, site, treatment) |
| `condition` | Exactly `A` (before) or `B` (after) |
| `parameter` | Numerical value used for weighting. It must be the same for all samples of a group; only the first value of each group is read |

**Note on notation.** In the metadata, `A` = before and `B` = after. Equation 1 of the article uses *A* and *B* the other way round (*A* for post-cultivation and *B* for pre-cultivation means). The calculation is the same; only the letters differ.

If no metadata file is found, `pws.py` writes a template, `metadata_example.tsv`.

---

## Usage

### Main pipeline

```bash
python main_pws.py
```

The pipeline runs `pws.py`, `filtering5.py` and `plotting_script.py` in sequence and then moves all outputs to `out_pwd/`. You will be asked:

1. Whether to clean `out_pwd/` if it already exists (`y`/`n`).
2. For a one-word name for your parameter, used in plot file names and labels:

```
My parameter for correlation is: Productivity
```

Use a single word without special characters (for example `Productivity`, `pH`, `Temperature`).

### Step by step

```bash
python pws.py               # 1. PWS calculation
python filtering5.py        # 2. Statistical tests and filtering
python plotting_script.py   # 3. Bar plots
```

When run separately, the scripts write their outputs to the current directory instead of `out_pwd/`.

---

## Reproducing the published analysis

The included `abundance_data.csv` and `metadata.tsv` are the species-level table and sample design of the study: 4 farms (F3–F6), 2 sampling times (before planting and after harvest), 3 samples per farm and time (24 samples), with garlic yield (kg/ha) as the parameter.

The article reports Welch's t-test for all comparisons. To reproduce it:

```bash
python pws.py                   # writes parameter_weighted_analysis.csv and top30.txt
python filtering_welch.py       # Welch's t-test for all taxa
python plot_weighted_scores.py  # bar plot of the Welch results
```

`filtering_welch.py` reads `parameter_weighted_analysis.csv` and `top30.txt` from the current directory. Run `pws.py` first, as above. **Do not** use `main_pws.py` for this step: it moves those files to `out_pwd/`, and `filtering_welch.py` will not find them.

The Welch outputs (`filtrado1.csv`, `welch_filtering_report.txt`, `parameter_weighted_barplot_welch.png`) are written to the current directory.

---

## Output files

### Main pipeline (written to `out_pwd/`)

| File | Content |
|---|---|
| `parameter_weighted_analysis.csv` | PWS, per-group scores, relative changes and means for every taxon, sorted by PWS |
| `top30.txt` | Summary with the 30 most positive and the 30 most negative PWS values |
| `sup_statistics.csv` | For every taxon and comparison: means, SDs, Shapiro–Wilk and Levene p-values, test selected and p-values |
| `filtrado.csv` | Taxa passing filters, per-taxon test selection (color column `Cor`) |
| `filtrado_conservativetest.csv` | Taxa passing filters, single test for all taxa (color column `Color`) |
| `filtering_report.txt` | Counts per significance class, tests used and comparison between the two filtered sets |
| `{Parameter}_weighted_barplot.png` | Bar plot, per-taxon tests |
| `{Parameter}_weighted_barplot2.png` | Bar plot, single test |
| `pipeline_log.txt` | Execution time stamp and list of generated files |

`parameter_weighted_analysis.csv` also has an `Impact_Classification` column based on PWS: High (> 0.1), Moderate (0.05–0.1) and Low (0–0.05), positive or negative.

### Welch workflow (written to the current directory)

| File | Content |
|---|---|
| `filtrado1.csv` | Taxa passing filters, all p-values from Welch's t-test |
| `welch_filtering_report.txt` | Counts per significance class and top 10 globally significant taxa |
| `parameter_weighted_barplot_welch.png` | Bar plot |

### Significance classes (bar plot colors)

| Color (column value) | Meaning |
|---|---|
| Red (`vermelho`) | Global test significant (p < 0.05) |
| Orange (`laranja`) | Significant within more than one group |
| Yellow (`amarelo`) | Significant within a single group |
| Gray (`cinza`) | Not significant |

---

## Statistical testing and filtering

For each taxon, two kinds of comparison are made:

- **Global:** all `A` samples vs all `B` samples
- **Per group:** `A` vs `B` within each group

All tests are two-sided, with α = 0.05 and no correction for multiple comparisons. The p-values are descriptive and are used for screening, not for formal inference.

### Per-taxon test selection (`filtrado.csv`)

For each comparison:

| Condition | Test |
|---|---|
| Both conditions normal (Shapiro–Wilk p > 0.05) and equal variances (Levene p > 0.05) | Student's t-test |
| Both normal, unequal variances | Welch's t-test |
| Any condition non-normal | Mann–Whitney U |

Zeros are removed before the Shapiro–Wilk test. If fewer than 3 non-zero values remain, normality cannot be tested and the data are treated as normal, so a t-test is used.

### Single test for all taxa (`filtrado_conservativetest.csv`)

One test is chosen for the whole dataset by screening up to 100 randomly selected taxa (global and per-group comparisons):

- Mann–Whitney U, if any non-normal comparison is found
- otherwise Welch's t-test, if any unequal variances are found
- otherwise Student's t-test

The test chosen is reported in `filtering_report.txt` and in `sup_statistics.csv` (`Conservative_Test_Type`).

**Small samples.** With 3 samples per condition, the smallest possible two-sided Mann–Whitney p-value is 0.1. If Mann–Whitney is chosen, no per-group comparison can reach p < 0.05; only the global comparison (12 vs 12 samples in the example data) can. This is one reason the published analysis used Welch's t-test.

### Welch's t-test for all taxa (`filtrado1.csv`)

`scipy.stats.ttest_ind(..., equal_var=False)` for every comparison. This is the approach used in the article.

### Filtering criteria (all workflows)

A taxon is kept when it meets both conditions:

1. Mean relative abundance ≥ 0.005 (0.5%) across all samples
2. At least one significant comparison (global or per group, p < 0.05), **or** listed in `top30.txt` (30 most positive or 30 most negative PWS)

---

## Troubleshooting

| Problem | What to check |
|---|---|
| Samples reported as missing | `sample-id` values in the metadata must match the sample column names in `abundance_data.csv` exactly |
| No taxa pass the abundance filter, or almost all do | Abundances must be proportions (0–1), not percentages |
| `filtering_welch.py` cannot find `parameter_weighted_analysis.csv` | Run `pws.py` in the same directory first (see [Reproducing the published analysis](#reproducing-the-published-analysis)) |
| Many `Insufficient_data` or p = 1.0 entries | Fewer than 2 samples in a condition for that comparison |
| Plot errors | The parameter name must be a single word without special characters |

---

## Citation

If you use PWS, please cite the article:

> Azevedo GPR, Fritsche Y, Giovanetti LK, da Silva D, Marques CEC, Mazzuco V, Welter LJ, Stefenon VM (2026). Changes in Soil Microbiome Linked to Management of *Allium sativum* Cultivation: Impacts on Soil Health and Agricultural Productivity. *Journal of Soil Science and Plant Nutrition*. https://doi.org/10.1007/s42729-026-03680-0

```bibtex
@article{azevedo2026garlic,
  author  = {Azevedo, Gustavo P. R. de and Fritsche, Yohan and Giovanetti, Leonardo Khaoê and da Silva, Danielle and Marques, Cassiano Eric de Carvalho and Mazzuco, Volni and Welter, Leocir José and Stefenon, Valdir M.},
  title   = {Changes in Soil Microbiome Linked to Management of \textit{Allium sativum} Cultivation: Impacts on Soil Health and Agricultural Productivity},
  journal = {Journal of Soil Science and Plant Nutrition},
  year    = {2026},
  doi     = {10.1007/s42729-026-03680-0}
}
```

---

## License

The code in this folder is released under the [MIT License](LICENSE).

---

## Contact

Developed at the Laboratory of Developmental Physiology and Plant Genetics (LFDGV), Department of Plant Sciences, Federal University of Santa Catarina (UFSC), Florianópolis, Brazil.

- Gustavo P. R. de Azevedo: gustavo.pitta@posgrad.ufsc.br
- Valdir M. Stefenon: valdir.stefenon@ufsc.br

Questions and bug reports can be opened as [issues](https://github.com/LFDGV-UFSC/LFDGV/issues) in this repository.
