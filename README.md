🧬 Peak Calling and Residual Regression Toolkit

This repository provides a complete workflow for **1 kb–resolution peak detection**, **peak refinement**, and **residual regression analysis** of chromatin accessibility or EdU incorporation datasets.  
It includes two Bash pipelines and one R analysis script.


1. 📁 Contents

| Script | Description |
|---------|--------------|
| `peak_calling.1kb.sh` | Performs genome-wide peak calling at 1 kb resolution directly from a bigWig coverage file. Iteratively refines peaks (three-pass detection) to identify narrow or broad regions of enrichment. |
| `refine_peaks_1kb_resolution.sh` | Refines any input BED file of peaks to 1 kb resolution using the corresponding bigWig. Identifies the most enriched 1 kb bin (summit) within each input region. |
| `residual_regression.R` | Performs residual regression analyses and generates annotated plots comparing ON/OFF conditions for ATAC-seq or EdU incorporation data. Can also be adjusted for other data types. |

---

2. ⚙️ Dependencies

🧩 Common Requirements

All scripts assume a **Unix/Linux** environment.

| Tool | Purpose |
|------|----------|
| `bedtools` | Genomic interval operations (binning, merging, intersection) |
| `bigWigAverageOverBed` | Extracts bigWig signal over genomic bins |
| `datamash` | Statistical summaries (percentiles, medians) |
| `bash` | Shell scripting environment |
| `R` | For statistical and visualization analysis (residual regression) |


📦 R Package Requirements

Install required R libraries before running `residual_regression.R`:

```r
install.packages(c("ggplot2", "dplyr", "broom"))
```


🧱 Installation of Command-Line Dependencies

### Ubuntu / Debian

```bash
sudo apt update
sudo apt install -y bedtools datamash wget
```

or install it from source

```
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make
```

Install bigWigAverageOverBed (UCSC binary)

```
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed
chmod +x bigWigAverageOverBed
sudo mv bigWigAverageOverBed /usr/local/bin/
```

macOS (Homebrew)

```bash
brew install bedtools datamash
brew install ucsc-bigwigaverageoverbed
```

---

3.🚀 Usage

#1️⃣ Peak Calling at 1 kb Resolution

Run interactively:

```bash
bash peak_calling.1kb.sh
```

You will be prompted for:
- Reference genome coordinates (BED format)
- bigWig file of signal
- Number of threads
- Output name

Output:  
`<output_name>.1kb_resolution.peaks.bed` — list of high-confidence peaks merged across all chromosomes.

> 🔧 TIPS

- Adjust sensitivity-specificity in line 37 by modifying the perc value of "datamash perc:95".  This is practically the background threshold. Values 1-100.  The lower the value the more peaks will be called. It is suggested to lower this value for broad peak calling. 
- To split broader peaks with multiple summits modify perc value in lines 55,78.  Values 1-100. Higher values split the broader peaks in multiple summits.  For peak calling of proader peaks low values are suggested.
- If peaks are too many and close enough (e.g from potential sequencing gaps or artifacts, or broad peak calling), please merge the detected peaks by bedtools merge '-d' values in lines 37,  60,  83, where '-d'  applies for the distance in bases (e.g. -d 10000 stands for merging peaks found in a distance of less than 10kb).


#2️⃣ Peak Refinement of Existing BED Files

Refine pre-called peaks to 1 kb bins and identify summits:

```bash
bash refine_peaks_1kb_resolution.sh
```

You will be prompted for:
- Input peaks BED file
- bigWig file
- Number of threads

Output:  
`<input_name>.summits.bed` — 1 kb peak summits per region.


#3️⃣ Residual Regression Analysis

The R script expects tidy TSV files with consistent column names.

Run the R script:

```bash
Rscript residual_regression.R
```

The script will:
- Load provided `.tsv` data files (`ATAC_OCT4_data.tsv`, `EdU_4h_OCT4_12h_data.tsv`, `EdU_4h_OCT4_18h_data.tsv`) or adjust the data accordingly.
- Perform residual regression per replication timing group (ES, MS, LS)
- Generate annotated regression plots for each dataset.

Output:  
Three plots (`p1`, `p2`, `p3`) showing relationships between ON/OFF conditions and regression statistics.

---

4. 👩‍🔬 Citation

If you use or adapt these scripts, please cite your dataset and this repository appropriately.
