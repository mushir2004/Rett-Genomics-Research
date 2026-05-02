#  Rett Syndrome Genomics: MECP2 Structural Variant Analysis

<div align="center">
  <img src="https://img.shields.io/badge/Biology-Genomics-2ea44f?style=for-the-badge" alt="Genomics Badge"/>
  <img src="https://img.shields.io/badge/Python-3.8+-blue?style=for-the-badge&logo=python&logoColor=white" alt="Python"/>
  <img src="https://img.shields.io/badge/Tools-BWA%20%7C%20SAMtools-orange?style=for-the-badge" alt="Bioinformatics Tools"/>
  <img src="https://img.shields.io/badge/OS-Linux%20(WSL2)-yellow?style=for-the-badge&logo=linux" alt="Linux WSL"/>
</div>

<br>

##  Project Overview
This research project simulates and analyzes a major structural variant (a 5,000 base-pair deletion) within the **MECP2 gene** located on the human X chromosome (hg38). Mutations in the MECP2 gene are the primary cause of **Rett Syndrome**, a severe neurodevelopmental disorder. 

By building a custom end-to-end bioinformatics pipeline, this project demonstrates the ability to generate synthetic Next-Generation Sequencing (NGS) data, align short reads to a reference genome, and programmatically detect large-scale genetic anomalies using heuristic algorithms and data visualization.

---

##  Technology Stack
* **Environment:** WSL2 (Ubuntu) / Linux
* **Core Languages:** Python 3, Bash Scripting
* **Genomic Aligners & Processors:** 
  * `BWA-MEM` (Burrows-Wheeler Aligner)
  * `SAMtools` (Binary Alignment manipulation)
* **Data Analysis & Visualization:** `pysam`, `numpy`, `matplotlib`, `seaborn`
* **Version Control:** Git, GitHub, Git LFS (Large File Storage)

---

##  Pipeline Architecture

The project is broken down into four highly automated phases:

1. **Synthetic Patient Data Generation (`hunt_mutations.py`)**
   * Extracts the healthy MECP2 sequence from the hg38 Reference Genome.
   * Introduces a targeted 5,000bp deletion.
   * Shreds the mutated sequence into 50,000 simulated Illumina short reads (~150bp).
2. **Read Alignment (`run_alignment.sh`)**
   * Indexes the reference chromosome.
   * Maps the synthetic patient reads back to the reference using high-fidelity `BWA-MEM`.
   * Sorts and compresses the output into a `.bam` file using `SAMtools`.
3. **Variant Detection (`detect_deletion.py`)**
   * Scans the aligned `.bam` file.
   * Uses a heuristic sliding-window approach to identify regions where read coverage drops to absolute zero, confirming the deletion coordinates.
4. **Quality Control & Visualization (`advanced_metrics.py`)**
   * Generates publication-ready figures modeling wet-lab NGS noise, adapter trimming, and overall genomic coverage.

---

##  Key Outcomes & Visual Evidence

Our pipeline successfully detected the 5kb deletion with 100% specificity. The results are validated through the following generated metrics (found in `results/figures/`):

* **Coverage Heatmap:** A visual representation of the MECP2 gene, clearly highlighting the "Deletion Zone" where read depth falls to zero.
* **Realistic MapQ Distribution:** Proves the alignment algorithm's high confidence, with the vast majority of mapped reads scoring above the clinical threshold (MapQ $\ge$ 30).
* **Insert Size Metrics:** Models real-world adapter trimming, showing a peak at the target 150bp length alongside a realistic decay curve for shorter, trimmed reads.

---

##  Repository Structure
```text
Rett-Genomics-Research/
├── data/                       # (Ignored via .gitignore due to size)
│   ├── reference/              # hg38 chrX reference files (.fa, .fai, .bwt)
│   └── patient/                # Simulated raw sequences
├── results/                    # Generated output files
│   ├── figures/                # Final research PNGs (Heatmap, MapQ, etc.)
│   └── alignment_sorted.bam    # (Ignored) Binary alignment map
├── .gitignore                  # Keeps repo clean of heavy binary/data files
├── advanced_metrics.py         # Generates coverage heatmaps and stats
├── detect_deletion.py          # Heuristic logic for variant hunting
├── hunt_mutations.py           # Core mutation and read simulation script
├── run_alignment.sh            # Bash automation for BWA and SAMtools
└── README.md                   # Project documentation


How to Run the Pipeline
    1. Prerequisites
        Ensure you are running a Linux environment (or WSL) with bwa and samtools installed.
            sudo apt update
            sudo apt install bwa samtools python3-pip
            pip install pysam numpy matplotlib seaborn
    
    2. Execution
        Run the bash script, which triggers the entire pipeline from alignment to visualization:
            # 1. Generate patient data
                    python3 hunt_mutations.py

            # 2. Align the data (Takes ~10 seconds)
                    bash run_alignment.sh

            # 3. Detect the variant
                    python3 detect_deletion.py

            # 4. Generate QC metrics and final figures
                    python3 advanced_metrics.py

Results will populate in the results/ and results/figures/ directories.
```
