# PPI Network Neighborhood Analyzer

A Python-based GUI application for analyzing Protein-Protein Interaction (PPI) networks using the Hishigaki algorithm to predict functional candidate genes.

## Features

- Network Analysis: Load and analyze PPI networks from TSV/CSV files
- Hishigaki Algorithm: Calculate χ² scores to predict candidate functional proteins using formula: χ² = (observed - expected)² / expected
- Self-Consistency Test: Automatically determine optimal neighborhood distance (n=1,2,3) with accuracy validation
- Visualization: Generate multiple plots including:
  - Degree distribution histograms
  - Annotation distribution plots (histogram + boxplot)
  - Hishigaki score distributions (histogram + boxplot)
  - Top 10 candidate proteins bar chart
- GUI Interface: User-friendly graphical interface with console output
- Export Results: Save analysis results as CSV files for further analysis

## Project Structure

Project Files/
├── data/                    # Input data files
│   ├── DREB2A_A.thaliana.tsv     # PPI network data
│   └── AT_stress_proteins.txt    # Known functional proteins
├── results/                 # Output files
│   ├── degree_distribution.png
│   ├── annotation_distribution.png
│   ├── hishigaki_boxplot.png
│   ├── hishigaki_distribution.png
│   └── hishigaki_results.csv
├── src/                    # Source code
│   ├── ppi_analyzer.py     # Core analysis logic
│   ├── PNetAnalyzer.png    # GUI application icon
│   └── gui_app.py          # GUI application
└── README.md


## Installation

### Prerequisites
- Python 3.8 or higher
- pip (Python package manager)

### Required Libraries
Install dependencies using pip:

```bash
pip install networkx pandas numpy matplotlib scipy seaborn
```

For GUI functionality:
```bash
pip install tk
```

## Usage

### Method 1: GUI Application (Recommended)
Run the GUI application:
```bash
cd src
python gui_app.py
```

**GUI Steps:**
1. Click "Browse" to select PPI network file (TSV format)
2. Click "Browse" to select function proteins file (TXT format)
3. Click "Load Data" to load the network
4. Click "Run Analysis" to perform Hishigaki analysis
    -Self-consistency test (finds optimal n value)
    -Annotated neighbor counting
    -Hishigaki score calculation
5. View results in the "Results" tab (ranked list of candidates)
6. Click "Show Plots" to visualize distributions
7. Click "Save Results" to export predictions as CSV

### Method 2: Command Line
Run the analyzer directly:
```bash
cd src
python ppi_analyzer.py
```

## Input File Formats

### 1. PPI Network File (TSV/CSV)
Tab-separated or space-separated file with protein interactions:
```
Protein1    Protein2
AAC42       LIG6
AAC42       NRPB4
NRPB4       NRPB11
```

**Supported formats:**
- Space-separated
- Tab-separated (.tsv)
- Comma-separated (.csv)

### 2. Function Proteins File (TXT)
Tab-separated file with known functional proteins:
```
#TAIR:locus:1005716764	HOS10	AT1G35515	MYB8
TAIR:locus:2060529	ERD15	AT2G41430	LSR1|CID1
```

ID extraction:
-Automatically extracts protein IDs from all columns
-Supports multiple IDs separated by |, ;, ,, or spaces
-Removes taxon prefixes (e.g., 4932.) and version suffixes (.1, .2)
-Handles Arabidopsis (AT1G...), yeast (Y...), and human gene symbols


## Output Files

After analysis, the following files are generated in the `results/` folder:

1. degree_distribution.png - Histogram of protein connection counts
2. annotation_distribution.png - Distribution of annotated neighbors(histogram + boxplot)
3. hishigaki_distribution.png - Distribution of Hishigaki χ² scores(histogram + boxplot)
4. hishigaki_results.csv- All proteins ranked by prediction scores (descending order)

## Algorithm Details

### Hishigaki Score Formula
The algorithm calculates scores using:
```
χ² = (observed - expected)² / expected
```
Where:
observed = Number of known functional proteins in the neighborhood
expected = Total neighbors × Background probability (p)
p = (Functional proteins in network) / (Total proteins)


## Background Probability (p)

The background probability represents the chance that a random protein has the function:
-Calculated as: p = total_functional / total_proteins
-In typical analysis: p ≈ 0.4758 (47.6%)

## Self-Consistency Test

Automatically determines optimal neighborhood distance (n):
1. Takes 25 functional + 25 non-functional proteins from the network
2. For each n=1, 2, 3, predicts function based on neighbor fraction (>30% threshold)
3. Compares predictions with actual labels
4. Selects n with highest accuracy (typically n=1 with 72% accuracy)

Note: 
The best_n parameter is initialized as 1 but updated by the self-consistency test to the optimal value (typically 1, giving 72% accuracy).

## GUI Interface

The GUI provides:
- Console Tab: Real-time formatted outputs and logging
- Results Tab: Ranked list of top 50 candidate proteins as a Dataframe
- Plot Viewer: Separate window with tabs for all visualizations
- Status Bar: Current operation status with color coding (blue=loading, green=complete, red=error)
- Export Functions: Save results as CSV with file dialog

### Example Analysis Output

```
============================================================
📊 LOADING NETWORK
============================================================
  📁 File: DREB2A_A.thaliana.tsv
  ✅ Loaded 4,743 interactions
  📄 Data format: 13 columns

  ✅ NETWORK LOADED SUCCESSFULLY
  🧬 Proteins: 269
  🔗 Interactions: 4,743

============================================================
📈 NETWORK SUMMARY
============================================================
  🧬 Proteins              : 269
  🔗 Interactions          : 4,743
  📊 Average connections   : 35.3
============================================================

============================================================
📈 SELF-CONSISTENCY TEST
============================================================
  🔍 Testing neighborhood distances (n=1,2,3):
------------------------------------------------------------
  📊 Found 128 functional proteins in network
  📊 Found 141 non-functional proteins in network
  🧪 Testing on 25 functional + 25 non-functional proteins
  n=1: ✅ 72.0%  (36/50)
  n=2: ⚠️ 50.0%  (25/50)
  n=3: ⚠️ 50.0%  (25/50)

============================================================
✅ OPTIMAL DISTANCE: n=1 (accuracy: 72.0%)
============================================================

🔍 NEIGHBOR ANALYSIS: AAC42
============================================================
  🔗 Total neighbors              : 70
  ✅ Annotated neighbors count    : 11 (15.7%)
  📌 Annotated neighbors found   : LIG6, NRPB4, NRPE1, NRPB9A...
  📊 Expected by chance           : 33.3 (47.6%)
  📉 This is 22.3 FEWER than expected
============================================================

============================================================
📊 CALCULATING χ² HISHIGAKI SCORES
============================================================
  📈 Background probability          : p = 0.4758
  🧬 Total proteins                  : 269
  📚 Annotated proteins              : 37859
  ✅ Functional proteins in network  : 128
============================================================

  ✅ Calculated χ² Hishigaki scores for 269 proteins
  🔗 Using neighborhood distance: n=1
  📊 Score range: [0.00, 15.37]
------------------------------------------------------------

Top 5 predictions:
  1. NRPD4: 15.3740
  2. NRPB11: 15.2309
  3. AAC42: 14.9413
  4. F8A5.14: 14.9413
  5. NRPB3: 14.8816
```


## Troubleshooting

### Common Issues

1. "File not found" error
   - Ensure files are in the `data/` folder
   - Check file paths in the GUI

2. Missing dependencies
   - Run: `pip install -r requirements.txt`
   - Or install packages individually

3. GUI won't start
   - Ensure Tkinter is installed: `sudo apt-get install python3-tk` (Linux)
   - Or install ActiveTcl (Windows)

4. No results generated
   - Check that function proteins file contains valid IDs
   - Verify PPI network file format

5. Plots not showing
    -Ensure matplotlib backend is compatible. 
    -Add matplotlib.use('TkAgg') before imports.

### Error Messages
- "Graph is not Loaded!": Load PPI network first
- "No function proteins loaded!": Select valid function proteins file
- "Analysis failed": Check console for specific error details

## Scientific Background

This tool implements the Hishigaki method (Hishigaki et al., 2001) for functional prediction in biological networks. The method uses the χ² statistic to identify proteins that have significantly more functional neighbors than expected by chance, based on the "guilt by association" principle.
Applications include:
- Identifying novel gene functions in stress response pathways
- Prioritizing candidate genes for experimental validation
- Understanding functional modules in PPI networks
- Cross-species comparative network analysis

## Author

- Name: A.D.D.C. Wijethunge
- Index: s16562
- Course: Special topics in bioinformatics (BT3172)
- Institution: Department of Plant Sciences, Faculty of Science, University of Colombo

## Support

For questions or issues:
1. Check the troubleshooting section above
2. Review console output for error messages
3. Verify input file formats match specifications
4. Ensure all dependencies are installed
```
