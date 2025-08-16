### **Product Requirements Document: Glycoform Analysis Automation Tool**

**Author:** Gemini, Project Manager
**Last Updated:** June 7, 2025
**Status:** Version 1.0 Finalized

-----

### **1. Introduction**

Glycoproteomics research involves complex data analysis to understand which sites on a protein are glycosylated and with what types of glycans. The output from mass spectrometry search engines, like the FragPipe-MSFragger-Glyco workflow, is rich but requires significant post-processing to yield biological insights. Researchers often perform repetitive, manual, and error-prone data wrangling and plotting to analyze site-specific glycan occupancy and composition.

This document outlines the requirements for the **Glycoform Analysis Automation Tool**, a command-line Python script designed to automate the post-processing of Peptide-Spectrum Match (PSM) data from FragPipe. The tool will standardize the analysis pipeline, from data ingestion and cleaning to the generation of publication-quality plots and summary reports, enabling researchers to accelerate their discovery process.

### **2. Goals and Objectives**

  * **Primary Goal:** To provide researchers with a robust, automated tool for performing site-specific glycoform analysis on a target protein from FragPipe-generated data.
  * **Key Objectives:**
      * **Reduce Analysis Time:** Drastically cut down the manual effort required to process and visualize glycoproteomics data.
      * **Improve Reproducibility:** Ensure consistent, standardized analysis by encapsulating the workflow in a single, version-controlled script.
      * **Enhance Data Interpretation:** Generate a comprehensive set of plots and summary statistics that facilitate the clear interpretation of glycosylation patterns.
      * **Increase Accessibility:** Provide a straightforward command-line interface that is accessible to bioinformaticians and bench scientists with basic computational skills.

### **3. Target Audience & User Personas**

**Primary User Persona: Dr. Evelyn Reed, Glycobiology Researcher**

  * **Role:** Principal Investigator/Postdoctoral Researcher in a molecular biology lab.
  * **Expertise:** Deep knowledge of protein biochemistry and mass spectrometry. Comfortable with data analysis concepts but is not a professional programmer. Uses tools like Excel, Prism, and occasionally R or Python for basic tasks.
  * **Goals:**
      * Wants to quickly see the glycosylation profile of her protein of interest after a mass spectrometry experiment.
      * Needs to compare glycoform distribution across different sites on the protein.
      * Needs to create figures for lab meetings, publications, and grant proposals.
  * **Pain Points:**
      * Spends hours or days manually filtering and pivoting `psm.tsv` files in Excel.
      * Struggles to create consistent and high-quality plots for all her analyses.
      * Worries about making manual errors when calculating glycosite positions or summarizing data.

### **4. Functional Requirements**

The tool shall be a Python script executed from the command line.

| **Feature ID** | **User Story** | **Acceptance Criteria** |
| :--- | :--- | :--- |
| **REQ-101** | As a researcher, I want to provide my `psm.tsv` file, the name of my target protein, and its signal peptide length as inputs so that the analysis is specific to my experiment. | \<ul\>\<li\>The script must accept three command-line arguments: `--input` (file path, required), `--protein` (protein name, required), and `--signal` (integer, optional, default=0).\</li\>\<li\>The script must gracefully handle cases where the input file is not found.\</li\>\<li\>The script must inform the user if no PSMs are found for the specified protein.\</li\>\</ul\> |
| **REQ-201** | As a researcher, I want the tool to automatically calculate the correct glycosylation site positions on my mature protein. | \<ul\>\<li\>The script must read the `protein_start` and `best_positions` columns from the input file.\</li\>\<li\>It must correctly calculate the final glycosite position by adjusting for the protein start position and the user-provided signal peptide length.\</li\>\<li\>It must correctly propagate the calculated glycosite position to all PSMs belonging to the same peptide sequence.\</li\>\</ul\> |
| **REQ-202** | As a researcher, I want the tool to automatically parse complex glycan composition strings and classify them into understandable biological categories. | \<ul\>\<li\>The script must parse strings like `HexNAc(2)Hex(5)Fuc(1)` for monosaccharide counts.\</li\>\<li\>It must use a defined logic to classify glycans into primary categories: `Oligomannose`, `Complex`, `Hybrid`, `O-Glycan`, `EndoE`, `Unoccupied`, etc.\</li\>\<li\>It must assign subcategories where applicable (e.g., `M9`, `M8`).\</li\>\</ul\> |
| **REQ-301** | As a researcher, I want a high-level summary file so I can quickly see the key statistics from my analysis. | \<ul\>\<li\>The script must generate a `.txt` file containing: Total PSMs analyzed, Target Protein PSMs, Unoccupied PSM count, a list of identified glycosites, and a table of PSM counts per glycoform category.\</li\>\</ul\> |
| **REQ-302** | As a researcher, I want a detailed data table in CSV format so I can perform my own downstream analysis if needed. | \<ul\>\<li\>The script must generate a `_filtered_glycoPSMs.csv` file.\</li\>\<li\>This file must contain one row per unique peptide/glycosite/glycan combination and include columns for the parsed glycoform categories and the total PSM count for that combination.\</li\>\</ul\> |
| **REQ-401** | As a researcher, I want a series of publication-quality plots to be automatically generated so I can visualize my results. | The script must generate the following plots as separate `.png` files, all clearly titled and labeled: \<br\>\<ul\>\<li\>**Overall Distribution:** A bar chart showing the total PSM count for each major glycoform category.\</li\>\<li\>**Site Coverage:** A stacked bar chart showing the percentage of each glycoform category at each glycosite.\</li\>\<li\>**Site-Specific Pies:** A separate pie chart for each glycosite, showing the relative abundance of glycoform categories at that site.\</li\>\<li\>**Total Composition:** A bar chart showing the total PSM count for every unique glycan composition across all sites, ordered by mass.\</li\>\<li\>**Site-Category Bars:** A grouped bar chart showing absolute PSM counts for each glycoform category at each site.\</li\>\<li\>**Site-Specific Composition:** A separate, detailed bar chart for each glycosite showing the PSM count for all identified glycan compositions (including those with zero counts at that site).\</li\>\<li\>**Site-Specific Categories:** A separate, detailed bar chart for each glycosite showing the PSM count for all glycoform categories (including those with zero counts at that site).\</li\>\</ul\>|
| **REQ-501**| As a researcher, I want all output files to be saved in a new, clearly named directory so that they don't clutter my current working directory. | \<ul\>\<li\>The script must create a new directory named `{input_file_stem}_analysis_results`.\</li\>\<li\>All generated text files, CSVs, and plots must be saved inside this directory.\</li\>\</ul\> |

### **5. Non-Functional Requirements**

| **ID** | **Requirement** | **Description** |
| :--- | :--- | :--- |
| **NFR-1**| **Usability** | The tool must be run via a simple, single command. Help text detailing the arguments must be available (`-h` or `--help`). |
| **NFR-2**| **Performance** | The script should complete analysis on a typical `psm.tsv` file (e.g., up to 500,000 rows) in under 2 minutes on a standard laptop. |
| **NFR-3**| **Maintainability** | The Python code must be well-commented and organized into logical functions for readability and future updates. |
| **NFR-4**| **Dependencies** | The script will rely on a standard scientific Python environment. All required libraries (`pandas`, `numpy`, `matplotlib`, `seaborn`) must be documented in the script's header or a `requirements.txt` file. |

### **6. Out of Scope / Future Work**

The following features will not be included in version 1.0 but may be considered for future releases:

  * A Graphical User Interface (GUI).
  * Interactive plots (e.g., using Plotly or Bokeh).
  * Support for data from other search engines (e.g., Byonic, pGlyco).
  * Analysis of multiple samples or conditions at once to perform comparative/statistical analysis.
  * Advanced glycan structure representation beyond simple composition.
  * Integration into a larger bioinformatics platform like Galaxy or Nextflow.

### **7. Success Metrics**

The success of the Glycoform Analysis Automation Tool will be measured by:

  * **Adoption:** The tool is successfully used by at least 3 research groups within the first 6 months of release.
  * **Efficiency Gains:** Users report a \>90% reduction in time spent on routine glycoform analysis.
  * **Impact:** The tool is cited or mentioned in the methods section of peer-reviewed publications.
  * **User Feedback:** Positive feedback from users regarding ease of use, quality of output, and reliability.