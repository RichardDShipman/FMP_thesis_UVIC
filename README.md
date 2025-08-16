# FMP Thesis UVIC

**Author:** Richard Shipman
**University:** El Campus Vic de la Universitat de Vic - Universitat Central de Catalunya (UVic-UCC)

This repository contains the code and data for my Final Masters Project (FMP). The project focuses on the analysis of glycopeptides and includes three main tools:

1.  **Glycopeptide Sequence Finder**: A Python tool for identifying theoretical glycopeptide sequences from protein FASTA files.
2.  **MSA Plotter with Sequon Markers**: An R script for analyzing and visualizing multiple sequence alignments with highlighted N-sequons.
3.  **PSM Analysis**: A Python script for analyzing and visualizing experimental glycopeptide Peptide-Spectrum Match (PSM) data from mass spectrometry experiments.

## About This Project

This project aims to provide a suite of tools for the in-silico and experimental analysis of glycopeptides. The tools can be used to predict potential glycopeptides, analyze their conservation across different species, and process and visualize experimental data from mass spectrometry.

## Tools

### Glycopeptide Sequence Finder

This tool identifies glycopeptide sequences from protein FASTA files. It can be used to generate theoretical libraries of glycopeptides for further analysis.

**Key Features:**

*   Supports various proteases for in-silico digestion.
*   Identifies N-linked, O-linked, and C-linked glycosylation sequons.
*   Calculates peptide properties such as mass, hydrophobicity, and isoelectric point.
*   Provides batch processing capabilities.

For more information, please see the [Glycopeptide Sequence Finder README](Glycopeptide_Sequence_Finder/README.md).

### MSA Plotter with Sequon Markers

This R script aligns sequences from a FASTA file and generates a readable, text-like Multiple Sequence Alignment (MSA) plot using 'ggplot2'. The plot is wrapped, highlights N-sequons, and adds labeled markers for them.

**Key Features:**

*   Performs multiple sequence alignment using the `msa` package.
*   Can also use pre-aligned FASTA files.
*   Generates a wrapped MSA plot for readability.
*   Highlights N-sequons and adds markers with position labels.

For more information, please see the [MSA Plotter with Sequon Markers R script](MSA_plotter_with_sequon_markers/MSA_plotter_with_sequon_markers.R).

### PSM Analysis

This Python script analyzes and visualizes experimental glycopeptide Peptide-Spectrum Match (PSM) data from FragPipe. It is designed to perform site-specific analysis of glycoforms.

**Key Features:**

*   Parses `psm.tsv` files from FragPipe.
*   Categorizes glycoforms based on their composition.
*   Generates various plots, including glycan coverage, pie charts of glycoform distribution, and total glycan composition.
*   Creates summary files for N- and O-glycans.

For more information, please see the [PSM Analysis Python script](psm_analysis/glycopeptide_psm_data_analysis.py).

