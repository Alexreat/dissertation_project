Metabolic Network Analysis of Hepatocellular Carcinoma (HCC)

Overview

This repository contains the complete computational workflow for the dissertation titled:

“A Systems-Level Analysis of Metabolic Reprogramming in Hepatocellular Carcinoma using Genome-Scale Modelling and Network Theory”

The project investigates how the metabolic network of Hepatocellular Carcinoma (HCC) diverges from that of healthy liver tissue. RNA-seq data from the TCGA-LIHC cohort was integrated with the Human-GEM genome-scale model to build 423 patient-specific metabolic models.

The analysis combines:
	•	Flux Balance Analysis (FBA) to assess metabolic functional capacity.
 
	•	Graph theory to quantify network topology.
 
	•	Machine learning (XGBoost) to classify samples as tumour or healthy based on network- and flux-derived features.
 

Key finding: HCC metabolism is a streamlined version of the healthy network — it achieves comparable theoretical growth by channeling flux through a smaller, more focused set of reactions rather than through radical structural rewiring.

⸻

Repository Structure

dissertation_project

-notebooks/

│   ├── FBA_per_sample.ipynb         # Flux Balance Analysis per model

│   ├── Graph_extra_features.ipynb   # Additional graph metrics extraction

│   ├── Graph_per_sample.ipynb       # Build per-sample metabolic graphs

│   ├── Graph_stats.ipynb            # Statistical analysis of graph metrics

│   ├── ML_model.ipynb               # Machine learning classification

│   └── hier_complexity.py           # Hierarchical complexity calculation

│
├── imat_script.m                    # MATLAB iMAT workflow (tumour models)

├── imat_script_healthy.m            # MATLAB iMAT workflow (healthy models)

│
├── mart_export.txt                  # Gene ID mapping reference

├── .gitignore                       # Ignore intermediate/data files

└── README.md                        # This file


⸻

Methodology Workflow

	1.	Data Preprocessing
	•	Raw gene expression (TPM) from TCGA-LIHC is log-transformed, Z-score normalized, and binarized for iMAT.
 
	2.	Model Reconstruction
	•	Processed expression vectors are used in MATLAB (COBRA Toolbox) via imat_script.m and imat_script_healthy.m to generate patient-specific metabolic models from Human-GEM.
 
	3.	Flux Balance Analysis (FBA)
	•	Implemented in FBA_per_sample.ipynb to calculate:
	•	Maximal biomass production
	•	Total absolute flux
	•	Active reaction count
 
	4.	Graph Construction & Feature Extraction
	•	Each model is converted into a metabolite-metabolite graph.
	•	Topological metrics (degree, clustering, modularity, assortativity, etc.) are computed in Graph_per_sample.ipynb and Graph_extra_features.ipynb.
 
	5.	Statistical Analysis
	•	Feature distributions compared between tumour and healthy using Mann–Whitney U tests and effect sizes in Graph_stats.ipynb.
 
	6.	Machine Learning Classification
	•	XGBoost model (ML_model.ipynb) trained on combined flux and graph features to distinguish tumour vs. healthy samples.

⸻

Reproducibility

All scripts and workflows required to reproduce the analysis are available in this repository:

📂 https://github.com/Alexreat/dissertation_project/tree/main

The MATLAB scripts require:
	•	MATLAB R2023a+
	•	COBRA Toolbox v3.x
	•	Human-GEM model
	•	TCGA-LIHC RNA-seq (STAR counts), The dataset (424 samples) can be downloaded from **UCSC Xena portal



