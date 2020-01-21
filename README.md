# Utility of telomere length and chromosome rearrangements for predicting radiation late effects

This repository is a supplement to my upcoming paper. 
  
  
## One sentence summary
---
Here we explore the utility of telomeres and chromosome rearrangements in developing predictive tools that identify patients at high risk for developing radiation late-effects post-therapy to improve personalized treatment options.
  
  
## Background
---
Telomeres are critical features at the ends of linear chromosomes, which guard the genome (DNA) against degradation and inappropriate activation of the DNA damage response. Telomeres shorten with cell division and therefore aging, but also from oxidative stress, lifestyle factors (i.e nutrition), and environmental exposures (i.e pollutants, radiation).  Short telomeres are associated with a range of aging-related pathologies, including cardiovascular disease (CVD), reduced immune function and cancer risk; thus the rate at which telomeres shorten links aging and aging-related pathologies, providing a biomarker which integrates an individual’s lifestyle exposures. Measurements of telomeres in response to radiation thus could provide a general assessment of an individual’s health, providing inference on future health outcomes and risks. Chromosome rearrangements (dicentrics, inversions, translocations) are well-established markers of ionizing radiation (IR) exposure, frequently used in biodosimetry, and associated with virtually all cancers. A patient's short- and long-term tendencies to develop chromosome rearrangements post-radiation exposure could thus provide insight on their personal risks for developing cancer after radiation therapy. Radiation late-effects refer to a broad range of debilitating, often permanent side effects experienced by radiation-sensitive patients after radiation therapy. These late-effects range from cardiovascular disease, scarring of heart and/or lung tissue, cognitive deficits, and even secondary cancers; and are of particular concern for pediatric patients.
  
   
## Research problem and objective
---
At present, no efficacious clinical means exist for predicting a patient’s radiation-sensitivity using pre-therapy data, preventing proper guidance of treatment decisions and unnecessary patient harm. Here, we explore the utility of telomeres propose generation of an imaging dataset could enable development of a machine learning tool capable of predicting which patients are likely to develop radiation late-effects using pre-therapy telomere data. Such a tool, capable of accurately predicting radiation late-effects, would improve personalized treatment decisions for patients predicted to be high-risk.
  

## Data handling/workflow
---
The notebooks which document and execute the data analysis for this project are provided below in multiple formats. All code is written in python. The Nbviewer and HTML links provide static renderings of the jupyter notebooks. The Jupyter lab Binder link provides an interactive, fully functional rendering of the entire repo or individual notebooks. Each notebook is dedicated to a specific task in the workflow, detailed in the filename suffix.

### 01_radiation_therapy_patients_data_EXTRACTION_CLEANING.ipynb
* Extraction and cleaning of telomere length and chromosome rearrangement data 

### 2_radiation_therapy_patients_data_VISUALIZATION_STATISTICS.ipynb
* Multiple types of highly customized data visualizations
* Feature engineering short/long telomeres
* ANOVAs, regressions, statistical modeling

### 3_radiation_therapy_patients_data_MACHINE_LEARNING.ipynb
* Extensive data-cleaning class pipelines to prepare data for building machine learning (ML) models
* Building XGBoost models trained on telomere data to predict either mean telomere length or # short telomeres (success)
* Building XGBoost models to predict incidence of chromosome rearrangements post-therapy (fails)
* GridSearch / Bayesian optimization (don't run if your laptop is older like mine)
* Extensive looks at model metrics

### 4_radiation_therapy_patients_data_DATA_STORY.ipynb
* **IN PROGRESS**
* Cohesive data story in jupyter notebook format
* Highlights and major implications from the analysis

| Launch entire repo (Jupyter Lab) |
| ---                       |
| [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Jared-Luxton/radiation-therapy-machine-learning/master?urlpath=lab)|

| Nbviewer | Jupyter Lab | HTML |
| ---      |  ---        | ---  |
| [01_radiation_therapy_patients_data_EXTRACTION_CLEANING.ipynb](https://nbviewer.jupyter.org/github/Jared-Luxton/radiation-therapy-machine-learning/blob/master/notebooks/01_radiation_therapy_patients_data_EXTRACTION_CLEANING.ipynb) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Jared-Luxton/radiation-therapy-machine-learning/master?urlpath=lab/tree/notebooks%2F01_radiation_therapy_patients_data_EXTRACTION_CLEANING.ipynb) | [HTML](https://raw.githack.com/Jared-Luxton/radiation-therapy-machine-learning/master/notebooks/01_radiation_therapy_patients_data_EXTRACTION_CLEANING.ipynb) |
| [02_radiation_therapy_patients_data_VISUALIZATION_STATISTICS.ipynb](https://nbviewer.jupyter.org/github/Jared-Luxton/radiation-therapy-machine-learning/blob/master/notebooks/02_radiation_therapy_patients_data_VISUALIZATION_STATISTICS.ipynb)|[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Jared-Luxton/radiation-therapy-machine-learning/master?urlpath=lab/tree/notebooks%2F02_radiation_therapy_patients_data_VISUALIZATION_STATISTICS.ipynb)|[HTML](https://raw.githack.com/Jared-Luxton/radiation-therapy-machine-learning/master/notebooks/02_radiation_therapy_patients_data_VISUALIZATION_STATISTICS.ipynb) |
| [03_radiation_therapy_patients_data_MACHINE_LEARNING.ipynb](https://nbviewer.jupyter.org/github/Jared-Luxton/radiation-therapy-machine-learning/blob/master/notebooks/03_radiation_therapy_patients_data_MACHINE_LEARNING.ipynb)|[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Jared-Luxton/radiation-therapy-machine-learning/master?urlpath=lab/tree/03_radiation_therapy_patients_data_MACHINE_LEARNING.ipynb)|[HTML](https://raw.githack.com/Jared-Luxton/radiation-therapy-machine-learning/master/notebooks/03_radiation_therapy_patients_data_MACHINE_LEARNING.ipynb) |
| [04_radiation_therapy_patients_data_STORY_HIGHLIGHTS.ipynb](https://nbviewer.jupyter.org/github/Jared-Luxton/radiation-therapy-machine-learning/blob/master/notebooks/04_radiation_therapy_patients_data_STORY_HIGHLIGHTS.ipynb)|[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Jared-Luxton/radiation-therapy-machine-learning/master?urlpath=lab/tree/04_radiation_therapy_patients_data_STORY_LEARNING.ipynb)|[HTML](https://raw.githack.com/Jared-Luxton/radiation-therapy-machine-learning/master/notebooks/04_radiation_therapy_patients_data_STORY_HIGHLIGHTS.ipynb) |



TO DO:  
organize readme
