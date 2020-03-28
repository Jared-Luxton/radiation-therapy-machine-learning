# Telomere length and chromosomal instability for predicting individual radiosensitivity and risk via machine learning

This repository is a supplement to my upcoming paper. 
&nbsp;   

## One sentence summary
Here we explore the utility of telomeres and chromosome rearrangements for predicting cancer patient's risks of radiation late effects from radiotherapy, and present the first implementation of individual telomere length data in a machine learning model (XGBoost), providing a general framework for predicting individual outcomes from radiotherapy.
&nbsp;   
  
## Abstract
The ability to predict a cancer patient’s response to radiotherapy and risk of developing adverse late health effects would greatly improve personalized treatment regimens and individual outcomes. Telomeres represent a compelling biomarker of individual radiosensitivity and risk, as exposure can result in dysfunctional telomere pathologies that coincidentally overlap with many radiation-induced late effects, ranging from degenerative conditions like fibrosis and cardiovascular disease to proliferative pathologies like cancer. Here, telomere length was longitudinally assessed in a cohort of fifteen prostate cancer patients undergoing Intensity Modulated Radiation Therapy (IMRT) utilizing Telomere Fluorescence in situ Hybridization (Telo-FISH). To evaluate genome instability and enhance predictions for individual patient risk of secondary malignancy, chromosome aberrations were also assessed utilizing directional Genomic Hybridization (dGH) for high-resolution inversion detection. We present the first implementation of individual telomere length data in a machine learning model, XGBoost, trained on pre-radiotherapy (baseline) and in vitro exposed (4 Gy γ-rays) telomere length measures, to predict post-radiotherapy telomeric outcomes, which together with chromosomal instability provide insight into individual radiosensitivity and risk for radiation-induced late effects.


## Data handling/workflow
The notebooks which document and execute the data analysis for this project are provided below in multiple formats. All code is written in python. The Nbviewer and HTML links provide static renderings of the jupyter notebooks. The Jupyter lab Binder link provides an interactive, fully functional rendering of the entire repo or individual notebooks. Each notebook is dedicated to a specific task in the workflow, detailed in the filename suffix.
&nbsp;   

---
**01_radiation_therapy_patients_data_EXTRACTION_CLEANING.ipynb**
* Extraction and cleaning of telomere length and chromosome rearrangement data 

**02_radiation_therapy_patients_data_VISUALIZATION_STATISTICS.ipynb**
* Multiple types of highly customized data visualizations
* Feature engineering short/long telomeres
* ANOVAs, regressions, statistical modeling

**03_radiation_therapy_patients_data_MACHINE_LEARNING.ipynb**
* Extensive data-cleaning class pipelines to prepare data machine learning (ML) models
* GridSearch / Bayesian optimization (don't run if your laptop is older like mine)
* Training XGBoost models on pre-therapy data to predict post-therapy otucomes
* Extensive evaluation of model metrics
* Hierarchical clustering of patients using longitudinal data into high/low risk sub-groups

---
### To launch all  notebooks (may be slow):
| Jupyter Lab |
| ---                       |
| [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Jared-Luxton/radiation-therapy-machine-learning/master?urlpath=lab)|

### To launch individual notebooks:
| Nbviewer | Jupyter Lab | HTML |
| ---      |  ---        | ---  |
| [01_radiation_therapy_patients_data_EXTRACTION_CLEANING.ipynb](https://nbviewer.jupyter.org/github/Jared-Luxton/radiation-therapy-machine-learning/blob/master/notebooks/01_radiation_therapy_patients_data_EXTRACTION_CLEANING.ipynb) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Jared-Luxton/radiation-therapy-machine-learning/master?urlpath=lab/tree/notebooks%2F01_radiation_therapy_patients_data_EXTRACTION_CLEANING.ipynb) | [HTML](https://raw.githack.com/Jared-Luxton/radiation-therapy-machine-learning/master/notebooks/html_copy_notebooks/01_radiation_therapy_patients_data_EXTRACTION_CLEANING.html) |
| [02_radiation_therapy_patients_data_VISUALIZATION_STATISTICS.ipynb](https://nbviewer.jupyter.org/github/Jared-Luxton/radiation-therapy-machine-learning/blob/master/notebooks/02_radiation_therapy_patients_data_VISUALIZATION_STATISTICS.ipynb)|[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Jared-Luxton/radiation-therapy-machine-learning/master?urlpath=lab/tree/notebooks%2F02_radiation_therapy_patients_data_VISUALIZATION_STATISTICS.ipynb)|[HTML](https://raw.githack.com/Jared-Luxton/radiation-therapy-machine-learning/master/notebooks/html_copy_notebooks/02_radiation_therapy_patients_data_VISUALIZATION_STATISTICS.html) |
| [03_radiation_therapy_patients_data_MACHINE_LEARNING.ipynb](https://nbviewer.jupyter.org/github/Jared-Luxton/radiation-therapy-machine-learning/blob/master/notebooks/03_radiation_therapy_patients_data_MACHINE_LEARNING.ipynb)|[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Jared-Luxton/radiation-therapy-machine-learning/master?urlpath=lab/tree/notebooks%2F03_radiation_therapy_patients_data_MACHINE_LEARNING.ipynb)|[HTML](https://raw.githack.com/Jared-Luxton/radiation-therapy-machine-learning/master/notebooks/html_copy_notebooks/03_radiation_therapy_patients_data_MACHINE_LEARNING.html) |

&nbsp;    
