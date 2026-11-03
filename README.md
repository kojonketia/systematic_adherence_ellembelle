# Mass drug administration adherence patterns in the Ellembelle district, Ghana
This study addressed the adherence patterns of Lymphatic Fllariasis in the Ellembelle district, Ghana. This repository contains the R scripts used for the statistical analysis and simulation study of adherence to MDA for LF in Ellembelle District, Western Region, Ghana.

The study investigates patterns of adherence to LF MDA, with particular emphasis on identifying individuals who consistently participate in repeated MDA rounds (systematic adherence) and those whose participation varies between rounds (non-systematic adherence). The analyses use de-identified MDA register data collected from communities in Ellembelle District during 2019, 2021 and 2022. The repository provides the analytical code to support transparency, reproducibility, and further investigation of MDA adherence patterns.

## Study Background
Lymphatic filariasis is a neglected tropical disease targeted for elimination as a public health problem through preventive chemotherapy using mass drug administration. In Ghana, ivermectin and albendazole are administered annually to eligible populations in endemic communities. The effectiveness of MDA programmes depends not only on population-level treatment coverage but also on repeated participation by individuals across successive treatment rounds. Individuals who consistently miss treatment may constitute a persistent reservoir of infection, potentially undermining the effectiveness of elimination programmes. Understanding individual-level adherence patterns is therefore important for identifying gaps in treatment participation and informing strategies to improve MDA programme effectiveness. The MDA registers contain individual treatment records and demographic information used to assess participation across treatment rounds. The study also considers historical MDA implementation and treatment participation to provide context for interpreting adherence patterns.

## Study Objectives
This study aims to:
- characterise individual-level adherence patterns to LF MDA in Ellembelle District.
- examine demographic factors associated with MDA adherence, including age and sex.
- model district antigen prevalence from baseline using an estimated adherence parameter from the data
- use simulation analyses to explore the implications of different adherence patterns and coverage on reaching <2% antigen prevalence.


## Analytical Approach
The repository is organised into two main components:

### 1. Statistical Analysis (`analysis/`)
This directory contains R scripts for data preparation, statistical analysis, and presentation of study findings.
The analyses include:
- data cleaning and preparation of individual-level MDA register records.
- classification of individuals according to their treatment participation across MDA rounds.
- description and statistical comparison of adherence patterns using appropriate descriptive and inferential methods.
- assessment of demographic characteristics associated with adherence.


### 2. Simulation Analyses (`simulations/`)
This directory contains R scripts for simulation analyses related to MDA adherence. The simulation component is intended to simulate the district's antigen prevalence using the observed adherence behaviors and how different patterns of individual treatment participation may influence MDA programme outcomes. The scripts provide the computational framework for exploring alternative adherence scenarios and their potential implications for LF elimination.

## Repository Structure
The repository is organised as follows:
```text
systematic_adherence_ellembelle/
│
├── analysis/
│   └── R scripts for data preparation,
│       statistical analyses, tables and figures
│
├── simulations/
│   └── R scripts for simulation analyses
│
└── README.md
```

## Software Requirements
The analyses were developed using R. The following R packages may be required, depending on the specific script being executed:
- `tidyverse`
- `dplyr`
- `ggplot2`
- `readxl`
- `gtsummary`
- `Hmisc`
- `parallel`

Additional packages are documented in the individual R scripts.
To install the main packages, run:

```r
install.packages(c("tidyverse", "dplyr", "ggplot2", "readxl", "gtsummary", "Hmisc", "parallel"))
```

## Reproducing the Analyses

### Step 1: Clone the repository
Clone the repository to your local computer:

```bash
git clone https://github.com/kojonketia/systematic_adherence_ellembelle.git
```

Navigate to the project directory:
```bash
cd systematic_adherence_ellembelle
```

### Step 2: Prepare the data
Ensure that the anonymised MDA register data required by the analysis scripts are available in the appropriate directory. Where necessary, update the file paths in the R scripts to reflect the location of the data on your computer.

**Note:** Individual-level study data should only be accessed and used in accordance with the applicable ethical approvals and data-sharing requirements.

### Step 3: Run the statistical analyses
Open the relevant R scripts in the `analysis/` directory using RStudio. Execute the scripts in the required order, beginning with data preparation, followed by statistical analysis and the generation of tables and figures.

### Step 4: Run the simulation analyses
Open the scripts in the `simulations/` directory and execute them according to their specified dependencies and simulation settings. Simulation parameters can be modified to investigate alternative MDA adherence scenarios.

### Step 5: Generate study outputs
The scripts generate the statistical summaries, tables, figures, and simulation results used to report the study findings. The output locations and file formats are specified within the individual scripts.

## Data Availability
The repository contains the analytical and simulation code developed for this study. Individual-level MDA register data are de-identifirf and subject to applicable ethical and data-sharing restrictions. Access to the underlying data may be requested from the corresponding study investigators, subject to the relevant approvals. No personally identifiable participant information should be disclosed or shared.

## Reproducibility
This repository is intended to facilitate reproducibility and transparency in the analysis of individual-level MDA adherence patterns. Researchers interested in applying the analytical framework to other LF-endemic settings are encouraged to adapt the scripts to their respective datasets, taking into consideration differences in MDA implementation, treatment eligibility, data availability, and programme context.

## Citation
If you use the code or analytical approaches presented in this repository, please cite the associated study.
**Associated publication:** To be updated following publication.

## Author and Correspondence
For questions, comments, or further information regarding the analytical code and study, please contact the repository owner through GitHub:
- **GitHub:** https://github.com/kojonketia
- **Repository:** https://github.com/kojonketia/systematic_adherence_ellembelle

