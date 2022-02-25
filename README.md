# Analyzing experiment LR-035B

## Installation
Rstudio is required

## Authors
Lisa Rodenburg

## Running the code
Input data should be Excel data from organoid swelling assays where foldchange and AUC values are already calculated.
Thid code can be used to make the graphs.

## Acknowledgments
Sam
Course wiring reproducible code

## Requirements
Language: R 4.0.3
Operating system: Windows 10 x64
Packages: 
  ggplot2_3.3.3   
  readxl_1.3.1    
  dplyr_1.0.4 

## Project organization
- PG = project-generated
- HW = human-writable
- RO = read only
```
.
├── .gitignore
├── CITATION.md
├── LICENSE.md
├── README.md
├── requirements.txt
├── bin                <- Compiled and external code, ignored by git (PG)
│   └── external       <- Any external source code, ignored by git (RO)
├── config             <- Configuration files (HW)
├── data               <- All project data, ignored by git
│   ├── processed      <- The final, canonical data sets for modeling. (PG)
│   ├── raw            <- The original, immutable data dump. (RO)
│   └── temp           <- Intermediate data that has been transformed. (PG)
├── docs               <- Documentation notebook for users (HW)
│   ├── manuscript     <- Manuscript source, e.g., LaTeX, Markdown, etc. (HW)
│   └── reports        <- Other project reports and notebooks (e.g. Jupyter, .Rmd) (HW)
├── results
│   ├── figures        <- Figures for the manuscript or reports (PG)
│   └── output         <- Other output for the manuscript or reports (PG)
└── src                <- Source code for this project (HW)

```


## License

This project is licensed under the terms of the [MIT License](/LICENSE.md)
