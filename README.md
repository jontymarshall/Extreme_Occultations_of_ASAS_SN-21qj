# Extreme_Occultations_of_ASAS_SN-21qj
LCO photometry of ASAS-SN 21qj along with ancilliary data and scripts. 

Optool (https://github.com/cdominik/optool) is required to run the extinction analysis for the dust compositions and size distributions.

```
.
├── analysis            <- Scripts used to produce the figures in the manuscript
├── data
│   ├── external    <- Photometry taken from online repositories (ATLAS, ASASSN, NEOWISE)
│   ├── processed   <- Photometry extracted from LCOGT images/BANZAI catalogue
│   └── raw         <- Scripts to download the LCOGT data from the archive
├── .gitignore      <- Avoids uploading data, credentials, outputs, system files etc
├── LICENCE
├── models      <- Dust composition models for occultation and stellar photosphere used to interpret SED, etc. 
├── paper       <- Manuscript
│   ├── manuscript.pdf  <- Generated manuscript file, equivalent to the arXiv version
│   ├── manuscript.zip  <- Archive containing the latex files and figures used to produce manuscript.pdf
├── README.md         <- This file
└── requirements.txt  <- Python modules used to run scripts in analysis directory
```