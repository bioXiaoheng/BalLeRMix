###  BalLeRMix---*Bal*ancing selection *L*ik*e*lihood *R*atio *Mix*ture models

This repository hosts the software package for BalLeRMix and scripts used in the study "*Flexible mixture model approaches that accommodate footprint size variability for robust detection of balancing selection*" (Cheng &amp; DeGiorgio 2019, Submitted). 

- For the software, go to `BalLeRMix/software/`
- For scripts used for SLiM simulations, go to `BalLeRMix/Simulation_scripts/`
- For scripts used in empirical analyses,  go to `BalLeRMix/Empirical_analysis/`

Please cite the following manuscript if using this software:

X Cheng, M DeGiorgio (2019) [Flexible mixture model approaches that accommodate footprint size variability for robust detection of for localizing balancing selection.](https://www.biorxiv.org/content/10.1101/645887v2) *bioRxiv* doi.org/10.1101/645887

Complete user manual for BalLeRMix v2 is coming soon...

------

#### Quick guide 

##### 1. Input format
For B<sub>0</sub> and B<sub>2</sub> statistics, the user should first generate the tab-delimited site frequency spectrum file, without header,e.g.:

> 1 50  0.3572
>
> 2 50  0.2024
>
> ...\<k\>  \<sample size n\>  \<proportion in the genome\>...
>
> ...
  
or the configuration file with polymorphism/substitution ratio, without header:

> 50  0.7346 0.2654
>
> ... \<sample size n\> \<\% of substitutions\>  \<\% of polymorphisms\>
>
  
##### 2. Running the *B* statistics
To perform B<sub>2</sub> scans on your input data, use

`python BalLeRMix.py -i <input> --spect <derived allele frequency spectrum> -o <output>`
