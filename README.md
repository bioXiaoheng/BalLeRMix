###  BalLeRMix---*Bal*ancing selection *L*ik*e*lihood *R*atio *Mix*ture models

This repository hosts the software package for BalLeRMix and scripts used in the study "*Flexible mixture model approaches that accommodate footprint size variability for robust detection of balancing selection*" (Cheng &amp; DeGiorgio 2019, Submitted). 

- For the software, go to `BalLeRMix/software/`
- For scripts used for SLiM simulations, go to `BalLeRMix/Simulation_scripts/`
- For scripts used in empirical analyses,  go to `BalLeRMix/Empirical_analysis/`

Please cite the following manuscript if using this software:

~~Xiaoheng Cheng, Michael DeGiorgio (2019) [Flexible mixture model approaches that accommodate footprint size variability for robust detection of for localizing balancing selection.](https://www.biorxiv.org/content/10.1101/645887v2) *bioRxiv* doi.org/10.1101/645887~~

Xiaoheng Cheng, Michael DeGiorgio (2020) [Flexible mixture model approaches that accommodate footprint size variability for robust detection of balancing selection.](https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msaa134/5848012) *Molecular Biology and Evolution* msaa134, doi.org/10.1093/molbev/msaa134

------

In BalLeRMix v2, we introduce the `-m <m>` argument to customize the presumed number of alleles being balanced at the selected sites, in case you want to look for multi-allelic balancing selection. The default value is 2.

Complete user manual for BalLeRMix v2 is coming soon...

**2020.6.22-Update:** Updated the model for multi-allelic balancing selection in v2.2.
**2020.2.5-Update:** Fixed a minor bug in the initialization module.

------

## Quick Guide 

```
usage: BalLeRMix.py [-h] -i INFILE --spect SPECTFILE [-o OUTFILE] [-m M]
                      [--getSpect] [--getConfig] [--nofreq] [--nosub] [--MAF]
                      [--physPos] [--rec RRATE] [--fixSize] [-w R]
                      [--noCenter] [-s STEP] [--fixX X] [--rangeA SEQA]
                      [--listA LISTA]
```
You can use `python BalLeRMix.py -h` to see the more detailed help page.

### 1. Input format
For B<sub>0</sub> and B<sub>2</sub> statistics, the user should first generate the __*tab-delimited*__ site frequency spectrum file, __*without header*__, *e.g.*:
> \<k\>|\<sample size n\>|\<proportion in the genome\>    
> :-----:|:-----:|-----
> 1|50|0.03572
> 2|50|0.02024
> ...
  
or the configuration file with polymorphism/substitution ratio, __*without header*__, *e.g.*:

> \<sample size n\> | \<\% of substitutions\> | \<\% of polymorphisms\>  
> :-----:|:-----:|:-----:   
> 50  |  0.7346  |  0.2654    

The input files should have four columns, presenting physical positions, genetic positions, number of derived (or minor) alleles observed, and total number of alleles observed (*i.e.* sample size). This file should be tab-delimited and should have a header, *e.g.*:

> physPos|genPos|x|n    
> :-----:|:-----:|:-----:|:-----:    
> 16|0.000016|50|50    
> 35|0.000035|12|50   
> ...
  
### 2. Running the *B* statistics
To perform B<sub>2</sub> scans on your input data, use

    python BalLeRMix.py -i <input> --spect <derived allele frequency spectrum> -o <output>

To perform B<sub>2,MAF</sub> scans on your input data, use

    python BalLeRMix.py -i <input> --spect <minor allele frequency spectrum> -o <output> --MAF

To perform B<sub>1</sub> scans on your input data, use

    python BalLeRMix.py -i <input> --config <sub/poly configuration file> -o <output> --nofreq

To perform B<sub>0</sub> scans on your input data, use

    python BalLeRMix.py -i <input> --config <derived allele frequency spectrum> -o <output> --nosub 

To perform B<sub>0,MAF</sub> scans on your input data, use

    python BalLeRMix.py -i <input> --config <minor allele frequency spectrum> -o <output> --nosub --MAF

### 3. Generate helper files
To generate spectrum file for B<sub>2</sub>:

    python BalLeRMix.py -i <concatenated input> --getSpect --spect <spectrum file name>

To generate spectrum file for B<sub>2,MAF</sub>:

    python BalLeRMix.py -i <concatenated input> --getSpect --MAF --spect <spectrum file name>

To generate spectrum file for B<sub>1</sub>:

    python BalLeRMix.py -i <concatenated input> --getConfig --spect <config file name>

To generate spectrum file for B<sub>0</sub>:

    python BalLeRMix.py -i <concatenated input> --getSpect --nosbub --spect <spectrum file name>

To generate spectrum file for B<sub>0,MAF</sub>:

    python BalLeRMix.py -i <concatenated input> --getSpect --nosub --MAF --spect <spectrum file name>

### 4. Customizing the scan
All arguments besides the aforementioned ones are for customizing the scan.

- `[--physPos] [--rec RRATE] `:

   Because `BalLeRMix` uses genetic distances (in cM) to compute likelihood, to direct the software to use physical positions instead, you should use `--physPos`, and indicate the uniform recombination rate (cM/nt) in your species of interest with `--rec`. The default value is 10<sup>-6</sup> cM/nt.
   
   This argument will be automatically incurred if you choose to fix the window size (*e.g.*, 1000bp, 5kb, *etc.* ), in which case yuou want to make sure the software is correctly informed of the recombination rate. Using physical positions will also change how you define window sizes and step sizes, if you were to customize the scanning window.

- ` [--fixX X] [--rangeA SEQA] [--listA LISTA]`:

   These areguments allow you to specify the parameter space that the software optimizes over. The presumed equilibrium frequency is *x*, and the rate of decay in linkage disequilibrium is *A*. If you choose to look for multi-allelic balancing selection where more than two alleles are being balanced, *x* should be a vector of descending equilibrium frequencies, and should match the number of balanced alleles you chose (via `-m`) to scan for.

- ` [--fixSize] [-w R] [--noCenter] [-s STEP] [--physPos]`:

   These areguments are for customizing the scanning window. You probably won't need them because `BalLeRMix` is robust to window sizes. For more details on how these arguments work, check the v1 software manual.
