
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10360718.svg)](https://doi.org/10.5281/zenodo.13959115)

# Integrative docking benchmark 

This repository contains a benchmark of hetero-dimeric complexes, with component atomic structures and chemical crosslinks. It is used to compare integrative modeling methods: EASAL and IMP.

## **Benchmark**

There are 30 cases in benchmark, varying in crosslinker type, number of crosslinks, source of crosslinks, and source of monomer structures of the constituent proteins. The benchmark cases are divided based on the source of crosslinks in `benchmark` directory. Refer to the `benchmark/simulated/README.md` and `benchmark/experimental/README.md` for details.

## **Installing and running integrative modeling softwares**

See `installation.md` for instructions for installing and running IMP and EASAL 

## **Comparing IMP and EASAL generated model ensembles**

The model ensembles obtained from IMP and EASAL are compared based on fit to the input crosslinks, similarity with the native structure, and efficiency. See `scripts/compare_ensembles/` for the scripts to compare the results from IMP and EASAL. 


## **Information**
**Author(s):** Yichi Zhang, Muskaan Jindal, Shruthi Viswanath, Meera Sitharam  
**License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
International License.  
**Publications:** Zhang, Y., _et._ al.A new discrete-geometry approach for integrative docking of proteins using chemical crosslinks.
 
 
 
 
 
 
 
 
