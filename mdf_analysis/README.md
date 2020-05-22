# metabolites_LP
A simple Scipy-based linear program to calculate the max-min driving force (MDF) as described by [Noor et al. 2014](https://doi.org/10.1371/journal.pcbi.1003483)
and associated metabolite 
concentrations for a given flux distribution.

The analysis takes a metabolic model in [COBRApy](https://doi.org/10.1186/1752-0509-7-74) or [ScrumPy](http://mudshark.brookes.ac.uk/ScrumPy) format, a flux distribution 
and a set of Gibbs free-energy parameters, and returns an MDF value and the associated set of 
metabolite concentrations.

## Example
The example in the ipython notebook ('mdf_analysis') calculates the MDF and associated metabolite concentrations 
for fluxes computed with flux balance analysis and thermodynamic flux balance analysis of the *Escherichia coli* core 
model (*reference needed*).

The Gibbs free-energy parameters used in this example have been downloaded from 
[equilibrator](https://doi.org/10.1093/nar/gkr874), and are provided in the 'params' folder. 
