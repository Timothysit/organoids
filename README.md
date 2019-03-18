# The organoid project 

This repository contains the code used for the microelectrode array (MEA) data analysis in the paper: 

[**Cerebral organoids at the air-liquid interface generate diverse nerve tracts with functional output**](https://www.nature.com/articles/s41593-019-0350-2)

by 

Stefano L. Giandomenico (1), Susanna B. Mierau (2), George M. Gibbons (3), Lea M.D. Wenger (3), Laura Masullo (1), Timothy Sit (2), Magdalena Sutcliffe (1), Jerome Boulanger (1), Marco Tripodi (1), Emmanuel Derivery (1), Ole Paulsen (2), András Lakatos (3),(4), Madeline A. Lancaster (1)

(1) MRC Laboratory of Molecular Biology, Cambridge Biomedical Campus, Cambridge CB2 0QH, UK
(2) Department of Physiology, Development and Neuroscience, University of Cambridge, Cambridge, CB2 3EG, UK.
(3) John van Geest Centre for Brain Repair and Division of Stem Cell Neurobiology, Department of Clinical Neurosciences, University of Cambridge, E.D. Adrian Building, Forvie Site, Robinson Way, Cambridge, CB2 0PY, UK.
(4) Wellcome Trust-MRC Cambridge Stem Cell Institute, Cambridge Biomedical Campus, Cambridge CB2 0AH, UK

This includes: 

1. Signal filtering 
2. Spike detection 
3. Correlation analysis 
4. And their visualisations 

(Note that supplementary materials here are not peer-reviewed)

All analysis were performed using `MATLAB`, and the main scripts used for analysis are `ogranoidProject.m` for spike detection, and `organoidCorrelation` for analysis of functional connectivity. 

Most visualisations were performed using `MATLAB`, with the exception of the network plot, which was made using the `R` interface for `igraph`. The code for this is available in the folder `R_analysis`. Placement of scale bars and legends for some plots were adjusted in Adobe Illustrator.






