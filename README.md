# aayushm-masters-thesis
Code for my master's thesis -  dealing with response diversity and ecosystem stability in the Western Ghats

This repository will have all the code I write during the course of my master's thesis, in folders as listed below:

1.) gams 

a.) rd_calculations.rmd - This is a R Markdown file containing code to estimate species occurrence probability curves along the Western Ghats using Generalised Additive Models. This represents the stacked species distribution model (SSDM) approach to modelling. The file also contains the code to estimate the divergence aspect of response diversity (Ross et al., 2023). The data used is species presence-absence data from Krishnadas et al., 2021 - encompassing 183 species from 286 sites across the Western Ghats. 

b.) geb13350-sup-0002-tree-species-occ - data for the above

2.) HMSC

a.) hmsc_procedure - RMD file containing code for the Hierarchical Modelling of Species Communities workflow as laid out in Tikhonov et al., 2019. HMSC is used to 'jointly' model multispecies occurrence probabilities (a JSDM, as opposed to the GAM-derived SSDM), and this is further used to estimate divergence. The raw dataset for this was given to me by Dr. Meghna Krishnadas, and after filtering for minimal species presence (> 10 individuals across the Western Ghats), it is in the form of abundance data for 363 species from 286 sites. 
