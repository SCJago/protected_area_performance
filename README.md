# Trade-offs Between Nature and People Reveal Challenges in Translating Global Conservation Targets Into National Realities
Analysis of Ethiopia's protected area network

This repository contains the code, data outputs, and figures supporting the manuscript:

**_Trade-offs Between Nature and People Reveal Challenges in Translating Global Conservation Targets Into National Realities_**  
Sophie Jago, Gebremeskel Gizaw, Bezawit Genanaw, Joe Langley, Ermias Lulekal, Joseph D. M. White, Adèle N. Rowlands, Tariku Geda., Kumara Wakjira, Fekede Regassa, Sebsebe Demissew, Wendawek Abebe, Julia P. G. Jones, Robert J. Smith, James S. Borrell

(Submitted for publication 2025)
---
## Summary

Achieving global biodiversity targets requires translating international commitments into effective national actions. In 2022, 196 parties agreed to protect 30% of the planet by 2030 under the Global Biodiversity Framework (GBF). However, concerns remain over whether current protected areas (PAs) are ecologically effective and socially equitable.

This project focuses on Ethiopia—a country with globally important biodiversity but high levels of poverty and food insecurity. 

We provide a comprehensive national-scale evaluation of Ethiopia’s progress towards the multiple dimensions of the 30-by-30 target. We assess the extent of Ethiopia’s protected area network and how well it represents national ecoregions and species. We then apply a robust quasi-experimental approach to assess both environmental (forest, agriculture and grassland cover change) and human wellbeing (change in months of adequate food, dietary diversity and material wellbeing) impacts of Ethiopia’s protected areas. Considering protected areas individually, we then examine whether funding allocation is a predictor of performance across environmental and social outcomes. Finally, we explore the views of key national stakeholders in conservation policy and practice and consider the alignment of national priorities and global goals. The research highlights the very real challenges faced by those tasked with turning a global commitment into reality on the ground. 

## Scripts

1. PA_expansion_plot.R - produces Figure 1B showing the expansion of Ethiopia's protected area network over time
2. WWF_ecoregion_PA_coverage_plot.R - produces Figure 2A showing the proportion of each ecoregion in Ethiopia that is protected in the current protected area network
3. Species_Representation_Analysis.R - produces Figure 2B showing the proportion of each species range in Ethiopia that is protected in the curent protected area network and carries out statistical tests to identify differences between taxa. Code written by J Langley, and figure code written by S. Jago
4. PCA_environmental_representation.R - produces Supplementary Figure S6 using a principal component analysis to determine the representation of Ethiopia's environmental space in its protected area network
5. Proprocessing_data_for_matching.RMD - assigns sampling units as treatment, buffer or control units and generates covariate data for each unit, ready for statistical matching 
6. Statistical_matching_and_outputs.RMD - carries out statistical matching for gridcells and households, post-matching checks, statistical analysis on effectiveness outputs, and generates Figure 3B, Figure 3C, Supplmentary Figure S4, Supplementary Figure S5 and Supplmentary Figure S8
7. Matching_robustness_checks.R - This is the code for carrying out matching robustness checks shown in Supplementary Figure S9. This code was adapted from code used in Devenish et al. (2022)  (Devenish, K., Desbureaux, S., Willcock, S., & Jones, J. P. G. (2022). On track to achieve no net loss of forest at Madagascar’s biggest mine. doi:10.1038/s41893-022-00850-7) the figures used to plot the results from our robustness checks were produced using code developed by Ariel Ortiz-Bobea (Ortiz-Bobea, A., Ault, T. R., Carrillo, C. M., Chambers, R. G. & Lobell, D. B. Anthropogenic climate change has slowed global agricultural productivity growth. doi:10.1038/s41558-021-01000-1.)
