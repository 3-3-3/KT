# KT
This repository contains code and resources for my research during the Summer of 2023, which investigated the usefulness of using KT diagrams–analagous to Hertzsprung-Russell diagrams in electromagnetic astronomy–in studying populations of white dwarf binaries which will be visible to future space based gravitational wave detectors. We generated data with COSMIC ([https://cosmic-popsynth.github.io/](https://cosmic-popsynth.github.io/)); in this repository, I include many of the scripts we used to analyze this data, as well as build files for COSMIC. Here is the list of files and how they have been used:

-KT.py: this includes helper functions for the scripts I created. 

-bpp.ipynb: bpp refers to the complete stellar evolution history output by COSMIC. This script was used to generate diagrams which helped us investigate the evolutionairy paths systems which ended up in similar places on our KT diagram took to get there. 

-Evolve.ipynb: this script uses COSMIC to evolve individual double white dwarf systems to get an intuition for their evolution in a KT diagram.

-Counts.ipynb: this script creates histograms for the number of semidetached systems in a fixed population for comparison against Figure 2 in [https://arxiv.org/abs/0705.3272](https://arxiv.org/abs/0705.3272).

-fix_the_fix.ipynb: this script uses the fixed populations for the Tormaline A galaxy to create a KT diagram. Note: this will not be representative of the milky way because all systems have been evolved to 13000 Myr.

-Time_Sampled_Fix.ipynb: this script uses correctly time sampled fixed populations we generated with Cosmic-Pop to get a more representative sample of the milky way galaxy. 

-Tourmaline_B.ipynb: this script creates a KT diagram for the Tourmaline B galaxy, which is the most up to date (as of this writing) galaxy population created with COSMIC.
