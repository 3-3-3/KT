# KT
This repository contains code and resources for my research during the Summer of 2023, which investigated the usefulness of using KT diagrams–analagous to Hertzsprung-Russell diagrams in electromagnetic astronomy–in studying populations of white dwarf binaries which will be visible to future space based gravitational wave detectors. We generated data with COSMIC ([https://cosmic-popsynth.github.io/](https://cosmic-popsynth.github.io/)); in this repository, I include some of the scripts we used to analyze COSMIC data, as well as other important resources:

* KT.py: this is where most of the code I wrote is contained. It is purely functional, and contains functions for KT plotting, as well as for processing COSMIC data.
* Example_Creating_Final_Fix.ipynb: an example of using functions in KT.py to generate a csv file from a bunch of fixed population .h5 files generated with cosmic-pop. Additionally, provides a discussion of how these fixed populations are generated using cosmic-pop.
* Example_SNR: an example of using signal to noise ratio analysis (with the function in KT.py) to limit the fixed population to what will be visible with LISA. Also contains examples of producing some KT diagrams with functions in KT.py. 
* \Agate: directory for the Agate project, which is meant to generate a galaxy population containing mass transferring systems.
  * buildFiles: the build files I have been using for generating fixed populations with realistic star formation sampling. 
  * Galaxy_Maker.py: an adaption of cosmic-pop, which loops until a critical amount of mass has been sampled rather than until convergence is reached. This will allow us to make the Agate galaxy (although, it is going to take a few months). 
