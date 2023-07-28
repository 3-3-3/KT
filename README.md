# KT
This repository contains code and resources for my research during the Summer of 2023, which investigated the usefulness of using KT diagrams–analagous to Hertzsprung-Russell diagrams in electromagnetic astronomy–in studying populations of white dwarf binaries which will be visible to future space based gravitational wave detectors. We generated data with COSMIC ([https://cosmic-popsynth.github.io/](https://cosmic-popsynth.github.io/)); in this repository, I include some of the scripts we used to analyze COSMIC data, as well as other important resources:

* KT.py: this is where most of the code I wrote is contained. It is purely functional, and contains functions for KT plotting, as well as for processing COSMIC data.
* Example_Creating_Final_Fix.ipynb: an example of using functions in KT.py to generate a csv file from a bunch of fixed population .h5 files generated with cosmic-pop. Additionally, provides a discussion of how these fixed populations are generated using cosmic-pop.
* Example_SNR: an example of using signal to noise ratio analysis (with the function in KT.py) to limit the fixed population to what will be visible with LISA. Also contains examples of producing some KT diagrams with functions in KT.py. 
* \Agate: directory for the Agate project, which is meant to generate a galaxy population containing mass transferring systems.
  * buildFiles: the build files I have been using for generating fixed populations with realistic star formation sampling. 
  * Dumb_Galaxy_Maker.py: an adaption of cosmic-pop, which loops until a critical amount of mass has been sampled rather than until convergence is reached. This will allow us to make the Agate galaxy (although, it is going to take a few months).

Using Dumb_Galaxy_Maker.py is similar to using cosmic-pop, but its being run explicitely with python, whereas cosmic-pop is set up to run out of bin. Additionally, we need to tell Dumb_Galaxy_Maker.py the total mass it needs to sample using the --total_mass argument (which I have added. So, for example, if we were to build the bulge, we would navigate to the Agate directory in terminal and run:

```
ipython3 Dumb_Galaxy_Maker.py -- --inifile buildFiles/bulge_Params.ini --final-kstar1 10 13 --final-kstar2 10 13 --total_mass 8.9e9 -n 4
```

The first -- passes no arguments to ipython3 so that we can pass the remaining commands to Dumb_Galaxy_Maker.py's argument parser. The commands are as follows: --inifile tells Dumb_Galaxy_Maker.py what inifile to use, --final-kstar1 and --final-kstar2 tell Dumb_Galaxy_Maker.py which range of kstars to keep, --total_mass tells Dumb_Galaxy_Maker.py how much mass to sample, and -n 4 tells Dumb_Galaxy_Maker.py to use four threads when evolving the systems. I empirically found this to be the fastest choice on the iMacs in the CC computational physics lab–likely anauther choice is better for a different architecture!
