# dsg

Computational model of the direction selective ganglion cells (DS cells) of the mouse retina.

## Overview

These scripts implement the fitted model of OODS retinal ganglion cell populations from Zylberberg, Cafaro, Turner, et al. Neuron 2016
[Direction selective circuits shape noise to ensure a precise population code](http://jzlab.org/team_retina_oods_neuron2016.pdf)

## Usage

Included should be everything you need to run the DS population model, and to reproduce the theoretical Fisher information calculations, in addition to the spiking responses recorded simultaneously from pairs of DS cells.
Below is a description of each file, the hierarchy of files, etc.
Please cite our Neuron paper for any work derived from use of this code and/or data.
All materials are provided as-is, with no guarantees.
Please alert me to any issues / bugs that you encounter Also, if you have any questions, don't hesitate to ask.

### DS cell model fitted to experimental data

For computing the amount of information the modeled population conveys about the stimulus.


The script ooDS_model.m loads up the model parameters (ooDS_params.mat) and runs the model, returning the activity stats for different stimuli.
ooDS_Fisher.m calls ooDS_model.m, and then uses the activity stats to compute the Fisher info.
Results should be able to reproduce Fig. 6E

In matlab:

```matlab
>> ooDS_Fisher
```

### Theoretical Calculations

For theoretical calcuations about how stimulus-dependent correlations between neurons affect information transmission from the eyes to the brain.

File details:

- The script TC.m contains the tuning curve shapes. 
- TC.m is called by do_FI_calc_stimdep.m, which generates responses from the tuning curves, and computes the covariance matrices using Eq. 2 of the paper (assuming Poisson noise). It then computes the Fisher info for the case of stim dependent correlations, and for the case of “matched” constant correlations. 
- There are two "looper" scripts that repeat the calculation for different population sizes and correlation strengths. looper_FI_HOMOG.m specifies homogeneous tuning curves, whereas looper_FI_HETEROG.m specifies heterogeneous ones. 

To Run for heterogeneous tuning curves (note, this may take several minutes to run):

```matlab
>> looper_FI_HETEROG.m
```

To Run for homogeneous tuning curves (note, this may take several minutes to run):

```matlab
>> looper_FI_HOMOG.m
```

### Data

Extracellular spiking recordings used to generate Fig. 1 of the paper, and to fit the model in Fig. 6

- data/readme.rtf contains a description of the matlab files and their formatting
- data/DSpairs_spikeData.mat contains the data

## Contact

Please contact Joel Zylberberg (joel.zylberberg [at] ucdenver [dot] edu) to report any bugs.