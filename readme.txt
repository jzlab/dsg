These scripts implement the fitted model of OODS retinal ganglion cell populations from Zylberberg, Cafaro, Turner, et al. Neuron 2016 “Direction selective circuits shape noise to ensure a precise population code”

Please contact Joel Zylberberg (joel.zylberberg@ucdenver.edu) to report any bugs.

The script ooDS_model.m loads up the model parameters (ooDS_params.mat) and runs the model, returning the activity stats for different stimuli.

ooDS_Fisher.m calls ooDS_model.m, and then uses the activity stats to compute the Fisher info.

Results should be able to reproduce Fig. 6E

To Run:
>> ooDS_Fisher