# Information

This is the code for the paper <br>
*Bayesian inference for biophysical neuron models enables stimulus optimization for retinal neuroprosthetics* <br>
by Oesterle et al. 2020.<br>
The paper can be found here: https://elifesciences.org/articles/54997

This repo can be used to reproduce the data and figures shown in the paper.
It also provides a python wrapper for some of the NeuronC functionality.

# Code-structure

While the core python code can be found [here](pythoncode), the experiments performed in the paper are described in jupyter notebooks (see below).
The experients code is structured similar to the paper. Since most experiments depend on the previous steps, the should be run in order. However, the data can also be loaded from zenodo, which allows to skip experiments or single (time-consuming) steps in single experiments.
Either way you have to [download the data](step0a_download_data) to get full functionality.

## Experiments
The experiments are ordered as follows:

- [Download the data.](step0a_download_data)
- [Proprocess cell iGluSnFR target data and the stimuli.](step0b_preprocess_iGluSnFR_data)
- [Parameter inference for the cone cell.](step1a_optimize_cones)
- [Remove single ion channels from the optimized cone to test effect of active ion channels.](step1b_analyse_optimized_cones)
- [Parameter inference for the CBC cells.](step2a_optimize_cbc)
- [Remove single ion channels from the optimized CBCs to test effect of active ion channels.](step2b_analyse_optimized_cbcs)
- [Parameter inference of the electrical parameters of the retina.](step3a_optimize_electrical_params)
- [Estimate the thresholds of the CBC relative to the GC thresholds estimates.](step3b_thresholds)
- [Estimate stimulus waveforms for selective stimulating of the OFF- or ON-CBC.](step4_optimize_stimulus)

After running the experiments or downloading the data, you can generate the [figures](create_figures) and the animations [animations](create_animations).

## NeuronC

The NeuronC [[1]](#1) version (6.3.14) used [here](neuronc) was downloaded from http://retina.anatomy.upenn.edu/~rob/neuronc.html and used with minor modifications.
To use it you have to run the makefile in [neuronc](neuronc) and in [retsim](neuronc/models/retsim).

## Electrical stimulation

The simulation of the electrical stimulation in the paper involded the software COMSOL. If you don't have access to this software, you can download the extracellular voltages from zenodo and simulate the electrical stimulation. If you want to recompute the extracellular voltage, you need to download the COMSOL files and rerun. Note that the experiments prepare the environment for this simulation, and the notebooks clearly indicate when COMSOL needs to be run. 

### Stimulus optimization

To optimize the stimuli we created a pipeline such that COMSOL can be run on a machine where it is installed, while SNPE and NeuronC can be run on a different machine.
For this we created two directories [COMSOL2retsim_COMSOL](step4_optimize_stimulus/COMSOL2retsim_COMSOL) and [COMSOL2retsim_interface](step4_optimize_stimulus/COMSOL2retsim_interface). The former should be move to a computer that can run both COMSOL *and* jupyter notebook. The notebook run on this computer will communicate with COMSOL and save the COMSOL output to the [COMSOL2retsim_interface](step4_optimize_stimulus/COMSOL2retsim_interface), i.e. it needs read and write permission. Then you can run [step4_optimize_stimulus.ipynb](step4_optimize_stimulus/1_optimize_stimulus.ipynb) on the second machine, that also needs read and write permission for [COMSOL2retsim_interface](step4_optimize_stimulus/COMSOL2retsim_interface). It will load the COMSOL output, use it to run SNPE and NeuronC, and tell the other notebook when and how to create new COMSOL outputs.

Similarily, this can be done for the [optimization of the electrical parameters of the retina](step3a_optimize_electrical_params).

# References
<a id="1">[1]</a> 
Smith, Robert G. (1992). 
NeuronC: a computational language for investigating functional architecture of neural circuits.
J. Neurosci. Methods 43: 83-108.
