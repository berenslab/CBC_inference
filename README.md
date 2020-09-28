# CBC_inference
Code etc. for CBC paper of Oesterle el al. 2020.
Can be used to reproduce the data and figure shown in the paper.
Also provides a more python wrapper for some of the NeuronC functionality.

# Code-structure

The python core code used can be found [here](_python_code), must the experiments performed in the paper are described in jupyter notebooks that use this code.
The code is structured similar to the paper. Since most experiments depend on the previous steps, the should be run in order. However, the data can also be loaded from zenodo, which allows to skip experiments or single (time-consuming) steps in single experiments.

The experiments are ordered as follows:

- [Parameter inference for the cone cell](step1a_optimize_cones)
- [Remove single ion channels from the optimized cone to test effect of active ion channels](step1b_analyse_optimized_cones)
- [Parameter inference for the CBC cells](step2a_optimize_cbc)
- [Remove single ion channels from the optimized CBCs to test effect of active ion channels](step2b_analyse_optimized_cbcs)
- [Parameter inference of the electrical parameters of the retins](step3_optimize_COMSOL_params)
- [Estimate the thresholds of the CBC relative to the GC thresholds estimates](step3b_thresholds)
- [Estimate stimulus waveforms for selective stimulating of the OFF- or ON-CBC](step4_optimize_stimulus)

The NeuronC[[1]](#1) version (6.3.14) was downloaded from http://retina.anatomy.upenn.edu/~rob/neuronc.html and used with minor modifications.
The code can be found in NeuronC

## References
<a id="1">[1]</a> 
Smith, Robert G. (1992). 
NeuronC: a computational language for investigating functional architecture of neural circuits.
J. Neurosci. Methods 43: 83-108.
