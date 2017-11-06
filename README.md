# EAG channels in the pre-synaptic terminal 

This repository contains code for a simulation of a presynaptic terminal, where we model a Calcium microdomain using a simple multi-compartment model based off the Hodgkin-Huxley equations. 


## Installation 

Get this repo from within `MATLAB` using my package manager:

```
% copy and paste this code in your MATLAB prompt
urlwrite('http://srinivas.gs/install.m','install.m'); 
install sg-s/srinivas.gs_mtools % you'll need this
install sg-s/puppeteer % for manipulation
install sg-s/xolotl % C++ code for simulation
install marderlab/eag % this repository 
```

or use git if you plan to develop this further: 

```
git clone https://github.com/sg-s/srinivas.gs_mtools
git clone https://github.com/sg-s/puppeteer
git clone https://github.com/sg-s/xolotl
git clone https://github.com/marderlab/eag
```

Finally, make sure you [configure MATLAB so that it is set up to delete files permanently](https://www.mathworks.com/help/matlab/ref/delete.html). Otherwise you will end up with a very large number of temporary files in your trash!

## Usage

Reproduce the modelling figure by running

```matlab
paper_figures
```

in your `MATLAB` prompt.

## License 

This software is free software. GPL v3. 