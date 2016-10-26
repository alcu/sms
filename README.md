# Simultaneous multislice fMRI

This repository contains MATLAB scripts for simultaneous multislice (SMS)
fMRI reconstruction with concentric ring k-space trajectories. Coil
compression is included as an option in each reconstruction.

## Setup

Run the startup.m script to set the MATLAB path.

## SMS GRAPPA

To do an SMS GRAPPA recon, run the main_grappa.m script in the recon/
directory. The variable xc is the reconstructed image.

## SMS SENSE

To do an SMS SENSE recon, run the main_sense.m script in the recon/
directory. The variable xc is the reconstructed image.

## Non-SMS

To do a non-SMS recon, run the main_ss.m script in the recon/
directory. The variable xc is the reconstructed image.

## Concentric ring trajectory

The concen_ring_traj/ directory contains files for generating concentric
ring trajectories. The main script to run is concen_rings2.m. The file
data/concen_ring/concen_rings2.mat contains saved variables from running
concen_rings2.m, and is used in the main recon scripts.

## Generating sensitivity maps with ESPIRiT

Run the demo_nuESPIRiT_alanchu.m file in the ESPIRiT/ directory to generate
sensitivity maps. demo_nuESPIRiT_alanchu.m is a modified version of
Lustig's demo_nuESPIRiT.m (see below). The variable c has the sensitivity
maps.

## References

The irt/ directory contains Jeff Fessler's image reconstruction toolbox
from http://web.eecs.umich.edu/~fessler/code/index.html. The relevant
publication is at http://dx.doi.org/10.1109/TSP.2002.807005.

concen_rings2.m is based on Jim Pipe's spiral waveform generation code from
http://www.ismrm.org/mri_unbound/. The relevant publication is at
http://dx.doi.org/10.1002/mrm.24675.

The ESPIRiT/ directory contains Miki Lustig's ESPIRiT code from
http://www.eecs.berkeley.edu/~mlustig/Software.html (plus a few minor bug
fixes). The relevant publication is at http://dx.doi.org/10.1002/mrm.24751.
