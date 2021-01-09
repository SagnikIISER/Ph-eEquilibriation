# Ph-eEquilibriation
This repository contains the files and codes related to the phonon-electron equilibriation problem that I am pursuing presently as MS Thesis at TIFR, Colaba with Dr Rajdeep Sensarma in the COViD-19 affected Fall-2020, Spring-2021 semesters. We start with a lattice of coupled Electrons and Phonons , with some of the electrons initially excited to high energy states and study how it thermalises using Non-Equilibrium Field Theory.

## Proposal
Consider an interacting electron-phonon system with some of the electrons initially in very high momentum states. We are interested in studying various aspects of how this system equilibrates.

This models various physical systems. For example, when a detector receives a charged particle, it would knock some of the electrons off to very high momentum states. The electrons interacting with each-other, in turn would initiate a cascade, producing a bunch of hot electrons in a crystal. Again in pump-probe experiments, electrons in a crystal are essentially excited to high energy states by shining Lasers. Various results of these experiments are dependent on how these systems thermalize by redistributing the energy.

In general, the electrons would thermalize between themselves quickly at a high temperature and then gradually transfer their energy to the phonons. So at some point there would be an electron temperature and a phonon temperature, and in between there will be a non-equilibrium phenomena. Finally both the temperature are expected to reach the same value thermodynamically. The aim is to study how this process happens.

There are various physical questions of interest. For e.g. how does the temperature rise depend on detector size. For a thermodynamic system, no temperature rise would occur. However, all the basic detectors, (CCDs etc.) are of finite size and gives a finite temperature rise. We are primarily interested in the phonon equilibration as often the measurements are done on them and they hardly thermalizes instantaneously. The interesting questions are what is the time scale of equilibration? Do every phonon mode equilibrates in the same time scale? How much is the efficiency of energy transfer at finite time scales? Etc.

We wish to study some of these with the Swinger-Keldysh formulation of non-equilibrium field theory. The idea is to solve for the time dependant interacting Greens functions from the Dyson Equation and obtain the time dependence of various quantities from them, which can then be used to address various aforementioned questions we have set out to delve in.

## Files
  1. **DysonIteration.tex**: The .tex for Scribbled Ntotes. Contains the derivation for algorithim and pseudocodes for iterating the Dyson Equation for the Bosonic Sector by Euler method and Self-Consistently.
  2. **RajdeepEuler.c**: The code for iterating Dyson Equation by Euler Method.
  3. **RajdeepSCF.c**: The code for iterating Dyson Equation by Self Consistent Method.
  4. **green's.h**: Header file containing the function call for Self Energies, Bare Green's Function and Baths.
