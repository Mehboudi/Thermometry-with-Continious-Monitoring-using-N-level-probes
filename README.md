Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

# Thermometry with Continious Monitoring using N level probes
 These codes are attached to our paper "[Optimal limits of continuously monitored thermometers and their Hamiltonian structure
](https://arxiv.org/abs/2408.01313v1)" [arXiv:2408.01313]. They help finding the optimal energy structure of N-level probes that continiously monitor the temperature of a fermionic/bosonic bath


# HERE are further instructions:


These documents find the optimal Hamiltonian structure for thermometry of a continiously monitored D-level system.
The bath can be fermionic or bosonic, each of which require a different code.

In the fermionic case, the structure seems to be an effective 3-level system. Thus we also take this as an ansatz.

You run the code from 
run_ham_structure
This will do the following:

I) fermion=1; means that we run for fermionic case (with flat SD)
If you set the parameter ferm_3lev=0, the code finds global optimal.
This turns out to be very close to an effective 2-level system, with one gap equal to zero.

If you set the parameter ferm_3lev=1, the code uses the 2-level ansatz. It leaves the degeneracy and gap
as the open parameters, and finds the best ones.





II) fermion=0; means that we run for bosonic case (with ohmicity s_o, you shold chose it within the correponding code)
If you set the parameter ferm_3lev=0, the code finds global optimal.
This turns out to be very close to an effective 2-level system, with one gap equal to zero.

If you set the parameter bose_3lev=1, the code uses the 2-level ansatz. It leaves the degeneracy
of each level as an open parameter, and finds the best one.

## Robustness to energy structure
In a new update, we also include some codes that consider robustness of the FI to the Hamiltonian structure. In particular, instead of the optimal structure, which is two levels with energies that are sharply equal, we consider two -levels with energies that can fluctuate around the optimal values. These are modeled with some Gaussian distribution of the energies around the optimal values.

You can run the robustness only for the fermionic case, but bosonic should be similar.
The code is called run_robustness_ham_structure. The main new parameter to change is the standard deviation of the energy levels, sigma. 
