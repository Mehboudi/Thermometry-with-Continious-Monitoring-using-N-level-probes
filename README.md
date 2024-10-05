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
