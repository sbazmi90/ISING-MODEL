# ISING-MODEL


## Magnetic Materials:

Some atoms have an unpaired electron possessing intrinsic angular momentum such as Fe (iron)

According to the quantum mechanics, the atomic spin can only hold two valus of opposite sign. s and -s. Those also
can be considered as +1 and -1, respectively.

If t spin of material are aligned, the sum of all atomic spins will produce a net magnetic moment of the material: $\langle M \rangle$.

## Simulation:

Real magnetic materials have $\mathcal{O}(10^{23})$ spins. However, infinite systems can mimick with periodic boundary conditions.

Apparently, if there are more spins pointing up than down (or vice versa), the lattice will have
a net magnetization (behaves like a magnet).

In totally disordered state, spins are arranged randomly means no net magnetization.

Instead of calculating averages by averaging over every state, we can select a subset of states. 
States must be chosen without bias:
– Markov chain
– chosen through “random walk”
– satisfy detail balance criterion

Average properties of the system can be estimated as the average from the sampled configurations.
