# LAMMPS-SRD

This repository contains LAMMPS with the improvements from the [Palmer group](https://github.com/palmergroup/mpcd_polymer_colloid/), providing an implementation of Multi-Particle Collision Dynamics (MPCD) as Stochastic Rotation Dynamics (SRD) or with the Andersen Thermostat variant. That code was later ported to a newer LAMMPS version by [hirotori](https://github.com/hirotori/lammps/). The original SRD module has been extended by the Palmer group with

- an Andersen thermostat coupling scheme,
- collision coupling of solvent and solutes for the SRD scheme, i.e., including colloids in the collision and thermostat steps.

This repository contributes by

- fixing cell handling during shear flow,
- adding a profile-unbiased thermostat.
- On top, an implementation of Random MPC is contained.

See [this issue](https://github.com/lammps/lammps/pull/4590) for further background.

The modified SRD module traces back to:

> Chen, R., Poling-Skutvik, R., Nikoubashman, A., Howard, M. P., Conrad, J. C., & Palmer, J. C. (2018).  
> *Coupling of nanoparticle dynamics to polymer center-of-mass motion in semidilute polymer solutions.*  
> *Macromolecules*, 51(5), 1865–1872.

and was updated to `LAMMPS (2 Aug 2023 - Development)`.

## Usage
You can find an example input file in the [Palmer group repository](https://github.com/palmergroup/mpcd_polymer_colloid/blob/master/MPCD%20input%20files/in.mpcd_polycoll/). Make sure to name the group of particles you want to couple in the collision step `small`.

You can turn on the profile-unbiased thermostat with the `unbiased yes` keyword (default `no`).
```
group small type 1
fix 1 solvent srd 100 NULL 1.0 1.0 49894 &
    tstat yes &
    unbiased yes &
    shift yes 42 &
    rescale no &
    small ${mass_small}
```

For Random MPC, the arguments are: $N$, the period of the collision step, $T$, the temperature of the virtual SRD particles, $m_\mathrm{SRD}$, their mass, and a seed.
```
fix 1 all rsrd 100 1.0 10 1.0 42
```

## Known limitations
- The code has been tested, but I do not guarantee its correctness.
- I did not test the Andersen thermostat, but still modified its code to support the `unbiased true` functionality.
- Using the profile-biased thermostat (`unbiased no`) under shear flow may lead to wrong temperatures for the colloids. You can keep track of that with a `compute t1 small temp/deform` or `compute t2 small temp/partial 0 1 1`. More precisely, this happens as soon as `fix deform` is present. This has to do with the fact that colloids may dominate the center of mass velocity due to their mass, and the way the rescaling thermostat works.

## Build instructions
Clone the repository and move into the source folder
```
git clone https://github.com/taenzel/lammps-srd.git
cd lammps-srd/src
```
Then, build with `make`. More information can be found in the [LAMMPS documentation](https://docs.lammps.org/Build_make.html).
```
make yes-SRD
make mpi -j <number of procs>
```

## Citation
Please cite the appropriate LAMMPS references, the corresponding article from the Palmer group, their github repository, this repository if you used its new features, and whatever else you deem appropriate. For Random MPC, cite the related literature and this repository.

## Collision step and thermostat
The [LAMMPS documentation](https://docs.lammps.org/fix_srd.html) and the [research articles](https://docs.lammps.org/fix_srd.html#default) provide some information on their SRD collision and thermostat scheme. By the way, additional thermostats may be found in [this code](https://matsci.org/t/inducing-flow-in-an-srd-fluid/25810/7). The following provides more details about the SRD thermostat and what is meant with profile-biased and -unbiased.

In each SRD collision cell $c$,

- $N_c$ is the number of particles (solvent and `small` colloids) in the cell,
- $M_c = \sum_{i \in c} m_i$ is the total mass in the cell,
- $`\mathbf{v}_c = \frac{1}{M_c}\sum_{i \in c} m_i \mathbf{v}_i`$ is the cell center-of-mass velocity,
- $T_c$ is the cell temperature measured relative to $\mathbf{v}_c$,
- $T_0$ is the target temperature
- $\zeta_c = \sqrt{\frac{T_0}{T_c}}$ is the rescaling factor of the thermostat

### Profile-biased form (`unbiased no`)

In the standard SRD collision and thermostat step, the post-update velocity can be written as

$$
\mathbf{v}_i' = \mathbf{u}(\mathbf{r}_c) + \zeta_c \left(\mathbf{v}_i - \mathbf{v}_c\right),
$$

where $\mathbf{u}(\mathbf{r}_c)$ is the imposed streaming velocity at the cell center $\mathbf{r}_c$. For simple shear flow with flow in $x$ and gradient in $y$,

$$
\mathbf{u}(\mathbf{r}_c) = \dot{\gamma} y_c \hat{\mathbf{x}}.
$$

The velocity profile is thus not a consequence of the Lees-Edwards boudary conditions of `fix deform`, but artificially imposed during the thermostatting.

### Profile-unbiased form (`unbiased yes`)

In the profile-unbiased variant (`unbiased yes`), the post-update velocity is

$$
\mathbf{v}_i' = \mathbf{v}_c + \zeta_c \left(\mathbf{v}_i - \mathbf{v}_c\right).
$$

When a `fix deform` is given, the velocity profile may then build up due to the boundary conditions alone.

## Random MPC
Random MPC is a method that is close to MPCD with SRD collision step, but without hydrodynamic interactions. No MPC particles are needed. Instead, the idea is to couple each solute particle to an effective solvent momentum drawn from a Maxwell-Boltzmann distribution. Shear flow can be accounted for, a linear shear profile is assumed.

You can find more information for example in 
> Ripoll, M., Winkler, R. G., & Gompper, G. (2007).
> *Hydrodynamic screening of star polymers in shear flow.*
> *The European Physical Journal E*, 23(4), 349-354.

## Practical advice
The implementation in HOOMD is probably preferable, as soon as it is available.
