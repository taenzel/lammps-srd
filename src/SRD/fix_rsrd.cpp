/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Victor Tänzel (UFR)
------------------------------------------------------------------------- */

#include "fix_rsrd.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_deform.h"
#include "force.h"
#include "random_mars.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_fix_rsrd[] =
    "fix rsrd command: doi:10.1140/epje/i2006-10220-0\n\n"
    "@Article{Ripoll07,\n"
    " author = {Ripoll, M and Winkler, RG and Gompper, G},\n"
    " title =  {Hydrodynamic screening of star polymers in shear flow},\n"
    " journal = {The European Physical Journal E},\n"
    " year =    {2007},\n"
    " volume =  {23},\n"
    " pages =   {349--354}\n"
    " publisher = {Springer}\n"
    "}\n\n";


FixRSRD::FixRSRD(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_rsrd);

  if (narg < 8) error->all(FLERR, "Too few arguments.");

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  temperature_srd = utils::numeric(FLERR, arg[4], false, lmp);
  srd_per_cell = utils::numeric(FLERR, arg[5], false, lmp);
  mass_srd = utils::numeric(FLERR, arg[6], false, lmp);
  seed = utils::inumeric(FLERR, arg[7], false, lmp);
  
  // error check
  if (nevery <= 0) error->all(FLERR, "Requires non-negative nevery");
  if (temperature_srd <= 0.0) error->all(FLERR, "Requires non-negative temperature");
  if (srd_per_cell <= 0) error->all(FLERR, "Requires non-negative srd density");
  if (mass_srd <= 0) error->all(FLERR, "Requires non-negative srd mass");
  if (seed <= 0) error->all(FLERR, "Requires non-negative seed");

  // initialize Marsaglia RNG with processor-unique seed
  me = comm->me;
  random = new RanMars(lmp, seed + me);

  lebcflag = domain->deform_vremap;
}

FixRSRD::~FixRSRD()
{
  delete random;
}

int FixRSRD::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

void FixRSRD::init() 
{
  lebcflag = domain->deform_vremap;
}

void FixRSRD::setup(int vflag)
{
  // set up mass of group particles
  // assumes same mass for all particles
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass) {
          mass_particle = rmass[i];
          break;
      } else {
          mass_particle = mass[type[i]];
          break;
      }
    }
  }
  mass_solvent = srd_per_cell*mass_srd;
  inv_mass_part_solv = 1.0 / (mass_particle + mass_solvent);
}

void FixRSRD::post_force(int /*vflag*/) 
{
  if (update->ntimestep % nevery == 0) {
    reset_velocities();
  }
}

void FixRSRD::reset_velocities()
{
  int j;
  int irandom, axis, sign;

  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  int nlocal = atom->nlocal;

  double vave[3];
  double psolvent[3];
  double vstream[3];
  double u[3];
  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  // The standard deviation of the random momentum
  sigma = sqrt(srd_per_cell * mass_srd * force->boltz * temperature_srd);

  // Prepare that atom positions are in [0, 1] to compute vstream
  if (lebcflag) domain->x2lamda(nlocal);

  for (j = 0; j < nlocal; j++) {
    
    irandom = static_cast<int>(6.0 * random->uniform());
    sign = irandom % 2;
    axis = irandom / 2;

    // Draw a random solvent momentum to which we will couple
    psolvent[0] = random->gaussian() * sigma;
    psolvent[1] = random->gaussian() * sigma;
    psolvent[2] = random->gaussian() * sigma;

    // Calculate streaming velocity, add to momentum to random momentum
    if (lebcflag) {
      vstream[0] = h_rate[0] * x[j][0] + h_rate[5] * x[j][1] + h_rate[4] * x[j][2] + h_ratelo[0];
      vstream[1] = h_rate[1] * x[j][1] + h_rate[3] * x[j][2] + h_ratelo[1];
      vstream[2] = h_rate[2] * x[j][2] + h_ratelo[2];

      psolvent[0] += srd_per_cell * mass_srd * vstream[0];
      psolvent[1] += srd_per_cell * mass_srd * vstream[1];
      psolvent[2] += srd_per_cell * mass_srd * vstream[2];
    }

    // Calculate the randomSRD center of mass velocity
    vave[0] = (mass_particle*v[j][0] + psolvent[0]) *  inv_mass_part_solv;
    vave[1] = (mass_particle*v[j][1] + psolvent[1]) *  inv_mass_part_solv;
    vave[2] = (mass_particle*v[j][2] + psolvent[2]) *  inv_mass_part_solv;

    // Perform the collision step as in fix_srd.cpp
    if (axis == 0) {
      u[0] = v[j][0] - vave[0];
      u[1] = sign ? v[j][2] - vave[2] : vave[2] - v[j][2];
      u[2] = sign ? vave[1] - v[j][1] : v[j][1] - vave[1];
    } else if (axis == 1) {
      u[1] = v[j][1] - vave[1];
      u[0] = sign ? v[j][2] - vave[2] : vave[2] - v[j][2];
      u[2] = sign ? vave[0] - v[j][0] : v[j][0] - vave[0];
    } else {
      u[2] = v[j][2] - vave[2];
      u[1] = sign ? v[j][0] - vave[0] : vave[0] - v[j][0];
      u[0] = sign ? vave[1] - v[j][1] : v[j][1] - vave[1];
    }

    v[j][0] = u[0] + vave[0];
    v[j][1] = u[1] + vave[1];
    v[j][2] = u[2] + vave[2];
  }
  
  if (lebcflag) domain->lamda2x(nlocal);

}


