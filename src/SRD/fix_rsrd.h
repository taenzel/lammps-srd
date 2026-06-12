/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(rsrd,FixRSRD);
// clang-format on
#else

#ifndef LMP_FIX_RSRD_H
#define LMP_FIX_RSRD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRSRD : public Fix {
 public:
  FixRSRD(class LAMMPS *, int, char **);
  ~FixRSRD() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;

 private:
  int me;

  // input
  int seed;  // seed for the RNG of the random momentum
  double temperature_srd;  // temperature
  double srd_per_cell;   // number of SRD particles the big particle couples to
  double mass_srd;  // mass of the imaginary SRD particles

  int lebcflag;  // 1 if fix deform remap v present, 0 else
  double mass_particle;  // mass of the big particle
  double mass_solvent;  // mass of the solvent a big particle couples to
  double inv_mass_part_solv;  // 1 / sum of particle and solvent mass
  
  double sigma;

  class RanMars *random;

  void reset_velocities();

};

}    // namespace LAMMPS_NS

#endif
#endif
