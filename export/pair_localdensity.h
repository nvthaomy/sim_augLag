/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(localdensity,PairLOCALDENSITY)

#else

#ifndef LMP_PAIR_LOCALDENSITY_H
#define LMP_PAIR_LOCALDENSITY_H

#include "stdio.h"
#include "pair.h"

namespace LAMMPS_NS {

class PairLOCALDENSITY : public Pair {
 public:

  // public variables 

  double cutmax;

  // data parsed from input file

  int nLD, nrho;	
  int **a, **b; // central and neighbor atom filters respectively
  double *uppercut, *lowercut; 
  double **frho, **rho, *rho_min, *rho_max, *delta_rho;

  // potentials in spline form used for force computation

  double ***frho_spline;

  PairLOCALDENSITY(class LAMMPS *);
  virtual ~PairLOCALDENSITY();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);

  virtual int pack_comm(int, int *, double **, int, int *);
  virtual void unpack_comm(int, int, double **);
  int pack_reverse_comm(int, int, double **);
  void unpack_reverse_comm(int, int *, double **);
  double memory_usage();

 protected:
  int nmax;                   // allocated size of per-atom arrays
  double cutforcesq;          // square of global upper cutoff

  // per-atom arrays

  double **localrho;     // local density
  double **fp;           // derivative of embedding function 

  void allocate();     
  void parse_file(char *);	
  void array2spline();
  void interpolate(int, double, double *, double **);
  double get_phi(double, double, double, double, int);
  void display(); // for debugging
  int DEBUG;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Cannot open Local Density potential file %s

The specified EAM potential file cannot be opened.  Check that the
path and name are correct.

*/
