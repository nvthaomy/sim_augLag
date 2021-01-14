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

/* ----------------------------------------------------------------------
   Contributing authors: Tanmoy Sanyal, M.Scott Shell, UC Santa Barbara
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_localdensity.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairLOCALDENSITY::PairLOCALDENSITY(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  one_coeff = 1; 	
  single_enable = 0;
  DEBUG = 0; // turn this off if display of parsed and splined arrays is not required	

  // read from file
  nLD = 0;
  nrho = 0;
  rho_min = NULL;
  rho_max = NULL;
  a = NULL;
  b = NULL;
  uppercut = NULL;
  lowercut = NULL;
  frho = NULL;
  rho = NULL;
  
  // splined arrays
  frho_spline = NULL;
  
  // per-atom arrays
  nmax = 0;
  fp = NULL;
  localrho = NULL;  

  // set comm size needed by this Pair

  comm_forward = 1;
  comm_reverse = 1;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairLOCALDENSITY::~PairLOCALDENSITY()
{

  memory->destroy(localrho);
  memory->destroy(fp);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }

  memory->destroy(frho_spline);
  
  memory->destroy(rho_min);  
  memory->destroy(rho_max);
  memory->destroy(delta_rho);	
  memory->destroy(uppercut);
  memory->destroy(lowercut);
  memory->destroy(frho);
  memory->destroy(rho);

  delete [] a;
  delete [] b;
  
  
}

/* ---------------------------------------------------------------------- */

void PairLOCALDENSITY::compute(int eflag, int vflag)
{
  
  int i,j,ii,jj,m,k,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,r,rsq;
  double uLD, phip, evdwl,fpair;
  double p, *coeff;
  int *ilist,*jlist,*numneigh,**firstneigh;

  uLD = 0.0;
  evdwl = 0.0;
  fpair = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  // grow local-density and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(localrho);
    memory->destroy(fp);
    nmax = atom->nmax; 
    memory->create(localrho, nLD, nmax, "pairLD:localrho");
    memory->create(fp, nLD, nmax, "pairLD:fp");
  }

  double **x = atom->x; 
  double **f = atom->f;
  int *type = atom->type; 
  int nlocal = atom->nlocal; 
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density

  if (newton_pair) {
    m = nlocal + atom->nghost;
    for (k = 0; k < nLD; k++) { 
    	for (i = 0; i < m; i++)	 		
		localrho[k][i] = 0.0;
    }	
  } 
  else {
    for (k = 0; k < nLD; k++){
  	for (i = 0; i < nlocal; i++)
  		localrho[k][i] = 0.0;

  	}
   }	

  // localrho = local density at each atom
  // loop over neighbors of my atoms and types of local-densities
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];	

      // calculate distance-squared between i,j atom-types
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;	
      
      // global cutoff check

      if (rsq > cutforcesq) {
	continue;
      }  
      // calculating local densities based on central and neighbor filters

      for (k = 0; k < nLD; k++) {
      	if (a[k][itype] && b[k][jtype])
		localrho[k][i]      += get_phi(-1, rsq, uppercut[k], lowercut[k], 0);
	if (a[k][jtype] && b[k][itype])
		localrho[k][i]      += get_phi(-1, rsq, uppercut[k], lowercut[k], 0);
     }
    }
  }
	
  // communicate and sum densities
  if (newton_pair) comm->reverse_comm_pair(this);

  // u_LD = embedding energy of each atom due to each LD potential type 
  // fp = derivative of embedding energy at each atom for each LD potential type

  for (ii = 0; ii < inum; ii++) {
    for (k = 0; k < nLD; k++) {
    	i = ilist[ii];
	itype = type[i];
	if (a[k][itype]) {
    		p = localrho[k][i] / delta_rho[k] + 1.0;
		m = static_cast<int> (p);
		m = MAX(1,MIN(m, nrho-1));
		p -= m;
		p = MIN(p, 1.0);
		coeff = frho_spline[k][m];
		fp[k][i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    		if (eflag) {
      			uLD = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      			if (eflag_global) eng_vdwl += uLD;
      			if (eflag_atom) eatom[i] += uLD;
       		}	
	}

  }
}

  // communicate derivatives of embedding function and localdensity

  comm->forward_comm_pair(this);
  
  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      // calculate square of distance between i,j atoms

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      // global cutoff check

      if (rsq > cutforcesq) continue;
    	        
      // calculate force between two atoms  
      
      r = sqrt(rsq);
      for (k = 0; k < nLD; k++) {
	phip  +=  get_phi(r, rsq, uppercut[k], lowercut[k], 1); // derivative of localdensity
	fpair += -(a[k][itype]*b[k][jtype]*fp[k][i] + a[k][jtype]*b[k][itype]*fp[k][j]) * phip;
      }	
      fpair *= (1/r); 			
        
      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      if (newton_pair || j < nlocal) {
	f[j][0] -= delx*fpair;
        f[j][1] -= dely*fpair;
        f[j][2] -= delz*fpair;
      }

      if (eflag) evdwl = 0.0; // eng_vwdl has already been completely built
      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  
  if (vflag_fdotr) virial_fdotr_compute();
}


/* ----------------------------------------------------------------------
   calculate phi(r_ij) and its derivative
------------------------------------------------------------------------- */

double PairLOCALDENSITY::get_phi(double r, double rsq, double R2, double R1, int option){

 // option = 0 means calculate phi 
 // option = 1 means calculate derivative of phi
 // if option = 0, then the calling function needs to supply r = -1

 double uppercutsq, lowercutsq, uppercutfourth ,uppercutsixth, cut_ratio, denom;
 double c0, c2, c4, c6;
 double returnval = 0; 

 lowercutsq =  R1*R1;
 uppercutsq =  R2*R2;
 uppercutfourth = uppercutsq*uppercutsq;
 uppercutsixth   = uppercutfourth*uppercutsq;
 cut_ratio = lowercutsq/uppercutsq;
 denom = (1-cut_ratio)*(1-cut_ratio)*(1-cut_ratio);
 
 
 c2 =  (1/uppercutsq) * (6*cut_ratio) / denom;
 c4 = -(1/uppercutfourth) * (3 + 3*cut_ratio) / denom;
 c6 = (1/uppercutsixth) * 2.00 / denom;
 
 if (option == 0) {
  if (rsq <= lowercutsq) returnval = 1.0;
  if (rsq >= uppercutsq) returnval = 0.0;
  if (rsq > lowercutsq && rsq < uppercutsq){
	c0 = (1 - 3*cut_ratio) / denom;
 	returnval = c0 + rsq * (c2 + rsq * (c4 + c6*rsq));
  }
 }

 if (option == 1) {
  if (rsq <= lowercutsq || rsq >= uppercutsq ) returnval = 0.0;
  else {
  	if (r!=-1) returnval = (1/r) * (r + rsq * (2*c2 + rsq * (4*c4 + 6*c6*rsq)));
  }
 }
	
 return returnval;
 
 }

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLOCALDENSITY::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
 
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;


}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLOCALDENSITY::settings(int narg, char **arg)
{
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for all type pairs
   read LD file
------------------------------------------------------------------------- */

void PairLOCALDENSITY::coeff(int narg, char **arg)
{
  int i, j;
  if (!allocated) allocate();

  if (narg != 3) error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  
  // parse LD file

  parse_file(arg[2]);


 // clear setflag since coeff() called once with I,J = * *

  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++)
      setflag[i][j] = 0;

  // set setflag for all i,j type pairs

  int count = 0;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
        setflag[i][j] = 1;
        count++;
      }
    }
  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLOCALDENSITY::init_style()
{
  // spline rho and frho arrays
  // request half neighbor list

  array2spline();
  if (DEBUG)
	display();

  neighbor->request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLOCALDENSITY::init_one(int i, int j)
{
  // single global cutoff = max of all uppercuts read in from LD file

  cutmax = 0.0;
  for (int k = 0; k < nLD; k++)
  	cutmax = MAX(cutmax,uppercut[k]);
    
  cutforcesq = cutmax*cutmax;

  return cutmax;
}

/* ----------------------------------------------------------------------
   read potential values from a single Local Density file
------------------------------------------------------------------------- */

void PairLOCALDENSITY::parse_file(char *filename) {
  	
  int k, n;
  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];

  if (me == 0) {
    fptr = fopen(filename, "r");
    if (fptr == NULL) {
      char str[128];
      sprintf(str,"Cannot open Local Density potential file %s",filename);
      error->one(FLERR,str);
    }
  }


 double *ftmp;	   // temprary variable to extract the complete 2D frho array from file
   
 // broadcast number of LD potentials and number of (rho,frho) pairs
 if (me == 0) {
    
    // first 2 comment lines ignored	
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    
    // extract number of potentials and number of (frho, rho) points
    fgets(line,MAXLINE,fptr);	
    sscanf(line,"%d %d",&nLD, &nrho);
    fgets(line,MAXLINE,fptr);
  }

  MPI_Bcast(&nLD,1,MPI_INT,0,world);
  MPI_Bcast(&nrho,1,MPI_INT,0,world);
  
  // setting up all arrays to be read from files and broadcasted
  memory->create(uppercut, nLD, "pairLD:uppercut");
  memory->create(lowercut, nLD, "pairLD:lowercut");
  memory->create(rho_min,  nLD, "pairLD:rho_min"); 
  memory->create(rho_max,  nLD, "pairLD:rho_max");
  memory->create(delta_rho, nLD,"pairLD:delta_rho");
  memory->create(ftmp, (nrho+1)*nLD, "pairLD:ftmp");
  
  // setting up central and neighbor atom filters		
  memory->create(a, nLD, atom->ntypes , "pairLD:a");
  memory->create(b, nLD, atom->ntypes, "pairLD:b"); 	
  if (me == 0) {
  	for (n = 0; n < atom->ntypes; n++){
		for (k = 0; k < nLD; k++) {
			a[k][n] = 0;
			b[k][n] = 0;
		}
  	}
  }	
  
 // read file block by block
  
  if (me == 0) {
  	for (k = 0; k < nLD; k++) {
	
		// parse upper and lower cut values	
		if (fgets(line,MAXLINE,fptr)==NULL) break;
		sscanf(line, "%lf %lf", &lowercut[k], &uppercut[k]);
	
		// parse and broadcast central atom filter
		fgets(line, MAXLINE, fptr);
		char *tmp = strtok(line, " /t/n/r/f");
		while (tmp != NULL) {
			a[k][atoi(tmp)-1] = 1;
			tmp = strtok(NULL, " /t/n/r/f");
		}
                MPI_Bcast(&a[k][0], atom->ntypes, MPI_INT, 0, world);
		
		// parse neighbor atom filter
		fgets(line, MAXLINE, fptr);
		tmp = strtok(line, " /t/n/r/f");
		while (tmp != NULL) {			
			b[k][atoi(tmp)-1] = 1;
			tmp = strtok(NULL, " /t/n/r/f");
		}
 		MPI_Bcast(&b[k][0], atom->ntypes, MPI_INT, 0, world);
	
		// parse min, max and delta rho values
		fgets(line, MAXLINE, fptr);
		sscanf(line, "%lf %lf %lf", &rho_min[k], &rho_max[k], &delta_rho[k]);
	 	
		// parse tabulated frho values from each line into temporary array
		for (n = 0; n < nrho+1; n++) {
			fgets(line,MAXLINE,fptr);
			sscanf(line, "%lf", &ftmp[k*(nrho+1) + n]);
		}
		
		// ignore blank line at the end of every block
		fgets(line,MAXLINE,fptr);
		

  	}
  }
	
  // Broadcast all parsed arrays	
  MPI_Bcast(&lowercut[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&uppercut[0], nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&rho_min[0],  nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&rho_max[0],  nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&delta_rho[0],nLD, MPI_DOUBLE, 0, world);
  MPI_Bcast(&ftmp[0], nLD*(nrho+1), MPI_DOUBLE,0, world);


  if (me == 0) fclose(fptr);

  // set up rho and frho arrays
  memory->create(rho, nLD, nrho+1, "pairLD:rho");
  memory->create(frho,nLD, nrho+1, "pairLD:frho"); 
  
  for (k = 0; k < nLD; k++) {
	for (n = 0; n < nrho+1; n++) {
		rho[k][n] = rho_min[k] + n*delta_rho[k];
		frho[k][n] = ftmp[k*(nrho+1) + n];
	}
 }

}

 
/*--------------------------------------------------------------------
   Spline the array frho read in from the file to create
   frho_spline
---------------------------------------------------------------------- */

void PairLOCALDENSITY::array2spline() {


  memory->destroy(frho_spline);
  memory->create(frho_spline,nLD,nrho+1,7,"pairLD:frho_spline");

  for (int k = 0; k < nLD; k++)
    interpolate(nrho,delta_rho[k],frho[k],frho_spline[k]);

}

/* ---------------------------------------------------------------------- */

void PairLOCALDENSITY::interpolate(int n, double delta, double *f, double **spline) {
  for (int m = 1; m <= n; m++) spline[m][6] = f[m];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5 * (spline[3][6]-spline[1][6]);
  spline[n-1][5] = 0.5 * (spline[n][6]-spline[n-2][6]);
  spline[n][5] = spline[n][6] - spline[n-1][6];

  for (int m = 3; m <= n-2; m++)
    spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) +
                    8.0*(spline[m+1][6]-spline[m-1][6])) / 12.0;

  for (int m = 1; m <= n-1; m++) {
    spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) -
      2.0*spline[m][5] - spline[m+1][5];
    spline[m][3] = spline[m][5] + spline[m+1][5] -
      2.0*(spline[m+1][6]-spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m <= n; m++) {
    spline[m][2] = spline[m][5]/delta;
    spline[m][1] = 2.0*spline[m][4]/delta;
    spline[m][0] = 3.0*spline[m][3]/delta;
  }
}

/* ----------------------------------------------------------------------
   communication routines
------------------------------------------------------------------------- */


int PairLOCALDENSITY::pack_comm(int n, int *list, double **buf, int pbc_flag, int *pbc) {
  int i,j,k;
  int l, m; 	

  l = 0;
  m = 0;
  for (k = 0; k < nLD; k++) {
    for (i = 0; i < n; i++) {
		j = list[i]; 
    		buf[l++][m++] = fp[k][j];
    }		
  }
  
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairLOCALDENSITY::unpack_comm(int n, int first, double **buf) {

  int i,k,l,m,last;
  
  l = 0;
  m = 0;
  last = first + n;
  for (k = 0; k < nLD; k++) {
	for (i = first; i < last; i++)
		fp[k][i] = buf[l++][m++];
 }		
}

/* ---------------------------------------------------------------------- */

int PairLOCALDENSITY::pack_reverse_comm(int n, int first, double **buf) {

  int i,k,l,m,last;

  l = 0;
  m = 0;
  last = first + n;
  for (k = 0; k < nLD; k++){
  	for (i = first; i < last; i++) 
		buf[l++][m++] = rho[k][i];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairLOCALDENSITY::unpack_reverse_comm(int n, int *list, double **buf) {

  int i,j,k;
  int l, m;

  l = 0;
  m = 0;
  for (k = 0; k < nLD; k++){
  	for (i = 0; i < n; i++) {
    		j = list[i];
    		rho[k][j] += buf[l++][m++];
    	}	
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairLOCALDENSITY::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += 2 * (nmax*nLD) * sizeof(double);
  return bytes;
}


/* ----------------------------------------------------------------------
   displaying parsed and splined data (only for debugging)
------------------------------------------------------------------------- */

void PairLOCALDENSITY::display(){
 int i,j,k;
 printf("\n LD-FILE PARSED, FRHO ARRAY SPLINED\n");
 printf("\ DATA FLOW TILL NOW (FOR DEBUGGING PURPOSES)\n");
 printf("\n\n -------------------------------------------\n\n");
 printf("\n N_LD = %d\n", nLD);
 printf("\n N_RHO = %d\n", nrho);
 
 printf("\n UPPERCUT\n");	
 for (k = 0; k < nLD; k++)
	printf("\%lf\t", uppercut[k]);
 printf("\n LOWERCUT\n");	
 for (k = 0; k < nLD; k++)
	printf("\%lf\t", lowercut[k]);

 printf("\n CENTRAL ATOM FILTER \n");
 for (k = 0; k < nLD; k++){
	for (i = 0; i < atom->ntypes; i++)
		printf("\%d\t", a[k][i]);
	printf("\n");
 }

 printf("\n NEIGHBOR ATOM FILTER \n");
 for (k = 0; k < nLD; k++){
	for (i = 0; i < atom->ntypes; i++)
		printf("\%d\t", b[k][i]);
	printf("\n");
 }

 printf("\n RHO_DATA\n");
 printf("\n RHO_MIN\tRHO_MAX\tDELTA_RHO\n");
 for (k = 0; k < nLD; k++) 
	printf("\n%lf\t%lf\t%lf\n", rho_min[k], rho_max[k], delta_rho[k]);

 printf("\n\n FRHO AS READ FROM FILE\n\n");
 for (i = 0; i < nrho+1; i++){
	for (k = 0; k < nLD; k++)
		printf("%lf\t", frho[k][i]);
	printf("\n");
 }

 printf("\n\n FRHO SPLINED\n\n");
 for (i = 0; i < nrho+1; i++){
	for (k = 0; k < nLD; k++) {
		for (j = 0; j < 7; j++)
			printf("%\lf\t", frho_spline[k][i][j]);
		printf("\n");
	}
	printf("\n\n------ potential %d -----------\n\n", i+1);
  }


 }
 

