#include "calcsurf.h"
#include "dssp\\vector_dssp.h"
#include "dssp\\calcaccsurf.h"


NEIGHRES NeighbourRes[NMAX + 1];
long LastNeighbourRes;

void Listentry(CAS* cas, double *xx, double *neighbour)
{
    double sqdist=Distsq(xx, neighbour);
    if(sqdist < (xx[3]+neighbour[3])*(xx[3]+neighbour[3]))
	CASAddNeigbourAtom(cas,
			   neighbour[0],neighbour[1], neighbour[2],
			   neighbour[3],sqdist);
}  /* Listentry */

int InBox(double *v, double *vmin, double *vmax)
{
  if (v[0] - v[3] < vmax[0] && v[1] - v[3] < vmax[1] &&
      v[2] - v[3] < vmax[2] && v[0] + v[3] > vmin[0] &&
      v[1] + v[3] > vmin[1] && v[2] + v[3] > vmin[2])
    return 1;
  else
    return 0;
}





void MinMax(double *v, double *vmin, double *vmax)
{
  long i;

  for (i = 0; i <= 2; i++) {
    if (v[i] - v[3] < vmin[i])
      vmin[i] = v[i] - v[3];
    if (v[i] + v[3] > vmax[i])
      vmax[i] = v[i] + v[3];
  }
}  
void Atom2Pointer4(const Atom& atom,double *v)
{
	v[0]=atom.pt.x;
	v[1]=atom.pt.y;
	v[2]=atom.pt.z;
	v[3]=atom.radium;
}
//void Flagaccess(int order)
//{
//    long i, k;
//    double f;
//    long FORLIM;
//    Backbone *WITH;
//    long FORLIM1;
//    CAS *cas= CASCreate(order);
//    AddToAllAtomRadii(RWATER); /* add the water radius :-) */
//    CalcResidueCenters();
//    FORLIM = lchain;
//    for (i = 1; i <= FORLIM; i++) {
//	if (chain[i].aa != '!') {
//	    WITH = &chain[i];
//	    FindNeighbourRes(WITH->boxmin, WITH->boxmax);
//	    f = Surface(cas,WITH->n) + Surface(cas, WITH->ca) +
//		Surface(cas, WITH->c) + Surface(cas, WITH->o);
//	    if (WITH->nsideatoms > 0) {
//		FORLIM1 = WITH->nsideatoms;
//		for (k = 0; k < FORLIM1; k++)
//		    f += Surface(cas, sidechain[WITH->atompointer + k]);
//	    }
//	    WITH->access = (long)floor(f + 0.5);
//	}
//    }
//    AddToAllAtomRadii(-RWATER);   /* remove the water :-) */
//    CASDelete(cas);
//}  /* Flagaccess */
//
