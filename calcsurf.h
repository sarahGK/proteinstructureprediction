#ifndef CALCSURF_H
#define CALCSURF_H

#include <vector>
#include <string>

using namespace std;

//define radius for atoms from DSSP package
#define RN              1.65
#define RCA             1.87
#define RC              1.76
#define RO              1.4
#define RSIDEATOM       1.8
#define RWATER          1.4

#define YVERTEX         0.8506508
#define ZVERTEX         0.5257311

#define NMAX            20000

typedef double Vector[3];
typedef double VectorRad[4];

typedef struct Polyeder {
    Vector *p;
    double *wp;
    char *accept;
    int np;
    int iter;
} Polyeder;
typedef struct {
    Polyeder* polyeder;
    VectorRad atom;
    VectorRad *neighbours;
    int nneighbours;
    int maxneighbours;
    int wantdots;
} CAS;

typedef struct
{
    int chain;
    long res;
} NEIGHRES;

#include "rd_pdb.h"

extern NEIGHRES NeighbourRes[NMAX + 1];
extern long LastNeighbourRes;

Polyeder * PolyederCreate(long order);
void PolyederInit(long order, Polyeder *polyed);


CAS* CASCreate(int order);

void Listentry(CAS* cas, double *xx, double *neighbour);
int InBox(double *v, double *vmin, double *vmax);
double Surface(CAS* cas,double *xatom);
void MinMax(double *v, double *vmin, double *vmax);


#endif
