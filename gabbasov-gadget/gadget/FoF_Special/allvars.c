
#include "allvars.h"



double Time;


/*  particle numbers */

int NumPart;
int N_DM;
double Omega0;
double BoxSize;
double BoxHalf;

double M_DM;


int MaxNodes;


double GroupMinMass;
double SearchRadius;
double LinkL;


int *Head, *Next;		/* for link-lists of groups */
int *Tail, *Len;
int *Nearest;

int *NextFinal, *IdSort;

int *GroupLen, *GroupTag;
float *GroupMass;
int Ngroups;
int NgroupsAboveMinLen;

int Snapshot;
int Files;


/* Quantities for all particles */

struct particle_data *P, *PP;

int *GridNext;

struct NODE *Nodes, *Nodes_base;


float *R2list;
int *Ngblist;


struct io_header header;
