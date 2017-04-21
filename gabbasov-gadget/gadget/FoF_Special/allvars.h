
#include "proto.h"


#define  MAX_NGB          20000

#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#define  GRAVITY     6.672e-8
#define  HUBBLE          3.2407789e-18	/* in h/sec */

extern  double Time;


/*  particle numbers */

extern int     NumPart;
extern int     N_DM;
extern double  Omega0;
extern double  BoxSize;
extern double  BoxHalf;
extern double  M_DM;

extern int    MaxNodes;


extern double GroupMinMass;
extern double SearchRadius;
extern double LinkL;


extern int *Head, *Next;		/* for link-lists of groups */
extern int *Tail, *Len;
extern int *Nearest;

extern int *NextFinal, *IdSort;
 
extern int   *GroupLen, *GroupTag;
extern float *GroupMass;
extern int   Ngroups;
extern int   NgroupsAboveMinLen;

extern int Snapshot;
extern int Files;




/* Quantities for all particles */

extern struct particle_data 
{
  float  Pos[3];
#ifdef VELOCITYDISPERSION
  float  Vel[3];
#endif
  int    FileIndex;
  int    Type;
  int    Nextnode;
  float  Mass;
  float  Sfr;
  float  Metallicity;
} *P, *PP;



extern struct NODE
{
  union
  {
    int suns[8];		/* ancestor nodes */
    struct 
    {
      float len; 	        /* sidelength of treenode */
      float center[3];	        /* geometrical center */
      int   sibling;
      int   nextnode;
    } d;
  } u;
}
*Nodes, *Nodes_base;


extern float *R2list;

extern int   *Ngblist;




extern struct io_header
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[96];		/* fills to 256 Bytes */
}
header;




















