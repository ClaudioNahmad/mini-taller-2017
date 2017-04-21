#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/* PP[0] = P[1] !
 */



/* index convention for accessing tree nodes:
   the indices   0...NumPart  reference single particles
   the indices   All.MaxPart.... All.MaxPart+nodes-1   reference tree nodes
 
   `Nodes_base' points to the first tree node
   `nodes' is shifted, such that nodes[All.MaxPart] gives the first tree node
*/



static int pmin;
static double r2min;

static int last;		/* auxialiary variable used to set-up non-recursive walk */

double part_dens;

int force_treebuild(void)
{
  int i, j, subnode, numnodes;
  int nfree, th, nn, startnode;
  double xmin[3], xmax[3], len, center[3];
  double Center[3], Len;
  double lenhalf;
  struct NODE *nfreep;

  startnode = 0;
  numnodes = 0;
  Nodes = Nodes_base - NumPart;


/* find enclosing rectangle */
  for(j = 0; j < 3; j++)
    xmin[j] = xmax[j] = PP[0].Pos[j];

  for(i = 1; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      {
	if(PP[i].Pos[j] > xmax[j])
	  xmax[j] = PP[i].Pos[j];
	if(PP[i].Pos[j] < xmin[j])
	  xmin[j] = PP[i].Pos[j];
      }

  part_dens = NumPart / ((xmax[0] - xmin[0]) * (xmax[1] - xmin[1]) * (xmax[2] - xmin[2]));

  /* determine maxmimum extension */
  for(j = 1, len = xmax[0] - xmin[0]; j < 3; j++)
    if((xmax[j] - xmin[j]) > len)
      len = xmax[j] - xmin[j];

  len *= 1.00001;

  for(j = 0; j < 3; j++)
    Center[j] = 0.5 * (xmax[j] + xmin[j]);
  Len = len;



  /* create a root node and insert first particle as its leaf */
  nfree = NumPart + startnode;	/* index */
  nfreep = &Nodes[nfree];	/* select first node */

  for(i = 0; i < 8; i++)
    nfreep->u.suns[i] = -1;

  numnodes++;
  nfree++;
  nfreep++;

  for(i = 0; i < NumPart; i++)	/* insert all  particles */
    {
      if(PP[i].Type == 1)
	{
	  th = NumPart + startnode;	/* select index of first node in tree */

	  len = Len;
	  lenhalf = Len / 2;
	  for(j = 0; j < 3; j++)
	    center[j] = Center[j];

	  while(1)		/* "th" will always point to an internal node */
	    {
	      len *= 0.5;
	      lenhalf *= 0.5;

	      subnode = 0;
	      if(PP[i].Pos[0] > center[0])
		{
		  center[0] += lenhalf;
		  subnode += 1;
		}
	      else
		{
		  center[0] -= lenhalf;
		}
	      if(PP[i].Pos[1] > center[1])
		{
		  center[1] += lenhalf;
		  subnode += 2;
		}
	      else
		{
		  center[1] -= lenhalf;
		}
	      if(PP[i].Pos[2] > center[2])
		{
		  center[2] += lenhalf;
		  subnode += 4;
		}
	      else
		{
		  center[2] -= lenhalf;
		}

	      nn = Nodes[th].u.suns[subnode];

	      if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
		{
		  if(nn >= NumPart)	/* the daugther node is an internal node. */
		    {
		      th = nn;
		    }
		  else		/* the daughter node is a particle. need to generate new node */
		    {
		      /* "nn"  is the particle index */

		      Nodes[th].u.suns[subnode] = nfree;
		      nfreep->u.suns[0] = -1;
		      nfreep->u.suns[1] = -1;
		      nfreep->u.suns[2] = -1;
		      nfreep->u.suns[3] = -1;
		      nfreep->u.suns[4] = -1;
		      nfreep->u.suns[5] = -1;
		      nfreep->u.suns[6] = -1;
		      nfreep->u.suns[7] = -1;

		      subnode = 0;
		      if(PP[nn].Pos[0] > center[0])
			subnode += 1;
		      if(PP[nn].Pos[1] > center[1])
			subnode += 2;
		      if(PP[nn].Pos[2] > center[2])
			subnode += 4;
#ifdef PERIODIC
		      if(len < 1.0e-6 * BoxSize)
			{
			  /* seems like we're dealing with particles   
			   * at identical locations. randomize 
			   * subnode index (well below gravitational softening scale). 
			   */
			  subnode = (int) (8.0 * drand48());
			  if(subnode >= 8)
			    subnode = 7;

			  printf("len=%g Len=%g i=%d BoxSize=%g subnode=%d (%g|%g|%g)\n",
				 len, Len, i, BoxSize, subnode, PP[i].Pos[0], PP[i].Pos[1], PP[i].Pos[2]);

			}
#endif
		      nfreep->u.suns[subnode] = nn;

		      th = nfree;	/* resume trying to insert the new particle at 
					   the newly created internal node */

		      numnodes++;
		      nfree++;
		      nfreep++;

		      if((numnodes + startnode) >= MaxNodes)
			{
			  printf("maximum number %d of tree-nodes reached.\n", MaxNodes);
			  printf("for particle %d\n", i);
			  exit(1);
			}
		    }
		}
	      else
		{
		  /* here we have found an empty slot where we can 
		   * attach the new particle as a leaf 
		   */
		  Nodes[th].u.suns[subnode] = i;
		  break;	/* done for this particle */
		}
	    }
	}
    }
  /* now walk the tree ones to set-up centers and faster walk */

  last = -1;

  force_update_node_recursive(NumPart + startnode, -1, Len, Center[0], Center[1], Center[2]);

  if(last >= 0)
    {
      if(last < NumPart)
	PP[last].Nextnode = -1;
      else
	Nodes[last].u.d.nextnode = -1;
    }

  printf("used %d nodes out of allocated %d. (filled fraction %g)\n",
	 numnodes, MaxNodes, (double) numnodes / MaxNodes);


  return numnodes;
}


void force_update_node_recursive(int no, int sib, double len, double cx, double cy, double cz)
{
  int j, jj, p, pp = 0, nextsib, suns[8];
  double ccx, ccy, ccz;


  if(no >= NumPart)		/* internal node */
    {
      for(j = 0; j < 8; j++)
	suns[j] = Nodes[no].u.suns[j];	/* this "backup" is necessary because the nextnode entry will
					   overwrite one element (union!) */
      Nodes[no].u.d.center[0] = cx;
      Nodes[no].u.d.center[1] = cy;
      Nodes[no].u.d.center[2] = cz;
      Nodes[no].u.d.len = len;
      Nodes[no].u.d.sibling = sib;

      if(last >= 0)
	{
	  if(last >= NumPart)
	    Nodes[last].u.d.nextnode = no;
	  else
	    PP[last].Nextnode = no;
	}
      last = no;

      for(j = 0; j < 8; j++)
	{
	  if((p = suns[j]) >= 0)
	    {
	      ccx = cx;
	      ccy = cy;
	      ccz = cz;

	      if((j & 1))
		ccx += len / 4;
	      else
		ccx -= len / 4;
	      if((j & 2))
		ccy += len / 4;
	      else
		ccy -= len / 4;
	      if((j & 4))
		ccz += len / 4;
	      else
		ccz -= len / 4;


	      /* check if we have a sibling on the same level */
	      for(jj = j + 1; jj < 8; jj++)
		if((pp = suns[jj]) >= 0)
		  break;

	      if(jj < 8)	/* yes, we do */
		nextsib = pp;
	      else
		nextsib = sib;

	      force_update_node_recursive(p, nextsib, len / 2, ccx, ccy, ccz);
	    }
	}
    }
  else				/* single particle */
    {
      if(last >= 0)
	{
	  if(last >= NumPart)
	    Nodes[last].u.d.nextnode = no;
	  else
	    PP[last].Nextnode = no;
	}
      last = no;
    }
}











float ngb_treefind(float xyz[3], int desngb, int *pnearest, float hguess)
{
  int numngb;
  float h2max;

  if(hguess == 0)
    hguess = pow(3 * desngb / (4 * M_PI) / part_dens, 0.33);

  do
    {
      numngb = ngb_treefind_variable(xyz, hguess);

      if(numngb == MAX_NGB)
	{
	  hguess /= 1.1;
	  continue;
	}

      if(numngb < desngb)
	{
	  hguess *= 1.26;
	  continue;
	}

      if(numngb >= desngb)
	{
	  h2max = ngb_select_closest(desngb, numngb, R2list - 1);
	  break;
	}

      hguess *= 1.26;
    }
  while(1);

  *pnearest = pmin;
  return h2max;
}



float ngb_select_closest(int k, int n, float *arr)
{
#define SWAP(a,b)  temp =(a);(a)=(b);(b)=temp;
  int i, ir, j, l, mid;
  float a, temp;

  l = 1;
  ir = n;
  while(1)
    {
      if(ir <= l + 1)
	{
	  if(ir == l + 1 && arr[ir] < arr[l])
	    {
	      SWAP(arr[l], arr[ir]);
	    }
	  return arr[k];
	}
      else
	{
	  mid = (l + ir) >> 1;
	  SWAP(arr[mid], arr[l + 1]);

	  if(arr[l] > arr[ir])
	    {
	      SWAP(arr[l], arr[ir]);
	    }
	  if(arr[l + 1] > arr[ir])
	    {
	      SWAP(arr[l + 1], arr[ir]);
	    }
	  if(arr[l] > arr[l + 1])
	    {
	      SWAP(arr[l], arr[l + 1]);
	    }
	  i = l + 1;
	  j = ir;
	  a = arr[l + 1];

	  while(1)
	    {
	      do
		i++;
	      while(arr[i] < a);

	      do
		j--;
	      while(arr[j] > a);

	      if(j < i)
		break;
	      SWAP(arr[i], arr[j]);
	    }

	  arr[l + 1] = arr[j];
	  arr[j] = a;

	  if(j >= k)
	    ir = j - 1;
	  if(j <= k)
	    l = i;
	}
    }
#undef SWAP
}


/* these macros map a coordinate difference
 * to the nearest periodic image
 */

#define NGB_PERIODIC(x) (((x)>BoxHalf)?((x)-BoxSize):(((x)<-BoxHalf)?((x)+BoxSize):(x)))






/*  ngb_treefind_variable() returns all neighbours (and only those) with distance <= hguess 
 *  and returns them in ngblistback and r2listback
 */
int ngb_treefind_variable(float searchcenter[3], double hguess)
{
  int k, numngb;
  int no, p, no_save;
  double dx, dy, dz, r2, h2;
  struct NODE *this;
  double searchmin[3], searchmax[3];

  h2 = hguess * hguess;

  for(k = 0; k < 3; k++)	/* cube-box window */
    {
      searchmin[k] = searchcenter[k] - hguess;
      searchmax[k] = searchcenter[k] + hguess;
    }

  numngb = 0;
  no = NumPart;
  r2min = 1.0e30;

  while(no >= 0)
    {
      if(no < NumPart)		/* single particle */
	{
	  p = no;
	  no = PP[no].Nextnode;

#ifdef PERIODIC
	  if(NGB_PERIODIC(PP[p].Pos[0] - searchcenter[0]) < (searchmin[0] - searchcenter[0]))
	    continue;
	  if(NGB_PERIODIC(PP[p].Pos[0] - searchcenter[0]) > (searchmax[0] - searchcenter[0]))
	    continue;
	  if(NGB_PERIODIC(PP[p].Pos[1] - searchcenter[1]) < (searchmin[1] - searchcenter[1]))
	    continue;
	  if(NGB_PERIODIC(PP[p].Pos[1] - searchcenter[1]) > (searchmax[1] - searchcenter[1]))
	    continue;
	  if(NGB_PERIODIC(PP[p].Pos[2] - searchcenter[2]) < (searchmin[2] - searchcenter[2]))
	    continue;
	  if(NGB_PERIODIC(PP[p].Pos[2] - searchcenter[2]) > (searchmax[2] - searchcenter[2]))
	    continue;
#else
	  if(PP[p].Pos[0] < (searchmin[0]))
	    continue;
	  if(PP[p].Pos[0] > (searchmax[0]))
	    continue;
	  if(PP[p].Pos[1] < (searchmin[1]))
	    continue;
	  if(PP[p].Pos[1] > (searchmax[1]))
	    continue;
	  if(PP[p].Pos[2] < (searchmin[2]))
	    continue;
	  if(PP[p].Pos[2] > (searchmax[2]))
	    continue;
#endif

#ifdef PERIODIC
	  dx = NGB_PERIODIC(PP[p].Pos[0] - searchcenter[0]);
	  dy = NGB_PERIODIC(PP[p].Pos[1] - searchcenter[1]);
	  dz = NGB_PERIODIC(PP[p].Pos[2] - searchcenter[2]);
#else
	  dx = PP[p].Pos[0] - searchcenter[0];
	  dy = PP[p].Pos[1] - searchcenter[1];
	  dz = PP[p].Pos[2] - searchcenter[2];
#endif
	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 < h2)
	    {
	      if(r2 < r2min)
		{
		  r2min = r2;
		  pmin = p;
		}
	      if(numngb < MAX_NGB)
		{
		  Ngblist[numngb] = p;
		  R2list[numngb++] = r2;
		}
	      else
		return numngb;
	    }
	}
      else
	{
	  this = &Nodes[no];

	  no_save = no;
	  no = Nodes[no].u.d.sibling;	/* in case the node can be discarded */
#ifdef PERIODIC
	  if((NGB_PERIODIC(this->u.d.center[0] - searchcenter[0]) + 0.5 * this->u.d.len) <
	     (searchmin[0] - searchcenter[0]))
	    continue;
	  if((NGB_PERIODIC(this->u.d.center[0] - searchcenter[0]) - 0.5 * this->u.d.len) >
	     (searchmax[0] - searchcenter[0]))
	    continue;
	  if((NGB_PERIODIC(this->u.d.center[1] - searchcenter[1]) + 0.5 * this->u.d.len) <
	     (searchmin[1] - searchcenter[1]))
	    continue;
	  if((NGB_PERIODIC(this->u.d.center[1] - searchcenter[1]) - 0.5 * this->u.d.len) >
	     (searchmax[1] - searchcenter[1]))
	    continue;
	  if((NGB_PERIODIC(this->u.d.center[2] - searchcenter[2]) + 0.5 * this->u.d.len) <
	     (searchmin[2] - searchcenter[2]))
	    continue;
	  if((NGB_PERIODIC(this->u.d.center[2] - searchcenter[2]) - 0.5 * this->u.d.len) >
	     (searchmax[2] - searchcenter[2]))
	    continue;
#else
	  if((this->u.d.center[0] + 0.5 * this->u.d.len) < (searchmin[0]))
	    continue;
	  if((this->u.d.center[0] - 0.5 * this->u.d.len) > (searchmax[0]))
	    continue;
	  if((this->u.d.center[1] + 0.5 * this->u.d.len) < (searchmin[1]))
	    continue;
	  if((this->u.d.center[1] - 0.5 * this->u.d.len) > (searchmax[1]))
	    continue;
	  if((this->u.d.center[2] + 0.5 * this->u.d.len) < (searchmin[2]))
	    continue;
	  if((this->u.d.center[2] - 0.5 * this->u.d.len) > (searchmax[2]))
	    continue;
#endif
	  no = this->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  return numngb;
}





/* Allocate memory for the neighbour lists
 */
void ngb_treeallocate(int maxnodes)
{
  int totbytes = 0, bytes;

  MaxNodes = maxnodes;

  if(!(Nodes_base = malloc(bytes = (MaxNodes + 1) * sizeof(struct NODE))))
    {
      printf("failed to allocate memory for %d tree-nodes (%d bytes).\n", MaxNodes, bytes);
      exit(3);
    }
  totbytes += bytes;


  if(!(R2list = malloc(MAX_NGB * sizeof(float))))
    {
      printf("failed to allocate memory for R2list\n");
      exit(3);
    }
  totbytes += MAX_NGB * sizeof(float);

  if(!(Ngblist = malloc(MAX_NGB * sizeof(int))))
    {
      printf("failed to allocate memory for Ngblist\n");
      exit(3);
    }
  totbytes += MAX_NGB * sizeof(int);

  printf("allocated %f Mbyte for ngb search.\n", ((double) totbytes) / (1024.0 * 1024.0));
}


/* To construct the neighbour tree,
 * we actually need to construct the force-tree,
 * because we use it now for the neighbour search.
 * This routine is obsolute at this point.
 */
void ngb_treebuild(void)
{
  printf("Begin Ngb-tree construction.\n");

  force_treebuild();

  printf("Ngb-Tree contruction finished \n");
}
