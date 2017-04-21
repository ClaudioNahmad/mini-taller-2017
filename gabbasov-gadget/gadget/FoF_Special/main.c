#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>


#include "allvars.h"
#include "nrsrc/nrutil.h"


double LinkLength = 0.2;	/* in terms of mean interparticle
				   seperation */

int GroupMinLen = 32;		/* store only groups in the catalogue
				   with a mass at least as high as the
				   equivalent of this number of
				   dark matter particles */


/* units used for Gadget run */

double UnitLength_in_cm  = 3.085678e21;     /*  1.0 kpc */
double UnitMass_in_g = 1.989e43;            /*  1.0e10 solar masses */
double UnitVelocity_in_cm_per_s = 1e5;      /*  1 km/sec */

/* -----------------------------------------------------------------------*/




int main(int argc, char **argv)
{
  int gr;
  int i, n, p;
  int tot;
  char path[500];
  char input_fname[400];
  char indexlist_fname[400];
  char catalogue_fname[400];
  char particles_fname[400];
  char basename[400];
  char buf[400];
  double rhomean;
  double UnitTime_in_s, G, Hubble;

  if(argc != 4)
    {
      fprintf(stderr, "\n\nwrong argument(s).  Specify:\n\n");
      fprintf(stderr, "<path>      (path)\n");
      fprintf(stderr, "<basename>  (basename of snapshot files)\n");
      fprintf(stderr, "<num>       (number of snapshot)\n");
      fprintf(stderr, "\n\n");
      exit(0);
    }

  strcpy(path, argv[1]);
  strcpy(basename, argv[2]);
  Snapshot = atoi(argv[3]);


  sprintf(input_fname, "%s/%s_%03d", path, basename, Snapshot);
  sprintf(indexlist_fname, "%s/groups_indexlist/fof_special_indexlist_%03d", path, Snapshot);
  sprintf(catalogue_fname, "%s/groups_catalogue/fof_special_catalogue_%03d", path, Snapshot);
  sprintf(particles_fname, "%s/groups_particles/fof_special_particles_%03d", path, Snapshot);

  sprintf(buf, "%s/groups_indexlist", path);
  mkdir(buf, 02755);
  sprintf(buf, "%s/groups_catalogue", path);
  mkdir(buf, 02755);
  sprintf(buf, "%s/groups_particles", path);
  mkdir(buf, 02755);


  Files = find_files(input_fname);

  loadpositions(input_fname, Files);
  printf("BoxSize= %g\n", BoxSize);
  printf("NumPart= %d\n", NumPart);


  peano_hilbert_order();

  /* set-up gadget units  */
  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
  Hubble = HUBBLE * UnitTime_in_s;

  rhomean = Omega0 * 3 * Hubble * Hubble / (8 * M_PI * G);

  LinkL = LinkLength * pow(M_DM / rhomean, 1.0 / 3);
  SearchRadius = LinkL;

  printf("Comoving linking length: %g kpc/h \n\n", LinkL);


  GroupMinMass = GroupMinLen * M_DM;



  Len = ivector(1, NumPart);
  Head = ivector(1, NumPart);
  Next = ivector(1, NumPart);
  Tail = ivector(1, NumPart);


  for(i = 1; i <= NumPart; i++)
    {
      Head[i] = Tail[i] = i;
      Next[i] = 0;
      Len[i] = 1;
    }



  ngb_treeallocate(0.9 * N_DM);
  ngb_treebuild();


  find_groups();


  Nearest = ivector(1, NumPart);

  determine_nearest();

  link_nearest();

  free_ivector(Nearest, 1, NumPart);


  for(n = 1, Ngroups = 0, tot = 0; n <= NumPart; n++)
    {
      if(Head[n] == n)
	{
	  Ngroups++;
	}
      else
	{
	  Len[n] = 0;
	}

      tot += Len[n];
    }

  printf("Groups=%d (all sizes)  NumPart=%d   tot=%d\n", Ngroups, NumPart, tot);



  free_ivector(Tail, 1, NumPart);

  GroupLen = ivector(1, Ngroups);
  GroupTag = ivector(1, Ngroups);
  GroupMass = vector(1, Ngroups);

  for(n = 1, gr = 1; n <= NumPart; n++)
    {
      if(Len[n])
	{
	  GroupLen[gr] = Len[n];
	  GroupTag[gr] = Head[n];

	  GroupMass[gr] = 0;
	  p = Head[n];
	  do
	    {
	      GroupMass[gr] += P[p].Mass;
	    }
	  while((p = Next[p]));

	  gr++;
	}
    }

  sort2_fltint(Ngroups, GroupMass, GroupLen, GroupTag);	/* ok, sorted by mass */


  /* Head, Tail, Len can also be freed */

  free_ivector(Head, 1, NumPart);
  free_ivector(Len, 1, NumPart);


  NextFinal = ivector(1, NumPart);
  IdSort = ivector(1, NumPart);

  for(gr = 1; gr <= Ngroups; gr++)
    {
      n = GroupLen[gr];
      p = GroupTag[gr];

      i = 1;

      do
	{
	  IdSort[i++] = p;
	}
      while((p = Next[p]));

      sort_int(n, IdSort);

      GroupTag[gr] = IdSort[1];

      for(i = 1; i < n; i++)
	NextFinal[IdSort[i]] = IdSort[i + 1];

      NextFinal[IdSort[n]] = IdSort[1];
    }

  free_ivector(IdSort, 1, NumPart);


  NgroupsAboveMinLen = 0;

  for(gr = 1; gr <= Ngroups; gr++)
    {
      if(GroupMass[gr] >= GroupMinMass)
	NgroupsAboveMinLen++;
    }

  printf("\nNumber of groups with a mass equivalent to %d dm particles: %d\n", GroupMinLen,
	 NgroupsAboveMinLen);
  printf("Largest group has %d particles.\n\n", GroupLen[Ngroups]);


  save_groups(indexlist_fname, catalogue_fname, particles_fname);

  free_ivector(Next, 1, NumPart);
  free_ivector(NextFinal, 1, NumPart);
  free_ivector(GroupTag, 1, Ngroups);
  free_ivector(GroupLen, 1, Ngroups);
  P++;
  free(P);

  return 0;
}






void save_groups(char *indexlist_fname, char *catalogue_fname, char *particles_fname)
{
  FILE *fd;
  int i, k, p, gr, offset;
  double mass, s[3], d[3];
  float ss[3];
  int ngas, ndm, nstars, ntot;
  float mgas, mdm, mstars, sfr, metgas, metstars, mcold;

#ifdef VELOCITYDISPERSION
  float sigma, vx, vy, vz;
  int count;
#endif

  printf("writing group catalogue...\n");

  if(!(fd = fopen(catalogue_fname, "w")))
    {
      printf("can't open file `%s`\n", catalogue_fname);
      exit(0);
    }

  fwrite(&NgroupsAboveMinLen, sizeof(int), 1, fd);

  for(gr = Ngroups; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      fwrite(&GroupLen[gr], sizeof(int), 1, fd);

  for(gr = Ngroups, offset = 0; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      {
	fwrite(&offset, sizeof(int), 1, fd);
	offset += GroupLen[gr];
      }

  for(gr = Ngroups; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      {
	fwrite(&GroupMass[gr], sizeof(float), 1, fd);
      }

  for(gr = Ngroups; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      {
	mass = s[0] = s[1] = s[2] = 0;
	for(i = 0, p = GroupTag[gr]; i < GroupLen[gr]; i++, p = NextFinal[p])
	  {
	    for(k = 0; k < 3; k++)
	      {
		d[k] = periodic(P[p].Pos[k] - P[GroupTag[gr]].Pos[k]);
		s[k] += P[p].Mass * d[k];
	      }
	    mass += P[p].Mass;
	  }

	for(k = 0; k < 3; k++)
	  {
	    s[k] = s[k] / mass + P[GroupTag[gr]].Pos[k];
	    ss[k] = periodic_wrap(s[k]);
	  }

	fwrite(ss, sizeof(float), 3, fd);
      }

  for(gr = Ngroups; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      {
	ngas = ndm = nstars = 0;
	for(i = 0, p = GroupTag[gr]; i < GroupLen[gr]; i++, p = NextFinal[p])
	  {
	    if(P[p].Type == 0)
	      ngas++;
	    if(P[p].Type == 1)
	      ndm++;
	    if(P[p].Type == 4)
	      nstars++;
	  }
	fwrite(&ngas, sizeof(int), 1, fd);
	fwrite(&ndm, sizeof(int), 1, fd);
	fwrite(&nstars, sizeof(int), 1, fd);
      }

  for(gr = Ngroups; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      {
	mgas = mdm = mstars = 0;
	for(i = 0, p = GroupTag[gr]; i < GroupLen[gr]; i++, p = NextFinal[p])
	  {
	    if(P[p].Type == 0)
	      mgas += P[p].Mass;
	    if(P[p].Type == 1)
	      mdm += P[p].Mass;
	    if(P[p].Type == 4)
	      mstars += P[p].Mass;
	  }
	fwrite(&mgas, sizeof(float), 1, fd);
	fwrite(&mdm, sizeof(float), 1, fd);
	fwrite(&mstars, sizeof(float), 1, fd);
      }

  for(gr = Ngroups; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      {
	sfr = 0;
	for(i = 0, p = GroupTag[gr]; i < GroupLen[gr]; i++, p = NextFinal[p])
	  {
	    if(P[p].Type == 0)
	      sfr += P[p].Sfr;
	  }
	fwrite(&sfr, sizeof(float), 1, fd);
      }

  for(gr = Ngroups; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      {
	metgas = metstars = 0;
	mgas = mstars = 0;
	for(i = 0, p = GroupTag[gr]; i < GroupLen[gr]; i++, p = NextFinal[p])
	  {
	    if(P[p].Type == 0)
	      {
		mgas += P[p].Mass;
		metgas += P[p].Mass * P[p].Metallicity;
	      }

	    if(P[p].Type == 4)
	      {
		mstars += P[p].Mass;
		metstars += P[p].Mass * P[p].Metallicity;
	      }
	  }
	if(mgas > 0)
	  metgas /= mgas;
	if(mstars > 0)
	  metstars /= mstars;


	fwrite(&metgas, sizeof(float), 1, fd);
	fwrite(&metstars, sizeof(float), 1, fd);
      }

  for(gr = Ngroups; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      {
	mcold = 0;
	for(i = 0, p = GroupTag[gr]; i < GroupLen[gr]; i++, p = NextFinal[p])
	  {
	    if(P[p].Type == 0)
	      {
		if(P[p].Sfr > 0)
		  mcold += P[p].Mass;
	      }
	  }

	fwrite(&mcold, sizeof(float), 1, fd);
      }
#ifdef VELOCITYDISPERSION
  for(gr = Ngroups; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      {
	vx = vy = vz = 0;
	count = 0;
	for(i = 0, p = GroupTag[gr]; i < GroupLen[gr]; i++, p = NextFinal[p])
	  {
	    if(P[p].Type == 4)
	      {
		vx += P[p].Vel[0];
		vy += P[p].Vel[1];
		vz += P[p].Vel[2];
		count++;
	      }
	  }
	if(count == 0)
	  sigma = 0;
	else
	  {
	    sigma = 0;
	    vx /= count;
	    vy /= count;
	    vz /= count;
	    for(i = 0, p = GroupTag[gr]; i < GroupLen[gr]; i++, p = NextFinal[p])
	      {
		if(P[p].Type == 4)
		  {
		    sigma += (P[p].Vel[0] - vx) * (P[p].Vel[0] - vx);
		    sigma += (P[p].Vel[1] - vy) * (P[p].Vel[1] - vy);
		    sigma += (P[p].Vel[2] - vz) * (P[p].Vel[2] - vz);
		  }
	      }
	    sigma = sqrt(sigma / (3 * count));
	  }

	sigma *= sqrt(Time);

	fwrite(&sigma, sizeof(float), 1, fd);
      }


  for(gr = Ngroups; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      {
	vx = vy = vz = 0;
	count = 0;
	for(i = 0, p = GroupTag[gr]; i < GroupLen[gr]; i++, p = NextFinal[p])
	  {
	    if(P[p].Type == 1)
	      {
		vx += P[p].Vel[0];
		vy += P[p].Vel[1];
		vz += P[p].Vel[2];
		count++;
	      }
	  }
	if(count == 0)
	  sigma = 0;
	else
	  {
	    sigma = 0;
	    vx /= count;
	    vy /= count;
	    vz /= count;
	    for(i = 0, p = GroupTag[gr]; i < GroupLen[gr]; i++, p = NextFinal[p])
	      {
		if(P[p].Type == 1)
		  {
		    sigma += (P[p].Vel[0] - vx) * (P[p].Vel[0] - vx);
		    sigma += (P[p].Vel[1] - vy) * (P[p].Vel[1] - vy);
		    sigma += (P[p].Vel[2] - vz) * (P[p].Vel[2] - vz);
		  }
	      }
	    sigma = sqrt(sigma / (3 * count));
	  }

	sigma *= sqrt(Time);

	fwrite(&sigma, sizeof(float), 1, fd);
      }
#endif

  fclose(fd);




  printf("writing group indexlist...\n");


  if(!(fd = fopen(indexlist_fname, "w")))
    {
      printf("can't open file `%s`\n", indexlist_fname);
      exit(0);
    }

  for(gr = Ngroups, ntot = 0; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      ntot += GroupLen[gr];

  fwrite(&ntot, sizeof(int), 1, fd);

  for(gr = Ngroups; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      {
	for(i = 0, p = GroupTag[gr]; i < GroupLen[gr]; i++, p = NextFinal[p])
	  {
	    fwrite(&P[p].FileIndex, sizeof(int), 1, fd);
	  }
      }
  fclose(fd);



  printf("writing group particles...\n");

  if(!(fd = fopen(particles_fname, "w")))
    {
      printf("can't open file `%s`\n", particles_fname);
      exit(0);
    }

  for(gr = Ngroups, ntot = 0; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      ntot += GroupLen[gr];

  fwrite(&ntot, sizeof(int), 1, fd);

  for(gr = Ngroups; gr >= 1; gr--)
    if(GroupMass[gr] >= GroupMinMass)
      {
	for(i = 0, p = GroupTag[gr]; i < GroupLen[gr]; i++, p = NextFinal[p])
	  {
	    fwrite(&P[p].Pos[0], sizeof(float), 3, fd);
	  }
      }
  fclose(fd);


  printf("done.\n");

}






void find_groups(void)
{
  int i, j, p, s, pp, ss, numngb, signal = 0;

  printf("linking...");
  fflush(stdout);

  for(i = 1; i <= NumPart; i++)
    {
      if(P[i].Type == 1)
	{
	  if(i > signal + NumPart / 100.0)
	    {
	      printf("x");
	      fflush(stdout);
	      signal += (NumPart / 100.0 + 1);
	    }

	  numngb = ngb_treefind_variable(P[i].Pos, LinkL);

	  for(j = 0; j < numngb; j++)
	    {
	      s = i;
	      p = Ngblist[j] + 1;

	      if(Head[p] != Head[s])	/* only if not yet linked */
		{
		  if(Len[Head[p]] > Len[Head[s]])	/* p group is longer */
		    {
		      Next[Tail[Head[p]]] = Head[s];

		      Tail[Head[p]] = Tail[Head[s]];

		      Len[Head[p]] += Len[Head[s]];

		      ss = Head[s];
		      do
			{
			  Head[ss] = Head[p];
			}
		      while((ss = Next[ss]));
		    }
		  else
		    {
		      Next[Tail[Head[s]]] = Head[p];

		      Tail[Head[s]] = Tail[Head[p]];

		      Len[Head[s]] += Len[Head[p]];

		      pp = Head[p];
		      do
			{
			  Head[pp] = Head[s];
			}
		      while((pp = Next[pp]));
		    }
		}
	    }
	}
    }

  printf("\ndone.\n");
}







void determine_nearest(void)
{
  int i, signal, pnearest;
  double h = 0;


  printf("finding nearest dark matter neighbour\n");

  for(i = 1, signal = 0; i <= NumPart; i++)
    {
      Nearest[i] = 0;

      if(P[i].Type != 1)
	{
	  if(i > (signal / 100.0) * NumPart)
	    {
	      printf("x");
	      fflush(stdout);
	      signal++;
	    }

	  h = sqrt(ngb_treefind(P[i].Pos, 10, &pnearest, h * 1.1));

	  Nearest[i] = pnearest + 1;
	}
    }

  printf("\ndone with finding of nearest DM particle.\n");
}


void link_nearest(void)
{
  int i, p, s, ss, pp;


  printf("linking to nearest DM particle.\n");

  for(i = 1; i <= NumPart; i++)
    {
      if(P[i].Type != 1 && Nearest[i] > 0)
	{
	  p = i;
	  s = Nearest[i];

	  if(Head[p] != Head[s])	/* only if not yet linked */
	    {
	      if(Len[Head[p]] > Len[Head[s]])	/* p group is longer */
		{
		  Next[Tail[Head[p]]] = Head[s];

		  Tail[Head[p]] = Tail[Head[s]];

		  Len[Head[p]] += Len[Head[s]];

		  ss = Head[s];
		  do
		    {
		      Head[ss] = Head[p];
		    }
		  while((ss = Next[ss]));
		}
	      else
		{
		  Next[Tail[Head[s]]] = Head[p];

		  Tail[Head[s]] = Tail[Head[p]];

		  Len[Head[s]] += Len[Head[p]];

		  pp = Head[p];
		  do
		    {
		      Head[pp] = Head[s];
		    }
		  while((pp = Next[pp]));
		}
	    }
	}
    }

  printf("done.\n");
}






double periodic(double x)
{
#ifdef PERIODIC
  if(x > 0.5 * BoxSize)
    x -= BoxSize;

  if(x < -0.5 * BoxSize)
    x += BoxSize;
#endif
  return x;
}

double periodic_wrap(double x)
{
#ifdef PERIODIC
  while(x > BoxSize)
    x -= BoxSize;

  while(x < 0)
    x += BoxSize;
#endif
  return x;
}
