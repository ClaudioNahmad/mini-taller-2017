#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


void loadpositions(char *fname, int files)
{
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
  FILE *fd;
  float mas;
  char buf[200];
  int i, j, k, dummy, dummy2;
  int n, pc, pc_new, ntot_withmasses;
  float pos[3], sfr, met;


  for(i = 0, pc = 1; i < files; i++, pc = pc_new)
    {
      if(files > 1)
	sprintf(buf, "%s.%d", fname, i);
      else
	sprintf(buf, "%s", fname);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s`\n", buf);
	  exit(0);
	}

      printf("reading `%s' ...\n", buf);
      fflush(stdout);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header, sizeof(header), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);


      printf("header.time= %g\n", header.time);

      if(files == 1)
	{
	  NumPart = header.npart[0] + header.npart[1] + header.npart[4];
	}
      else
	{
	  NumPart = header.npartTotal[0] + header.npartTotal[1] + header.npartTotal[4];
	}

      for(k = 0, ntot_withmasses = 0; k < 6; k++)
	{
	  if(header.mass[k] == 0 && header.npartTotal[k] > 0)
	    ntot_withmasses += header.npartTotal[k];
	}

      if(i == 0)
	printf("ntot_withmasses= %d\n", ntot_withmasses);

      if(i == 0)
	{
	  allocate_memory();

	  for(j = 0; j < 6; j++)
	    printf("npart[%d]= %d\n", j, header.npartTotal[j]);
	}


      fread(&dummy, sizeof(dummy), 1, fd);
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header.npart[k]; n++)
	    {
	      fread(&pos[0], sizeof(float), 3, fd);

	      if(k == 0 || k == 1 || k == 4)
		{
		  P[pc_new].Pos[0] = pos[0];
		  P[pc_new].Pos[1] = pos[1];
		  P[pc_new].Pos[2] = pos[2];
		  P[pc_new].FileIndex = pc_new;
		  pc_new++;
		}
	    }
	}
      fread(&dummy, sizeof(dummy), 1, fd);

#ifndef VELOCITYDISPERSION
      fread(&dummy, sizeof(dummy), 1, fd);
      fseek(fd, dummy, SEEK_CUR);	/* skip velocities */
      fread(&dummy, sizeof(dummy), 1, fd);
#else
      fread(&dummy, sizeof(dummy), 1, fd);
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header.npart[k]; n++)
	    {
	      fread(&pos[0], sizeof(float), 3, fd);
	      if(k == 0 || k == 1 || k == 4)
		{
		  P[pc_new].Vel[0] = pos[0];
		  P[pc_new].Vel[1] = pos[1];
		  P[pc_new].Vel[2] = pos[2];
		  pc_new++;
		}
	    }
	}
      fread(&dummy, sizeof(dummy), 1, fd);
#endif
      fread(&dummy, sizeof(dummy), 1, fd);
      fseek(fd, dummy, SEEK_CUR);	/* skip ID */
      fread(&dummy, sizeof(dummy), 1, fd);


      if(ntot_withmasses)
	fread(&dummy, sizeof(dummy), 1, fd);
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header.npart[k]; n++)
	    {
	      if(header.mass[k] == 0)
		fread(&mas, sizeof(float), 1, fd);
	      else
		mas = header.mass[k];

	      if(k == 0 || k == 1 || k == 4)
		{
		  P[pc_new].Mass = mas;
		  pc_new++;
		}
	    }
	}
      if(ntot_withmasses)
	fread(&dummy, sizeof(dummy), 1, fd);


      fread(&dummy, sizeof(dummy), 1, fd);	/* skip Temp */
      fseek(fd, sizeof(float) * header.npart[0], SEEK_CUR);
      fread(&dummy, sizeof(dummy), 1, fd);

      fread(&dummy, sizeof(dummy), 1, fd);	/* skip rho */
      fseek(fd, sizeof(float) * header.npart[0], SEEK_CUR);
      fread(&dummy, sizeof(dummy), 1, fd);


      /* now electron abundance */
      if(header.flag_cooling)
	{
	  fread(&dummy, sizeof(dummy), 1, fd);	/* skip Ne */
	  fseek(fd, sizeof(float) * header.npart[0], SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fd);

	  fread(&dummy, sizeof(dummy), 1, fd);	/* skip NH0 */
	  fseek(fd, sizeof(float) * header.npart[0], SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fd);
	}

      fread(&dummy, sizeof(dummy), 1, fd);	/* skip hsml */
      fseek(fd, sizeof(float) * header.npart[0], SEEK_CUR);
      fread(&dummy, sizeof(dummy), 1, fd);


      if(header.flag_sfr)
	{
	  fread(&dummy, sizeof(dummy), 1, fd);
	  for(k = 0, pc_new = pc; k < 6; k++)
	    {
	      for(n = 0; n < header.npart[k]; n++)
		{
		  if(k == 0 || k == 1 || k == 4)
		    {
		      if(k == 0)
			{
			  fread(&sfr, sizeof(float), 1, fd);	/* SFR */
			  P[pc_new].Sfr = sfr;
			}
		      else
			P[pc_new].Sfr = 0;
		      pc_new++;
		    }
		}
	    }
	  fread(&dummy2, sizeof(dummy2), 1, fd);
	  if(dummy2 != dummy)
	    {
	      printf("dummy=%d  dummy2=%d\n", dummy, dummy2);
	      exit(30);
	    }
	}
      else
	{
	  for(k = 0, pc_new = pc; k < 6; k++)
	    {
	      for(n = 0; n < header.npart[k]; n++)
		{
		  if(k == 0 || k == 1 || k == 4)
		    {
		      P[pc_new].Sfr = 0;
		      pc_new++;
		    }
		}
	    }
	}


      if(header.flag_sfr)
	{
	  if(header.npart[4])
	    {
	      fread(&dummy, sizeof(dummy), 1, fd);	/* skip stellar ages */
	      fseek(fd, sizeof(float) * header.npart[4], SEEK_CUR);
	      fread(&dummy2, sizeof(dummy2), 1, fd);
	      if(dummy2 != dummy)
		{
		  printf("problem: dummy=%d  dummy2=%d  %d\n", dummy, dummy2, header.npart[4]);
		  exit(30);
		}
	    }

	  fread(&dummy, sizeof(dummy), 1, fd);
	  for(k = 0, pc_new = pc; k < 6; k++)
	    {
	      for(n = 0; n < header.npart[k]; n++)
		{
		  if(k == 0 || k == 1 || k == 4)
		    {
		      if(k == 0 || k == 4)
			{
			  fread(&met, sizeof(float), 1, fd);	/* metallicity */
			  P[pc_new].Metallicity = met;
			}
		      else
			P[pc_new].Metallicity = 0;
		      pc_new++;
		    }
		}
	    }
	  fread(&dummy2, sizeof(dummy2), 1, fd);

	  if(dummy2 != dummy)
	    {
	      printf("dummy=%d  dummy2!=%d\n", dummy, dummy2);
	      exit(30);
	    }
	}
      else
	{
	  for(k = 0, pc_new = pc; k < 6; k++)
	    {
	      for(n = 0; n < header.npart[k]; n++)
		{
		  if(k == 0 || k == 1 || k == 4)
		    {
		      P[pc_new].Metallicity = 0;
		      pc_new++;
		    }
		}
	    }
	}


      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header.npart[k]; n++)
	    {
	      if(k == 0 || k == 1 || k == 4)
		{
		  P[pc_new].Type = k;
		  pc_new++;
		}
	    }
	}


      fclose(fd);
    }


  Time = header.time;

  Omega0 = header.Omega0;
  BoxSize = header.BoxSize;
  BoxHalf = BoxSize / 2;

  if(files == 1)
    N_DM = header.npart[1];
  else
    N_DM = header.npartTotal[1];


  for(n = 1; n <= NumPart; n++)
    if(P[n].Type == 1)
      {
	M_DM = P[n].Mass;
	break;
      }
  if(n > NumPart)
    exit(0);


}


int find_files(char *fname)
{
  FILE *fd;
  char buf[200], buf1[200];
  int dummy;


  sprintf(buf, "%s.%d", fname, 0);
  sprintf(buf1, "%s", fname);

  if((fd = fopen(buf, "r")))
    {
      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header, sizeof(header), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      fclose(fd);

      return header.num_files;
    }

  if((fd = fopen(buf1, "r")))
    {
      fclose(fd);

      return 1;
    }

  printf("Error. Can't find snapshot!\nneither as `%s'\nnor as `%s'\n\n", buf, buf1);

  exit(1);
  return 0;
}
