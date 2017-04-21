#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



#include "allvars.h"






void allocate_memory(void)
{
  fprintf(stderr, "allocating memory...\n");

  if(NumPart > 0)
    {
      if(!(P = malloc(NumPart * sizeof(struct particle_data))))
	{
	  fprintf(stderr, "failed to allocate memory. (A)\n");
	  exit(0);
	}

      PP = P;

      P--;			/* start with offset 1 */

      /*
         if(!(Id = malloc(NumPart * sizeof(int))))
         {
         fprintf(stderr, "failed to allocate memory. (B)\n");
         exit(0);
         }

         Id--;  */
      /* start with offset 1 */
    }

  fprintf(stderr, "allocating memory...done\n");
}
