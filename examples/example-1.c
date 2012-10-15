/**
    libgapmis: a library for pairwise sequence aligment with a single gap.
    Copyright (C) 2012 Nikolaos Alachiotis, Simon Berger, Tomas Flouri, and
    Solon P. Pissis. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <errno.h>
#include <sys/time.h>
#include <libgapmis/gapmis.h>
#include <libgapmis/errors.h>

#define MAIN_MAX 	             500

int main ( int argc, char * argv [] )
 {
   struct gapmis_params         in;
   struct gapmis_align        * out;
   
   FILE                       * fd_query;   /* File descriptors for the reads*/
   FILE                       * fd_target;  /* File descriptors for the reads*/
   FILE                       * fd_errors;  /* File descriptors for the errors*/
   
   char                         query[MAIN_MAX];     /* Buffer for storing the query */
   char                         queryId[MAIN_MAX];   /* Buffer for storing the Id of the query */
   char                         target[MAIN_MAX];     /* Buffer for storing the target */
   char                         targetId[MAIN_MAX];   /* Buffer for storing the Id of the target */
   char                         errors[MAIN_MAX];    /* Buffer for storing the errors */

   unsigned int			valid = 0, correct = 0, total = 0, invalid = 0, incorrect = 0;
   unsigned int i;
   unsigned int			gaps_len = 0, num_gaps = 0, num_mis = 0;
   unsigned int			real_gaps_len = 0, real_num_gaps = 0, real_num_mis = 0;

   /* assign the input parameters */
   in . scoring_matrix = 0;
   in . gap_open_pen   = atof( argv[5] );
   in . gap_extend_pen = atof( argv[6] );
   in . max_gap        = atoi( argv[7] ); 
   in . num_frags      = atoi( argv[8] );


   /* open file 1 (query sequence) */
   if ( ! ( fd_query = fopen ( argv[1], "r") ) ) 
    {
      fprintf ( stderr, "Cannot open file %s\n", argv[1] );
      return ( 0 ); 
    }

   /* open file 2 (target sequence) */
   if ( ! ( fd_target = fopen ( argv[2], "r") ) ) 
    {
      fprintf ( stderr, "Cannot open file %s\n", argv[2] );
      return ( 0 );
    }

   /* open file 3 ( errors ) */
   if ( ! ( fd_errors = fopen ( argv[3], "r") ) ) 
    {
      fprintf ( stderr, "Cannot open file %s\n", argv[3] );
      return ( 0 );
    }

   while ( ! feof ( fd_query ) && ! feof ( fd_target ) && ! feof ( fd_errors ) )
    {

      /* allocate the space for the output */
      out = ( struct gapmis_align * ) calloc ( in . num_frags , sizeof ( struct gapmis_align ) );

      /* read file 1 (query sequence) */
      if ( fgetc ( fd_query ) && fgets ( queryId, MAIN_MAX, fd_query ) && fgets ( query, MAIN_MAX, fd_query ) )
       {
          query[ strlen ( query ) - 1] = 0;
          queryId[ strlen ( queryId ) - 1] = 0;
       }
      else
		break; 

      /* read file 2 (target sequence) */
      if ( fgetc ( fd_target ) && fgets ( targetId, MAIN_MAX, fd_target ) && fgets ( target, MAIN_MAX, fd_target ) )
       {
         target[ strlen ( target ) - 1] = 0;
         targetId[ strlen ( targetId ) - 1] = 0;
       }
      else
		break; 

      /* read file 3 (errors) */
      if ( fgets ( errors, MAIN_MAX, fd_errors ) )
       {
		char * pch;
		pch = strtok ( errors, " ,.-" );
    		real_num_gaps = atoi( pch );
    		pch = strtok ( NULL, " ,.-");
    		real_gaps_len = atoi( pch );
    		pch = strtok ( NULL, " ,.-");
    		real_num_mis = atoi( pch );
       }
      else
		break; 

   	gaps_len = 0, num_gaps = 0, num_mis = 0;

   	gapsmis_one_to_one_onf ( query, target, &in, out ) ;

        for ( i = 0; i < in . num_frags; i ++ )
	{
		if ( out[i] . min_gap > 0 ) 
		{
			num_gaps += 1;
			gaps_len += out[i] . min_gap;
		}
		num_mis += out [i] . num_mis;	
	}	
		
	if ( num_gaps <= real_num_gaps )
	{

		valid ++;
		if ( gaps_len <= real_gaps_len && num_mis <= real_num_mis )
			correct ++;
		else
		{
			incorrect++;
		#if 0
                fprintf ( stdout, "%s\n", query );
                fprintf ( stdout, "%s\n", target );
                fprintf ( stdout, "\nnum_gaps: %d gaps_len: %d num_mis: %d\n", num_gaps, gaps_len, num_mis );
                fprintf ( stdout, "\nreal_num_gaps: %d real_gaps_len: %d real_num_mis: %d\n", real_num_gaps, real_gaps_len, real_num_mis );
   		gapsmis_results_one_to_one ( argv[4], query, queryId, target, targetId, &in, out );
		exit ( 0 );
		#endif
		}
	}
	else
	{
		invalid ++;
	}
	total ++;
        free ( out );
    }

   fprintf ( stdout, "\nTotal: %d Valid: %d Correct: %d Invalid: %d Incorrect: %d\n", total, valid, correct, invalid, incorrect );

   fclose ( fd_target );
   fclose ( fd_query );
   fclose ( fd_errors );

   return ( 0 );
 }
