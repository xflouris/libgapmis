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
#include <time.h>
#include "gapmis.h"
#include "errors.h"
#include "EDNAFULL.h"
#include "EBLOSUM62.h"

/* Computes the optimal semi-global alignment between a number of fragments of pattern p and text t */
unsigned int gapsmis_one_to_one ( const char * p, const char * t, const struct gapmis_params * in, struct gapmis_align * out )
{
   unsigned int		i;
   unsigned int		m;
   unsigned int		pb;
   unsigned int		mm;
   unsigned int		num_frags = in -> num_frags;
   
   m = strlen ( p );
 
   if ( num_frags == 0 || num_frags > m )
    {
      errno = FRAGS; //Error: the num of frags should be less or equal to the length of p!!!
      return ( 0 );
    }

   mm = m / num_frags;
   pb = 0;		//processed bytes

   for ( i = 0; i < num_frags; ++ i )
    {
         char		* frag;

	/* calculate the length of the pattern's fragment */
	if ( i == num_frags - 1 )  mm = m - pb;  

	/* allocate and copy the pattern's fragment */
        if ( ! ( frag = ( char * ) calloc ( mm + 1 , sizeof ( char ) ) ) )
         {
           errno = MALLOC; //Error: frag could not be allocated!!!
           return ( 0 );
         } 
        memcpy ( frag, p, mm );	
 	frag[mm] = '\0';

	/* perform the single-gap alignment */
        if ( ! ( gapmis_one_to_one ( frag, t, in, &out[i] ) ) )
           return ( 0 );
	
	/* update the pointers */
	p  += mm;
        switch ( out[i] . where )
        {
          case 0:	//no gap
          t  += mm; break;
          case 1:	//gap is in the text
          t  += mm - out[i] . min_gap; break;
          case 2:	//gap is in the pattern
          t  += mm + out[i] . min_gap; break;
	  default:
          break;
	}
	pb += mm;

        free ( frag );
    }

   return ( 1 );
}

/* Computes the optimal semi-global alignment between the optimal number of fragments of pattern p and text t */
unsigned int gapsmis_one_to_one_onf ( const char * p, const char * t, const struct gapmis_params * in, struct gapmis_align * out )
{
   unsigned int		i;
   unsigned int		m;
   unsigned int		pb;
   unsigned int		mm;
   unsigned int		num_frags     = in -> num_frags;
   unsigned int		opt_num_frags = in -> num_frags;
   struct gapmis_params in2;
   double		scr, max_scr;

   m = strlen ( p );
   memcpy ( &in2, in, sizeof ( struct gapmis_params ) );
   
   if ( num_frags == 0 || num_frags > m )
    {
      errno = FRAGS; //Error: the num of frags should be less or equal to the length of p!!!
      return ( 0 );
    }

   max_scr = 0;
   for ( i = 0; i < num_frags; ++ i )
    {
      in2 . num_frags = i + 1;
      scr = 0;
      gapsmis_one_to_one_scr ( p, t, &in2, &scr );
      if ( scr > max_scr )
      {
        opt_num_frags = i + 1;
        max_scr = scr;
      }
    }

   in2 . num_frags = opt_num_frags;
  
   mm = m / opt_num_frags;
   pb = 0;		//processed bytes

   for ( i = 0; i < opt_num_frags; ++ i )
    {
         char		* frag;

	/* calculate the length of the pattern's fragment */
	if ( i == opt_num_frags - 1 )  mm = m - pb;  

	/* allocate and copy the pattern's fragment */
        if ( ! ( frag = ( char * ) calloc ( mm + 1 , sizeof ( char ) ) ) )
         {
           errno = MALLOC; //Error: frag could not be allocated!!!
           return ( 0 );
         } 
        memcpy ( frag, p, mm );	
 	frag[mm] = '\0';

	/* perform the single-gap alignment */
        if ( ! ( gapmis_one_to_one ( frag, t, &in2, &out[i] ) ) )
           return ( 0 );
	
	/* update the pointers */
	p  += mm;
        switch ( out[i] . where )
        {
          case 0:	//no gap
          t  += mm; break;
          case 1:	//gap is in the text
          t  += mm - out[i] . min_gap; break;
          case 2:	//gap is in the pattern
          t  += mm + out[i] . min_gap; break;
	  default:
          break;
	}
	pb += mm;

        free ( frag );
    }

   return ( 1 );
}

/* Computes the score of the optimal semi-global alignment between a number of fragments of pattern p and text t */
unsigned int gapsmis_one_to_one_scr ( const char * p, const char * t, const struct gapmis_params * in, double * scr )
{
   unsigned int		i;
   unsigned int		m;
   unsigned int		pb;
   unsigned int		mm;
   unsigned int		num_frags = in -> num_frags; 

   m = strlen ( p );
   if ( num_frags == 0 || num_frags > m )
    {
      errno = FRAGS; //Error: the num offrags should be less or equal to the length of p!!!
      return ( 0 );
    }

   mm = m / num_frags;
   pb = 0;		//processed bytes
   ( *scr ) = 0;

   for ( i = 0; i < num_frags; ++ i )
    {
        char                 * frag;
        struct gapmis_align  out;

	/* calculate the length of the fragment */
	if ( i == num_frags - 1 )  mm = m - pb;  

	/* allocate and copy the fragment */
        if ( ! ( frag = ( char * ) calloc ( mm + 1 , sizeof ( char ) ) ) )
         {
           errno = MALLOC; //Error: frag could not be allocated!!!
           return ( 0 );
         } 
        memcpy ( frag, p, mm );	
	frag[mm] = '\0';

	/* do the single-gap alignment */
        /* Note: unfortunately here we cannot call gapmis_one_to_one_scr() because we need out . min_gap later on*/
        if ( ! ( gapmis_one_to_one ( frag, t, in, &out ) ) )
          return ( 0 );

	( *scr ) += out . max_score;

	/* update the pointers */
	p  += mm;
        switch ( out . where )
        {
          case 0:	//no gap
          t  += mm; break;
          case 1:	//gap is in the text
          t  += mm - out . min_gap; break;
          case 2:	//gap is in the pattern
          t  += mm + out . min_gap; break;
	  default:
          break;
	}
	pb += mm;

        free ( frag );
    }

   return ( 1 );
}

/* Computes the optimal semi-global alignment with the maximum score between each pattern p and all texts t */
unsigned int gapmis_many_to_many_opt ( const char const ** p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out )
 {
   for ( ; *p; ++ p, ++out )
    {
      if ( ! ( gapmis_one_to_many_opt ( *p, t, in, out ) ) )
        return ( 0 );  
    }

   return ( 1 );
 }

/* Computes the optimal semi-global alignment with the maximum score between a pattern p and all texts t */
unsigned int gapmis_one_to_many_opt ( const char * p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out )
 {
   double       tmp_scr     = -DBL_MAX;
   double       scr         = 0;
   unsigned int i           = 0;
   unsigned int max_t       = 0;
   const char   ** Tmp      = t;	
	
   for ( ; *Tmp; ++Tmp, ++ i )	//computing the alignment with the maximum score
    {
      if ( ! ( gapmis_one_to_one_scr ( p, *Tmp, in, &scr ) ) )
        return ( 0 );
      if ( scr > tmp_scr )
       {
         max_t = i;
         tmp_scr = scr;
       } 
    } 
   
   if ( ! ( gapmis_one_to_one ( p, t[ max_t ], in, out ) ) ) //computing the rest details of the alignment with the maximum score
     return ( 0 );

   return ( 1 );
 }

/* Computes only the maximum score of the optimal semi-global alignment between pattern p and text t */
unsigned int gapmis_one_to_one_scr ( const char * p, const char * t, const struct gapmis_params * in, double * scr )
 {
   int               ** G; 		//dynamic programming matrix
   unsigned int         n;
   unsigned int         m;
   unsigned int         i;

   /* Checks the input parameters */
   n = strlen ( t );
   m = strlen ( p );
   if ( m > n )
    {
      errno = LENGTH; //Error: the length of p should be less or equal to the length of t!!!
      return ( 0 );
    }
   
   if ( in -> scoring_matrix > 1 )
    {
      errno = MATRIX; //Error: the value of the scoring matrix parameter should be either 0 (nucleotide sequences) or 1 (protein sequences)!!!
      return ( 0 );
    }

   if ( in -> max_gap >= n )
    {
      errno = MAXGAP; //Error: the value of the max gap parameter should be less than the length of t!!!
      return ( 0 );
    }

   /* 2d dynamic memory allocation for matrices G and H*/
   if ( ! ( G = ( int ** ) malloc ( ( n + 1 ) * sizeof ( int * ) ) ) )
    {
      errno = MALLOC; //Error: DP matrix could not be allocated!!!
      return ( 0 );
    } 
   
   if ( ! ( G[0] = ( int * ) calloc ( ( n + 1 ) * ( m + 1 ), sizeof ( int ) ) ) )
    {
      errno = MALLOC; //Error: DP matrix could not be allocated!!!
      return ( 0 );
    } 
   
   for ( i = 1; i < n + 1; ++ i )
     G[i] = ( void * ) G[0] + i * ( m + 1 ) * sizeof ( int );
     
   /* dynamic programming algorithm */
   if ( ! ( dp_algorithm_scr( G, t, n, p, m, in ) ) )
    {
      //Error: dp_algorithm_scr() failed due to bad character!!!
      return ( 0 );	
    }
   
   /* computes the optimal alignment based on the matrix score and the affine gap penalty function */
   opt_solution_scr ( G, n, m, in, scr );
   
   free ( G[0] );
   free ( G );	
   return ( 1 );
   
 } 

#if 1
 
/* Computes the optimal semi-global alignment between each pattern p and all texts t */
unsigned int gapmis_many_to_many ( const char const ** p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out )
 {
   
   unsigned int         stride = 0;
   const char        ** Tmp;
   for ( Tmp = t; *Tmp; ++ Tmp, ++ stride );  //counting the total number of texts

   for ( ; *p; ++ p, out += stride )
    {
      if ( ! ( gapmis_one_to_many ( *p, t, in, out ) ) )
        return ( 0 );
    }

   return ( 1 );
 }

/* Computes the optimal semi-global alignment between a pattern p and all texts t */
unsigned int gapmis_one_to_many ( const char * p, const char const ** t, const struct gapmis_params * in, struct gapmis_align * out )
 {
   for ( ; *t; ++ t, ++ out )
    {
      if ( ! ( gapmis_one_to_one ( p, *t, in, out ) ) )
        return ( 0 );
    }

   return ( 1 );
 }

/* Computes the optimal semi-global alignment between pattern p and text t */
unsigned int gapmis_one_to_one ( const char * p, const char * t, const struct gapmis_params * in, struct gapmis_align * out )
 {
   int               ** G; 		//dynamic programming matrix
   unsigned int      ** H; 		//backtracing matrix
   unsigned int         n;
   unsigned int         m;
   unsigned int         i;
   unsigned int         start = 0;

   /* Checks the input parameters */
   n = strlen ( t );
   m = strlen ( p );
   if ( m > n )
    {
      errno = LENGTH; //Error: the length of p should be less or equal to the length of t!!!
      return ( 0 );
    }
   
   if ( in -> scoring_matrix > 1 )
    {
      errno = MATRIX; //Error: the value of the scoring matrix parameter should be either 0 (nucleotide sequences) or 1 (protein sequences)!!!
      return ( 0 );
    }

   if ( in -> max_gap >= n )
    {
      errno = MAXGAP; //Error: the value of the max gap parameter should be less than the length of t!!!
      return ( 0 );
    }

   /* 2d dynamic memory allocation for matrices G and H*/
   if ( ! ( G = ( int ** ) malloc ( ( n + 1 ) * sizeof ( int * ) ) ) )
    {
      errno = MALLOC; //Error: DP matrix could not be allocated!!!
      return ( 0 );
    } 
   
   if ( ! ( G[0] = ( int * ) calloc ( ( n + 1 ) * ( m + 1 ), sizeof ( int ) ) ) )
    {
      errno = MALLOC; //Error: DP matrix could not be allocated!!!
      return ( 0 );
    } 
   
   for ( i = 1; i < n + 1; ++ i )
     G[i] = ( void * ) G[0] + i * ( m + 1 ) * sizeof ( int );
   
   if ( ! ( H = ( unsigned int ** ) malloc ( ( n + 1 ) * sizeof ( unsigned int * ) ) ) )
    {
      errno = MALLOC; //Error: DP matrix could not be allocated!!!
      return ( 0 );
    } 
   
   if ( ! ( H[0] = ( unsigned int * ) calloc ( ( n + 1 ) * ( m + 1 ) , sizeof ( unsigned int ) ) ) )
    {
      errno = MALLOC; //Error: DP matrix could not be allocated!!!
      return ( 0 );
    }
   
   for ( i = 1 ; i < n + 1 ; ++ i )
     H[i] = ( void * ) H[0] + i * ( m + 1 ) * sizeof ( unsigned int );
   
   /* dynamic programming algorithm */
   if ( ! ( dp_algorithm( G, H, t, n, p, m, in ) ) )
    {
      //Error: dp_algorithm_scr() failed due to bad character!!!
      return ( 0 );	
    }
   
   /* computes the optimal alignment based on the matrix score and the affine gap penalty function */
   opt_solution ( G, n, m, in, out, &start );
   
   /* computes the position of the gap */
   if ( out -> min_gap > 0 ) 
        backtracing ( H, m, n, start, out );

   /* computes the number of mismatches in the alignment */
   if ( out -> where == 1 )     //gap is the text
     num_mismatch ( t, n, p, m, out );
   else                         //gap is in the pattern or there is no gap
     num_mismatch ( p, m, t, n, out );
   
   free ( G[0] );
   free ( H[0] );
   free ( G );
   free ( H );

   return ( 1 );
   
 } 
#endif
/* The dynamic programming algorithm for calculating matrices G and H */
unsigned int dp_algorithm ( int ** G, unsigned int ** H, const char * t, unsigned int n, const char * p, unsigned int m, const struct gapmis_params * in )
 {
   int                  gap;
   int                  mis;
   unsigned int         i;
   unsigned int         j;
   int                  matching_score;
   unsigned int 	j_min;
   unsigned int 	j_max;
   unsigned int 	valM;
   unsigned int 	i_max;

   for( i = 0; i < n + 1 ; i++ )      H[i][0] = i;
   for( j = 0; j < m + 1 ; j++ )      H[0][j] = j;

   i_max = GM_min ( n, m + in -> max_gap );

   for( i = 1; i < i_max + 1; i++)
     {
       j_min = GM_max ( 1, (int) ( i - in -> max_gap ));
       j_max = GM_min ( m, (int) ( i + in -> max_gap ));
       for( j = j_min; j <= j_max; j++ )
        {
           matching_score = ( in -> scoring_matrix ? ( int ) pro_delta( t[i - 1], p[j - 1] ) : ( int ) nuc_delta( t[i - 1], p[j - 1] ) );
           if ( matching_score == BADCHAR )
             {  
                errno = BADCHAR; //Error: unrecognizable character!!!
                return ( 0 );
             }	

           mis = G[i - 1][j - 1] + matching_score;
	   gap = G[j][j];
	   valM = i - j;

	   if( j > i )	
	     {
	       gap = G[i][i];
               valM = j - i;
	     }

           if( gap > mis )	H[i][j] = valM;
	   if( i == j )		gap = mis - 1;

           G[i][j] = GM_max ( mis, gap );
         }
      }
   return ( 1 );
 }

/* The dynamic programming algorithm for calculating only matrix G */
unsigned int dp_algorithm_scr ( int ** G, const char * t, unsigned int n, const char * p, unsigned int m, const struct gapmis_params * in )
 {
   int                  gap;
   int                  mis;
   unsigned int         i;
   unsigned int         j;
   int                  matching_score;
   unsigned int 	j_min;
   unsigned int 	j_max;
   unsigned int 	i_max;

   i_max = GM_min ( n, m + in -> max_gap );

   for( i = 1; i < i_max + 1; i++)
     {
       j_min = GM_max ( 1, (int) ( i - in -> max_gap ));
       j_max = GM_min ( m, (int) ( i + in -> max_gap ));
       for( j = j_min; j <= j_max; j++ )
        {
           matching_score = ( in -> scoring_matrix ? ( int ) pro_delta( t[i - 1], p[j - 1] ) : ( int ) nuc_delta( t[i - 1], p[j - 1] ) );
           if ( matching_score == BADCHAR )
            {  
               errno = BADCHAR; //Error: unrecognizable character!!!
               return ( 0 );
            }	

           mis = G[i - 1][j - 1] + matching_score;
	   gap = G[j][j];
	   if( j > i )		gap = G[i][i];
	   if( i == j )		gap = mis - 1;
           G[i][j] = GM_max ( mis, gap );
         }
      }
   return ( 1 );
 }

/* Returns the score for matching character a and b based on EDNAFULL matrix */
int nuc_delta ( char a, char b )
 {
   unsigned int index_a = nuc_char_to_index ( a );
   unsigned int index_b = nuc_char_to_index ( b );

   if ( ( index_a < GM_NUC_SCORING_MATRIX_SIZE ) && ( index_b < GM_NUC_SCORING_MATRIX_SIZE ) )
     return ( EDNAFULL_matrix[ index_a ][ index_b ] );
   else //Error: unrecognizable character!!!
     return ( BADCHAR );
 }

/* Returns the score for matching character a and b based on EBLOSUM62 matrix */
int pro_delta ( char a, char b )
 {
   unsigned int index_a = pro_char_to_index ( a );
   unsigned int index_b = pro_char_to_index ( b );

   if ( ( index_a < GM_PRO_SCORING_MATRIX_SIZE ) && ( index_b < GM_PRO_SCORING_MATRIX_SIZE ) )
     return ( EBLOSUM62_matrix[ index_a ][ index_b ] );
   else //Error: unrecognizable character!!!
     return ( BADCHAR );
 }

/* Returns the index of char a in EDNAFULL matrix */
unsigned int nuc_char_to_index ( char a )
 {
   unsigned int         index; 

   switch ( a )
    {
      case 'A':
        index = 0; break;

      case 'T':
        index = 1; break;

      case 'G':
        index = 2; break;

      case 'C':
        index = 3; break;

      case 'S':
        index = 4; break;

      case 'W':
        index = 5; break;

      case 'R':
        index = 6; break;

      case 'Y':
        index = 7; break;

      case 'K':
        index = 8; break;

      case 'M':
        index = 9; break;

      case 'B':
        index = 10; break;

      case 'V':
        index = 11; break;

      case 'H':
        index = 12; break;

      case 'D':
        index = 13; break;

      case 'N':
        index = 14; break;

      default:
        index = BADCHAR; break;
    }
   
   return ( index );
 }

/* Returns the index of char a in EBLOSUM62 matrix */
unsigned int pro_char_to_index ( char a )
 {
   unsigned int         index; 

   switch ( a )
    {
      case 'A':
        index = 0; break;

      case 'R':
        index = 1; break;

      case 'N':
        index = 2; break;

      case 'D':
        index = 3; break;

      case 'C':
        index = 4; break;

      case 'Q':
        index = 5; break;

      case 'E':
        index = 6; break;

      case 'G':
        index = 7; break;

      case 'H':
        index = 8; break;

      case 'I':
        index = 9; break;

      case 'L':
        index = 10; break;

      case 'K':
        index = 11; break;

      case 'M':
        index = 12; break;

      case 'F':
        index = 13; break;

      case 'P':
        index = 14; break;

      case 'S':
        index = 15; break;

      case 'T':
        index = 16; break;

      case 'W':
        index = 17; break;

      case 'Y':
        index = 18; break;

      case 'V':
        index = 19; break;

      case 'B':
        index = 20; break;

      case 'Z':
        index = 21; break;

      case 'X':
        index = 22; break;

      case '*':
        index = 23; break;

      default:
        index = BADCHAR; break;
    }
   return ( index );
 }

/* Computes the optimal alignment using matrix G */
unsigned int opt_solution ( int ** G, unsigned int n, unsigned int m, const struct gapmis_params * in, struct gapmis_align * out, unsigned int * start )
 {
   unsigned int         i;
   double               score = -DBL_MAX;
   unsigned int         up    = 0;
   unsigned int         down  = 0;
   double               temp_score;
   
   i_limits ( n, m, &up, &down, in -> max_gap );			// computes the i coordinates for matrix G for the last column

   for ( i = up ; i <= down ; i++ )
    {
      temp_score = 0.0;
      if ( i < m )
       {
         if ( m - i <= in -> max_gap )
          {
            temp_score = total_scoring ( m - i, G[i][m], - in -> gap_open_pen, - in -> gap_extend_pen );
            if ( temp_score > score )
             {
               score            = temp_score;
               out -> max_score = score; 
               out -> min_gap   = m - i;
               out -> where     = 1;		//where: gap is in the text and start backtracing from the last column
               ( *start )       = i;		//backtrace from cell G[start,m]
             }
          }
       }
      else if ( i > m )
       {
         if ( i - m <= in -> max_gap )
          {
            temp_score = total_scoring( i - m, G[i][m], - in -> gap_open_pen, - in -> gap_extend_pen );
            if (  temp_score > score )
             {
               score            = temp_score;
               out -> max_score = score; 
               out -> min_gap   = i - m;
               out -> where     = 2;		//where: gap is in the pattern and start backtracing from last column
               ( *start )       = i;		//backtrace from cell G[start,m]
             }
          }
       }
      else if ( i == m )
       {
         temp_score = total_scoring( 0, G[i][m], - in -> gap_open_pen, - in -> gap_extend_pen );
         if (  temp_score > score ) // mgap = 0 
          {
            score            = temp_score;
            out -> max_score = score; 
            out -> min_gap   = 0;
            out -> where     = 0;		//there is no gap
            ( *start )       = m;		//no need to backtrace
          }
       }
    }

   return ( 1 );
 }

/* Computes only the maximum score of the optimal alignment using matrix G */
unsigned int opt_solution_scr ( int ** G, unsigned int n, unsigned int m, const struct gapmis_params * in, double * scr )
 {
   unsigned int         i;
   double               score = -DBL_MAX;
   unsigned int         up    = 0;
   unsigned int         down  = 0;
   double               temp_score;
   
   i_limits ( n, m, &up, &down, in -> max_gap );			

   for ( i = up ; i <= down ; i++ )
    {
      temp_score = 0.0;
      if ( i < m )
       {
         if ( m - i <= in -> max_gap )
          {
            temp_score = total_scoring ( m - i, G[i][m], - in -> gap_open_pen, - in -> gap_extend_pen );
          }
       }
      else if ( i > m )
       {
         if ( i - m <= in -> max_gap )
          {
            temp_score = total_scoring( i - m, G[i][m], - in -> gap_open_pen, - in -> gap_extend_pen );        
          }
       }
      else if ( i == m )
       {
         temp_score = total_scoring( 0, G[i][m], - in -> gap_open_pen, - in -> gap_extend_pen );
       }
      if (  temp_score > score )
       {
         score     = temp_score;
         ( * scr ) = score;
       }
    }

   return ( 1 );
 }

/* Computes the limits of the i-th coordinate for the matrix G in constant time */
unsigned int i_limits ( unsigned int n, unsigned int m, unsigned int * up, unsigned int * down, unsigned int max_gap )
 {
   (* up )   = ( ( int ) m - ( int ) max_gap < 0 ) ?  0 : m - max_gap;
   (* down ) = ( m + max_gap > n )                 ?  n : m + max_gap;
   return ( 0 );
 }


/*
Gives the total score of an alignment in constant time
Note: double matrix_score is only the value of G[i][m], i.e. the score of an alignment WITHOUT the affine gap penalty
*/
double total_scoring( unsigned int gap, double matrix_score, double gap_open_penalty, double gap_extend_penalty )
 {
   return ( matrix_score + ( ( gap > 0 ) ? ( gap - 1 ) * gap_extend_penalty + gap_open_penalty : 0 ) );
 }

/* Computes the position of the gap */
unsigned int backtracing ( unsigned int ** H, unsigned int m, unsigned int n, unsigned int start, struct gapmis_align * out )
 {
   unsigned int         i, j;

   out -> gap_pos = 0;

   if ( out -> where == 1 || out -> where == 2 )
    {
      i = start; j = m; 	//we start backtracing from the last column
    }
   else
    {
      i = n; j = start;	//we start backtracing from the last row
    }
   while ( i >= 0 && j >= 0 )
    {
      if ( H[i][j] == 0 )
       {
         --i; --j;
       }
      else				
       {
         out -> gap_pos = ( i > j ) ? j : i;
         break;	
       }
    }

   return ( 1 );
 }

/* Computes the number of mismatches between seqa and seqb*/
unsigned int num_mismatch ( const char * seqa, unsigned int seqa_len, const char * seqb, unsigned int seqb_len, struct gapmis_align * out )
 {
   unsigned int i;
   unsigned int min_mis = 0;

   if ( out -> min_gap > 0 )
    {
      for ( i = 0; i < out -> gap_pos; ++ i )
        if ( seqa[i] != seqb[i] ) ++ min_mis;

      for ( ; i < seqb_len - out -> min_gap && i < seqa_len ; ++ i )
        if ( seqa[i] != seqb[i + out -> min_gap] ) ++ min_mis;
    }
   else
    {
      for ( i = 0; i < seqa_len; ++ i )
        if ( seqa[i] != seqb[i] ) ++ min_mis;
    }	
   out -> num_mis = min_mis;

   return ( 1 );
 }

/* Prints the header in the output file */
void print_header ( FILE * out, const char * filename, const struct gapmis_params* in )
 {
   time_t               t;
   time ( &t );

   fprintf ( out, "####################################\n" );
   fprintf ( out, "# Program: GapMis\n" );
   fprintf ( out, "# Rundate: %s", ctime ( &t ) );
   fprintf ( out, "# Report file: %s\n", filename );
   fprintf ( out, "# Matrix: %s\n", ( in -> scoring_matrix ? "BLOSUM62" : "EDNAFULL" ) );
   fprintf ( out, "# Gap penalty: %.3f\n", in -> gap_open_pen );
   fprintf ( out, "# Extend penalty: %.3f\n", in -> gap_extend_pen );
   fprintf ( out, "####################################\n\n" );
 }

/* Creates seq_gap and mark_mis, and computes min_mis */
unsigned int print_alignment ( const char * seqa, unsigned int seqa_len, const char * seqb, unsigned int seqb_len, char* seq_gap, char* mark_mis, const struct gapmis_params * in, struct gapmis_align * out )
 {
   unsigned int i, j;

   if ( out -> min_gap > 0 )
    {

      for ( i = 0; i < out -> gap_pos; ++ i )
       {
         seq_gap[i] = seqa[i];
         if ( seqa[i] != seqb[i] )	
          {
            mark_mis[i] = '.';
          }
         else				
           mark_mis[i] = '|';
       }

      for ( j = 0; j < out -> min_gap; ++ j )
       {
         seq_gap[ j + i ] = '-'; 
         mark_mis[ j + i ] = ' ';
       }

      for ( ; i < seqb_len - out -> min_gap && i < seqa_len ; ++ i )
       {
         seq_gap[j + i] = seqa[i];
         if ( seqa[i] != seqb[i + out -> min_gap] )	
          {
            mark_mis[ j + i ] = '.';
          }
         else
           mark_mis[j + i] = '|';
       }
      
      for ( ; i < seqa_len; ++ i )
       {
         seq_gap[j + i] = seqa[i];
         mark_mis[ j + i ] = '|';
       }
    }
   else
    {
      for ( i = 0; i < seqa_len; ++ i )
       {
         seq_gap[i] = seqa[i];
         if ( seqa[i] != seqb[i] )
          {
            mark_mis[i] = '.';
          }
         else			
           mark_mis[i] = '|';
       }
    }	
   
   return ( 1 );
 }

void print_line ( const char * s, int start, int offset, int stop, int * nr_gaps, int end, int diff, FILE * output, const char* header )
 {
   int                  k;

   if ( start == stop ) return;

   if ( diff )
    {
      fprintf ( output, "%25s", "" );
    }
   else
    {
     if ( header )
       fprintf ( output, "%-13.13s %10d ", header, start + 1 - *nr_gaps + offset );
     else
       fprintf ( output, "%-13.13s %10d ", "", start + 1 - *nr_gaps + offset );
    }

   for ( ; start < stop; ++ start )
    {
      fputc ( s[start], output );
      if ( s[start] == '-' && ! diff ) ++ ( *nr_gaps );
    }

   if ( stop != end )
    {
      for ( k = stop; k < end; ++ k )
       {
         fputc ( ' ', output );
       }
    }
   if ( ! diff )  fprintf ( output, " %-10d", start - *nr_gaps + offset );
   fprintf ( output, "\n" );
 }

/* Wrap two sequences s1 and s2 including the differences (diff) so that the line width is at most len */
void wrap ( const char * s1, const char * s1_header, const char * s2, const char * s2_header, const char * diff, int len, int s1_offset, int s2_offset, FILE * output )
 {
   int                  m, n, i, j;
   int                  nr_gaps_a;
   int                  nr_gaps_b;
   int                  nr_lines;

   if ( ! len ) 
     return;

   m = strlen ( s1 );
   n = strlen ( s2 );

   if ( ! n && ! m ) 
     return;

   i         = 0;
   j         = 0;
   nr_gaps_a = 0;
   nr_gaps_b = 0;

   //nr_lines = m / len;
   nr_lines = ( n > m ? m : n ) / len;
   for ( i = 0; i < nr_lines; ++ i )
    {
      /* Sequence s1 */
      print_line ( s1, i * len, s1_offset, ( i + 1 ) * len, &nr_gaps_a, ( i + 1 ) * len, 0, output, s1_header );

      /* Difference */
      print_line ( diff, i * len, 0, ( i + 1 ) * len, NULL, ( i + 1 ) * len, 1, output, NULL );

      /* Sequence s2 */
      print_line ( s2 , i * len, s2_offset, ( i + 1 ) * len, &nr_gaps_b, ( i + 1 ) * len, 0, output, s2_header );
      fprintf ( output, "\n" );
    }

   /* Last line of first sequence and difference */
   j = i * len;
   if ( j < m || j < n ) 
    {
      print_line ( s1,  i * len , s1_offset, GM_min ( m, n ), &nr_gaps_a, ( i + 1 ) * len, 0, output, s1_header );
      print_line ( diff,  i * len , 0, ( m < n ) ? m : n, NULL, ( i + 1 ) * len, 1, output, NULL );
      print_line ( s2,  i * len , s2_offset, GM_min ( n, m), &nr_gaps_b, ( i + 1 ) * len, 0, output, s2_header );
    }

 }

/* Creates the output file with the one to one alignment */
unsigned int gapmis_results_one_to_one ( const char * filename, const char * p, const char * p_header, const char * t, const char * t_header, const struct gapmis_params* in, struct gapmis_align* out )
 {

   FILE          * output;
   char          * seq_gap;            //the sequence with the inserted gap 
   char          * mark_mis;           //a string with the mismatches marked as '|' (and the matches as ' ')
   unsigned int    n;
   unsigned int    m;
 
   n = strlen ( t );
   m = strlen ( p );
   if ( m > n )
    {
      errno = LENGTH; //Error: the length of p should be less or equal to the length of t!!!
      return ( 0 );
    }
   
   /* Dynamic memory allocation for seq_gap */
   if ( ! ( seq_gap = ( char * ) calloc ( n + 1 + ( out -> min_gap ) + 1, sizeof ( char ) ) ) )
    {
      errno = MALLOC; //Error: seq_gap could not be allocated!!!
      return ( 0 );
    } 
   
   if ( ! ( mark_mis = ( char* ) calloc ( n + ( out -> min_gap ) + 1, sizeof( char ) ) ) )
    {
      errno = MALLOC; //Error: mark_mis could not be allocated!!!
      return ( 0 );
    } 
   
   if ( ! ( output = fopen ( filename, "w" ) ) )
    {
      errno = IO; 
      return ( 0 );
    }
   
   print_header ( output, filename, in );
   
   if ( out -> where == 1 ) //gap is in the text
    {
      print_alignment ( t, n, p, m, seq_gap, mark_mis, in, out );
      wrap ( p, p_header, seq_gap, t_header, mark_mis, GM_LINE_LNG, 0, 0, output ); 
    }
   else                     //gap is in the pattern
    {
      print_alignment ( p, m, t, n, seq_gap, mark_mis, in, out );
      wrap ( seq_gap, p_header, t, t_header, mark_mis, GM_LINE_LNG, 0, 0, output ); 
    }
   
   fprintf ( output, "\n" );
   fprintf ( output, "Alignment score: %lf\n", out -> max_score );
   fprintf ( output, "Number of mismatches: %d\n", out -> num_mis );
   fprintf ( output, "Length of gap: %d\n", out -> min_gap );
   
   if( out -> min_gap > 0 )
    {
      fprintf ( output, "The gap was inserted in %.13s after position %d\n", ( out -> where == 1 ) ? t_header : p_header, out -> gap_pos );
    } 
   
   fprintf ( output, "\n\n" );
   
   if ( fclose ( output ) ) 
    {
      errno = IO; 
      return ( 0 );
    }
   
   free ( mark_mis );
   free ( seq_gap );
   
   return ( 1 );	
 }

/* Creates the output file with the one to many alignments */
unsigned int gapmis_results_one_to_many ( const char * filename, const char * p, const char const * p_header, const char const ** t, const char const ** t_header, const struct gapmis_params* in, struct gapmis_align* out )
 {

   FILE          * output;
   char          * seq_gap;            //the sequence with the inserted gap 
   char          * mark_mis;           //a string with the mismatches marked as '|' (and the matches as ' ')
   unsigned int    n;
   unsigned int    m;
  
   m = strlen ( p );
   if ( ! ( output = fopen ( filename, "w" ) ) )
    {
      errno = IO; 
      return ( 0 );
    }
   
   for ( ; *t && *t_header; ++ t, ++ t_header, ++ out )
    {
       n = strlen ( *t );
       
       if ( m > n )
        {
          errno = LENGTH; //Error: the length of p should be less or equal to the length of t!!!
          return ( 0 );
        }
   
      /* Dynamic memory allocation for seq_gap */
      if ( ! ( seq_gap = ( char * ) calloc ( n + 1 + ( out -> min_gap ) + 1, sizeof ( char ) ) ) )
       {
         errno = MALLOC; //Error: seq_gap could not be allocated!!!
         return ( 0 );
       } 
   
      if ( ! ( mark_mis = ( char* ) calloc ( n + ( out -> min_gap ) + 1, sizeof( char ) ) ) )
       {
         errno = MALLOC; //Error: mark_mis could not be allocated!!!
         return ( 0 );
       } 
      print_header ( output, filename, in );
   
      if ( out -> where == 1 ) //gap is in the text
       {
         print_alignment ( *t, n, p, m, seq_gap, mark_mis, in, out );
         wrap ( p, p_header, seq_gap, *t_header, mark_mis, GM_LINE_LNG, 0, 0, output ); 
       }
      else                     //gap is in the pattern
       {
         print_alignment ( p, m, *t, n, seq_gap, mark_mis, in, out );
         wrap ( seq_gap, p_header, *t, *t_header, mark_mis, GM_LINE_LNG, 0, 0, output ); 
       }
   
      fprintf ( output, "\n" );
      fprintf ( output, "Alignment score: %lf\n", out -> max_score );
      fprintf ( output, "Number of mismatches: %d\n", out -> num_mis );
      fprintf ( output, "Length of gap: %d\n", out -> min_gap );
   
      if( out -> min_gap > 0 )
       {
         fprintf ( output, "The gap was inserted in %.13s after position %d\n", ( out -> where == 1 ) ? *t_header : p_header, out -> gap_pos );
       } 
   
      fprintf ( output, "\n\n" );
      free ( mark_mis );
      free ( seq_gap ); 
    }

   if ( fclose ( output ) ) 
    {
      errno = IO; 
      return ( 0 );
    }
     
   return ( 1 );	
 }

/* Creates the output file with the many to many alignments */
unsigned int gapmis_results_many_to_many ( const char * filename, const char const ** p, const char const ** p_header, const char const ** t, const char const ** t_header, const struct gapmis_params* in, struct gapmis_align* out )
 {

   FILE          * output;
   char          * seq_gap;            //the sequence with the inserted gap 
   char          * mark_mis;           //a string with the mismatches marked as '|' (and the matches as ' ')
   unsigned int    n;
   unsigned int    m;
   const char   ** Tmp;
   const char   ** Tmp_header;

 
   if ( ! ( output = fopen ( filename, "w" ) ) )
    {
      errno = IO; 
      return ( 0 );
    }

   for ( ; *p && *p_header; ++ p, ++ p_header )
    {
      Tmp = t;
      Tmp_header = t_header;
      for ( ; *Tmp && *Tmp_header; ++ Tmp, ++ Tmp_header, ++ out )
       {
         n = strlen ( *Tmp );
         m = strlen ( *p );
         if ( m > n )
          {
            errno = LENGTH; //Error: the length of p should be less or equal to the length of t!!!
            return ( 0 );
          }
   
         /* Dynamic memory allocation for seq_gap */
         if ( ! ( seq_gap = ( char * ) calloc ( n + 1 + ( out -> min_gap ) + 1, sizeof ( char ) ) ) )
          {
            errno = MALLOC; //Error: seq_gap could not be allocated!!!
            return ( 0 );
          } 
   
         if ( ! ( mark_mis = ( char* ) calloc ( n + ( out -> min_gap ) + 1, sizeof( char ) ) ) )
          {
            errno = MALLOC; //Error: mark_mis could not be allocated!!!
            return ( 0 );
          } 
         print_header ( output, filename, in );
   
         if ( out -> where == 1 ) //gap is in the text
          {
            print_alignment ( *Tmp, n, *p, m, seq_gap, mark_mis, in, out );
            wrap ( *p, *p_header, seq_gap, *Tmp_header, mark_mis, GM_LINE_LNG, 0, 0, output ); 
          }
         else                     //gap is in the pattern
          {
            print_alignment ( *p, m, *Tmp, n, seq_gap, mark_mis, in, out );
            wrap ( seq_gap, *p_header, *Tmp, *Tmp_header, mark_mis, GM_LINE_LNG, 0, 0, output ); 
          }
   
         fprintf ( output, "\n" );
         fprintf ( output, "Alignment score: %lf\n", out -> max_score );
         fprintf ( output, "Number of mismatches: %d\n", out -> num_mis );
         fprintf ( output, "Length of gap: %d\n", out -> min_gap );
   
         if( out -> min_gap > 0 )
          {
            fprintf ( output, "The gap was inserted in %.13s after position %d\n", ( out -> where == 1 ) ? *Tmp_header : *p_header, out -> gap_pos );
          } 
   
         fprintf ( output, "\n\n" );
         free ( mark_mis );
         free ( seq_gap ); 
    }
  }
  if ( fclose ( output ) ) 
   {
     errno = IO; 
     return ( 0 );
   }
     
   return ( 1 );	
 }

/* Creates the output file with the one to one alignment */
unsigned int gapsmis_results_one_to_one ( const char * filename, const char * p, const char * p_header, const char * t, const char * t_header, const struct gapmis_params* in, struct gapmis_align* out )
 {

   FILE          * output;
   unsigned int  m;
   unsigned int  n;

   unsigned int  p_offset = 0;
   unsigned int  t_offset = 0;
   unsigned int	 i;
   unsigned int	 pb;
   unsigned int	 mm;
   double        total_max_score = 0;
   unsigned int  total_num_mis = 0;
   unsigned int  total_min_gap = 0;
   unsigned int  total_num_gaps = 0;
   unsigned int	 num_frags = in -> num_frags;

   m = strlen ( p );
   if ( num_frags == 0 || num_frags > m )
    {
      errno = FRAGS; //Error: the num of frags should be less or equal to the length of p!!!
      return ( 0 );
    }

   n = strlen ( t );
   if ( m > n )
    {
      errno = LENGTH; //Error: the length of p should be less or equal to the length of t!!!
      return ( 0 );
    }

   if ( ! ( output = fopen ( filename, "w" ) ) )
    {
      errno = IO; 
      return ( 0 );
    }
   
   print_header ( output, filename, in );

   mm = m / num_frags;
   pb = 0;		//processed bytes

   for ( i = 0; i < num_frags; ++ i, ++ out )
    { 
        char          * seq_gap;            // the sequence with the inserted gap 
        char          * mark_mis;           // a string with the matches marked as '|' and the mismatches as '.'
        char          * frag;		    // the fragment

	/* calculate the length of the current fragment */
	if ( i == num_frags - 1 )  mm = m - pb;  

	/* allocate and copy the fragment */
        if ( ! ( frag = ( char * ) calloc ( mm + 1 , sizeof ( char ) ) ) )
         {
           errno = MALLOC; //Error: seq_gap could not be allocated!!!
           return ( 0 );
         } 
        memcpy ( frag, p, mm );	
	frag[mm] = '\0';

        /* allocation for seq_gap and mark_mis */
        if ( ! ( seq_gap = ( char * ) calloc ( mm + 1 + ( out -> min_gap ) + 1, sizeof ( char ) ) ) )
         {
           errno = MALLOC; //Error: seq_gap could not be allocated!!!
           return ( 0 );
         } 
   
        if ( ! ( mark_mis = ( char* ) calloc ( mm + ( out -> min_gap ) + 1, sizeof( char ) ) ) )
         {
           errno = MALLOC; //Error: mark_mis could not be allocated!!!
           return ( 0 );
         }

	/* write the results for this fragment */
        if ( out -> where == 1 ) //gap is in the text
          {
            print_alignment ( t, mm, frag, mm, seq_gap, mark_mis, in, out );
            wrap ( frag, p_header, seq_gap, t_header, mark_mis, GM_LINE_LNG, p_offset, t_offset, output ); 
            t_offset += mm;
            p_offset += mm + out -> min_gap; //gap is in the text so we have to increase pattern's index
          }
        else                     //gap is in the pattern
          {
            print_alignment ( frag, mm, t, mm + out -> min_gap, seq_gap, mark_mis, in, out );
            wrap ( seq_gap, p_header, t, t_header, mark_mis, GM_LINE_LNG, p_offset, t_offset, output ); 
            t_offset += mm + out -> min_gap;  //gap is in the pattern so we have to increase text's index
            p_offset += mm;
          } 
	
	/* update the pointers */
	p  += mm;
        switch ( out -> where )
        {
          case 0:	//no gap
          t  += mm; break;
          case 1:	//gap is in the text
          t  += mm - out -> min_gap; break;
          case 2:	//gap is in the pattern
          t  += mm + out -> min_gap; break;
	  default:
          break;
	}
	pb += mm;

        total_max_score += out -> max_score;
        total_num_mis   += out -> num_mis;
        total_min_gap   += out -> min_gap;
        total_num_gaps  += ( out -> min_gap > 0 ? 1 : 0 );
   
        free ( frag );
        free ( mark_mis );
        free ( seq_gap );
    }

   fprintf ( output, "\n" );
   fprintf ( output, "Alignment score: %lf\n", total_max_score );
   fprintf ( output, "Number of mismatches: %d\n", total_num_mis );
   fprintf ( output, "Number of gap(s): %d\n", total_num_gaps );
   fprintf ( output, "Total length of gap(s): %d\n", total_min_gap );
   fprintf ( output, "Number of fragment(s): %d\n", num_frags );
   
   fprintf ( output, "\n\n" );
   
   if ( fclose ( output ) ) 
    {
      errno = IO; 
      return ( 0 );
    }
   
   return ( 1 );	
 }
