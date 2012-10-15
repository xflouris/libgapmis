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

#ifndef         GAPMIS_H
#define         GAPMIS_H

#ifdef _USE_GPU
#include "CL/opencl.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

//#ifndef _HACK_CLEAN_NAMESPACE /* these defines partly collide with my implementation and/or c++ in general. */

#ifndef __cplusplus
#define GM_max(a,b) ((a) > (b)) ? (a) : (b)
#define GM_min(a,b) ((a) < (b)) ? (a) : (b)
#endif

#define GM_NUC_SCORING_MATRIX_SIZE	15		
#define GM_PRO_SCORING_MATRIX_SIZE      24
#define GM_LINE_LNG 			53

//#endif

struct gapmis_params
 {
   unsigned int         max_gap;
   unsigned int         scoring_matrix;
   double               gap_open_pen;
   double               gap_extend_pen;
   unsigned int		num_frags;
 };

struct gapmis_align
 {
   double               max_score;
   unsigned int         min_gap;
   unsigned int         where;
   unsigned int         gap_pos;
   unsigned int         num_mis;
 };

/* Computes the optimal semi-global alignment between a number of fragments of pattern p and text t */
unsigned int gapsmis_one_to_one ( const char * p, const char * t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment between the optimal number of fragments of pattern p and text t */
unsigned int gapsmis_one_to_one_onf ( const char * p, const char * t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes only the maximum score of the optimal semi-global alignment between a number of fragments of pattern p and text t */
unsigned int gapsmis_one_to_one_scr ( const char * p, const char * t, const struct gapmis_params * in, double* scr );

/* Computes only the maximum score of the optimal semi-global alignment between pattern p and text t */
unsigned int gapmis_one_to_one_scr ( const char * p, const char * t, const struct gapmis_params * in, double* scr );

/* Computes the optimal semi-global alignment with the maximum score between a pattern p and all texts t */
unsigned int gapmis_one_to_many_opt ( const char * p, const char  ** t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment with the maximum score between each pattern p and all texts t */
unsigned int gapmis_many_to_many_opt ( const char ** p, const char ** t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment between pattern p and text t */
unsigned int gapmis_one_to_one ( const char * p, const char * t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment between a pattern p and all texts t */
unsigned int gapmis_one_to_many ( const char * p, const char  ** t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment between each pattern p and all texts t */
unsigned int gapmis_many_to_many ( const char ** p, const char ** t, const struct gapmis_params * in, struct gapmis_align * out );

/* Creates the output file with the many to many alignments */
unsigned int gapmis_results_many_to_many ( const char * filename, const char ** p, const char ** p_header, const char ** t, const char ** t_header, const struct gapmis_params * in, struct gapmis_align * out );

/* Creates the output file with the one to many alignments */
unsigned int gapmis_results_one_to_many ( const char * filename, const char * p, const char * p_header, const char ** t, const char ** t_header, const struct gapmis_params * in, struct gapmis_align * out );

/* Creates the output file with the one to one alignment */
unsigned int gapmis_results_one_to_one ( const char * filename, const char * p, const char * p_header, const char * t, const char * t_header, const struct gapmis_params * in, struct gapmis_align * out );

/* Creates the output file with the one to one alignment */
unsigned int gapsmis_results_one_to_one ( const char * filename, const char * p, const char * p_header, const char * t, const char * t_header, const struct gapmis_params * in, struct gapmis_align * out );

unsigned int dp_algorithm ( int ** G, unsigned int** H, const char* t, unsigned int n, const char * p, unsigned int m, const struct gapmis_params * in );
unsigned int dp_algorithm_scr ( int ** G, const char* t, unsigned int n, const char * p, unsigned int m, const struct gapmis_params * in );
int nuc_delta ( char a, char b );
int pro_delta ( char a, char b );
unsigned int nuc_char_to_index ( char a );
unsigned int pro_char_to_index ( char a );
unsigned int i_limits( unsigned int n, unsigned int m, unsigned int * up, unsigned int * down, unsigned int max_gap );
unsigned int opt_solution ( int ** G, unsigned int n, unsigned int m, const struct gapmis_params * in, struct gapmis_align * out, unsigned int * start );
unsigned int opt_solution_scr ( int ** G, unsigned int n, unsigned int m, const struct gapmis_params * in, double * scr );
unsigned int backtracing ( unsigned int ** H, unsigned int m, unsigned int n, unsigned int start, struct gapmis_align * out );
double total_scoring( unsigned int gap, double matrix_score, double gap_open_penalty, double gap_extend_penalty );
unsigned int num_mismatch ( const char * seqa, unsigned int seqa_len, const char * seqb, unsigned int seqb_len, struct gapmis_align * out );

#ifdef _USE_GPU
/* Computes the optimal semi-global alignment with the maximum score between a pattern p and all texts t */
unsigned int gapmis_one_to_many_opt_gpu ( const char * p, const char ** t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment with the maximum score between each pattern p and all texts t */
unsigned int gapmis_many_to_many_opt_gpu (const char ** p, const char ** t, const struct gapmis_params * in, struct gapmis_align * out );

cl_platform_id get_gpu_id(int * error);
cl_device_id get_dev_id(cl_platform_id gpu_id, int * error);
cl_context create_context(cl_device_id dev_id, int * error);
cl_command_queue create_cmd_queue (cl_device_id dev_id, cl_context context, int * error);
cl_kernel load_kernel (char * name, char * kernel_name, cl_device_id dev_id, cl_context context, int * error);
cl_mem malloc_device (cl_context context, size_t size, int * error);
void init_device_mem_int (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, int * mem, size_t size, int * error);
void init_device_mem_uint (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, unsigned int * mem, size_t size, int * error);
void init_device_mem_float (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, float * mem, size_t size, int * error);
void set_kernel_arguments (cl_kernel kernel, cl_command_queue cmd_queue, cl_mem cl_mem0, cl_mem cl_mem1, cl_mem cl_mem2, cl_mem cl_mem3, cl_mem cl_mem4, cl_mem cl_mem5, cl_mem cl_mem6, cl_mem cl_mem7);
void read_device_mem_float (cl_command_queue cmd_queue, size_t size, float * mem, cl_mem dev_mem, int * error);
void fill_txtsLenVec (unsigned int totalTxts, const char ** in, int * outVec);
void fill_argsVec (unsigned int totalPats, unsigned int totalTxts, const char ** in, unsigned int max_gap, unsigned int pBlockSize, unsigned int maxPatLen, unsigned int maxTxtLen, int * outVec);
void fill_patsVec (unsigned int total, unsigned int blockSize, const char ** in, unsigned int * outVec, int matrix);
void fill_txtsVec (unsigned int total, unsigned int blockSize, const char ** in, unsigned int * outVec, int matrix);
unsigned int get_pblock_size (unsigned int input, unsigned int mult);
unsigned int get_max_length (unsigned int total, const char ** input);
unsigned int get_min_length (unsigned int total, const char ** input);
unsigned int get_number_of_sequences (const char ** input);
unsigned int get_number_of_groups (int elements, int groupSize);
unsigned int kernel_launch (cl_kernel kernel, cl_context context, cl_command_queue cmd_queue, const char ** p, const char ** t, const struct gapmis_params * in, float * scores);
unsigned int kernel_launch_l (cl_kernel kernel, cl_context context, cl_command_queue cmd_queue, const char ** p, const char ** t, const struct gapmis_params * in, float * scores, struct gapmis_align * out);
void update_group_match (float * groupScores, int * groupMatch, float * groupMatchScores, unsigned int patGroupSize, unsigned int txtGroupSize, int pats, int txts, int patGroupIndex, int txtGroupIndex);
void set_invalid(int * groupMatch, int groupSize);
void set_minimum(float * groupMatchScores, int groupSize);
void set_null (const char ** input, int size);
void initialize_pointers (const char * groupPatterns[], int groupIndex, int groupSize, const char ** source, int sourceSize);
int check_sequences(unsigned int total, const char ** sequences, unsigned int matrix);

#endif


#ifdef _USE_SSE
#if 0 
/*
 * the non *_opt functions for SSE are kind of a hack. Don't use them.
 */

/* Computes the optimal semi-global alignment between pattern p and text t */
unsigned int gapmis_one_to_one_sse ( const char * p, const char * t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment between a pattern p and all texts t */
unsigned int gapmis_one_to_many_sse ( const char * p, const char ** t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment between each pattern p and all texts t */
unsigned int gapmis_many_to_many_sse ( const char ** p, const char ** t, const struct gapmis_params * in, struct gapmis_align * out );
#endif

/* Computes the optimal semi-global alignment with the maximum score between a pattern p and all texts t */
unsigned int gapmis_one_to_many_opt_sse ( const char * p, const char * const * t, const struct gapmis_params * in, struct gapmis_align * out );

/* Computes the optimal semi-global alignment between a set of factors and a set of patterns */
unsigned int gapmis_many_to_many_opt_sse ( const char * const * p, const char * const * t, const struct gapmis_params * in, struct gapmis_align * out );

/* Sets the number of threads to be used */
void gapmis_sse_hint_num_threads( size_t num );

#endif

#ifdef __cplusplus
}
#endif



#endif
