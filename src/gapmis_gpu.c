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

#ifdef _USE_GPU

unsigned int gapmis_one_to_many_opt_gpu ( const char * p1, const char ** t, const struct gapmis_params * in, struct gapmis_align * out )
{
	const char * p[] = { p1, NULL};

	if ( in -> scoring_matrix > 1 )
	{
		errno = MATRIX;
		return ( 0 );
	}

	unsigned int 	pats = get_number_of_sequences (p);
	unsigned int 	txts = get_number_of_sequences (t);
	unsigned int	maxPatLen = get_max_length (pats, p);
	unsigned int	minTxtLen = get_min_length (txts, t);

	if (check_sequences(pats,p,in->scoring_matrix)==0)
	{
		errno = BADCHAR;
      		return ( 0 );
	}

	if (check_sequences(txts,t,in->scoring_matrix)==0)
	{
		errno = BADCHAR;
      		return ( 0 );
	}

	if(maxPatLen > minTxtLen)
	{
		errno = LENGTH;
      		return ( 0 );
	}

	if ( in -> max_gap >= minTxtLen )
	{
		errno = MAXGAP; 
		return ( 0 );
	}

	int err = -1;

	/* get the GPU id */
	cl_platform_id gpu_id = get_gpu_id(&err);	
	if(err)
	{	
	 	errno = NOGPU;
      		return ( 0 );
	}

        /* get the device id */
	cl_device_id dev_id = get_dev_id(gpu_id, &err);
	if(err)
	{	
	 	errno = NOGPU;
      		return ( 0 );
	}

	/* create the context using dev_id */
	cl_context context = create_context(dev_id, &err);
	if(err)
	{	
	 	errno = GPUERROR;
      		return ( 0 );
	}

	/* create a list with the commands to be executed by GPU */
	cl_command_queue cmd_queue = create_cmd_queue (dev_id, context, &err);
	if(err)
	{	
	 	errno = GPUERROR;
      		return ( 0 );
	}

	/* create a kernel */
	cl_kernel kernel;

	/* load the kernel ``kernel_dna.cl'' with name ``gapmis_kernel''*/
	if(in->scoring_matrix==0)
		kernel = load_kernel ("kernel_dna.cl", "gapmis_kernel", dev_id, context, &err);
	else
		kernel = load_kernel ("kernel_pro.cl", "gapmis_kernel", dev_id, context, &err);

	if(err)
	{	
	 	errno = KERNEL;
      		return ( 0 );
	}

	const unsigned int patGroupSize = 1;
	const unsigned int txtGroupSize = 768;
	unsigned int i, j;	
	unsigned int patGroups = get_number_of_groups (pats, patGroupSize);
	unsigned int txtGroups = get_number_of_groups (txts, txtGroupSize);	

	const char * groupPatterns[patGroupSize+1];
	set_null (groupPatterns, patGroupSize+1);

	const char * groupTexts[txtGroupSize+1];
	set_null (groupTexts, txtGroupSize+1);

	float * groupScores;
        groupScores = calloc (patGroupSize*txtGroupSize, sizeof(float) );

	int groupMatch [patGroupSize];
	float groupMatchScores [patGroupSize];
	set_invalid(groupMatch,patGroupSize);
	set_minimum(groupMatchScores,patGroupSize);	

	for(i=0;i<patGroups;i++)
	{
		set_null (groupPatterns, patGroupSize+1);
		initialize_pointers (groupPatterns,i,patGroupSize,p,pats);
		set_invalid(groupMatch,patGroupSize);
		set_minimum(groupMatchScores,patGroupSize);
		
		for(j=0;j<txtGroups;j++)
		{			
			set_null (groupTexts, txtGroupSize+1);
			initialize_pointers (groupTexts,j,txtGroupSize,t,txts);

			if( ! ( kernel_launch (kernel, context, cmd_queue, groupPatterns, groupTexts, in, groupScores) ))
				return ( 0 );			

			update_group_match (groupScores,groupMatch,groupMatchScores,patGroupSize,txtGroupSize, pats, txts, i, j);
		
		}

		for(j=0;j<patGroupSize;j++)
		{
			if(i*patGroupSize+j<pats)
			{
				groupPatterns[0] = p[i*patGroupSize+j];
				groupPatterns[1] = NULL;

				groupTexts[0] = t[groupMatch[j]];
				groupTexts[1] = NULL;
				
				if( !( kernel_launch_l (kernel, context, cmd_queue, groupPatterns, groupTexts, in, groupScores,&out[i*patGroupSize+j] ) ) )
					return ( 0 );				
			}
		}
	}

        free ( groupScores );
        clReleaseContext ( context );
	clReleaseCommandQueue ( cmd_queue );
        clReleaseKernel(kernel);

	return ( 1 );
 }

unsigned int gapmis_many_to_many_opt_gpu ( const char ** p, const char ** t, const struct gapmis_params * in, struct gapmis_align * out )
{

	if ( in -> scoring_matrix > 1 )
	{
		errno = MATRIX;
		return ( 0 );
	}

	unsigned int pats = get_number_of_sequences (p);
	unsigned int txts = get_number_of_sequences (t);
	unsigned int	maxPatLen = get_max_length (pats, p);
	unsigned int	minTxtLen = get_min_length (txts, t);

	if (check_sequences(pats,p,in->scoring_matrix)==0)
	{
		errno = BADCHAR;
      		return ( 0 );
	}

	if (check_sequences(txts,t,in->scoring_matrix)==0)
	{
		errno = BADCHAR;
      		return ( 0 );
	}

	if(maxPatLen > minTxtLen)
	{
		errno = LENGTH;
      		return ( 0 );
	}

	if ( in -> max_gap >= minTxtLen )
	{
		errno = MAXGAP; 
		return ( 0 );
	}


	int err = -1;

	cl_platform_id  gpu_id = get_gpu_id(&err);
	
	if(err)
	{	
	 	errno = NOGPU;
      		return ( 0 );
	}

	cl_device_id  dev_id = get_dev_id(gpu_id, &err);

	if(err)
	{	
	 	errno = NOGPU;
      		return ( 0 );
	}

	cl_context context = create_context(dev_id, &err);
	if(err)
	{	
	 	errno = GPUERROR;
      		return ( 0 );
	}

	cl_command_queue cmd_queue = create_cmd_queue (dev_id, context, &err);
	if(err)
	{	
	 	errno = GPUERROR;
      		return ( 0 );
	}

	cl_kernel kernel;

	if(in->scoring_matrix==0)
		kernel = load_kernel ("kernel_dna.cl", "gapmis_kernel", dev_id, context, &err);
	else
		kernel = load_kernel ("kernel_pro.cl", "gapmis_kernel", dev_id, context, &err);

	if(err)
	{	
	 	errno = KERNEL;
      		return ( 0 );
	}

	const unsigned int patGroupSize = 1024;
	const unsigned int txtGroupSize = 32;
	unsigned int i, j;	
	unsigned int patGroups = get_number_of_groups (pats, patGroupSize);
	unsigned int txtGroups = get_number_of_groups (txts, txtGroupSize);	

	const char * groupPatterns[patGroupSize+1];
	set_null (groupPatterns, patGroupSize+1);

	const char * groupTexts[txtGroupSize+1];
	set_null (groupTexts, txtGroupSize+1);

	float * groupScores;
        groupScores = calloc (patGroupSize*txtGroupSize, sizeof(float) );

	int groupMatch [patGroupSize];
	float groupMatchScores [patGroupSize];
	set_invalid(groupMatch,patGroupSize);
	set_minimum(groupMatchScores,patGroupSize);	
      
	for(i=0;i<patGroups;i++)
	{
		set_null (groupPatterns, patGroupSize+1);
		initialize_pointers (groupPatterns,i,patGroupSize,p,pats);
		set_invalid(groupMatch,patGroupSize);
		set_minimum(groupMatchScores,patGroupSize);
		
		for(j=0;j<txtGroups;j++)
		{			
			set_null (groupTexts, txtGroupSize+1);
			initialize_pointers (groupTexts,j,txtGroupSize,t,txts);

			if( ! ( kernel_launch (kernel, context, cmd_queue, groupPatterns, groupTexts, in, groupScores) ) )
				return (0);			
			
			update_group_match (groupScores, groupMatch, groupMatchScores, patGroupSize, txtGroupSize, pats, txts, i, j);
		
		}

		for(j=0;j<patGroupSize;j++)
		{
			if(i*patGroupSize+j<pats)
			{
				groupPatterns[0] = p[i*patGroupSize+j];
				groupPatterns[1] = NULL;

				groupTexts[0] = t[groupMatch[j]];
				groupTexts[1] = NULL;
				
				if( !( kernel_launch_l (kernel, context, cmd_queue, groupPatterns, groupTexts, in, groupScores, &out[i*patGroupSize+j]) ) )
					return (0);				
			}
		}
	}

        free ( groupScores );
        clReleaseContext ( context );
	clReleaseCommandQueue ( cmd_queue );
        clReleaseKernel(kernel);
	
        return ( 1 );
 }

/* runs on the CPU */
void kernel_one_to_one_dna (unsigned int groupID, unsigned int localID, unsigned int * patsVec, unsigned int * txtsVec, int * argsVec, int * txtsLenVec, float * pensVec, int * hproVec, int * dproVec, float * scrsVec, unsigned int ** H, struct gapmis_align * out, int * start)
{
	int i, j, cur_diag_nxt, mis, gap,
	    max_gap = argsVec[2], 
	    blockSize = argsVec[3],
	    maxPatLen = argsVec[4],
	    dproVecGsize = argsVec[5],
	    hproVecGsize = argsVec[6],
	    m = argsVec[groupID + 7],
	    n = txtsLenVec[localID],
	    doffset = 1,gapmismax;

	unsigned int j_min, j_max, abs_ij, patChar, txtChar;

	float temp_score, score = -1000000.0;		
	
	unsigned int min_gap = 0;
        unsigned int where;

	for( i = 0; i < n + 1 ; i++ )      H[i][0] = i;
	for( j = 0; j < m + 1 ; j++ )      H[0][j] = j;

	if(GM_max ( 1,  m - max_gap )==1)
	{
		score = ( m - 1 ) * pensVec[1] + pensVec[0];		 
               	min_gap   = m;
               	where     = 1;		
              	(*start)       = 0;	
	}	

	for( i = 1; i <= m; i++ )
	{
		patChar = patsVec [groupID * maxPatLen + i - 1];

		j_min = GM_max ( 1,  i - max_gap );

		j_max = GM_min ( n,  i + max_gap );

		if(i<= max_gap+1)
		{
			cur_diag_nxt = 0;
		}
		else
		{
			cur_diag_nxt = hproVec[hproVecGsize*groupID + doffset*blockSize + localID];		
			doffset++;		
		}

		for( j = j_min; j <= j_max; j++)
		{
			txtChar = txtsVec[(j-1)*blockSize + localID]; 

			mis = cur_diag_nxt + EDNAFULL_matrix [txtChar][patChar];
			
			if(j<i)
			{	
				abs_ij = i-j;
				gap = dproVec[dproVecGsize*groupID + j*blockSize + localID];
			}
			else
			{	
				abs_ij = j-i;
				gap = dproVec[dproVecGsize*groupID + i*blockSize + localID]; 
			}

			if(i==j)
			{
				gap = mis - 1;
				dproVec[dproVecGsize*groupID + j*blockSize + localID]= mis;
			}
	
			if( gap > mis )	H[j][i] = abs_ij;

			cur_diag_nxt = hproVec[hproVecGsize*groupID + j*blockSize + localID];
	
			gapmismax = GM_max ( mis, gap );
		
			hproVec[hproVecGsize*groupID + j*blockSize + localID] = gapmismax;

			if(i==m)
			{				
				if(abs_ij<=max_gap)
				{
					temp_score = (float)gapmismax;			
				
					if(abs_ij>0)
						temp_score += ( abs_ij - 1 ) * pensVec[1] + pensVec[0];			

					if(temp_score>score)
					{
						score = temp_score;

						if(i>j)
						{							
               						 min_gap   = m-j;
               						 where     = 1;		
              						 (*start)  = j;		
						}

						if(i==j)
						{
							
            						min_gap   = 0;
            						where     = 0;		
            						(*start)  = m;		
						}

						if(i<j)
						{						  
						       min_gap   = j-m;
						       where     = 2;		
						       ( *start )   = j;		
						}								
					}
				}
			}								
		}
	}

	scrsVec[groupID*blockSize+localID]=score;
	out->max_score = score;
	out->min_gap = min_gap;
	out->where = where;	
}

/* runs on the CPU */
void kernel_one_to_one_pro (unsigned int groupID, unsigned int localID, unsigned int * patsVec, unsigned int * txtsVec, int * argsVec, int * txtsLenVec, float * pensVec, int * hproVec, int * dproVec, float * scrsVec, unsigned int ** H, struct gapmis_align * out, int * start)
{
	int i, j, cur_diag_nxt, mis, gap,
	    max_gap = argsVec[2], 
	    blockSize = argsVec[3],
	    maxPatLen = argsVec[4],
	    dproVecGsize = argsVec[5],
	    hproVecGsize = argsVec[6],
	    m = argsVec[groupID + 7],
	    n = txtsLenVec[localID],
	    doffset = 1,gapmismax;

	unsigned int j_min, j_max, abs_ij, patChar, txtChar;

	float temp_score, score = -1000000.0;		
	
	unsigned int min_gap = 0;
        unsigned int where;

	for( i = 0; i < n + 1 ; i++ )      H[i][0] = i;
	for( j = 0; j < m + 1 ; j++ )      H[0][j] = j;

	if(GM_max ( 1,  m - max_gap )==1)
	{
		score = ( m - 1 ) * pensVec[1] + pensVec[0];		 
               	min_gap   = m;
               	where     = 1;		
              	(*start)       = 0;	
	}	

	for( i = 1; i <= m; i++ )
	{
		patChar = patsVec [groupID * maxPatLen + i - 1];

		j_min = GM_max ( 1,  i - max_gap );

		j_max = GM_min ( n,  i + max_gap );

		if(i<= max_gap+1)
		{
			cur_diag_nxt = 0;
		}
		else
		{
			cur_diag_nxt = hproVec[hproVecGsize*groupID + doffset*blockSize + localID];		
			doffset++;		
		}

		for( j = j_min; j <= j_max; j++)
		{
			txtChar = txtsVec[(j-1)*blockSize + localID]; 

			mis = cur_diag_nxt + EBLOSUM62_matrix [txtChar][patChar];
			
			if(j<i)
			{	
				abs_ij = i-j;
				gap = dproVec[dproVecGsize*groupID + j*blockSize + localID];
			}
			else
			{	
				abs_ij = j-i;
				gap = dproVec[dproVecGsize*groupID + i*blockSize + localID]; 
			}

			if(i==j)
			{
				gap = mis - 1;
				dproVec[dproVecGsize*groupID + j*blockSize + localID]= mis;
			}
	
			if( gap > mis )	H[j][i] = abs_ij;

			cur_diag_nxt = hproVec[hproVecGsize*groupID + j*blockSize + localID];
	
			gapmismax = GM_max ( mis, gap );
		
			hproVec[hproVecGsize*groupID + j*blockSize + localID] = gapmismax;

			if(i==m)
			{				
				if(abs_ij<=max_gap)
				{
					temp_score = (float)gapmismax;				
		
				
					if(abs_ij>0)
						temp_score += ( abs_ij - 1 ) * pensVec[1] + pensVec[0];			

					if(temp_score>score)
					{
						score = temp_score;

						if(i>j)
						{							
               						 min_gap   = m-j;
               						 where     = 1;		
              						 (*start)  = j;		
						}

						if(i==j)
						{							
            						min_gap   = 0;
            						where     = 0;		
            						(*start)  = m;		
						}

						if(i<j)
						{						  
						       min_gap   = j-m;
						       where     = 2;		
						       ( *start ) = j;		
						}								
					}
				}
			}								
		}
	}

	scrsVec[groupID*blockSize+localID]=score;
	out->max_score = score;
	out->min_gap = min_gap;
	out->where = where;
}

unsigned int kernel_launch (cl_kernel kernel, cl_context context, cl_command_queue cmd_queue, const char ** p, const char ** t, const struct gapmis_params * in, float * scores)
{
        int error=1;
	unsigned int	pats = get_number_of_sequences (p);
	unsigned int	txts = get_number_of_sequences (t);

	unsigned int	maxTxtLen = get_max_length (txts, t);
	unsigned int	maxPatLen = get_max_length (pats, p);
	unsigned int	pBlockSize = get_pblock_size (txts,32); 
	unsigned int 	hproVecLen = pats * pBlockSize * (maxTxtLen + 1);
	unsigned int 	dproVecLen = pats * pBlockSize * (maxPatLen + 1);
	
	unsigned int * txtsVec = calloc(maxTxtLen*pBlockSize, sizeof(unsigned int));
	unsigned int * patsVec = calloc(maxPatLen*pats, sizeof(unsigned int));
	         int * argsVec = malloc(sizeof(int)*(pats+7));		
		 int * txtsLenVec = calloc(pBlockSize,sizeof(int));
	       float * pensVec = malloc(sizeof(float)*2);
	         int * hproVec = calloc(hproVecLen,sizeof(int));
	         int * dproVec = calloc(dproVecLen,sizeof(int));

	cl_int err;	

	if(patsVec==NULL   || txtsVec==NULL      || argsVec==NULL || 
           pensVec == NULL || txtsLenVec == NULL || hproVec==NULL || 
	   dproVec==NULL)
	{	
	 	errno = MALLOC;
      		return ( 0 );
	}


	fill_txtsVec (txts, pBlockSize, t, txtsVec, in->scoring_matrix);
	fill_patsVec (pats, maxPatLen, p, patsVec, in->scoring_matrix);
	fill_argsVec (pats, txts, p, in->max_gap, pBlockSize, maxPatLen, maxTxtLen, argsVec);
	fill_txtsLenVec (txts, t, txtsLenVec);

	pensVec[0] = - in -> gap_open_pen;
	pensVec[1] = - in -> gap_extend_pen;	
	
	/* GPU malloc */
	cl_mem txtsVec_device = malloc_device (context, (maxTxtLen*pBlockSize)*sizeof(unsigned int), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}
	
	/* copy from CPU to GPU mem */
	init_device_mem_uint (context, cmd_queue, txtsVec_device, txtsVec, maxTxtLen*pBlockSize, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}
	
	cl_mem patsVec_device = malloc_device (context, (maxPatLen*pats)*sizeof(unsigned int), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}

	init_device_mem_uint (context, cmd_queue, patsVec_device, patsVec,maxPatLen*pats, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}

	cl_mem argsVec_device = malloc_device (context, (pats+7)*sizeof(int), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}

	init_device_mem_int (context, cmd_queue, argsVec_device, argsVec, pats+7, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}

	cl_mem txtsLenVec_device = malloc_device (context, pBlockSize*sizeof(int), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}

	init_device_mem_int (context, cmd_queue, txtsLenVec_device, txtsLenVec, pBlockSize, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}

	cl_mem pensVec_device = malloc_device (context, 2*sizeof(float), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}

	init_device_mem_float (context, cmd_queue, pensVec_device, pensVec, 2, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}

	cl_mem hproVec_device = malloc_device (context, hproVecLen*sizeof(int), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}

	init_device_mem_int (context, cmd_queue, hproVec_device, hproVec, hproVecLen, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}

	cl_mem dproVec_device = malloc_device (context, dproVecLen*sizeof(int), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}

	init_device_mem_int (context, cmd_queue, dproVec_device, dproVec, dproVecLen, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}

	cl_mem scrsVec_device = malloc_device (context, (pats*pBlockSize)*sizeof(float), &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}

	err = clFinish(cmd_queue);
	if(err != CL_SUCCESS)
	{
	 	errno = GPUERROR;
      		return ( 0 );	
	}

	/* connect the input arguments of the kernel with the corresponding mem */
	set_kernel_arguments (kernel, cmd_queue, patsVec_device, txtsVec_device, argsVec_device, txtsLenVec_device, pensVec_device, hproVec_device, dproVec_device, scrsVec_device);

	/* synchronisation */
	err = clFinish(cmd_queue);
	if(err != CL_SUCCESS)
	{
	 	errno = GPUERROR;
      		return ( 0 );	
	}	


	/* WorkSizeGlobal is the total number of threads of the device*/
	size_t WorkSizeGlobal[] = {pBlockSize * pats};
	/* WorkSizeLocal is the number of threads per group*/
	size_t WorkSizeLocal[] = {pBlockSize};

	/* kernel enters the command queue using WorkSizeGlobal and WorkSizeLocal */
	err = clEnqueueNDRangeKernel(cmd_queue, kernel, 1, NULL, WorkSizeGlobal, WorkSizeLocal, 0, NULL, NULL);
	if(error)
	{
	 	errno = KERNEL;
      		return ( 0 );	
	}

	/* finalise the kernel */
	err = clFinish(cmd_queue);
	if(err != CL_SUCCESS)
	{
	 	errno = GPUERROR;
      		return ( 0 );	
	}	

	/* return the results from the GPU to the CPU */
	read_device_mem_float (cmd_queue, pats*pBlockSize, scores, scrsVec_device, &error);
	if(error)
	{
	 	errno = GPUMALLOC;
      		return ( 0 );	
	}

	/* deallocation */
	free (txtsVec);
	free (patsVec);
	free (argsVec);
	free (txtsLenVec);
	free (pensVec);
	free (hproVec);
	free (dproVec);

	clReleaseMemObject(patsVec_device);
	clReleaseMemObject(txtsVec_device);
	clReleaseMemObject(argsVec_device);
	clReleaseMemObject(txtsLenVec_device);
	clReleaseMemObject(pensVec_device);
	clReleaseMemObject(hproVec_device);
	clReleaseMemObject(dproVec_device);
	clReleaseMemObject(scrsVec_device);

	return ( 1 );
}

unsigned int kernel_launch_l (cl_kernel kernel, cl_context context, cl_command_queue cmd_queue, const char ** p, const char ** t, const struct gapmis_params * in, float * scores, struct gapmis_align * out)
{
	unsigned int 	i;
		 int	start;
	unsigned int	pats = 1;
	unsigned int	txts = 1;

	unsigned int	maxTxtLen = strlen(t[0]);
	unsigned int	maxPatLen = strlen(p[0]);
	unsigned int	pBlockSize = 1;
	unsigned int 	hproVecLen = pats * pBlockSize * (maxTxtLen + 1);
	unsigned int 	dproVecLen = pats * pBlockSize * (maxPatLen + 1);
	
	unsigned int * txtsVec = calloc(maxTxtLen*pBlockSize, sizeof(unsigned int));
	unsigned int * patsVec = calloc(maxPatLen*pats, sizeof(unsigned int));
	         int * argsVec = malloc(sizeof(int)*(pats+7));		
		 int * txtsLenVec = calloc(pBlockSize,sizeof(int));
	       float * pensVec = malloc(sizeof(float)*2);
	         int * hproVec = calloc(hproVecLen,sizeof(int));
	         int * dproVec = calloc(dproVecLen,sizeof(int));
       unsigned int ** H;

	H = malloc((maxTxtLen+1)*sizeof(unsigned int*));
	H[0] = calloc ((maxTxtLen+1)*(maxPatLen+1), sizeof(unsigned int));
	for ( i = 1 ; i < maxTxtLen + 1 ; ++ i )
     		H[i] = (void*)H[0] + i*(maxPatLen+1)* sizeof(unsigned int);


	if(patsVec == NULL || txtsVec == NULL    || argsVec == NULL || 
           pensVec == NULL || txtsLenVec == NULL || hproVec == NULL || 
	   dproVec == NULL || H == NULL)
	{
	 	errno = MALLOC;
      		return ( 0 );	
	}

	fill_txtsVec (txts, pBlockSize, t, txtsVec,in->scoring_matrix);
	fill_patsVec (pats, maxPatLen, p, patsVec,in->scoring_matrix);
	fill_argsVec (pats, txts, p, in->max_gap, pBlockSize, maxPatLen, maxTxtLen, argsVec);
	fill_txtsLenVec (txts, t, txtsLenVec);

	pensVec[0] = - in -> gap_open_pen;
	pensVec[1] = - in -> gap_extend_pen;


	if(in->scoring_matrix==0)
		kernel_one_to_one_dna (0, 0, patsVec, txtsVec, argsVec, txtsLenVec, pensVec, hproVec, dproVec, scores, H, out, &start);
	else
		kernel_one_to_one_pro (0, 0, patsVec, txtsVec, argsVec, txtsLenVec, pensVec, hproVec, dproVec, scores, H, out, &start);  

   	if ( out -> min_gap > 0 ) 
        	backtracing ( H, maxPatLen, maxTxtLen, start, out );


   	if ( out -> where == 1 ) 
     		num_mismatch ( t[0], maxTxtLen, p[0], maxPatLen, out );
   	else                   
     		num_mismatch ( p[0], maxPatLen, t[0], maxTxtLen, out );

	free (txtsVec);
	free (patsVec);
	free (argsVec);
	free (txtsLenVec);
	free (pensVec);
	free (hproVec);
	free (dproVec);
        free ( H[0] );
	free (H);
	
	return 1;
}

cl_platform_id get_gpu_id(int * error)
{
	cl_platform_id gpu_id = NULL;

	cl_uint platforms_total=0;
	
	cl_int err;

	err = clGetPlatformIDs (0, NULL, &platforms_total);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}
	if(platforms_total<=0)
	{
		*error = 1;
		return NULL;
	}
	cl_platform_id * gpu_id_vec = malloc(sizeof(cl_platform_id) * platforms_total);

	err = clGetPlatformIDs (platforms_total, gpu_id_vec, NULL);

   // printf( "platforms: %d\n", platforms_total );
    int i;
    int use_platform = -1;
    
    
    // choose the nvidia platform
    for( i = 0; i < platforms_total; ++i ) {
        char str[256];
        
        clGetPlatformInfo( gpu_id_vec[i], CL_PLATFORM_VENDOR, sizeof(str), str, NULL );
        
     //   printf( "vendor: %s\n", str );
        
        if( strstr( str, "NVIDIA" ) != NULL ) {
            use_platform = i;
            break;
        }
    }
    
	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return NULL;
	}
	
	if( use_platform == -1 ) {
        *error = 1;
        return NULL;
    }
	gpu_id = gpu_id_vec[use_platform];

	free( gpu_id_vec );
	*error = 0;

	return gpu_id;
}

cl_device_id get_dev_id(cl_platform_id  gpu_id, int * error)
{
	cl_device_id dev_id = NULL;

	cl_uint devices_total;

	cl_int err;
	
	err = clGetDeviceIDs(gpu_id, CL_DEVICE_TYPE_GPU, 0, NULL, &devices_total);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	if(devices_total<=0)
	{
		*error = 1;
		return NULL;
	}

	cl_device_id * dev_id_vec = malloc(sizeof(cl_device_id) * devices_total);

	err = clGetDeviceIDs(gpu_id, CL_DEVICE_TYPE_ALL, devices_total, dev_id_vec, NULL);

    
    
    
	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	dev_id = dev_id_vec[0];

	free( dev_id_vec );
	*error = 0;
	return dev_id;
}

cl_context create_context(cl_device_id dev_id, int * error)
{


	cl_context context;

	cl_int err;
    
	context = clCreateContext (0,1, &dev_id, NULL,NULL, &err);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}
	
	*error = 0;
	return context;
}

cl_command_queue create_cmd_queue (cl_device_id dev_id, cl_context context, int * error)
{
	cl_int err;

	cl_command_queue cmd_queue;

	cmd_queue = clCreateCommandQueue(context, dev_id, 0, &err);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	return cmd_queue;
}

cl_kernel load_kernel (char * name, char * kernel_name, cl_device_id dev_id, cl_context context, int * error)
{
	cl_kernel kernel;

	FILE * fp = fopen(name, "r");

	if(fp==NULL)
	{
		*error = 1;
		return NULL;
	}

	char * source;
	size_t size;

	fseek(fp, 0, SEEK_END);
	size = ftell(fp);
	fseek(fp, 0, SEEK_SET);
		
	source = malloc((size+1)*sizeof(char));
	
	size = fread(source, 1, size, fp);

	fclose(fp);
	
	source[size] = '\0';

	cl_int err;

	cl_program program = clCreateProgramWithSource (context, 1, (const char **) &source, &size, &err);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	clBuildProgram(program, 0, NULL, NULL, NULL, NULL);

	cl_build_status status;

	clGetProgramBuildInfo(program, dev_id, CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status, NULL);

	if(status!=CL_BUILD_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	kernel = clCreateKernel(program, kernel_name , &err);

	if(err!=CL_SUCCESS)
	{
		*error = 1;
		return NULL;
	}

	*error = 0;
        free ( source );
	clReleaseProgram(program);
	return kernel;
}

cl_mem malloc_device (cl_context context, size_t size, int * error)
{
	cl_mem mem = NULL;

	cl_int err;

	mem = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &err);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return NULL;
	}

	*error=0;
	return mem;	
}

void init_device_mem_int (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, int * mem, size_t size, int * error)
{
	cl_int err;

	err = clEnqueueWriteBuffer(cmd_queue, dev_mem, CL_FALSE, 0, size * sizeof(int), mem, 0, NULL, NULL);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	err = clFinish(cmd_queue);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	*error = 0;
	return;
}

void init_device_mem_uint (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, unsigned int * mem, size_t size, int * error)
{
	cl_int err;

	err = clEnqueueWriteBuffer(cmd_queue, dev_mem, CL_FALSE, 0, size * sizeof(unsigned int), mem, 0, NULL, NULL);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	err = clFinish(cmd_queue);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	*error = 0;
	return;
}

void init_device_mem_float (cl_context context, cl_command_queue cmd_queue, cl_mem dev_mem, float * mem, size_t size, int * error)
{
	cl_int err;

	err = clEnqueueWriteBuffer(cmd_queue, dev_mem, CL_FALSE, 0, size * sizeof(float), mem, 0, NULL, NULL);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	err = clFinish(cmd_queue);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	*error = 0;
	return;
}

void set_kernel_arguments (cl_kernel kernel, cl_command_queue cmd_queue, cl_mem cl_mem0, cl_mem cl_mem1, cl_mem cl_mem2, cl_mem cl_mem3, cl_mem cl_mem4, cl_mem cl_mem5, cl_mem cl_mem6, cl_mem cl_mem7)
{
	/* connect each allocated vector to the corresponding argument */
	clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &cl_mem0);
	clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &cl_mem1);
	clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &cl_mem2);
	clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &cl_mem3);
	clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &cl_mem4);
	clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &cl_mem5);
	clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &cl_mem6);
	clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &cl_mem7);
	clFinish(cmd_queue);
}

void read_device_mem_float (cl_command_queue cmd_queue, size_t size, float * mem, cl_mem dev_mem, int * error)
{
	cl_int err;

	err = clEnqueueReadBuffer(cmd_queue, dev_mem, CL_FALSE, 0, size * sizeof(float), mem, 0, NULL, NULL);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}

	err = clFinish(cmd_queue);

	if(err!=CL_SUCCESS)
	{	
		*error = 1;
		return;
	}
	
	*error = 0;
	return;
}

void update_group_match (float * groupScores, int * groupMatch, float * groupMatchScores, unsigned int patGroupSize, unsigned int txtGroupSize, int pats, int txts, int patGroupIndex, int txtGroupIndex)
{
	int i,j;

	float max_score;
	int position=-1;

	for(i=0;i<patGroupSize;i++)
	{
		if(patGroupIndex*patGroupSize+i<pats)
		{		
			max_score=groupMatchScores[i];
			position=groupMatch[i];
	
			for(j=0;j<txtGroupSize;j++)
			{		
				if( groupScores[i*txtGroupSize+j] > max_score && txtGroupIndex*txtGroupSize + j < txts)
				{
					max_score = groupScores[i*txtGroupSize+j];
					position = txtGroupIndex*txtGroupSize+j;	
				}
			}

			groupMatch[i]=position;
			groupMatchScores[i]=max_score;
		}
	}
}

void set_invalid(int * groupMatch, int groupSize)
{
	int i;
	for(i=0;i<groupSize;i++)
		groupMatch[i]=-1;
}

void set_minimum(float * groupMatchScores, int groupSize)
{
	int i;
	for(i=0;i<groupSize;i++)
		groupMatchScores[i]=-DBL_MAX;
}

unsigned int get_number_of_groups (int elements, int groupSize)
{	
	unsigned int div32 = elements / groupSize;
	unsigned int mod32 = elements % groupSize;
	
	unsigned int groups = mod32!=0?(div32+1):div32;	
	
	return groups;
}

void set_null (const char ** input, int size)
{
	int i;
	for(i=0;i<size;i++)
		input[i]=NULL;
}

void initialize_pointers (const char * groupPatterns[], int groupIndex, int groupSize, const char ** source, int sourceSize)
{
	int i, elements = sourceSize - groupIndex * groupSize;

	if(elements>=groupSize)
	{
		for(i=0;i<groupSize;i++)
			groupPatterns[i]=source[i + groupIndex * groupSize];
	}
	else
	{
		for(i=0;i<elements;i++)
			groupPatterns[i]=source[i + groupIndex * groupSize];

		for(i=elements;i<groupSize;i++)
			groupPatterns[i]=NULL;
	}
}

unsigned int get_number_of_sequences (const char ** input)
{
	unsigned int 	total = 0;
	const char 	** Tmp;

	for ( Tmp = input; *Tmp; ++ Tmp, ++ total );

	return total;  
}

unsigned int get_max_length (unsigned int total, const char ** input)
{
	unsigned int	i, curLen, maxLen=0;
   
	for (i=0;i<total;i++)
	{
		curLen = strlen(input[i]);
		if(curLen>maxLen)
			maxLen = curLen;
	} 

	return maxLen;
}

unsigned int get_min_length (unsigned int total, const char ** input)
{
	unsigned int	i, curLen, minLen=4294967295u;
   
	for (i=0;i<total;i++)
	{
		curLen = strlen(input[i]);
		if(curLen<minLen)
			minLen = curLen;
	} 

	return minLen;
}

unsigned int get_pblock_size (unsigned int input, unsigned int mult)
{
	unsigned int div32 = input / mult;
	unsigned int mod32 = input % mult;
	
	unsigned int result = mod32!=0?(div32+1)*mult:div32*mult;	
	
	return result;
}

void fill_txtsVec (unsigned int total, unsigned int blockSize, const char ** in, unsigned int * outVec, int matrix)
{
	unsigned int i,j,len;
	
	if(matrix==0)
		for(i=0;i<total;i++)
		{		
			len = strlen(in[i]);
	
			for(j=0;j<len;j++)
				outVec[j*blockSize + i] = nuc_char_to_index ( in[i][j]);
		}
	else
		for(i=0;i<total;i++)
		{		
			len = strlen(in[i]);
	
			for(j=0;j<len;j++)
				outVec[j*blockSize + i] = pro_char_to_index ( in[i][j]);
		}
}

void fill_patsVec (unsigned int total, unsigned int blockSize, const char ** in, unsigned int * outVec, int matrix)
{
	unsigned int i,j,len,spos;

	if(matrix==0)
		for(i=0;i<total;i++)
		{		
			len = strlen(in[i]);
			spos = i * blockSize;

			for(j=0;j<len;j++)
				outVec[spos + j] = nuc_char_to_index ( in[i][j]);
		}
	else
		for(i=0;i<total;i++)
		{		
			len = strlen(in[i]);
			spos = i * blockSize;

			for(j=0;j<len;j++)
				outVec[spos + j] = pro_char_to_index ( in[i][j]);
		}	
}

void fill_argsVec (unsigned int totalPats, unsigned int totalTxts, const char ** in, unsigned int max_gap, unsigned int pBlockSize, unsigned int maxPatLen, unsigned int maxTxtLen, int * outVec)
{
	unsigned int i;

	outVec[0] = totalPats;
	outVec[1] = totalTxts;
	outVec[2] = max_gap;
	outVec[3] = pBlockSize;
	outVec[4] = maxPatLen;
	outVec[5] = pBlockSize * (maxPatLen + 1);
	outVec[6] = pBlockSize * (maxTxtLen + 1);

	for(i=7;i<totalPats+7;i++)		
		outVec[i] = strlen(in[i-7]);
}

void fill_txtsLenVec (unsigned int totalTxts, const char ** in, int * outVec)
{
	unsigned int i;

	for(i=0;i<totalTxts;i++)		
		outVec[i] = strlen(in[i]);
}

int check_sequences(unsigned int total, const char ** sequences, unsigned int matrix)
{
	unsigned int i,j,length, index;

	if(matrix==0)
		for(i=0;i<total;i++)
		{
			length = strlen(sequences[i]);
		
			for(j=0;j<length;j++)
			{
				index = nuc_char_to_index (sequences[i][j]);

				if(index >= GM_NUC_SCORING_MATRIX_SIZE)
					return 0;			
			}
		}
	else
		for(i=0;i<total;i++)
		{
			length = strlen(sequences[i]);
		
			for(j=0;j<length;j++)
			{
				index = pro_char_to_index (sequences[i][j]);

				if(index >= GM_PRO_SCORING_MATRIX_SIZE)
					return 0;
			}
		}
	return 1;
}
#endif
