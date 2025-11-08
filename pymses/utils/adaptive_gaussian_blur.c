/******************************************************************************************
*   This file is part of PyMSES.                                                          *
*                                                                                         *
*   PyMSES is free software: you can redistribute it and/or modify                        *
*   it under the terms of the GNU General Public License as published by                  *
*   the Free Software Foundation, either version 3 of the License, or                     *
*   (at your option) any later version.                                                   *
*                                                                                         *
*   PyMSES is distributed in the hope that it will be useful,                             *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of                        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                         *
*   GNU General Public License for more details.                                          *
*                                                                                         *
*   You should have received a copy of the GNU General Public License                     *
*   along with PyMSES.  If not, see <http://www.gnu.org/licenses/>.                       *
******************************************************************************************/
//include <stdio.h>
//include <math.h>
//include "adaptive_gaussian_blur.h"

static void compute_same_value_pixel_size_map_C(int * result_map, double * map, int map_max_size_i, int map_max_size_j){
	printf("Computing same value pixel size map...\n");	
	double val, map_min, map_max, max_val_diff, equal_test_precision;
	int i, j, ii, jj, same_value_pixel_size;
	int map_size = map_max_size_i * map_max_size_j;
	int set_size_limit = 300;
	char test;
	
	// find map min and max to get an adapted equal_test_precision value:
	map_min = map[0];
	map_max = map[0];
	for (i = 0 ; i < map_size ; i++){
		if (map[i] > map_max)
			map_max = map[i];
		if (map[i] < map_min)
			map_min = map[i];
	}
	max_val_diff = map_max - map_min;
	equal_test_precision = 0.000001 * max_val_diff;
	
	// loop on map value
	for (i = 0 ; i < map_max_size_i ; i++){
		for (j = 0 ; j < map_max_size_j ; j++){
			if (result_map[i*map_max_size_j + j] == 0){
				val = map[i * map_max_size_j + j];
				
				// first loop to get the same_value_pixel_size
				same_value_pixel_size = 1;
				ii = i;
				test = 1;
				// loop j++
				jj = j;
				while(test && jj < map_max_size_j && jj - j < set_size_limit){
					jj++;						
					test = fabs(map[ii * map_max_size_j + jj] - val)<=equal_test_precision;
				}
				if (jj - j > same_value_pixel_size)
					same_value_pixel_size = jj - j;
					
				result_map[i * map_max_size_j + j] = same_value_pixel_size;
				// second exact same loop to report result inside result_map
				test = 1;
				// loop j++
				jj = j;
				while(test && jj < map_max_size_j && jj - j < set_size_limit){
					jj++;						
					test = fabs(map[ii * map_max_size_j + jj] - val)<=equal_test_precision;
					if (test)
						result_map[ii * map_max_size_j + jj] = same_value_pixel_size;
				}
			}
		}
	}

	// Last loop on map value : extend gauss filter size area by one pixel to reduce the difference where there is a change of pixel size
	int * save_result_map = malloc(map_size * sizeof(int));
	for (i = 0 ; i < map_size ; i++)
		save_result_map[i] = result_map[i];
	char higher_value_neighbor;
	for (i = 0 ; i < map_max_size_i ; i++){
		for (j = 0 ; j < map_max_size_j ; j++){
			val = save_result_map[i*map_max_size_j + j];
			// this following code does :
			// if one neighbor pixel in save_result_map is higher
			// then we add 1
			higher_value_neighbor = 0;
			if (i > 0){
				if (j >0 && save_result_map[(i-1)*map_max_size_j + j-1] > val)
					higher_value_neighbor = 1;
				if (save_result_map[(i-1)*map_max_size_j + j] > val)
					higher_value_neighbor = 1;
				if (j+1 < map_max_size_j && save_result_map[(i-1)*map_max_size_j + j+1] > val)
					higher_value_neighbor = 1;
			}
			if (j >0 && save_result_map[i*map_max_size_j + j-1] > val)
				higher_value_neighbor = 1;
			if (j+1 < map_max_size_j && save_result_map[i*map_max_size_j + j+1] > val)
				higher_value_neighbor = 1;
			if (i+1 < map_max_size_i){
				if (j >0 && save_result_map[(i+1)*map_max_size_j + j-1] > val)
					higher_value_neighbor = 1;
				if (save_result_map[(i+1)*map_max_size_j + j] > val)
					higher_value_neighbor = 1;
				if (j+1 < map_max_size_j && save_result_map[(i+1)*map_max_size_j + j+1] > val)
					higher_value_neighbor = 1;
			}
			if (higher_value_neighbor)
				result_map[i*map_max_size_j + j] += 1;
		}
	}
	free(save_result_map);
}

static void adaptive_gaussian_blur_C(double * map_filtered, double * original_map, int * same_value_pixel_size_map, int map_max_size_i, int map_max_size_j){
	double distance_to_mask_center_square, gaus_val, mask_val, norm_sum;
	printf("Computing adaptive gaussian blur...\n");
	int max_size = 0;
	int mask_size, size, two_sigma_square, max_mask_size_array, index, mask_size_index;
	int different_size_number = 0;
	int i, j, k, ii, jj;
	char new_value;
	//for (j =0;j<16;j++){
	//	printf(" %i ",j);
	//	printf("same_value_pixel_size_map = %i \n",same_value_pixel_size_map[j]);
	//	printf("map_filtered = %f \n",map_filtered[j]);
	//	printf("original_map = %f \n",original_map[j]);
	//}
	int * different_size = malloc(map_max_size_i * sizeof(int));
	// create mask dict : find every different size mask needed 
	for (i = 0 ; i < map_max_size_i ; i++){
		for (j = 0 ; j < map_max_size_j ; j++){
			size = same_value_pixel_size_map[i * map_max_size_j + j];
			// try to see if it is a new size value
			if (size > 0) {
				new_value = 1;
				for (k = 0 ; k < different_size_number ; k++){
					if (size == different_size[k])
						new_value = 0;
				}
				if (new_value) {
					different_size[different_size_number] = size;
					different_size_number++;
					if (size > max_size)
						max_size = size;
				}
			}
		}
	}
	
	mask_size = max_size + 1;
	max_mask_size_array = mask_size * mask_size;
	double* mask_size_array = malloc(different_size_number * max_mask_size_array * sizeof(double));
	// create mask dict : compute gaussian values for every different size mask 
	for (k = 0 ; k < different_size_number ; k++){
		size = different_size[k];
		mask_size = size + 1;
		norm_sum = 0;
		// we try here with sigma = size but it is probably too big :
		// what about double all mask size and reduce sigma by 2?
		two_sigma_square = 2 * size*size;
		// Use symmetry to reduce thoose loops size by 4 here
		for (i = 0 ; i < mask_size ; i++){
			for (j = 0 ; j < mask_size ; j++){
				distance_to_mask_center_square = (double)(i*i + j*j);
				gaus_val = exp((double)(-distance_to_mask_center_square/two_sigma_square));
				mask_size_array[k * max_mask_size_array + i * mask_size + j] = gaus_val;
				if (i == 0){
					if (j == 0)
						norm_sum += gaus_val;
					else
						norm_sum += 2 * gaus_val;
				}
				else{
					if (j == 0)
						norm_sum += 2 * gaus_val;
					else
						norm_sum += 4 * gaus_val;
				}
			}
		}
		// Norm:
		for (i = 0 ; i < mask_size * mask_size ; i++){
			mask_size_array[k * max_mask_size_array + i] /= norm_sum;
			//if (size == 1){
			//	printf("\n mask_size_array[k * max_mask_size_array + i] : %f",mask_size_array[k * max_mask_size_array + i]);
			//}
		}
	}
	//printf(" different_size_number : %i ",different_size_number);
	// we compute the resulting map
	for (i = 0 ; i < map_max_size_i ; i++){
		for (j = 0 ; j < map_max_size_j ; j++){
			index = i * map_max_size_j + j;
			size = same_value_pixel_size_map[index];
			//if (i == 3 && j == 2){
			//	printf("\noriginal_map[index] : %f",original_map[index]);
			//	printf(" size : %i",size);
			//}
			if (size == 0)
				map_filtered[index] += original_map[index];
			else {
				// find mask size index :
				for (k = 0 ; k < different_size_number ; k++){
					// we should always find a value as it is previously computed
					if (different_size[k] == size)
						mask_size_index = k;
				}
				// eventually compute the gaussian filtered pixel value:
				// Like in mask computation, use symmetry to reduce thoose loops size by 4 here
				// middle mask value :
				map_filtered[index] += original_map[index] * mask_size_array[mask_size_index * max_mask_size_array];
				//if (i == 3 && j == 2){
				//	printf("\nmask_mid : %f ",mask_size_array[mask_size_index * max_mask_size_array]);
				//	printf(" * %f \n",original_map[index]);
				//}
				mask_size = size + 1;
				for (ii = 1 ; ii < mask_size ; ii++){
					// middle cross mask values :
					mask_val = mask_size_array[mask_size_index * max_mask_size_array + ii];
					//if (i == 3 && j == 2){
					//	printf("\nmask_mid : %f ",mask_val);
					//	printf(" * %f \n",original_map[(i-ii) * map_max_size_j + j]);
					//	printf(" * %f \n",original_map[(i+ii) * map_max_size_j + j]);
					//	printf(" * %f \n",original_map[i * map_max_size_j + j-ii]);
					//	printf(" * %f \n",original_map[i * map_max_size_j + j+ii]);
					//}
					if ((i-ii) >= 0)
						map_filtered[index] += original_map[(i-ii) * map_max_size_j + j] * mask_val;
					else
						map_filtered[index] += original_map[j] * mask_val;
					if ((i+ii) < map_max_size_i)
						map_filtered[index] += original_map[(i+ii) * map_max_size_j + j] * mask_val;
					else
						map_filtered[index] += original_map[(map_max_size_i - 1) * map_max_size_j + j] * mask_val;
					if ((j-ii) >= 0)
						map_filtered[index] += original_map[i * map_max_size_j + j - ii] * mask_val;
					else
						map_filtered[index] += original_map[i * map_max_size_j] * mask_val;
					if ((j+ii) < map_max_size_j)
						map_filtered[index] += original_map[i * map_max_size_j + j + ii] * mask_val;
					else
						map_filtered[index] += original_map[i * map_max_size_j + map_max_size_j - 1] * mask_val;
					// others values :
					for (jj = 1 ; jj < mask_size ; jj++){
						mask_val = mask_size_array[mask_size_index * max_mask_size_array + ii * mask_size + jj];
						//if (i == 3 && j == 2){
						//	printf("\nmask_val : %f ",mask_val);
						//	printf(" * %f \n",original_map[(i-ii) * map_max_size_j + j-jj]);
						//	printf(" * %f \n",original_map[(i-ii) * map_max_size_j + j+jj]);
						//	printf(" * %f \n",original_map[(i+ii) * map_max_size_j + j-jj]);
						//	printf(" * %f \n",original_map[(i+ii) * map_max_size_j + j+jj]);
						//}
						if ((i-ii) >= 0) {
							if ((j-jj) >= 0)
								map_filtered[index] += original_map[(i-ii) * map_max_size_j + j-jj] * mask_val;
							else
								map_filtered[index] += original_map[(i-ii) * map_max_size_j] * mask_val;
							if ((j+jj) < map_max_size_j)
								map_filtered[index] += original_map[(i-ii) * map_max_size_j + j+jj] * mask_val;
							else
								map_filtered[index] += original_map[(i-ii) * map_max_size_j + map_max_size_j - 1] * mask_val;
						}
						else {
							if ((j-jj) >= 0)
								map_filtered[index] += original_map[j-jj] * mask_val;
							else
								map_filtered[index] += original_map[0] * mask_val;
							if ((j+jj) < map_max_size_j)
								map_filtered[index] += original_map[j+jj] * mask_val;
							else
								map_filtered[index] += original_map[map_max_size_j - 1] * mask_val;
						}
						if ((i+ii) < map_max_size_i) {
							if ((j-jj) >= 0)
								map_filtered[index] += original_map[(i+ii) * map_max_size_j + j-jj] * mask_val;
							else
								map_filtered[index] += original_map[(i+ii) * map_max_size_j] * mask_val;
							if ((j+jj) < map_max_size_j)
								map_filtered[index] += original_map[(i+ii) * map_max_size_j + j+jj] * mask_val;
							else
								map_filtered[index] += original_map[(i+ii) * map_max_size_j + map_max_size_j - 1] * mask_val;
						}
						else {
							if ((j-jj) >= 0)
								map_filtered[index] += original_map[(map_max_size_i - 1) * map_max_size_j + j-jj] * mask_val;
							else
								map_filtered[index] += original_map[(map_max_size_i - 1) * map_max_size_j] * mask_val;
							if ((j+jj) < map_max_size_j)
								map_filtered[index] += original_map[(map_max_size_i - 1) * map_max_size_j + j+jj] * mask_val;
							else
								map_filtered[index] += original_map[(map_max_size_i - 1) * map_max_size_j + map_max_size_j - 1] * mask_val;
						}
					}
				}
			}
		}
	}
	free(mask_size_array);
	free(different_size);
}

	
	
	
	
	
	
	
	
	
	
	