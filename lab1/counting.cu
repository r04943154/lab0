#include "counting.h"
#include <cstdio>
#include <cassert>
#include <thrust/scan.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/inner_product.h>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

using namespace std;


__device__ __host__ int CeilDiv(int a, int b) { return (a-1)/b + 1; }
__device__ __host__ int CeilAlign(int a, int b) { return CeilDiv(a, b) * b; }

// determines whether the character is alphabetical
__host__ __device__
bool is_alpha(const char c)
{
	return (c >= 'a' && c <= 'z');
}

// determines whether the right character begins a new word
struct is_word_start
	: public thrust::binary_function<const char&, const char&, bool>
{
	__host__ __device__
		bool operator()(const char& left, const char& right) const
	{
		return is_alpha(right) && !is_alpha(left);
	}
};
int word_count(const thrust::device_vector<int>& input)
{
	// check for empty string
	if (input.empty())
		return 0;

	// compute the number characters that start a new word
	int wc = thrust::inner_product(input.begin(), input.end() - 1,  // sequence of left characters
		input.begin() + 1,               // sequence of right characters
		0,                               // initialize sum to 0
		thrust::plus<int>(),             // sum values together
		is_word_start());       // how to compare the left and right characters

	// if the first character is alphabetical, then it also begins a word
	if (is_alpha(input.front()))
		wc++;

	return wc;
}

__global__ void toOne(const char* A, int* B, int C)
{
	int i =  blockIdx.x + threadIdx.x;
	if (A[i] != '\n'){
		//printf(".@%d.",i);
		B[i] = 1;
	}
	else
	{
		B[i] = 0; 
		//printf(".~%d.", i);
	}
	/*int j = 1;
	for (int hi = 0; hi < C; hi++){
		if (A[hi] != '\n'){
			B[hi] = j;
			j++;
		}
		else{
			B[hi] = 0;
			j = 1;
		}
	}*/
}
__global__ void STadd(int A, int* B, int C)
{
	int i = blockIdx.x + threadIdx.x;
	if (B[i*2]==1&&B[i*2+1]==1)
	{
		B[C + i] = 2;
		if (i > 0){
			if (B[i * 2 - 2] == 0 && B[i * 2 - 1] == A && B[i * 2 + 2] == 0 && B[i * 2 + 3] == 0)B[C + i] = 2+A;
			if (B[i * 2 - 2] == 0 && B[i * 2 - 1] == A && B[i * 2 + 2] == A && B[i * 2 + 3] == 0)B[C + i] = 2+A+A;
			if (B[i * 2 - 2] == 0 && B[i * 2 - 1] == 0 && B[i * 2 + 2] == A && B[i * 2 + 3] == 0)B[C + i] = 2+A;
		}
	}
}

void CountPosition(const char *text, int *pos, int text_size)
{
	cout << endl << pos << endl;

	int numofblock = text_size / 256 + 1;
	toOne << <1, 1 >> >(text, pos, text_size);
	for (int i = 1; i < 11; i++){
		cudaDeviceSynchronize();
		STadd << <numofblock/pow(2,i), 256 >> >( i, pos, text_size);
	}
	int j = 1;


	thrust::device_ptr<const char> text_d(text);
	thrust::device_vector<char> input(text_d, text_d + 1000);

	// count words
	//int wc = word_count(input);

	//std::cout << "Text sample contains " << wc << " words" << std::endl;
	
}




int ExtractHead(const int *pos, int *head, int text_size)
{
	int *buffer;
	int nhead=1;
	cudaMalloc(&buffer, sizeof(int)*text_size*2); // this is enough
	thrust::device_ptr<const int> pos_d(pos);
	thrust::device_ptr<int> head_d(head), flag_d(buffer), cumsum_d(buffer+text_size);

	// TODO
	cout << endl << pos << endl;
	//thrust::inclusive_scan(thrust::host, pos_d, pos_d + text_size, head_d);
	
	cudaFree(buffer);
	return nhead;
}

void Part3(char *text, int *pos, int *head, int text_size, int n_head)
{
}
