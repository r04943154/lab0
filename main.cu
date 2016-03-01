#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include <cuda.h>
#include "SyncedMemory.h"
#include "device_launch_parameters.h"
#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include <fstream>

using namespace std;



__global__ void SomeTransform(char *input_gpu, int fsize) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < fsize && input_gpu[idx] != '\n') {
		input_gpu[idx] = input_gpu[idx]+2;//let every char ASCII code plus two
	}
}

int main()
{
	// init, and check
	char file[] = "test.txt";
	
	FILE *fp = fopen("test.txt", "r");
	if (!fp) {
		printf("Cannot open %s", file);
		//abort();
	}
	// get file size
	fseek(fp, 0, SEEK_END);
	size_t fsize = ftell(fp);
	fseek(fp, 0, SEEK_SET);

	// read files
	MemoryBuffer<char> text(fsize + 1);
	auto text_smem = text.CreateSync(fsize);

	fread(text_smem.get_cpu_wo(), 1, fsize, fp);
	text_smem.get_cpu_wo()[fsize] = '\0';
	fclose(fp);

	// TODO: do your transform here
	char *input_gpu = text_smem.get_gpu_rw();
	// An example: transform the first 64 characters to '!'
	// Don't transform over the tail
	// And don't transform the line breaks
	SomeTransform << <2, fsize >> >(input_gpu, fsize);

	puts(text_smem.get_cpu_ro());
	printf("%d" ,text_smem );
	system("pause");
	return 0;
}
