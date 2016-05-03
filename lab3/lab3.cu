#include "lab3.h"
#include <cstdio>

__device__ __host__ int CeilDiv(int a, int b) { return (a-1)/b + 1; }
__device__ __host__ int CeilAlign(int a, int b) { return CeilDiv(a, b) * b; }

__global__ void SimpleClone(
	const float *background,
	const float *target,
	const float *mask,
	float *output,
	const int wb, const int hb, const int wt, const int ht,
	const int oy, const int ox
)
{
	const int yt = blockIdx.y * blockDim.y + threadIdx.y;
	const int xt = blockIdx.x * blockDim.x + threadIdx.x;
	const int curt = wt*yt+xt;
	if (yt < ht && xt < wt && mask[curt] > 127.0f) {
		const int yb = oy+yt, xb = ox+xt;
		const int curb = wb*yb+xb;
		if (0 <= yb && yb < hb && 0 <= xb && xb < wb) {
			output[curb*3+0] = target[curt*3+0];
			output[curb*3+1] = target[curt*3+1];
			output[curb*3+2] = target[curt*3+2];
		}
	}
}

__global__ void CalculateFixed(const float *background, const float *target, const float *mask, float *fixed,
	                           const int wb, const int hb, const int wt, const int ht, const int oy, const int ox){
	float nb[3] = { 0.0, 0.0, 0.0 }, sb[3] = { 0.0, 0.0, 0.0 }, eb[3] = { 0.0, 0.0, 0.0 }, wwb[3] = { 0.0, 0.0, 0.0 };
	float nt[3] = { 0.0, 0.0, 0.0 }, st[3] = { 0.0, 0.0, 0.0 }, et[3] = { 0.0, 0.0, 0.0 }, wwt[3] = { 0.0, 0.0, 0.0 };
	const int yt = blockIdx.y * blockDim.y + threadIdx.y;
	const int xt = blockIdx.x * blockDim.x + threadIdx.x;
	const int curt = wt*yt + xt;
	const int yb = oy + yt, xb = ox + xt;
	const int curb = wb*yb + xb;

	if (mask[curt] == 0){
		//fixed[curt * 3 + 0] = background[curb * 3 + 0];
		//fixed[curt * 3 + 1] = background[curb * 3 + 1];
		//fixed[curt * 3 + 2] = background[curb * 3 + 2];
	}
	else{
		/////N
		if (curt >= wt){
			if (mask[curt - wt] == 255){
				//nb[0] = background[(curb - wb) * 3 + 0];
				//nb[1] = background[(curb - wb) * 3 + 1];
				//nb[2] = background[(curb - wb) * 3 + 2];
			}
			else{
				nb[0] = background[(curb - wb) * 3 + 0];
				nb[1] = background[(curb - wb) * 3 + 1];
				nb[2] = background[(curb - wb) * 3 + 2];
				//nb[0] = fixed[(curt - wt) * 3 + 0];
				//nb[1] = fixed[(curt - wt) * 3 + 1];
				//nb[2] = fixed[(curt - wt) * 3 + 2];
				//mask[curt] = 0;
			}
			nt[0] = target[(curt - wt) * 3 + 0];
			nt[1] = target[(curt - wt) * 3 + 1];
			nt[2] = target[(curt - wt) * 3 + 2];
		}
		else{
			nt[0] = target[curt * 3 + 0];
			nt[1] = target[curt * 3 + 1];
			nt[2] = target[curt * 3 + 2];
			nb[0] = background[(curb - wb) * 3 + 0];
			nb[1] = background[(curb - wb) * 3 + 1];
			nb[2] = background[(curb - wb) * 3 + 2];
		}
		///////////////////////S
		if (curt + wt<wt*ht){
			if (mask[curt + wt] == 255){
				//sb[0] = background[(curb + wb) * 3 + 0];
				//sb[1] = background[(curb + wb) * 3 + 1];
				//sb[2] = background[(curb + wb) * 3 + 2];
			}
			else{
				sb[0] = background[(curb + wb) * 3 + 0];
				sb[1] = background[(curb + wb) * 3 + 1];
				sb[2] = background[(curb + wb) * 3 + 2];
				//sb[0] = fixed[(curt + wt) * 3 + 0];
				//sb[1] = fixed[(curt + wt) * 3 + 1];
				//sb[2] = fixed[(curt + wt) * 3 + 2];
			}
			st[0] = target[(curt + wt) * 3 + 0];
			st[1] = target[(curt + wt) * 3 + 1];
			st[2] = target[(curt + wt) * 3 + 2];
		}
		else{
			st[0] = target[curt * 3 + 0];
			st[1] = target[curt * 3 + 1];
			st[2] = target[curt * 3 + 2];
			sb[0] = background[(curb + wb) * 3 + 0];
			sb[1] = background[(curb + wb) * 3 + 1];
			sb[2] = background[(curb + wb) * 3 + 2];
			//sb[0] = fixed[curt * 3 + 0];
			//sb[1] = fixed[curt * 3 + 1];
			//sb[2] = fixed[curt * 3 + 2];
		}
		///////////////////W
		if (curt%wt != 0){
			if (mask[curt - 1] == 255){
				//wwb[0] = background[(curb - 1) * 3 + 0];
				//wwb[1] = background[(curb - 1) * 3 + 1];
				//wwb[2] = background[(curb - 1) * 3 + 2];
			}
			else{
				wwb[0] = background[(curb - 1) * 3 + 0];
				wwb[1] = background[(curb - 1) * 3 + 1];
				wwb[2] = background[(curb - 1) * 3 + 2];
				//wwb[0] = fixed[(curt - 1) * 3 + 0];
				//wwb[1] = fixed[(curt - 1) * 3 + 1];
				//wwb[2] = fixed[(curt - 1) * 3 + 2];
			}
			wwt[0] = target[(curt - 1) * 3 + 0];
			wwt[1] = target[(curt - 1) * 3 + 1];
			wwt[2] = target[(curt - 1) * 3 + 2];
		}
		else{
			wwt[0] = target[curt * 3 + 0];
			wwt[1] = target[curt * 3 + 1];
			wwt[2] = target[curt * 3 + 2];
			wwb[0] = background[(curb - 1) * 3 + 0];
			wwb[1] = background[(curb - 1) * 3 + 1];
			wwb[2] = background[(curb - 1) * 3 + 2];
			//wwb[0] = fixed[curt * 3 + 0];
			//wwb[1] = fixed[curt * 3 + 1];
			//wwb[2] = fixed[curt * 3 + 2];
		}
		///////////////E
		if ((curt + 1) % wt != 0){
			if (mask[curt + 1] == 255){
				//eb[0] = background[(curb + 1) * 3 + 0];
				//eb[1] = background[(curb + 1) * 3 + 1];
				//eb[2] = background[(curb + 1) * 3 + 2];
			}
			else{
				eb[0] = background[(curb + 1) * 3 + 0];
				eb[1] = background[(curb + 1) * 3 + 1];
				eb[2] = background[(curb + 1) * 3 + 2];
				//eb[0] = fixed[(curt + 1) * 3 + 0];
				//eb[1] = fixed[(curt + 1) * 3 + 1];
				//eb[2] = fixed[(curt + 1) * 3 + 2];
			}
			et[0] = target[(curt + 1) * 3 + 0];
			et[1] = target[(curt + 1) * 3 + 1];
			et[2] = target[(curt + 1) * 3 + 2];
		}
		else{
			et[0] = target[curt * 3 + 0];
			et[1] = target[curt * 3 + 1];
			et[2] = target[curt * 3 + 2];
			eb[0] = background[(curb + 1) * 3 + 0];
			eb[1] = background[(curb + 1) * 3 + 1];
			eb[2] = background[(curb + 1) * 3 + 2];
			//eb[0] = fixed[curt * 3 + 0];
			//eb[1] = fixed[curt * 3 + 1];
			//eb[2] = fixed[curt * 3 + 2];
		}
		fixed[curt * 3 + 0] = 4 * target[curt * 3 + 0] - nt[0] - st[0] - wwt[0] - et[0] + nb[0] + sb[0] + wwb[0] + eb[0];
		fixed[curt * 3 + 1] = 4 * target[curt * 3 + 1] - nt[1] - st[1] - wwt[1] - et[1] + nb[1] + sb[1] + wwb[1] + eb[1];
		fixed[curt * 3 + 2] = 4 * target[curt * 3 + 2] - nt[2] - st[2] - wwt[2] - et[2] + nb[2] + sb[2] + wwb[2] + eb[2];
	}
}


__global__ void PoissonImageCloneing(float *fixed, const float *mask, const float *buf1, float *buf2, const int wt, const int ht)
{
	const int yt = blockIdx.y * blockDim.y + threadIdx.y;
	const int xt = blockIdx.x * blockDim.x + threadIdx.x;
	const int curt = wt*yt + xt;
	//float nb[3] = { 255.0, 255.0, 255.0 }, sb[3] = { 255.0, 255.0, 255.0 }, eb[3] = { 255.0, 255.0, 255.0 }, wwb[3] = { 255.0, 255.0, 255.0 };
	//float nt[3] = { 255.0, 255.0, 255.0 }, st[3] = { 255.0, 255.0, 255.0 }, et[3] = { 255.0, 255.0, 255.0 }, wwt[3] = { 255.0, 255.0, 255.0 };
	float nb[3] = { 0.0, 0.0, 0.0 }, sb[3] = { 0.0, 0.0, 0.0 }, eb[3] = { 0.0, 0.0, 0.0 }, wwb[3] = { 0.0, 0.0, 0.0 };
	float nt[3] = { 0.0, 0.0, 0.0 }, st[3] = { 0.0, 0.0, 0.0 }, et[3] = { 0.0, 0.0, 0.0 }, wwt[3] = { 0.0, 0.0, 0.0 };

	/////N
	if (curt >= wt){
		if (mask[curt - wt] == 255){
			nb[0] = buf1[(curt - wt) * 3 + 0];
			nb[1] = buf1[(curt - wt) * 3 + 1];
			nb[2] = buf1[(curt - wt) * 3 + 2];
		}
		else{
			//nb[0] = fixed[(curt - wt) * 3 + 0];
			//nb[1] = fixed[(curt - wt) * 3 + 1];
			//nb[2] = fixed[(curt - wt) * 3 + 2];
			//mask[curt] = 0;
		}		
	}
	else{
		//nb[0] = buf1[curt * 3 + 0];
		//nb[1] = buf1[curt * 3 + 1];
		//nb[2] = buf1[curt * 3 + 2];
	}
	///////////////////////S
	if (curt + wt<wt*ht){
		if (mask[curt + wt] == 255){
			sb[0] = buf1[(curt + wt) * 3 + 0];
			sb[1] = buf1[(curt + wt) * 3 + 1];
			sb[2] = buf1[(curt + wt) * 3 + 2];
		}
		else{

		}
		
	}
	else{
	
		//sb[0] = buf1[curt * 3 + 0];
		//sb[1] = buf1[curt * 3 + 1];
		//sb[2] = buf1[curt * 3 + 2];
	}
	///////////////////W
	if (curt%wt != 0){
		if (mask[curt - 1] == 255){
			wwb[0] = buf1[(curt - 1) * 3 + 0];
			wwb[1] = buf1[(curt - 1) * 3 + 1];
			wwb[2] = buf1[(curt - 1) * 3 + 2];
		}
		else{
		
		}
	
	}
	else{
	
		//wwb[0] = buf1[curt * 3 + 0];
		//wwb[1] = buf1[curt * 3 + 1];
		//wwb[2] = buf1[curt * 3 + 2];
	}
	///////////////E
	if ((curt + 1) % wt != 0){
		if (mask[curt + 1] == 255){
			eb[0] = buf1[(curt + 1) * 3 + 0];
			eb[1] = buf1[(curt + 1) * 3 + 1];
			eb[2] = buf1[(curt + 1) * 3 + 2];
		}
		else{
			
		}
		
	}
	else{
		
		//eb[0] = buf1[curt * 3 + 0];
		//eb[1] = buf1[curt * 3 + 1];
		//eb[2] = buf1[curt * 3 + 2];
	}
	if (mask[curt] == 255){
		buf2[curt * 3 + 0] = (fixed[curt * 3 + 0] + nb[0] + sb[0] + wwb[0] + eb[0]) / 4;
		buf2[curt * 3 + 1] = (fixed[curt * 3 + 1] + nb[1] + sb[1] + wwb[1] + eb[1]) / 4;
		buf2[curt * 3 + 2] = (fixed[curt * 3 + 2] + nb[2] + sb[2] + wwb[2] + eb[2]) / 4;
	}
	/*
	/////N
	if (curt>=wt){
		if (mask[curt - wt]==255){
			nb[0] = buf1[(curt - wt) * 3 + 0];
			nb[1] = buf1[(curt - wt) * 3 + 1];
			nb[2] = buf1[(curt - wt) * 3 + 2];
		}
		else{
			nb[0] = fixed[(curt - wt) * 3 + 0];
			nb[1] = fixed[(curt - wt) * 3 + 1];
			nb[2] = fixed[(curt - wt) * 3 + 2];
			//mask[curt] = 0;
		}
		nt[0] = fixed[(curt - wt) * 3 + 0];
		nt[1] = fixed[(curt - wt) * 3 + 1];
		nt[2] = fixed[(curt - wt) * 3 + 2];
	}
	else{
		nt[0] = fixed[curt * 3 + 0];
		nt[1] = fixed[curt * 3 + 1];
		nt[2] = fixed[curt * 3 + 2];
		nb[0] = fixed[curt * 3 + 0];
		nb[1] = fixed[curt * 3 + 1];
		nb[2] = fixed[curt * 3 + 2];
	}
	///////////////////////S
	if (curt + wt<wt*ht){
		if (mask[curt + wt]==255){
			sb[0] = buf1[(curt + wt) * 3 + 0];
			sb[1] = buf1[(curt + wt) * 3 + 1];
			sb[2] = buf1[(curt + wt) * 3 + 2];
		}
		else{
			sb[0] = fixed[(curt + wt) * 3 + 0];
			sb[1] = fixed[(curt + wt) * 3 + 1];
			sb[2] = fixed[(curt + wt) * 3 + 2];
		}
		st[0] = fixed[(curt + wt) * 3 + 0];
		st[1] = fixed[(curt + wt) * 3 + 1];
		st[2] = fixed[(curt + wt) * 3 + 2];
	}
	else{
		st[0] = fixed[curt * 3 + 0];
		st[1] = fixed[curt * 3 + 1];
		st[2] = fixed[curt * 3 + 2];
		sb[0] = fixed[curt * 3 + 0];
		sb[1] = fixed[curt * 3 + 1];
		sb[2] = fixed[curt * 3 + 2];
	}
	///////////////////W
	if (curt%wt != 0){
		if (mask[curt - 1]==255){
			wwb[0] = buf1[(curt - 1) * 3 + 0];
			wwb[1] = buf1[(curt - 1) * 3 + 1];
			wwb[2] = buf1[(curt - 1) * 3 + 2];
		}
		else{
			wwb[0] = fixed[(curt - 1) * 3 + 0];
			wwb[1] = fixed[(curt - 1) * 3 + 1];
			wwb[2] = fixed[(curt - 1) * 3 + 2];
		}
		wwt[0] = fixed[(curt - 1) * 3 + 0];
		wwt[1] = fixed[(curt - 1) * 3 + 1];
		wwt[2] = fixed[(curt - 1) * 3 + 2];
	}
	else{
		wwt[0] = fixed[curt * 3 + 0];
		wwt[1] = fixed[curt * 3 + 1];
		wwt[2] = fixed[curt * 3 + 2];
		wwb[0] = fixed[curt * 3 + 0];
		wwb[1] = fixed[curt * 3 + 1];
		wwb[2] = fixed[curt * 3 + 2];
	}
	///////////////E
	if ((curt + 1) % wt != 0){
		if (mask[curt + 1]==255){
			eb[0] = buf1[(curt + 1) * 3 + 0];
			eb[1] = buf1[(curt + 1) * 3 + 1];
			eb[2] = buf1[(curt + 1) * 3 + 2];
		}
		else{
			eb[0] = fixed[(curt + 1) * 3 + 0];
			eb[1] = fixed[(curt + 1) * 3 + 1];
			eb[2] = fixed[(curt + 1) * 3 + 2];
		}
		et[0] = fixed[(curt + 1) * 3 + 0];
		et[1] = fixed[(curt + 1) * 3 + 1];
		et[2] = fixed[(curt + 1) * 3 + 2];
	}
	else{
		et[0] = fixed[curt * 3 + 0];
		et[1] = fixed[curt * 3 + 1];
		et[2] = fixed[curt * 3 + 2];
		eb[0] = fixed[curt * 3 + 0];
		eb[1] = fixed[curt * 3 + 1];
		eb[2] = fixed[curt * 3 + 2];
	}
	if (mask[curt] == 255){
		buf2[curt * 3 + 0] = (4 * buf1[curt * 3 + 0] - nt[0] - st[0] - wwt[0] - et[0] + nb[0] + sb[0] + wwb[0] + eb[0]) / 4;
		buf2[curt * 3 + 1] = (4 * buf1[curt * 3 + 1] - nt[1] - st[1] - wwt[1] - et[1] + nb[1] + sb[1] + wwb[1] + eb[1]) / 4;
		buf2[curt * 3 + 2] = (4 * buf1[curt * 3 + 2] - nt[2] - st[2] - wwt[2] - et[2] + nb[2] + sb[2] + wwb[2] + eb[2]) / 4;
	}
	*/
}

void PoissonImageCloning(
	const float *background,
	const float *target,
	const float *mask,
	float *output,
	const int wb, const int hb, const int wt, const int ht,
	const int oy, const int ox
)
{
	float *fixed, *buf1, *buf2;
	cudaMalloc(&fixed, 3 * wt*ht*sizeof(float));
	cudaMalloc(&buf1, 3 * wt*ht*sizeof(float));
	cudaMalloc(&buf2, 3 * wt*ht*sizeof(float));

	// intialize
	dim3 gdim(CeilDiv(wt, 32), CeilDiv(ht, 16)), bdim(32, 16);
	CalculateFixed<<<gdim, bdim>>>(background, target, mask, fixed, wb, hb, wt, ht, oy, ox);
	cudaMemcpy(buf1, target, wt*ht*sizeof(float) * 3, cudaMemcpyDeviceToDevice);

	//iterate
	for (int i = 0; i < 10000; i++){
		PoissonImageCloneing<<<gdim, bdim>>>(fixed, mask, buf1, buf2, wt, ht);
		PoissonImageCloneing<<<gdim, bdim>>>(fixed, mask, buf2, buf1, wt, ht);
	}


	cudaMemcpy(output, background, wb*hb*sizeof(float)*3, cudaMemcpyDeviceToDevice);
	SimpleClone<<<dim3(CeilDiv(wt,32), CeilDiv(ht,16)), dim3(32,16)>>>(
		background, buf1, mask, output,
		wb, hb, wt, ht, oy, ox
	);
	
}
