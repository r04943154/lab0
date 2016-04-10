#include "timer.h"   

#include "lab2.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"


static const unsigned W = 640;
static const unsigned H = 480;
static const unsigned NFRAME = 240;


__global__ void simple_kernel(uint8_t *pos, unsigned int width, unsigned int height, float time)
{
	unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int y = blockIdx.y*blockDim.y + threadIdx.y;

	
	float u = x / (float)width;
	float v = y / (float)height;
	u = u*2.0f - 1.0f;
	v = v*2.0f - 1.0f;

	
	float freq = 4.0f;
	float w = sinf(u*freq + time/30) * cosf(v*freq + time/30) * 2.5f;

	if (time>120)w = sinf(u*freq + time/60) * cosf(v*freq + time/60) * 40.5f;

	
	//pos[y*width + x] = make_float4(u, w, v, 1.0f);
	pos[y*width + x] = (uint8_t)(u + v + w )*1000 % 256;
	if (time<40||time>200)pos[y*width + x] = (uint8_t)(x+y+time) % 256;
	//printf("=%f %f %f=\n",u , v, w);
}


struct Lab2VideoGenerator::Impl {
	int t = 0;
};

Lab2VideoGenerator::Lab2VideoGenerator(): impl(new Impl) {
}

Lab2VideoGenerator::~Lab2VideoGenerator() {}

void Lab2VideoGenerator::get_info(Lab2VideoInfo &info) {
	info.w = W;
	info.h = H;
	info.n_frame = NFRAME;
	// fps = 24/1 = 24
	info.fps_n = 24;
	info.fps_d = 1;
};


void Lab2VideoGenerator::Generate(uint8_t *yuv) {

	
	//float t = (float)time();

	dim3 block(8, 8, 1);
	dim3 grid(W / block.x, H*1.5 / block.y, 1);
	simple_kernel << < grid, block >> >(yuv, W, H, impl->t);

	//cudaMemset(yuv, (impl->t), W*H/2);
	//cudaMemset(yuv + W*H / 2, ((impl->t)+128)%256, W*H/2);
	//cudaMemset(yuv + W*H, rand(), W*H / 2);
	//impl->t = rand() % 256;
	impl->t += 1;
}
