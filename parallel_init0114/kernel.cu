
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>
#include <thrust\host_vector.h>
#include <thrust\device_vector.h>
#include <thrust\transform_reduce.h>
#include <utility>
#include <ctime>
#include <device_functions.h>
#include <stdio.h>
#include <vector>
#include <thrust/random.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/tuple.h>
#include <fstream>

using namespace std;

using namespace thrust;

template <typename T>
struct square
{
	__host__ __device__
		T operator()(const T& x) const {
		return x * x;
	}
};

struct prg
{
	float a, b;

	__host__ __device__
		prg(float _a, float _b) : a(_a), b(_b) {};

	__host__ __device__
		float operator()(const unsigned int n) const
	{
		thrust::default_random_engine rng;
		thrust::uniform_real_distribution<float> dist(a, b);
		rng.discard(n);

		return dist(rng);
	}
};


struct Point
{
	float x;
	float y;
	Point() :x(0.0f), y(0.0f) {}

	__host__ __device__
		Point(float h, float w) : x(h), y(w) {}
};



//__shared__ PointTuple *raw_ptr;

__constant__ int wNum[1];
__constant__ int hNum[1];
__constant__ float width[1];
__constant__ float radius[1];//R
__constant__ int cinit[1];//³õÊ¼»¯µãÊýÁ¿
#define M_PI 3.141592654f

cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);
cudaError_t init_point(thrust::device_ptr<Point>, unsigned int);
cudaError_t randomInsertPoints(unsigned int);


__global__ void addKernel(int *c, const int *a, const int *b)
{
	int i = blockIdx.x;
	c[i] = a[i] + b[i];
}

__global__ void testPointer(Point *raw_ptr) {
	//	cout << "pointer[0]" << *a<<endl;

	extern __shared__ Point dev_zip[];
	extern __shared__ int sum[];
	sum[0] = 0;
	int j = threadIdx.x;
	int i = threadIdx.x + blockDim.x*blockIdx.x;

	if (i < 300000) dev_zip[j] = raw_ptr[i];

	__syncthreads();



	atomicAdd(&sum[0], 1);
}



__global__ void testNum(int *raw_ptr) {

	/*extern __shared__ int share_nums[];
	extern __shared__ int sum[];
	sum[0] = 0;
	int j = threadIdx.x;
	for (int i = j; i < 30; i += (512))
	{
	share_nums[i] = raw_ptr[i];
	}

	__syncthreads();


	int num = 0;
	for (; j < 30; j += (512))
	{
	if (share_nums[j])
	num++;
	}
	atomicAdd(&sum[0], num);*/

	extern __shared__ int sm[];

	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int j = threadIdx.x;
	if (i<300000) sm[j] = raw_ptr[i];

	__syncthreads();


}




__global__ void addWithSM(int *c, const int *a) {
	int i = threadIdx.x;
	extern __shared__ int seme[];
	seme[i] = a[i];
	__syncthreads();

	if (i == 0) {
		c[0] = 0;
		for (int d = 0; d < 5; d++) {
			c[0] += seme[d] * seme[d];
		}
	}
	if (i == 1) {
		c[1] = 0;
		for (int d = 0; d < 5; d++) {
			c[i] += seme[d];
		}
	}
	if (i == 2)
	{
		c[i] = 1;
		for (int d = 0; d < 5; d++) {
			c[2] *= seme[d];
		}
	}

}

__global__ void addKernel_thd(int *c, int *a, int *b) {
	int i = threadIdx.x;
	c[i] = a[i] + b[i];
}


__global__ void initpoints(float *h_ptr, float *w_ptr, Point *init_ptr) {
	int i = blockIdx.x;
	int j = threadIdx.x;
	int index = i*blockDim.x + j;
	if (index < cinit[0]) {
		init_ptr[index] = Point(h_ptr[index], w_ptr[index]);
	}
}




__global__ void insertPointsErase(Point *raw_ptr, bool *erase_ptr, float *ran_ptr, Point *points_ptr, bool *flag_ptr, int *neis_ptr) {


	extern __shared__ Point p[];
	extern __shared__ int save[];
	extern __shared__ int n[];
	//extern __shared__ int value[];

	//thrust::default_random_engine rng;
	//thrust::uniform_real_distribution<float> range(0.0f, 1000.0f);
	//rng.discard(i);

	int i = blockIdx.x;
	int j = threadIdx.x;



	//if (i == j){
	//	ran_ptr[i] = range(rng);
	//	//int index = ran_ptr[i] / width[0];
	//}
	p[j] = raw_ptr[j];

	__syncthreads();

	if (i != j) {
		float dis = sqrt(pow(p[j].x - p[i].x, 2) + pow(p[j].y - p[i].y, 2));
		if (dis < radius[0]) {
			if (ran_ptr[i] < ran_ptr[j] && !erase_ptr[j]) {
				erase_ptr[i] = true;
			}
			//erase_ptr[i] = true;
		}
	}

	if (j < 6) {
		if (j == 0) n[0] = 0;
		save[j] = -1;
		//value[j] = -1;
	}
	__syncthreads();



	if (!erase_ptr[i]) {
		/*int x = p[i].x / width[0];
		int y = p[i].y / width[0];

		int index = x*wNum[0]+y;*/

		for (size_t k = j; k < cinit[0]; k += 512)
		{
			/*int hi = k / wNum[0];
			int wi = k - hi*wNum[0];
			if ((hi - x <= 2 && hi - x>-2) && (wi - y <= 2&&wi-y>=-2)){

			}*/
			float dX = points_ptr[k].x - p[i].x;
			float dY = points_ptr[k].y - p[i].y;
			float dis = sqrt(pow(dX, 2) + pow(dY, 2));

			if (dis<radius[0] && !flag_ptr[k]) {


				atomicAdd(&n[0], 1);

				float vectorA[] = { -1, 0 };
				float vectorB[] = { dX, dY };


				double dotProduct = -1 * dX;
				double cos = dotProduct / (1 * sqrt(pow(dX, 2) + pow(dY, 2)));



				double res = acos(cos);
				res = dY > 0 ? res : 2 * M_PI - res;

				int s = floor(res / (M_PI / 3));//·ÖÇø
				s = s > 5 ? 5 : s;
				save[s] = k;
				flag_ptr[k] = true;
			}
		}
	}
	__syncthreads();
	if (j<6) {
		if (save[j] != -1) flag_ptr[save[j]] = false;
	}
	if (i == j) {
		neis_ptr[j] = n[0];
		if (n[0] == 1) {
			for (size_t t = 0; t < 6; t++) {
				if (save[t] != -1) points_ptr[save[t]] = raw_ptr[i];
			}
		}
	}

}




float window_width = 350.0;
float window_height = 350.0;

const size_t N = 1000;

const unsigned int numRandom = 300000;
int init[1];

float ratio = 0.75;
float R[1];//¶¨Òå±ê×¼µã¾à




		   //float u = 1.0;//¶¨Òåµ¥Î»¾àÀë



float d[1]; //¸ñ×Ó±ß³¤


int h[1];
int w[1];

int main()
{
	ratio = 0.75;
	float D = 2 * sqrt(window_width*window_height / (2 * N*sqrt(3)));
	R[0] = ratio * D;
	cout << R[0] << endl;
	d[0] = sqrt(2)*R[0] / 2;

	h[0] = ceil(window_height / d[0]);
	w[0] = ceil(window_width / d[0]);
	init[0] = numRandom;
	//copy host  data to constant memory
	cudaMemcpyToSymbol(hNum, h, sizeof(int));
	cudaMemcpyToSymbol(wNum, w, sizeof(int));
	cudaMemcpyToSymbol(width, d, sizeof(float));
	cudaMemcpyToSymbol(radius, R, sizeof(float));
	cudaMemcpyToSymbol(cinit, init, sizeof(int));


	cudaError_t cudaStatus;
	int num = 0;
	cudaDeviceProp prop;
	cudaStatus = cudaGetDeviceCount(&num);

	for (int i = 0; i<num; i++)
	{
		cudaGetDeviceProperties(&prop, i);
	}




	cudaStatus = randomInsertPoints(512);
	// Add vectors in parallel.
	// cudaStatus = addWithCuda(c, a, b, arraySize);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "compute failed!");
		return 1;
	}


	// cudaDeviceReset must be called before exiting in order for profiling and
	// tracing tools such as Nsight and Visual Profiler to show complete traces.
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceReset failed!");
		return 1;
	}

	/*for (int i = 0; i<arraySize; i++)
	{
	if (c[i] != (a[i] + b[i]))
	{
	printf("Error in %d\n", i);
	}
	}*/
	return 0;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size)
{
	int *dev_a = 0;
	int *dev_b = 0;
	int *dev_c = 0;
	cudaError_t cudaStatus;


	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	// Allocate GPU buffers for three vectors (two input, one output)    .
	cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

	cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}



	//time
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	for (int i = 0; i < 1000; i++) {
		addKernel_thd << <1, 512 >> >(dev_c, dev_a, dev_b);

	}


	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float tm;
	cudaEventElapsedTime(&tm, start, stop);
	printf("GPU ll time:%.6f ms.\n", tm);


	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		goto Error;
	}

	// Copy output vector from GPU buffer to host memory.
	cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

Error:
	cudaFree(dev_c);
	cudaFree(dev_a);
	cudaFree(dev_b);

	return cudaStatus;
}





cudaError_t randomInsertPoints(const unsigned int size)
{

	cudaError_t cudaStatus;
	bool *erase_ptr;


	thrust::device_ptr<Point> dev_init = thrust::device_malloc<Point>(numRandom);

	cudaStatus = init_point(dev_init, numRandom);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "init points error!\n", cudaStatus);
		return cudaStatus;
		//goto Errors;
	}


	Point *points_ptr = thrust::raw_pointer_cast(dev_init);

	/*for (size_t i = 0; i < numRandom; i+=10000)
	{
	cout << points_ptr[i].x << " " << points_ptr[i].y << endl;
	}*/
	//random insert points set




	cudaMalloc((void**)&erase_ptr, size * sizeof(bool));
	thrust::device_ptr<bool> dev_erase_ptr(erase_ptr);



	thrust::device_ptr<bool>dev_flag_ptr = thrust::device_malloc<bool>(numRandom);
	thrust::fill(dev_flag_ptr, dev_flag_ptr + numRandom, false);
	bool *flag_ptr = thrust::raw_pointer_cast(dev_flag_ptr);



	thrust::device_vector<float> x(size);
	thrust::device_vector<float> y(size);

	//threads'random value sets
	thrust::device_vector<float> ran(size);

	//threads' neighbor points to erase
	//thrust::device_vector<int> neis(size);
	thrust::device_ptr<int>dev_nei_ptr = thrust::device_malloc<int>(size);
	int *neis_ptr = thrust::raw_pointer_cast(dev_nei_ptr);


	int times = 0;

	int flagnum = 0;
	bool flag[numRandom];
	Point *points;
	points = (Point*)malloc(numRandom * sizeof(Point));

	string file = "parallel";
	

	ofstream outfile;


	while (true)
	{
		times++;

		thrust::fill(dev_erase_ptr, dev_erase_ptr + size, false);

		thrust::fill(dev_nei_ptr, dev_nei_ptr + size, -1);

		prg rx = prg(0.0f, window_height);
		unsigned int offset = time(NULL);
		thrust::counting_iterator<unsigned int> index_sequence_beginx(offset), index_sequence_beginy(offset + times);
		thrust::transform(index_sequence_beginx,
			index_sequence_beginx + size,
			x.begin(),
			rx);

		prg ry = prg(0.0f, window_width);
		thrust::transform(index_sequence_beginy,
			index_sequence_beginy + size,
			y.begin(),
			ry);




		thrust::counting_iterator<unsigned int> index_sequence_begin_ran(offset + 1000 + times);
		prg rran = prg(0.0f, 1000.0f);
		thrust::transform(index_sequence_begin_ran,
			index_sequence_begin_ran + size,
			ran.begin(),
			rran);


		float *ran_ptr = thrust::raw_pointer_cast(ran.data());

		thrust::device_vector<Point> dev_tuvec;

		for (int i = 0; i < size; i++) {
			//cout << x[i] << " " << y[i]<<endl;
			int xIndex = x[i] / d[0];
			int yIndex = y[i] / d[0];
			if (xIndex >= h[0] || yIndex >= w[0]) {
				cout << h[0] << " " << w[0] << " " << d[0] << endl;
				std::system("pause");
			}
			dev_tuvec.push_back(Point(x[i], y[i]));

		}


		Point *raw_ptr = thrust::raw_pointer_cast(dev_tuvec.data());





		//time
		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start, 0);

		insertPointsErase << < size, size, size * sizeof(Point) + 7 * sizeof(int), 0 >> >(raw_ptr, erase_ptr, ran_ptr, points_ptr, flag_ptr, neis_ptr);


		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		float tm;
		cudaEventElapsedTime(&tm, start, stop);
		printf("GPU Elapsed time:%.6f ms.\n", tm);


		// Check for any errors launching the kernel
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			goto Errors;
		}

		// cudaDeviceSynchronize waits for the kernel to finish, and returns
		// any errors encountered during the launch.
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
			goto Errors;
		}



		bool erase[512];
		cudaStatus = cudaMemcpy(erase, erase_ptr, size * sizeof(bool), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Errors;
		}


		cudaStatus = cudaMemcpy(flag, flag_ptr, numRandom * sizeof(bool), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Errors;
		}



		cudaStatus = cudaMemcpy(points, points_ptr, numRandom * sizeof(Point), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Errors;
		}
		// copy memory to a new device_vector (which automatically allocates memory)
		//thrust::device_vector<bool> vec(dev_erase_ptr, dev_erase_ptr + size);
		/*int eraseNum = 0;
		for (int i = 0; i < size; i++){
		if (erase[i]) eraseNum++;
		}
		cout << "erase£º" << eraseNum << endl;
		*/


		int eraseNum = 0;
		int neis = 0;
		for (int i = 0; i < size; i++) {
			if (dev_erase_ptr[i]) eraseNum++;
			if (dev_nei_ptr[i] == 1) neis++;
		}
		cout << "erase" << eraseNum << endl;
		cout << "only one neighbors" << neis << endl;


		flagnum = 0;
		for (size_t i = 0; i < numRandom; i++)
		{
			if (flag[i]) {
				flagnum++;
				//cout << i << endl;
			}
		}
		cout << times << " " << flagnum << endl;
		if (flagnum >= numRandom - 1000) break;

	}
	file.append(to_string(flagnum));
	file.append(".txt");

	outfile.open(file);
	outfile << numRandom - flagnum << endl;
	for (size_t i = 0; i < numRandom; i++)
	{
		if (!flag[i]) {
			//cout << i << endl;
			outfile << points[i].x / window_height << " " << points[i].y / window_width << endl;
		}
	}
	outfile.close();

Errors:
	cudaFree(erase_ptr);


	return cudaStatus;
}



cudaError_t init_point(thrust::device_ptr<Point>dev_init, unsigned int size) {



	cudaError_t cudaStatus;




	//random points init
	//thrust::device_vector<Point> initpoints;
	unsigned int off = time(NULL);
	thrust::counting_iterator<unsigned int> index_sequence_begin_h(off), index_sequence_begin_w(off + 1);
	thrust::device_vector<float> h(size);
	thrust::device_vector<float> w(size);

	prg rh = prg(0.0f, window_height);
	thrust::transform(index_sequence_begin_h,
		index_sequence_begin_h + size,
		h.begin(),
		rh);

	prg rw = prg(0.0f, window_width);
	thrust::transform(index_sequence_begin_w,
		index_sequence_begin_w + size,
		w.begin(),
		rw);


	cudaEvent_t start_host, stop_host;
	cudaEventCreate(&start_host);
	cudaEventCreate(&stop_host);
	cudaEventRecord(start_host, 0);




	float *h_ptr = thrust::raw_pointer_cast(h.data());
	float *w_ptr = thrust::raw_pointer_cast(w.data());



	//thrust::device_ptr<Point> dev_init = thrust::device_malloc<Point>(numRandom);
	Point *init_ptr = thrust::raw_pointer_cast(dev_init);

	size_t block = ceil(size / 512);
	initpoints << <block, 512 >> >(h_ptr, w_ptr, init_ptr);



	cudaEventRecord(stop_host, 0);
	cudaEventSynchronize(stop_host);
	float tm_host;
	cudaEventElapsedTime(&tm_host, start_host, stop_host);
	printf("Cpu Elapsed time:%.6f ms.\n", tm_host);




	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "init error! cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		return cudaStatus;
	}


	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "init error! launch failed: %s\n", cudaGetErrorString(cudaStatus));
		return cudaStatus;
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "init error! cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);

		return cudaStatus;
	}

	return cudaStatus;
}