
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <omp.h>
using namespace std;
#define CUT_OFF 4

bool isSorted(int ref[], int data[], const size_t size) {
	std::sort(ref, ref + size);
	for (size_t idx = 0; idx < size; ++idx) {
		if (ref[idx] != data[idx]) {
			return false;
		}
	}
	return true;
}


/**
  * sequential merge step (straight-forward implementation)
  */
  // TODO: cut-off could also apply here (extra parameter?)
  // TODO: optional: we can also break merge in two halves
void MsMergeSequential(int* out, int* in, long begin1, long end1, long begin2, long end2, long outBegin,int depth) {
	long left = begin1;
	long right = begin2;
	long idx = outBegin;
	while (left < end1 && right < end2) {
		if (in[left] <= in[right]) {
			out[idx] = in[left];
			left++;
		}
		else {
			out[idx] = in[right];
			right++;
		}
		idx++;
	}
	while (left < end1) {
		out[idx] = in[left];
		left++, idx++;
	}

	while (right < end2) {
		out[idx] = in[right];
		right++, idx++;
	}
}
int num = 0;

long binary_search(int* arr, long low, long high, int x) {
	while (low < high) {
		long mid = (low + high) / 2;
		if (arr[mid] <= x) {
			low = mid + 1;
		}
		else {
			high = mid;
		}
	}
	return low;
}
void MsMergeParallel(int* out, int* in, long begin1, long end1, long begin2, long end2, long outBegin,int depth) {
	if (depth < CUT_OFF) {
		long len1 = end1 - begin1;
		long len2 = end2 - begin2;
		if (len1 < len2) {
			std::swap(begin1, begin2);
			std::swap(end1, end2);
			std::swap(len1, len2);
		}
		if (len1 == 0) {
			return;
		}
		long mid1 = (begin1 + end1) / 2;
		long mid2 = binary_search(in, begin2, end2, in[mid1]);
		long outMid = outBegin + (mid1 - begin1) + (mid2 - begin2);
		out[outMid] = in[mid1];
		#pragma omp task
		MsMergeParallel(out, in, begin1, mid1, begin2, mid2, outBegin,depth);  // ×ó
		#pragma omp task
		MsMergeParallel(out, in, mid1 + 1, end1, mid2, end2, outMid + 1,depth);  // ÓÒ
		#pragma omp taskwait
	}
	else {
		MsMergeSequential(out, in, begin1, end1, begin2, end2, outBegin, depth);
	}
}
/**
  * sequential MergeSort
  */
  // TODO: remember one additional parameter (depth)
  // TODO: recursive calls could be taskyfied
  // TODO: task synchronization also is required
void MsSequential(int* array, int* tmp, bool inplace, long begin, long end,long depth) {
	if (begin < (end - 1)) 
	{
		depth++;
		const long half = (begin + end) / 2;
		if (depth <= CUT_OFF) {
			#pragma omp task
			{
				//int execute_thread = omp_get_thread_num(); 
				MsSequential(array, tmp, !inplace, begin, half, depth);
			}
			#pragma omp task
			{
				MsSequential(array, tmp, !inplace, half, end, depth);
			}
			#pragma omp taskwait

			if (inplace) {
				MsMergeSequential(array, tmp, begin, half, half, end, begin, depth);
			}
			else {
				MsMergeSequential(tmp, array, begin, half, half, end, begin, depth);
			}
		}
		else{
			MsSequential(array, tmp, !inplace, begin, half, depth);
			MsSequential(array, tmp, !inplace, half, end, depth);
			if (inplace) {
				MsMergeSequential(array, tmp, begin, half, half, end, begin, depth);
			}
			else {
				MsMergeSequential(tmp, array, begin, half, half, end, begin, depth);
			}
		}
	}
	else if (!inplace) {
		tmp[begin] = array[begin];
	}
}


/**
  * Serial MergeSort
  */
  // TODO: this function should create the parallel region
  // TODO: good point to compute a good depth level (cut-off)
void MsSerial(int* array, int* tmp, const size_t size) {
	#pragma omp parallel
	{
		#pragma omp single			// always using the same thread to create 
		{
			std::cout << "Initial task created by thread " << omp_get_thread_num() << std::endl;
			MsSequential(array, tmp, true, 0, size, 0);
		}
	}// TODO: parallel version of MsSequential will receive one more parameter: 'depth' (used as cut-off)
	
}


int main() {
	int* array, * tmp, * ref;
	size_t size = 100000000;
	array = new int[size];
	tmp = new int[size];
	ref = new int[size];
	srand(997);			// array entries always be same
	for (int i = 0; i < size; i++) {
		array[i] = rand() % 23;
		ref[i] = array[i];
	}
	auto start = std::chrono::high_resolution_clock::now();
	MsSerial(array, tmp, size);
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	std::cout << "compute time: " << duration.count() << " seconds" << std::endl;
	if (isSorted(ref, array, size)) {
		printf("Successful.\n");
	}
	else {
		printf(" FAILED.\n");
	}
	return 0;
}