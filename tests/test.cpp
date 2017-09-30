/*
 * This file demonstrates how TENSOR4X as a template class can be used.
 */

#include <iostream>
#include <vector>
#include <cassert>
#include "../include/tensor4x.hpp"

int main(int argc, char *argv[])
{
    const int N = 4;

	// Constructe a 3*3*3*3 tensor with all elements zeros
    TENSOR4X<double> V(N);
    // Set some elements
    V(0,0,0,0) = 0.5;
    V(0,1,0,1) = 1.;    // set 0101, 1001, 0110 and 1010
    V(0,1,2,3) = 1.5;   // set 0123, 0132, 1023, 1032,
                        //     2301, 2310, 3201 and 3210
    // check symmetry is indeed enabled
    assert(V(0,1,0,1) == V(1,0,0,1));

	// or print all non-zero elements
	std::cout << "V:" << std::endl << V;

	// Alternatively, we can construct it from std:vector
	// 55 is the number of unique elements in a 4*4*4*4
	// four-tensor of 8-fold symmetry.
	std::vector<double> a(55, 0.);
	for (size_t i = 0, ij = 0; i < N; i++)
	for (size_t j = 0; j <= i; j++, ij++)
	for (size_t k = 0, kl = 0; k < N; k++)
	for (size_t l = 0; l <= k; l++, kl++)
		if (ij <= kl)
			a[cpind(ij,kl)] = i+j+k+l;
	TENSOR4X<double> A(N, a);
	
	std::cout << "A:" << std::endl << A;

	// You can also perform elementwise add/subtraction
	TENSOR4X<double> B = A - V;

    return 0;
}
