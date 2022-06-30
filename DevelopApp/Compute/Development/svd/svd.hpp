//#define PRINT_DEBUGGING_OUTPUT
#define USE_SCALAR_IMPLEMENTATION
//#define USE_ACCURATE_RSQRT_IN_JACOBI_CONJUGATION
//#define USE_SSE_IMPLEMENTATION
//#define USE_AVX_IMPLEMENTATION

#define COMPUTE_V_AS_MATRIX
//#define COMPUTE_V_AS_QUATERNION
#define COMPUTE_U_AS_MATRIX
//#define COMPUTE_U_AS_QUATERNION

#include "Singular_Value_Decomposition_Preamble.hpp"

inline void computeMat3SVD(double* M, double* U, double* V, double* A)
{
#include "Singular_Value_Decomposition_Kernel_Declarations.hpp"
	ENABLE_SCALAR_IMPLEMENTATION(Sa11.f = M[0];)                                        ENABLE_SSE_IMPLEMENTATION(Va11 = _mm_set1_ps(M[0]);)                                    ENABLE_AVX_IMPLEMENTATION(Va11 = _mm256_set1_ps(M[0]);)                                 ENABLE_AVX512_IMPLEMENTATION(Va11 = _mm512_set1_ps(M[0]);)
	ENABLE_SCALAR_IMPLEMENTATION(Sa21.f = M[1];)                                        ENABLE_SSE_IMPLEMENTATION(Va21 = _mm_set1_ps(M[1]);)                                    ENABLE_AVX_IMPLEMENTATION(Va21 = _mm256_set1_ps(M[1]);)                                 ENABLE_AVX512_IMPLEMENTATION(Va21 = _mm512_set1_ps(M[1]);)
	ENABLE_SCALAR_IMPLEMENTATION(Sa31.f = M[2];)                                        ENABLE_SSE_IMPLEMENTATION(Va31 = _mm_set1_ps(M[2]);)                                    ENABLE_AVX_IMPLEMENTATION(Va31 = _mm256_set1_ps(M[2]);)                                 ENABLE_AVX512_IMPLEMENTATION(Va31 = _mm512_set1_ps(M[2]);)
	ENABLE_SCALAR_IMPLEMENTATION(Sa12.f = M[3];)                                        ENABLE_SSE_IMPLEMENTATION(Va12 = _mm_set1_ps(M[3]);)                                    ENABLE_AVX_IMPLEMENTATION(Va12 = _mm256_set1_ps(M[3]);)                                 ENABLE_AVX512_IMPLEMENTATION(Va12 = _mm512_set1_ps(M[3]);)
	ENABLE_SCALAR_IMPLEMENTATION(Sa22.f = M[4];)                                        ENABLE_SSE_IMPLEMENTATION(Va22 = _mm_set1_ps(M[4]);)                                    ENABLE_AVX_IMPLEMENTATION(Va22 = _mm256_set1_ps(M[4]);)                                 ENABLE_AVX512_IMPLEMENTATION(Va22 = _mm512_set1_ps(M[4]);)
	ENABLE_SCALAR_IMPLEMENTATION(Sa32.f = M[5];)                                        ENABLE_SSE_IMPLEMENTATION(Va32 = _mm_set1_ps(M[5]);)                                    ENABLE_AVX_IMPLEMENTATION(Va32 = _mm256_set1_ps(M[5]);)                                 ENABLE_AVX512_IMPLEMENTATION(Va32 = _mm512_set1_ps(M[5]);)
	ENABLE_SCALAR_IMPLEMENTATION(Sa13.f = M[6];)                                        ENABLE_SSE_IMPLEMENTATION(Va13 = _mm_set1_ps(M[6]);)                                    ENABLE_AVX_IMPLEMENTATION(Va13 = _mm256_set1_ps(M[6]);)                                 ENABLE_AVX512_IMPLEMENTATION(Va13 = _mm512_set1_ps(M[6]);)
	ENABLE_SCALAR_IMPLEMENTATION(Sa23.f = M[7];)                                        ENABLE_SSE_IMPLEMENTATION(Va23 = _mm_set1_ps(M[7]);)                                    ENABLE_AVX_IMPLEMENTATION(Va23 = _mm256_set1_ps(M[7]);)                                 ENABLE_AVX512_IMPLEMENTATION(Va23 = _mm512_set1_ps(M[7]);)
	ENABLE_SCALAR_IMPLEMENTATION(Sa33.f = M[8];)                                        ENABLE_SSE_IMPLEMENTATION(Va33 = _mm_set1_ps(M[8]);)                                    ENABLE_AVX_IMPLEMENTATION(Va33 = _mm256_set1_ps(M[8]);)                                 ENABLE_AVX512_IMPLEMENTATION(Va33 = _mm512_set1_ps(M[8]);)

#include "Singular_Value_Decomposition_Main_Kernel_Body.hpp"
	ENABLE_SCALAR_IMPLEMENTATION(U[0] = Su11.f;)
	ENABLE_SCALAR_IMPLEMENTATION(U[1] = Su21.f;)
	ENABLE_SCALAR_IMPLEMENTATION(U[2] = Su31.f;)
	ENABLE_SCALAR_IMPLEMENTATION(U[3] = Su12.f;)
	ENABLE_SCALAR_IMPLEMENTATION(U[4] = Su22.f;)
	ENABLE_SCALAR_IMPLEMENTATION(U[5] = Su32.f;)
	ENABLE_SCALAR_IMPLEMENTATION(U[6] = Su13.f;)
	ENABLE_SCALAR_IMPLEMENTATION(U[7] = Su23.f;)
	ENABLE_SCALAR_IMPLEMENTATION(U[8] = Su33.f;)

	ENABLE_SCALAR_IMPLEMENTATION(V[0] = Sv11.f;)
	ENABLE_SCALAR_IMPLEMENTATION(V[1] = Sv21.f;)
	ENABLE_SCALAR_IMPLEMENTATION(V[2] = Sv31.f;)
	ENABLE_SCALAR_IMPLEMENTATION(V[3] = Sv12.f;)
	ENABLE_SCALAR_IMPLEMENTATION(V[4] = Sv22.f;)
	ENABLE_SCALAR_IMPLEMENTATION(V[5] = Sv32.f;)
	ENABLE_SCALAR_IMPLEMENTATION(V[6] = Sv13.f;)
	ENABLE_SCALAR_IMPLEMENTATION(V[7] = Sv23.f;)
	ENABLE_SCALAR_IMPLEMENTATION(V[8] = Sv33.f;)

	ENABLE_SCALAR_IMPLEMENTATION(A[0] = Sa11.f;)
	ENABLE_SCALAR_IMPLEMENTATION(A[1] = Sa22.f;)
	ENABLE_SCALAR_IMPLEMENTATION(A[2] = Sa33.f;)
}