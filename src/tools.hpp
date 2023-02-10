#ifndef MATRIX_MUL_MOD_TOOLS_H
#define MATRIX_MUL_MOD_TOOLS_H

#include <vector>
#include <cstdint>
#include <cassert>

typedef unsigned __int128 uint128_t;

void mul_matrix_matrix_mod(uint64_t* A, uint64_t A_rows, uint64_t A_cols, uint64_t A_sub_rows, uint64_t A_sub_cols,
							uint64_t* B, uint64_t B_rows, uint64_t B_cols, uint64_t B_sub_rows, uint64_t B_sub_cols,
							uint64_t* result, uint64_t N, std::vector<uint64_t> coeff_moduli, std::vector<uint128_t> m);

uint64_t mulmod(uint64_t a, uint64_t b, uint64_t n, uint128_t m);

uint64_t addmod(uint64_t a, uint64_t b, uint64_t n);

#endif