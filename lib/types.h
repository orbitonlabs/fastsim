#pragma once

#include <utility>
#include <valarray>
#include <functional>
#include <iostream>
#include <cxxplot/cxxplot>
#include "range.h"

using std::function;
using std::valarray;
using util::lang::indices;

typedef std::function<double(valarray<double>, double)> vector_function_type;
typedef valarray<vector_function_type> functional_vector;
typedef valarray<double> vector_type;
typedef std::function<int(vector_type&, double)> custom_parametric_function;
typedef valarray<valarray<double>> matrix_type;

static double zero(vector_type v, double t) {
    return 0;
}

static vector_type eval(const functional_vector& f, const vector_type& x, double t) {
	vector_type v(0.0, f.size());
	for (size_t i = 0; i != f.size(); i++) v[i] = f[i](x, t);
	return v;
}

namespace vector {
	static double dot(const vector_type& v1, const vector_type& v2) {
		vector_type v = v1 * v2;
		double dp = 0.0;
		for (auto e : v) dp += e;
		return dp;
	}

	static void print(const vector_type& v) {
		std::cout << "[";
		for (auto e : v) std::cout << e << ", ";
		std::cout << "]\n";
	}
}

namespace matrix {

	class Matrix {
	private:
		matrix_type m;
		size_t rows;
		size_t cols;

	public:
		Matrix(matrix_type matrix, const size_t r, const size_t c) {
			m = std::move(matrix);
			rows = r;
			cols = c;
		}

		explicit Matrix(vector_type v) {
			m = matrix_type({0.0}, v.size());
			for (size_t i = 0; i < v.size(); i++) m[i] = v[i];
			rows = v.size();
			cols = 1;
		}

		Matrix(const size_t r, const size_t c) {
			auto vec = vector_type(0.0, c);
			m = matrix_type(vec, r);
			rows = r;
			cols = c;
		}

		Matrix operator + (const Matrix& matrix) {
			return {m + matrix.m, rows, cols};
		}

		Matrix operator - (const Matrix& matrix) {
			return {m - matrix.m, rows, cols};
		}

		Matrix transpose() {
			Matrix t(cols, rows);
			for (size_t ri = 0; ri < rows; ri++) {
				for (size_t ci = 0; ci < cols; ci++) {
					t.m[ci][ri] = m[ri][ci];
				}
			}
			return t;
		}

		Matrix operator * (Matrix& matrix) {
			Matrix r(rows, matrix.cols);
			auto t = matrix.transpose();
			for (size_t ri = 0; ri < rows; ri++) {
				for (size_t ci = 0; ci < matrix.cols; ci++) {
					r.m[ri][ci] = vector::dot(m[ri], t.m[ci]);
				}
			}
			return r;
		}

		vector_type dot(const vector_type& v) {
			vector_type r(v.size());
			for (size_t i = 0; i < rows; i++) {
				r[i] = vector::dot(m[i], v);
			}
			return r;
		}

		void print() {
			for (auto a : indices(m)) {
				auto vec = m[a];
				for (auto b : indices(vec))
					std::cout << vec[b] << "\t";
				std::cout << "\n";
			}
		}

		// http://dx.doi.org/10.4018/jtd.2010010102
		static double inverse(const Matrix& a, Matrix& res, size_t size) {
			double pivot, det = 1.0;
			size_t i, j, p;
			res.m = a.m;
			for (p = 0; p < size; p++) {
				pivot = res.m[p][p];
				det = det * pivot;
				if (fabs(pivot) < 1e-5) return 0;
				for (i = 0; i < size; i++)
					res.m[i][p] = -res.m[i][p] / pivot;
				for (i = 0; i < size; i++)
					if (i != p)
						for (j = 0; j < size; j++)
							if (j != p)
								res.m[i][j] = res.m[i][j] + res.m[p][j] * res.m[i][p];
				for (j = 0; j < size; j++)
					res.m[p][j] = res.m[p][j] / pivot;
				res.m[p][p] = 1 / pivot;
			}
			return det;
		}

		static double geninv(Matrix& a, Matrix& geninv) {
			auto at = a.transpose();
			auto sq = a * at;
			const size_t size = sq.rows;
			Matrix inv(size, size);
			double res = inverse(sq, inv, size);
			geninv = at * inv;
			return res;
		}
	};

}