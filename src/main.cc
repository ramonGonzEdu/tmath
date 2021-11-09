#include <array>
#include <cmath>
#include <concepts>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

template <typename T>
concept Matrixable = requires(T a, T b) {
  { a *b } -> std::convertible_to<T>;
  { a + b } -> std::convertible_to<T>;
};

template <int M = 1, int N = 1, Matrixable Numeric = double> struct Matrix {

  std::array<std::array<Numeric, N>, M> data{};

#pragma region matrix_operations
  template <int P>
  constexpr Matrix<M, P, Numeric> operator*(Matrix<N, P, Numeric> other) {
    Matrix<M, P, Numeric> output;

    for (std::size_t m = 0; m < M; m++) {
      for (std::size_t p = 0; p < P; p++) {
        // Do cell m,p of output matrix
        for (std::size_t n = 0; n < N; n++) {
          output.data[m][p] += data[m][n] * other.data[n][p];
        }
      }
    }

    return output;
  }

  constexpr Matrix<M, N, Numeric> operator+(Matrix<M, N, Numeric> other) {
    Matrix<M, N, Numeric> output;

    for (std::size_t m = 0; m < M; m++) {
      for (std::size_t n = 0; n < N; n++) {
        // Do cell m,n of output matrix
        output.data[m][n] = data[m][n] + other.data[m][n];
      }
    }

    return output;
  }

  constexpr Matrix<M, N, Numeric> operator-(Matrix<M, N, Numeric> other) {
    Matrix<M, N, Numeric> output;

    for (std::size_t m = 0; m < M; m++) {
      for (std::size_t n = 0; n < N; n++) {
        // Do cell m,n of output matrix
        output.data[m][n] = data[m][n] - other.data[m][n];
      }
    }

    return output;
  }

  template <int P> void operator*=(Matrix<N, N, Numeric> other) {
    Matrix<M, N, Numeric> output;

    for (std::size_t m = 0; m < M; m++) {
      for (std::size_t p = 0; p < P; p++) {
        // Do cell m,p of output matrix
        output.data[m][p] = 0;
        for (std::size_t n = 0; n < N; n++) {
          output.data[m][p] += data[m][n] * other.data[n][p];
        }
      }
    }

    data = output.data;
  }

  void operator+=(Matrix<M, N, Numeric> other) {
    for (std::size_t m = 0; m < M; m++) {
      for (std::size_t n = 0; n < N; n++) {
        // Do cell m,n of output matrix
        data[m][n] += other.data[m][n];
      }
    }
  }

  void operator-=(Matrix<M, N, Numeric> other) {
    for (std::size_t m = 0; m < M; m++) {
      for (std::size_t n = 0; n < N; n++) {
        // Do cell m,n of output matrix
        data[m][n] -= other.data[m][n];
      }
    }
  }

#pragma endregion matrix_operations

#pragma region scalar_operations
  constexpr Matrix<M, N, Numeric> operator*(Numeric other) {
    Matrix<M, N, Numeric> output;

    for (std::size_t m = 0; m < M; m++) {
      for (std::size_t n = 0; n < N; n++) {
        // Do cell m,p of output matrix
        output.data[m][n] = data[m][n] * other;
      }
    }

    return output;
  }

  constexpr Matrix<M, N, Numeric> operator+(Numeric other) {
    Matrix<M, N, Numeric> output;

    for (std::size_t m = 0; m < M; m++) {
      for (std::size_t n = 0; n < N; n++) {
        // Do cell m,n of output matrix
        output.data[m][n] = data[m][n] + other;
      }
    }

    return output;
  }

  constexpr Matrix<M, N, Numeric> operator-(Numeric other) {
    Matrix<M, N, Numeric> output;

    for (std::size_t m = 0; m < M; m++) {
      for (std::size_t n = 0; n < N; n++) {
        // Do cell m,n of output matrix
        output.data[m][n] = data[m][n] - other;
      }
    }

    return output;
  }

  void operator*=(Numeric other) {

    for (std::size_t m = 0; m < M; m++) {
      for (std::size_t n = 0; n < N; n++) {
        // Do cell m,p of output matrix
        data[m][n] *= other;
      }
    }
  }

  void operator+=(Numeric other) {

    for (std::size_t m = 0; m < M; m++) {
      for (std::size_t n = 0; n < N; n++) {
        // Do cell m,n of output matrix
        data[m][n] += other;
      }
    }
  }

  void operator-=(Numeric other) {

    for (std::size_t m = 0; m < M; m++) {
      for (std::size_t n = 0; n < N; n++) {
        // Do cell m,n of output matrix
        data[m][n] -= other;
      }
    }
  }

#pragma endregion scalar_operations

#pragma region matrix_generators
  constexpr Matrix<N, N, Numeric> identity() {
    Matrix<N, N, Numeric> output;
    for (std::size_t m = 0; m < N; m++)
      for (std::size_t p = 0; p < N; p++)
        if (m == p)
          output.data[m][p] = 1;
        else
          output.data[m][p] = 0;

    return output;
  }

  template <int P> constexpr Matrix<P, P, Numeric> inverse() = delete;

#pragma endregion matrix_generators

  std::string toString() const {
    std::stringstream out;
    for (std::size_t m = 0; m < M; m++) {
      for (std::size_t n = 0; n < N; n++) {
        out << data[m][n] << '\t';
      }
      out << '\n';
    }
    return out.str();
  }
};

namespace staticGenerators {

template <int N, Matrixable Numeric = double>
constexpr Matrix<N + 1, N + 1, Numeric>
transformTranslate(const Matrix<1, N, Numeric> &v) {
  Matrix<N + 1, N + 1, Numeric> output;
  for (std::size_t m = 0; m < N + 1; m++)
    for (std::size_t p = 0; p < N + 1; p++)
      if (m == p)
        output.data[m][p] = 1;
      else if (m == N)
        output.data[m][p] = v.data[0][p];
      else
        output.data[m][p] = 0;

  return output;
}

template <int N, Matrixable Numeric = double>
constexpr Matrix<N + 1, N + 1, Numeric>
transformScale(const Matrix<1, N, Numeric> &v) {
  Matrix<N + 1, N + 1, Numeric> output;
  for (std::size_t m = 0; m < N + 1; m++)
    for (std::size_t p = 0; p < N + 1; p++)
      if (m == p) {
        if (m == N)
          output.data[m][p] = 1;
        else
          output.data[m][p] = v.data[0][p];
      } else
        output.data[m][p] = 0;

  return output;
}

template <int N, Matrixable Numeric = double>
constexpr Matrix<4, 4, Numeric>
transformRotate(const Matrix<1, 3, Numeric> &axis, double degAngle) {

  double radAngle = degAngle * M_PI / 180;
  double c = std::cos(radAngle);
  double s = std::sin(radAngle);
  double C = 1 - c;

  // Matrix<4, 4, Numeric> output
  return {
      c + axis.data[0][0] * axis.data[0][0] * C,                   //
      axis.data[0][0] * axis.data[0][1] * C - axis.data[0][2] * s, //
      axis.data[0][0] * axis.data[0][2] * C + axis.data[0][1] * s, //
      0,                                                           //
                                                                   //
      axis.data[0][1] * axis.data[0][0] * C + axis.data[0][2] * s, //
      c + axis.data[0][1] * axis.data[0][1] * C,                   //
      axis.data[0][1] * axis.data[0][2] * C - axis.data[0][0] * s, //
      0,                                                           //
                                                                   //
      axis.data[0][2] * axis.data[0][0] * C - axis.data[0][1] * s, //
      axis.data[0][2] * axis.data[0][1] * C + axis.data[0][0] * s, //
      c + axis.data[0][2] * axis.data[0][2] * C,                   //
      0,                                                           //
                                                                   //
      0,                                                           //
      0,                                                           //
      0,                                                           //
      1                                                            //
  };

  // return output;
}

}; // namespace staticGenerators

template <int M, int N, Matrixable Numeric = double>
std::ostream &operator<<(std::ostream &os, const Matrix<M, N, Numeric> &m) {
  return os << m.toString();
}

template <int N, Matrixable Numeric = double>
using Vector = Matrix<1, N + 1, Numeric>;

template <int N, Matrixable Numeric = double>
Vector<N, Numeric> createVector(const Matrix<N, 1, Numeric> &m) {
  Vector<N, Numeric> v;
  for (std::size_t i = 0; i < N; i++)
    v.data[0][i] = m.data[0][i];
  v.data[0][N] = 1; // Set Translate Variable
  return v;
}

int main() {
  //
  auto a = createVector<3>({1, 2, 3});
  auto b = staticGenerators::transformScale<3>({2, 2, 2}) *
           staticGenerators::transformRotate<3>({1, 0, 0}, 90) *
           staticGenerators::transformTranslate<3>({0, 0, 0});

  std::cout << a << '\n';
  std::cout << b << '\n';
  std::cout << a * b << '\n';
}
