#pragma once

using namespace System;

namespace MKLSharp {
  public ref class Blas1 {
  public:
    static float sasum(int n, array<float>^ x, int incX);
    static double dasum(int n, array<double>^ x, int incX);

    static void saxpy(int n, float a, array<float>^ x, int incX, array<float>^ y, int incY);
    static void daxpy(int n, double a, array<double>^ x, int incX, array<double>^ y, int incY);
  };
}