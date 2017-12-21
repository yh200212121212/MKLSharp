#include "stdafx.h"

#include "Blas1.h"

namespace MKLSharp {
  float Blas1::sasum(int n, array<float>^ x, int incX) {
    pin_ptr<float> ptr_x = &x[0];
    auto res = cblas_sasum(n, ptr_x, incX);
    ptr_x = nullptr;
    return res;
  }
  double Blas1::dasum(int n, array<double>^ x, int incX) {
    pin_ptr<double> ptr_x = &x[0];
    auto res = cblas_dasum(n, ptr_x, incX);
    ptr_x = nullptr;
    return res;
  }

  void Blas1::saxpy(int n, float a, array<float>^ x, int incX, array<float>^ y, int incY) {
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_y = &y[0];
    cblas_saxpy(n, a, ptr_x, incX, ptr_y, incY);
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
  void Blas1::daxpy(int n, double a, array<double>^ x, int incX, array<double>^ y, int incY) {
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_y = &y[0];
    cblas_daxpy(n, a, ptr_x, incX, ptr_y, incY);
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
}