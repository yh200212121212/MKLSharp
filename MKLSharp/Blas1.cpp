#include "stdafx.h"

#include "Blas1.h"

namespace MKLSharp {
  float Blas1::sasum(int n, array<float>^ x, int incX) {
    pin_ptr<float> ptr_x = &x[0];
    auto res = cblas_sasum(n, ptr_x, incX);
    ptr_x = nullptr;
    return res;
  }
}