#include "stdafx.h"

#include "Blas2.h"

namespace MKLSharp {
  void Blas2::sgbmv(CBlasLayout Layout, CBlasTranspose Trans,
                    long m, long n, long kl, long ku,
                    float alpha, array<float>^ a, long lda,
                    array<float>^ x, long incX, float beta, array<float>^ y, long incY) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_y = &y[0];
    cblas_sgbmv((CBLAS_LAYOUT)Layout, (CBLAS_TRANSPOSE)Trans,
                m, n, kl, ku, alpha, ptr_a, lda, ptr_x, incX, beta, ptr_y, incY);
    ptr_a = nullptr;
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
  void Blas2::dgbmv(CBlasLayout Layout, CBlasTranspose Trans,
                    long m, long n, long kl, long ku,
                    double alpha, array<double>^ a, long lda,
                    array<double>^ x, long incX, double beta, array<double>^ y, long incY) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_y = &y[0];
    cblas_dgbmv((CBLAS_LAYOUT)Layout, (CBLAS_TRANSPOSE)Trans,
                m, n, kl, ku, alpha, ptr_a, lda, ptr_x, incX, beta, ptr_y, incY);
    ptr_a = nullptr;
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
}