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

  void Blas2::sgemv(CBlasLayout Layout, CBlasTranspose Trans,
                    long m, long n, float alpha, array<float>^ a, long lda,
                    array<float>^ x, long incX, float beta, array<float>^ y, long incY) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_y = &y[0];
    cblas_sgemv((CBLAS_LAYOUT)Layout, (CBLAS_TRANSPOSE)Trans,
                m, n, alpha, ptr_a, lda, ptr_x, incX, beta, ptr_y, incY);
    ptr_a = nullptr;
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
  void Blas2::dgemv(CBlasLayout Layout, CBlasTranspose Trans,
                    long m, long n, double alpha, array<double>^ a, long lda,
                    array<double>^ x, long incX, double beta, array<double>^ y, long incY) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_y = &y[0];
    cblas_dgemv((CBLAS_LAYOUT)Layout, (CBLAS_TRANSPOSE)Trans,
                m, n, alpha, ptr_a, lda, ptr_x, incX, beta, ptr_y, incY);
    ptr_a = nullptr;
    ptr_x = nullptr;
    ptr_y = nullptr;
  }

  void Blas2::sger(CBlasLayout Layout, long m, long n,
                   float alpha, array<float>^ x, long incX,
                   array<float>^ y, long incY,
                   array<float>^ a, long lda) {
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_y = &y[0];
    pin_ptr<float> ptr_a = &a[0];
    cblas_sger((CBLAS_LAYOUT)Layout, m, n, alpha, ptr_x, incX, ptr_y, incY, ptr_a, lda);
    ptr_x = nullptr;
    ptr_y = nullptr;
    ptr_a = nullptr;
  }
  void Blas2::dger(CBlasLayout Layout, long m, long n,
                   double alpha, array<double>^ x, long incX,
                   array<double>^ y, long incY,
                   array<double>^ a, long lda) {
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_y = &y[0];
    pin_ptr<double> ptr_a = &a[0];
    cblas_dger((CBLAS_LAYOUT)Layout, m, n, alpha, ptr_x, incX, ptr_y, incY, ptr_a, lda);
    ptr_x = nullptr;
    ptr_y = nullptr;
    ptr_a = nullptr;
  }

  void Blas2::ssbmv(CBlasLayout Layout, CBlasUpLo UpLo,
                    long n, long k, float alpha, array<float>^ a, long lda,
                    array<float>^ x, long incX,
                    float beta, array<float>^ y, long incY) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_y = &y[0];
    cblas_ssbmv((CBLAS_LAYOUT)Layout, (CBLAS_UPLO)UpLo,
                n, k, alpha, ptr_a, lda, ptr_x, incX, beta, ptr_y, incY);
    ptr_a = nullptr;
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
  void Blas2::dsbmv(CBlasLayout Layout, CBlasUpLo UpLo,
                    long n, long k, double alpha, array<double>^ a, long lda,
                    array<double>^ x, long incX,
                    double beta, array<double>^ y, long incY) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_y = &y[0];
    cblas_dsbmv((CBLAS_LAYOUT)Layout, (CBLAS_UPLO)UpLo,
                n, k, alpha, ptr_a, lda, ptr_x, incX, beta, ptr_y, incY);
    ptr_a = nullptr;
    ptr_x = nullptr;
    ptr_y = nullptr;
  }

  void Blas2::sspmv(CBlasLayout Layout, CBlasUpLo UpLo,
                    long n, float alpha, array<float>^ ap,
                    array<float>^ x, long incX,
                    float beta, array<float>^ y, long incY) {
    pin_ptr<float> ptr_a = &ap[0];
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_y = &y[0];
    cblas_sspmv((CBLAS_LAYOUT)Layout, (CBLAS_UPLO)UpLo,
                n, alpha, ptr_a, ptr_x, incX, beta, ptr_y, incY);
    ptr_a = nullptr;
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
  void Blas2::dspmv(CBlasLayout Layout, CBlasUpLo UpLo,
                    long n, double alpha, array<double>^ ap,
                    array<double>^ x, long incX,
                    double beta, array<double>^ y, long incY) {
    pin_ptr<double> ptr_a = &ap[0];
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_y = &y[0];
    cblas_dspmv((CBLAS_LAYOUT)Layout, (CBLAS_UPLO)UpLo,
                n, alpha, ptr_a, ptr_x, incX, beta, ptr_y, incY);
    ptr_a = nullptr;
    ptr_x = nullptr;
    ptr_y = nullptr;
  }

  void Blas2::sspr(CBlasLayout Layout, CBlasUpLo UpLo,
                   long n, float alpha, array<float>^ x, long incX, array<float>^ ap) {
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_a = &ap[0];
    cblas_sspr((CBLAS_LAYOUT)Layout, (CBLAS_UPLO)UpLo, n, alpha, ptr_x, incX, ptr_a);
    ptr_x = nullptr;
    ptr_a = nullptr;
  }
  void Blas2::dspr(CBlasLayout Layout, CBlasUpLo UpLo,
                   long n, double alpha, array<double>^ x, long incX, array<double>^ ap) {
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_a = &ap[0];
    cblas_dspr((CBLAS_LAYOUT)Layout, (CBLAS_UPLO)UpLo, n, alpha, ptr_x, incX, ptr_a);
    ptr_x = nullptr;
    ptr_a = nullptr;
  }
}