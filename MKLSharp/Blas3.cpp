#include "stdafx.h"

#include "Blas3.h"

namespace MKLSharp {
  void Blas3::sgemm(CBlasLayout Layout, CBlasTranspose TransA, CBlasTranspose TransB,
                    long m, long n, long k,
                    float alpha, array<float>^ a, long lda,
                    array<float>^ b, long ldb,
                    float beta, array<float>^ c, long ldc) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_c = &c[0];
    cblas_sgemm((CBLAS_LAYOUT)Layout, (CBLAS_TRANSPOSE)TransA, (CBLAS_TRANSPOSE)TransB,
                m, n, k, alpha, ptr_a, lda, ptr_b, ldb, beta, ptr_c, ldc);
    ptr_a = nullptr;
    ptr_b = nullptr;
    ptr_c = nullptr;
  }
  void Blas3::dgemm(CBlasLayout Layout, CBlasTranspose TransA, CBlasTranspose TransB,
                    long m, long n, long k,
                    double alpha, array<double>^ a, long lda,
                    array<double>^ b, long ldb,
                    double beta, array<double>^ c, long ldc) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_c = &c[0];
    cblas_dgemm((CBLAS_LAYOUT)Layout, (CBLAS_TRANSPOSE)TransA, (CBLAS_TRANSPOSE)TransB,
                m, n, k, alpha, ptr_a, lda, ptr_b, ldb, beta, ptr_c, ldc);
    ptr_a = nullptr;
    ptr_b = nullptr;
    ptr_c = nullptr;
  }

  void Blas3::ssymm(CBlasLayout Layout, CBlasSide Side, CBlasUpLo UpLo,
                    long m, long n, float alpha, array<float>^ a, long lda,
                    array<float>^ b, long ldb,
                    float beta, array<float>^ c, long ldc) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_b = &b[0];
    pin_ptr<float> ptr_c = &c[0];
    cblas_ssymm((CBLAS_LAYOUT)Layout, (CBLAS_SIDE)Side, (CBLAS_UPLO)UpLo,
                m, n, alpha, ptr_a, lda, ptr_b, ldb, beta, ptr_c, ldc);
    ptr_a = nullptr;
    ptr_b = nullptr;
    ptr_c = nullptr;
  }
  void Blas3::dsymm(CBlasLayout Layout, CBlasSide Side, CBlasUpLo UpLo,
                    long m, long n, double alpha, array<double>^ a, long lda,
                    array<double>^ b, long ldb,
                    double beta, array<double>^ c, long ldc) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_b = &b[0];
    pin_ptr<double> ptr_c = &c[0];
    cblas_dsymm((CBLAS_LAYOUT)Layout, (CBLAS_SIDE)Side, (CBLAS_UPLO)UpLo,
                m, n, alpha, ptr_a, lda, ptr_b, ldb, beta, ptr_c, ldc);
    ptr_a = nullptr;
    ptr_b = nullptr;
    ptr_c = nullptr;
  }

  void Blas3::ssyrk(CBlasLayout Layout, CBlasUpLo UpLo, CBlasTranspose Trans,
                    long n, long k, float alpha, array<float>^ a, long lda,
                    float beta, array<float>^ c, long ldc) {
    pin_ptr<float> ptr_a = &a[0];
    pin_ptr<float> ptr_c = &c[0];
    cblas_ssyrk((CBLAS_LAYOUT)Layout, (CBLAS_UPLO)UpLo, (CBLAS_TRANSPOSE)Trans,
                n, k, alpha, ptr_a, lda, beta, ptr_c, ldc);
    ptr_a = nullptr;
    ptr_c = nullptr;
  }
  void Blas3::dsyrk(CBlasLayout Layout, CBlasUpLo UpLo, CBlasTranspose Trans,
                    long n, long k, double alpha, array<double>^ a, long lda,
                    double beta, array<double>^ c, long ldc) {
    pin_ptr<double> ptr_a = &a[0];
    pin_ptr<double> ptr_c = &c[0];
    cblas_dsyrk((CBLAS_LAYOUT)Layout, (CBLAS_UPLO)UpLo, (CBLAS_TRANSPOSE)Trans,
                n, k, alpha, ptr_a, lda, beta, ptr_c, ldc);
    ptr_a = nullptr;
    ptr_c = nullptr;
  }
}