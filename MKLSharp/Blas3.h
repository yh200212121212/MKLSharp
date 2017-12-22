#pragma once

using namespace System::Runtime::InteropServices;

namespace MKLSharp {
  public ref class Blas3 {
  public:
    static void sgemm(CBlasLayout Layout, CBlasTranspose TransA, CBlasTranspose TransB,
                      long m, long n, long k,
                      float alpha, array<float>^ a, long lda,
                      array<float>^ b, long ldb,
                      float beta, array<float>^ c, long ldc);
    static void dgemm(CBlasLayout Layout, CBlasTranspose TransA, CBlasTranspose TransB,
                      long m, long n, long k,
                      double alpha, array<double>^ a, long lda,
                      array<double>^ b, long ldb,
                      double beta, array<double>^ c, long ldc);

    static void ssymm(CBlasLayout Layout, CBlasSide Side, CBlasUpLo UpLo,
                      long m, long n, float alpha, array<float>^ a, long lda,
                      array<float>^ b, long ldb,
                      float beta, array<float>^ c, long ldc);
    static void dsymm(CBlasLayout Layout, CBlasSide Side, CBlasUpLo UpLo,
                      long m, long n, double alpha, array<double>^ a, long lda,
                      array<double>^ b, long ldb,
                      double beta, array<double>^ c, long ldc);

    static void ssyrk(CBlasLayout Layout, CBlasUpLo UpLo, CBlasTranspose Trans,
                      long n, long k, float alpha, array<float>^ a, long lda,
                      float beta, array<float>^ c, long ldc);
    static void dsyrk(CBlasLayout Layout, CBlasUpLo UpLo, CBlasTranspose Trans,
                      long n, long k, double alpha, array<double>^ a, long lda,
                      double beta, array<double>^ c, long ldc);

    static void ssyr2k(CBlasLayout Layout, CBlasUpLo UpLo, CBlasTranspose Trans,
                       long n, long k, float alpha, array<float>^ a, long lda,
                       array<float>^ b, long ldb,
                       float beta, array<float>^ c, long ldc);
    static void dsyr2k(CBlasLayout Layout, CBlasUpLo UpLo, CBlasTranspose Trans,
                       long n, long k, double alpha, array<double>^ a, long lda,
                       array<double>^ b, long ldb,
                       double beta, array<double>^ c, long ldc);

    static void strmm(CBlasLayout Layout, CBlasSide Side, CBlasUpLo UpLo,
                      CBlasTranspose TransA, CBlasDiag Diag, long m, long n,
                      float alpha, array<float>^ a, long lda,
                      array<float>^ b, long ldb);
    static void dtrmm(CBlasLayout Layout, CBlasSide Side, CBlasUpLo UpLo,
                      CBlasTranspose TransA, CBlasDiag Diag, long m, long n,
                      double alpha, array<double>^ a, long lda,
                      array<double>^ b, long ldb);

    static void strsm(CBlasLayout Layout, CBlasSide Side, CBlasUpLo UpLo,
                      CBlasTranspose TransA, CBlasDiag Diag, long m, long n,
                      float alpha, array<float>^ a, long lda,
                      array<float>^ b, long ldb);
    static void dtrsm(CBlasLayout Layout, CBlasSide Side, CBlasUpLo UpLo,
                      CBlasTranspose TransA, CBlasDiag Diag, long m, long n,
                      double alpha, array<double>^ a, long lda,
                      array<double>^ b, long ldb);
  };
}