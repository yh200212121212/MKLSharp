#pragma once

using namespace System::Runtime::InteropServices;

namespace MKLSharp {
  public ref class Blas2 {
    static void sgbmv(CBlasLayout Layout, CBlasTranspose Trans,
                      long m, long n, long kl, long ku,
                      float alpha, array<float>^ a, long lda,
                      array<float>^ x, long incX, float beta, array<float>^ y, long incY);
    static void dgbmv(CBlasLayout Layout, CBlasTranspose Trans,
                      long m, long n, long kl, long ku,
                      double alpha, array<double>^ a, long lda,
                      array<double>^ x, long incX, double beta, array<double>^ y, long incY);

    static void sgemv(CBlasLayout Layout, CBlasTranspose Trans,
                      long m, long n, float alpha, array<float>^ a, long lda,
                      array<float>^ x, long incX, float beta, array<float>^ y, long incY);
    static void dgemv(CBlasLayout Layout, CBlasTranspose Trans,
                      long m, long n, double alpha, array<double>^ a, long lda,
                      array<double>^ x, long incX, double beta, array<double>^ y, long incY);
  };
}