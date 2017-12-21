#include "stdafx.h"

#include "Blas1.h"

namespace MKLSharp {
  float Blas1::sasum(long n, array<float>^ x, long incX) {
    pin_ptr<float> ptr_x = &x[0];
    auto res = cblas_sasum(n, ptr_x, incX);
    ptr_x = nullptr;
    return res;
  }
  double Blas1::dasum(long n, array<double>^ x, long incX) {
    pin_ptr<double> ptr_x = &x[0];
    auto res = cblas_dasum(n, ptr_x, incX);
    ptr_x = nullptr;
    return res;
  }

  void Blas1::saxpy(long n, float a, array<float>^ x, long incX, array<float>^ y, long incY) {
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_y = &y[0];
    cblas_saxpy(n, a, ptr_x, incX, ptr_y, incY);
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
  void Blas1::daxpy(long n, double a, array<double>^ x, long incX, array<double>^ y, long incY) {
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_y = &y[0];
    cblas_daxpy(n, a, ptr_x, incX, ptr_y, incY);
    ptr_x = nullptr;
    ptr_y = nullptr;
  }

  void Blas1::scopy(long n, array<float>^ x, long incX, array<float>^ y, long incY) {
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_y = &y[0];
    cblas_scopy(n, ptr_x, incX, ptr_y, incY);
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
  void Blas1::scopy(long n, array<float>^ x, long incX, [Out]array<float>^% y, long incY) {
    pin_ptr<float> ptr_x = &x[0];
    y = gcnew array<float>(n);
    pin_ptr<float> ptr_y = &y[0];
    cblas_scopy(n, ptr_x, incX, ptr_y, incY);
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
  void Blas1::dcopy(long n, array<double>^ x, long incX, array<double>^ y, long incY) {
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_y = &y[0];
    cblas_dcopy(n, ptr_x, incX, ptr_y, incY);
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
  void Blas1::dcopy(long n, array<double>^ x, long incX, [Out]array<double>^% y, long incY) {
    pin_ptr<double> ptr_x = &x[0];
    y = gcnew array<double>(n);
    pin_ptr<double> ptr_y = &y[0];
    cblas_dcopy(n, ptr_x, incX, ptr_y, incY);
    ptr_x = nullptr;
    ptr_y = nullptr;
  }

  float Blas1::sdot(long n, array<float>^ x, long incX, array<float>^ y, long incY) {
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_y = &y[0];
    auto res = cblas_sdot(n, ptr_x, incX, ptr_y, incY);
    ptr_x = nullptr;
    ptr_y = nullptr;
    return res;
  }
  double Blas1::ddot(long n, array<double>^ x, long incX, array<double>^ y, long incY) {
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_y = &y[0];
    auto res = cblas_ddot(n, ptr_x, incX, ptr_y, incY);
    ptr_x = nullptr;
    ptr_y = nullptr;
    return res;
  }

  float Blas1::sdsdot(long n, float sb, array<float>^ sx, long incX, array<float>^ sy, long incY) {
    pin_ptr<float> ptr_sx = &sx[0];
    pin_ptr<float> ptr_sy = &sy[0];
    auto res = cblas_sdsdot(n, sb, ptr_sx, incX, ptr_sy, incY);
    ptr_sx = nullptr;
    ptr_sy = nullptr;
    return res;
  }
  double Blas1::dsdot(long n, array<float>^ sx, long incX, array<float>^ sy, long incY) {
    pin_ptr<float> ptr_sx = &sx[0];
    pin_ptr<float> ptr_sy = &sy[0];
    auto res = cblas_dsdot(n, ptr_sx, incX, ptr_sy, incY);
    ptr_sx = nullptr;
    ptr_sy = nullptr;
    return res;
  }

  float Blas1::snrm2(long n, array<float>^ x, long incX) {
    pin_ptr<float> ptr_x = &x[0];
    auto res = cblas_snrm2(n, ptr_x, incX);
    ptr_x = nullptr;
    return res;
  }
  double Blas1::dnrm2(long n, array<double>^ x, long incX) {
    pin_ptr<double> ptr_x = &x[0];
    auto res = cblas_dnrm2(n, ptr_x, incX);
    ptr_x = nullptr;
    return res;
  }

  void Blas1::srot(long n, array<float>^ x, long incX, array<float>^ y, long incY, float c, float s) {
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_y = &y[0];
    cblas_srot(n, ptr_x, incX, ptr_y, incY, c, s);
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
  void Blas1::drot(long n, array<double>^ x, long incX, array<double>^ y, long incY, double c, double s) {
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_y = &y[0];
    cblas_drot(n, ptr_x, incX, ptr_y, incY, c, s);
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
}