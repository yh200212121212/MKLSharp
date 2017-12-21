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

  void Blas1::srotg(float% a, float% b, [Out]float% c, [Out]float% s) {
    pin_ptr<float> ptr_a = &a;
    pin_ptr<float> ptr_b = &b;
    pin_ptr<float> ptr_c = &c;
    pin_ptr<float> ptr_s = &s;
    cblas_srotg(ptr_a, ptr_b, ptr_c, ptr_s);
    ptr_a = nullptr;
    ptr_b = nullptr;
    ptr_c = nullptr;
    ptr_s = nullptr;
  }
  void Blas1::drotg(double% a, double% b, [Out]double% c, [Out]double% s) {
    pin_ptr<double> ptr_a = &a;
    pin_ptr<double> ptr_b = &b;
    pin_ptr<double> ptr_c = &c;
    pin_ptr<double> ptr_s = &s;
    cblas_drotg(ptr_a, ptr_b, ptr_c, ptr_s);
    ptr_a = nullptr;
    ptr_b = nullptr;
    ptr_c = nullptr;
    ptr_s = nullptr;
  }

  void Blas1::srotm(long n, array<float>^ x, long incX, array<float>^ y, long incY, array<float>^ param) {
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_y = &y[0];
    pin_ptr<float> ptr_p = &param[0];
    cblas_srotm(n, ptr_x, incX, ptr_y, incY, ptr_p);
    ptr_x = nullptr;
    ptr_y = nullptr;
    ptr_p = nullptr;
  }
  void Blas1::drotm(long n, array<double>^ x, long incX, array<double>^ y, long incY, array<double>^ param) {
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_y = &y[0];
    pin_ptr<double> ptr_p = &param[0];
    cblas_drotm(n, ptr_x, incX, ptr_y, incY, ptr_p);
    ptr_x = nullptr;
    ptr_y = nullptr;
    ptr_p = nullptr;
  }

  void Blas1::srotmg(float% d1, float% d2, float% x1, float y1, [Out]array<float>^% param) {
    pin_ptr<float> ptr_d1 = &d1;
    pin_ptr<float> ptr_d2 = &d2;
    pin_ptr<float> ptr_x1 = &x1;
    param = gcnew array<float>(5);
    pin_ptr<float> ptr_p = &param[0];
    cblas_srotmg(ptr_d1, ptr_d2, ptr_x1, y1, ptr_p);
    ptr_d1 = nullptr;
    ptr_d2 = nullptr;
    ptr_x1 = nullptr;
    ptr_p = nullptr;
  }
  void Blas1::drotmg(double% d1, double% d2, double% x1, double y1, [Out]array<double>^% param) {
    pin_ptr<double> ptr_d1 = &d1;
    pin_ptr<double> ptr_d2 = &d2;
    pin_ptr<double> ptr_x1 = &x1;
    param = gcnew array<double>(5);
    pin_ptr<double> ptr_p = &param[0];
    cblas_drotmg(ptr_d1, ptr_d2, ptr_x1, y1, ptr_p);
    ptr_d1 = nullptr;
    ptr_d2 = nullptr;
    ptr_x1 = nullptr;
    ptr_p = nullptr;
  }

  void Blas1::sscal(long n, float a, array<float>^ x, long incX) {
    pin_ptr<float> ptr_x = &x[0];
    cblas_sscal(n, a, ptr_x, incX);
    ptr_x = nullptr;
  }
  void Blas1::dscal(long n, double a, array<double>^ x, long incX) {
    pin_ptr<double> ptr_x = &x[0];
    cblas_dscal(n, a, ptr_x, incX);
    ptr_x = nullptr;
  }

  void Blas1::sswap(long n, array<float>^ x, long incX, array<float>^ y, long incY) {
    pin_ptr<float> ptr_x = &x[0];
    pin_ptr<float> ptr_y = &y[0];
    cblas_sswap(n, ptr_x, incX, ptr_y, incY);
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
  void Blas1::dswap(long n, array<double>^ x, long incX, array<double>^ y, long incY) {
    pin_ptr<double> ptr_x = &x[0];
    pin_ptr<double> ptr_y = &y[0];
    cblas_dswap(n, ptr_x, incX, ptr_y, incY);
    ptr_x = nullptr;
    ptr_y = nullptr;
  }
}