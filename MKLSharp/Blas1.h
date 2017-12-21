#pragma once

using namespace System::Runtime::InteropServices;

namespace MKLSharp {
  public ref class Blas1 {
  public:
    static float sasum(long n, array<float>^ x, long incX);
    static double dasum(long n, array<double>^ x, long incX);

    static void saxpy(long n, float a, array<float>^ x, long incX, array<float>^ y, long incY);
    static void daxpy(long n, double a, array<double>^ x, long incX, array<double>^ y, long incY);

    static void scopy(long n, array<float>^ x, long incX, array<float>^ y, long incY);
    static void scopy(long n, array<float>^ x, long incX, [Out]array<float>^% y, long incY);
    static void dcopy(long n, array<double>^ x, long incX, array<double>^ y, long incY);
    static void dcopy(long n, array<double>^ x, long incX, [Out]array<double>^% y, long incY);

    static float sdot(long n, array<float>^ x, long incX, array<float>^ y, long incY);
    static double ddot(long n, array<double>^ x, long incX, array<double>^ y, long incY);

    static float sdsdot(long n, float sb, array<float>^ sx, long incX, array<float>^ sy, long incY);
    static double dsdot(long n, array<float>^ sx, long incX, array<float>^ sy, long incY);

    static float snrm2(long n, array<float>^ x, long incX);
    static double dnrm2(long n, array<double>^ x, long incX);

    static void srot(long n, array<float>^ x, long incX, array<float>^ y, long incY, float c, float s);
    static void drot(long n, array<double>^ x, long incX, array<double>^ y, long incY, double c, double s);

    static void srotg(float% a, float% b, [Out]float% c, [Out]float% s);
    static void drotg(double% a, double% b, [Out]double% c, [Out]double% s);

    static void srotm(long n, array<float>^ x, long incX, array<float>^ y, long incY, array<float>^ param);
    static void drotm(long n, array<double>^ x, long incX, array<double>^ y, long incY, array<double>^ param);

    static void srotmg(float% d1, float% d2, float% x1, float y1, [Out]array<float>^% param);
    static void drotmg(double% d1, double% d2, double% x1, double y1, [Out]array<double>^% param);

    static void sscal(long n, float a, array<float>^ x, long incX);
    static void dscal(long n, double a, array<double>^ x, long incX);
  };
}