#pragma once

using namespace System::Runtime::InteropServices;

namespace MKLSharp {
  public ref class Blas1 {
  public:
    static float sasum(int n, array<float>^ x, int incX);
    static double dasum(int n, array<double>^ x, int incX);

    static void saxpy(int n, float a, array<float>^ x, int incX, array<float>^ y, int incY);
    static void daxpy(int n, double a, array<double>^ x, int incX, array<double>^ y, int incY);

    static void scopy(int n, array<float>^ x, int incX, array<float>^ y, int incY);
    static void scopy(int n, array<float>^ x, int incX, [Out]array<float>^% y, int incY);
    static void dcopy(int n, array<double>^ x, int incX, array<double>^ y, int incY);
    static void dcopy(int n, array<double>^ x, int incX, [Out]array<double>^% y, int incY);

    static float sdot(int n, array<float>^ x, int incX, array<float>^ y, int incY);
    static double ddot(int n, array<double>^ x, int incX, array<double>^ y, int incY);

    static float sdsdot(int n, float sb, array<float>^ sx, int incX, array<float>^ sy, int incY);
    static double dsdot(int n, array<float>^ sx, int incX, array<float>^ sy, int incY);
  };
}