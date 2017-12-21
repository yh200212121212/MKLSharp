#pragma once

using namespace System;

namespace MKLSharp {

  public ref class Blas1 {
  public:
    static float sasum(int n, array<float>^ x, int incX);
  };
}