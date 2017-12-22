#pragma once

#include <mkl.h>

namespace MKLSharp {
  public enum class CBlasLayout { 
    RowMajor = CblasRowMajor,
    ColMajor = CblasColMajor
  };
  public enum class CBlasTranspose {
    NoTrans = CblasNoTrans,
    Trans = CblasTrans,
    ConjTrans = CblasConjTrans
  };
  public enum class CBlasUpLo {
    Upper = CblasUpper,
    Lower = CblasLower
  };
  public enum class CBlasDiag {
    Unit = CblasUnit,
    NonUnit = CblasNonUnit
  };
  public enum class CBlasSide {
    Left = CblasLeft,
    Right = CblasRight
  };
}