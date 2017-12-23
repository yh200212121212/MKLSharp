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
  public enum class LapackLayout {
    RowMajor = LAPACK_ROW_MAJOR,
    ColMajor = LAPACK_COL_MAJOR
  };
  public enum class LapackTranspose {
    N = 'N',
    T = 'T',
    C = 'C'
  };
  public enum class LapackNorm {
    One = '1',
    O = 'O',
    I = 'I'
  };
  public enum class LapackEquil {
    N = 'N',
    R = 'R',
    C = 'C',
    B = 'B'
  };
  public enum class LapackUpLo {
    U = 'U',
    L = 'L'
  };
}