#pragma once

#include <mkl.h>

namespace MKLSharp {
  public enum class CBlasLayout { RowMajor = CblasRowMajor, ColMajor = CblasColMajor };
  public enum class CBlasTranspose { NoTrans = CblasNoTrans, Trans = CblasTrans, ConjTrans = CblasConjTrans };
}