using MKLSharp;
using static System.Console;

namespace Samples {
  class Program {
    static void Main(string[] args) {
      var x = new float[] { 1.0f, 1.0f, 1.0f };
      WriteLine("Level1 BLAS sasum call test.");
      WriteLine(Blas1.sasum(x.Length, x, 1));

      WriteLine("Level1 BLAS scopy call test.");
      Blas1.scopy(x.Length, x, 1, out var y, 1);
      for (var i = 0; i < y.Length; i++)
        Write(y[i] + " ");
      WriteLine("\n");

      WriteLine("Level2 BLAS dgemv call test.");
      var ad = new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
      var xd = new double[] { 1.0, 1.0, 1.0 };
      var yd = new double[] { 0.0, 0.0, 0.0 };
      Blas2.dgemv(CBlasLayout.RowMajor, CBlasTranspose.NoTrans, 3, 3, 1.0, ad, 3, xd, 1, 1.0, yd, 1);
      for (var i = 0; i < yd.Length; i++)
        Write(yd[i] + " ");
      WriteLine("\n");

      WriteLine("LAPACK General Matrix call test.");
      var ag = new double[] { 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0 };
      var bg = new double[] { 6.0, 7.0, 12.0, 15.0 };
      Lapack.dgetrf(LapackLayout.RowMajor, 4, 4, ag, 4, out var ipiv);
      Lapack.dgetrs(LapackLayout.RowMajor, (sbyte)'N', 4, 1, ag, 4, ipiv, bg, 1);
      for (var i = 0; i < bg.Length; i++)
        Write(bg[i] + " ");
      WriteLine();

      WriteLine("Please press Enter key...");
      ReadLine();
    }
  }
}
