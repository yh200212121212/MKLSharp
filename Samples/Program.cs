using MKLSharp;
using static System.Console;

namespace Samples {
  class Program {
    static void Main(string[] args) {
      var x = new float[] { 1.0f, 1.0f, 1.0f };
      WriteLine("Level1 BLAS sasum call test.");
      WriteLine(Blas1.sasum(x.Length, x, 1));

      WriteLine("Please press Enter key...");
      ReadLine();
    }
  }
}
