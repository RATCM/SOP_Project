using System.Runtime.InteropServices;
using System.Diagnostics;
using LinAlg.Matrix;
using LinAlg.Complex;
using System.Windows.Markup;
using BenchmarkDotNet;
using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Order;
using BenchmarkDotNet.Running;

namespace MatrixCalc
{
    [MemoryDiagnoser]
    [Orderer(SummaryOrderPolicy.FastestToSlowest)]
    [RankColumn]
    public class MatrixMulBenchmark
    {
        static float[][] InitJaggedArray(int length)
        {
            var f = new float[length][];

            for(int i = 0; i < length; i++)
            {
                f[i] = new float[length];
            }

            return f;
        }

        private Matrix matrix = new Matrix(5,5);
        private float[,] f1 = new float[100,100];
        private float[][] f2 = InitJaggedArray(100);

        [Benchmark]
        public void NewMatrix1()
        {
            new Matrix(f1);
        }

        [Benchmark]
        public void NewMatrix2()
        {
            new Matrix(f2);
        }
    }
    class Program
    {
        static void Main()
        {
            BenchmarkRunner.Run<MatrixMulBenchmark>();
        }
    }
}