using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Order;
using LinAlg.Matricies;

namespace MatrixCalc
{
    [MemoryDiagnoser]
    [Orderer(SummaryOrderPolicy.FastestToSlowest)]
    [RankColumn]
    public class MatrixMulBenchmark
    {
        private Matrix matrix = new Matrix(100, 100);

        [Benchmark]
        public void ParallelMultiply()
        {
            Matrix.ParallelMultiply(matrix, matrix);
        }

        [Benchmark]
        public void RegularMultiply()
        {
            Matrix.RegularMultiply(matrix, matrix);
        }
    }
    class Program
    {
        static void Main()
        {
            Matrix m1 = new float[,]
            {
                { 5, 2, 6 },
                { 3, 6, 3 },
                { 1, 8, -2 }
            };

            Vector v1 = new float[]
            { 4, 2, 8 };

            Console.WriteLine(m1.Inverse * v1);
        }
    }
    //BenchmarkRunner.Run<MatrixMulBenchmark>();
}