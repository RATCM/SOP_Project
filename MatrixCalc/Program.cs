using System.Runtime.InteropServices;
using System.Diagnostics;
using LinAlg.Matrix;
using LinAlg.Complex;
using System.Windows.Markup;

namespace MatrixCalc
{
    class Program
    {
        static void Main()
        {
            Random rand = new Random();

            //int dim = rand.Next(10, 1000);
            int dim = 4;
            float[,] a = new float[dim, dim];
            float[,] b = new float[dim, dim];

            Console.WriteLine($"Generating matricies with dimensions {dim}x{dim}");
            for (int i = 0; i < dim; i++)
            {
                for (int j = 0; j < dim; j++)
                {
                    a[i, j] = (float)(rand.Next(1, 1000) / 100f);
                    b[i, j] = (float)(rand.Next(1, 1000) / 100f);
                }
            }
            Console.WriteLine();

            Matrix m1 = new Matrix(a);
            Matrix m2 = new Matrix(b);


            //Matrix m = new Matrix(new float[,] { { 8.16f, 9.11f }, { 2.58f, 1.5f } });


            Stopwatch stopwatch = new Stopwatch();

            Console.WriteLine("Calculating matrix 1,000,000 multiplications...");
            stopwatch.Start();
            for (int i = 0; i < 1_000_000; i++)
            {
                Matrix m3 = m1 * m2;
            }
            stopwatch.Stop();
            Console.WriteLine($"Elapsed time: {stopwatch.ElapsedMilliseconds} ms");
            //stopwatch = new Stopwatch();

            //Console.WriteLine();


            //Console.WriteLine($"Calculating inverse of the matrix");
            //stopwatch.Start();
            //Matrix m3_inv = m3.Inverse;
            //stopwatch.Stop();
            //Console.WriteLine($"Elapsed time: {stopwatch.ElapsedMilliseconds} ms");

            //Console.WriteLine();

            //stopwatch = new Stopwatch();

            //Console.WriteLine($"Calculating determinant of the matrix");
            //stopwatch.Start();
            //var det = m3.Determinant;
            //stopwatch.Stop();
            //Console.WriteLine($"Elapsed time: {stopwatch.ElapsedMilliseconds} ms");

            //Console.WriteLine();

            Console.WriteLine("Done!");

            //Console.WriteLine("Matrix:");
            //Console.WriteLine(m3);

            //Console.WriteLine("Eigenvalues:");
            //foreach (var value in values)
            //{
            //    Console.WriteLine(value);
            //    Console.WriteLine();
            //}



            //Console.WriteLine("Eigenvectors");


            //var vectors = m.eigenvectors;

            //for (int i = 0; i < vectors.Length; i++)
            //{
            //    Console.WriteLine(vectors[i]);
            //}

            //Console.WriteLine();
            //Console.WriteLine("Check eigenvalues");

            //for (int i = 0; i < vectors.Length; i++)
            //{
            //    Console.WriteLine((m - Matrix.Identity(m.rowLength) * values[i]).Determinant);
            //}
            //var newVec = (m * vector).normalized - vector.normalized;
        }
    }
}