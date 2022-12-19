using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;
using LinAlg.Complex;

namespace LinAlg.Matrix
{
    public class Matrix
    {

        private (Matrix, Matrix) QRDecomposition()
        {
            Matrix[] H = new Matrix[rowLength + 1];
            Matrix[] A = new Matrix[rowLength + 1];
            A[0] = new Matrix(this.rows);

            // The output matricies
            Matrix Q = Matrix.Identity(rowLength);
            Matrix R = Matrix.Identity(rowLength);

            // Go through rowLength iterations
            for (int i = 0; i < rowLength; i++)
            {
                // Dont worry about it
                Matrix choiceOfMatrix = Matrix.Identity(rowLength - i);

                for (int j = i; j < rowLength; j++)
                {
                    for (int k = i; k < colLength; k++)
                    {
                        choiceOfMatrix[j - i, k - i] = A[i][j, k];
                    }
                }
    

                Vector vec = new Vector(rowLength - i);

                // Just a column vector with the first value as 1
                Vector colVec = new float[rowLength - i];
                colVec[0] = new Complex.Complex(1,0);

                float maxNum = 0;
                for (int j = 0; j < rowLength - i; j++)
                {
                    vec[j] = choiceOfMatrix[j, 0];
                    if (maxNum < MathF.Abs(vec[j].Real))
                        maxNum = MathF.Abs(vec[j].Real);
                }
                float sign = vec[0].Real >= 0 ? 1 : -1;

                var u = vec + sign * vec.norm * colVec;
                var v = 1f / u[0].Real * u;

                float multiplier = 2 / ((v * v).Real);

                var vecAsMatrix = v.ToMatrix();

                var HMatrix = Matrix.Identity(rowLength - i) - multiplier * vecAsMatrix * vecAsMatrix.Transposed;

                Vector HOnVec = -sign * vec.norm * colVec;

                var subMatrixA = HMatrix * A[i].SubMatrix(i, i);

                H[i] = Matrix.Identity(rowLength);
                A[i + 1] = new Matrix(A[i].rows);

                for (int j = i; j < rowLength; j++)
                {
                    for (int k = i; k < colLength; k++)
                    {
                        H[i][j, k] = HMatrix[j - i, k - i];
                        A[i + 1][j, k] = subMatrixA[j - i, k - i];
                    }
                }
                Q = Q * H[i];
                A[i + 1] = H[i] * A[i];
                R = A[i];
            }
            //Console.WriteLine(this);

            return (Q, (this.Transposed * Q).Transposed);
        }

        private Complex.Complex[] Calc2x2EigenValues()
        {
            Complex.Complex[] complex = new Complex.Complex[2];

            var det = this[0, 0] * this[1, 1] - this[0, 1] * this[1, 0];

            float a = 1;
            float b = (-this[0, 0] - this[1, 1]).Real;
            float c = (this[0, 0] * this[1, 1] - this[1, 0] * this[0, 1]).Real;

            float dis = b * b - 4 * a * c;

            if(dis < 0)
            {
                complex[0] = new Complex.Complex(-b / (2 * a), MathF.Sqrt(-dis) / (2 * a));
                complex[1] = new Complex.Complex(-b / (2 * a),-MathF.Sqrt(-dis) / (2 * a));
            }
            else
            {
                complex[0] = new Complex.Complex((-b + MathF.Sqrt(dis)) / (2 * a),0);
                complex[1] = new Complex.Complex((-b - MathF.Sqrt(dis)) / (2 * a),0);
            }

            return complex;
        }

        private Vector power_iteration()
        {
            //float value = 0;
            int iterations = 10;
            Random rand = new Random();
            Vector bk = new float[rowLength];
            for (int i = 0; i < rowLength; i++)
            { // Give random numbers to vector
                bk[i] = new Complex.Complex(rand.Next(1000) / 100f,0);
            }
            for (int i = 0; i < iterations; i++)
            {
                Vector bk_1 = this * bk;

                var bk_1_norm = MathF.Sqrt((bk_1 * bk_1).Real);

                bk = bk_1 * (1f / bk_1_norm);
            }
            return bk;
        }
        
        public Complex.Complex this[int index1, int index2]
        {
            get => rows[index1].nums[index2];
            set => rows[index1].nums[index2] = value;
        }
        
        /// <summary>
        /// The eigenvalues of the matrix
        /// </summary>
        public Complex.Complex[] eigenvalues
        {
            get
            {
                //Complex.Complex[] values = new Complex.Complex[rowLength];
                List<Complex.Complex> values = new();
                Matrix A = new Matrix(this.rows);
                //(Matrix, Matrix) QR = new (Matrix.Identity(rowLength), Matrix.Identity(rowLength));
                for (int i = 0; i < 20 * A.colLength * A.colLength * A.colLength; i++)
                {
                    var QR = A.QRDecomposition();

                    A = QR.Item2 * QR.Item1;
                }

                for (int i = 0; i < rowLength; i++)
                {
                    if (MathF.Abs((this - Matrix.Identity(this.rowLength) * A[i, i]).Determinant) > 0.3)
                        continue;

                    values.Add(new Complex.Complex(A[i, i].Real, 0));
                    //if (i > 0 && MathF.Abs(A[i,i-1].Real) > 1E-20)
                    //{

                    //    if(i == rowLength - 1)
                    //    {
                    //        values.RemoveAt(values.Count - 1);
                    //    }
                    //    continue;
                    //}
                    //else
                    //{
                    //    //values[i] = new Complex.Complex(A[i, i].Real,0);
                    //}
                }

                return values.ToArray();
            }
        }

        /// <summary>
        /// The eigenvectors of the matrix
        /// </summary>
        public Vector[] eigenvectors
        {
            get
            {
                var values = this.eigenvalues;
                Vector[] vecotrs = new Vector[values.Length];
                for (int i = 0; i < vecotrs.Length; i++)
                {
                    var m = this - (Matrix.Identity(rowLength) * values[i]);

                    Vector vec = new Vector(rowLength);
                    // We define the last element to be 1 to avoid trivial solutions
                    vec[0] = new Complex.Complex(1,0);

                    // solve equation
                    var v2 = m.Inverse * vec;

                    //Vector v2 = m ^ vec;

                    // normalize vector by making last value 1
                    vecotrs[i] = 1f / v2[rowLength - 1].Real * v2;

                    if(v2.arr.Select(r => r.Real).Contains(float.NaN))
                    {
                        // idk what to do at this point, this is probably the result of a floating point error
                        throw new Exception();
                    }    
                }

                return vecotrs;
            }
        }

        /// <summary>
        /// The number of rows in the matrix
        /// </summary>
        public int rowLength => rows.Length;

        /// <summary>
        /// The number of columns in the matrix
        /// </summary>
        public int colLength => rows[0].nums.Length;


        /// <summary>
        /// The inverse of the matrix
        /// </summary>
        public Matrix Inverse => (this ^ Identity(rows.Length)).GetMatricies().Item2;

        /// <summary>
        /// The determinant of the matrix
        /// </summary>
        public float Determinant => determinant();

        /// <summary>
        /// The matrix transposed
        /// </summary>
        public Matrix Transposed => Transpose();

        protected Row[] rows;

        #region Constructors

        /// <summary>
        /// Creates a matrix
        /// </summary>
        /// <param name="arr"></param>
        /// <exception cref="InvalidSizeException"></exception>
        public Matrix(float[,] arr)
        {
            int row_count = arr.GetLength(0);
            int col_count = arr.GetLength(1);
            if(row_count == 0 || col_count == 0)
            {
                throw new InvalidSizeException("The size of the arrays rows and columns must be greater than zero");
            }

            this.rows = new Row[row_count];

            for (int i = 0; i < row_count; i++)
            {
                float[] temp = new float[col_count];

                for (int j = 0; j < col_count; j++)
                    temp[j] = arr[i, j];

                this.rows[i] = new Row(temp);
            }
        }

        /// <summary>
        /// Creates a matrix
        /// </summary>
        /// <param name="arr"></param>
        /// <exception cref="InvalidSizeException"></exception>
        public Matrix(float[][] arr)
        {
            //if(arr.Length == 0 || arr[0].Length == 0)
            //{
            //    throw new InvalidSizeException("The size of the input must be greater than zero");
            //}

            int row_count = arr.Length;
            int col_count = arr[0].Length;

            // Its faster to iterate through a span than a array
            var temp_arr = arr.AsSpan();
            this.rows = new Row[row_count];

            for (int i = 0; i < row_count; i++)
            {
                //if (temp_arr[i].Length != col_count)
                //    throw new InvalidSizeException("The size of the columns must all be equal");

                this.rows[i] = new Row(temp_arr[i].AsSpan());
            }
        }

        /// <summary>
        /// Creates a 1xn matrix
        /// </summary>
        /// <param name="arr"></param>
        /// <exception cref="InvalidSizeException"></exception>
        public Matrix(IEnumerable<float> arr)
        {
            if(arr == null || !arr.Any())
                throw new InvalidSizeException("The array must be non null and of size 0 or greater");

            rows = new Row[] { new Row(arr) };
        }

        /// <summary>
        /// Creates a matrix
        /// </summary>
        /// <param name="rows"></param>
        /// <exception cref="InvalidSizeException"></exception>
        public Matrix(IEnumerable<Row> rows)
        {
            if(rows == null || !rows.Any())
                throw new InvalidSizeException("The array must be non null and of size 0 or greater");

            this.rows = (rows.Select(x => new Row(x.nums))).ToArray();
        }

        /// <summary>
        /// Initializes the matrix and creates a zero matrix
        /// </summary>
        /// <param name="rows"></param>
        /// <param name="cols"></param>
        /// <exception cref="InvalidSizeException"></exception>
        public Matrix(int rows, int cols)
        {
            if(rows <= 0 || cols <= 0)
                throw new InvalidSizeException("The rows and column size must 0 or greater");

            this.rows = new Row[rows];

            for (int i = 0; i < rows; i++)
                this.rows[i] = new Row(cols);
        }

        #endregion

        /// <summary>
        /// Checks if the lengths of the matricies rows and columns are equal
        /// </summary>
        /// <param name="mat2"></param>
        /// <returns>True or False</returns>
        public bool EqualSize(Matrix mat2) =>
            this.rows.Length == mat2.rows.Length &&
            this.rows[0].nums.Length == mat2.rows[0].nums.Length;

        /// <summary>
        /// Adds 2 matricies of equal size
        /// </summary>
        /// <param name="mat1"></param>
        /// <param name="mat2"></param>
        /// <returns>The resulting matrix</returns>
        /// <exception cref="InvalidDimensionsException"></exception>
        public static Matrix operator +(Matrix mat1, Matrix mat2)
        {
            if (!mat1.EqualSize(mat2))
                throw new InvalidDimensionsException(
                    "The dimensions of the 2 matricies must be equal, expected " +
                    $"{mat1.rowLength}x{mat1.colLength} + {mat1.rowLength}x{mat1.colLength} but got " +
                    $"{mat1.rowLength}x{mat1.colLength} + {mat1.rowLength}x{mat1.colLength}");

            return new Matrix(mat1.rows.Zip(mat2.rows, (x, y) => x + y));
        }

        /// <summary>
        /// Subtracts 2 matricies of equal size
        /// </summary>
        /// <param name="mat1"></param>
        /// <param name="mat2"></param>
        /// <returns>The resulting matrix</returns>
        /// <exception cref="InvalidDimensionsException"></exception>
        public static Matrix operator -(Matrix mat1, Matrix mat2)
        {
            if (!mat1.EqualSize(mat2))
                throw new InvalidDimensionsException(
                        "The dimensions of the 2 matricies must be equal, " +
                        $"expected {mat1.rowLength}x{mat1.colLength} - {mat1.rowLength}x{mat1.colLength} " +
                        $"but got {mat1.rowLength}x{mat1.colLength} - {mat2.rowLength}x{mat2.colLength}");

            return new Matrix(mat1.rows.Zip(mat2.rows, (x, y) => x - y));
        }

        /// <summary>
        /// Multiplies a matrix by a real scalar
        /// </summary>
        /// <param name="mat"></param>
        /// <param name="multiplier"></param>
        /// <returns>The resulting matrix</returns>
        /// <exception cref="ArgumentException"></exception>
        public static Matrix operator *(Matrix mat, float multiplier)
        {
            if (multiplier == float.NaN)
                throw new ArgumentException("The multiplier must not be NaN");

            return new Matrix(mat.rows.Select(x => x * multiplier));
        }

        /// <summary>
        /// Multiplies a matrix by a real scalar
        /// </summary>
        /// <param name="mat"></param>
        /// <param name="multiplier"></param>
        /// <returns>The resulting matrix</returns>
        /// <exception cref="ArgumentException"></exception>
        public static Matrix operator *(float multiplier, Matrix mat) =>
            mat * multiplier;

        /// <summary>
        /// Multiplies a matrix by a complex scalar
        /// </summary>
        /// <param name="mat"></param>
        /// <param name="multiplier"></param>
        /// <returns>The resulting matrix</returns>
        /// <exception cref="ArgumentException"></exception>
        public static Matrix operator *(Matrix mat, Complex.Complex multiplier)
        {
            if(multiplier.Real == float.NaN || (float)multiplier.Imaginary == float.NaN)
                throw new ArgumentException("The multiplier must not be NaN");

            return new Matrix(mat.rows.Select(x => x * multiplier));
        }

        /// <summary>
        /// Multiplies a matrix by a complex scalar
        /// </summary>
        /// <param name="mat"></param>
        /// <param name="multiplier"></param>
        /// <returns>The resulting matrix</returns>
        /// <exception cref="ArgumentException"></exception>
        public static Matrix operator *(Complex.Complex multiplier, Matrix mat) =>
            mat * multiplier;

        /// <summary>
        /// Multiplies 2 matricies
        /// </summary>
        /// <param name="mat1"></param>
        /// <param name="mat2"></param>
        /// <returns>The resulting matrix</returns>
        /// <exception cref="InvalidDimensionsException"></exception>
        public static Matrix operator *(Matrix mat1, Matrix mat2)
        {
            if(mat1.colLength != mat2.rowLength)
            {
                throw new InvalidDimensionsException(
                    "The number of columns in the first matrix must be equal to the number of rows in the second matrix, " +
                    $"expected {mat1.rowLength}x{mat1.colLength} * {mat1.colLength}x{mat2.colLength} " +
                    $"but got {mat1.rowLength}x{mat1.colLength} * {mat2.rowLength}x{mat2.colLength}");
            }

            var rowCount = mat1.rows.Length;
            var colCount = mat2.rows[0].nums.Length;

            //var rect = new Matrix(rowCount, colCount);

            float[,] ret = new float[rowCount, colCount];

            for (int i = 0; i < mat1.rowLength; i++)
            {
                for (int j = 0; j < colCount; j++)
                {
                    for (int k = 0; k < mat2.rows.Length; k++)
                    {
                        ret[i, j] += (mat1[i, k] * mat2[k, j]).Real;
                    }
                }
            }
            return new Matrix(ret);
        }

        public static Matrix ParallelMultiply(Matrix m1, Matrix m2)
        {
            if (m1.colLength != m2.rowLength)
            {
                throw new InvalidDimensionsException(
                    "The number of columns in the first matrix must be equal to the number of rows in the second matrix, " +
                    $"expected {m1.rowLength}x{m1.colLength} * {m1.colLength}x{m2.colLength} " +
                    $"but got {m1.rowLength}x{m1.colLength} * {m2.rowLength}x{m2.colLength}");
            }

            int row_count = m1.rowLength;
            int col_count = m2.colLength;

            float[,] ret = new float[row_count, col_count];

            Parallel.For(0, m1.rowLength, i =>
            {
                for (int j = 0; j < col_count; j++)
                {
                    for (int k = 0; k < m2.rows.Length; k++)
                    {
                        ret[i,j] += (m1[i, k] * m2[k, j]).Real;
                    }
                }
            });

            //float[][] ret = new float[row_count][];

            //Parallel.For(0, m1.rowLength, i =>
            //{
            //    ret[i] = new float[col_count];
            //    for (int j = 0; j < col_count; j++)
            //    {
            //        for (int k = 0; k < m2.rows.Length; k++)
            //        {
            //            ret[i][j] += (m1[i, k] * m2[k, j]).Real;
            //        }
            //    }
            //});

            return new Matrix(ret);
        }

        /// <summary>
        /// Multiplies a matrix with a vector
        /// </summary>
        /// <param name="mat"></param>
        /// <param name="vec"></param>
        /// <returns>The resulting vector</returns>
        /// <exception cref="InvalidDimensionsException"></exception>
        public static Vector operator *(Matrix mat, Vector vec)
        {
            if(mat.colLength != vec.size)
            {
                throw new InvalidDimensionsException(
                    "The number of columns in the first matrix must be equal to the number of rows in the second matrix, " +
                    $"expected {mat.rowLength}x{mat.colLength} * {mat.colLength}x1 " +
                    $"but got {mat.rowLength}x{mat.colLength} * {vec.size}x1");
            }

            Vector vecOut = new Vector(vec.size);

            int rowCount = mat.rows.Length;

            for (int i = 0; i < rowCount; i++)
            {
                for (int k = 0; k < vecOut.size; k++)
                {
                    vecOut[i] += mat[i, k] * vec[k];
                }
            }
            return vecOut;
        }

        // Implicitly converts arrays to 1xn matricies
        public static implicit operator Matrix(float[] array) =>
            new Matrix(array);


        // Implicitly converts 2d arrays to matricies
        public static implicit operator Matrix(float[,] array) =>
            new Matrix(array);

        /// <summary>
        /// Creates the identity matrix
        /// </summary>
        /// <param name="dimensions"></param>
        /// <returns>A matrix with 1's on the diagonal</returns>
        /// <exception cref="ArgumentException"></exception>
        public static Matrix Identity(int dimensions)
        {
            if(dimensions <= 0)
            {
                throw new ArgumentException("The number of dimensions cant be less or equal to zero");
            }

            Matrix ret = new Matrix(dimensions, dimensions);
            for (int i = 0; i < dimensions; i++)
                ret.rows[i].nums[i] = new Complex.Complex(1,0);
            return ret;
        }

        /// <summary>
        /// Creates an augmented matrix of the 2 matricies
        /// </summary>
        /// <param name="mat1"></param>
        /// <param name="mat2"></param>
        /// <returns>The resulting augmented matrix</returns>
        /// <exception cref="InvalidDimensionsException"></exception>
        public static AugmentedMatrix operator |(Matrix mat1, Matrix mat2) // Agumented matrix
        {
            if(mat1.rowLength != mat2.rowLength)
                throw new InvalidDimensionsException("The row size of the 2 matricies must be equal");

            int rowLen = mat1.rowLength;
            int colLen = mat1.colLength + mat2.colLength;
            Matrix ret = new Matrix(rowLen, colLen);

            for (int i = 0; i < rowLen; i++)
            {
                Row row = new Row(mat1.rows[i].nums.Concat(mat2.rows[i].nums));
                ret.rows[i] = row;
            }

            return new AugmentedMatrix(ret.rows, mat1.rows[0].nums.Length);
        }

        public static bool operator ==(Matrix m1, Matrix m2) =>
            m1.rows.SequenceEqual(m2.rows);

        public static bool operator !=(Matrix m1, Matrix m2) =>
            !m1.rows.SequenceEqual(m2.rows);



        // This swaps the rows where the number on the diagonal is zero
        // The method is only used for gaussian elimination and shouldn't be used anywhere else
        private static (int[] pos, Matrix matrix) RealignMatrix(Matrix mat1)
        {
            Matrix ret = new Matrix(mat1.rows);

            int[] pos = new int[mat1.rowLength];
            var zero = new Complex.Complex(0, 0);


            // ensures that the numbers on diagonal isnt 0
            for (int row = 0; row < ret.rowLength; row++)
            {
                if (ret[row, row] != zero)
                {
                    pos[row] = row;
                    continue;
                }

                Row? NonZeroRow = ret.rows.Skip(row).FirstOrDefault(x => x.nums[row] != zero);
                int colIndex = Array.FindIndex(ret.rows, row, ret.rowLength, x => x.Equals(NonZeroRow));


                if (NonZeroRow is null)
                    throw new Exception("Couldn't realign the matrix");


                ret.rows[colIndex] = ret.rows[row];
                ret.rows[row] = NonZeroRow;

                pos[row] = colIndex;
                pos[colIndex] = row;
              
            }

            return (pos, ret);
        }

        /// <summary>
        /// This does gaussian elimination on the matrix 
        /// </summary>
        /// <param name="mat1"></param>
        /// <param name="mat2"></param>
        /// <returns>The resulting augmented matrix</returns>
        /// <exception cref="InvalidDimensionsException"></exception>
        /// <exception cref="Exception"></exception>
        public static AugmentedMatrix operator ^(Matrix mat1, Matrix mat2) // Gaussian elimination
        {

            (int[] pos, mat1) = RealignMatrix(mat1);
            AugmentedMatrix matrix = mat1 | mat2;


            for (int col = 0; col < matrix.linePos + 1; col++)
            {
                var pivot = matrix[col, col];

                for (int row = col + 1; row < matrix.rows.Length; row++)
                {
                    var multiplier = matrix[row, col] / pivot;
                    matrix.rows[row] -= matrix.rows[col] * multiplier;
                }
                matrix.rows[col] *= 1f / (matrix[col,col]);
            }

            for (int col = matrix.linePos; col > 0; col--)
            {
                var pivot = matrix.rows[col].nums[col];

                for (int row = col - 1; row >= 0; row--)
                {
                    var multiplier = matrix.rows[row].nums[col] / pivot;
                    matrix.rows[row] = matrix.rows[row] - matrix.rows[col] * multiplier;
                }
            }

            for(int i = 0; i < pos.Length; i++)
            {
                matrix.SwapCols(i, pos[i]);
            }

            return matrix;
        }

        // Swaps the columns of a matrix
        private void SwapCols(int c1, int c2)
        {
            for(int i = 0; i < rowLength; i++)
            {
                var temp = this[i, c1];

                this[i, c1] = this[i, c2];
                this[i, c2] = temp;
            }
        }

        // Calculates the determinant
        private float determinant()
        {
            var mat = LUDecomposition();

            float det = 1;
            for (int i = 0; i < rows.Length; i++)
                det *= mat.Item2.rows[i].nums[i].Real;
            return det;
        }

        // Calculates the resulting L,U matrices from LU decomposition
        private (Matrix, Matrix) LUDecomposition()
        {
            Matrix L = Matrix.Identity(rows.Length);
            Matrix U = new Matrix(rows);

            for (int col = 0; col < rows.Length; col++)
            {
                float pivot = U.rows[col].nums[col].Real;
                for (int row = col + 1; row < L.rows.Length; row++)
                {
                    float multiplier = U.rows[row].nums[col].Real / pivot;
                    U.rows[row] = U.rows[row] - U.rows[col] * multiplier;
                    L.rows[row].nums[col] = new Complex.Complex(multiplier,0);
                }
                //U.rows[col] *= 1 / (U.rows[col].nums[col]);
            }
            return (L, U);
        }

        // Transposes the matrix
        private Matrix Transpose()
        {
            Matrix mT = new Matrix(rows[0].nums.Length, rows.Length);

            for (int row = 0; row < rows.Length; row++)
            {
                for (int col = 0; col < rows[row].nums.Length; col++)
                {
                    mT.rows[col].nums[row] = rows[row].nums[col];
                }
            }

            return mT;
        }

        /// <summary>
        /// Calculates the Submatrix of a given matrix
        /// </summary>
        /// <param name="row"></param>
        /// <param name="col"></param>
        /// <returns>The resulting submatrix</returns>
        /// <exception cref="ArgumentException"></exception>
        public Matrix SubMatrix(int row, int col)
        {
            if (row < 0 || col < 0 || row >= rowLength || col >= colLength)
                throw new ArgumentException("The row and column index cant be less than 0 or greater than the size");

            Matrix sub = new Matrix(rowLength - row, colLength - col);

            for (int r = row; r < rowLength; r++)
            {
                for (int c = col; c < colLength; c++)
                {
                    sub[r - row, c - col] = this[r, c];
                }
            }
            return sub;
        }

        /// <summary>
        /// Calculates the Submatrix of a given matrix
        /// </summary>
        /// <param name="m1"></param>
        /// <param name="m2"></param>
        /// <returns>The resulting submatrix</returns>
        /// <exception cref="ArgumentException"></exception>
        public Matrix SubMatrix((int row, int col) m1, (int col, int row) m2)
        {
            if(m2.row < m1.row || m2.col < m1.col)
                throw new ArgumentException("the indexes in the first parameter must be less or equal to the indexes in the second parameter");
            else if (m1.row < 0 || m1.col < 0 || m2.row >= rowLength || m2.col >= colLength)
                throw new ArgumentException("The row and column index cant be less than 0 or greater than the size");

            Matrix sub = new Matrix(m2.row - m1.row+1, m2.col - m1.col+1);

            for(int r = m1.row; r < m2.row+1; r++)
            {
                for(int c = m1.col; c < m2.col+1; c++)
                {
                    sub[r - m1.row, c - m1.col] = this[r, c];
                }
            }

            return sub;
        }
        public override string ToString()
        {
            string output = "";
            for (int i = 0; i < rows.Length; i++)
                output += $"({rows[i]}), \n";
            return output;
        }
    }
}