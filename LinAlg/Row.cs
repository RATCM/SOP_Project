using LinAlg.Complex;


namespace LinAlg.Matricies
{
    public class Row
    {
        public Complex.ComplexNumber[] nums;

        #region Constructors

        /// <summary>
        /// Creates a row
        /// </summary>
        /// <param name="arr"></param>
        /// <exception cref="ArgumentException"></exception>
        public Row(float[] arr)
        {
            if (arr == null || arr.Length == 0)
                throw new ArgumentException("The collection can't be null or have a size of 0");

            nums = arr.Select(x => new ComplexNumber(x, 0)).ToArray();
        }

        public Row(Span<float> arr)
        {
            nums = new ComplexNumber[arr.Length];
            for (int i = 0; i < arr.Length; i++)
            {
                nums[i] = new ComplexNumber(arr[i], 0);
            }
        }

        /// <summary>
        /// Creates a row
        /// </summary>
        /// <param name="arr"></param>
        /// <exception cref="ArgumentException"></exception>
        public Row(IEnumerable<float> arr)
        {
            if (arr == null || !arr.Any())
                throw new ArgumentException("The collection can't be null or have a size of 0");

            nums = arr.Select(x => new ComplexNumber(x, 0)).ToArray();
        }

        /// <summary>
        /// Creates a row
        /// </summary>
        /// <param name="arr"></param>
        /// <exception cref="ArgumentException"></exception>
        public Row(ComplexNumber[] arr)
        {
            if (arr == null || arr.Length == 0)
                throw new ArgumentException("The collection can't be null or have a size of 0");

            nums = arr.Select(x => x).ToArray();
        }

        /// <summary>
        /// Creates a row
        /// </summary>
        /// <param name="arr"></param>
        /// <exception cref="ArgumentException"></exception>
        public Row(IEnumerable<ComplexNumber> arr)
        {
            if (arr == null || !arr.Any())
                throw new ArgumentException("The collection can't be null or have a size of 0");

            nums = arr.Select(x => x).ToArray();
        }

        /// <summary>
        /// Creates a row with all elements being 0
        /// </summary>
        /// <param name="length"></param>
        /// <exception cref="ArgumentException"></exception>
        public Row(int length)
        {
            if (length <= 0)
                throw new ArgumentException("The length of a row must be greater or equal to 0");

            nums = new Complex.ComplexNumber[length];
        }
        #endregion

        /// <summary>
        /// Adds the 2 rows together
        /// </summary>
        /// <param name="r1"></param>
        /// <param name="r2"></param>
        /// <returns>The resulting row</returns>
        /// <exception cref="ArgumentException"></exception>
        public static Row operator +(Row r1, Row r2)
        {
            if (r1.nums.Length != r2.nums.Length)
                throw new ArgumentException("The length of the rows needs to be equal");

            return new Row(r1.nums.Zip(r2.nums, (x, y) => x + y).ToArray());
        }

        /// <summary>
        /// Subtracts the 2 rows
        /// </summary>
        /// <param name="r1"></param>
        /// <param name="r2"></param>
        /// <returns>The resulting row</returns>
        /// <exception cref="ArgumentException"></exception>
        public static Row operator -(Row r1, Row r2)
        {
            if (r1.nums.Length != r2.nums.Length)
                throw new ArgumentException("The length of the rows needs to be equal");

            return new Row(r1.nums.Zip(r2.nums, (x, y) => x - y).ToArray());
        }

        /// <summary>
        /// Multiplies the row with a number
        /// </summary>
        /// <param name="r1"></param>
        /// <param name="multiplier"></param>
        /// <returns>The resulting row</returns>
        /// <exception cref="ArgumentException"></exception>
        public static Row operator *(Row r1, float multiplier)
        {
            if (multiplier == float.NaN)
                throw new ArgumentException("Multiplier cannot be NaN");

            return new Row(r1.nums.Select(x => x * multiplier).ToArray());
        }

        /// <summary>
        /// Multiplies the row with a number
        /// </summary>
        /// <param name="r1"></param>
        /// <param name="multiplier"></param>
        /// <returns>The resulting row</returns>
        /// <exception cref="ArgumentException"></exception>
        public static Row operator *(Row r1, Complex.ComplexNumber multiplier)
        {
            if (multiplier.Real == float.NaN || (float)multiplier.Imaginary == float.NaN)
                throw new ArgumentException("Multiplier cannot be NaN");

            return new Row(r1.nums.Select(x => x * multiplier).ToArray());
        }

        #region Comparison
        public static bool operator ==(Row r1, Row r2) =>
            r1.nums.SequenceEqual(r2.nums);

        public static bool operator !=(Row r1, Row r2) =>
            !r1.nums.SequenceEqual(r2.nums);

        public override bool Equals(object? obj) =>
            obj is Row && obj != null
            ? this == (Row)obj
            : false;
        #endregion

        public override string ToString()
        {
            string output = "";
            foreach (var num in nums)
                output += $" {(num)}";
            return output;
        }

        public override int GetHashCode()
        {
            return nums.GetHashCode();
        }
    }
}