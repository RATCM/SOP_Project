using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace LinAlg.Matrix
{
    public class Row
    {
        public Complex.Complex[] nums;
        public Row(float[] arr) { nums = arr.Select(x => new Complex.Complex(x,0)).ToArray(); }
        public Row(IEnumerable<float> arr) { nums = arr.Select(x => new Complex.Complex(x, 0)).ToArray(); }
        public Row(Complex.Complex[] arr) { nums = arr.Select(x => x).ToArray(); }
        public Row(IEnumerable<Complex.Complex> arr) { nums = arr.Select(x => x).ToArray(); }

        public Row(int length) { nums = new Complex.Complex[length]; }


        public static Row operator +(Row r1, Row r2) =>
            new Row(r1.nums.Zip(r2.nums, (x, y) => x + y).ToArray());

        public static Row operator -(Row r1, Row r2) =>
            new Row(r1.nums.Zip(r2.nums, (x, y) => x - y).ToArray());
        public static Row operator *(Row r1, float multiplier) =>
            new Row(r1.nums.Select(x => x * multiplier).ToArray());

        public static Row operator *(Row r1, Complex.Complex multiplier) =>
            new Row(r1.nums.Select(x => x * multiplier).ToArray());


        #region Comparison
        public static bool operator ==(Row r1, Row r2) =>
            r1.nums.SequenceEqual(r2.nums);
        public static bool operator !=(Row r1, Row r2) =>
            !r1.nums.SequenceEqual(r2.nums);
        #endregion

        public override string ToString()
        {
            string output = "";
            foreach (var num in nums)
                output += $" {(num)}";
            return output;
        }

        public override bool Equals(object? obj) =>
            obj is Row && obj != null
            ? this == (Row)obj
            : false;

        public override int GetHashCode()
        {
            return nums.GetHashCode();
        }
    }
}