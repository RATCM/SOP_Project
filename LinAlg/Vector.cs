using LinAlg.Extensions;


namespace LinAlg.Matricies
{
    public class Vector
    {
        public Complex.ComplexNumber[] arr;
        public int size => arr.Length;
        public Vector normalized =>
                this * (1f / magnitude);
        public float magnitude =>
            MathF.Sqrt(arr.Select(x => (x.Real, x.Imaginary))
                .Sum(x => x.Real * x.Real + (float)x.Imaginary * (float)x.Imaginary));
        //MathF.Sqrt(arr.Aggregate((x,y) =>
        //x.Real*x.Real + (float)x.Imaginary + (float)x.Imaginary +
        //y.Real*y.Real + (float)y.Imaginary + (float)y.Imaginary));

        public Complex.ComplexNumber norm =>
            (this * this).Pow(0.5f);

        public static Vector operator +(Vector vec1, Vector vec2)
        {
            return new Vector(vec1.arr.Zip(vec2.arr, (x, y) => x + y).ToArray());
        }

        public static Vector operator -(Vector vec1, Vector vec2) =>
            new Vector(vec1.arr.Zip(vec2.arr, (x, y) => x - y).ToArray());

        public static Vector operator *(Vector vec, float scalar) =>
            new Vector(vec.arr.Select(x => x * scalar));
        public static Vector operator *(float scalar, Vector vec) =>
            vec * scalar;

        public static Vector operator *(Vector vec, Complex.ComplexNumber scalar) =>
            new Vector(vec.arr.Select(x => x * scalar));
        public static Vector operator *(Complex.ComplexNumber scalar, Vector vec) =>
            new Vector(vec.arr.Select(x => x * scalar));

        private Vector(IEnumerable<float> array) =>
            arr = array.Select(x => new Complex.ComplexNumber(x, 0)).ToArray();

        private Vector(IEnumerable<Complex.ComplexNumber> array) =>
            arr = array.Select(x => x).ToArray();

        public Vector(int len) =>
            arr = new Complex.ComplexNumber[len];

        public static implicit operator Vector(float[] array) =>
            new Vector(array);

        public Matrix ToMatrix()
        {
            var mOut = new Matrix(size, 1);
            for (int i = 0; i < size; i++)
            {
                mOut[i, 0] = this[i];
            }
            return mOut;
        }

        public Complex.ComplexNumber this[int index]
        {
            get => arr[index];
            set => arr[index] = value;
        }

        public override string ToString()
        {
            string output = "(";
            for (int i = 0; i < arr.Length; i++)
            {
                output += $"{arr[i]}, ";
            }
            output = output.RemoveLast(2);
            output += ")";
            return output;
        }
        public static Complex.ComplexNumber operator *(Vector v1, Vector v2) =>
            v1.arr.Zip(v2.arr, (x, y) => x * y).Aggregate((x, y) => x + y);
    }
}