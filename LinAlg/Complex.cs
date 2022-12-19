using System;
using System.Collections.Generic;
using System.Data;
using System.Data.Common;
using System.Linq;
using System.Numerics;
using System.Reflection.Metadata;
using System.Runtime.CompilerServices;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;

namespace LinAlg.Complex
{
    public struct Imaginary
    {
        float Val;

        private Imaginary(float val)
        {
            Val = val;
        }

        public static implicit operator Imaginary(float val) =>
            new Imaginary(val);

        public static explicit operator float(Imaginary imag) =>
            imag.Val;


        public static Imaginary operator *(Imaginary imag, float real) =>
            new Imaginary(imag.Val * real);

        public static Imaginary operator *(float real, Imaginary imag) =>
            imag * real;

        public static Imaginary operator +(Imaginary imag1, Imaginary imag2) =>
            new Imaginary(imag1.Val + imag2.Val);

        public static Imaginary operator -(Imaginary imag1, Imaginary imag2) =>
            imag1 + (-1 * imag2);

        public static float operator *(Imaginary imag1, Imaginary imag2) =>
            -(imag1.Val * imag2.Val);


        public static bool operator ==(Imaginary c1, Imaginary c2) =>
            c1.Val == c2.Val;

        public static bool operator !=(Imaginary c1, Imaginary c2) =>
            c1.Val != c2.Val;

        public static Complex operator ^(Imaginary imag, int pow)
        {
            if (pow == 0)
                return (1, 0);

            Complex temp = (0, imag.Val);
            Complex ret = (0, imag.Val);

            for (int i = 1; i < pow; i++)
                ret = ret * temp;

            return ret;
        }

        // Empty string if Value is 0
        public override string ToString() =>
            (Val >= 0 ? "":"-") + (Val != 1 ? MathF.Abs(Val)+"i" : Val != 0 ? "i":"");

        public override bool Equals(object? obj)
        {
            return obj is Imaginary && ((Imaginary)obj).Val == this.Val;
        }

        public override int GetHashCode()
        {
            return Val.GetHashCode();
        }

    }

    public struct Complex
    {
        public float Real;
        public Imaginary Imaginary;

        public Complex()
        {
            Real = 0;
            Imaginary = 0;
        }

        public Complex(float real, float imaginary)
        {
            Real = real;
            Imaginary = imaginary;
        }

        public static implicit operator Complex((float, float) complex) =>
            new Complex(complex.Item1, complex.Item2);

        #region Complex-Real opeartors
        public static Complex operator +(Complex complex, float real) =>
            new Complex(complex.Real + real, (float)complex.Imaginary);

        public static Complex operator +(float real, Complex complex) =>
            complex + real;

        public static Complex operator -(Complex complex, float real) =>
            complex + (real * -1);

        public static Complex operator -(float real, Complex complex) =>
            real + (complex * -1);

        public static Complex operator *(Complex complex, float real) =>
            new Complex(complex.Real * real, (float)complex.Imaginary * real);

        public static Complex operator *(float real, Complex complex) =>
            complex * real;

        public static Complex operator /(float real, Complex complex) =>
            new Complex(real, 0) / complex;

        public static Complex operator /(Complex complex, float real) =>
            complex * (1/real);
        #endregion

        #region Complex-Imaginary operators
        public static Complex operator +(Complex complex, Imaginary imag) =>
            new Complex(complex.Real, (float)complex.Imaginary + (float)imag);

        public static Complex operator +(Imaginary imag, Complex complex) =>
            complex + imag;

        public static Complex operator -(Complex complex, Imaginary imag) =>
            complex + (-1 * imag);

        public static Complex operator -(Imaginary imag, Complex complex) =>
            (-1*complex) + imag;
        #endregion

        #region Complex-Complex operators
        public static Complex operator +(Complex c1, Complex c2) =>
            new Complex(c1.Real + c2.Real, (float)c1.Imaginary + (float)c2.Imaginary);

        public static Complex operator -(Complex c1, Complex c2) =>
            c1 + c2 * -1;

        public static Complex operator *(Complex c1, Complex c2) =>
            new Complex(c1.Real * c2.Real + c1.Imaginary * c2.Imaginary, (float)(c1.Real * c2.Imaginary + c2.Real * c1.Imaginary));

        public static Complex operator /(Complex c1, Complex c2) =>
            new Complex(
                (c1.Real * c2.Real + (float)c1.Imaginary * (float)c2.Imaginary)/ // Real part
                (c2.Real * c2.Real + (float)c2.Imaginary * (float)c2.Imaginary),
                (c2.Real * (float)c1.Imaginary - c1.Real * (float)c2.Imaginary)/ // Imaginary part
                (c2.Real * c2.Real + (float)c2.Imaginary * (float)c2.Imaginary));

        public static bool operator ==(Complex c1, Complex c2) =>
            c1.Real == c2.Real && c1.Imaginary == c2.Imaginary;

        public static bool operator !=(Complex c1, Complex c2) =>
            c1.Real != c2.Real || c1.Imaginary != c2.Imaginary;
        #endregion

        #region Singular operators
        public static Complex operator -(Complex c1) =>
            c1 * -1;
        #endregion

        public override string ToString()
        {
            string real_part = "";
            string opeartor = "";
            string imaginary_part = "";

            string total = "";

            if (Real != 0)
            {
                real_part = Real.ToString();
                if ((float)Imaginary != 0)
                    opeartor = (float)Imaginary > 0 ? " + " : " - "; 
            }
            else if (Imaginary != 0)
                opeartor = (float)Imaginary < 0 ? "-" : "";

            if (Imaginary != 0)
            {
                if (MathF.Abs((float)Imaginary) != 1)
                    imaginary_part += MathF.Abs((float)Imaginary).ToString();

                imaginary_part +=  "i";
            }

            total = real_part + opeartor + imaginary_part;

            return total.Length != 0 ? total : "0";
        }

        public override bool Equals(object? obj)
        {
            return (obj is Complex) && (Complex)(obj) == this;
        }

        public override int GetHashCode()
        {
            return Real.GetHashCode() ^ Imaginary.GetHashCode();
        }

        (float r, float theta) ToPolar()
        {
            (float r, float theta) polar;

            polar.r = MathF.Sqrt((Real * Real) + ((float)Imaginary * (float)Imaginary));
            polar.theta = MathF.Acos(this.Real / polar.r);

            return polar;
        }

        static Complex PolarToComplex((float r, float theta) polar) =>
            new Complex(polar.r * MathF.Cos(polar.theta), polar.r * MathF.Sin(polar.theta));

        public Complex Pow(float n)
        {
            var p1 = ToPolar();
            var p2 = (MathF.Pow(p1.r, n), p1.theta * n);
            return PolarToComplex(p2);
        }
    }
}