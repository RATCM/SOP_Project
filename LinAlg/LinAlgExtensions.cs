using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace LinAlg.Extensions
{
    public static class Extensions
    {
        public static string RemoveLast(this string str, int count = 1)
        {
            return str.Substring(0, str.Length - count);
        }
    }
}