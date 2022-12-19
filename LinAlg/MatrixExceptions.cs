using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace LinAlg
{
    public class InvalidDimensionsException : Exception
    {
        public InvalidDimensionsException(string message)
            : base(message) { }
    }

    public class InvalidSizeException : Exception
    {
        public InvalidSizeException(string message)
            : base(message) { }
    }
}
