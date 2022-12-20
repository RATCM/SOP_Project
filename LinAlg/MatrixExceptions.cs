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
