namespace LinAlg.Matricies
{
    public class AugmentedMatrix : Matrix
    {
        public int linePos;
        public AugmentedMatrix(IEnumerable<Row> rows, int first_matrix_length) : base(rows)
        {
            linePos = first_matrix_length - 1;
        }
        public override string ToString()
        {
            string output = "";

            for (int i = 0; i < rows.Length; i++)
            {
                output += "[";
                for (int j = 0; j < rows[0].nums.Length; j++)
                {
                    output += $" ({rows[i].nums[j]}) ";
                    if (j == linePos)
                    {
                        output += "] | [";
                    }
                }
                output += "]\n";
            }

            return output;
        }

        public (Matrix, Matrix) GetMatricies()
        {

            Matrix mat1 = new(rows.Select(x => new Row(x.nums.Take(linePos + 1))));
            Matrix mat2 = new(rows.Select(x => new Row(x.nums.Skip(linePos + 1))));

            return (mat1, mat2);
        }

    }
}