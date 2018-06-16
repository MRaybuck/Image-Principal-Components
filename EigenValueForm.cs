using System;
using System.Windows.Forms;

namespace MattRaybuckPCA
{
    public partial class EigenValueForm : Form
    {
        public EigenValueForm(double[] eigenvalues, double[,] covarianceMatrix)
        {
            InitializeComponent();

            double eigensum = Math.Floor(eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);
            double varianceSumValue = Math.Floor(covarianceMatrix[0, 0] + covarianceMatrix[1, 1] + covarianceMatrix[2, 2]);

            // Eigenvalues
            eigenvalue1_Label.Text = Math.Round(eigenvalues[0], 5).ToString();
            eigenvalue2_Label.Text = Math.Round(eigenvalues[1], 5).ToString();
            eigenvalue3_Label.Text = Math.Round(eigenvalues[2], 5).ToString();

            // Eigenvalue Sum
            eigenvalueSum_Label.Text = eigensum.ToString();

            // Variance Sum
            varianceSum_Label.Text = varianceSumValue.ToString();

            // If the variance and the sum of the eigenvalues are the same then the correct eigenvalues have been found.
        }
    }
}
