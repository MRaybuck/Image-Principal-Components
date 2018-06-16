using System;
using System.Windows.Forms;

namespace MattRaybuckPCA
{
    public partial class PrincipalComponentsForm : Form
    {
        public PrincipalComponentsForm(double[] first, double[] second, double[] third)
        {
            InitializeComponent();

            // Eigenvectors
            eigenvector1_Label.Text = Math.Round(first[0], 5).ToString();
            eigenvector2_Label.Text = Math.Round(first[1], 5).ToString();
            eigenvector3_Label.Text = Math.Round(first[2], 5).ToString();

            eigenvector4_Label.Text = Math.Round(second[0], 5).ToString();
            eigenvector5_Label.Text = Math.Round(second[1], 5).ToString();
            eigenvector6_Label.Text = Math.Round(second[2], 5).ToString();

            eigenvector7_Label.Text = Math.Round(third[0], 5).ToString();
            eigenvector8_Label.Text = Math.Round(third[1], 5).ToString();
            eigenvector9_Label.Text = Math.Round(third[2], 5).ToString();
        }
    }
}
