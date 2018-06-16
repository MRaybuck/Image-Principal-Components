using System;
using System.Windows.Forms;

namespace MattRaybuckPCA
{
    public partial class EnergyForm : Form
    {
        public EnergyForm(double[] eigenvalues)
        {
            InitializeComponent();

            // Sum
            double eigensum = Math.Floor(eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);

            // Eigenvalues
            eigenvalue1_Label.Text = Math.Round(eigenvalues[0], 5).ToString();
            eigenvalue2_Label.Text = Math.Round(eigenvalues[1], 5).ToString();
            eigenvalue3_Label.Text = Math.Round(eigenvalues[2], 5).ToString();

            // Eigenvalue Energies
            eigEnergy1.Text = ( (Math.Round(eigenvalues[0], 5) / eigensum) * 100).ToString();
            eigEnergy2.Text = ( (Math.Round(eigenvalues[1], 5) / eigensum) * 100).ToString();
            eigEnergy3.Text = ( (Math.Round(eigenvalues[2], 5) / eigensum) * 100).ToString();
        }
    }
}
