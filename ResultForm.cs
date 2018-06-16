using System.Drawing;
using System.Windows.Forms;

namespace MattRaybuckPCA
{
    public partial class ResultForm : Form
    {
        private Bitmap image;

        public ResultForm(Bitmap inputimage)
        {
            image = inputimage;

            InitializeComponent();
        }

        public PictureBox ResultPictureBox
        {
            get
            {
                return this.pictureBox;
            }
        }
    }
}
