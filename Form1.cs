using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;
using System.Windows.Forms;

namespace MattRaybuckPCA
{
    public partial class Form1 : Form
    {

        Bitmap image;
        private int imageStrideWidth, imageHeight;
        private byte[] greyscalePixelBuffer, CrPixelBuffer, CbPixelBuffer, eigenYPixelBuffer, eigenCrPixelBuffer, eigenCbPixelBuffer;
        private byte[] PixelBuffer;
        private double[,] covariance_matrix;
        private double[] eigenValues, firstEigenvector, secondEigenvector, thirdEigenvector;

        public Form1()
        {
            InitializeComponent();
            this.Text = "Principal Component Analysis";
            this.Name = "Matt Raybuck";

            energyButton.Enabled = false;
            PC_Button.Enabled = false;
            compareY.Enabled = false;
            compareCr.Enabled = false;
            compareCb.Enabled = false;
        }


        #region Form Internal Functions

        private void SeperateImageComponents(Bitmap inputImage)
        {
            Bitmap resultBitmap = new Bitmap(inputImage.Width, inputImage.Height);

            // First create a byte array for the greyscale image. Then I don't have to worry about pointers or mark the function unsafe.
            BitmapData inputImageData = inputImage.LockBits(new Rectangle(0, 0, inputImage.Width, inputImage.Height),
                                                            ImageLockMode.ReadOnly,
                                                            PixelFormat.Format32bppArgb);

            // Set global variables for other function.
            imageStrideWidth = inputImageData.Stride;
            imageHeight = inputImageData.Height;

            PixelBuffer = new byte[imageStrideWidth * imageHeight];
            greyscalePixelBuffer = new byte[imageStrideWidth * imageHeight];
            CrPixelBuffer = new byte[imageStrideWidth * imageHeight];
            CbPixelBuffer = new byte[imageStrideWidth * imageHeight];

            Marshal.Copy(inputImageData.Scan0, PixelBuffer, 0, PixelBuffer.Length);

            inputImage.UnlockBits(inputImageData);

            float YComponent, CrComponent, CbComponent;
            int j = 0;
            double[] Ybuffer = new double[(imageStrideWidth / 4) * imageHeight];
            double[] CrBuffer = new double[(imageStrideWidth / 4) * imageHeight];
            double[] CbBuffer = new double[(imageStrideWidth / 4) * imageHeight];

            // Calculate the intensity for each pixel and replace the original pixel data with the greyscale data
            for (int i = 0; i < PixelBuffer.Length; i += 4) // each pixel is of the form BGRA so iterate i + 4 to skip to the next pixel in buffer
            {
                j = i / 4;
                YComponent = (PixelBuffer[i] * 0.114f)     // Blue 
                             + (PixelBuffer[i + 1] * 0.587f) // Green
                             + (PixelBuffer[i + 2] * 0.299f);// Red

                CrComponent = (PixelBuffer[i] * -0.081f)     // Blue 
                             + (PixelBuffer[i + 1] * 0.419f) // Green
                             + (PixelBuffer[i + 2] * 0.5f);// Red

                CbComponent = (PixelBuffer[i] * 0.5f)     // Blue 
                             + (PixelBuffer[i + 1] * -0.331f) // Green
                             + (PixelBuffer[i + 2] * -0.169f);// Red

                // Set RGB values all equal to intensity
                Ybuffer[j] = YComponent;
                CrBuffer[j] = CrComponent;
                CbBuffer[j] = CbComponent;

            }

            double maxY = 0, minCr = 0, minCb = 0;

            // Find max value in eigenY for normalization, min values in Cr and Cb for Offset.
            for (int i = 0; i < Ybuffer.Length; i++)
            {
                if (maxY < Ybuffer[i])
                {
                    maxY = Ybuffer[i];
                }
                if (minCr > CrBuffer[i])
                {
                    minCr = CrBuffer[i];
                }
                if (minCb > CbBuffer[i])
                {
                    minCb = CbBuffer[i];
                }
            }

            // Normalize eigen Y, Offset eigen Cr, Offset eigen Cb
            byte YByte, CrByte, CbByte;
            j = 0;
            for (int i = 0; i < Ybuffer.Length; i++)
            {
                if (maxY > 255)
                {
                    // Normalize eigen Y value
                    YByte = (byte)((Ybuffer[i] / maxY) * 255);
                }
                else
                {
                    // Cast eigen Y value to byte value
                    YByte = (byte)Ybuffer[i];
                }

                if (minCr < 0)
                {
                    CrByte = (byte)(CrBuffer[i] + Math.Abs(minCr));
                }
                else
                {
                    CrByte = (byte)CrBuffer[i];
                }

                if (minCb < 0)
                {
                    CbByte = (byte)(CbBuffer[i] + Math.Abs(minCb));
                }
                else
                {
                    CbByte = (byte)CbBuffer[i];
                }

                greyscalePixelBuffer[j] = YByte;
                greyscalePixelBuffer[j + 1] = YByte;
                greyscalePixelBuffer[j + 2] = YByte;
                greyscalePixelBuffer[j + 3] = 255;

                CrPixelBuffer[j] = CrByte;
                CrPixelBuffer[j + 1] = CrByte;
                CrPixelBuffer[j + 2] = CrByte;
                CrPixelBuffer[j + 3] = 255;

                CbPixelBuffer[j] = CbByte;
                CbPixelBuffer[j + 1] = CbByte;
                CbPixelBuffer[j + 2] = CbByte;
                CbPixelBuffer[j + 3] = 255;

                j += 4;
            }


        }

        private void SeperateImageComponentsWithEigenvectors()
        {
            double[] eigenYbuffer = new double[(imageStrideWidth / 4) * imageHeight];
            double[] eigenCrBuffer = new double[(imageStrideWidth / 4) * imageHeight];
            double[] eigenCbBuffer = new double[(imageStrideWidth / 4) * imageHeight];

            double YComponent, CrComponent, CbComponent;

            int j = 0;
            for (int i = 0; i < PixelBuffer.Length; i += 4) // each pixel is of the form BGRA so iterate i + 4 to skip to the next pixel in buffer
            {
                j = i / 4;
                YComponent = (PixelBuffer[i] * firstEigenvector[2])     // Blue 
                             + (PixelBuffer[i + 1] * firstEigenvector[1]) // Green
                             + (PixelBuffer[i + 2] * firstEigenvector[0]);// Red

                CrComponent = (PixelBuffer[i] * secondEigenvector[2])     // Blue 
                             + (PixelBuffer[i + 1] * secondEigenvector[1]) // Green
                             + (PixelBuffer[i + 2] * secondEigenvector[0]);// Red

                CbComponent = (PixelBuffer[i] * thirdEigenvector[2])     // Blue 
                             + (PixelBuffer[i + 1] * thirdEigenvector[1]) // Green
                             + (PixelBuffer[i + 2] * thirdEigenvector[0]);// Red

                eigenYbuffer[j] = YComponent;
                eigenCrBuffer[j] = CrComponent;
                eigenCbBuffer[j] = CbComponent;

            }

            double maxY = 0, minCr = 0, minCb = 0;

            // Find max value in eigenY for normalization, min values in Cr and Cb for Offset.
            for (int i = 0; i < eigenYbuffer.Length; i++)
            {
                if (maxY < eigenYbuffer[i])
                {
                    maxY = eigenYbuffer[i];
                }
                if (minCr > eigenCrBuffer[i])
                {
                    minCr = eigenCrBuffer[i];
                }
                if (minCb > eigenCbBuffer[i])
                {
                    minCb = eigenCbBuffer[i];
                }
            }


            eigenYPixelBuffer = new byte[imageStrideWidth * imageHeight];
            eigenCrPixelBuffer = new byte[imageStrideWidth * imageHeight];
            eigenCbPixelBuffer = new byte[imageStrideWidth * imageHeight];

            // Normalize eigen Y, Offset eigen Cr, Offset eigen Cb
            byte eigenYByte, eigenCrByte, eigenCbByte;
            j = 0;
            for (int i = 0; i < eigenYbuffer.Length; i++)
            {
                if (maxY > 255)
                {
                    // Normalize eigen Y value
                    eigenYByte = (byte)((eigenYbuffer[i] / maxY) * 255);
                }
                else
                {
                    // Cast eigen Y value to byte value
                    eigenYByte = (byte)eigenYbuffer[i];
                }

                if (minCr < 0)
                {
                    eigenCrByte = (byte)(eigenCrBuffer[i] + Math.Abs(minCr));
                }
                else
                {
                    eigenCrByte = (byte)eigenCrBuffer[i];
                }

                if (minCb < 0)
                {
                    eigenCbByte = (byte)(eigenCbBuffer[i] + Math.Abs(minCb));
                }
                else
                {
                    eigenCbByte = (byte)eigenCbBuffer[i];
                }

                eigenYPixelBuffer[j] = eigenYByte;
                eigenYPixelBuffer[j + 1] = eigenYByte;
                eigenYPixelBuffer[j + 2] = eigenYByte;
                eigenYPixelBuffer[j + 3] = 255;

                eigenCrPixelBuffer[j] = eigenCrByte;
                eigenCrPixelBuffer[j + 1] = eigenCrByte;
                eigenCrPixelBuffer[j + 2] = eigenCrByte;
                eigenCrPixelBuffer[j + 3] = 255;

                eigenCbPixelBuffer[j] = eigenCbByte;
                eigenCbPixelBuffer[j + 1] = eigenCbByte;
                eigenCbPixelBuffer[j + 2] = eigenCbByte;
                eigenCbPixelBuffer[j + 3] = 255;

                j += 4;
            }

        }

        private Bitmap ByteArrayToBitmap(byte[] inputByteArray)
        {
            Bitmap resultBitmap = new Bitmap(image.Width, image.Height);

            BitmapData resultBitmapData = resultBitmap.LockBits(new Rectangle(0, 0, resultBitmap.Width, resultBitmap.Height), ImageLockMode.WriteOnly, PixelFormat.Format32bppArgb);

            Marshal.Copy(inputByteArray, 0, resultBitmapData.Scan0, inputByteArray.Length);

            resultBitmap.UnlockBits(resultBitmapData);

            return resultBitmap;
        }

        private double[,] CalculateCoVarianceMatrix(byte[] buffer)
        {
            double[,] covarianceMatrix = new double[3, 3];

            double rBar, gBar, bBar;
            int rSum = 0, gSum = 0, bSum = 0;

            for (int i = 0; i < buffer.Length; i+=4)
            {
                rSum += buffer[i + 2]; // Red Pixel Value
                gSum += buffer[i + 1]; // Green Pixel Value
                bSum += buffer[i];     // Blue Pixel Value
            }

            rBar = rSum / (image.Width * image.Height);
            gBar = gSum / (image.Width * image.Height);
            bBar = bSum / (image.Width * image.Height);

            double Red, Green, Blue, res;

            for (int i = 0; i < buffer.Length; i+=4)
            {
                Red = buffer[i + 2] - rBar;
                Green = buffer[i + 1] - gBar;
                Blue = buffer[i] - bBar;

                // Row 1
                covarianceMatrix[0, 0] += (Red) * (Red);
                covarianceMatrix[0, 1] += (Red) * (Green);
                covarianceMatrix[0, 2] += (Red) * (Blue);

                // Row 2
                covarianceMatrix[1, 1] += (Green) * (Green);
                covarianceMatrix[1, 2] += (Green) * (Blue);

                // Row 3
                covarianceMatrix[2, 2] += (Blue) * (Blue);
            }
            res = (1 / ((image.Width * image.Height) - 1.0));

            // Row 1
            covarianceMatrix[0, 0] = covarianceMatrix[0, 0] * res;
            covarianceMatrix[0, 1] = covarianceMatrix[0, 1] * res;
            covarianceMatrix[0, 2] = covarianceMatrix[0, 2] * res;

            // Row 2
            covarianceMatrix[1, 0] = covarianceMatrix[0, 1];
            covarianceMatrix[1, 1] = covarianceMatrix[1, 1] * res;
            covarianceMatrix[1, 2] = covarianceMatrix[1, 2] * res;

            // Row 3
            covarianceMatrix[2, 0] = covarianceMatrix[0, 2];
            covarianceMatrix[2, 2] = covarianceMatrix[2, 2] * res;
            covarianceMatrix[2, 1] = covarianceMatrix[1, 2];

            return covarianceMatrix;
        }

        private double[] CalculateCharacteristicPolynomialCoefs(double[,] covariance_matrix)
        {
            double[] polyCoefs = new double[4];
            double traceOfSquaredMatrix = 0;

            // coef a: set equal to 1
            polyCoefs[0] = 1;

            // coef b: tr(covariance matrix)
            polyCoefs[1] = covariance_matrix[0, 0] + covariance_matrix[1, 1] + covariance_matrix[2, 2];

            // coef c: tr(covariance matrix)^2 - tr(covariance matrix ^2)
            traceOfSquaredMatrix = ((covariance_matrix[0,0] * covariance_matrix[0,0]) + (covariance_matrix[0, 1] * covariance_matrix[1, 0]) + (covariance_matrix[0, 2] * covariance_matrix[2, 0]))      // Row 1 Column 1
                                   + ((covariance_matrix[1, 0] * covariance_matrix[0, 1]) + (covariance_matrix[1, 1] * covariance_matrix[1, 1]) + (covariance_matrix[1, 2] * covariance_matrix[2, 1]))  // Row 2 Column 2
                                   + ((covariance_matrix[2, 0] * covariance_matrix[0, 2]) + (covariance_matrix[2, 1] * covariance_matrix[1, 2]) + (covariance_matrix[2, 2] * covariance_matrix[2, 2])); // Row 3 Column 3

            polyCoefs[2] = (1 / 2.0) * (Math.Pow(covariance_matrix[0, 0] + covariance_matrix[1, 1] + covariance_matrix[2, 2], 2) - traceOfSquaredMatrix);

            // coef d: find determinant of covariance matrix
            polyCoefs[3] = (covariance_matrix[0, 0] * covariance_matrix[1, 1] * covariance_matrix[2, 2])
                           + (covariance_matrix[0, 1] * covariance_matrix[1, 2] * covariance_matrix[2, 0])
                           + (covariance_matrix[0, 2] * covariance_matrix[1, 0] * covariance_matrix[2, 1])
                           - (covariance_matrix[0, 1] * covariance_matrix[1, 0] * covariance_matrix[2, 2])
                           - (covariance_matrix[0, 0] * covariance_matrix[1, 2] * covariance_matrix[2, 1])
                           - (covariance_matrix[0, 2] * covariance_matrix[1, 1] * covariance_matrix[2, 0]);

            return polyCoefs;
        }

        private double[] CalculatePolynomialRoots(double[] polyCoefs)
        {
            double[] x = new double[6], roots = new double[3];
            double deltaT = 0, stepValue = 100;
            double Ft, prevFt = (polyCoefs[0] * Math.Pow(deltaT, 3)) - (polyCoefs[1] * Math.Pow(deltaT, 2)) + (polyCoefs[2] * deltaT) - polyCoefs[3];

            int i = 0, j = 0;

            // First time through, either get lucky and find the roots
            // Or find the bounds which the roots lay between.
            while (i < 3)
            {

                Ft = (polyCoefs[0] * Math.Pow(deltaT, 3)) - (polyCoefs[1] * Math.Pow(deltaT, 2)) + (polyCoefs[2] * deltaT) - polyCoefs[3];

                if ((Math.Sign(Ft) != Math.Sign(prevFt) && Ft != 0))
                {
                    j++;
                    x[j] = deltaT;
                    j++;

                    if (j > 5)
                    {
                        break;
                    }
                }

                if (Ft == 0)
                {
                    roots[i] = deltaT;
                    i++;
                }
                else
                {
                        x[j] = deltaT;
                }

                // Increase deltaT by stepValue
                deltaT = deltaT + stepValue;

                if (Ft != 0)
                {
                    prevFt = Ft;
                }
            }

            switch (i)
            {
                case 0:
                    // No roots have been found
                    roots[2] = FindRootWithinBounds(x[0], x[1], stepValue, polyCoefs);
                    roots[1] = FindRootWithinBounds(x[2], x[3], stepValue, polyCoefs);
                    roots[0] = FindRootWithinBounds(x[4], x[5], stepValue, polyCoefs);
                    break;
                case 1:
                    // One root has been found
                    roots[1] = FindRootWithinBounds(x[2], x[3], stepValue, polyCoefs);
                    roots[2] = FindRootWithinBounds(x[4], x[5], stepValue, polyCoefs);
                    break;
                case 2:
                    // Two roots have been found
                    roots[2] = FindRootWithinBounds(x[4], x[5], stepValue, polyCoefs);
                    break;
                default:
                    break;
            }

            return roots;

        }

        private double FindRootWithinBounds(double lower, double upper, double stepValue, double[] polyCoefs)
        {
            double Ft, prevFt, newLower = lower;

            // previous calculation intialization
            prevFt = (polyCoefs[0] * Math.Pow(lower, 3)) - (polyCoefs[1] * Math.Pow(lower, 2)) + (polyCoefs[2] * lower) - polyCoefs[3];

            while (true)
            {

                // Iterate between the given bounds to try and find the root
                for (double deltaT = lower; deltaT < upper; deltaT = deltaT + stepValue)
                {
                    Ft = ((polyCoefs[0] * Math.Pow(deltaT, 3)) - (polyCoefs[1] * Math.Pow(deltaT, 2)) + (polyCoefs[2] * deltaT) - polyCoefs[3]);

                    // The function of t crossed over the root. So the new upper bound is deltaT
                    if ((Math.Sign(Ft) != Math.Sign(prevFt) && Ft != 0))
                    {
                        upper = deltaT;
                        break;
                    }

                    if ( Math.Floor( Math.Abs(Ft) ) == 0)
                    {
                        // The root was found!
                        return deltaT;
                    }
                    else
                    {
                        // Update the lower bound so the next time this loop runs it will have fewer iterations.
                        lower = deltaT;
                    }

                    if (Ft != 0)
                    {
                        prevFt = Ft;
                    }
                }

                stepValue = stepValue / 10;
            }
        }

        // FindEigenvectors: To do this repeatedly apply Cramers rule to solve the equations produced by multiplying the (varianceCovarianceMatrix - lambda*I) by a 3x1 vector x.
        //                   The third x value in the vector is always taken to be equal to one.
        private void FindEigenvectors(double[] eigenValues, double[,] matrix)
        {
            // Find x11
            double x11 = ApplyCramersRule(eigenValues[0], matrix);

            // Find x12
            double x12 = ApplyCramersRule(eigenValues[0], matrix, false);

            // Normalize to get eigenvector
            firstEigenvector = NormalizeEigenvector(x11, x12);

            // Find x21
            double x21 = ApplyCramersRule(eigenValues[1], matrix);

            // Find x22
            double x22 = ApplyCramersRule(eigenValues[1], matrix, false);

            // Normalize to get eigenvector
            secondEigenvector = NormalizeEigenvector(x21, x22);

            // Find x31
            double x31 = ApplyCramersRule(eigenValues[2], matrix);

            // Find x32
            double x32 = ApplyCramersRule(eigenValues[2], matrix, false);

            // Normalize to get eigenvector
            thirdEigenvector = NormalizeEigenvector(x31, x32);

        }

        private double ApplyCramersRule(double eigenvalue, double[,] matrix, bool firstColumn = true)
        {
            double xValue, denom = (((matrix[0, 0] - eigenvalue) * (matrix[1, 1] - eigenvalue)) - (matrix[1, 0] * matrix[1, 0]));

            if (firstColumn)
            {
                xValue = ( ((-matrix[0, 2]) * (matrix[1, 1] - eigenvalue)) - ((-matrix[1, 2]) * matrix[0, 1]) ) / denom;
            }
            else
            {
                xValue = (((-matrix[1, 2]) * (matrix[0, 0] - eigenvalue)) - ((-matrix[0, 2]) * matrix[0, 1])) / denom;
            }

            return xValue;
        }

        // NormalizeEigenvector: Normalizes the values found by Cramers Rule and returns a unit eigenvector.
        private double[] NormalizeEigenvector(double x11, double x12)
        {
            double[] eigenvector = new double[3];
            double denom;

            denom = Math.Sqrt(Math.Pow(x11, 2) + Math.Pow(x12, 2) + 1);

            eigenvector[0] = x11 / denom;
            eigenvector[1] = x12 / denom;
            eigenvector[2] = 1 / denom;

            return eigenvector;
        }

        // CreateResultWindow: Programmatically creates a new form for a given image and sets the forms text to the given string.
        private void CreateResultWindow(Bitmap resultImage, String windowText)
        {
            ResultForm resultForm = new ResultForm(resultImage);
            resultForm.ResultPictureBox.Image = resultImage;

            // Resize window to fit the image
            resultForm.ResultPictureBox.Size = resultImage.Size;
            resultForm.Width = resultImage.Width + 40;
            resultForm.Height = resultImage.Height + 62;

            resultForm.Text = windowText;
            resultForm.Show();
        }

        #endregion

        #region Form GUI Control Functions

        private void loadPictureButton_Click(object sender, EventArgs e)
        {
            OpenFileDialog LoadedImage = new OpenFileDialog();

            LoadedImage.Filter = "Image Files(*.jpg;*.jpeg;*.gif;*.bmp;*.png)|*.jpg;*.jpeg;*.gif;*.bmp;*.png";

            if (LoadedImage.ShowDialog() == DialogResult.OK)
            {
                image = new Bitmap(LoadedImage.FileName);
                pictureBox.Image = image;
                pictureBox.Size = image.Size;
                this.Width = (image.Width + 100) < 800 ? 800 : (image.Width + 100);
                this.Height = image.Height < 600 ? 600 : image.Height;

                SeperateImageComponents(image);
                energyButton.Enabled = false;
                PC_Button.Enabled = false;
                compareY.Enabled = false;
                compareCr.Enabled = false;
                compareCb.Enabled = false;
            }
        }

        private void exitButton_Click(object sender, EventArgs e)
        {
            Application.Exit();
        }

        private void PCA_Button_Click(object sender, EventArgs e)
        {
            if (image != null)
            {

                covariance_matrix = CalculateCoVarianceMatrix(PixelBuffer);
                double[] polynomialCoefs = CalculateCharacteristicPolynomialCoefs(covariance_matrix);
                eigenValues = CalculatePolynomialRoots(polynomialCoefs);

                EigenValueForm eigenvaluesWindow = new EigenValueForm(eigenValues, covariance_matrix);
                eigenvaluesWindow.Show();

                energyButton.Enabled = true;

            }

        }

        private void energyButton_Click(object sender, EventArgs e)
        {
            // Display the energy of the eigenvalues
            EnergyForm energyWindow = new EnergyForm(eigenValues);
            energyWindow.Show();

            // Find the eigenvectors
            FindEigenvectors(eigenValues, covariance_matrix);

            PC_Button.Enabled = true;
        }


        private void PC_Button_Click(object sender, EventArgs e)
        {
            // Display the eigenvectors
            PrincipalComponentsForm PC_Window = new PrincipalComponentsForm(firstEigenvector, secondEigenvector, thirdEigenvector);
            PC_Window.Show();

            // Calculate the Y Cr Cb with the eigenvectors
            SeperateImageComponentsWithEigenvectors();

            compareY.Enabled = true;
            compareCr.Enabled = true;
            compareCb.Enabled = true;
        }

        private void compareY_Click(object sender, EventArgs e)
        {
            CreateResultWindow(ByteArrayToBitmap(greyscalePixelBuffer), "Original Y");
            CreateResultWindow(ByteArrayToBitmap(eigenYPixelBuffer), "PC Y");
        }

        private void compareCr_Click(object sender, EventArgs e)
        {
            CreateResultWindow(ByteArrayToBitmap(CrPixelBuffer), "Original Cr");
            CreateResultWindow(ByteArrayToBitmap(eigenCrPixelBuffer), "PC Cr");
        }

        private void compareCb_Click(object sender, EventArgs e)
        {
            CreateResultWindow(ByteArrayToBitmap(CbPixelBuffer), "Original Cb");
            CreateResultWindow(ByteArrayToBitmap(eigenCbPixelBuffer), "PC Cb");
        }

        #endregion


    }
}
