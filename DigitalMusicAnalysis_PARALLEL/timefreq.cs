using System;
using System.Numerics;
using System.Threading.Tasks;
using System.Diagnostics;

namespace DigitalMusicAnalysis
{
    public class timefreq
    {
        public float[][] timeFreqData;
        public int wSamp;
        public Complex[] twiddles;
        private Stopwatch stopwatch3 = new Stopwatch();
        public static double timeinstft;



        public timefreq(float[] x, int windowSamp)
        {
            double pi = 3.14159265;
            Complex i = Complex.ImaginaryOne;
            this.wSamp = windowSamp;
            twiddles = new Complex[wSamp];

            Parallel.For(0, (wSamp), (int ii, ParallelLoopState state) =>
            {
                double a = 2 * pi * ii / (double)wSamp;
                twiddles[ii] = Complex.Pow(Complex.Exp(-i), (float)a);
            });

            timeFreqData = new float[wSamp / 2][];

            int nearest = (int)Math.Ceiling((double)x.Length / (double)wSamp);
            nearest = nearest * wSamp;

            Complex[] compX = new Complex[nearest];

            for (int kk = 0; kk < nearest; kk++)
            {
                if (kk < x.Length)
                {
                    compX[kk] = x[kk];
                }
                else
                {
                    compX[kk] = Complex.Zero;
                }
            }

            int cols = 2 * nearest / wSamp;

            for (int jj = 0; jj < wSamp / 2; jj++)
            {
                timeFreqData[jj] = new float[cols];
            }
            timeFreqData = stft(compX, wSamp);
        }

        float[][] stft(Complex[] x, int wSamp)
        {
            stopwatch3.Start();
            int ii = 0;
            int jj = 0;
            int kk = 0;
            // int ll = 0;
            int N = x.Length;
            float fftMax = 0;

            float[][] Y = new float[wSamp / 2][];
            /*
            for (ll = 0; ll < wSamp / 2; ll++)
            {
                Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
            }
            */
            Parallel.For(0, wSamp / 2, new ParallelOptions { MaxDegreeOfParallelism = MainWindow.threadPool }, ll =>
            {
                Y[ll] = new float[2 * (int)Math.Floor((double)N / (double)wSamp)];
            });

            Complex[] temp = new Complex[wSamp];
            Complex[] tempFFT = new Complex[wSamp];

            for (ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            {

                for (jj = 0; jj < wSamp; jj++)
                {
                    temp[jj] = x[ii * (wSamp / 2) + jj];
                }

                tempFFT = Iterative(temp, wSamp, twiddles);

                for (kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] = (float)Complex.Abs(tempFFT[kk]);

                    if (Y[kk][ii] > fftMax)
                    {
                        fftMax = Y[kk][ii];
                    }
                }
            }

            for (ii = 0; ii < 2 * Math.Floor((double)N / (double)wSamp) - 1; ii++)
            {
                for (kk = 0; kk < wSamp / 2; kk++)
                {
                    Y[kk][ii] /= fftMax;
                }
            }
            stopwatch3.Stop();
            timeinstft = stopwatch3.Elapsed.TotalMilliseconds;
            return Y;
        }

        static int ReverseBits(int bit, int n)
        {
            int reverse = n;
            int counter = bit - 1;
            n >>= 1;
            while (n > 0)
            {
                reverse = (reverse << 1) | (n & 1);
                counter--;
                n >>= 1;
            }
            return ((reverse << counter) & ((1 << bit) - 1));
        }


        public static  Complex[] Iterative(Complex[] x, int L, Complex[] twiddles) //public static Complex[]
        {
        int N = x.Length;
            Complex[] Y = new Complex[N];
            int bits = (int)Math.Log(N, 2);

            for (int i = 0; i < N; i++)
            {
                int pos = ReverseBits(bits, i);
                Y[i] = x[pos];
            }

            for (int ii = 2; ii <= N; ii <<= 1)
            {
                for (int jj = 0; jj < N; jj += ii)
                {
                    for (int kk = 0; kk < ii / 2; kk++)
                    {
                        int e = jj + kk;
                        int o = jj + kk + (ii / 2);
                        Complex even = Y[e];
                        Complex odd = Y[o];

                        Y[e] = even + odd * twiddles[kk * (L / ii)];
                        Y[o] = even + odd * twiddles[(kk + (ii / 2)) * (L / ii)];
                    }
                }
            }
            return Y;
        }
    }
}


