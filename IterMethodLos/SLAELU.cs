using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Runtime.InteropServices.ComTypes;
using System.Text;
using System.Threading.Tasks;


namespace IterMethodLos
{
    public class SLAELU //класс для решения СЛАУ
    {
        public static int N;
        public static int max_iter;
        public static int[] ig;
        public static int[] jg;
        public static double eps;
        public static double[] di;
        public static double[] ggl;
        public static double[] ggu;

        // for Solver
        public static double[] l, u, d, d1;
        public static double[] F, temp, temp0;
        public static double[] r, z, p;

        public double[] q { get; set; }
        // LU факторизация
        static void CalcLU()
        {
            //for (int i = 0; i < di.Length; i++)
            //    d[i] = Math.Sqrt(di[i]);
            //for (int i = 0; i < di.Length; i++)
            //    d1[i] = 1 / d[i];
            //for (int i = 0; i < ggl.Length; i++)
            //    l[i] = 0;
            //for (int i = 0; i < ggu.Length; i++)
            //    u[i] = 0;


            for (int i = 0; i < di.Length; i++)
                d[i] = di[i];
            for (int i = 0; i < ggl.Length; i++)
                l[i] = ggl[i];
            for (int i = 0; i < ggu.Length; i++)
                u[i] = ggu[i];

            double sumU, sumL, sumD;
            int n = N;

            for (int i = 0; i < n; i++)
            {
                sumD = 0;

                int begI = ig[i];
                int endI = ig[i + 1];
                for (int igi = begI; igi < endI; igi++)
                {
                    sumU = 0;
                    sumL = 0;

                    int Jindex = jg[igi];

                    for (int igj = begI; igj < igi; igj++)
                    {
                        int begJ = ig[Jindex];
                        int endJ = ig[Jindex + 1];

                        for (int jgi = begJ; jgi < endJ; jgi++)
                        {
                            if (jg[igj] == jg[jgi])
                            {
                                sumL += l[igj] * u[jgi];
                                sumU += l[jgi] * u[igj];
                            }
                        }
                    }
                    l[igi] -= sumL;
                    u[igi] -= sumU;
                    u[igi] /= d[Jindex];
                    sumD += l[igi] * u[igi];
                }

                d[i] -= sumD;
            }
        }
        // Прямой ход Ly = F
        static void CalcDir(double[] y, double[] F)
        {
            double sum, buf;
            int n = N;

            for (int i = 0; i < n; i++)
            {
                y[i] = F[i];
            }

            for (int i = 0; i < n; i++)
            {
                sum = 0;

                int begI = ig[i];
                int endI = ig[i + 1];

                for (int igi = begI; igi < endI; igi++)
                {
                    sum += y[jg[igi]] * l[igi];
                }

                buf = y[i] - sum;

                y[i] = buf / d[i];
            }
            //for (int i = 0; i < n; i++)
            //    y[i] = F[i] * d1[i];
        }
        // Обратный ход Ux = y
        static void CalcRev(double[] x, double[] y)
        {
            int n = N;

            for (int i = 0; i < n; i++)
            {
                x[i] = y[i];
            }

            for (int i = n - 1; i >= 0; i--)
            {
                int begI = ig[i];
                int endI = ig[i + 1];

                for (int igi = begI; igi < endI; igi++)
                {
                    x[jg[igi]] -= x[i] * u[igi];
                }
            }
        }
        // Процедура умножения матрицы на вектор Ax = res
        static void MultMV(int[] ig, int[] jg, double[] x, double[] res)
        {
            int n = x.Length;

            for (int i = 0; i < n; i++)
            {
                res[i] = di[i] * x[i];

                int begI = ig[i];
                int endI = ig[i + 1];

                for (int igi = begI; igi < endI; igi++)
                {
                    int Jindex = jg[igi];

                    res[i] += ggl[igi] * x[Jindex];
                    res[Jindex] += ggu[igi] * x[i];
                }
            }
        }
        static void MultMV1(int[] ig, int[] jg, double[] x, double[] res)
        {
            int n = x.Length;

            for (int i = 0; i < n; i++)
            {
                res[i] = d[i] * x[i];

                int begI = ig[i];
                int endI = ig[i + 1];

                for (int igi = begI; igi < endI; igi++)
                {
                    int Jindex = jg[igi];

                    res[i] += l[igi] * x[Jindex];
                    //res[Jindex] += ggu[igi] * x[i];
                }
            }
        }
        // Функция скалярного произведение двух векторов
        static double ScalarProd(double[] x, double[] y)
        {
            int n = x.Length;

            double result = 0;

            for (int i = 0; i < n; i++)
            {
                result += x[i] * y[i];
            }

            return result;
        }
        // Локально-оптимальная схема c факторизацией LU
        public void LOS_LU()
        {
            double alpha, beta, norm;


            int n = N, maxiter = max_iter;
            double epsilon = eps;

            CalcLU();
            // A * x0
            MultMV(ig, jg, q, temp);

            // f - A * x0
            for (int i = 0; i < n; i++)
            {
                temp[i] = F[i] - temp[i];
            }
            double norm_r0 = Math.Sqrt(ScalarProd(temp, temp));

            // L * r0 = f - A * x0
            CalcDir(r, temp);

            // U * z0 = r0
            
            CalcRev(z, r);

            // A * z0
            MultMV(ig, jg, z, temp);

            // L * p0 = A * z0
            CalcDir(p, temp);

            int k;

            double[] r_k = new double[n];

            for (k = 0; k < maxiter; k++)
            {
                double sp = ScalarProd(p, p);


                alpha = ScalarProd(p, r) / sp;

                for (int i = 0; i < n; i++)
                {
                    q[i] = q[i] + alpha * z[i];
                    r[i] = r[i] - alpha * p[i];
                }

                // U * temp = r
                
                CalcRev(temp, r);

                // A * U-1 * r = temp0
                MultMV(ig, jg, temp, temp0);

                // L * temp = A * U-1 * r 
                CalcDir(temp, temp0);

                beta = -1 * ScalarProd(p, temp) / sp;

                // U * temp0 = r
                CalcRev(temp0, r);

                MultMV1(ig, jg, r, r_k);

                norm = Math.Sqrt(ScalarProd(r_k, r_k))/norm_r0;

                if(k%10 == 0)
                    Console.WriteLine((k+1).ToString() + " " + (norm).ToString());

                if (norm < epsilon)
                {
                    break;
                }

                for (int i = 0; i < n; i++)
                {
                    z[i] = temp0[i] + beta * z[i];
                    p[i] = temp[i] + beta * p[i];
                }

            }
            Console.WriteLine("\niter: {0}: ", k - 1);
        }
        static void ReadKuslau()
        {
            string[] filetext = File.ReadAllLines("kuslau");
            N = int.Parse(filetext[0]);
            eps = double.Parse(filetext[1].Replace(",", "."), CultureInfo.InvariantCulture);
            max_iter = int.Parse(filetext[2]);
        }
        static void ReadBinOfLong(string filename, int[] mas)
        {
            using (FileStream filestream = new FileStream(filename, FileMode.Open, FileAccess.Read))
            {
                using (var reader = new BinaryReader(filestream))
                {
                    for (int i = 0; i < mas.Length; i++)
                    {
                        mas[i] = reader.ReadInt32();

                    }
                }
            }
        }
        static void ReadBinOfDouble(string filename, double[] mas)
        {
            using (FileStream filestream = new FileStream(filename, FileMode.Open, FileAccess.Read))
            {
                using (var reader = new BinaryReader(filestream))
                {
                    for (int i = 0; i < mas.Length; i++)
                    {
                        mas[i] = reader.ReadDouble();
                    }
                }
            }
        }
        public SLAELU()
        {
            ReadKuslau();

            ig = new int[N + 1];
            ReadBinOfLong("ig", ig);
            for (int i = 0; i < N + 1; i++)
            {
                ig[i]--;
            }

            jg = new int[ig[N]];
            ReadBinOfLong("jg",jg);
            for (int i = 0; i < ig[N]; i++)
            {
                jg[i]--;
            }
            
            di = new double[N];
            ReadBinOfDouble("di", di);
            
            ggl = new double[ig[N]];
            ReadBinOfDouble("gg", ggl);

            ggu = ggl;

            F = new double[N];
            ReadBinOfDouble("pr", F);

            int size = ig[N];

            l = new double[size];
            u = new double[size];

            d = new double[N];
            d1 = new double[N];

            q = new double[N];
            temp = new double[N];
            temp0 = new double[N];
            r = new double[N];
            z = new double[N];
            p = new double[N];
        }
    }

}
