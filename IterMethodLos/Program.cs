using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace IterMethodLos
{
    internal class Program
    {
        static double NormDown(double[] x)
        {
            int n = x.Length/2;

            double result = 0;

            for (int i = 0; i < n; i++)
            {
                result += x[i] * x[i];
            }

            return Math.Sqrt(result);
        }
        static double NormUp(double[] x, double[] y)
        {
            int n = x.Length/2;

            double result = 0;

            for (int i = 0; i < n; i++)
            {
                double tmp = x[i] - y[i];
                result += tmp * tmp;
            }

            return Math.Sqrt(result);
        }
        static double[] ReadDataFromPardiso(string filename)
        {
            string[] datafromfile = File.ReadAllLines(filename);
            double[] data = new double[datafromfile.Length];

            for (int i =0; i < datafromfile.Length; i++)
            {
                data[i] = double.Parse(datafromfile[i].Replace(",", "."), CultureInfo.InvariantCulture);
            }
            return data;
        }
        static void Main(string[] args)
        {
            SLAEdiag solver = new SLAEdiag();

            Stopwatch stopwatch = new Stopwatch();

            stopwatch.Start();
            solver.LOS_LU();
            stopwatch.Stop();

            Console.WriteLine("time: " + stopwatch.ElapsedMilliseconds);

            double[] x = ReadDataFromPardiso("x4.txt");

            Console.WriteLine("Error: " + (NormUp(x, solver.q) / NormDown(x) * 100).ToString() + " %");

            Console.ReadLine();
        }
    }

}
