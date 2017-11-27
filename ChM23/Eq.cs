using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ChM23
{
    public class Eq
    {
        public int n;
        public double norm;
        public double[,] matr;
        public double[] b, x;

        public Eq(double[,] m, double[] ba)
        {
            n = ba.Length;
            x = new double[n];
            b = new double[ba.Length];
            matr = new double[ba.Length, ba.Length];

            norm = 0;
            for (int i = 0; i < ba.Length; i++)
            {
                b[i] = ba[i];
                double sum = 0;
                for (int j = 0; j < ba.Length; j++)
                {
                    matr[i, j] = m[i, j];
                    sum += Math.Abs(m[j, i]);
                }
                if (sum > norm)
                    norm = sum;
            }
        }

        public bool IsDP()
        {
            bool isDP = true;
            for (int i = 0; i < n; i++)
            {
                double sum = -Math.Abs(matr[i, i]);
                for (int j = 0; j < n; j++)
                {
                    sum += Math.Abs(matr[i, j]);
                }
                if (matr[i, i] > sum)
                    continue;
                else
                {
                    isDP = false;
                    break;
                }
            }
            return isDP;
        }

        public bool IsPO()
        {
            bool isPO = true;
            for (int i = 1; i <= n; i++)
            {
                double[,] temp = new double[i, i];
                for (int j = 0; j < i; j++)
                {
                    for (int k = 0; k < i; k++)
                    {
                        temp[j, k] = matr[j, k];
                    }
                }
                double det = Eq.Det(temp);
                if (det > 0)
                    continue;
                else
                {
                    isPO = false;
                    break;
                }
            }
            return isPO;
        }

        public void Diag()
        {
            for (int i = 1; i < this.n; i++)
            {
                this.Sort(i);
                for (int j = i; j < this.n; j++)
                {
                    double a = this.matr[j, i - 1];

                    this.b[j] = this.b[j] - this.b[i - 1] * (a / this.matr[i - 1, i - 1]);

                    for (int m = i - 1; m < this.n; m++)
                    {
                        this.matr[j, m] = this.matr[j, m] - this.matr[i - 1, m] * (a / this.matr[i - 1, i - 1]);
                    }
                }
            }
        }

         public void Sort(int l)
        {
            int k = this.n - 1;
            if (this.matr[l - 1, l - 1] == 0)
            {
                while (this.matr[k, l - 1] == 0)
                {
                    if (k != 0)
                    {
                        k--;
                    }
                    else
                        break;
                }

                if (this.matr[k, l - 1] != 0)
                {
                    double temp = this.b[l - 1];
                    this.b[l - 1] = this.b[k];
                    this.b[k] = temp;

                    for (int j = l - 1; j < this.n; j++)
                    {
                        temp = this.matr[l - 1, j];
                        this.matr[l - 1, j] = this.matr[k, j];
                        this.matr[k, j] = temp;
                    }
                }
            }
        }

        public static double Det(double[,] ma)
        {
            double det = 1;
            double[] b = new double[(int)Math.Sqrt(ma.Length)];

            Eq q = new Eq(ma, b);

            q.Diag();

            for (int i = 0; i < q.n; i++)
            {
                det *= q.matr[i, i];
            }

            return det;
        }        
    }

}