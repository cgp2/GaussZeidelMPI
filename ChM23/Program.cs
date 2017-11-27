using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ChM23
{
    class Program
    {
        static void Main(string[] args)
        {
            int k = 4;
            Eq equat;
            Console.WriteLine("Введите размерность матрицы");
            int n = int.Parse(Console.ReadLine());
            double[,] matr = new double[n, n];
            double[] b = new double[n];

            Console.WriteLine("Меню:\n 1) Задать произвольную матрицу \n 2) Задать матрицу из примера");
            int s = int.Parse(Console.ReadLine());

            switch (s)
            {
                case 1:
                    for (int i = 0; i < n; i++)
                    {
                        for (int j = 0; j < n; j++)
                        {
                            Console.WriteLine("Введите элемент матрицы {0}-{1} :", i, j);
                            matr[i, j] = double.Parse(Console.ReadLine());
                        }
                    }
                    for (int i = 0; i < n; i++)
                    {
                        Console.WriteLine("Введите элемент столбца {0} :", i);
                        b[i] = double.Parse(Console.ReadLine());
                    }
                    break;
                case 2:
                    k = 4;

                    for (int i = 0; i < n; i++)
                    {
                        b[i] = 15 + k;
                        k += 2;
                    }

                    k = 2;

                    for (int i = 0; i < n; i++)
                    {
                        for (int j = 0; j < n; j++)
                        {
                            if (i == j)
                            {
                                matr[i, j] = 15 + k;
                                k += 2;
                            }
                            else
                                matr[i, j] = 1;
                        }
                    }

                    break;
            }

            equat = new Eq(matr, b);

            Console.WriteLine("Матрица:");
            for (int i = 0; i < equat.n; i++)
            {
                for (int j = 0; j < equat.n; j++)
                {
                    Console.Write("{0}\t", equat.matr[i, j]);
                }
                Console.WriteLine();
            }

            Console.WriteLine("Cтолбец:");
            foreach (double i in equat.b)
                Console.WriteLine(i);

            Console.WriteLine("*************************************************");
            Console.WriteLine("Методы: \n1)Гаусс \n2)Зейдель \n3)МПИ");
            s = int.Parse(Console.ReadLine());  

            switch (s)
            {
                case 1:
                    equat.x = Gauss(equat);
                    break;
                case 2:
                    equat.x = Zeid(equat);
                    break;
                case 3:
                    equat.x =MPI(equat);
                    break;
            }

            Console.WriteLine("Решение:");
            foreach (double i in equat.x)
                Console.WriteLine(i);

            Console.ReadKey();
        }

        static double[] Gauss(Eq equat)
        {
            double[,] m = new double[equat.n, equat.n];
            for (int i = 0; i < equat.n; i++)
            {
                for (int j = 0; j < equat.n; j++)
                {
                    m[i, j] = equat.matr[i, j];
                }
            }
            double d = Eq.Det(m);
            if (d != 0)
            {
                equat.Diag();
                for (int i = equat.n - 1; i >= 0; i--)
                {
                    double sum = 0;
                    for (int j = equat.n - 1; j > i; j--)
                    {
                        sum += equat.x[j] * equat.matr[i, j];
                    }
                    equat.x[i] = (equat.b[i] - sum) / equat.matr[i, i];
                }
            }
            return equat.x;
        }

        static double[] Zeid(Eq equat)
        {
            if (equat.IsDP() == true && equat.IsPO() == true)
            {
                double[] x1 = new double[equat.n];
                double eps = 0, eps0 = Math.Abs((1 - equat.norm) / equat.norm) * 0.0001;

                do
                {
                    for (int i = 0; i < equat.n; i++)
                        x1[i] = equat.x[i];

                    for (int i = 0; i < equat.n; i++)
                    {
                        double sum = 0;
                        for (int j = 0; j < i; j++)
                            sum += equat.matr[i, j] * equat.x[j];
                        for (int j = i + 1; j < equat.n; j++)
                            sum += equat.matr[i, j] * x1[j];
                        equat.x[i] = (equat.b[i] - sum) / equat.matr[i, i];
                    }
                    eps = 0;
                    for (int i = 0; i < equat.n; i++)
                    {
                        eps += Math.Abs(equat.x[i] - x1[i]);
                    }
                } while (eps > eps0);
                
            }
            return equat.x;
            else
            {
                Console.WriteLine("Матрциа не ПО или не ДП");
                Console.ReadKey();
                Environment.Exit(0);
                return null;
            }
        }

        static double[] MPI(Eq equat)
        {
            Eq equat1;

            if (equat.IsDP() == false)
            {
                if (equat.IsPO() == true)
                {
      

                    double[,] matr1 = new double[equat.n, equat.n];
                    double[] b1 = new double[equat.n];

                    double m = 1 / (equat.norm);

                    for (int i = 0; i < equat.n; i++)
                    {
                        b1[i] = m * equat.b[i];
                        for (int j = 0; j < equat.n; j++)
                        {
                            matr1[i, j] = 1 - m * equat.matr[i, j];
                        }
                    }
                    equat1 = new Eq(matr1, b1);
                    equat1.Diag();
                    double[] x1 = new double[equat1.n], x2 = new double[equat1.n];
                    double eps = 0, eps0 = 0;

                    for (int i = 0; i < equat1.n; i++)
                    {
                        equat1.x[i] = equat1.b[i];
                    }

                    do
                    {
                        for (int i = 0; i < equat1.n; i++)
                            x2[i] = x1[i];
                        for (int i = 0; i < equat1.n; i++)
                            x1[i] = equat1.x[i];
                        for (int i = 0; i < equat1.n; i++)
                        {
                            double var = 0;
                            for (int j = 0; j < equat1.n; j++)
                            {
                                if (i != j)
                                {
                                    var += (equat1.matr[i, j] * x1[j]);
                                }
                            }

                            equat1.x[i] = (equat1.b[i] - var) / equat1.matr[i, i];
                        }
                        eps = 0;

                        for (int i = 0; i < equat1.n; i++)
                        {
                            eps += Math.Abs(equat1.x[i] - x1[i]);
                        }

                        eps0 = 0;
                        for (int i = 0; i < equat1.n; i++)
                        {
                            eps0 += Math.Abs(x1[i] - x2[i]);
                        }
                        eps0 *= equat1.norm / (1 - equat1.norm);
                    } while (eps > eps0);
                    return equat1.x;
                }
                else
                {
                    double[,] matr1 = new double[equat.n, equat.n];
                    double[] b1 = new double[equat.n];

                    for (int i = 0; i < equat.n; i++)
                    {
                        for (int j = 0; j < equat.n; j++)
                        {
                            b1[i] += equat.matr[j, i] * equat.b[j];
                            for (int l = 0; l < equat.n; l++)
                            {
                                matr1[i, j] += equat.matr[l, i] * equat.matr[l, j];
                            }
                        }
                    }

                    equat1 = new Eq(matr1, b1);
                    equat1.Diag();

                    double[] x1 = new double[equat1.n], x2 = new double[equat1.n];
                    double eps = 0, eps0 = Math.Abs((1 - equat.norm) / equat.norm) * 0.0001;

                    for (int i = 0; i < equat1.n; i++)
                    {
                        equat1.x[i] = equat1.b[i];
                    }

                    do
                    {
                        for (int i = 0; i < equat1.n; i++)
                            x2[i] = x1[i];
                        for (int i = 0; i < equat1.n; i++)
                            x1[i] = equat1.x[i];
                        for (int i = 0; i < equat1.n; i++)
                        {
                            double var = 0;
                            for (int j = 0; j < equat1.n; j++)
                            {
                                if (i != j)
                                {
                                    var += (equat1.matr[i, j] * x1[j]);
                                }
                            }

                            equat1.x[i] = (equat1.b[i] - var) / equat1.matr[i, i];
                        }
                        eps = 0;

                        for (int i = 0; i < equat1.n; i++)
                        {
                            eps += Math.Abs(equat1.x[i] - x1[i]);
                        }

                        eps0 = 0;
                        for (int i = 0; i < equat1.n; i++)
                        {
                            eps0 += Math.Abs(x1[i] - x2[i]);
                        }
                        eps0 *= equat1.norm / (1 - equat1.norm);
                    } while (eps > eps0);
                    return equat1.x;
                    }


            }
            else
            {
                double[] x1 = new double[equat.n];
                double eps = 0, eps0 = Math.Abs((1 - equat.norm) / equat.norm) * 0.0001;

                for (int i = 0; i < equat.n; i++)
                {
                    equat.x[i] = equat.b[i];
                }

                do
                {
                    for (int i = 0; i < equat.n; i++)
                        x1[i] = equat.x[i];
                    for (int i = 0; i < equat.n; i++)
                    {
                        double var = 0;
                        for (int j = 0; j < equat.n; j++)
                        {
                            if (i != j)
                            {
                                var += (equat.matr[i, j] * x1[j]);
                            }
                        }

                        equat.x[i] = (equat.b[i] - var) / equat.matr[i, i];
                    }
                    eps = 0;

                    for (int i = 0; i < equat.n; i++)
                    {
                        eps += Math.Abs(equat.x[i] - x1[i]);
                    }
                } while (eps > eps0);
                return equat.x;
            }
        }
    }
}
