using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection.Metadata.Ecma335;
using System.Text;

namespace LightningPCA
{
    class ProgramLowMem
    {
        static string[] HELP = { "--help", "/?", "-h" };

        static void MainLowMem(string[] args)
        {
            if (args.Length == 0 || args.Any(o => HELP.Contains(o)))
            {
                //show help

                return;
            }

            var path = args[0];

            var isPathRooted = Path.IsPathRooted(path);
            if (!isPathRooted)
                path = Path.Combine(Environment.CurrentDirectory, path);

            var dtBegin = DateTime.UtcNow;

            int m, n;
            if (ScanFile(path, out m, out n)) return;

            //goto CLOSE;

            var binCachePath = path + "_bin.tmp";
            if (CreateBinCache(path, binCachePath)) return;

            var binTCachePath = path + "_binT.tmp";
            if (CreateBinTCache(binCachePath, binTCachePath, m, n)) return;

            var meanCachePath = path + "_mean.tmp";
            CalculateMeans(binCachePath, meanCachePath, m, n);

            Log("Calculating 1st eigenvector...");
            var x0 = new double[n];
            FillRandoms(x0);
            Normalize(x0);

            var fsBin = File.OpenRead(binCachePath);
            var brBin = new BinaryReader(fsBin);
            var fsBinT = File.OpenRead(binTCachePath);
            var brBinT = new BinaryReader(fsBinT);
            var fsMean = File.OpenRead(meanCachePath);
            var brMean = new BinaryReader(fsMean);
            double sum_r;
            var tmp1 = new double[m];
            var x1 = new double[n];
            for (int k = 1; ; k++)
            {
                Log("eiv 1 itr " + k + "...");

                //!m1=1/m C'C x0

                //tmp1 = C x0
                brBin.BaseStream.Position = 0;
                for (int i = 0; i < m; i++)
                {
                    sum_r = 0;
                    for (int j = 0; j < n; j++)
                    {
                        sum_r += brBin.ReadByte() * x0[j];
                    }

                    tmp1[i] = sum_r;
                }

                //tmp2 = C' tmp1
                brBinT.BaseStream.Position = 0;
                for (int j = 0; j < n; j++)
                {
                    sum_r = 0;
                    for (int i = 0; i < m; i++)
                    {
                        //brBin.BaseStream.Position = i * n + j;
                        sum_r += brBinT.ReadByte() * tmp1[i];
                    }

                    x1[j] = sum_r;
                }

                //m1 = tmp2 / m
                for (int j = 0; j < n; j++)
                {
                    x1[j] = x1[j] / m;
                }

                //m2=uu'x0

                //sum_r = u'x0
                sum_r = 0;
                brMean.BaseStream.Position = 0;
                for (int j = 0; j < n; j++)
                {
                    sum_r += x0[j] * brMean.ReadDouble();
                }

                //m2 = u sum_r
                //result=m1-m2
                brMean.BaseStream.Position = 0;
                for (int j = 0; j < n; j++)
                {
                    x1[j] = x1[j] - brMean.ReadDouble() * sum_r;
                }

                // lamda = x1x0/x0x0
                sum_r = VectorDotProduct(x0, x1);

                //error=|A x - lamda x|/lamda
                double error = 0;
                for (int j = 0; j < n; j++)
                {
                    var a = x1[j] - sum_r * x0[j];
                    error += a * a;
                }
                error = Math.Sqrt(error) / sum_r;

                Normalize(x1);

                x1.CopyTo(x0, 0);

                Log("error " + error);

                if (error < 5e-2)
                    break;
            }

            var lamda1 = sum_r;

            var alpha1CachePath = path + "_alpha1.tmp";
            if (File.Exists(alpha1CachePath))
            {
                File.Delete(alpha1CachePath);
            }
            var fsAlpha1 = File.Create(alpha1CachePath);
            var bwAlpha1 = new BinaryWriter(fsAlpha1);
            for (int i = 0; i < x1.Length; i++)
            {
                bwAlpha1.Write(x1[i]);
            }
            bwAlpha1.Close();
            fsAlpha1.Close();

            Log("Calculating 2nd eigenvector...");
            fsAlpha1 = File.OpenRead(alpha1CachePath);
            var brAlpha1 = new BinaryReader(fsAlpha1);
            FillRandoms(x0);
            Normalize(x0);
            for (int k = 1; ; k++)
            {
                Log("eiv 2 itr " + k + "...");

                //!m1=1/m C'C x0

                //tmp1 = C x0
                brBin.BaseStream.Position = 0;
                for (int i = 0; i < m; i++)
                {
                    sum_r = 0;
                    for (int j = 0; j < n; j++)
                    {
                        sum_r += brBin.ReadByte() * x0[j];
                    }

                    tmp1[i] = sum_r;
                }

                //tmp2 = C' tmp1
                brBinT.BaseStream.Position = 0;
                for (int j = 0; j < n; j++)
                {
                    sum_r = 0;
                    for (int i = 0; i < m; i++)
                    {
                        //brBin.BaseStream.Position = i * n + j;
                        sum_r += brBinT.ReadByte() * tmp1[i];
                    }

                    x1[j] = sum_r;
                }

                //m1 = tmp2 / m
                for (int j = 0; j < n; j++)
                {
                    x1[j] = x1[j] / m;
                }

                //m2=uu'x0

                //sum_r = u'x0
                sum_r = 0;
                brMean.BaseStream.Position = 0;
                for (int j = 0; j < n; j++)
                {
                    sum_r += x0[j] * brMean.ReadDouble();
                }

                //m2 = u sum_r
                brMean.BaseStream.Position = 0;
                for (int j = 0; j < n; j++)
                {
                    x1[j] = x1[j] - brMean.ReadDouble() * sum_r;
                }

                //m3=lamda1 alpha1 alpha1' x0
                //alpha1' x0
                sum_r = 0;
                brAlpha1.BaseStream.Position = 0;
                for (int j = 0; j < n; j++)
                {
                    sum_r += brAlpha1.ReadDouble() * x0[j];
                }
                //m3
                brAlpha1.BaseStream.Position = 0;
                for (int j = 0; j < n; j++)
                {
                    x1[j] = x1[j] - lamda1 * sum_r * brAlpha1.ReadDouble();
                }

                //result=m1-m2-m3

                //2nd eigen value
                sum_r = VectorDotProduct(x1, x0);

                //error=|A x - lamda x|/lamda
                double error = 0;
                for (int j = 0; j < n; j++)
                {
                    var a = x1[j] - sum_r * x0[j];
                    error += a * a;
                }
                error = Math.Sqrt(error) / sum_r;

                Normalize(x1);

                x1.CopyTo(x0, 0);

                Log("error " + error);

                if (error < 5e-2)
                    break;
            }

            //get 2D data
            Log("Mapping to 2D...");
            var resultPath = path + "_2D_PCA_RESULT.txt";
            if (File.Exists(resultPath))
            {
                File.Delete(resultPath);
            }
            var swResult = File.CreateText(resultPath);

            //m mapper = (c-eu')mapper = c mapper - e u' mapper

            //u' mapper
            double tmp_r1 = 0;
            double tmp_r2 = 0;
            brMean.BaseStream.Position = 0;
            brAlpha1.BaseStream.Position = 0;
            for (int j = 0; j < n; j++)
            {
                sum_r = brMean.ReadDouble();
                tmp_r1 += sum_r * brAlpha1.ReadDouble();
                tmp_r2 += sum_r * x1[j];
            }

            //c mapper
            var data2d1 = new double[m];
            var data2d2 = new double[m];
            brBin.BaseStream.Position = 0;
            for (int i = 0; i < m; i++)
            {
                brAlpha1.BaseStream.Position = 0;
                for (int j = 0; j < n; j++)
                {
                    var cij = brBin.ReadByte();
                    data2d1[i] += cij * brAlpha1.ReadDouble();
                    data2d2[i] += cij * x1[j];
                }

                data2d1[i] = data2d1[i] - tmp_r1;
                data2d2[i] = data2d2[i] - tmp_r2;

                swResult.WriteLine(data2d1[i] + "\t" + data2d2[i]);
            }

            //m mapper = c mapper - e u' mapper
            //for (int i = 0; i < m; i++)
            //{
            //    data2d1[i] = data2d1[i] - tmp_r1;
            //    data2d2[i] = data2d2[i] - tmp_r2;
            //}

            swResult.Close();

            //CLOSE:

            Log("Finished. Time spent: " + (DateTime.UtcNow - dtBegin).TotalMinutes + " minute(s).");
            Log("Press any key to continue...");
            Console.ReadKey();
        }

        private static void Log(string text)
        {
            Console.WriteLine(DateTime.Now.ToString("O") + ": " + text);
        }

        private static void FillRandoms(double[] x0)
        {
            var r = new Random();
            for (int j = 0; j < x0.Length; j++)
            {
                x0[j] = r.NextDouble();
            }
        }

        private static double VectorDotProduct(double[] x0, double[] x1)
        {
            double product = 0;
            for (int i = 0; i < x0.Length; i++)
            {
                product += x0[i] * x1[i];
            }

            return product;
        }

        private static void Normalize(double[] x0)
        {
            double mode = 0;
            for (int j = 0; j < x0.Length; j++)
            {
                mode += x0[j] * x0[j];
            }

            mode = Math.Sqrt(mode);
            for (int j = 0; j < x0.Length; j++)
            {
                x0[j] = x0[j] / mode;
            }
        }

        private static void CalculateMeans(string binCachePath, string meanCachePath, int m, int n)
        {
            Log("Calculating means...");
            var x0 = new double[n];
            var fsBin = File.OpenRead(binCachePath);
            var br = new BinaryReader(fsBin);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    x0[j] += br.ReadByte();
                }
            }

            if (File.Exists(meanCachePath))
            {
                File.Delete(meanCachePath);
            }
            var fsMean = File.Create(meanCachePath);
            var bw = new BinaryWriter(fsMean);
            for (int j = 0; j < n; j++)
            {
                bw.Write(x0[j] / m);
            }

            bw.Close();
            fsMean.Close();

            br.Close();
            fsBin.Close();
        }

        private static bool CreateBinCache(string path, string binCachePath)
        {
            //Log("Preparing disk cache...");

            Log("Building matrix cache file " + binCachePath);
            if (File.Exists(binCachePath))
            {
                File.Delete(binCachePath);
            }
            var fsBin = File.Create(binCachePath);
            var bw = new BinaryWriter(fsBin);

            var sr = File.OpenText(path);
            var lineNum = 1;
            var colNum = 1;
            Char c = (char)0;
            StringBuilder sb = new StringBuilder();
            byte element;
            while (!sr.EndOfStream)
            {
                c = (char)sr.Read();
                switch (c)
                {
                    case '\r':
                    case '\n':
                        if (c == '\r')
                        {
                            Char peek = (char)sr.Peek();
                            if (peek == '\n')
                                sr.Read();
                        }

                        if (sb.Length != 0)
                        {
                            //colNum++;

                            var tryParse = byte.TryParse(sb.ToString(), out element);
                            if (!tryParse)
                            {
                                Log(typeof(byte) + " conversion failed. Line: " + lineNum + " Column: " +
                                                  colNum);
                                return true;
                            }

                            bw.Write(element);

                            sb.Clear();
                        }

                        //bw.Write();

                        lineNum++;
                        colNum = 1;
                        break;

                    case '\t':
                    case ' ':
                        if (sb.Length != 0)
                        {
                            var tryParse = byte.TryParse(sb.ToString(), out element);
                            if (!tryParse)
                            {
                                Log(typeof(byte) + "Integer conversion failed. Line: " + lineNum +
                                                  " Column: " +
                                                  colNum);
                                return true;
                            }

                            bw.Write(element);

                            sb.Clear();

                            colNum++;
                        }

                        break;

                    default:
                        sb.Append(c);
                        break;
                }
            }

            if (sb.Length != 0)
            {
                var tryParse = byte.TryParse(sb.ToString(), out element);
                if (!tryParse)
                {
                    Log(typeof(byte) + "Integer conversion failed. Line: " + lineNum + " Column: " +
                                      colNum);
                    return true;
                }

                bw.Write(element);
            }

            bw.Close();
            sr.Close();
            fsBin.Close();
            return false;
        }

        private static bool CreateBinTCache(string binCachePath, string binTCachePath, int m, int n)
        {
            //Log("Preparing disk cache...");

            Log("Building transposed matrix cache file " + binTCachePath);
            if (File.Exists(binTCachePath))
            {
                File.Delete(binTCachePath);
            }
            var fsBinT = File.Create(binTCachePath);
            var bw = new BinaryWriter(fsBinT);

            var fs = File.OpenRead(binCachePath);
            var brBin = new BinaryReader(fs);

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < m; i++)
                {
                    brBin.BaseStream.Position = i * n + j;
                    bw.Write(brBin.ReadByte());
                }
            }

            bw.Close();
            fsBinT.Close();
            brBin.Close();
            fs.Close();
            return false;
        }

        private static bool ScanFile(string path, out int m, out int n)
        {
            m = 0;
            n = 0;

            FileStream fs;
            try
            {
                fs = File.OpenRead(path);
            }
            catch (Exception e)
            {
                Log(e.Message);
                Log(e.StackTrace);
                return true;
            }

            Log("Scaning " + path);

            var lineElementCounts = new List<int>();
            var elementCount = 0;
            var isElement = false;
            Char c = (char)0;

            var sr = new StreamReader(fs);
            //long totalCount = 0;
            while (!sr.EndOfStream)
            {
                c = (char)sr.Read();
                switch (c)
                {
                    case '\r':
                    case '\n':
                        if (c == '\r')
                        {
                            Char peek = (char)sr.Peek();
                            if (peek == '\n')
                                sr.Read();
                        }

                        isElement = false;

                        lineElementCounts.Add(elementCount);
                        elementCount = 0;

                        break;

                    case '\t':
                    case ' ':
                        isElement = false;
                        break;

                    default:
                        if (!isElement)
                        {
                            elementCount++;
                            //totalCount++;
                        }

                        isElement = true;
                        break;
                }
            }
            //Log(totalCount);

            sr.Close();
            fs.Close();

            if (c != 0 && c != '\r' && c != '\n')
                lineElementCounts.Add(elementCount);

            if (lineElementCounts.Count == 0)
            {
                Log("File is empty.");
                return true;
            }

            var firstLineElementCount = lineElementCounts[0];

            if (firstLineElementCount == 0)
            {
                Log("First line has zero element.");
                return true;
            }

            for (var i = 1; i < lineElementCounts.Count; i++)
            {
                if (lineElementCounts[i] != firstLineElementCount)
                {
                    Log("Line " + (i + 1) + " has " + lineElementCounts[i] + " elements but line 1 has " +
                                      firstLineElementCount);
                    return true;
                }
            }

            m = lineElementCounts.Count;
            n = firstLineElementCount;

            Log(m + " sample(s) and " + n + " metric(s) found.");

            return false;
        }
    }
}