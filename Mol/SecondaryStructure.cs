/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System;
using System.IO;
using System.Text;

namespace MaxQuant.Mol
{
    [Serializable]
    public class SecondaryStructure
    {
        private const int pssPerLong = 12;
        private const int basis = 32;

        private ulong[] prediction;
        private int len;

        public SecondaryStructure()
        {
        }

        public SecondaryStructure(string predict)
        {
            len = predict.Length;
            int a = len / pssPerLong;
            int b = len % pssPerLong;
            int n = (b == 0) ? a : a + 1;
            prediction = new ulong[n];

            for (int i = 0; i < a; i++)
            {
                prediction[i] = Encode(predict.Substring(pssPerLong * i, pssPerLong));
            }

            if (b > 0)
            {
                prediction[a] = Encode(predict.Substring(pssPerLong * a, len - pssPerLong * a));
            }
        }

        public static SecondaryStructure Read(BinaryReader reader)
        {
            try
            {
                SecondaryStructure result = new SecondaryStructure();
                result.len = reader.ReadInt32();
                result.prediction = new ulong[reader.ReadInt32()];

                for (int i = 0; i < result.prediction.Length; i++)
                {
                    result.prediction[i] = reader.ReadUInt16();
                }
                reader.Close();

                return result;
            }
            catch (Exception)
            {
                return null;
            }
        }

        public void Write(BinaryWriter writer)
        {
            writer.Write(len);
            writer.Write(prediction.Length);
            foreach (ulong u in prediction)
            {
                writer.Write(u);
            }
        }

        public override string ToString()
        {
            StringBuilder b = new StringBuilder();
            for (int i = 0; i < prediction.Length - 1; i++)
            {
                b.Append(Decode(prediction[i], pssPerLong));
            }

            int l = len - (prediction.Length - 1) * pssPerLong;
            b.Append(Decode(prediction[prediction.Length - 1], l));

            return b.ToString();
        }

        private static ulong Encode(string predict)
        {
            ulong result = 0;
            for (int i = 0; i < predict.Length; i++)
            {
                char c = predict[i];
                int x = SecondaryStructureElement.singleLetter.IndexOf(c);
                if (x == -1)
                {
                    x = 3;
                }
                result = (result << 5) + (ulong)x;
            }

            return result;
        }

        private static string Decode(ulong a, int n)
        {
            char[] result = new char[n];
            for (int i = 0; i < n; i++)
            {
                int x = (int)(a % basis);
                a = a >> 5;
                if (x == 3)
                {
                    result[n - 1 - i] = 'X';
                }
                else
                {
                    result[n - 1 - i] = SecondaryStructureElement.singleLetter[x];
                }
            }

            return new string(result);
        }
    }
}