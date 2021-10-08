﻿
// This class outputs the information generated by CalcNeutralLosses into its
// own separate results file (rootname_NLI.txt)
//
// -------------------------------------------------------------------------------
// Written by Ken Auberry for the Department of Energy (PNNL, Richland, WA)
// Coding for this class began on 19 November 2004
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0

using System.Collections.Generic;
using System.IO;

namespace SequestResultsProcessor
{
    internal class OutputNLIFile
    {
        private readonly string m_RootName;
        private readonly string m_OutputDirectory;
        private FileStream m_FileStream;
        private StreamWriter m_StreamWriter;
        private readonly List<NLIEntry> m_DataList;
        private readonly string m_Ext = "_NLI.txt";

        public OutputNLIFile(string rootFileName, string OutputDirectory)
        {
            m_RootName = rootFileName;
            m_OutputDirectory = OutputDirectory;
            m_DataList = new List<NLIEntry>();
        }

        private void WriteEntries(List<NLIEntry> dataList)
        {
            dataList.Sort(new NLIIndexComparer());
            var fi = new FileInfo(OutputFilePath);
            if (fi.Exists)
                fi.Delete();
            m_FileStream = new FileStream(OutputFilePath, FileMode.CreateNew);
            m_StreamWriter = new StreamWriter(m_FileStream);
            string headerLine = "Scannum" + '\t' + "NL1_Intensity" + '\t' + "NL2_Intensity" + '\t' + "NL3_Intensity" + '\t';
            m_StreamWriter.WriteLine(headerLine);
            foreach (var entry in dataList)
            {
                string outputLine = entry.ScanNumber.ToString() + '\t' + entry.NeutralLossTxt;
                m_StreamWriter.WriteLine(outputLine);
            }

            m_StreamWriter.Close();
            m_StreamWriter = null;
            m_FileStream.Close();
            m_FileStream = null;
        }

        public void MakeNLIEntry(int ScanNum, NeutralLossList NeutralLosses)
        {
            m_DataList.Add(new NLIEntry(ScanNum, NeutralLosses));
        }

        // ReSharper disable once UnusedMember.Global
        public void CloseNLIWriter()
        {
            WriteEntries(m_DataList);
        }

        public string OutputFilePath => Path.Combine(m_OutputDirectory, m_RootName + m_Ext);

        public struct NLIEntry
        {
            public NLIEntry(int ScanNumber, NeutralLossList NeutralLosses)
            {
                this.ScanNumber = ScanNumber;
                NeutralLossTxt = NeutralLosses.NL1Intensity.ToString() + '\t' + NeutralLosses.NL2Intensity.ToString() + '\t' + NeutralLosses.NL3Intensity.ToString() + '\t';
            }

            public int ScanNumber;
            public string NeutralLossTxt;
        }

        private class NLIIndexComparer : IComparer<NLIEntry>
        {
            public int Compare(NLIEntry x, NLIEntry y)
            {
                if (x.ScanNumber > y.ScanNumber)
                {
                    return 1;
                }

                if (x.ScanNumber < y.ScanNumber)
                {
                    return -1;
                }
                return 0;
            }
        }
    }
}