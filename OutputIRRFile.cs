﻿
// This class extracts observed vs theoretical ion ratio information from sequest
// .out files (the n/m value in the peptide hit lines) and places it into its own
// separate results file (rootname_IRR.txt)
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
    internal class OutputIRRFile
    {
        private readonly string m_RootName;
        private readonly string m_OutputDirectory;
        private FileStream m_FileStream;
        private StreamWriter m_StreamWriter;
        private int m_CachedScanNum;
        private int m_CachedCS;
        private int m_cachedRankXC;
        private readonly List<IRREntry> m_DataList;
        private readonly string m_Ext = "_IRR.txt";

        public OutputIRRFile(string rootFileName, string OutputDirectory)
        {
            m_RootName = rootFileName;
            m_OutputDirectory = OutputDirectory;
            m_DataList = new List<IRREntry>();
        }

        public void MakeIRREntry(int ScanNum, int ChargeState, int RankXC, int ObsIonCount, int TheoIonCount)
        {
            if (ChargeState != m_CachedCS)
            {
                m_cachedRankXC = 0;
            }
            else if (ScanNum != m_CachedScanNum)
            {
                m_cachedRankXC = 0;
            }

            if (RankXC != m_cachedRankXC)
            {
                m_DataList.Add(new IRREntry(ScanNum, ChargeState, RankXC, ObsIonCount, TheoIonCount));
                m_CachedCS = ChargeState;
                m_CachedScanNum = ScanNum;
                m_cachedRankXC = RankXC;
            }
        }

        protected void WriteEntries(List<IRREntry> dataList)
        {
            dataList.Sort(new RecordIndexComparer());
            var fi = new FileInfo(OutputFilePath);
            if (fi.Exists)
                fi.Delete();
            m_FileStream = new FileStream(OutputFilePath, FileMode.CreateNew);
            m_StreamWriter = new StreamWriter(m_FileStream);
            string headerLine = "Scannum" + ControlChars.Tab.ToString() + "CS" + ControlChars.Tab.ToString() + "RankXc" + ControlChars.Tab.ToString() + "ObservedIons" + ControlChars.Tab.ToString() + "PossibleIons" + ControlChars.Tab.ToString();
            m_StreamWriter.WriteLine(headerLine);
            foreach (var entry in dataList)
            {
                string outputLine = entry.ScanNumber.ToString() + ControlChars.Tab.ToString() + entry.ChargeState.ToString() + ControlChars.Tab.ToString() + entry.RankXc.ToString() + ControlChars.Tab.ToString() + entry.ObsIons.ToString() + ControlChars.Tab.ToString() + entry.PossIons.ToString() + ControlChars.Tab.ToString();
                m_StreamWriter.WriteLine(outputLine);
            }

            m_StreamWriter.Close();
            m_StreamWriter = null;
            m_FileStream.Close();
            m_FileStream = null;
        }

        public void CloseIRRWriter()
        {
            WriteEntries(m_DataList);
        }

        public string OutputFilePath => Path.Combine(m_OutputDirectory, m_RootName + m_Ext);

        public struct IRREntry
        {
            public IRREntry(int ScanNumber, int ChargeState, int RankXc, int ObsIons, int PossIons)
            {
                this.ScanNumber = ScanNumber;
                this.ChargeState = ChargeState;
                this.RankXc = RankXc;
                this.ObsIons = ObsIons;
                this.PossIons = PossIons;
            }

            public int ScanNumber;
            public int ChargeState;
            public int RankXc;
            public int ObsIons;
            public int PossIons;
        }

        private class RecordIndexComparer : IComparer<IRREntry>
        {
            public int Compare(IRREntry x, IRREntry y)
            {
                if (x.ScanNumber > y.ScanNumber)
                {
                    return 1;
                }
                else if (x.ScanNumber < y.ScanNumber)
                {
                    return -1;
                }
                else if (x.ChargeState > y.ChargeState)
                {
                    return 1;
                }
                else if (x.ChargeState < y.ChargeState)
                {
                    return -1;
                }
                else if (x.RankXc > y.RankXc)
                {
                    return 1;
                }
                else if (x.RankXc < y.RankXc)
                {
                    return -1;
                }
                else
                {
                    return 0;
                }
            }
        }
    }
}