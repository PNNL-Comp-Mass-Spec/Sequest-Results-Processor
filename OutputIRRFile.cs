﻿
// This class extracts observed vs theoretical ion ratio information from sequest
// .out files (the n/m value in the peptide hit lines) and places it into its own
// separate results file (RootName_IRR.txt)
//
// -------------------------------------------------------------------------------
// Written by Ken Auberry for the Department of Energy (PNNL, Richland, WA)
// Coding for this class began in November 2004
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
        private readonly string mRootName;
        private readonly string mOutputDirectory;
        private FileStream mFileStream;
        private StreamWriter mStreamWriter;
        private int mCachedScanNum;
        private int mCachedCS;
        private int mCachedRankXC;
        private readonly List<IRREntry> mDataList;
        private readonly string mExt = "_IRR.txt";

        public OutputIRRFile(string rootFileName, string OutputDirectory)
        {
            mRootName = rootFileName;
            mOutputDirectory = OutputDirectory;
            mDataList = new List<IRREntry>();
        }

        public void MakeIRREntry(int ScanNum, int ChargeState, int RankXC, int ObsIonCount, int TheoreticalIonCount)
        {
            if (ChargeState != mCachedCS)
            {
                mCachedRankXC = 0;
            }
            else if (ScanNum != mCachedScanNum)
            {
                mCachedRankXC = 0;
            }

            if (RankXC != mCachedRankXC)
            {
                mDataList.Add(new IRREntry(ScanNum, ChargeState, RankXC, ObsIonCount, TheoreticalIonCount));
                mCachedCS = ChargeState;
                mCachedScanNum = ScanNum;
                mCachedRankXC = RankXC;
            }
        }

        private void WriteEntries(List<IRREntry> dataList)
        {
            dataList.Sort(new RecordIndexComparer());
            var fi = new FileInfo(OutputFilePath);
            if (fi.Exists)
                fi.Delete();

            mFileStream = new FileStream(OutputFilePath, FileMode.CreateNew);
            mStreamWriter = new StreamWriter(mFileStream);
            var headerLine = "Scannum" + '\t' + "CS" + '\t' + "RankXc" + '\t' + "ObservedIons" + '\t' + "PossibleIons" + '\t';
            mStreamWriter.WriteLine(headerLine);

            foreach (var entry in dataList)
            {
                var outputLine = entry.ScanNumber.ToString() + '\t' + entry.ChargeState + '\t' + entry.RankXc + '\t' + entry.ObsIons + '\t' + entry.PossIons + '\t';
                mStreamWriter.WriteLine(outputLine);
            }

            mStreamWriter.Close();
            mStreamWriter = null;
            mFileStream.Close();
            mFileStream = null;
        }

        public void CloseIRRWriter()
        {
            WriteEntries(mDataList);
        }

        public string OutputFilePath => Path.Combine(mOutputDirectory, mRootName + mExt);

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

                if (x.ScanNumber < y.ScanNumber)
                {
                    return -1;
                }

                if (x.ChargeState > y.ChargeState)
                {
                    return 1;
                }

                if (x.ChargeState < y.ChargeState)
                {
                    return -1;
                }

                if (x.RankXc > y.RankXc)
                {
                    return 1;
                }

                if (x.RankXc < y.RankXc)
                {
                    return -1;
                }

                return 0;
            }
        }
    }
}