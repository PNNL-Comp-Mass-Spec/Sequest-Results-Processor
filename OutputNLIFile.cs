﻿
// This class outputs the information generated by CalcNeutralLosses into its
// own separate results file (RootName_NLI.txt)
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
using System.Globalization;
using System.IO;

namespace SequestResultsProcessor
{
    internal class OutputNLIFile
    {
        private readonly string mRootName;
        private readonly string mOutputDirectory;
        private FileStream mFileStream;
        private StreamWriter mStreamWriter;
        private readonly List<NLIEntry> mDataList;
        private readonly string mExt = "_NLI.txt";

        public OutputNLIFile(string rootFileName, string OutputDirectory)
        {
            mRootName = rootFileName;
            mOutputDirectory = OutputDirectory;
            mDataList = new List<NLIEntry>();
        }

        private void WriteEntries(List<NLIEntry> dataList)
        {
            dataList.Sort(new NLIIndexComparer());
            var fi = new FileInfo(OutputFilePath);
            if (fi.Exists)
                fi.Delete();
            mFileStream = new FileStream(OutputFilePath, FileMode.CreateNew);
            mStreamWriter = new StreamWriter(mFileStream);
            var headerLine = "Scannum" + '\t' + "NL1_Intensity" + '\t' + "NL2_Intensity" + '\t' + "NL3_Intensity" + '\t';
            mStreamWriter.WriteLine(headerLine);
            foreach (var entry in dataList)
            {
                var outputLine = entry.ScanNumber.ToString() + '\t' + entry.NeutralLossTxt;
                mStreamWriter.WriteLine(outputLine);
            }

            mStreamWriter.Close();
            mStreamWriter = null;
            mFileStream.Close();
            mFileStream = null;
        }

        public void MakeNLIEntry(int ScanNum, NeutralLossList NeutralLosses)
        {
            mDataList.Add(new NLIEntry(ScanNum, NeutralLosses));
        }

        // ReSharper disable once UnusedMember.Global
        public void CloseNLIWriter()
        {
            WriteEntries(mDataList);
        }

        public string OutputFilePath => Path.Combine(mOutputDirectory, mRootName + mExt);

        public struct NLIEntry
        {
            public NLIEntry(int ScanNumber, NeutralLossList NeutralLosses)
            {
                this.ScanNumber = ScanNumber;
                NeutralLossTxt = NeutralLosses.NL1Intensity.ToString(CultureInfo.InvariantCulture) + '\t' + NeutralLosses.NL2Intensity + '\t' + NeutralLosses.NL3Intensity + '\t';
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