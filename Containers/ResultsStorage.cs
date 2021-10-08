// This class stages and stores information extracted by the SequestFileExtractor
// until it is output to the final synopsis and first hits files
//
// -------------------------------------------------------------------------------
// Written by Ken Auberry for the Department of Energy (PNNL, Richland, WA)
// Coding for this class began on 19 November 2004
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//

using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;

namespace SequestResultsProcessor.Containers
{
    /// <summary>
    /// Tracks peptide results from several .Out files
    /// </summary>
    /// <remarks></remarks>
    internal class ResultsStorage
    {
        public enum OutputTypeList
        {
            Syn,
            FHT
        }

        // Keys in this dictionary are StartScan.EndScan.Charge
        private Dictionary<string, ResultsFileEntry> m_Results;
        private ResultsFileEntry m_CachedResultsFileEntry;
        private const string HEADER = "HitNum\tScanNum\tScanCount\tChargeState\tMH\tXCorr\tDelCn\tSp\tReference\tMultiProtein\tPeptide\tDelCn2\tRankSp\tRankXc\tDelM\tXcRatio\tIons_Observed\tIons_Expected\tNumTrypticEnds\tDelM_PPM";

        public int Count
        {
            get
            {
                if (m_Results is null)
                {
                    return 0;
                }

                return m_Results.Count;
            }
        }

        public ResultsStorage()
        {
            // Key/value collection ->
            // key: "StartScan.EndScan.ChargeState"
            // value: ResultsFileEntry entity
            m_Results = new Dictionary<string, ResultsFileEntry>();
        }

        public void AddPeptideResults(double HeaderMass, PeptideHitEntry peptideResults)
        {
            var tmpKey = GetKey(peptideResults.StartScanNum, peptideResults.EndScanNum, peptideResults.ChargeState);
            if (Math.Abs(peptideResults.XCorr) < float.Epsilon)
            {
                return;
            }

            if (m_Results.TryGetValue(tmpKey, out var cachedResultsFileEntry))
            {
                // Make sure m_CachedResultsFileEntry is up-to-date
                if (!ReferenceEquals(m_CachedResultsFileEntry, cachedResultsFileEntry))
                {
                    m_CachedResultsFileEntry = cachedResultsFileEntry;
                }
            }
            else
            {
                m_CachedResultsFileEntry = new ResultsFileEntry(peptideResults.StartScanNum, peptideResults.EndScanNum, peptideResults.ChargeState);
                m_CachedResultsFileEntry.AddHeaderMass(HeaderMass);
                m_Results.Add(tmpKey, m_CachedResultsFileEntry);
            }

            m_CachedResultsFileEntry.AddPeptideResults(peptideResults);
        }

        public void ClearResults()
        {
            m_Results.Clear();
        }

        public void ExportContents(OutputTypeList outputType, double XCorrCutoff, bool ExpandMultiProtein, string ExportFilePath, List<OutputRecordIndex> outputRecordIndexList)
        {
            var proteinXrefFilePath = Path.Combine(Path.GetDirectoryName(ExportFilePath), Path.GetFileNameWithoutExtension(ExportFilePath) + "_prot" + Path.GetExtension(ExportFilePath));
            // two loops
            // outer: each result file in the set
            // inner: each peptide result in a result file

            SortedList<string, string> proteinExportList = null;
            StreamWriter xrefWriter = null;
            var outputWriter = new StreamWriter(ExportFilePath, true, Encoding.ASCII);
            if (outputType == OutputTypeList.FHT)
            {
                xrefWriter = new StreamWriter(proteinXrefFilePath, true, Encoding.ASCII);
                var proteinOutputFileInfo = new FileInfo(proteinXrefFilePath);
                if (proteinOutputFileInfo.Length == 0L)
                {
                    xrefWriter.WriteLine("RankXc" + '\t' + "ScanNum" + '\t' + "ChargeState" + '\t' + "MultiProteinID" + '\t' + "Reference" + '\t');
                }
            }

            var outputFileInfo = new FileInfo(ExportFilePath);
            var currentPosition = outputFileInfo.Length;
            foreach (var resultsFile in m_Results.Values)
            {
                foreach (var peptideHit in resultsFile.PeptideHits.Values)
                {
                    if (outputType == OutputTypeList.FHT & peptideHit.HitNum > 1)
                    {
                        break;
                    }

                    var exportList = peptideHit.ExportContents(ExpandMultiProtein);
                    if (outputType == OutputTypeList.FHT)
                    {
                        proteinExportList = peptideHit.ExportMultiProteinXref(outputType);
                    }

                    var currentMultiProteinID = 0;

                    // Old test (deprecated in Fall 2011)
                    // If peptideHit.XCorr > XCorrCutoff Or peptideHit.PassFilt > 0 Then

                    if (peptideHit.XCorr > XCorrCutoff)
                    {
                        foreach (var peptideLine in exportList.Values)
                        {
                            outputWriter.WriteLine(peptideLine);
                            var tmpLength = peptideLine.Length + outputWriter.NewLine.Length;
                            outputRecordIndexList.Add(new OutputRecordIndex(peptideHit.XCorr, peptideHit.StartScanNum, peptideHit.EndScanNum, peptideHit.ChargeState, peptideHit.HitNum, currentMultiProteinID, currentPosition, tmpLength));
                            currentPosition += tmpLength;
                            currentMultiProteinID += 1;
                        }

                        if (outputType == OutputTypeList.FHT && proteinExportList is object)
                        {
                            foreach (var proteinLine in proteinExportList.Values)
                                xrefWriter.WriteLine(proteinLine);
                        }
                    }
                }
            }

            if (outputWriter is object)
            {
                outputWriter.Close();
            }

            if (xrefWriter is object)
            {
                xrefWriter.Close();
            }
        }

        /// <summary>
        /// Sorts the peptides and writes them to disk
        /// </summary>
        /// <param name="outputFilePath"></param>
        /// <param name="outputRecordList"></param>
        /// <param name="finalOutputPath"></param>
        /// <returns>Dictionary where keys are XCorr threshold and values are the number of peptides with an XCorr over the threshold</returns>
        public Dictionary<int, int> SortPeptides(string outputFilePath, List<OutputRecordIndex> outputRecordList, string finalOutputPath)
        {
            var proteinOutputPath = MakeProteinXrefOutputPath(outputFilePath);
            var finalProteinOutputPath = MakeProteinXrefOutputPath(finalOutputPath);
            var proteinOutputFileInfo = new FileInfo(proteinOutputPath);
            var outputFileInfo = new FileInfo(outputFilePath);
            if (proteinOutputFileInfo.Exists)
            {
                proteinOutputFileInfo.CopyTo(finalProteinOutputPath, true);
            }

            byte[] buffer;
            var fi = new FileInfo(outputFilePath);
            FileStream fs;
            if (fi.Length == 0L)
            {
                fs = File.Create(finalOutputPath);
                fs.Close();
            }
            else
            {
                var rowCount = 1;
                var reUpdateHitNum = new Regex(@"^\d+");
                string inputString;
                string outputString;
                outputRecordList.Sort(new OutputRecordIndexComparer());
                using (var fsInFile = new FileStream(outputFilePath, FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    using (var swOutFile = new StreamWriter(new FileStream(finalOutputPath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                    {
                        swOutFile.WriteLine(HEADER);
                        foreach (var entry in outputRecordList)
                        {
                            buffer = new byte[entry.RecordLength + 1];
                            fsInFile.Seek(entry.StartOffset, SeekOrigin.Begin);
                            fsInFile.Read(buffer, 0, entry.RecordLength);
                            inputString = Encoding.Default.GetString(buffer);

                            // Update the the row number (the first number on the line)
                            outputString = reUpdateHitNum.Replace(inputString, rowCount.ToString());
                            swOutFile.Write(outputString.Trim('\0'));
                            rowCount += 1;
                        }
                    }
                }
            }

            var XCorrStatsGenerator = new XCorrSummaryGenerator(outputRecordList);
            XCorrStatsGenerator.GenerateXCorrSummaries();
            if (outputFileInfo.Exists)
            {
                outputFileInfo.Delete();
            }

            if (proteinOutputFileInfo.Exists)
            {
                proteinOutputFileInfo.Delete();
            }

            return XCorrStatsGenerator.StatsTable;
        }

        private string MakeProteinXrefOutputPath(string outputFilePath)
        {
            return Path.Combine(Path.GetDirectoryName(outputFilePath), Path.GetFileNameWithoutExtension(outputFilePath) + "_prot" + Path.GetExtension(outputFilePath));
        }

        private string GetKey(int startScan, int endScan, int chargeState)
        {
            return startScan.ToString("000000") + "." + endScan.ToString("000000") + "." + chargeState.ToString("00");
        }

        internal class OutputRecordIndex
        {
            public OutputRecordIndex(double XCorr, int StartScanNum, int EndScanNum, int ChargeState, int HitNum, int MultiProteinID, long StartOffset, int RecordLength)
            {
                this.RecordLength = RecordLength;
                this.StartOffset = StartOffset;
                this.XCorr = XCorr;
                this.StartScanNum = StartScanNum;
                this.EndScanNum = EndScanNum;
                this.ChargeState = ChargeState;
                this.HitNum = HitNum;
                this.MultiProteinID = MultiProteinID;
            }

            public double XCorr { get; set; }
            public int StartScanNum { get; set; }
            public int EndScanNum { get; set; }
            public int ChargeState { get; set; }
            public long StartOffset { get; set; }
            public int RecordLength { get; set; }
            public int HitNum { get; set; }
            public int MultiProteinID { get; set; }
        }

        private class OutputRecordIndexComparer : IComparer<OutputRecordIndex>
        {
            public int Compare(OutputRecordIndex x, OutputRecordIndex y)
            {
                if (x.XCorr > y.XCorr)
                {
                    return -1;
                }

                if (x.XCorr < y.XCorr)
                {
                    return 1;
                }

                if (x.StartScanNum > y.StartScanNum)
                {
                    return 1;
                }

                if (x.StartScanNum < y.StartScanNum)
                {
                    return -1;
                }

                if (x.EndScanNum > y.EndScanNum)
                {
                    return 1;
                }

                if (x.EndScanNum < y.EndScanNum)
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

                if (x.HitNum > y.HitNum)
                {
                    return 1;
                }

                if (x.HitNum < y.HitNum)
                {
                    return -1;
                }

                if (x.MultiProteinID > y.MultiProteinID)
                {
                    return 1;
                }

                if (x.MultiProteinID < y.MultiProteinID)
                {
                    return -1;
                }

                return 0;
            }
        }

        internal class XCorrSummaryGenerator
        {
            private readonly List<OutputRecordIndex> m_OutputIndexList;

            public XCorrSummaryGenerator(List<OutputRecordIndex> outputIndexList)
            {
                m_OutputIndexList = outputIndexList;
                StatsTable = new Dictionary<int, int>();
            }

            public Dictionary<int, int> GenerateXCorrSummaries()
            {
                var GT5 = default(int);
                var GT4 = default(int);
                var GT3 = default(int);
                var GT2 = default(int);
                var GT1 = default(int);
                var GT0 = default(int);
                foreach (var entry in m_OutputIndexList)
                {
                    switch (entry.XCorr)
                    {
                        case var @case when @case > 5.0d:
                            {
                                GT5 += 1;
                                GT4 += 1;
                                GT3 += 1;
                                GT2 += 1;
                                GT1 += 1;
                                GT0 += 1;
                                break;
                            }

                        case var case1 when case1 > 4.0d:
                            {
                                GT4 += 1;
                                GT3 += 1;
                                GT2 += 1;
                                GT1 += 1;
                                GT0 += 1;
                                break;
                            }

                        case var case2 when case2 > 3.0d:
                            {
                                GT3 += 1;
                                GT2 += 1;
                                GT1 += 1;
                                GT0 += 1;
                                break;
                            }

                        case var case3 when case3 > 2.0d:
                            {
                                GT2 += 1;
                                GT1 += 1;
                                GT0 += 1;
                                break;
                            }

                        case var case4 when case4 > 1.0d:
                            {
                                GT1 += 1;
                                GT0 += 1;
                                break;
                            }

                        default:
                            {
                                GT0 += 1;
                                break;
                            }
                    }
                }

                {
                    var withBlock = StatsTable;
                    withBlock.Add(5, GT5);
                    withBlock.Add(4, GT4);
                    withBlock.Add(3, GT3);
                    withBlock.Add(2, GT2);
                    withBlock.Add(1, GT1);
                    withBlock.Add(0, GT0);
                }

                return StatsTable;
            }

            private int GetSummaryTableValue(int xcorrThreshold)
            {
                if (StatsTable.TryGetValue(xcorrThreshold, out var peptideCount))
                {
                    return peptideCount;
                }

                return 0;
            }

            /// <summary>
            /// Keys in this dictionary are XCorr threshold; values are the number of peptides with an XCorr over the threshold
            /// </summary>
            /// <returns></returns>
            public Dictionary<int, int> StatsTable { get; private set; }
        }
    }
}