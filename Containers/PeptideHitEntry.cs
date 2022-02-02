using System;
using System.Collections.Generic;
using System.Text;

namespace SequestResultsProcessor.Containers
{
    // This class stores information on individual peptide hits and is used in conjunction
    // with the SequestFileExtractor
    //
    // -------------------------------------------------------------------------------
    // Written by Ken Auberry for the Department of Energy (PNNL, Richland, WA)
    // Coding for this class began on 19 November 2004
    // -------------------------------------------------------------------------------
    //
    // Licensed under the Apache License, Version 2.0; you may not use this file except
    // in compliance with the License.  You may obtain a copy of the License at
    // http://www.apache.org/licenses/LICENSE-2.0

    internal class PeptideHitEntry
    {
        private string mReference;
        private SortedList<int, string> mMultiProteinEntries;
        private int mMultiProteinID;
        private static readonly CleavageStateCalculator mCleavageStateCalculator = new();

        // ReSharper disable once UnusedMember.Global
        public void Clear()
        {
            HitNum = 0;
            MH = 0.0d;
            XCorr = 0.0d;
            DelCn = 0.0d;
            Sp = 0.0d;
            mReference = "";
            MultiProteinCount = 0;
            Peptide = "";
            DelCn2 = 0.0d;
            RankSp = 0;
            RankXc = 0;
            DelM = 0.0d;
            XcRatio = 0.0d;
            // mPassFilt = 0
            // mMScore = 0.0
            NumTrypticEnds = 0;
            ObsIons = 0;
            PossIons = 0;
            ChargeState = 0;
            ScanCount = 0;
            StartScanNum = 0;
            EndScanNum = 0;
            // mObsMH = 0.0

            mMultiProteinEntries.Clear();
        }

        // ReSharper disable once UnusedMember.Global
        public string IDKey => GetIDKey();

        private string GetIDKey()
        {
            var idKey = new StringBuilder();

            idKey.Append(StartScanNum.ToString("000000"));
            idKey.Append(".");
            idKey.Append(EndScanNum.ToString("000000"));
            idKey.Append(".");
            idKey.Append(ChargeState.ToString("00"));
            idKey.Append(".");
            idKey.Append(HitNum.ToString("000"));

            return idKey.ToString();
        }

        public void AddMultiProteinRef(string refName)
        {
            mMultiProteinEntries ??= new SortedList<int, string>();

            mMultiProteinID++;
            mMultiProteinEntries.Add(mMultiProteinID, refName);
        }

        public void CalculateScoreComponents()
        {
            NumTrypticEnds = CountTrypticEnds(Peptide);
        }

        // ReSharper disable once IdentifierTypo
        // ReSharper disable once UnusedMember.Global
        public int CalculatePassFilt()
        {
            var myScore = CalculateFilterScore(XCorr, DelCn2, NumTrypticEnds, ChargeState);
            var optFilterScore = 0.1d + 0.03d * RankXc;
            return myScore >= optFilterScore ? 1 : 0;
        }

        private double CalculateFilterScore(double xcorr, double deltaCn, int numTrypticEnds, int chargeState)
        {
            double[] trypRelScore;

            switch (chargeState)
            {
                case 1:
                    trypRelScore = new[] { 0.35d, 0.9d, 1.04d };
                    var myScore = trypRelScore[numTrypticEnds] * (1.03d / (1d + Math.Exp((1.49d - xcorr) / 0.25d))) * (0.98d / (1d + Math.Exp((0.07d - deltaCn) / 0.085d)));
                    return myScore - 0.1d;

                case 2:
                    trypRelScore = new[] { 0.31d, 0.76d, 0.98d };
                    return trypRelScore[numTrypticEnds] * (1.03d / (1d + Math.Exp((2.44d - xcorr) / 0.524d))) * (1.02d / (1d + Math.Exp((0.09d - deltaCn) / 0.03d)));

                case 3:
                    trypRelScore = new[] { 0.21d, 0.62d, 1.04d };
                    return trypRelScore[numTrypticEnds] * (1.01d / (1d + Math.Exp((3.2d - xcorr) / 0.48d))) * (1.04d / (1d + Math.Exp((0.11d - deltaCn) / 0.06d)));
            }

            return 0;
        }

        private int CountTrypticEnds(string peptideSeq)
        {
            return mCleavageStateCalculator.CountTrypticEnds(peptideSeq);
        }

        public SortedList<string, string> ExportMultiProteinXref(ResultsStorage.OutputTypeList ResultsType)
        {
            SortedList<string, string> proteinList = null;
            var identifierSB = new StringBuilder();
            var lineSB = new StringBuilder();
            if (MultiProteinCount > 0)
            {
                if (ResultsType == ResultsStorage.OutputTypeList.Syn || ResultsType == ResultsStorage.OutputTypeList.FHT && RankXc == 1)
                {
                    proteinList = new SortedList<string, string>();
                    foreach (var objEntry in mMultiProteinEntries)
                    {
                        identifierSB.Clear();
                        identifierSB.Append(StartScanNum.ToString("##0000"));
                        identifierSB.Append(".");
                        identifierSB.Append(EndScanNum.ToString("##0000"));
                        identifierSB.Append(".");
                        identifierSB.Append(ChargeState.ToString());
                        identifierSB.Append(".");
                        identifierSB.Append(RankXc.ToString());
                        identifierSB.Append(".");
                        identifierSB.Append(objEntry.Key.ToString());
                        lineSB.Clear();
                        lineSB.Append(RankXc.ToString());
                        lineSB.Append('\t');
                        lineSB.Append(StartScanNum.ToString());
                        lineSB.Append('\t');
                        lineSB.Append(ChargeState.ToString());
                        lineSB.Append('\t');
                        lineSB.Append(objEntry.Key.ToString());
                        lineSB.Append('\t');
                        lineSB.Append(objEntry.Value);
                        proteinList.Add(identifierSB.ToString(), lineSB.ToString());
                    }
                }
            }

            return proteinList;
        }

        public SortedList<int, string> ExportContents(bool expandMultiProteinEntries)
        {
            var startColumns = new StringBuilder(150);
            var endColumns = new StringBuilder(150);

            const char delimiter = '\t';

            // Tracks protein names:
            // keys are 0, 1, 2, etc.
            // Values are protein names
            var exportList = new SortedList<int, string>();

            CalculateScoreComponents();

            // startColumns and endColumns contain the constant parts of the peptide result
            // each multiorf reference will be stuck in the middle during the integration
            // phase at the end of the function

            startColumns.Append(HitNum.ToString());
            startColumns.Append(delimiter);
            startColumns.Append(StartScanNum.ToString("##0000"));
            startColumns.Append(delimiter);
            startColumns.Append(ScanCount.ToString());
            startColumns.Append(delimiter);
            startColumns.Append(ChargeState.ToString());
            startColumns.Append(delimiter);
            startColumns.Append(Math.Round(MH, 5).ToString("#####0.00000"));
            startColumns.Append(delimiter);
            startColumns.Append(Math.Round(XCorr, 4).ToString("##0.0000"));
            startColumns.Append(delimiter);
            startColumns.Append(Math.Round(DelCn, 4).ToString("##0.0000"));
            startColumns.Append(delimiter);
            startColumns.Append(Math.Round(Sp, 1).ToString("######0.0"));
            startColumns.Append(delimiter);

            if (MultiProteinCount > 0)
            {
                endColumns.Append(MultiProteinCount.ToString("+0"));
            }
            else
            {
                endColumns.Append("0");
            }

            endColumns.Append(delimiter);
            endColumns.Append(Peptide);
            endColumns.Append(delimiter);
            endColumns.Append(Math.Round(DelCn2, 4).ToString("##0.0000"));
            endColumns.Append(delimiter);
            endColumns.Append(RankSp.ToString());
            endColumns.Append(delimiter);
            endColumns.Append(RankXc.ToString());
            endColumns.Append(delimiter);
            endColumns.Append(Math.Round(DelM, 5).ToString("##0.00000"));
            endColumns.Append(delimiter);
            endColumns.Append(Math.Round(XcRatio, 3).ToString("0.000"));
            endColumns.Append(delimiter);

            // PassFilt and MScore are Legacy columns
            // endColumns.Append(PassFilt.ToString)
            // endColumns.Append(delimiter)
            // If .MScore = 10 OrElse .MScore = 0 Then
            // endColumns.Append(MScore.ToString("0"))
            // Else
            // endColumns.Append(Math.Round(.MScore, 2).ToString("##0.00"))
            // End If
            // endColumns.Append(delimiter)

            endColumns.Append(ObsIons.ToString());
            endColumns.Append(delimiter);
            endColumns.Append(PossIons.ToString());
            endColumns.Append(delimiter);
            endColumns.Append(NumTrypticEnds.ToString());
            endColumns.Append(delimiter);
            endColumns.Append(DelMPPM.ToString("##0.0000"));

            exportList.Add(0, startColumns + Reference + delimiter + endColumns);
            if (expandMultiProteinEntries && mMultiProteinEntries != null)
            {
                // Keys in mMultiProteinEntries are MultiProteinID #, values are protein names
                foreach (var objEntry in mMultiProteinEntries)
                {
                    if ((Reference ?? "") != (objEntry.Value ?? ""))
                    {
                        exportList.Add(objEntry.Key, startColumns + objEntry.Value + delimiter + endColumns);
                    }
                }
            }

            return exportList;
        }

        public int HitNum { get; set; }
        public int EndScanNum { get; set; }
        public int StartScanNum { get; set; }
        public int ScanCount { get; set; }
        public int ChargeState { get; set; }
        public double MH { get; set; }
        public double XCorr { get; set; }
        public double DelCn { get; set; }
        public double Sp { get; set; }

        public string Reference
        {
            get => mReference;

            set
            {
                mReference = value;
                AddMultiProteinRef(value);
            }
        }

        public int MultiProteinCount { get; set; }
        public string Peptide { get; set; }
        public double DelCn2 { get; set; }
        public int RankSp { get; set; }
        public int RankXc { get; set; }
        public double DelM { get; set; }
        public double XcRatio { get; set; } = 1.0d;

        public int ObsIons { get; set; }
        public int PossIons { get; set; }
        public int NumTrypticEnds { get; set; }
        public double DelMPPM { get; set; }
    }
}