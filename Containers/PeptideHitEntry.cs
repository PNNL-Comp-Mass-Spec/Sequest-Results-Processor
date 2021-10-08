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
        private string m_Reference;
        private SortedList<int, string> m_MultiProteinEntries;
        private int m_MultiProteinID = 0;
        private static CleavageStateCalculator m_CleavageStateCalculator = new CleavageStateCalculator();
        private static DiscriminantCalc _s_DiscriminantCalc;

        public static event dtaLoadProgressEventHandler dtaLoadProgress;

        public delegate void dtaLoadProgressEventHandler(string taskDescription, double fractionDone);

        public static void OnDTALoadUpdate(string taskDescription, double fractiondone)
        {
            dtaLoadProgress?.Invoke(taskDescription, fractiondone);
        }

        public void Clear()
        {
            HitNum = 0;
            MH = 0.0d;
            XCorr = 0.0d;
            DelCn = 0.0d;
            Sp = 0.0d;
            m_Reference = "";
            MultiProteinCount = 0;
            Peptide = "";
            DelCn2 = 0.0d;
            RankSp = 0;
            RankXc = 0;
            DelM = 0.0d;
            XcRatio = 0.0d;
            // m_PassFilt = 0
            // m_MScore = 0.0
            NumTrypticEnds = 0;
            ObsIons = 0;
            PossIons = 0;
            ChargeState = 0;
            ScanCount = 0;
            StartScanNum = 0;
            EndScanNum = 0;
            // m_ObsMH = 0.0

            m_MultiProteinEntries.Clear();
        }

        public string IDKey => GetIDKey();

        private string GetIDKey()
        {
            var sb = new StringBuilder();
            sb.Append(StartScanNum.ToString("000000"));
            sb.Append(".");
            sb.Append(EndScanNum.ToString("000000"));
            sb.Append(".");
            sb.Append(ChargeState.ToString("00"));
            sb.Append(".");
            sb.Append(HitNum.ToString("000"));
            return sb.ToString();
        }

        public void AddMultiProteinRef(string refName)
        {
            if (m_MultiProteinEntries is null)
            {
                m_MultiProteinEntries = new SortedList<int, string>();
            }

            m_MultiProteinID += 1;
            m_MultiProteinEntries.Add(m_MultiProteinID, refName);
        }

        public void CalculateScoreComponents()
        {
            NumTrypticEnds = CountTrypticEnds(Peptide);
        }

        public int CalculatePassFilt()
        {
            double myScore = CalculateFilterScore(XCorr, DelCn2, NumTrypticEnds, ChargeState);
            double optFilterScore = 0.1d + 0.03d * RankXc;
            if (myScore >= optFilterScore)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }

        private double CalculateFilterScore(double dblXCorr, double deltCn, int intNumTrypticEnds, int intChargeState)
        {
            double myScore = 0.0d;
            double[] trypRelScore;
            switch (intChargeState)
            {
                case 1:
                    {
                        trypRelScore = new double[] { 0.35d, 0.9d, 1.04d };
                        myScore = trypRelScore[intNumTrypticEnds] * (1.03d / (1d + Math.Exp((1.49d - dblXCorr) / 0.25d))) * (0.98d / (1d + Math.Exp((0.07d - deltCn) / 0.085d)));
                        myScore = myScore - 0.1d;
                        break;
                    }

                case 2:
                    {
                        trypRelScore = new double[] { 0.31d, 0.76d, 0.98d };
                        myScore = trypRelScore[intNumTrypticEnds] * (1.03d / (1d + Math.Exp((2.44d - dblXCorr) / 0.524d))) * (1.02d / (1d + Math.Exp((0.09d - deltCn) / 0.03d)));
                        break;
                    }

                case 3:
                    {
                        trypRelScore = new double[] { 0.21d, 0.62d, 1.04d };
                        myScore = trypRelScore[intNumTrypticEnds] * (1.01d / (1d + Math.Exp((3.2d - dblXCorr) / 0.48d))) * (1.04d / (1d + Math.Exp((0.11d - deltCn) / 0.06d)));
                        break;
                    }
            }

            return myScore;
        }

        private int CountTrypticEnds(string peptideSeq)
        {
            return m_CleavageStateCalculator.CountTrypticEnds(peptideSeq);
        }

        public int ExportLength { get; private set; } = 0;

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
                    foreach (KeyValuePair<int, string> objEntry in m_MultiProteinEntries)
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

        public SortedList<int, string> ExportContents(bool ExpandMultiProteinEntries)
        {
            StringBuilder sbFront;
            StringBuilder sbRear;
            sbFront = new StringBuilder(150);
            sbRear = new StringBuilder(150);
            const char delim = '\t';

            // Tracks protein names:
            // keys are 0, 1, 2, etc.
            // Values are protein names
            SortedList<int, string> exportList;
            exportList = new SortedList<int, string>();
            CalculateScoreComponents();

            // sbFront and srRear contain the constant parts of the peptide result
            // each multiorf reference will be stuck in the middle during the integration
            // phase at the end of the function

            {
                var withBlock = this;
                sbFront.Append(withBlock.HitNum.ToString());
                sbFront.Append(delim);
                sbFront.Append(withBlock.StartScanNum.ToString("##0000"));
                sbFront.Append(delim);
                sbFront.Append(withBlock.ScanCount.ToString());
                sbFront.Append(delim);
                sbFront.Append(withBlock.ChargeState.ToString());
                sbFront.Append(delim);
                sbFront.Append(Math.Round(withBlock.MH, 5).ToString("#####0.00000"));
                sbFront.Append(delim);
                sbFront.Append(Math.Round(withBlock.XCorr, 4).ToString("##0.0000"));
                sbFront.Append(delim);
                sbFront.Append(Math.Round(withBlock.DelCn, 4).ToString("##0.0000"));
                sbFront.Append(delim);
                sbFront.Append(Math.Round(withBlock.Sp, 1).ToString("######0.0"));
                sbFront.Append(delim);
                if (withBlock.MultiProteinCount > 0)
                {
                    sbRear.Append(withBlock.MultiProteinCount.ToString("+0"));
                }
                else
                {
                    sbRear.Append("0");
                }

                sbRear.Append(delim);
                sbRear.Append(withBlock.Peptide);
                sbRear.Append(delim);
                sbRear.Append(Math.Round(withBlock.DelCn2, 4).ToString("##0.0000"));
                sbRear.Append(delim);
                sbRear.Append(withBlock.RankSp.ToString());
                sbRear.Append(delim);
                sbRear.Append(withBlock.RankXc.ToString());
                sbRear.Append(delim);
                sbRear.Append(Math.Round(withBlock.DelM, 5).ToString("##0.00000"));
                sbRear.Append(delim);
                sbRear.Append(Math.Round(withBlock.XcRatio, 3).ToString("0.000"));
                sbRear.Append(delim);

                // PassFilt and MScore are Legacy columns
                // sbRear.Append(.PassFilt.ToString)
                // sbRear.Append(delim)
                // If .MScore = 10 OrElse .MScore = 0 Then
                // sbRear.Append(.MScore.ToString("0"))
                // Else
                // sbRear.Append(Math.Round(.MScore, 2).ToString("##0.00"))
                // End If
                // sbRear.Append(delim)

                sbRear.Append(withBlock.ObsIons.ToString());
                sbRear.Append(delim);
                sbRear.Append(withBlock.PossIons.ToString());
                sbRear.Append(delim);
                sbRear.Append(withBlock.NumTrypticEnds.ToString());
                sbRear.Append(delim);
                sbRear.Append(withBlock.DelMPPM.ToString("##0.0000"));
            }

            exportList.Add(0, sbFront.ToString() + Reference + delim.ToString() + sbRear.ToString());
            if (ExpandMultiProteinEntries && m_MultiProteinEntries is object)
            {
                // Keys in m_MultiProteinEntries are MultiProteinID #, values are protein names
                foreach (KeyValuePair<int, string> objEntry in m_MultiProteinEntries)
                {
                    if ((Reference ?? "") != (objEntry.Value ?? ""))
                    {
                        exportList.Add(objEntry.Key, sbFront.ToString() + objEntry.Value.ToString() + delim.ToString() + sbRear.ToString());
                    }
                }
            }

            return exportList;
        }

        #region  Properties

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
            get => m_Reference;

            set
            {
                m_Reference = value;
                AddMultiProteinRef(value);
            }
        }

        public int MultiProteinCount { get; set; }
        public string Peptide { get; set; }
        public double DelCn2 { get; set; } = 0.0d;
        public int RankSp { get; set; }
        public int RankXc { get; set; }
        public double DelM { get; set; } = 0.0d;
        public double XcRatio { get; set; } = 1.0d;
        
        public int ObsIons { get; set; }
        public int PossIons { get; set; }
        public int NumTrypticEnds { get; set; }
        public double DelMPPM { get; set; } = 0d;

        #endregion

    }
}