
// This set of classes stages and stores information to be included in the generated
// synopsis and first hits files
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

namespace SequestResultsProcessor.Containers
{
    /// <summary>
    /// Tracks the peptides listed in a single .Out file
    /// </summary>
    /// <remarks></remarks>
    internal class ResultsFileEntry
    {
        private double mHeaderMass;
        private PeptideHitEntry mCachedPreviousHit;
        private double mCachedHighestXcorr;
        private SortedList<int, PeptideHitEntry> mPeptideHits;
        private readonly ComputeDelMPPM mComputeDelMPPM;

        public ResultsFileEntry(int startScanNum, int endScanNum, int chargeState)
        {
            StartScanNumber = startScanNum;
            EndScanNumber = endScanNum;
            ChargeState = chargeState;
            mComputeDelMPPM = new ComputeDelMPPM();
        }

        public void AddHeaderMass(double headerMass)
        {
            mHeaderMass = headerMass;
        }

        public void AddPeptideResults(PeptideHitEntry peptideResults)
        {
            mPeptideHits ??= new SortedList<int, PeptideHitEntry>();

            // See if peptideResults is already present in mPeptideHits
            foreach (var objItem in mPeptideHits)
            {
                if (peptideResults.StartScanNum == objItem.Value.StartScanNum && peptideResults.ChargeState == objItem.Value.ChargeState && (peptideResults.Peptide ?? "") == (objItem.Value.Peptide ?? ""))
                {
                    // Duplicate entry; update XCorr if higher
                    if (peptideResults.XCorr > objItem.Value.XCorr)
                    {
                        objItem.Value.XCorr = peptideResults.XCorr;
                        if (objItem.Key == 0)
                        {
                            mCachedHighestXcorr = objItem.Value.XCorr;
                        }
                    }

                    return;
                }
            }

            var newIndex = mPeptideHits.Count + 1;
            if (newIndex > 1)
            {
                mCachedPreviousHit = mPeptideHits[newIndex - 1];
                mCachedPreviousHit.DelCn2 = CalculateDelCN2(mCachedPreviousHit.XCorr, peptideResults.XCorr);
            }
            else
            {
                mCachedHighestXcorr = peptideResults.XCorr;
            }

            peptideResults.DelM = CalculateDelM(mHeaderMass, peptideResults.MH);
            peptideResults.DelMPPM = CalculateDelMPPM(mHeaderMass, peptideResults.MH);
            peptideResults.XcRatio = CalculateXCRatio(peptideResults.XCorr, mCachedHighestXcorr);
            mPeptideHits.Add(newIndex, peptideResults);
        }

        #region Calculate Differences
        private double CalculateDelCN2(double xcorrCurrent, double xcorrNextLowest)
        {
            return (xcorrCurrent - xcorrNextLowest) / xcorrCurrent;
        }

        private double CalculateXCRatio(double xcorrCurrent, double xcorrFirstHit)
        {
            return xcorrCurrent / xcorrFirstHit;
        }

        private double CalculateDelM(double headerMass, double peptideObsMass)
        {
            return peptideObsMass - headerMass;
        }

        private double CalculateDelMPPM(double precursorMassMH, double peptideTheoreticalMH)
        {
            return mComputeDelMPPM.ComputeDelMCorrected(precursorMassMH, peptideTheoreticalMH);
        }

        #endregion

        #region  Properties

        public int StartScanNumber { get; }
        public int EndScanNumber { get; }
        public int ChargeState { get; }

        public SortedList<int, PeptideHitEntry> PeptideHits => mPeptideHits;

        #endregion
    }
}