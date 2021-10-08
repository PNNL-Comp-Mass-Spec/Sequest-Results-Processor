
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
        private double m_HeaderMass;
        private PeptideHitEntry m_CachedPreviousHit;
        private double m_CachedHighestXcorr;
        private SortedList<int, PeptideHitEntry> m_PeptideHits;
        private readonly ComputeDelMPPM m_ComputeDelMPPM;

        public ResultsFileEntry(int startScanNum, int endScanNum, int chargeState)
        {
            StartScanNumber = startScanNum;
            EndScanNumber = endScanNum;
            ChargeState = chargeState;
            m_ComputeDelMPPM = new ComputeDelMPPM();
        }

        public void AddHeaderMass(double HeaderMass)
        {
            m_HeaderMass = HeaderMass;
        }

        public void AddTotalIntensity(double HeaderIntensity)
        {
        }

        public void AddPeptideResults(PeptideHitEntry peptideResults)
        {
            if (m_PeptideHits is null)
            {
                m_PeptideHits = new SortedList<int, PeptideHitEntry>();
            }

            int newIndex;

            // See if peptideResults is already present in m_PeptideHits
            foreach (var objItem in m_PeptideHits)
            {
                if (peptideResults.StartScanNum == objItem.Value.StartScanNum && peptideResults.ChargeState == objItem.Value.ChargeState && (peptideResults.Peptide ?? "") == (objItem.Value.Peptide ?? ""))
                {
                    // Duplicate entry; update XCorr if higher
                    if (peptideResults.XCorr > objItem.Value.XCorr)
                    {
                        objItem.Value.XCorr = peptideResults.XCorr;
                        if (objItem.Key == 0)
                        {
                            m_CachedHighestXcorr = objItem.Value.XCorr;
                        }
                    }

                    return;
                }
            }

            newIndex = m_PeptideHits.Count + 1;
            if (newIndex > 1)
            {
                m_CachedPreviousHit = m_PeptideHits[newIndex - 1];
                m_CachedPreviousHit.DelCn2 = CalculateDelCN2(m_CachedPreviousHit.XCorr, peptideResults.XCorr);
            }
            else
            {
                m_CachedHighestXcorr = peptideResults.XCorr;
            }

            peptideResults.DelM = CalculateDelM(m_HeaderMass, peptideResults.MH);
            peptideResults.DelMPPM = CalculateDelMPPM(m_HeaderMass, peptideResults.MH);
            peptideResults.XcRatio = CalculateXCRatio(peptideResults.XCorr, m_CachedHighestXcorr);
            m_PeptideHits.Add(newIndex, peptideResults);
        }

        #region Calculate Differences

        private double CalculateDelCN2(double XCorrCurrent, double XCorrNextLowest)
        {
            double tmpDelCN2;
            tmpDelCN2 = (XCorrCurrent - XCorrNextLowest) / XCorrCurrent;
            return tmpDelCN2;
        }

        private double CalculateXCRatio(double XCorrCurrent, double XCorrFirstHit)
        {
            double tmpXCR;
            tmpXCR = XCorrCurrent / XCorrFirstHit;
            return tmpXCR;
        }

        private double CalculateDelM(double HeaderMass, double PeptideObsMass)
        {
            return PeptideObsMass - HeaderMass;
        }

        private double CalculateDelMPPM(double dblPrecursorMassMH, double dblPeptideTheoreticalMH)
        {
            return m_ComputeDelMPPM.ComputeDelMCorrected(dblPrecursorMassMH, dblPeptideTheoreticalMH);
        }

        #endregion

        #region  Properties

        public int StartScanNumber { get; private set; }
        public int EndScanNumber { get; private set; }
        public int ChargeState { get; private set; }

        public SortedList<int, PeptideHitEntry> PeptideHits => m_PeptideHits;

        #endregion
    }
}