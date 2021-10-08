
namespace SequestResultsProcessor
{
    public class ComputeDelMPPM
    {
        #region Constants
        public const double MASS_C13 = 1.00335483d;
        public const double MASS_PROTON = 1.00727649d;               // Note that this is the mass of hydrogen minus the mass of one electron
        #endregion

        public double ComputeDelMCorrected(double dblPrecursorMassMH, double dblPeptideTheoreticalMH)
        {
            double dblDelM;
            double dblPrecursorMonoMass;
            double dblPeptideMonoisotopicMass;
            dblDelM = dblPrecursorMassMH - dblPeptideTheoreticalMH;
            dblPrecursorMonoMass = dblPrecursorMassMH - MASS_PROTON;
            dblPeptideMonoisotopicMass = dblPeptideTheoreticalMH - MASS_PROTON;
            return ComputeDelMCorrected(dblDelM, dblPrecursorMonoMass, true, dblPeptideMonoisotopicMass);
        }

        public double ComputeDelMCorrected(double dblDelM, double dblPrecursorMonoMass, bool blnAdjustPrecursorMassForC13, double dblPeptideMonoisotopicMass)
        {
            int intCorrectionCount = 0;

            // Examine dblDelM to determine which isotope was chosen
            if (dblDelM >= -0.5d)
            {
                // This is the typical case
                while (dblDelM > 0.5d)
                {
                    dblDelM -= MASS_C13;
                    intCorrectionCount += 1;
                }
            }
            else
            {
                // This happens less often; but we'll still account for it
                // In this case, intCorrectionCount will be negative
                while (dblDelM < -0.5d)
                {
                    dblDelM += MASS_C13;
                    intCorrectionCount -= 1;
                }
            }

            if (intCorrectionCount != 0)
            {
                if (blnAdjustPrecursorMassForC13)
                {
                    // Adjust the precursor mono mass based on intCorrectionCount
                    dblPrecursorMonoMass -= intCorrectionCount * MASS_C13;
                }

                // Compute a new DelM value
                dblDelM = dblPrecursorMonoMass - dblPeptideMonoisotopicMass;
            }

            return MassToPPM(dblDelM, dblPeptideMonoisotopicMass);
        }

        public double MassToPPM(double dblMassToConvert, double dblCurrentMZ)
        {
            // Converts dblMassToConvert to ppm, based on the value of dblCurrentMZ

            return dblMassToConvert * 1000000.0d / dblCurrentMZ;
        }
    }
}