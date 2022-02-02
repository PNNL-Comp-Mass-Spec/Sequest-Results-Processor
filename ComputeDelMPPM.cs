
namespace SequestResultsProcessor
{
    public class ComputeDelMPPM
    {
        public const double MASS_C13 = 1.00335483d;
        public const double MASS_PROTON = 1.00727649d;               // Note that this is the mass of hydrogen minus the mass of one electron

        public double ComputeDelMCorrected(double precursorMassMH, double peptideTheoreticalMH)
        {
            var delM = precursorMassMH - peptideTheoreticalMH;
            var precursorMonoMass = precursorMassMH - MASS_PROTON;
            var peptideMonoisotopicMass = peptideTheoreticalMH - MASS_PROTON;
            return ComputeDelMCorrected(delM, precursorMonoMass, true, peptideMonoisotopicMass);
        }

        public double ComputeDelMCorrected(double delM, double precursorMonoMass, bool adjustPrecursorMassForC13, double peptideMonoisotopicMass)
        {
            var correctionCount = 0;

            // Examine delM to determine which isotope was chosen
            if (delM >= -0.5d)
            {
                // This is the typical case
                while (delM > 0.5d)
                {
                    delM -= MASS_C13;
                    correctionCount++;
                }
            }
            else
            {
                // This happens less often; but we'll still account for it
                // In this case, correctionCount will be negative
                while (delM < -0.5d)
                {
                    delM += MASS_C13;
                    correctionCount--;
                }
            }

            if (correctionCount != 0)
            {
                if (adjustPrecursorMassForC13)
                {
                    // Adjust the precursor mono mass based on correctionCount
                    precursorMonoMass -= correctionCount * MASS_C13;
                }

                // Compute a new DelM value
                delM = precursorMonoMass - peptideMonoisotopicMass;
            }

            return MassToPPM(delM, peptideMonoisotopicMass);
        }

        public double MassToPPM(double massToConvert, double currentMZ)
        {
            // Converts massToConvert to ppm, based on the value of currentMZ

            return massToConvert * 1000000.0d / currentMZ;
        }
    }
}