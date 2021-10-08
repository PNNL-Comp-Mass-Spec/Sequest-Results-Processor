
// This class calculates intensities for observed neutral loss species in
// phosphoproteomic studies
//
// -------------------------------------------------------------------------------
// Written by Ken Auberry for the Department of Energy (PNNL, Richland, WA)
// Coding for this class began on 19 November 2004
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0

using System;

namespace SequestResultsProcessor
{
    internal class CalcNeutralLosses
    {
        public NeutralLossList CalculateNeutralLosses(DiscriminantCalc.PeptideIntensities scan, double massTolerance)
        {
            var nll = new NeutralLossList();
            var tmpNeuLoss1 = default(double);
            var tmpNeuLoss2 = default(double);
            var tmpNeuLoss3 = default(double);
            double tmpMass;
            double tmpAbun;
            var parentMZ = scan.ParentMZ;
            var NLPeak1MZ = parentMZ - 32.7d;
            var NLPeak2MZ = parentMZ - 49.0d;
            var NLPeak3MZ = parentMZ - 98.0d;
            foreach (var f in scan.FragmentList)
            {
                tmpMass = f.MZ;
                tmpAbun = scan.GetNormalizedIntensity(f);
                if (Math.Abs(tmpMass - NLPeak1MZ) <= massTolerance & tmpAbun > tmpNeuLoss1)
                {
                    nll.NL1Intensity = Math.Round(tmpAbun, 2);
                }
                else if (Math.Abs(tmpMass - NLPeak2MZ) <= massTolerance & tmpAbun > tmpNeuLoss2)
                {
                    nll.NL2Intensity = Math.Round(tmpAbun, 2);
                }
                else if (Math.Abs(tmpMass - NLPeak3MZ) <= massTolerance & tmpAbun > tmpNeuLoss3)
                {
                    nll.NL3Intensity = Math.Round(tmpAbun, 2);
                }
                else if (tmpMass > parentMZ + massTolerance)
                {
                    break;
                }
            }

            return nll;
        }

        public NeutralLossList GetFillerNeutralLosses()
        {
            var nll = new NeutralLossList();
            return nll;
        }
    }

    internal class NeutralLossList
    {
        public double NL1Intensity { get; set; } = 0.0d;
        public double NL2Intensity { get; set; } = 0.0d;
        public double NL3Intensity { get; set; } = 0.0d;
    }
}