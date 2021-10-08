using System;
using System.Collections.Generic;
using System.Text.RegularExpressions;

namespace SequestResultsProcessor.Containers
{

    // This class stores information about theoretical peptide fragments and is used by
    // clsDiscriminantCalc
    //
    // -------------------------------------------------------------------------------
    // Written by Ken Auberry for the Department of Energy (PNNL, Richland, WA)
    // Coding for this class began on 10 October 2005
    // -------------------------------------------------------------------------------
    //
    // Licensed under the Apache License, Version 2.0; you may not use this file except
    // in compliance with the License.  You may obtain a copy of the License at
    // http://www.apache.org/licenses/LICENSE-2.0

    internal class TheoreticalFragmentInfo
    {
        private static Dictionary<string, ResidueInfo> s_Intensities = new Dictionary<string, ResidueInfo>();
        private List<Fragment> m_TheoYIons;
        private List<Fragment> m_TheoBIons;
        private string m_Sequence;

        public struct Fragment
        {
            public double Mass;
            public double Intensity;
            public FragmentTypes FragmentType;
            public int FragmentPosition;

            public Fragment(double Mass, double Intensity, FragmentTypes FragmentType, int FragmentPosition)
            {
                this.Mass = Mass;
                this.Intensity = Intensity;
                this.FragmentType = FragmentType;
                this.FragmentPosition = FragmentPosition;
            }
        }

        public struct ResidueInfo
        {
            public Residues Name;
            public double Mass;
            public double LeftIntensity;
            public double RightIntensity;
        }

        public enum Residues
        {
            A,
            C,
            D,
            E,
            F,
            G,
            H,
            I,
            K,
            L,
            M,
            N,
            P,
            Q,
            R,
            S,
            T,
            V,
            W,
            Y
        }

        public enum FragmentTypes
        {
            // a
            b,
            // c
            // x
            y
            // z
        }

        private const string KNOWN_RESIDUES = "ACDEFGHIKLMNPQRSTVWY";

        public List<Fragment> YIons => m_TheoYIons;

        public List<Fragment> BIons => m_TheoBIons;

        public TheoreticalFragmentInfo(string peptideSequence, int chargeState)
        {
            if (s_Intensities is null)
            {
                SetupResidueHashes();
            }

            m_Sequence = peptideSequence;
            CalculateTheoreticalIonHashes(peptideSequence, chargeState);
        }

        private void SetupResidueHashes()
        {
            char[] peptideResidues;
            double[] peptideLeft;
            double[] peptideRight;
            double[] peptideMass;

            peptideResidues = KNOWN_RESIDUES.ToCharArray();
            peptideLeft = new double[] { -0.2d, -0.75d, 0.45d, -0.05d, 0.4d, -0.8d, 0.35d, -0.4d, -0.25d, -0.15d, -0.2d, -0.5d, -1.15d, -0.35d, -0.75d, -0.6d, -0.65d, 0d, 0.25d, -0.4d };
            peptideRight = new double[] { 0.35d, -0.2d, -0.45d, -0.15d, 0.45d, 0.5d, 0.25d, 0.35d, 0.3d, 0.05d, 0.1d, 0.1d, 1.15d, -0.05d, -0.35d, 0.5d, 0.45d, 0.2d, 0.45d, 0.4d };
            peptideMass = new double[] { 71.037d, 103.009d, 115.027d, 129.043d, 147.068d, 57.022d, 137.059d, 113.084d, 128.095d, 113.084d, 131.04d, 114.043d, 97.053d, 128.059d, 156.101d, 87.032d, 101.048d, 99.068d, 186.08d, 163.063d };
            int maxIndex = peptideResidues.Length - 1;
            int counter;
            s_Intensities = new Dictionary<string, ResidueInfo>();
            var loopTo = maxIndex;
            for (counter = 0; counter <= loopTo; counter++)
            {
                var residueSymbol = peptideResidues[counter].ToString();
                var r = new ResidueInfo();
                if (!Enum.TryParse(residueSymbol, out Residues symbol))
                {
                    throw new Exception("Unrecognized residue: " + residueSymbol);
                }

                r.Name = symbol;
                r.Mass = peptideMass[counter];
                r.LeftIntensity = peptideLeft[counter];
                r.RightIntensity = peptideRight[counter];
                s_Intensities.Add(residueSymbol, r);
            }
        }

        private void CalculateTheoreticalIonHashes(string dirtySeq, int chargeState)
        {
            string cleanSeq = GetCleanSequence(dirtySeq);
            int peptideLength = cleanSeq.Length;
            string IntensityFrag;
            string BFragment;
            string YFragment;
            double tmpYMass;
            double tmpBMass;
            double tmpPepMass;
            double tmpInt;
            m_TheoYIons = new List<Fragment>();
            m_TheoBIons = new List<Fragment>();
            int counter;
            tmpPepMass = GetMass(cleanSeq, chargeState);
            tmpBMass = 1.01d;
            var loopTo = peptideLength - 1;
            for (counter = 1; counter <= loopTo; counter++)
            {
                IntensityFrag = Strings.Mid(cleanSeq, counter, 2);
                BFragment = Strings.Left(IntensityFrag, 1);
                YFragment = Strings.Right(IntensityFrag, 1);
                tmpBMass += GetMass(BFragment, chargeState);
                tmpYMass = tmpPepMass - tmpBMass + 20.02d;
                tmpInt = GetIntensity(BFragment, YFragment);
                m_TheoBIons.Add(new Fragment(Math.Round(tmpBMass, 4), tmpInt, FragmentTypes.b, counter));
                m_TheoYIons.Add(new Fragment(Math.Round(tmpYMass, 4), tmpInt, FragmentTypes.y, peptideLength - counter));
            }

            var sorter = new TheoIonComparer();
            m_TheoBIons.Sort(sorter);
            m_TheoYIons.Sort(sorter);
            if (m_TheoBIons.Count > 0)
                m_TheoBIons.RemoveAt(m_TheoBIons.Count - 1);
            if (m_TheoYIons.Count > 0)
                m_TheoYIons.RemoveAt(m_TheoYIons.Count - 1);
        }

        private string GetCleanSequence(string rawPeptideSeq)
        {
            var r = new Regex(@"^*.\.(?<cleanseq>\S+)\..*$");
            var m = r.Match(rawPeptideSeq);
            string cleanSeq = m.Groups["cleanseq"].Value.ToString();
            cleanSeq = Regex.Replace(cleanSeq, "[^" + KNOWN_RESIDUES + "]", "");
            return cleanSeq.ToUpper();
        }

        private double GetMass(string sequence, int cs)
        {
            var seq = sequence.ToCharArray();
            ResidueInfo residue;
            var tmpMass = default(double);
            double tmpMZ;
            foreach (var resChar in seq)
            {
                if (s_Intensities.TryGetValue(resChar.ToString(), out residue))
                {
                    tmpMass += residue.Mass;
                }
            }

            tmpMZ = tmpMass / cs;
            return tmpMZ;
        }

        private double GetIntensity(string LeftResidue, string RightResidue)
        {
            ResidueInfo resLeft;
            ResidueInfo resRight;
            if (s_Intensities.TryGetValue(LeftResidue, out resLeft))
            {
                if (s_Intensities.TryGetValue(RightResidue, out resRight))
                {
                    double tmpIntens = Math.Round(resLeft.LeftIntensity + resRight.RightIntensity, 2);
                    return tmpIntens;
                }
            }

            return 0d;
        }

        private class TheoIonComparer : IComparer<Fragment>
        {
            public int Compare(Fragment x, Fragment y)
            {
                if (x.Mass > y.Mass)
                {
                    return 1;
                }
                else if (x.Mass < y.Mass)
                {
                    return -1;
                }
                else
                {
                    return 0;
                }
            }
        }
    }
}