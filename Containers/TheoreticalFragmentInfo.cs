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
        private static Dictionary<string, ResidueInfo> s_Intensities = new();
        private List<Fragment> mTheoreticalYIons;
        private List<Fragment> mTheoreticalBIons;

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
            // ReSharper disable UnusedMember.Global
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
            // ReSharper restore UnusedMember.Global
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

        public List<Fragment> YIons => mTheoreticalYIons;

        public List<Fragment> BIons => mTheoreticalBIons;

        public TheoreticalFragmentInfo(string peptideSequence, int chargeState)
        {
            if (s_Intensities is null)
            {
                SetupResidueHashes();
            }

            CalculateTheoreticalIonHashes(peptideSequence, chargeState);
        }

        private void SetupResidueHashes()
        {
            var peptideResidues = KNOWN_RESIDUES.ToCharArray();
            var peptideLeft = new[] { -0.2d, -0.75d, 0.45d, -0.05d, 0.4d, -0.8d, 0.35d, -0.4d, -0.25d, -0.15d, -0.2d, -0.5d, -1.15d, -0.35d, -0.75d, -0.6d, -0.65d, 0d, 0.25d, -0.4d };
            var peptideRight = new[] { 0.35d, -0.2d, -0.45d, -0.15d, 0.45d, 0.5d, 0.25d, 0.35d, 0.3d, 0.05d, 0.1d, 0.1d, 1.15d, -0.05d, -0.35d, 0.5d, 0.45d, 0.2d, 0.45d, 0.4d };
            var peptideMass = new[] { 71.037d, 103.009d, 115.027d, 129.043d, 147.068d, 57.022d, 137.059d, 113.084d, 128.095d, 113.084d, 131.04d, 114.043d, 97.053d, 128.059d, 156.101d, 87.032d, 101.048d, 99.068d, 186.08d, 163.063d };

            var maxIndex = peptideResidues.Length - 1;
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
            var cleanSeq = GetCleanSequence(dirtySeq);
            var peptideLength = cleanSeq.Length;
            mTheoreticalYIons = new List<Fragment>();
            mTheoreticalBIons = new List<Fragment>();
            int counter;

            var tmpPepMass = GetMass(cleanSeq, chargeState);
            var tmpBMass = 1.01d;

            for (var i = 0; i < peptideLength - 1; i++)
            {
                var intensityFrag = cleanSeq.Substring(i, 2);
                var bFragment = intensityFrag.Substring(0, 1);
                var yFragment = intensityFrag.Substring(1, 1);

                tmpBMass += GetMass(bFragment, chargeState);
                var tmpYMass = tmpPepMass - tmpBMass + 20.02d;
                var tmpInt = GetIntensity(bFragment, yFragment);
                mTheoreticalBIons.Add(new Fragment(Math.Round(tmpBMass, 4), tmpInt, FragmentTypes.b, i + 1));
                mTheoreticalYIons.Add(new Fragment(Math.Round(tmpYMass, 4), tmpInt, FragmentTypes.y, peptideLength - (i + 1)));
            }

            var sorter = new TheoreticalIonComparer();
            mTheoreticalBIons.Sort(sorter);
            mTheoreticalYIons.Sort(sorter);
            if (mTheoreticalBIons.Count > 0)
                mTheoreticalBIons.RemoveAt(mTheoreticalBIons.Count - 1);
            if (mTheoreticalYIons.Count > 0)
                mTheoreticalYIons.RemoveAt(mTheoreticalYIons.Count - 1);
        }

        private string GetCleanSequence(string rawPeptideSeq)
        {
            var r = new Regex(@"^.*\.(?<cleanseq>\S+)\..*$");
            var m = r.Match(rawPeptideSeq);
            var cleanSeq = m.Groups["cleanseq"].Value;
            cleanSeq = Regex.Replace(cleanSeq, "[^" + KNOWN_RESIDUES + "]", "");
            return cleanSeq.ToUpper();
        }

        private double GetMass(string sequence, int cs)
        {
            var seq = sequence.ToCharArray();
            var tmpMass = 0.0;
            double tmpMZ;
            foreach (var resChar in seq)
            {
                if (s_Intensities.TryGetValue(resChar.ToString(), out var residue))
                {
                    tmpMass += residue.Mass;
                }
            }

            return tmpMass / cs;
        }

        private double GetIntensity(string LeftResidue, string RightResidue)
        {
            if (s_Intensities.TryGetValue(LeftResidue, out var resLeft) && s_Intensities.TryGetValue(RightResidue, out var resRight))
            {
                return Math.Round(resLeft.LeftIntensity + resRight.RightIntensity, 2);
            }

            return 0d;
        }

        private class TheoreticalIonComparer : IComparer<Fragment>
        {
            public int Compare(Fragment x, Fragment y)
            {
                if (x.Mass > y.Mass)
                {
                    return 1;
                }

                if (x.Mass < y.Mass)
                {
                    return -1;
                }

                return 0;
            }
        }
    }
}