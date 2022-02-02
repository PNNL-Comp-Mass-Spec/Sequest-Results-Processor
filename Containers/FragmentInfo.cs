
// This class stores information about observed peptide fragments and is used
// by class DiscriminantCalc
//
// -------------------------------------------------------------------------------
// Written by Ken Auberry for the Department of Energy (PNNL, Richland, WA) in 2005
// Converted to C# in 2021
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0

using System;
using System.Collections.Generic;

namespace SequestResultsProcessor.Containers
{
    public class FragmentInfo
    {
        private double mMaxIntensity;

        public struct Fragment
        {
            public double MZ;
            public double Intensity;
            public double NormIntensity;
        }

        public FragmentInfo()
        {
            FragmentList = new List<Fragment>();
        }

        public List<Fragment> FragmentList { get; }

        public virtual void Add(double mZ, double intensity)
        {
            var f = new Fragment
            {
                MZ = mZ,
                Intensity = intensity
            };

            FragmentList.Add(f);
            if (intensity > mMaxIntensity)
            {
                mMaxIntensity = intensity;
            }
        }

        public void Clear()
        {
            FragmentList.Clear();
            mMaxIntensity = 0d;
        }

        public double GetMass(int index)
        {
            return FragmentList[index].MZ;
        }

        // ReSharper disable once UnusedMember.Global
        public double GetIntensity(int index)
        {
            return FragmentList[index].Intensity;
        }

        public double GetNormalizedIntensity(int index)
        {
            if (Math.Abs(mMaxIntensity) < float.Epsilon)
                return 0d;
            return FragmentList[index].Intensity / mMaxIntensity;
        }

        public double GetNormalizedIntensity(Fragment f)
        {
            if (Math.Abs(mMaxIntensity) < float.Epsilon)
                return 0d;
            return f.Intensity / mMaxIntensity;
        }

        protected void NormalizeIntensities()
        {
            mMaxIntensity = 0d;
            foreach (var f in FragmentList)
            {
                if (f.Intensity > mMaxIntensity)
                {
                    mMaxIntensity = f.Intensity;
                }
            }

            if (Math.Abs(mMaxIntensity) < float.Epsilon)
                return;

            foreach (var item in FragmentList)
            {
                var fragment = item;
                fragment.NormIntensity = fragment.Intensity / mMaxIntensity;
            }
        }
    }
}