
// This class stores information about observed peptide fragments and is used
// by class DiscriminantCalc
//
// -------------------------------------------------------------------------------
// Written by Ken Auberry for the Department of Energy (PNNL, Richland, WA)
// Coding for this class began on 10 October 2005
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0

using System;
using System.Collections.Generic;

namespace SequestResultsProcessor.Containers
{

    // ReSharper disable once InconsistentNaming
    public class FragmentInfo
    {
        protected List<Fragment> m_FragmentList;
        protected double m_MaxIntensity;

        public struct Fragment
        {
            public double MZ;
            public double Intensity;
            public double NormIntensity;
        }

        public FragmentInfo()
        {
            m_FragmentList = new List<Fragment>();
        }

        public List<Fragment> FragmentList => m_FragmentList;

        public int maxIndex => m_FragmentList.Count - 1;

        public virtual void Add(double mZ, double intensity)
        {
            var f = new Fragment();
            f.MZ = mZ;
            f.Intensity = intensity;
            m_FragmentList.Add(f);
            if (intensity > m_MaxIntensity)
            {
                m_MaxIntensity = intensity;
            }
        }

        public void Clear()
        {
            m_FragmentList.Clear();
            m_MaxIntensity = 0d;
        }

        public double GetMass(int index)
        {
            return m_FragmentList[index].MZ;
        }

        public double GetIntensity(int index)
        {
            return m_FragmentList[index].Intensity;
        }

        public double GetNormalizedIntensity(int index)
        {
            if (Math.Abs(m_MaxIntensity) < float.Epsilon)
                return 0d;
            return m_FragmentList[index].Intensity / m_MaxIntensity;
        }

        public double GetNormalizedIntensity(Fragment f)
        {
            if (Math.Abs(m_MaxIntensity) < float.Epsilon)
                return 0d;
            return f.Intensity / m_MaxIntensity;
        }

        protected virtual void NormalizeIntensities()
        {
            m_MaxIntensity = 0d;
            foreach (var f in m_FragmentList)
            {
                if (f.Intensity > m_MaxIntensity)
                {
                    m_MaxIntensity = f.Intensity;
                }
            }

            if (Math.Abs(m_MaxIntensity) < float.Epsilon)
                return;
            foreach (var f in m_FragmentList)
                f.NormIntensity = f.Intensity / m_MaxIntensity;
        }
    }
}