using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Text.RegularExpressions;
using SequestResultsProcessor.Containers;

namespace SequestResultsProcessor
{
    // This class calculates the M-Score component of the discriminant score that was
    // originally designed to support the AMT Tag approach to proteomics studies
    //
    // -------------------------------------------------------------------------------
    // Originally written in Perl by Eric F. Strittmatter
    // Ported to VB.NET by Ken Auberry for the Department of Energy (PNNL, Richland, WA)
    // Coding for this class began on 19 November 2004
    // -------------------------------------------------------------------------------
    //
    // Licensed under the Apache License, Version 2.0; you may not use this file except
    // in compliance with the License.  You may obtain a copy of the License at
    // http://www.apache.org/licenses/LICENSE-2.0

    internal class DiscriminantCalc
    {
        private DTAFileInformation _m_dtaFileInfo;

        private DTAFileInformation m_dtaFileInfo
        {
            [MethodImpl(MethodImplOptions.Synchronized)]
            get => _m_dtaFileInfo;

            [MethodImpl(MethodImplOptions.Synchronized)]
            set
            {
                if (_m_dtaFileInfo != null)
                {
                    _m_dtaFileInfo.OffsetProgress -= OffsetLoadingProgressHandler;
                }

                _m_dtaFileInfo = value;
                if (_m_dtaFileInfo != null)
                {
                    _m_dtaFileInfo.OffsetProgress += OffsetLoadingProgressHandler;
                }
            }
        }

        private PeptideIntensities m_CachedScanInfo;
        private string m_CachedScanKey;
        private int m_CachedCS;
        private readonly OutputNLIFile m_NLIDumper;
        private double m_MassTol;
        private readonly string m_dtaFilepath;
        private readonly string m_Version;
        private readonly bool m_NoDTAs;

        public event ProgressUpdateEventHandler ProgressUpdate;

        public delegate void ProgressUpdateEventHandler(string TaskDescription, double fractionDone);

        private void OnProgressUpdate(string taskDescription, double fractionDone)
        {
            ProgressUpdate?.Invoke(taskDescription, fractionDone);
        }

        private void OffsetLoadingProgressHandler(double fractionDone)
        {
            OnProgressUpdate("Loading .dta File Locations (" + (fractionDone * 100d).ToString("0.0") + "% Completed)", fractionDone);
        }

        public DiscriminantCalc(string dtaFilePath)
        {
            var fi = new FileInfo(dtaFilePath);
            var rootFileName = fi.Name;
            rootFileName = Regex.Replace(rootFileName, "_dta.txt", "");
            m_dtaFilepath = dtaFilePath;
            if (fi.Exists)
            {
                m_dtaFileInfo = new DTAFileInformation(dtaFilePath);
                m_NoDTAs = false;
                m_NLIDumper = new OutputNLIFile(rootFileName, fi.DirectoryName);
            }
            else
            {
                m_NoDTAs = true;
            }

            m_Version = FileVersionInfo.GetVersionInfo(Assembly.GetExecutingAssembly().Location).FileVersion;
        }

        /// <summary>
        /// Version of the module
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public string Version()
        {
            return m_Version;
        }

        // --------------------
        // Calculates MScore for given peptide record
        // using fragment intensities from
        // concatenated dta file
        [Obsolete("Unused")]
        private double GetPeptideMScore(string PeptideSeq, int StartScanNumber, int EndScanNumber, int ChargeState)
        {
            var ScanKey = StartScanNumber + "." + EndScanNumber;
            if (m_NoDTAs)
                return 10d;

            m_dtaFileInfo ??= new DTAFileInformation(m_dtaFilepath);

            m_MassTol = 0.7d;
            if (ChargeState != m_CachedCS || ScanKey != (m_CachedScanKey ?? string.Empty))
            {
                m_dtaFileInfo.GetDTAFileIntensities(StartScanNumber, EndScanNumber, ChargeState);
                if (m_dtaFileInfo.DTAScanInfo.FragmentList.Count == 0)
                {
                    return 10d;
                }

                if (ScanKey != (m_CachedScanKey ?? string.Empty))
                {
                    m_dtaFileInfo.DTAScanInfo.GetNeutralLosses(m_MassTol);
                    m_NLIDumper.MakeNLIEntry(StartScanNumber, m_dtaFileInfo.DTAScanInfo.NeutralLosses);
                }

                m_CachedScanInfo = m_dtaFileInfo.DTAScanInfo;
                m_CachedScanKey = ScanKey;
                m_CachedCS = ChargeState;
            }

            return CalculateMScore(PeptideSeq, ChargeState, m_CachedScanInfo);
        }

        private double CalculateMScore(string peptideSequence, int peptideChargeState, FragmentInfo scanIntensities)
        {
            var match = 0.0;

            var mTol = m_MassTol;
            if (scanIntensities.FragmentList.Count == 0)
            {
                return 10d;
            }

            if (Regex.IsMatch(peptideSequence.ToUpper(), "[BZXOJU*@#!$%^&]"))
            {
                return 10d;
            }

            if (peptideSequence.Length <= 5)
            {
                return 10d;
            }
            var noModsInSequence = true;

            // Get values for +1 charge state
            var theoreticalFragments = new TheoreticalFragmentInfo(peptideSequence, 1);

            // Check and score +1 possibilities for B-Ions
            match += HashScanner(scanIntensities, theoreticalFragments.BIons, mTol, 1);

            // Check and score +1 possibilities for Y-Ions
            match += HashScanner(scanIntensities, theoreticalFragments.YIons, mTol, 1);
            if (peptideChargeState > 2)
            {
                // get values for +2 charge state
                theoreticalFragments = new TheoreticalFragmentInfo(peptideSequence, 2);

                // Check and score +2 possibilities for B-Ions
                match += HashScanner(scanIntensities, theoreticalFragments.BIons, mTol, 2);

                // Check and score +2 possibilities for Y-Ions
                match += HashScanner(scanIntensities, theoreticalFragments.YIons, mTol, 2);
            }

            if (noModsInSequence)
            {
                return Math.Round(10d + match, 2);
            }

            return 10d;
        }

        private double HashScanner(FragmentInfo peptideRecord, List<TheoreticalFragmentInfo.Fragment> theoreticalFragments, double massTol, int CSToCheck)
        {
            var maxObsRecord = peptideRecord.FragmentList.Count;
            if (maxObsRecord == 0)
            {
                return 0.0d;
            }

            var maxObsIndex = maxObsRecord - 1;

            var obsCount = 0;
            var match = 0d;
            var tmpObsMass = peptideRecord.GetMass(obsCount);
            foreach (var fragment in theoreticalFragments)
            {
                var theoreticalMass = fragment.Mass;
                if (tmpObsMass > theoreticalMass + massTol && obsCount <= maxObsIndex)
                {
                }
                else
                {
                    while (tmpObsMass < theoreticalMass - massTol && obsCount < maxObsIndex)
                    {
                        obsCount++;
                        tmpObsMass = peptideRecord.GetMass(obsCount);
                    }

                    if (tmpObsMass < theoreticalMass + massTol && !(theoreticalMass > tmpObsMass + massTol) && obsCount <= maxObsIndex)

                    {
                        if (CSToCheck >= 2)
                        {
                            match += peptideRecord.GetNormalizedIntensity(obsCount);
                        }
                        else
                        {
                            match += peptideRecord.GetNormalizedIntensity(obsCount) * fragment.Intensity;
                        }
                    }
                }
            }

            return match;
        }

        private class DTAFileInformation
        {
            private dtaFileOffsets _m_Offsets;

            private dtaFileOffsets m_Offsets
            {
                [MethodImpl(MethodImplOptions.Synchronized)]
                get => _m_Offsets;

                [MethodImpl(MethodImplOptions.Synchronized)]
                set
                {
                    if (_m_Offsets != null)
                    {
                        _m_Offsets.dtaScanProgress -= OffsetLoadingProgressHandler;
                    }

                    _m_Offsets = value;
                    if (_m_Offsets != null)
                    {
                        _m_Offsets.dtaScanProgress += OffsetLoadingProgressHandler;
                    }
                }
            }

            private PeptideIntensities m_DTAFileIntensities;
            private FileStream m_dtaStream;
            private readonly string m_dtaFilePath;

            public event OffsetProgressEventHandler OffsetProgress;

            public delegate void OffsetProgressEventHandler(double fractionDone);

            public DTAFileInformation(string dtaFilePath)
            {
                // Scan the dta file for dta file boundaries
                m_Offsets = new dtaFileOffsets();
                m_dtaFilePath = dtaFilePath;
            }

            // ReSharper disable once UnusedMember.Local
            public void Configure()
            {
                m_Offsets.LoadOffsetsFromDTAFile(m_dtaFilePath);
                m_dtaStream = new FileStream(m_dtaFilePath, FileMode.Open);
                m_DTAFileIntensities = new PeptideIntensities(m_dtaStream);
            }

            // ReSharper disable once UnusedMember.Local
            public void Close()
            {
                m_dtaStream.Close();
            }

            public PeptideIntensities DTAScanInfo => m_DTAFileIntensities;

            private void OffsetLoadingProgressHandler(double fractionDone)
            {
                OffsetProgress?.Invoke(fractionDone);
            }

            public void GetDTAFileIntensities(int StartScanNumber, int EndScanNumber, int ChargeState)
            {
                var fileOffset = m_Offsets.get_GetOffset(StartScanNumber, EndScanNumber, ChargeState);
                if (fileOffset > 0L)
                {
                    {
                        var withBlock = m_DTAFileIntensities;
                        withBlock.GetIntensitiesFromDTA(fileOffset);
                    }
                }
                else if (m_DTAFileIntensities != null)
                {
                    m_DTAFileIntensities.Clear();
                }
                else
                {
                    m_DTAFileIntensities = new PeptideIntensities(m_dtaStream);
                }
            }

            private class dtaFileOffsets
            {
                private readonly Dictionary<string, long> mOffsets = new();

                public long get_GetOffset(int StartScanNumber, int EndScanNumber, int ChargeState)
                {
                    var keyName = StartScanNumber + "." + EndScanNumber + "." + ChargeState;
                    if (mOffsets.TryGetValue(keyName, out var offset))
                    {
                        return offset;
                    }

                    throw new KeyNotFoundException("Offsets dictionary does not have key " + keyName);
                }

                private void AddOffset(int StartScanNumber, int EndScanNumber, int ChargeState, long StartOffset, string ChargeExtra = "")
                {
                    var keyName = StartScanNumber + "." + EndScanNumber + "." + ChargeState;
                    if (ChargeExtra.Length > 0)
                    {
                        keyName += "_" + ChargeExtra;
                    }

                    mOffsets.Add(keyName, StartOffset);
                }

                public void LoadOffsetsFromDTAFile(string dtaFilePath)
                {
                    var fi = new FileInfo(dtaFilePath);
                    long currentPosition = 0;
                    var dtaCount = 0;
                    var r = new Regex(@"^\s*[=]{5,}\s+\""(?<rootname>.+)\.(?<StartScan>\d+)\.(?<EndScan>\d+)\.(?<ChargeBlock>(?<ChargeState>\d+)[^0-9]?(?<ChargeExtra>\S*))\.(?<FileType>.+)\""\s+[=]{5,}\s*$", RegexOptions.CultureInvariant | RegexOptions.Compiled);

                    var lineMatch = new Regex("^===*");
                    if (!fi.Exists)
                    {
                        return;
                    }

                    var lineEndCharCount = SequestFileExtractor.LineEndCharacterCount(fi);
                    var fileLength = fi.Length;
                    TextReader tr = fi.OpenText();

                    var dataLine = tr.ReadLine();
                    while (dataLine != null)
                    {
                        currentPosition += dataLine.Length + lineEndCharCount;
                        if (lineMatch.IsMatch(dataLine))
                        {
                            var dtaStartPos = currentPosition - dataLine.Length - lineEndCharCount;
                            var m = r.Match(dataLine);
                            string chargeExtra;
                            if (m.Groups["ChargeExtra"].Length > 0)
                            {
                                chargeExtra = m.Groups["ChargeExtra"].Value;
                            }
                            else
                            {
                                chargeExtra = "";
                            }

                            AddOffset(
                                int.Parse(m.Groups["StartScan"].Value),
                                int.Parse(m.Groups["EndScan"].Value),
                                int.Parse(m.Groups["ChargeState"].Value), dtaStartPos, chargeExtra);
                            dtaCount++;
                            if (dtaCount % 3000 == 0)
                            {
                                OnDTAScanUpdate(currentPosition / (double)fileLength);
                            }
                        }

                        dataLine = tr.ReadLine();
                    }

                    tr.Close();
                }

                public event dtaScanProgressEventHandler dtaScanProgress;

                public delegate void dtaScanProgressEventHandler(double fractionDone);

                private void OnDTAScanUpdate(double fractionDone)
                {
                    dtaScanProgress?.Invoke(fractionDone);
                }
            }

            ~DTAFileInformation()
            {
                m_DTAFileIntensities = null;
            }
        }

        internal class PeptideIntensities : FragmentInfo
        {
            private static FileStream m_FileStream;
            private NeutralLossList m_NeutralLoss;

            private double m_ParentMH;
            private int m_ParentCS;
            private static CalcNeutralLosses NLCalc;

            public PeptideIntensities(FileStream dtaFileStream)
            {
                if (m_FileStream is null)
                {
                    m_FileStream = dtaFileStream;
                }
                else if ((m_FileStream.Name ?? "") != (dtaFileStream.Name ?? ""))
                {
                    m_FileStream = dtaFileStream;
                }

                NLCalc ??= new CalcNeutralLosses();
            }

            public NeutralLossList NeutralLosses => m_NeutralLoss;

            public double ParentMZ => (m_ParentMH - 1.0d + m_ParentCS) / m_ParentCS;

            public void GetIntensitiesFromDTA(long startOffset)
            {
                var maxIntensity = 0.0;
                var sr = new StreamReader(m_FileStream);
                Clear();
                sr.BaseStream.Seek(startOffset, SeekOrigin.Begin);

                // Read the header line
                sr.ReadLine();

                //var HeaderLine = new Regex(@"^=+\s+\""\S+\.(?<Scan>\d+)\.\d+\.\d+\.dta");
                //if (HeaderLine.IsMatch(s))
                //{
                //    var headerLineMatch = HeaderLine.Match(s);
                //    m_ScanNum = int.Parse(headerLineMatch.Groups["Scan"].Value);
                //}

                var precursorLine = sr.ReadLine();
                var ParentLine = new Regex(@"^(?<ParentMass>\d+\.*\d*)\s+(?<ChargeState>\d+)");
                if (Regex.IsMatch(precursorLine, @"^\S+"))
                {
                    var parentLineMatch = ParentLine.Match(precursorLine);
                    m_ParentMH = double.Parse(parentLineMatch.Groups["ParentMass"].Value);
                    m_ParentCS = int.Parse(parentLineMatch.Groups["ChargeState"].Value);
                }

                var LineMatch = new Regex(@"^(?<Mass>\d+\.\d+)\s(?<Intensity>\d+\.\d+)");

                var dataLine = sr.ReadLine();
                while (dataLine != null && LineMatch.IsMatch(dataLine))
                {
                    var m = LineMatch.Match(dataLine);
                    var tmpMass = double.Parse(m.Groups["Mass"].Value);
                    var tmpIntensity = double.Parse(m.Groups["Intensity"].Value);
                    if (tmpIntensity > maxIntensity)
                    {
                        maxIntensity = tmpIntensity;
                    }

                    Add(tmpMass, tmpIntensity);

                    dataLine = sr.ReadLine();
                }

                NormalizeIntensities();
            }

            public void GetNeutralLosses(double massTol = 0.7)
            {
                m_NeutralLoss = NLCalc.CalculateNeutralLosses(this, massTol);
            }
        }
    }
}