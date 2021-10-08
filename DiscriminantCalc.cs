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
                    _m_dtaFileInfo.OffsetProgress -= OffsetLoadingProgHandler;
                }

                _m_dtaFileInfo = value;
                if (_m_dtaFileInfo != null)
                {
                    _m_dtaFileInfo.OffsetProgress += OffsetLoadingProgHandler;
                }
            }
        }

        private PeptideIntensities m_CachedScanInfo;
        private int m_CachedScanNum;
        private string m_CachedScanKey;
        private int m_CachedCS;
        private OutputNLIFile m_NLIDumper;
        private double m_MassTol;
        private string m_dtaFilepath;
        private string m_Version;
        private bool m_NoDTAs;

        public event ProgressUpdateEventHandler ProgressUpdate;

        public delegate void ProgressUpdateEventHandler(string TaskDescription, double fractionDone);

        private void OnProgressUpdate(string taskDescription, double fractionDone)
        {
            ProgressUpdate?.Invoke(taskDescription, fractionDone);
        }

        private void OffsetLoadingProgHandler(double fractionDone)
        {
            OnProgressUpdate("Loading .dta File Locations (" + (fractionDone * 100d).ToString("0.0") + "% Completed)", fractionDone);
        }

        public DiscriminantCalc(string dtaFilePath)
        {
            var fi = new FileInfo(dtaFilePath);
            string rootFileName;
            rootFileName = fi.Name;
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

        private void ConfigureDiscriminantCalc()
        {
            if (!m_NoDTAs)
            {
                m_dtaFileInfo.Configure();
            }
        }

        private string GetVersionString()
        {
            return m_Version;
        }

        private void CloseOut()
        {
            if (m_dtaFileInfo is object)
            {
                m_dtaFileInfo.Close();
                // m_dtaFileInfo = Nothing
            }

            if (m_NLIDumper is object)
            {
                m_NLIDumper.CloseNLIWriter();
            }
        }

        // --------------------
        // Version of the module (be sure to update)
        //
        public string version()
        {
            return "1.2";
        }

        // --------------------
        // Calculates MScore for given peptide record
        // using fragment intensities from
        // concatenated dta file
        //
        private double GetPeptideMScore(string PeptideSeq, int StartScanNumber, int EndScanNumber, int ChargeState)
        {
            string ScanKey = StartScanNumber.ToString() + "." + EndScanNumber.ToString();
            if (m_NoDTAs)
                return 10d;
            if (m_dtaFileInfo is null)
            {
                m_dtaFileInfo = new DTAFileInformation(m_dtaFilepath);
            }

            m_MassTol = 0.7d;
            if (ChargeState != m_CachedCS | (ScanKey ?? "") != (m_CachedScanKey ?? ""))
            {
                m_dtaFileInfo.GetDTAFileIntensities(StartScanNumber, EndScanNumber, ChargeState);
                if (m_dtaFileInfo.DTAScanInfo.FragmentList.Count == 0)
                {
                    return 10d;
                }

                if ((ScanKey ?? "") != (m_CachedScanKey ?? ""))
                {
                    m_dtaFileInfo.DTAScanInfo.GetNeutralLosses(m_MassTol);
                    m_NLIDumper.MakeNLIEntry(StartScanNumber, m_dtaFileInfo.DTAScanInfo.NeutralLosses);
                }

                m_CachedScanInfo = m_dtaFileInfo.DTAScanInfo;
                m_CachedScanNum = StartScanNumber;
                m_CachedScanKey = ScanKey;
                m_CachedCS = ChargeState;
            }

            double mscore = CalculateMScore(PeptideSeq, ChargeState, m_CachedScanInfo);
            return mscore;
        }

        private double CalculateMScore(string peptideSequence, int peptideChargeState, PeptideIntensities scanIntensities)
        {
            var match = default(double);
            double mScore;
            bool noModsInSequence;
            TheoreticalFragmentInfo theoFrags;
            double mTol = m_MassTol;
            if (scanIntensities.FragmentList.Count == 0)
            {
                return 10d;
            }

            if (Regex.IsMatch(peptideSequence.ToUpper(), "[BZXOJU*@#!$%^&]"))
            {
                return 10d;
            }
            else if (peptideSequence.Length <= 5)
            {
                return 10d;
            }
            else
            {
                noModsInSequence = true;
            }

            // Get values for +1 charge state
            theoFrags = new TheoreticalFragmentInfo(peptideSequence, 1);

            // Check and score +1 possibilities for B-Ions
            match += HashScanner(scanIntensities, theoFrags.BIons, mTol, 1);

            // Check and score +1 possibilities for Y-Ions
            match += HashScanner(scanIntensities, theoFrags.YIons, mTol, 1);
            if (peptideChargeState > 2)
            {
                // get values for +2 charge state
                theoFrags = new TheoreticalFragmentInfo(peptideSequence, 2);

                // Check and score +2 possibilities for B-Ions
                match += HashScanner(scanIntensities, theoFrags.BIons, mTol, 2);

                // Check and score +2 possibilities for Y-Ions
                match += HashScanner(scanIntensities, theoFrags.YIons, mTol, 2);
            }

            if (noModsInSequence)
            {
                mScore = Math.Round(10d + match, 2);
            }
            else
            {
                mScore = 10d;
            }

            return mScore;
        }

        private double HashScanner(PeptideIntensities peptideRecord, List<TheoreticalFragmentInfo.Fragment> theoFrags, double massTol, int CSToCheck)
        {
            double tmpObsMass;
            int maxObsRecord = peptideRecord.FragmentList.Count;
            if (maxObsRecord == 0)
            {
                return 0.0d;
            }

            int maxObsIndex = maxObsRecord - 1;
            double match;

            int obsCount = 0;
            match = 0d;
            tmpObsMass = peptideRecord.GetMass(obsCount);
            foreach (var theoFrag in theoFrags)
            {
                double tmpTheoMass = theoFrag.Mass;
                if (tmpObsMass > tmpTheoMass + massTol & obsCount <= maxObsIndex)
                {
                }
                else
                {
                    while (tmpObsMass < tmpTheoMass - massTol & obsCount < maxObsIndex)
                    {
                        obsCount += 1;
                        tmpObsMass = peptideRecord.GetMass(obsCount);
                    }

                    if (tmpObsMass < tmpTheoMass + massTol & !(tmpTheoMass > tmpObsMass + massTol) & obsCount <= maxObsIndex)

                    {
                        if (CSToCheck >= 2)
                        {
                            match += peptideRecord.GetNormalizedIntensity(obsCount);
                        }
                        else
                        {
                            match += peptideRecord.GetNormalizedIntensity(obsCount) * theoFrag.Intensity;
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
            private string m_dtaFilePath;

            public event OffsetProgressEventHandler OffsetProgress;

            public delegate void OffsetProgressEventHandler(double fractionDone);

            public DTAFileInformation(string dtaFilePath)
            {
                // Scan the dta file for dta file boundaries
                m_Offsets = new dtaFileOffsets();
                m_dtaFilePath = dtaFilePath;
            }

            public void Configure()
            {
                m_Offsets.LoadOffsetsFromDTAFile(m_dtaFilePath);
                m_dtaStream = new FileStream(m_dtaFilePath, FileMode.Open);
                m_DTAFileIntensities = new PeptideIntensities(m_dtaStream);
            }

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
                long fileOffset = m_Offsets.get_GetOffset(StartScanNumber, EndScanNumber, ChargeState);
                if (fileOffset > 0L)
                {
                    {
                        var withBlock = m_DTAFileIntensities;
                        withBlock.GetIntensitiesFromDTA(fileOffset);
                    }
                }
                else if (m_DTAFileIntensities is object)
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
                private Dictionary<string, long> mOffsets = new Dictionary<string, long>();
                private string m_FilePath;

                public long get_GetOffset(int StartScanNumber, int EndScanNumber, int ChargeState)
                {
                    string keyName = StartScanNumber.ToString() + "." + EndScanNumber.ToString() + "." + ChargeState.ToString();
                    long offset;
                    if (mOffsets.TryGetValue(keyName, out offset))
                    {
                        return offset;
                    }

                    throw new KeyNotFoundException("Offsets dictionary does not have key " + keyName);
                }

                public void AddOffset(int StartScanNumber, int EndScanNumber, int ChargeState, long StartOffset, string ChargeExtra = "")
                {
                    string keyName = StartScanNumber.ToString() + "." + EndScanNumber.ToString() + "." + ChargeState.ToString();
                    if (ChargeExtra.Length > 0)
                    {
                        keyName += "_" + ChargeExtra;
                    }

                    mOffsets.Add(keyName, StartOffset);
                }

                public int LoadOffsetsFromDTAFile(string dtaFilePath)
                {
                    var fi = new FileInfo(dtaFilePath);
                    TextReader tr;
                    string s;
                    int lineEndCharCount;
                    var currPos = default(long);
                    long dtaStartPos;
                    string chargeExtra;
                    var dtaCount = default(int);
                    long fileLength;
                    var r = new Regex(@"^\s*[=]{5,}\s+\""(?<rootname>.+)\.(?<StartScan>\d+)\.(?<EndScan>\d+)\.(?<ChargeBlock>(?<ChargeState>\d+)[^0-9]?(?<ChargeExtra>\S*))\.(?<FileType>.+)\""\s+[=]{5,}\s*$", RegexOptions.CultureInvariant | RegexOptions.Compiled);

                    var lineMatch = new Regex("^===*");
                    Match m;
                    if (fi.Exists)
                    {
                        lineEndCharCount = SequestFileExtractor.LineEndCharacterCount(fi);
                        fileLength = fi.Length;
                        tr = fi.OpenText();
                        s = tr.ReadLine();
                        while (s is object)
                        {
                            currPos += s.Length + lineEndCharCount;
                            if (lineMatch.IsMatch(s))
                            {
                                dtaStartPos = currPos - s.Length - lineEndCharCount;
                                m = r.Match(s);
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
                                dtaCount += 1;
                                if (dtaCount % 3000 == 0)
                                {
                                    Debug.WriteLine(dtaCount);
                                    OnDTAScanUpdate(currPos / (double)fileLength);
                                }
                            }

                            s = tr.ReadLine();
                        }

                        tr.Close();
                    }

                    return dtaCount;
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
            private double m_BFNLT;
            private int m_IsPoorSpec;
            private double m_ParentMH;
            private int m_ScanNum;
            private int m_ParentCS;
            private static CalcNeutralLosses NLCalc;

            public PeptideIntensities(FileStream dtaFileStream) : base()
            {
                if (m_FileStream is null)
                {
                    m_FileStream = dtaFileStream;
                }
                else if ((m_FileStream.Name ?? "") != (dtaFileStream.Name ?? ""))
                {
                    m_FileStream = dtaFileStream;
                }

                if (NLCalc is null)
                {
                    NLCalc = new CalcNeutralLosses();
                }
            }

            public NeutralLossList NeutralLosses => m_NeutralLoss;


            public double ParentMZ => (m_ParentMH - 1.0d + m_ParentCS) / m_ParentCS;

            public void GetIntensitiesFromDTA(long startOffset)
            {
                string s;
                double tmpMass;
                double tmpIntensity;
                var maxIntensity = default(double);
                var sr = new StreamReader(m_FileStream);
                Clear();
                sr.BaseStream.Seek(startOffset, SeekOrigin.Begin);
                s = sr.ReadLine();
                var HeaderLine = new Regex(@"^=+\s+\""\S+\.(?<Scan>\d+)\.\d+\.\d+\.dta");
                if (HeaderLine.IsMatch(s))
                {
                    var headerLineMatch = HeaderLine.Match(s);
                    m_ScanNum = int.Parse(headerLineMatch.Groups["Scan"].Value);
                }

                s = sr.ReadLine();
                var ParentLine = new Regex(@"^(?<ParentMass>\d+\.*\d*)\s+(?<ChargeState>\d+)");
                Match parentLineMatch;
                if (Regex.IsMatch(s, @"^\S+"))
                {
                    parentLineMatch = ParentLine.Match(s);
                    m_ParentMH = double.Parse(parentLineMatch.Groups["ParentMass"].Value);
                    m_ParentCS = int.Parse(parentLineMatch.Groups["ChargeState"].Value);
                }

                var LineMatch = new Regex(@"^(?<Mass>\d+\.\d+)\s(?<Intensity>\d+\.\d+)");
                Match m;
                s = sr.ReadLine();
                while (LineMatch.IsMatch(s))
                {
                    m = LineMatch.Match(s);
                    tmpMass = double.Parse(m.Groups["Mass"].Value);
                    tmpIntensity = double.Parse(m.Groups["Intensity"].Value);
                    if (tmpIntensity > maxIntensity)
                    {
                        maxIntensity = tmpIntensity;
                    }

                    Add(tmpMass, tmpIntensity);
                    s = sr.ReadLine();
                    if (s is null)
                        break;
                }

                NormalizeIntensities();
            }

            public void GetNeutralLosses(double MassTol)
            {
                m_NeutralLoss = NLCalc.CalculateNeutralLosses(this, 0.7d);
            }
            }
    }
}