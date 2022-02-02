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
        private DTAFileInformation mDtaFileInfo;

        private DTAFileInformation DtaFileInfo
        {
            [MethodImpl(MethodImplOptions.Synchronized)]
            get => mDtaFileInfo;

            [MethodImpl(MethodImplOptions.Synchronized)]
            set
            {
                if (mDtaFileInfo != null)
                {
                    mDtaFileInfo.OffsetProgress -= OffsetLoadingProgressHandler;
                }

                mDtaFileInfo = value;
                if (mDtaFileInfo != null)
                {
                    mDtaFileInfo.OffsetProgress += OffsetLoadingProgressHandler;
                }
            }
        }

        private PeptideIntensities mCachedScanInfo;
        private string mCachedScanKey;
        private int mCachedCS;
        private readonly OutputNLIFile mNLIDumper;
        private double mMassTol;
        private readonly string mDtaFilepath;
        private readonly string mVersion;
        private readonly bool mNoDTAs;

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
            mDtaFilepath = dtaFilePath;
            if (fi.Exists)
            {
                DtaFileInfo = new DTAFileInformation(dtaFilePath);
                mNoDTAs = false;
                mNLIDumper = new OutputNLIFile(rootFileName, fi.DirectoryName);
            }
            else
            {
                mNoDTAs = true;
            }

            mVersion = FileVersionInfo.GetVersionInfo(Assembly.GetExecutingAssembly().Location).FileVersion;
        }

        /// <summary>
        /// Version of the module
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public string Version()
        {
            return mVersion;
        }

        // --------------------
        // Calculates MScore for given peptide record
        // using fragment intensities from
        // concatenated dta file
        [Obsolete("Unused")]
        private double GetPeptideMScore(string PeptideSeq, int StartScanNumber, int EndScanNumber, int ChargeState)
        {
            var ScanKey = StartScanNumber + "." + EndScanNumber;
            if (mNoDTAs)
                return 10d;

            DtaFileInfo ??= new DTAFileInformation(mDtaFilepath);

            mMassTol = 0.7d;
            if (ChargeState != mCachedCS || ScanKey != (mCachedScanKey ?? string.Empty))
            {
                DtaFileInfo.GetDTAFileIntensities(StartScanNumber, EndScanNumber, ChargeState);
                if (DtaFileInfo.DTAScanInfo.FragmentList.Count == 0)
                {
                    return 10d;
                }

                if (ScanKey != (mCachedScanKey ?? string.Empty))
                {
                    DtaFileInfo.DTAScanInfo.GetNeutralLosses(mMassTol);
                    mNLIDumper.MakeNLIEntry(StartScanNumber, DtaFileInfo.DTAScanInfo.NeutralLosses);
                }

                mCachedScanInfo = DtaFileInfo.DTAScanInfo;
                mCachedScanKey = ScanKey;
                mCachedCS = ChargeState;
            }

            return CalculateMScore(PeptideSeq, ChargeState, mCachedScanInfo);
        }

        private double CalculateMScore(string peptideSequence, int peptideChargeState, FragmentInfo scanIntensities)
        {
            var match = 0.0;

            var mTol = mMassTol;
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

            const bool noModsInSequence = true;

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

        private double HashScanner(FragmentInfo peptideRecord, List<TheoreticalFragmentInfo.Fragment> theoreticalFragments, double massTol, int chargeStateToCheck)
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
                    continue;
                }

                while (tmpObsMass < theoreticalMass - massTol && obsCount < maxObsIndex)
                {
                    obsCount++;
                    tmpObsMass = peptideRecord.GetMass(obsCount);
                }

                if (tmpObsMass >= theoreticalMass + massTol || theoreticalMass > tmpObsMass + massTol || obsCount > maxObsIndex)
                {
                    continue;
                }

                if (chargeStateToCheck >= 2)
                {
                    match += peptideRecord.GetNormalizedIntensity(obsCount);
                }
                else
                {
                    match += peptideRecord.GetNormalizedIntensity(obsCount) * fragment.Intensity;
                }
            }

            return match;
        }

        private class DTAFileInformation
        {
            private DTAFileOffsets mOffsets;

            private DTAFileOffsets Offsets
            {
                [MethodImpl(MethodImplOptions.Synchronized)]
                get => mOffsets;

                [MethodImpl(MethodImplOptions.Synchronized)]
                set
                {
                    if (mOffsets != null)
                    {
                        mOffsets.DtaScanProgress -= OffsetLoadingProgressHandler;
                    }

                    mOffsets = value;
                    if (mOffsets != null)
                    {
                        mOffsets.DtaScanProgress += OffsetLoadingProgressHandler;
                    }
                }
            }

            private PeptideIntensities mDTAFileIntensities;
            private FileStream mDtaStream;
            private readonly string mDtaFilePath;

            public event OffsetProgressEventHandler OffsetProgress;

            public delegate void OffsetProgressEventHandler(double fractionDone);

            public DTAFileInformation(string dtaFilePath)
            {
                // Scan the dta file for dta file boundaries
                Offsets = new DTAFileOffsets();
                mDtaFilePath = dtaFilePath;
            }

            // ReSharper disable once UnusedMember.Local
            public void Configure()
            {
                Offsets.LoadOffsetsFromDTAFile(mDtaFilePath);
                mDtaStream = new FileStream(mDtaFilePath, FileMode.Open);
                mDTAFileIntensities = new PeptideIntensities(mDtaStream);
            }

            // ReSharper disable once UnusedMember.Local
            public void Close()
            {
                mDtaStream.Close();
            }

            public PeptideIntensities DTAScanInfo => mDTAFileIntensities;

            private void OffsetLoadingProgressHandler(double fractionDone)
            {
                OffsetProgress?.Invoke(fractionDone);
            }

            public void GetDTAFileIntensities(int StartScanNumber, int EndScanNumber, int ChargeState)
            {
                var fileOffset = Offsets.GetOffset(StartScanNumber, EndScanNumber, ChargeState);
                if (fileOffset > 0L)
                {
                    mDTAFileIntensities.GetIntensitiesFromDTA(fileOffset);
                }
                else if (mDTAFileIntensities != null)
                {
                    mDTAFileIntensities.Clear();
                }
                else
                {
                    mDTAFileIntensities = new PeptideIntensities(mDtaStream);
                }
            }

            private class DTAFileOffsets
            {
                private readonly Dictionary<string, long> mOffsets = new();

                public long GetOffset(int StartScanNumber, int EndScanNumber, int ChargeState)
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

                public event dtaScanProgressEventHandler DtaScanProgress;

                public delegate void dtaScanProgressEventHandler(double fractionDone);

                private void OnDTAScanUpdate(double fractionDone)
                {
                    DtaScanProgress?.Invoke(fractionDone);
                }
            }

            ~DTAFileInformation()
            {
                mDTAFileIntensities = null;
            }
        }

        internal class PeptideIntensities : FragmentInfo
        {
            private static FileStream mFileStream;
            private NeutralLossList mNeutralLoss;

            private double mParentMH;
            private int mParentCS;
            private static CalcNeutralLosses NLCalc;

            public PeptideIntensities(FileStream dtaFileStream)
            {
                if (mFileStream is null)
                {
                    mFileStream = dtaFileStream;
                }
                else if ((mFileStream.Name ?? "") != (dtaFileStream.Name ?? ""))
                {
                    mFileStream = dtaFileStream;
                }

                NLCalc ??= new CalcNeutralLosses();
            }

            public NeutralLossList NeutralLosses => mNeutralLoss;

            public double ParentMZ => (mParentMH - 1.0d + mParentCS) / mParentCS;

            public void GetIntensitiesFromDTA(long startOffset)
            {
                var maxIntensity = 0.0;
                var reader = new StreamReader(mFileStream);
                Clear();
                reader.BaseStream.Seek(startOffset, SeekOrigin.Begin);

                // Read the header line
                reader.ReadLine();

                //var HeaderLine = new Regex(@"^=+\s+\""\S+\.(?<Scan>\d+)\.\d+\.\d+\.dta");
                //if (HeaderLine.IsMatch(s))
                //{
                //    var headerLineMatch = HeaderLine.Match(s);
                //    mScanNum = int.Parse(headerLineMatch.Groups["Scan"].Value);
                //}

                var precursorLine = reader.ReadLine();
                var ParentLine = new Regex(@"^(?<ParentMass>\d+\.*\d*)\s+(?<ChargeState>\d+)");
                if (Regex.IsMatch(precursorLine, @"^\S+"))
                {
                    var parentLineMatch = ParentLine.Match(precursorLine);
                    mParentMH = double.Parse(parentLineMatch.Groups["ParentMass"].Value);
                    mParentCS = int.Parse(parentLineMatch.Groups["ChargeState"].Value);
                }

                var LineMatch = new Regex(@"^(?<Mass>\d+\.\d+)\s(?<Intensity>\d+\.\d+)");

                var dataLine = reader.ReadLine();
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

                    dataLine = reader.ReadLine();
                }

                NormalizeIntensities();
            }

            public void GetNeutralLosses(double massTol = 0.7)
            {
                mNeutralLoss = NLCalc.CalculateNeutralLosses(this, massTol);
            }
        }
    }
}