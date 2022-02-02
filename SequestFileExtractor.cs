using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Text.RegularExpressions;
using System.Threading;
using PRISM.Logging;
using SequestResultsProcessor.Containers;

namespace SequestResultsProcessor
{
    // This class extracts peptide identifications from SEQUEST .out files and stores
    // them in synopsis files (_syn.txt) and first hits files (_fht.txt)
    //
    // -------------------------------------------------------------------------------
    // Originally written in Perl by Gary Kiebel
    // Ported to VB.NET by Ken Auberry for the Department of Energy (PNNL, Richland, WA)
    // Coding for this class began in November 2004
    // -------------------------------------------------------------------------------
    //
    // Licensed under the Apache License, Version 2.0; you may not use this file except
    // in compliance with the License.  You may obtain a copy of the License at
    // http://www.apache.org/licenses/LICENSE-2.0

    // ReSharper disable once CheckNamespace
    public class SequestFileExtractor
    {
        // Const DEBUG_FLAG As Boolean = True

        private ConcatenatedOutFileProcessor mParserInstance;

        private ConcatenatedOutFileProcessor OutFileParser
        {
            [MethodImpl(MethodImplOptions.Synchronized)]
            get => mParserInstance;

            [MethodImpl(MethodImplOptions.Synchronized)]
            set
            {
                if (mParserInstance != null)
                {
                    mParserInstance.EndingTask -= TaskEndHandler;
                    mParserInstance.ProgressReport -= ImportProgressHandler;
                }

                mParserInstance = value;
                if (mParserInstance != null)
                {
                    mParserInstance.EndingTask += TaskEndHandler;
                    mParserInstance.ProgressReport += ImportProgressHandler;
                }
            }
        }

        private static Regex mExtraProteinLineMatcher;
        private static Regex mHitLineMatcher;
        private static Regex mHitLineMatcherNoReference;
        private static Regex mTopProteinsMatcher;
        private static Regex mHeaderMassMatcher;
        private static Regex mDataBlockDelimiterMatcher;

        private enum eHitMatchType
        {
            NoMatch = 0,
            MatchWithProtein = 1,
            MatchWithoutProtein = 2
        }

        public SequestFileExtractor(StartupArguments startupArgs)
        {
            OutFileParser = new ConcatenatedOutFileProcessor(startupArgs);
            InitializeMatchers();
        }

        private void InitializeMatchers()
        {
            mExtraProteinLineMatcher ??= new Regex(@"^\s+\d*\s+(?<reference>\S+)\s*(?<description>.*)",
                RegexOptions.Compiled | RegexOptions.IgnoreCase | RegexOptions.Singleline);

            // Notes:
            // The ID column is not always present; thus the use of * in (?<idblock>...)*
            // The sf column is present in Bioworks 3.3.x, which has SEQUEST v.28 (rev. 13)
            // Sf stands for "Final Score" and is a number meant to reflect the strength of the
            // SEQUEST hit on a scale of 0 to 1. The number is created by running
            // SEQUEST Utilities -> Spectra (DTA) Tools -> Final Score.
            // Numbers above 0.7 are considered "good"
            // The multiorf column is only present when a peptide is in multiple proteins; thus the use of * in <?<multiorfblock>...)*
            mHitLineMatcher ??= new Regex(
                @"^\s*(?<hitnum>\d+)\.\s+" +
                @"(?<rankxc>\d+)\s*\/\s*" + @"(?<ranksp>\d+)\s+" +
                @"(?<idblock>(?<id>\d+)\s+)*" +
                @"(?<mhmass>\d+\.\d+)\s+" +
                @"(?<delcn>\d+\.\d+)\s+" +
                @"(?<xcorr>\d+\.\d+)\s+" +
                @"(?<sp>\d+\.\d+)\s+" +
                @"(?<sfblock>(?<sf>[0-9.]+)\s+)*" +
                @"(?<obsions>\d+)\s*\/\s*" +
                @"(?<theoions>\d+)\s+" +
                @"(?<reference>\S+)\s+" +
                @"(?<multiorfblock>\+(?<multiorf>\d+)\s+)*" +
                @"(?<sequence>\S+)",
                RegexOptions.Compiled | RegexOptions.IgnoreCase | RegexOptions.Singleline);

            // The following is used to match lines that do not have any text in the Reference (protein name) column
            mHitLineMatcherNoReference ??= new Regex(
                @"^\s*(?<hitnum>\d+)\.\s+" +
                @"(?<rankxc>\d+)\s*\/\s*" +
                @"(?<ranksp>\d+)\s+" +
                @"(?<idblock>(?<id>\d+)\s+)*" +
                @"(?<mhmass>\d+\.\d+)\s+" +
                @"(?<delcn>\d+\.\d+)\s+" +
                @"(?<xcorr>\d+\.\d+)\s+" +
                @"(?<sp>\d+\.\d+)\s+" +
                @"(?<sfblock>(?<sf>[0-9.]+)\s+)*" +
                @"(?<obsions>\d+)\s*\/\s*" +
                @"(?<theoions>\d+)\s+" +
                @"(?<sequence>\S+)",
                RegexOptions.Compiled | RegexOptions.IgnoreCase | RegexOptions.Singleline);

            mTopProteinsMatcher ??= new Regex(@"^\s+\d+\.\s+\d*\s+(?<reference>\S+)\s+(?<description>.+)$",
                RegexOptions.Compiled | RegexOptions.IgnoreCase | RegexOptions.Singleline);

            mHeaderMassMatcher ??= new Regex(@"mass\s+=\s+(?<HeaderMass>\d+\.\d+)", RegexOptions.Compiled);

            mDataBlockDelimiterMatcher ??= new Regex(@"^\s+---\s+-{3,}", RegexOptions.Compiled);
        }

        public void ProcessInputFile()
        {
            OutFileParser.ProcessInputFile();
        }

        public event ProgressReportEventHandler ProgressReport;

        public delegate void ProgressReportEventHandler(double fractionDone);

        public event StatusReportEventHandler StatusReport;

        public delegate void StatusReportEventHandler(string taskString);

        public event EndTaskEventHandler EndTask;

        public delegate void EndTaskEventHandler();

        private void UpdateProgress(string currentTask, double fractionDone)
        {
            StatusReport?.Invoke(currentTask);
            ProgressReport?.Invoke(fractionDone);
        }

        private void UpdateProgress(string currentTask, long currentFileCount, long totalFileCount)
        {
            double fractionDone;
            if (totalFileCount > 0L)
            {
                fractionDone = currentFileCount / (double)totalFileCount;
            }
            else
            {
                fractionDone = 0.0d;
            }

            UpdateProgress(currentTask, fractionDone);
        }

        private void TaskEndHandler()
        {
            EndTask?.Invoke();
        }

        private void ImportProgressHandler(string currentTask, long currentPosition, long totalSize)
        {
            UpdateProgress(currentTask, currentPosition, totalSize);
        }

        private class ConcatenatedOutFileProcessor
        {
            private StartupArguments mStartupArguments;
            private string mSourceFileFullPath;
            private ResultsStorage mResults;
            private OutputIRRFile mIRRDumper;
            private PRISM.Logging.FileLogger mLogger;
            private bool mStopProcessing;
            private const int RESULTS_DUMPING_INTERVAL = 200;

            #region  Progress Update Events
            internal event ProgressReportEventHandler ProgressReport;

            internal delegate void ProgressReportEventHandler(string currentTask, long currentPosition, long totalSize);

            internal event EndingTaskEventHandler EndingTask;

            internal delegate void EndingTaskEventHandler();

            private void EndTask()
            {
                EndingTask?.Invoke();
            }

            private void StartingNewTask(string taskName)
            {
                ProgressReport?.Invoke(taskName, 0L, 0L);
            }

            private void UpdateProgressCountingOuts(long currentFilePos, long totalFileSize)
            {
                const string statusString = "Counting .out Files... ";
                ProgressReport?.Invoke(statusString, currentFilePos, totalFileSize);
            }

            private void UpdateProgressExtracting(int currentOutFileCount, int totalOutFileCount)
            {
                var statusString = "Extracting from " + Path.GetFileNameWithoutExtension(mSourceFileFullPath) + " (File " + currentOutFileCount + " of " + totalOutFileCount + ")";
                ProgressReport?.Invoke(statusString, currentOutFileCount, totalOutFileCount);
            }

            #endregion

            public ConcatenatedOutFileProcessor(StartupArguments ProcessSettings)
            {
                InitializeSettings(ProcessSettings);
            }

            public void AbortProcessing()
            {
                mStopProcessing = true;
            }

            private bool AdvanceReaderUntilMatch(StreamReader reader, Regex reMatcher, string dataLine, out string matchingLine)
            {
                if (reMatcher.IsMatch(dataLine))
                {
                    matchingLine = dataLine;
                    return true;
                }

                while (true)
                {
                    if (!reader.EndOfStream)
                    {
                        var nextLine = reader.ReadLine();
                        var matchFound = reMatcher.IsMatch(nextLine ?? string.Empty);
                        if (matchFound)
                        {
                            matchingLine = nextLine;
                            return true;
                        }
                    }
                    else
                    {
                        matchingLine = string.Empty;
                        return false;
                    }
                }
            }

            private void DumpCachedResults(string tmpFHTPath, string tmpSynPath, List<ResultsStorage.OutputRecordIndex> FHTOutputIndexList, List<ResultsStorage.OutputRecordIndex> SynOutputIndexList)
            {
                mResults.ExportContents(ResultsStorage.OutputTypeList.FHT, mStartupArguments.FHTXCorrThreshold, false, tmpFHTPath, FHTOutputIndexList);
                mResults.ExportContents(ResultsStorage.OutputTypeList.Syn, mStartupArguments.SynXCorrThreshold, mStartupArguments.ExpandMultiORF, tmpSynPath, SynOutputIndexList);
                mResults.ClearResults();
            }

            private void InitializeSettings(StartupArguments ProcessSettings)
            {
                mStartupArguments = ProcessSettings;
                mSourceFileFullPath = mStartupArguments.InputFileFullPath;
                mStartupArguments = ProcessSettings;
                mLogger = new FileLogger(Path.Combine(mStartupArguments.DestinationDirectory, Path.GetFileNameWithoutExtension(mStartupArguments.LogFileName)));
            }

            /// <summary>
            /// Open a _out.txt file and extract the peptide info
            /// </summary>
            /// <remarks></remarks>
            public void ProcessInputFile()
            {
                var makeIRRFile = false;
                var extractorVersion = Assembly.GetExecutingAssembly().GetName().Version.ToString();
                mLogger.LogMessage(BaseLogger.LogLevels.INFO, "------ Peptide File Extractor, " + Path.GetFileName(Assembly.GetExecutingAssembly().Location) + " v" + extractorVersion);
                ProgressReport?.Invoke("Initializing", 0L, 0L);
                mResults = new ResultsStorage();
                var lstOutFilesProcessed = new SortedSet<string>();
                var FHTOutputIndexList = new List<ResultsStorage.OutputRecordIndex>();
                var SynOutputIndexList = new List<ResultsStorage.OutputRecordIndex>();
                var r_FileDelimiterMatcher = new Regex(@"^\s*[=]{5,}\s+\""(?<rootname>.+)\.(?<StartScan>\d+)\.(?<EndScan>\d+)\.(?<ChargeBlock>(?<ChargeState>\d+)[^0-9]?(?<ChargeExtra>\S*))\.(?<FileType>.+)\""\s+[=]{5,}\s*$", RegexOptions.CultureInvariant | RegexOptions.Compiled | RegexOptions.IgnoreCase);
                if (mStartupArguments.MakeIRRFile)
                {
                    mIRRDumper = new OutputIRRFile(mStartupArguments.RootFileName, mStartupArguments.DestinationDirectory);
                    makeIRRFile = true;
                }

                var fi = new FileInfo(mSourceFileFullPath);
                var removeDupMultiProteinRefs = mStartupArguments.RemoveDuplicatedMultiProteinRefs;
                var tmpFHTPath = Path.Combine(Path.GetDirectoryName(mSourceFileFullPath), "Tmp_FHT.txt");
                var tmpSynPath = Path.Combine(Path.GetDirectoryName(mSourceFileFullPath), "Tmp_Syn.txt");
                var tmpFHTProtPath = Path.Combine(Path.GetDirectoryName(mSourceFileFullPath), "Tmp_FHT_Prot.txt");
                var tmpSynProtPath = Path.Combine(Path.GetDirectoryName(mSourceFileFullPath), "Tmp_Syn_Prot.txt");
                var tmpFHTFI = new FileInfo(tmpFHTPath);
                var tmpSynFI = new FileInfo(tmpSynPath);
                var tmpFHTProtFI = new FileInfo(tmpFHTProtPath);
                var tmpSynProtFI = new FileInfo(tmpSynProtPath);
                if (tmpFHTFI.Exists)
                {
                    tmpFHTFI.Delete();
                }

                if (tmpSynFI.Exists)
                {
                    tmpSynFI.Delete();
                }

                if (tmpFHTProtFI.Exists)
                {
                    tmpFHTProtFI.Delete();
                }

                if (tmpSynProtFI.Exists)
                {
                    tmpSynProtFI.Delete();
                }

                if (!fi.Exists)
                {
                    mLogger.LogMessage(BaseLogger.LogLevels.ERROR, "-------- ERROR: '" + mSourceFileFullPath + "' apparently doesn't exist");
                    mLogger.LogMessage(BaseLogger.LogLevels.ERROR, "-------- Exiting program ---------");
                    EndTask();
                    return;
                }

                var totalOutFileCount = CountOutFiles();
                var currentOutFileCount = 0;
                if (totalOutFileCount == 0)
                {
                    mLogger.LogMessage(BaseLogger.LogLevels.ERROR, "-------- ERROR: '" + mStartupArguments.InputFileName + "' contained no concatenated .out files");
                    mLogger.LogMessage(BaseLogger.LogLevels.ERROR, "-------- Exiting program ---------");
                    EndTask();
                    return;
                }

                using (var reader = new StreamReader(new FileStream(fi.FullName, FileMode.Open, FileAccess.Read, FileShare.Read)))
                {
                    mLogger.LogMessage(BaseLogger.LogLevels.INFO, "Processing '" + mStartupArguments.InputFileName + "'");
                    while (!reader.EndOfStream)
                    {
                        var dataLine = reader.ReadLine();

                        var outFileFound = AdvanceReaderUntilMatch(reader, r_FileDelimiterMatcher, dataLine, out var matchingLine);
                        if (outFileFound)
                        {
                            // Check whether this .out file has already been processed
                            if (lstOutFilesProcessed.Contains(matchingLine))
                            {
                                // Skip this .out file
                                mLogger.LogMessage(BaseLogger.LogLevels.WARN, "Skipping duplicate .out file: " + matchingLine.Trim('='));
                                Console.WriteLine("Warning, skipping duplicate .out file: " + matchingLine.Trim('='));
                                outFileFound = false;
                            }
                            else
                            {
                                lstOutFilesProcessed.Add(matchingLine);
                            }
                        }

                        if (outFileFound)
                        {
                            // Increment the currentOutFile counter
                            currentOutFileCount++;
                            if (currentOutFileCount % 1000 == 0)
                            {
                                mLogger.LogMessage(BaseLogger.LogLevels.INFO, " ... " + currentOutFileCount.ToString().PadLeft(5) + " spectra processed");
                            }

                            ReadAndStoreOutFileResults(reader, matchingLine, r_FileDelimiterMatcher, makeIRRFile, removeDupMultiProteinRefs);
                        }

                        if (currentOutFileCount % RESULTS_DUMPING_INTERVAL == 0 || currentOutFileCount >= totalOutFileCount)
                        {
                            DumpCachedResults(tmpFHTPath, tmpSynPath, FHTOutputIndexList, SynOutputIndexList);
                            UpdateProgressExtracting(currentOutFileCount, totalOutFileCount);
                            GC.Collect();
                            GC.WaitForPendingFinalizers();
                            Thread.Sleep(100);
                        }
                    }
                }

                if (mResults.Count > 0)
                {
                    DumpCachedResults(tmpFHTPath, tmpSynPath, FHTOutputIndexList, SynOutputIndexList);
                }

                mLogger.LogMessage(BaseLogger.LogLevels.INFO, "Peak file '" + Path.GetFileName(tmpSynPath) + "' was generated");
                mLogger.LogMessage(BaseLogger.LogLevels.INFO, "Peak file '" + Path.GetFileName(tmpFHTPath) + "' was generated");
                mLogger.LogMessage(BaseLogger.LogLevels.INFO, "Peak file '" + Path.GetFileName(tmpSynPath) + "' contains " + SynOutputIndexList.Count.ToString().PadLeft(7) + " peptides " + "(XCorr Threshold was " + mStartupArguments.SynXCorrThreshold + " and subsequently filtered)");
                mLogger.LogMessage(BaseLogger.LogLevels.INFO, "Peak file '" + Path.GetFileName(tmpFHTPath) + "' contains " + FHTOutputIndexList.Count.ToString().PadLeft(7) + " peptides " + "(XCorr Threshold was " + mStartupArguments.FHTXCorrThreshold + ")");

                // Keys in these dictionaries are XCorr threshold; values are the number of peptides with an XCorr over the threshold
                if (!mStopProcessing)
                {
                    StartingNewTask("Sorting results...");
                    mLogger.LogMessage(BaseLogger.LogLevels.INFO, "Sorting peptides in Syn file");
                    var SynSummaryResults = mResults.SortPeptides(tmpSynPath, SynOutputIndexList, mStartupArguments.SynopsisFileFullPath);
                    mLogger.LogMessage(BaseLogger.LogLevels.INFO, "Sorting peptides in Fht file");
                    var FHTSummaryResults = mResults.SortPeptides(tmpFHTPath, FHTOutputIndexList, mStartupArguments.FirstHitsFullPath);
                    mLogger.LogMessage(BaseLogger.LogLevels.INFO, GetScoreSummary(SynSummaryResults, "all peptides"));
                    mLogger.LogMessage(BaseLogger.LogLevels.INFO, GetScoreSummary(FHTSummaryResults, "first hits only"));
                    mLogger.LogMessage(BaseLogger.LogLevels.INFO, "Synopsis File   '" + Path.GetFileName(mStartupArguments.SynopsisFileFullPath) + "' was generated");
                    mLogger.LogMessage(BaseLogger.LogLevels.INFO, "First Hits File '" + Path.GetFileName(mStartupArguments.FirstHitsFullPath) + "' was generated");
                }

                mLogger = null;
                if (mStartupArguments.MakeIRRFile)
                {
                    mIRRDumper.CloseIRRWriter();
                }

                mResults = null;
                EndTask();
            }

            private string GetScoreSummary(Dictionary<int, int> peptideCountsByXCorr, string statsDescription)
            {
                // "Scores (all peptides)     -> "
                // "Scores (first hits only)  -> "

                var msg = string.Format("Scores {0,-18} -> ", "(" + statsDescription + ")");
                var query = from item in peptideCountsByXCorr.Keys
                            orderby item
                            select item;

                foreach (var scoreThreshold in query)
                {
                    msg += string.Format("{0,7} peptides above {1}, ", peptideCountsByXCorr[scoreThreshold], scoreThreshold);
                }

                // Remove the trailing comma and space
                var trimmedMsg = msg.Trim();
                if (trimmedMsg.EndsWith(","))
                {
                    return trimmedMsg.Substring(0, trimmedMsg.Length - 1);
                }

                return trimmedMsg;
            }

            private int CountOutFiles()
            {
                var fi = new FileInfo(mSourceFileFullPath);
                var outFileCount = 0;
                long currentPosition = 0;
                long lineCount = 0;
                var lineEndCharCount = LineEndCharacterCount(fi);
                var r = new Regex("^===*", RegexOptions.Compiled);
                if (fi.Exists)
                {
                    TextReader tr = fi.OpenText();
                    var s = tr.ReadLine();
                    while (s != null)
                    {
                        lineCount++;
                        currentPosition += s.Length + lineEndCharCount;
                        if (r.IsMatch(s))
                            outFileCount++;
                        s = tr.ReadLine();
                        if (lineCount % 500L == 0L)
                        {
                            UpdateProgressCountingOuts(currentPosition, fi.Length);
                            if (mStopProcessing)
                                break;
                        }
                    }

                    UpdateProgressCountingOuts(currentPosition, fi.Length);
                    tr.Close();
                    return outFileCount;
                }

                return default;
            }

            private void ReadAndStoreOutFileResults(StreamReader reader, string fileHeaderLine, Regex fileDelimiterMatcher, bool makeIRRFile, bool removeDupMultiProteinRefs)
            {
                int currentStartScan;
                int currentEndScan;
                int currentCS;

                // Extract File Header info
                var fileHeaderMatch = fileDelimiterMatcher.Match(fileHeaderLine);
                if (fileHeaderMatch.Success)
                {
                    currentStartScan = int.Parse(fileHeaderMatch.Groups["StartScan"].Value);
                    currentEndScan = int.Parse(fileHeaderMatch.Groups["EndScan"].Value);
                    currentCS = int.Parse(fileHeaderMatch.Groups["ChargeState"].Value);
                }
                else
                {
                    currentStartScan = 0;
                    currentEndScan = 0;
                    currentCS = 0;
                }

                var tmpMultiProteinRefs = new List<string>();

                // Wait until we see the measured mass show up in the header
                var foundHeaderMass = AdvanceReaderUntilMatch(reader, mHeaderMassMatcher, "", out var matchingLine);
                if (!foundHeaderMass)
                    return;

                // grab the header mass value
                var headerMassMatch = mHeaderMassMatcher.Match(matchingLine);
                var currentHeaderMass = double.Parse(headerMassMatch.Groups["HeaderMass"].Value);

                // Wait until we see the dashed line underneath the headings for the data block

                var foundDataBlock = AdvanceReaderUntilMatch(reader, mDataBlockDelimiterMatcher, "", out _);
                if (!foundDataBlock || reader.EndOfStream)
                    return;

                // Read the first line of the data block
                var dataLine = reader.ReadLine();

                // As long as we keep seeing hit lines, keep grabbing them (also, allow one blank line)
                // (to separate that VERY last hit that creeps in)

                var matchType = eHitMatchType.MatchWithProtein;
                while (matchType != eHitMatchType.NoMatch && dataLine != null)
                {
                    var dataLineMatch = mHitLineMatcher.Match(dataLine);
                    if (dataLineMatch.Success)
                    {
                        matchType = eHitMatchType.MatchWithProtein;
                    }
                    else
                    {
                        dataLineMatch = mHitLineMatcherNoReference.Match(dataLine);
                        if (dataLineMatch.Success)
                        {
                            matchType = eHitMatchType.MatchWithoutProtein;
                        }
                        else
                        {
                            matchType = eHitMatchType.NoMatch;
                        }
                    }

                    if (matchType != eHitMatchType.NoMatch)
                    {
                        var currentPeptide = new PeptideHitEntry
                        {
                            HitNum = int.Parse(dataLineMatch.Groups["hitnum"].Value),
                            MH = Math.Round(double.Parse(dataLineMatch.Groups["mhmass"].Value), 5),
                            DelCn = double.Parse(dataLineMatch.Groups["delcn"].Value),
                            XCorr = double.Parse(dataLineMatch.Groups["xcorr"].Value),
                            Sp = double.Parse(dataLineMatch.Groups["sp"].Value)
                        };

                        if (matchType == eHitMatchType.MatchWithProtein)
                        {
                            currentPeptide.Reference = dataLineMatch.Groups["reference"].Value;
                            if (dataLineMatch.Groups["multiorf"].Length > 0)
                            {
                                currentPeptide.MultiProteinCount = int.Parse(dataLineMatch.Groups["multiorf"].Value);
                            }
                            else
                            {
                                currentPeptide.MultiProteinCount = 0;
                            }
                        }
                        else
                        {
                            currentPeptide.Reference = string.Empty;
                            currentPeptide.MultiProteinCount = 0;
                        }

                        currentPeptide.Peptide = dataLineMatch.Groups["sequence"].Value;
                        currentPeptide.RankSp = int.Parse(dataLineMatch.Groups["ranksp"].Value);
                        currentPeptide.RankXc = int.Parse(dataLineMatch.Groups["rankxc"].Value);

                        currentPeptide.ObsIons = int.Parse(dataLineMatch.Groups["obsions"].Value);
                        currentPeptide.PossIons = int.Parse(dataLineMatch.Groups["theoions"].Value);
                        currentPeptide.XcRatio = 1d;

                        currentPeptide.StartScanNum = currentStartScan;
                        currentPeptide.EndScanNum = currentEndScan;
                        currentPeptide.ScanCount = currentEndScan - currentStartScan + 1;
                        currentPeptide.ChargeState = currentCS;

                        if (makeIRRFile)
                        {
                            mIRRDumper.MakeIRREntry(currentStartScan, currentCS, currentPeptide.RankXc, currentPeptide.ObsIons, currentPeptide.PossIons);
                        }

                        // Read the next line
                        if (!reader.EndOfStream)
                        {
                            dataLine = reader.ReadLine();

                            // Look for multi-protein hit lines, but make sure to exclude those top scoring protein lines underneath them
                            while (true)
                            {
                                var extraProteinLineMatch = mExtraProteinLineMatcher.Match(dataLine ?? string.Empty);

                                if (extraProteinLineMatch.Success && !mTopProteinsMatcher.IsMatch(dataLine ?? string.Empty))
                                {
                                    var tmpMultiProteinRef = extraProteinLineMatch.Groups["reference"].Value;

                                    if (tmpMultiProteinRefs.Contains(tmpMultiProteinRef) && removeDupMultiProteinRefs && !string.Equals(tmpMultiProteinRef, currentPeptide.Reference, StringComparison.OrdinalIgnoreCase))
                                    {
                                        // Skip this protein name
                                    }
                                    else
                                    {
                                        // Store this protein name
                                        currentPeptide.AddMultiProteinRef(tmpMultiProteinRef);
                                        tmpMultiProteinRefs.Add(tmpMultiProteinRef);
                                    }

                                    if (!reader.EndOfStream)
                                    {
                                        dataLine = reader.ReadLine();
                                    }
                                    else
                                    {
                                        dataLine = null;
                                        break;
                                    }
                                }
                                else
                                {
                                    break;
                                }
                            }
                        }

                        currentPeptide.MultiProteinCount = tmpMultiProteinRefs.Count;

                        // Store the results for this peptide
                        mResults.AddPeptideResults(currentHeaderMass, currentPeptide);
                        tmpMultiProteinRefs.Clear();
                    }
                }
            }
        }

        /// <summary>
        /// This method reads the input file one byte at a time, looking for the first occurrence of character code 10 or 13 (aka CR or LF)
        /// When found, the next byte is examined
        /// If the next byte is also character code 10 or 13, the line terminator is assumed to be 2 bytes; if not found, it is assumed to be one byte
        /// </summary>
        /// <param name="fi"></param>
        /// <returns>1 if a one-byte line terminator; 2 if a two-byte line terminator</returns>
        internal static int LineEndCharacterCount(FileInfo fi)
        {
            var endCount = 1; // Initially assume a one-byte line terminator

            if (!fi.Exists)
                return endCount;

            TextReader tr = fi.OpenText();
            for (var counter = 1; counter <= fi.Length; counter++)
            {
                var testCode = tr.Read();
                if (testCode is 10 or 13)
                {
                    var testCode2 = tr.Read();
                    if (testCode2 is 10 or 13)
                    {
                        endCount = 2;
                        break;
                    }

                    endCount = 1;
                    break;
                }
            }

            tr.Close();

            return endCount;
        }

        public class StartupArguments
        {
            private readonly string mSourceDirectory;
            private string mDestinationDirectory;
            private string mInputFileName;
            private string mSynopsisFileName;
            private string mFirstHitsFileName;
            private string mDTAFileName;
            private string mLogFileName;

            /// <summary>
            /// Parameterless constructor
            /// </summary>
            // ReSharper disable once UnusedMember.Global
            public StartupArguments()
            {
            }

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="SourceDirectory"></param>
            /// <param name="RootFileName"></param>
            // ReSharper disable once UnusedMember.Global
            public StartupArguments(string SourceDirectory, string RootFileName)
            {
                mSourceDirectory = SourceDirectory;
                this.RootFileName = RootFileName;
            }

            private bool CheckFileExists(string filePath)
            {
                var fi = new FileInfo(filePath);
                return fi.Exists;
            }

            public string DestinationDirectory
            {
                get
                {
                    if (mDestinationDirectory is null)
                    {
                        return mSourceDirectory;
                    }

                    return mDestinationDirectory;
                }

                set => mDestinationDirectory = value;
            }

            // ReSharper disable once UnusedMember.Global
            public bool CatOutFileExists => CheckFileExists(InputFileFullPath);

            // ReSharper disable once UnusedMember.Global
            public bool CatDTAFileExists => CheckFileExists(DTAFileFullPath);

            public string RootFileName { get; set; }

            public string InputFileName
            {
                get
                {
                    if (mInputFileName is null)
                    {
                        return RootFileName + "_out.txt";
                    }

                    return mInputFileName;
                }

                set => mInputFileName = value;
            }

            public string SynopsisFileName
            {
                get
                {
                    if (mSynopsisFileName is null)
                    {
                        return RootFileName + "_syn.txt";
                    }

                    return mSynopsisFileName;
                }

                set => mSynopsisFileName = value;
            }

            public string FirstHitsFileName
            {
                get
                {
                    if (mFirstHitsFileName is null)
                    {
                        return RootFileName + "_fht.txt";
                    }

                    return mFirstHitsFileName;
                }

                set => mFirstHitsFileName = value;
            }

            public string DTAFileName
            {
                get
                {
                    if (mDTAFileName is null)
                    {
                        return RootFileName + "_dta.txt";
                    }

                    return mDTAFileName;
                }

                set => mDTAFileName = value;
            }

            public string LogFileName
            {
                get
                {
                    if (mLogFileName is null)
                    {
                        return RootFileName + "_log.txt";
                    }

                    return mLogFileName;
                }

                set => mLogFileName = value;
            }

            public bool MakeIRRFile { get; set; }
            public bool ExpandMultiORF { get; set; } = true;
            public double FHTXCorrThreshold { get; set; }
            public double SynXCorrThreshold { get; set; } = 1.5d;
            public bool RemoveDuplicatedMultiProteinRefs { get; set; }

            public string InputFileFullPath => Path.Combine(DestinationDirectory, InputFileName);

            public string DTAFileFullPath => Path.Combine(DestinationDirectory, DTAFileName);

            public string FirstHitsFullPath => Path.Combine(DestinationDirectory, FirstHitsFileName);

            public string SynopsisFileFullPath => Path.Combine(DestinationDirectory, SynopsisFileName);

            public string LogFileFullPath => Path.Combine(DestinationDirectory, LogFileName);
        }
    }
}