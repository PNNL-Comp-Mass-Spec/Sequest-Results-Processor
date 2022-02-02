using System.Text.RegularExpressions;

namespace SequestResultsProcessor
{
    public class CleavageStateCalculator
    {
        private readonly Regex mFullyTryptic;
        private readonly Regex mFullyTrypticProlineCheck;
        private readonly Regex mFullyTrypticModKR;
        private readonly Regex mPeptideSpanningProtein;
        private readonly Regex mFullyTrypticAtNTerm;
        private readonly Regex mFullyTrypticAtNTermModKR;
        private readonly Regex mFullyTrypticAtCTerm;
        private readonly Regex mPartiallyTryptic1;
        private readonly Regex mPartiallyTryptic2;
        private readonly Regex mPartiallyTrypticModKR;

        public CleavageStateCalculator()
        {
            mFullyTryptic = new Regex(@"[KR]\..+[KR]\.[^P]", RegexOptions.Compiled | RegexOptions.Singleline);
            mFullyTrypticModKR = new Regex(@"[KR]\..+[KR][^A-Z]\.[^P]", RegexOptions.Compiled | RegexOptions.Singleline);
            mFullyTrypticProlineCheck = new Regex(@"\.P.*", RegexOptions.Compiled | RegexOptions.Singleline);
            mPeptideSpanningProtein = new Regex(@"-\..+\.-", RegexOptions.Compiled | RegexOptions.Singleline);
            mFullyTrypticAtNTerm = new Regex(@"-\..+[KR]\.[^P]", RegexOptions.Compiled | RegexOptions.Singleline);
            mFullyTrypticAtNTermModKR = new Regex(@"-\..+[KR][^A-Z]\.[^P]", RegexOptions.Compiled | RegexOptions.Singleline);
            mFullyTrypticAtCTerm = new Regex(@"[KR]\.[^P].+\.-", RegexOptions.Compiled | RegexOptions.Singleline);
            mPartiallyTryptic1 = new Regex(@"[KR]\.[^P].*\.\S+", RegexOptions.Compiled | RegexOptions.Singleline);
            mPartiallyTryptic2 = new Regex(@"\S+\..*[KR]\.[^P-]", RegexOptions.Compiled | RegexOptions.Singleline);
            mPartiallyTrypticModKR = new Regex(@"\S+\.[^KR].*[KR][^A-Z]\.[^P-]", RegexOptions.Compiled | RegexOptions.Singleline);
        }

        public int CountTrypticEnds(string peptideSeq)
        {
            // Fully Tryptic
            if (mFullyTryptic.IsMatch(peptideSeq) && !mFullyTrypticProlineCheck.IsMatch(peptideSeq))
            {
                return 2;
            }

            // Fully Tryptic, allowing modified K or R

            if (mFullyTrypticModKR.IsMatch(peptideSeq) && !mFullyTrypticProlineCheck.IsMatch(peptideSeq))
            {
                return 2;
            }

            // Label sequences spanning the entire protein as fully tryptic
            if (mPeptideSpanningProtein.IsMatch(peptideSeq))
            {
                return 2;
            }

            // Fully tryptic at N-Terminus
            if (mFullyTrypticAtNTerm.IsMatch(peptideSeq))
            {
                return 2;
            }

            // Fully tryptic at N-Terminus, allowing modified K or R
            if (mFullyTrypticAtNTermModKR.IsMatch(peptideSeq))
            {
                return 2;
            }

            // Fully tryptic at C-Terminus
            if (mFullyTrypticAtCTerm.IsMatch(peptideSeq))
            {
                return 2;
            }

            // Partially tryptic
            if (mPartiallyTryptic1.IsMatch(peptideSeq))
            {
                return 1;
            }

            // Partially tryptic
            if (mPartiallyTryptic2.IsMatch(peptideSeq))
            {
                return 1;
            }

            // Partially tryptic, allowing modified K or R
            if (mPartiallyTrypticModKR.IsMatch(peptideSeq))
            {
                return 1;
            }
            return 0;
        }
    }
}