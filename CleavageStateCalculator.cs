using System.Text.RegularExpressions;

namespace SequestResultsProcessor
{
    public class CleavageStateCalculator
    {
        private readonly Regex m_reFullyTryptic;
        private readonly Regex m_reFullyTrypticProlineCheck;
        private readonly Regex m_reFullyTrypticModKR;
        private readonly Regex m_rePeptideSpanningProtein;
        private readonly Regex m_reFullyTrypticAtNTerm;
        private readonly Regex m_reFullyTrypticAtNTermModKR;
        private readonly Regex m_reFullyTrypticAtCTerm;
        private readonly Regex m_rePartiallyTryptic1;
        private readonly Regex m_rePartiallyTryptic2;
        private readonly Regex m_rePartiallyTrypticModKR;

        public CleavageStateCalculator()
        {
            m_reFullyTryptic = new Regex(@"[KR]\..+[KR]\.[^P]", RegexOptions.Compiled | RegexOptions.Singleline);
            m_reFullyTrypticModKR = new Regex(@"[KR]\..+[KR][^A-Z]\.[^P]", RegexOptions.Compiled | RegexOptions.Singleline);
            m_reFullyTrypticProlineCheck = new Regex(@"\.P.*", RegexOptions.Compiled | RegexOptions.Singleline);
            m_rePeptideSpanningProtein = new Regex(@"-\..+\.-", RegexOptions.Compiled | RegexOptions.Singleline);
            m_reFullyTrypticAtNTerm = new Regex(@"-\..+[KR]\.[^P]", RegexOptions.Compiled | RegexOptions.Singleline);
            m_reFullyTrypticAtNTermModKR = new Regex(@"-\..+[KR][^A-Z]\.[^P]", RegexOptions.Compiled | RegexOptions.Singleline);
            m_reFullyTrypticAtCTerm = new Regex(@"[KR]\.[^P].+\.-", RegexOptions.Compiled | RegexOptions.Singleline);
            m_rePartiallyTryptic1 = new Regex(@"[KR]\.[^P].*\.\S+", RegexOptions.Compiled | RegexOptions.Singleline);
            m_rePartiallyTryptic2 = new Regex(@"\S+\..*[KR]\.[^P-]", RegexOptions.Compiled | RegexOptions.Singleline);
            m_rePartiallyTrypticModKR = new Regex(@"\S+\.[^KR].*[KR][^A-Z]\.[^P-]", RegexOptions.Compiled | RegexOptions.Singleline);
        }

        public int CountTrypticEnds(string peptideSeq)
        {
            // Fully Tryptic
            if (m_reFullyTryptic.IsMatch(peptideSeq) && !m_reFullyTrypticProlineCheck.IsMatch(peptideSeq))
            {
                return 2;
            }

            // Fully Tryptic, allowing modified K or R

            if (m_reFullyTrypticModKR.IsMatch(peptideSeq) && !m_reFullyTrypticProlineCheck.IsMatch(peptideSeq))
            {
                return 2;
            }

            // Label sequences spanning the entire protein as fully tryptic
            if (m_rePeptideSpanningProtein.IsMatch(peptideSeq))
            {
                return 2;
            }

            // Fully tryptic at N-Terminus
            if (m_reFullyTrypticAtNTerm.IsMatch(peptideSeq))
            {
                return 2;
            }

            // Fully tryptic at N-Terminus, allowing modified K or R
            if (m_reFullyTrypticAtNTermModKR.IsMatch(peptideSeq))
            {
                return 2;
            }

            // Fully tryptic at C-Terminus
            if (m_reFullyTrypticAtCTerm.IsMatch(peptideSeq))
            {
                return 2;
            }

            // Partially tryptic
            if (m_rePartiallyTryptic1.IsMatch(peptideSeq))
            {
                return 1;
            }

            // Partially tryptic
            if (m_rePartiallyTryptic2.IsMatch(peptideSeq))
            {
                return 1;
            }

            // Partially tryptic, allowing modified K or R
            if (m_rePartiallyTrypticModKR.IsMatch(peptideSeq))
            {
                return 1;
            }
            return 0;
        }
    }
}