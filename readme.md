# Preliminary
To run the _notebook.ipynb_, please download the DS-CRISPR bank of 26-nt at: 
- https://drive.google.com/file/d/1-VhZMjlX7u2edGCqUXeRCsqK4ESOWu7c/view?usp=sharing, 
- https://drive.google.com/file/d/1L9PPqJevngIQPWQpb0kIXK8uSB0_1rnt/view?usp=sharing, 
- https://drive.google.com/file/d/1LGVQ4Gff1_kr2tiEEBKLSBdOyIMGbt-I/view?usp=sharing
- Copy the downloaded datasets into folder _crisprbank_.

Expected Output:
- Result_TargetedSequences.csv: Extracted targeting DNA Codes.
- Result_PAM_Analysis.pdf: Plot for PAM Distribution Analysis.
- Result_mRNA.csv: Optimized gRNA (small) Molecule.
- AlphaFold.csv: Meta-data file to test on AlphaFold 2 using ColabFold API.

**trainset** (DNA Codes)
- CTLA-4:     https://www.ncbi.nlm.nih.gov/nuccore/M74363.1?report=fasta
- CD80:       https://www.ncbi.nlm.nih.gov/nuccore/AH002809.2?report=fasta
- CD9:        https://www.ncbi.nlm.nih.gov/nuccore/AH006868.3?report=fasta
- CD74:       https://www.ncbi.nlm.nih.gov/nuccore/AH001484.2?report=fasta

**testset** (Protein
Uniprot
>sp|P16410|CTLA4_HUMAN Cytotoxic T-lymphocyte protein 4 OS=Homo sapiens OX=9606 GN=CTLA4 PE=1 SV=3
MACLGFQRHKAQLNLATRTWPCTLLFFLLFIPVFCKAMHVAQPAVVLASSRGIASFVCEY
ASPGKATEVRVTVLRQADSQVTEVCAATYMMGNELTFLDDSICTGTSSGNQVNLTIQGLR
AMDTGLYICKVELMYPPPYYLGIGNGTQIYVIDPEPCPDSDFLLWILAAVSSGLFFYSFL
LTAVSLSKMLKKRSPLTTGVYVKMPPTEPECEKQFQPYFIPIN


>sp|P21926|CD9_HUMAN CD9 antigen OS=Homo sapiens OX=9606 GN=CD9 PE=1 SV=4
MPVKGGTKCIKYLLFGFNFIFWLAGIAVLAIGLWLRFDSQTKSIFEQETNNNNSSFYTGV
YILIGAGALMMLVGFLGCCGAVQESQCMLGLFFGFLLVIFAIEIAAAIWGYSHKDEVIKE
VQEFYKDTYNKLKTKDEPQRETLKAIHYALNCCGLAGGVEQFISDICPKKDVLETFTVKS
CPDAIKEVFDNKFHIIGAVGIGIAVVMIFGMIFSMILCCAIRRNREMV

>sp|P33681|CD80_HUMAN T-lymphocyte activation antigen CD80 OS=Homo sapiens OX=9606 GN=CD80 PE=1 SV=1
MGHTRRQGTSPSKCPYLNFFQLLVLAGLSHFCSGVIHVTKEVKEVATLSCGHNVSVEELA
QTRIYWQKEKKMVLTMMSGDMNIWPEYKNRTIFDITNNLSIVILALRPSDEGTYECVVLK
YEKDAFKREHLAEVTLSVKADFPTPSISDFEIPTSNIRRIICSTSGGFPEPHLSWLENGE
ELNAINTTVSQDPETELYAVSSKLDFNMTTNHSFMCLIKYGHLRVNQTFNWNTTKQEHFP
DNLLPSWAITLISVNGIFVICCLTYCFAPRCRERRRNERLRRESVRPV

>sp|P04233|HG2A_HUMAN HLA class II histocompatibility antigen gamma chain OS=Homo sapiens OX=9606 GN=CD74 PE=1 SV=3
MHRRRSRSCREDQKPVMDDQRDLISNNEQLPMLGRRPGAPESKCSRGALYTGFSILVTLL
LAGQATTAYFLYQQQGRLDKLTVTSQNLQLENLRMKLPKPPKPVSKMRMATPLLMQALPM
GALPQGPMQNATKYGNMTEDHVMHLLQNADPLKVYPPLKGSFPENLRHLKNTMETIDWKV
FESWMHHWLLFEMSRHSLEQKPTDAPPKVLTKCQEEVSHIPAVHPGSFRPKCDENGNYLP
LQCYGSIGYCWCVFPNGTEVPNTRSRGHHNCSESLELEDPSSGLGVTKQDLGPVPM

