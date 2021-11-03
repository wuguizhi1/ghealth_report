#16srDNA pcr primres
* forward               CCTACGGGNGGCWGCAG
* reverse               GACTACHVGGGTATCTAATCC

# second barcode pairs info;
* F1R1             ACTCCT;GACTAC
* F2R2             CGAGCC;TTGACT
* F3R3             GTGATC;CGAGAC
* F4R4             TACTGG;ACTCTG

## to ensure reads pair have complete pcr primers
* F1R1
* read1 cut first 3nt, ACT       reminds: CCTACGGGNGGCW
* read2 cut       0nt,           reminds: GACTACHVGGGTA
>-----------
* F2R2
* read1 cut first 4nt, CGAG      reminds: CCTACGGGNGGCW
* read2 cut first 2nt, TT        reminds: GACTACHVGGGTA
>-----------
* F3R3
* read1 cut first 5nt, GTGAT     reminds: CCTACGGGNGGCW
* read2 cut first 3nt, CGA       reminds: GACTACHVGGGTA
>-----------
* F4R4
* read1 cut first 6nt, TACTGG    reminds: CCTACGGGNGGCW
* read2 cut first 6nt, ACTCTG    reminds: GACTACHVGGGTA
