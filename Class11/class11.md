Untitled
================

GitHub Documents
----------------

This is an R Markdown format used for publishing markdown documents to GitHub. When you click the **Knit** button all R code chunks are run and a markdown file (.md) suitable for publishing to GitHub is generated.

Including Code
--------------

``` r
stats<-read.csv("Data Export Summary.csv",row.names=1)
stats
```

    ##                     Proteins Nucleic.Acids Protein.NA.Complex Other  Total
    ## X-Ray                 124770          1993               6451    10 133224
    ## NMR                    10988          1273                257     8  12526
    ## Electron Microscopy     2057            31                723     0   2811
    ## Other                    250             4                  6    13    273
    ## Multi Method             127             5                  2     1    135

### Q1 determine the percentage of structures solved by X-Ray and Electron Microscopy.

``` r
pre.by.method<-stats$Total/sum(stats$Total) *100
names(pre.by.method) <-rownames(stats)
pre.by.method
```

    ##               X-Ray                 NMR Electron Microscopy 
    ##         89.43068692          8.40846082          1.88696977 
    ##               Other        Multi Method 
    ##          0.18325960          0.09062288

Also can you determine what proportion of structures are protein?
=================================================================

``` r
protein<-sum(stats$Proteins)
round((protein/sum(stats$Total)*100),1)
```

    ## [1] 92.8

``` r
library("bio3d")
```

``` r
pdb<-read.pdb("1hsg")
```

    ##   Note: Accessing on-line PDB file

``` r
pdb
```

    ## 
    ##  Call:  read.pdb(file = "1hsg")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
library(devtools)
devtools::install_bitbucket("Grantlab/bio3d-view")
```

``` r
library("bio3d.view")
view(pdb,"overview",col="sse")
```

    ## Computing connectivity from coordinates...

extract the protein only portion of this PDB structure and write it out to a new PDB file
=========================================================================================

Extract the ligand(i.e) drug
============================

``` r
ca.inds <- atom.select(pdb, "calpha")
ca.inds
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb, string = "calpha")
    ## 
    ##    Atom Indices#: 198  ($atom)
    ##    XYZ  Indices#: 594  ($xyz)
    ## 
    ## + attr: atom, xyz, call

``` r
inds.ligand<-atom.select(pdb,"ligand")
inds.ligand
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb, string = "ligand")
    ## 
    ##    Atom Indices#: 45  ($atom)
    ##    XYZ  Indices#: 135  ($xyz)
    ## 
    ## + attr: atom, xyz, call

``` r
# 45 atoms in our ligand
```

``` r
inds.ligand$atom
```

    ##  [1] 1515 1516 1517 1518 1519 1520 1521 1522 1523 1524 1525 1526 1527 1528
    ## [15] 1529 1530 1531 1532 1533 1534 1535 1536 1537 1538 1539 1540 1541 1542
    ## [29] 1543 1544 1545 1546 1547 1548 1549 1550 1551 1552 1553 1554 1555 1556
    ## [43] 1557 1558 1559

``` r
pdb$atom$resid
```

    ##    [1] "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "GLN" "GLN" "GLN" "GLN"
    ##   [12] "GLN" "GLN" "GLN" "GLN" "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##   [23] "ILE" "ILE" "THR" "THR" "THR" "THR" "THR" "THR" "THR" "LEU" "LEU"
    ##   [34] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "TRP" "TRP" "TRP" "TRP" "TRP"
    ##   [45] "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "GLN" "GLN"
    ##   [56] "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "ARG" "ARG" "ARG" "ARG"
    ##   [67] "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "PRO" "PRO" "PRO" "PRO"
    ##   [78] "PRO" "PRO" "PRO" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##   [89] "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "THR" "THR" "THR" "THR"
    ##  [100] "THR" "THR" "THR" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [111] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "ILE" "ILE"
    ##  [122] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY" "GLY"
    ##  [133] "GLY" "GLY" "GLY" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN"
    ##  [144] "GLN" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LYS" "LYS"
    ##  [155] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "GLU" "GLU" "GLU" "GLU"
    ##  [166] "GLU" "GLU" "GLU" "GLU" "GLU" "ALA" "ALA" "ALA" "ALA" "ALA" "LEU"
    ##  [177] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [188] "LEU" "LEU" "LEU" "LEU" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP"
    ##  [199] "ASP" "THR" "THR" "THR" "THR" "THR" "THR" "THR" "GLY" "GLY" "GLY"
    ##  [210] "GLY" "ALA" "ALA" "ALA" "ALA" "ALA" "ASP" "ASP" "ASP" "ASP" "ASP"
    ##  [221] "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP"
    ##  [232] "THR" "THR" "THR" "THR" "THR" "THR" "THR" "VAL" "VAL" "VAL" "VAL"
    ##  [243] "VAL" "VAL" "VAL" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [254] "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU"
    ##  [265] "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "MET" "MET" "MET" "MET"
    ##  [276] "MET" "MET" "MET" "MET" "SER" "SER" "SER" "SER" "SER" "SER" "LEU"
    ##  [287] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "PRO" "PRO" "PRO" "PRO"
    ##  [298] "PRO" "PRO" "PRO" "GLY" "GLY" "GLY" "GLY" "ARG" "ARG" "ARG" "ARG"
    ##  [309] "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "TRP" "TRP" "TRP" "TRP"
    ##  [320] "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "LYS"
    ##  [331] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "PRO" "PRO" "PRO"
    ##  [342] "PRO" "PRO" "PRO" "PRO" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS"
    ##  [353] "LYS" "LYS" "MET" "MET" "MET" "MET" "MET" "MET" "MET" "MET" "ILE"
    ##  [364] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY"
    ##  [375] "GLY" "GLY" "GLY" "GLY" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [386] "ILE" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "PHE" "PHE"
    ##  [397] "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "ILE" "ILE"
    ##  [408] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "LYS" "LYS" "LYS" "LYS" "LYS"
    ##  [419] "LYS" "LYS" "LYS" "LYS" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL"
    ##  [430] "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ##  [441] "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "TYR" "TYR"
    ##  [452] "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "ASP"
    ##  [463] "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "GLN" "GLN" "GLN" "GLN"
    ##  [474] "GLN" "GLN" "GLN" "GLN" "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [485] "ILE" "ILE" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "ILE"
    ##  [496] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLU" "GLU" "GLU" "GLU"
    ##  [507] "GLU" "GLU" "GLU" "GLU" "GLU" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [518] "ILE" "ILE" "CYS" "CYS" "CYS" "CYS" "CYS" "CYS" "GLY" "GLY" "GLY"
    ##  [529] "GLY" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS"
    ##  [540] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "ALA" "ALA"
    ##  [551] "ALA" "ALA" "ALA" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [562] "GLY" "GLY" "GLY" "GLY" "THR" "THR" "THR" "THR" "THR" "THR" "THR"
    ##  [573] "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "LEU" "LEU" "LEU" "LEU"
    ##  [584] "LEU" "LEU" "LEU" "LEU" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL"
    ##  [595] "GLY" "GLY" "GLY" "GLY" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO"
    ##  [606] "THR" "THR" "THR" "THR" "THR" "THR" "THR" "PRO" "PRO" "PRO" "PRO"
    ##  [617] "PRO" "PRO" "PRO" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "ASN"
    ##  [628] "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ILE" "ILE" "ILE" "ILE"
    ##  [639] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [650] "ILE" "GLY" "GLY" "GLY" "GLY" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ##  [661] "ARG" "ARG" "ARG" "ARG" "ARG" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN"
    ##  [672] "ASN" "ASN" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [683] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "THR" "THR" "THR" "THR"
    ##  [694] "THR" "THR" "THR" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN"
    ##  [705] "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY"
    ##  [716] "GLY" "GLY" "CYS" "CYS" "CYS" "CYS" "CYS" "CYS" "THR" "THR" "THR"
    ##  [727] "THR" "THR" "THR" "THR" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [738] "LEU" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "PHE" "PHE"
    ##  [749] "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PRO" "PRO"
    ##  [760] "PRO" "PRO" "PRO" "PRO" "PRO" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN"
    ##  [771] "GLN" "GLN" "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [782] "THR" "THR" "THR" "THR" "THR" "THR" "THR" "LEU" "LEU" "LEU" "LEU"
    ##  [793] "LEU" "LEU" "LEU" "LEU" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP"
    ##  [804] "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "GLN" "GLN" "GLN" "GLN"
    ##  [815] "GLN" "GLN" "GLN" "GLN" "GLN" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ##  [826] "ARG" "ARG" "ARG" "ARG" "ARG" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO"
    ##  [837] "PRO" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "VAL" "VAL"
    ##  [848] "VAL" "VAL" "VAL" "VAL" "VAL" "THR" "THR" "THR" "THR" "THR" "THR"
    ##  [859] "THR" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "LYS" "LYS"
    ##  [870] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "ILE" "ILE" "ILE" "ILE"
    ##  [881] "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY"
    ##  [892] "GLY" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "LEU"
    ##  [903] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LYS" "LYS" "LYS" "LYS"
    ##  [914] "LYS" "LYS" "LYS" "LYS" "LYS" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU"
    ##  [925] "GLU" "GLU" "GLU" "ALA" "ALA" "ALA" "ALA" "ALA" "LEU" "LEU" "LEU"
    ##  [936] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [947] "LEU" "LEU" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "THR"
    ##  [958] "THR" "THR" "THR" "THR" "THR" "THR" "GLY" "GLY" "GLY" "GLY" "ALA"
    ##  [969] "ALA" "ALA" "ALA" "ALA" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP"
    ##  [980] "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "THR" "THR"
    ##  [991] "THR" "THR" "THR" "THR" "THR" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL"
    ## [1002] "VAL" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "GLU" "GLU"
    ## [1013] "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU"
    ## [1024] "GLU" "GLU" "GLU" "GLU" "GLU" "MET" "MET" "MET" "MET" "MET" "MET"
    ## [1035] "MET" "MET" "SER" "SER" "SER" "SER" "SER" "SER" "LEU" "LEU" "LEU"
    ## [1046] "LEU" "LEU" "LEU" "LEU" "LEU" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO"
    ## [1057] "PRO" "GLY" "GLY" "GLY" "GLY" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ## [1068] "ARG" "ARG" "ARG" "ARG" "ARG" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP"
    ## [1079] "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "LYS" "LYS" "LYS"
    ## [1090] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "PRO" "PRO" "PRO" "PRO" "PRO"
    ## [1101] "PRO" "PRO" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS"
    ## [1112] "MET" "MET" "MET" "MET" "MET" "MET" "MET" "MET" "ILE" "ILE" "ILE"
    ## [1123] "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY"
    ## [1134] "GLY" "GLY" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY"
    ## [1145] "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "PHE" "PHE" "PHE" "PHE"
    ## [1156] "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "ILE" "ILE" "ILE" "ILE"
    ## [1167] "ILE" "ILE" "ILE" "ILE" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS"
    ## [1178] "LYS" "LYS" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "ARG" "ARG"
    ## [1189] "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "GLN" "GLN"
    ## [1200] "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "TYR" "TYR" "TYR" "TYR"
    ## [1211] "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "ASP" "ASP" "ASP"
    ## [1222] "ASP" "ASP" "ASP" "ASP" "ASP" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN"
    ## [1233] "GLN" "GLN" "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ## [1244] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "ILE" "ILE" "ILE"
    ## [1255] "ILE" "ILE" "ILE" "ILE" "ILE" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU"
    ## [1266] "GLU" "GLU" "GLU" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ## [1277] "CYS" "CYS" "CYS" "CYS" "CYS" "CYS" "GLY" "GLY" "GLY" "GLY" "HIS"
    ## [1288] "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "LYS" "LYS"
    ## [1299] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "ALA" "ALA" "ALA" "ALA"
    ## [1310] "ALA" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY"
    ## [1321] "GLY" "GLY" "THR" "THR" "THR" "THR" "THR" "THR" "THR" "VAL" "VAL"
    ## [1332] "VAL" "VAL" "VAL" "VAL" "VAL" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ## [1343] "LEU" "LEU" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "GLY" "GLY"
    ## [1354] "GLY" "GLY" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "THR" "THR"
    ## [1365] "THR" "THR" "THR" "THR" "THR" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO"
    ## [1376] "PRO" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "ASN" "ASN" "ASN"
    ## [1387] "ASN" "ASN" "ASN" "ASN" "ASN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ## [1398] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY"
    ## [1409] "GLY" "GLY" "GLY" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ## [1420] "ARG" "ARG" "ARG" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN"
    ## [1431] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ## [1442] "LEU" "LEU" "LEU" "LEU" "LEU" "THR" "THR" "THR" "THR" "THR" "THR"
    ## [1453] "THR" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "ILE"
    ## [1464] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY"
    ## [1475] "CYS" "CYS" "CYS" "CYS" "CYS" "CYS" "THR" "THR" "THR" "THR" "THR"
    ## [1486] "THR" "THR" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "ASN"
    ## [1497] "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "PHE" "PHE" "PHE" "PHE"
    ## [1508] "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "MK1" "MK1" "MK1" "MK1"
    ## [1519] "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1"
    ## [1530] "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1"
    ## [1541] "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1"
    ## [1552] "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "HOH" "HOH" "HOH"
    ## [1563] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1574] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1585] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1596] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1607] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1618] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1629] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1640] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1651] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1662] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1673] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1684] "HOH" "HOH" "HOH"

``` r
ligand.pdb<-trim.pdb(pdb,inds.ligand)
ligand.pdb
```

    ## 
    ##  Call:  trim.pdb(pdb = pdb, inds.ligand)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 45,  XYZs#: 135  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 45  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
view(ligand.pdb)
```

    ## Computing connectivity from coordinates...

``` r
write.pdb(ligand.pdb,file="ligand.pdb")
```

``` r
pdb.2<-read.pdb("1Hel")
```

    ##   Note: Accessing on-line PDB file

``` r
#normal mode analysis
modes<-nma(pdb.2)
```

    ##  Building Hessian...     Done in 0.017 seconds.
    ##  Diagonalizing Hessian...    Done in 0.094 seconds.

``` r
plot(modes)
```

![](class11_files/figure-markdown_github/unnamed-chunk-15-1.png)

``` r
m7 <- mktrj(modes,
 mode=7,
 file="mode_7.pdb")
view(m7, col=vec2color( rmsf(m7) ))
```

    ## Potential all C-alpha atom structure(s) detected: Using calpha.connectivity()
