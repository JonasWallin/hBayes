Explanation of the data follows below.   I will convert the files to 
QTL Cartographer format at some time in the future.  -Chris Basten




X-Sender: zeng@statgen.ncsu.edu
Date: Wed, 15 Sep 1999 16:24:54 -0500
To: basten
From: Zhao-Bang Zeng <zeng@brooks.statgen.ncsu.edu>
Subject: QTL data
Cc: zeng

Chris:

I include four files of the data for 

Zeng, Z.-B., J. Liu, L. F. Stam, C.-H. Kao, J. M. Mercer and C. C. Laurie.
Genetic architecture of a morphological shape difference between 
two Drosophila species. Genetics (in press)

The data set includes two backcrosses between Drosophila simulans and 
D. mauritiana. Each backcross has two independent samples. The sample
size is:
184 for bs4zb.out
287 for bs6zb.out
192 for bm4zb.out
299 for bm6zb.out

There are 45 markers with names
ewg, w, rp, v, sd, run, gl, pgk, cg, gpdh, 
ninaC, glt, prd, mhc, ddc, duc, eve, sli, plu, egfr,
twi, zip, lsp, ve, acr, dbi, h, cyc, fz, eip, 
tra, rdg, ht, ant, ninaE, fas, mst, odh, tub, hb, 
rox, ald, mlc, jan, efi                                   

The linkage map of the markers is
5 marker intervals (cM) for the first 6 markers on chromosom X:
 3.60, 10.60,  9.20, 17.20, 18.70 
15 marker intervals (cM) for the next 16 markers on chromosom 2:
 6.98, 10.10,  4.94,  6.51,  6.19, 20.46, 12.78,  3.90,  4.55,  7.49, 
30.02, 16.85,  4.34,  3.71,  7.03
22 marker intervals (cM) for the next 23 markers on chromosom 2:
 4.99,   9.34,  6.97,  7.44, 14.46,  6.79,  3.55,  6.32, 11.86,  4.58,  
 6.85,  6.35, 11.79, 12.88,  9.15,  3.30, 7.98, 13.09, 10.04,  3.70,  
 9.79,  3.43

There are 5 traits with names
pc1, adjpc1, area, areat, tibia
But the paper only reports the results on the first trait.

The code for marker data in chromosome 2 and 3 is 
2==MM (homozygotes for mauritiana alleles)
1==MS (heterozygote)
0==SS (homozygotes for simulans alleles)

However, for markers in chromosome X, the code is
2==M in bm4zb.out and bm6zb.out,  or S in bs4zb.out and bs6zb.out (hemizygote
allele)
0==null allele

Missing data are indicated by 9 for markers or 999.9999 for traits.

Zhao-Bang Zeng 

