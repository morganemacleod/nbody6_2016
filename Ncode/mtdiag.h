c
c  common block for diagnostic subroutine Nbody4/6
c  written by Michele Trenti 2005 
c  this program comes with no warranty

       integer mtNMAX,mtNMASS
       parameter(mtNMAX=230000,mtNMASS=27)

       double precision  mtx,mtvx,mtm,mttime,mtdr
       integer mtNTOT,mtNS,mtNB,mtIFIRST,mtN,mtNZERO
       double precision mtMASSdis,mtMASSvalue,mtZMBAR
       
       common /mtd/ mtZMBAR,mtNTOT,mtNS,mtNB,mtIFIRST,mtN
     &    ,mttime,mtdr(3),       
     &    mtNZERO,mtx(3,mtNMAX), mtvx(3,mtNMAX), mtm(mtNMAX)
     &    ,mtMASSdis(mtNMASS),mtMASSvalue(mtNMASS)
