c common data file for writing snapshots


      integer LTBM,Ltable,Ltablemax
      parameter(LTBM=20000)

      real table
      common /tbl/ Ltable,Ltablemax,table(LTBM)
