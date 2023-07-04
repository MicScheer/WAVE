*CMZ :  4.00/14 30/12/2021  11.01.32  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
*KEEP,mhbook_module.
      include 'mhbook_module.cmn'
*KEND.
      program mhbook_main

      use mhbook_mod

      implicit none

      print*,"Running mhbook_main"

      call mh_path('//MHBOOK')
      call mh_fopen('mh_book.mhb')

      call mh_limit(100)
      chvar_mh(1)='var1'
      call mh_bookn(-1,'Ntitle',1,chvar_mh,2)
      var_mh(1)=1.0d0
      var_mh(2)=2.0d0
      call mh_filln(-1,var_mh)

      call mh_book1(10,"h10",1,-10.0d0,10.0d0)
      call mh_book1(1,"h1",1,-10.0d0,10.0d0)
      call mh_fill1(1,1.0d0,4.0d0)
      call mh_fill1(1,1.0d0,0.0d0)
      call mh_type(1)

      call mh_book2(2,"h2",1,-10.0d0,10.0d0,1,-3.0d0,3.0d0)
      call mh_fill2(2,1.0d0,2.0d0,10.0d0)
      call mh_info(2)

      call mh_opera(2,'/',2,22,1.0d0,1.0d0)
      call mh_info(22)

      call mh_type(22)

      call mh_delete(1)
      call mh_list1
      call mh_list2

      call mh_rout(0)
      call mh_end

      print*,"End of mhbook_main"
      end
