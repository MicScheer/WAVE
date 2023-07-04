*CMZ :  4.00/14 31/12/2021  11.02.33  by  Michael Scheer
*-- Author :    Michael Scheer   07/12/2021
      subroutine mh_kind(id,ihkind)

      use mhbook_mod
      implicit none

      integer mh_exists,id,ind,ihkind,kind(32)


      ind=mh_exists(id,ihkind)
      kind=0

      if (ihkind.eq.1) then
        kind(1)=1
        lastid_mh=id
        lastind_mh=ind
      else if (ihkind.eq.2) then
        kind(2)=1
        lastid_mh=id
        lastind_mh=ind
      else if (ihkind.eq.3) then
        kind(3)=1
        lastnid_mh=id
        lastnind_mh=ind
      endif

      return
      end
