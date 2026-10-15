C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     BLDPT06
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     Build the inverted OPAL 2006 EOS table:
C
C        rlpt(j,i,k) = log10( rho )   at   log10 P = plgr(j)  [P in Mbar,
C                                          TOTAL pressure incl. radiation]
C                                          T6       = t6a(i)  [OPAL temps]
C                                          X        = xa(k)   [OPAL X's]
C
C     Every node is computed ONCE with the existing iterative inversion
C     RHOOFP06 (irad=1), so each stored value is exactly what the current
C     code would return at that node.  Nodes that RHOOFP06 cannot invert,
C     or that invert to a density outside OPAL's tabulated range at that
C     temperature, are stored as -999 (no extrapolated values in the table).
C
C     Pressure is the FIRST index of rlpt: Fortran stores the first index
C     contiguously, and the runtime lookup steps through pressure most.
C
C     Called once, from RHOPT06, on its first call.  Sets IPTB06=12345678
C     when done (same "already done" convention as ITIME in ESAC06).
C
C     Grid: log10 P [Mbar] from PLMIN in steps of DLP, NPT points.
C           Measured OPAL coverage (check_pspacing.py, EOSOPAL06Z0.016492):
C           log10 P ~ -13.2 to 12.3 Mbar; finest native gap ~0.06 dex.
C           -13.5 + 520*0.05 = 12.5, so NPT=521 covers it with margin.
C
C--------------------------------------------------------------------------
      subroutine bldpt06

      parameter (mx=5,mv=10,nr=169,nt=197)
      parameter (npt=521)
      IMPLICIT REAL*8 (A-H,O-Z)
      common/aeos06/  xz(mx,mv,nt,nr),
     . t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),esk2(nt,nr),dfsx(mx)
     . ,dfs(nt),dfsr(nr),m,mf,xa(mx)
      common/beos06/ zz(mx),iri(10),index(10),nta(nr),nra(nt)
      common/ptab06/ rlpt(npt,nt,mx),plgr(npt),tlgr(nt),xgr(mx),
     .               plmin,dlp,iptb06
      COMMON/LUOUT/ILAST,IDEBUG,ITRACK,ISHORT,IMILNE,IMODPT,ISTOR,IOWR
      real*4 tcpu0,tcpu1
      save

C     ---- 1. make sure the OPAL tables are in memory (same dummy call
C             RHOOFP06 makes on its first use)
      xdum   = 0.5D0
      t6dum  = 1.0D0
      rdum   = 0.001D0
      idum   = 1
      irad   = 1
      CALL ESAC06 (xdum,t6dum,rdum,idum,irad,*10)
   10 continue

C     ---- 2. grids
      plmin = -13.5D0
      dlp   = 0.05D0
      do j=1,npt
         plgr(j) = plmin + dble(j-1)*dlp
      enddo
      do i=1,nt
         tlgr(i) = dlog10(t6a(i))
      enddo
      do k=1,mx
         xgr(k) = xa(k)
      enddo

C     ---- 3. one RHOOFP06 inversion per node
      nok  = 0
      nbad = 0
      call cpu_time(tcpu0)
      do k=1,mx
         do i=1,nt
            rhomin = rho(1)*(1.0D0-1.0D-10)
            rhomax = rho(nra(i))*(1.0D0+1.0D-10)
            do j=1,npt
               p12 = 10.0D0**plgr(j)
               r   = rhoofp06(xgr(k),t6a(i),p12,irad)
               if (r .ge. rhomin .and. r .le. rhomax) then
                  rlpt(j,i,k) = dlog10(r)
                  nok = nok + 1
               else
                  rlpt(j,i,k) = -999.0D0
                  nbad = nbad + 1
               endif
            enddo
         enddo
      enddo
      call cpu_time(tcpu1)

C     ---- 4. report, and mark the table as built
      iptb06 = 12345678
      write(*,100)      nok, nbad, tcpu1-tcpu0
      write(ISHORT,100) nok, nbad, tcpu1-tcpu0
  100 format(' BLDPT06: inverted OPAL06 EOS table built:',i9,
     .       ' nodes ok,',i9,' outside OPAL,',f9.2,' s CPU')

      return
      end
