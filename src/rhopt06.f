C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     RHOPT06
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     Density from (X, T6, total P) for the OPAL 2006 EOS, by lookup in
C     the precomputed inverted table built by BLDPT06.
C
C     Drop-in replacement for RHOOFP06(X,T6,P12,IRAD) in OEQOS06: same
C     arguments, same output (rho in g/cc), same failure value (-999)
C     because every doubtful case is handed to RHOOFP06 itself.
C
C     Inputs:  x    hydrogen mass fraction
C              t6   temperature / 1e6 K
C              p12  pressure in Mbar (1e12 dyn/cm2); TOTAL (gas+radiation)
C                   when irad=1, gas only when irad=0 -- as in RHOOFP06
C              irad 1 = pressure includes radiation, 0 = it does not
C
C     The table is built with irad=1 only (total pressure), so it is used
C     only when irad=1.  Falls back to the original RHOOFP06 when:
C       - irad is not 1
C       - radiation pressure > 50% of total (inversion poorly conditioned)
C       - X, log T6 or log P outside the table grid
C       - any of the 27 table nodes used is -999 (near OPAL's coverage edge)
C
C     Otherwise: 3-point (quadratic) Lagrange interpolation, first of
C     log10(rho) in log10 P, then in log10 T6, then -- in X -- of the
C     SPECIFIC VOLUME 1/rho.  At fixed P and T an ionized gas has
C     1/rho proportional to 1/mu ~ (3+5X), i.e. LINEAR in X, so quadratic
C     interpolation of 1/rho in X is exact there.  (Interpolating log rho
C     in X instead gave a systematic ~2e-4 dex error at X~0.7, which
C     shifted the solar model's L by ~0.45%.)  OPAL's ESAC06 is exact in
C     the same limit because it interpolates P/P0, also linear in 1/mu.
C
C     A running count of calls and fallbacks is written to the .short
C     file every NREP calls.
C
C     DIAGNOSTIC MODE (ldbg=.true.): every time the table is used, the
C     original RHOOFP06 is ALSO called at the same point and the two
C     log10(rho) values are compared.  The model still uses the TABLE
C     value; the comparison only observes.  Counts of differences above
C     1e-5, 1e-4, 1e-3 dex and the location of the largest difference are
C     added to the periodic report.  Set ldbg=.false. for production.
C
C--------------------------------------------------------------------------
      function rhopt06(x,t6,p12,irad)

      parameter (mx=5,nt=197)
      parameter (npt=521)
      IMPLICIT REAL*8 (A-H,O-Z)
      common/ptab06/ rlpt(npt,nt,mx),plgr(npt),tlgr(nt),xgr(mx),
     .               plmin,dlp,iptb06
      COMMON/LUOUT/ILAST,IDEBUG,ITRACK,ISHORT,IMILNE,IMODPT,ISTOR,IOWR
      dimension ft(3),fx(3)
      logical ldbg
      parameter (ldbg=.false.)
      parameter (nrep=50000)
      data sigmacc/1.8914785D-3/
      data ncall,nfall/0,0/
      data ncmp,n5,n4,n3,nrf/0,0,0,0,0/
      data dmax,xmax,tmax,pmax,fmax/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/
      save

C     ---- build the table on the first call
      if (iptb06 .ne. 12345678) call bldpt06

      ncall = ncall + 1
      if (mod(ncall,nrep) .eq. 0) then
         write(ISHORT,'(" RHOPT06: calls=",i12,"  fallbacks to",
     .        " RHOOFP06=",i12)') ncall, nfall
         if (ldbg) then
            write(ISHORT,'(" RHOPT06 DBG: compared=",i11,
     .        "  |dlogrho|>1e-5:",i10,"  >1e-4:",i10,"  >1e-3:",i10,
     .        "  ref failed:",i8)') ncmp, n5, n4, n3, nrf
            write(ISHORT,'(" RHOPT06 DBG: max |dlogrho|=",1pe10.3,
     .        " at X=",0pf7.4," logT=",f7.4," logP(cgs)=",f8.4,
     .        " Prad/P=",1pe10.3)') dmax, xmax, tmax, pmax, fmax
         endif
      endif

C     ---- cases handed to the original routine
C     table holds rho(total P) only, built with irad=1
      if (irad .ne. 1) goto 900
      if (p12 .le. 0.0D0) goto 900
C     radiation fraction, same expression as RHOOFP06 (in Mbar)
      pr = 4.0D0/3.0D0*sigmacc*t6**4
      if (pr/p12 .gt. 0.5D0) goto 900
      if (x .lt. xgr(1) .or. x .gt. xgr(mx)) goto 900
      tl = dlog10(t6)
      if (tl .gt. tlgr(1) .or. tl .lt. tlgr(nt)) goto 900
      pl = dlog10(p12)
      fp = (pl-plmin)/dlp + 1.0D0
      if (fp .lt. 1.0D0 .or. fp .gt. dble(npt)) goto 900

C     ---- pressure stencil (uniform grid): j0, j0+1, j0+2 centred on pl
      j0 = nint(fp) - 1
      j0 = max(1, min(j0, npt-2))

C     ---- temperature stencil: tlgr DEcreases with index.
C          bisection for tlgr(ilo) >= tl >= tlgr(ihi), ihi = ilo+1
      ilo = 1
      ihi = nt
   10 if (ihi-ilo .gt. 1) then
         imid = (ilo+ihi)/2
         if (tl .le. tlgr(imid)) then
            ilo = imid
         else
            ihi = imid
         endif
         goto 10
      endif
C     third point on whichever side tl is closer to
      if (tl-tlgr(ihi) .lt. tlgr(ilo)-tl) then
         i0 = ilo
      else
         i0 = ilo - 1
      endif
      i0 = max(1, min(i0, nt-2))

C     ---- composition stencil: xgr INcreases with index
      k = 1
   20 if (k .lt. mx-1 .and. x .gt. xgr(k+1)) then
         k = k + 1
         goto 20
      endif
      if (x-xgr(k) .lt. xgr(k+1)-x) then
         k0 = k - 1
      else
         k0 = k
      endif
      k0 = max(1, min(k0, mx-2))

C     ---- interpolate: log rho in log P, then in log T; 1/rho in X
      do kk=1,3
         k = k0 + kk - 1
         do ii=1,3
            i = i0 + ii - 1
            do jj=1,3
               if (rlpt(j0+jj-1,i,k) .lt. -998.0D0) goto 900
            enddo
            ft(ii) = q3pt06(pl, plgr(j0),plgr(j0+1),plgr(j0+2),
     .               rlpt(j0,i,k),rlpt(j0+1,i,k),rlpt(j0+2,i,k))
         enddo
         fx(kk) = q3pt06(tl, tlgr(i0),tlgr(i0+1),tlgr(i0+2),
     .            ft(1),ft(2),ft(3))
C        log rho -> specific volume 1/rho at this X node
         fx(kk) = 10.0D0**(-fx(kk))
      enddo
      v  = q3pt06(x, xgr(k0),xgr(k0+1),xgr(k0+2), fx(1),fx(2),fx(3))
      if (v .le. 0.0D0) goto 900
      rl = -dlog10(v)

      rhopt06 = 10.0D0**rl

C     ---- diagnostic: compare with the original inversion (observe only)
      if (ldbg) then
         rref = rhoofp06(x,t6,p12,irad)
         ncmp = ncmp + 1
         if (rref .gt. 0.0D0) then
            d = abs(rl - dlog10(rref))
            if (d .gt. 1.0D-5) n5 = n5 + 1
            if (d .gt. 1.0D-4) n4 = n4 + 1
            if (d .gt. 1.0D-3) n3 = n3 + 1
            if (d .gt. dmax) then
               dmax = d
               xmax = x
               tmax = tl + 6.0D0
               pmax = pl + 12.0D0
               fmax = pr/p12
            endif
         else
            nrf = nrf + 1
         endif
      endif
      return

C     ---- original iterative inversion
  900 continue
      nfall   = nfall + 1
      rhopt06 = rhoofp06(x,t6,p12,irad)
      return
      end

C--------------------------------------------------------------------------
C     3-point Lagrange (quadratic) interpolation through (x1,f1),(x2,f2),
C     (x3,f3), evaluated at x.  Nodes need not be equally spaced.
C--------------------------------------------------------------------------
      function q3pt06(x,x1,x2,x3,f1,f2,f3)
      IMPLICIT REAL*8 (A-H,O-Z)
      q3pt06 = f1*(x-x2)*(x-x3)/((x1-x2)*(x1-x3))
     .       + f2*(x-x1)*(x-x3)/((x2-x1)*(x2-x3))
     .       + f3*(x-x1)*(x-x2)/((x3-x1)*(x3-x2))
      return
      end
