c------------------------------------------------------------------------c
      program isolf
c     Aaron Dotter, Steve Bjork, & Brian Chaboyer | updated 6/24/08      c
c     Contact aaron.dotter@dartmouth.edu with questions or comments      c
c     Program reads in tabulated isochrone files with any number of      c
c     ages and writes out Luminosity Functions based on user-supplied    c
c     input.                                                             c
c     To compile:      g77 -o isolf isolf.f                              c
c     The program must be run in the directory where the isochrones      c
c     are located.  The program assumes either:                          c
c        - a power-law  IMF where  dN/dM ~ M^x                           c
c        - a log-normal IMF where  dN/dM ~ exp(-log(M/Mc)^2/2sigma^2)    c
c                                                                        c
c     Output file includes a header similar to that of the input         c
c     isochrone, followed by sequential lines of LF data. Including:     c
c     1. bin number                                                      c
c     2. absolute magnitude of the bin in the chosen filter              c
c     3. log_10 cumulative number of stars, total number is nstars       c
c        (default is nstars=10^5, to change, edit nstars below)          c
c     4. log_10 number of stars in each bin                              c
c                                                                        c
c     Output from this program can be split into individual ages with    c
c     the isolf_split program.                                           c
c                                                                        c
c     This program may be run quietly from the command line (if 'IArgC'  c
c     and 'getarg' are supported as in g77) by uncommenting the lines    c
c     under 'run in quiet mode' below. To run from the command line      c
c     simply supply the arguments in the same order as requested by      c
c     the interactive mode, for example:                                 c
c             ./isolf fehm10afep0.cmd test.lf 8 4 0.1 1 -2.35            c
c     will generate a new file called 'test.lf' using empirical colors   c
c     at [Fe/H]=-1, [alpha/Fe]=0, LF in F606W, bin size 0.1 mag, and     c
c     power-law Salpeter IMF (x=-2.35). Note that the lognormal IMF      c
c     requires one additional parameter on the command line (sigma).     c
c                                                                        c
c     Note: in interactive mode it is best if both input and output      c
c     files are in the same directory as the isolf executable.           c
c------------------------------------------------------------------------c
      implicit none
      logical lquiet
      integer m0,m1,m2,mdl,num,jjj,i,j,k,jj,nn,n0,j0,n1,mgm,nages,max,
     *     iage,npts,ier,kc,junk,nstars,imf,iclr,ncol
      parameter(max=300,mgm=600,nstars=100000,ncol=10)
      integer cols(ncol), corr
      double precision dnt(mgm),dntl,ns,nsl,msave(max),tsave(max),logg,
     *     tau,teff,te(max),logl,mag(max),smt,sml(max),dmag,mb0,sigma,
     *     mb,mtot,dm,dn(mgm,3),mgrid(mgm),mx,mxmag,mlast,mnext,slast,
     *     snext,blue,v(77),v1,m,lognormal,tiny,mk1,mk2
      character infile*50,outfile*50,hdr1*60,hdr2*60,hdr3*60,
     *     sep*60,nagestr*16,arg*10,phot(ncol)*12,filter*7
      data lquiet/.false./
      data tiny/1.0d-12/
      data cols/4,4,8,8,8,5,7,7,4,77/
      data phot/'BVRI        ','uvby        ','UBVRIJHKs   ',
     *          'HST/ACS-WFC ','HST/WFPC2   ','SDSS ugriz  ',
     *          'ACS GGCT-Emp','ACS GGCT-Syn','Spitzer IRAC',
     *          'HST WFC3'/

      lquiet=.false.
c run in quiet mode - remove comments below
      if(IArgC().gt.0) then
         lquiet=.true.
         call getarg(1,infile)
         call getarg(2,outfile)
         call getarg(3,arg)
         read(arg,'(i2)') iclr
         call getarg(4,arg)
         read(arg,'(i2)') m0
         call getarg(5,arg)
         read(arg,'(f10.5)') dmag
         call getarg(6,arg)
         read(arg,'(i1)') imf
         call getarg(7,arg)
         read(arg,'(f10.5)') mx
         if(imf.eq.2) call getarg(8,arg)
         read(arg,'(f10.5)') sigma
         goto 1
      endif

c run interactively - skip this section if run from the command line
c print header info and get user-supplied input
      !Input file
      print *, '-------------------------------------------'
      print *, '     DSEP Luminosity Function Program'
      print *, '-------------------------------------------'
      print *, "  Enter isochrone filename:"
      read (*,'(a)') infile
      !Output file
      print *, "  Enter output filename:"
      read (*,'(a)') outfile
      !Choice of photometric system and filter     
      print *, '-------------------------------------------'
      print *, '  Choose photometric system: (1-8)'
      do i=1,ncol
         write(*,'(a3,i2,a2,a12)') '  [',i,'] ',phot(i)
      enddo
      print *, '  Enter choice of phot. system:'
      read *, iclr
      !fudge factor for filter names and numbers for HST/WFC3
      if(iclr.eq.10) then
         corr = 5
      else
         corr = 0
      endif
      print *, '  Choose filter to build LF from:'
      do i=1,cols(iclr)/2
         write(*,'(2(8x,i2,4x,a7))') i+corr, filter(iclr,i), 
     *      i+cols(iclr)/2+corr, filter(iclr,i+cols(iclr)/2)
      enddo
      if(mod(cols(iclr),2).eq.1)then ! if cols(iclr) is odd
         write(*,'(8x,i2,4x,a7)') cols(iclr)+5, filter(iclr,cols(iclr))
      endif
      !write(*,'(8(4x,i2,4x))') (i,i=1,cols(iclr)), 
      !write(*,'(8(1x,a7,1x))') (filter(iclr,i),i=1,cols(iclr))
      print *, "  Enter choice of filter:"
      read *, m0
      !Bin size
      print *, '-------------------------------------------'
      print *, "  Enter width of LF bins (in mag):"
      read *, dmag
      !IMF
      print *, '-------------------------------------------'
      print *, '  Enter choice of IMF:'
      print *, '  1. Power-law      2. Log-normal'
      read *, imf
      if(imf.eq.1) then
      !Power law IMF version:
         print *, '-------------------------------------------'
         print *, "  IMF power law: dN/dM ~ M^x"
         print *, "  Enter value of x:"
         read *, mx
      else if(imf.eq.2) then
      !Log-normal IMF version
         print *, '-------------------------------------------'
         print *, "  Log-normal IMF"
         print *, "  Enter value of Mc and Sigma:"
         read*, mx, sigma
      endif

 1    continue

      !redundant... but necessary to catch corr in quiet mode
      !fudge factor for filter names and numbers for HST/WFC3
      if(iclr.eq.10) then
         corr = 5
      else
         corr = 0
      endif

      !correct m0, see above
      m0 = m0 - corr
      m1 = 0
      m2 = 0

c set up color based on choice of filter
      if(iclr.eq.7.or.iclr.eq.8)then                      
c ACS GGC includes filters in 3 different photometric systems 
         if(m0.eq.1) then
            m1=1
            m2=2
         elseif(m0.eq.2.or.m0.eq.3)then
            m1=2
            m2=3
         elseif(m0.eq.4.or.m0.eq.5)then
            m1=4
            m2=5
         elseif(m0.eq.6.or.m0.eq.7)then
            m1=6
            m2=7
         endif
      else if(iclr.eq.10)then
         if(m0.eq.39)then !m0=F814W
            m1=25 !F606W
            m2=39
         else if(m0.eq.75)then !m0=F160W
            m1=65 !F110W
            m2=75
         else if(m0.lt.63)then !m0 in UVIS, use F814W
            m1=m0
            m2=39
         else  !m0 in IR, use F160W
            m1=m0
            m2=75
         endif
      else
         if(m0.lt.cols(iclr)) then
            m1=m0
            m2=m0+1
         else
            m1=m0-1
            m2=m0
         endif
      endif

      open(unit=1,file=infile,status='old',err=900)
      open(unit=2,file=outfile,status='unknown',err=901)
      goto 2
 900  if(.not.lquiet) print *,'Input file ',infile,' does not exist.'
 901  if(.not.lquiet) print *,'Output file ',outfile,' does not exist.'
      stop

c read composition, etc from input file
 2    read(1,'(a16,i2)') nagestr, nages
      read(1,'(a)') sep
      read(1,'(a)') hdr1
      read(1,'(a)') hdr2
      read(1,*)
      read(1,'(a)') hdr3
      read(1,*)
c write header to output file
      write(2,'(a16,i2)') nagestr, nages
      write(2,'(a60)') sep
      write(2,'(a60)') hdr1
      write(2,'(a60)') hdr2
      write(2,'(a60)') sep
      write(2,'(a60)') hdr3
      write(2,'(a60)') sep

c loop back here to read in each new age;
c calculate smt from nstar, for selected s, to get greater precision
      do jjj=1,nages
         read(1,'(5x,f6.3,6x,i3)') tau, npts
         read(1,*)
         iage= nint(tau*1.0d3)
         do i=1,npts
            if (i.gt.max) stop 'i exceeded dimensions'
            read(1,'(1x,i3,f10.6,99f8.4)')mdl,smt,logg,teff,logl,
     *           (v(k),k=1,cols(iclr))
            tsave(i)=v(m1)-v(m2)
            msave(i)=v(m0)
            if(i.eq.1) then
               v1=v(1)
               mb0=msave(1)
            endif
            sml(i)=log10(smt)
         enddo
         if(jjj.lt.nages) read(1,*)
         if(jjj.lt.nages) read(1,*)

c store mags and colors in arrays to be used in calculating the LF
         do j=1,npts
            mag(j)=msave(j)
            te(j)=tsave(j)
         enddo

c locate bluest point
         blue=te(1)
         do i=2,npts
            if(te(i).lt.blue) then
               blue=te(i)
               mb=mag(i)
            endif
         enddo
         m=(mb0-mb)/dmag
         mb=mb+m*dmag

c find mag(j) next fainter than mb and set mag(j0)=mb
         do j=2,npts
            if(mag(j).lt.mb) go to 150
         enddo
         if(.not.lquiet) write(*,*), 'all isochrone pts for',tau,
     *        ' are fainter than lf ', 'lower bound'
         go to 100
 150     j0=j-1
         n0=j-2
         n1=3
         nn=4
         if(n0.eq.0.or.(mag(j)-mag(j0))*(mag(j0)-mag(n0)).le.0.0d0)then
            n0=j0
            n1=2
            nn=3
         endif
         if(j+1.gt.npts.or.(mag(j+1)-mag(j))*(mag(j)-mag(n0)).le.0.0d0) 
     *        nn=nn-1
         call parrot (mag(n0),sml(n0),n1,nn,mb,sml(j0),ier)
         mag(j0)=mb
c define mgrid to be faint boundary of mag intervals
      mxmag = 1.0d2
      do i=j0,npts
        if(mag(i).lt.mxmag) mxmag=mag(i)
      enddo
      mgrid(1)=mb
      do j=2,mgm
        mgrid(j)=mgrid(1)-(j-1)*dmag
        if(mgrid(j).lt.mxmag) go to 170
      enddo
      stop 'exceeded mgrid dimension'
170   jj=j
      if(iage.eq.0) jj=jj-1
c zero out arrays which will contain sums of stars, etc in each interval
c and column
        do j=1,jj+1
          do k=1,3
            dn(j,k)=0.0d0
          enddo
        enddo
c set up parameters for initial point;
c kc is number of output columns in table; normally 1, but 2 or 3
c if mags decrease between turnoff and base of rgb
        i=j0
        j=1
        kc=1
        mlast=mag(i)
        slast=sml(i)
c follow isochrone, dividing nstar into mag strips, and columns (kc)
c check to see that points and grid intervals don't get out of phase
 200    if(mlast.lt.mgrid(j+1).or.mlast.gt.mgrid(j)) 
     1       stop 'error - mag grid out of phase'
        if(mag(i+1)-mlast) 210,220,230
c next mag is brighter than mlast; check kc and 
c move up one box in grid if necessary
 210    if(kc.eq.2.and.mag(i-1).lt.mag(i)) kc=3
        if(mlast.eq.mgrid(j+1)) j=j+1
 220    mnext=mag(i+1)
        snext=sml(i+1)
        if(mag(i+1).ge.mgrid(j+1)) go to 250
c next mag is not in same mgrid interval; interpolate on grid boundary
        mnext=mgrid(j+1)
        n0=i-1
        n1=3
        nn=4
        if(n0.eq.0.or.mag(i-1).le.mag(i)) then
           n0=i
           n1=2
           nn=3
        endif
        if(i+2.gt.npts.or.mag(i+1).le.mag(i+2)) nn=nn-1
        call parrot (mag(n0),sml(n0),n1,nn,mnext,snext,ier)
        if(ier.ne.0) then
           if(.not.lquiet) print *, '2. parrot error ', ier
           stop
        endif
        go to 250

c next mag is fainter than mlast; check kc and 
c move down one box in grid if necessary
 230    if(kc.eq.1.and.mag(i).lt.mag(i-1)) kc=2
        if(mlast.eq.mgrid(j)) j=j-1
        mnext=mag(i+1)
        snext=sml(i+1)
        if(mag(i+1).le.mgrid(j)) go to 250
c next mag is not in same mgrid interval
        mnext=mgrid(j)
        n0=i-1
        n1=3
        nn=4
        if(n0.eq.0.or.mag(i).le.mag(i-1)) then
           n0=i
           n1=2
           nn=3
        endif
        if(i+2.gt.npts.or.mag(i+2).le.mag(i+1)) nn=nn-1
        call parrot (mag(n0),sml(n0),n1,nn,mnext,snext,ier)
        if(ier.ne.0) then
           if(.not.lquiet) print *, '3. parrot error ', ier
           stop
        endif
c bc allow linear extrapolation     
c sum fraction of dn in mag strip
 250    if(snext.eq.slast) go to 280
        if(imf.eq.1.and.mx.eq.-1.0d0) then
           dm=snext-slast
        else if(imf.eq.1.and.mx.ne.-1.0d0)then
           dm=10.d0**((mx+1.0d0)*snext)-10.0d0**((mx+1.0d0)*slast) 
        else if(imf.eq.2) then
           dm=lognormal(snext,mx,sigma)-lognormal(slast,mx,sigma)
        else if(imf.eq.3) then  !Kroupa,Tout,&Gilmore(1993)IMF
           if(slast.gt.1.0d0) mk1=10.0d0**(-1.7d0*slast)
           if(slast.le.1.0d0) mk1=10.0d0**(-1.2d0*slast)
           if(slast.le.0.5d0) mk1=10.0d0**(-0.3d0*slast)
           if(snext.gt.1.0d0) mk2=10.0d0**(-1.7d0*snext)
           if(snext.le.1.0d0) mk2=10.0d0**(-1.2d0*snext)
           if(snext.le.0.5d0) mk2=10.0d0**(-0.3d0*snext)
           dm=mk2-mk1
        endif
        dn(j,kc)=dn(j,kc)+dm
c move to next segment
 280    mlast=mnext
        slast=snext
        if(mnext.eq.mag(i+1)) then
           i=i+1
           if(i.eq.npts) go to 300
        endif
        go to 200

c go back to one column if mag decrease spans less than one interval
 300    if(kc.gt.1) then
           do j=2,jj
              if(dn(j,2).gt.0.0d0.and.dn(j-1,2).eq.0.0d0
     *             .and.dn(j+1,2).eq.0.0d0) kc=1
              if(kc.eq.1) then
                 dn(j,1)=dn(j,1)+dn(j,2)+dn(j,3)
                 dn(j,2)=0.0d0
                 dn(j,3)=0.0d0
              endif
           enddo
        endif
c write individual LF header
        junk=jj-2
        if(imf.eq.1) write(2,1010) tau, junk, mx 
        if(imf.eq.2) write(2,1012) tau, junk, mx, sigma
 1010   format('#AGE=',f6.3,' BINS=',i3,' Power-law IMF: x=',f6.3)
 1012   format('#AGE=',f6.3,' BINS=',i3,' Log-normal IMF: Mc=',f5.3,
     *       ' sigma=',f5.3)
        write(2,1011) filter(iclr,m0) 
 1011   format('#',1x,'N',2x,a7,6x,'Log10(N)',7x,'Log10(dN)')

c find total dn
        ns = 0.0d0
        do j=jj-2, 1, -1
           dnt(j)=dn(j,1)+dn(j,2)+dn(j,3)
c srb 7/03 calculate also the integrated number of stars n
           ns = ns + dnt(j)
        enddo

c break into two loops to rescale the total number of stars to a constant
        mtot=dble(nstars)/ns
        ns=0.0d0
        num=0
        do j=jj-2,1,-1
           num=num+1
           dnt(j)=mtot*dnt(j) + tiny !so that the number is never zero
           ns = ns + dnt(j)          !for the next couple of lines
           if(ns.ne.0.0d0) nsl = log10(ns)
c srb 2/03 added output of log(dn), and the integrated number of stars n
           if(dnt(j).gt.0.0d0) dntl = log10(dnt(j))
c write table
           write(2,'(i3,f7.3,2f15.4)') num,mgrid(j),nsl,dntl
        enddo
        if(jjj.lt.nages) write(2,*)
        if(jjj.lt.nages) write(2,*)
c read more ages
 100    continue
      enddo
      
c close files and quit
      close(1)
      close(2)
      if(.not.lquiet) then
         print *, '-------------------------------------------'
         print *, "  LF generation completed."
         print *, "  Output written to ", outfile
         print *, '-------------------------------------------'
      endif
      
c All done!      
      end
c-----------------------------------------------------------------------
      subroutine parrot (x,y,j,n,x0,y0,ierr)
c parrot interpolates between a set of points in arrays x and y. (the x
c array must be monotonic, increasing or decreasing.)  the curve produced
c by parrot follows the data very closely: polynomial "wiggles" between 
c points are almost entirely eliminated.  the first and second derivatives 
c are everywhere continuous.
c parrot finds one rotated parabola passing through the points j-2, j-1,
c and j, with vertex exactly at j-1, and a second parabola passing through
c j-1, j, and j+1, with vertex at j, then determines the y values
c corresponding to x0 for each curve.  the returned value y0 is the
c weighted sum of the y values on the two parabolas, with the weights
c determined by the distance of the interpolated point from points j and
c j-1, respectively. (in the first and last intervals, the result is a
c weighted average of a parabola and a straight line.) 
c
c  input parameters (unchanged by parrot) :
c    x(i),y(i) - arrays of known points to be interpolated; x monotonic.
c            n - largest defined index of x(i) and y(i)
c           x0 - point for which interpolated value is to be returned
c            j - x(j) and x(j-1) must bracket x0
c
c  output parameters :
c           y0 - returned value corresponding to x0
c         ierr - error flag = 0 if no errors; = 1 if x0 is outside the
c                interval between x(1) and x(n); = 2 if x(i) are not
c                monotonic; = 3 if x0 is outside the interval between
c                x(j) and x(j-1).  if ierr = 1, y0 is a linear 
c                extrapolation.  y0 is not set for ierr > 1.
      implicit none
      integer ierr, n, i, j, l, j0, j1, j3, iline
      double precision x(n),y(n),x0,y0,x1,y1,fy,fx,x3, y3,res(2),d1,d3
      double precision r, d, ca, sa, rca, p, p2, q, q3, root,theta, cc
      double precision alpha, alt, phi, test, w, cw, sw, xp,yp,xp2,cw2
      double precision a, b, c, aa, xa, dx, dy, wt1, wt2, xo, qr, save

c check for error conditions
      ierr=0
c  ierr=1 ?
      if((x(n)-x0)*(x0-x(1)).lt.0.0d0) go to 400
c  ierr=2 ?
      do l=2,n-1
         if((x(l-1)-x(l))*(x(l)-x(l+1)).le.0.0d0) go to 410
      enddo
c  ierr=3 ?
      if((n-j)*(j-1).lt.0.) go to 420
      if((x(j)-x0)*(x0-x(j-1)).lt.0.0d0) go to 420

      do 300 i=1,2
         j0=i+j-2
         xo=x0-x(j0)
c translate coordinates to vertex of each future parabola, with x1
c always chosen to be closest to x0 (+ linear extrapolation past ends)
c (since j0 is either j or j-1, j1 is chosen to be the other)
         j1=j-i+1
         x1=x(j1)-x(j0)
         y1=y(j1)-y(j0)
         if(j0.gt.1) go to 110
         x3=x(1)-x(2)
         y3=y(1)-y(2)
         go to 125
 110     if(j0.lt.n) go to 120
         x3=x(n)-x(n-1)
         y3=y(n)-y(n-1)
         go to 125
 120     j3=j+3*i-5
         x3=x(j3)-x(j0)
         y3=y(j3)-y(j0)
 125     res(i)=0.0d0
         if(y1.eq.0.0d0.and.y3.eq.0.0d0) go to 195
c scale x coordinates to approximately same size as y coordinates
         fx=max(abs(x3),abs(x1))
         fy=max(abs(y3),abs(y1))
         if(fy.eq.0.0d0) fy=1.0d0
         x1=x1/fx
         x3=x3/fx
         xo=xo/fx
         y1=y1/fy
         y3=y3/fy
c determine angle of rotation w such that a parabola of the form y=a*x**2
c passes thru the three points (if theta is the angle from the origin
c and (x1,y1) to the new (rotated) y axis; need to solve a cubic equation
c in y = tan(theta):  y**3 + p*y**2 + 2*r*ca*y - r*sa = 0, where
c p = ca*(1-r*ca)/sa, ca = cos(alpha), sa = sin(alpha), and alpha is
c the angle between the lines from the origin to (x1,y1) and (x3,y3).
c if y=x-p/3, can get equation of form x**3 + a*x + b = 0.
         d1=x1*x1+y1*y1
         d3=x3*x3+y3*y3
         r=sqrt(d3/d1)
         d=sqrt(d1*d3)
         ca=(x1*x3+y1*y3)/d
         sa=(x1*y3-x3*y1)/d
         if(abs(sa).lt.0.001d0) go to 150
         rca=r*ca
         p=ca*(1.0d0-rca)/sa
         p2=p*p
         q=2.0d0*rca-p2/3.0d0
         r=2.0d0*(p2*p/27.0d0-p*rca/3.0d0)-r*sa
         q3=q*q*q/27.0d0
         root=0.25d0*r*r+q3
         if(root.gt.0.0d0) then ! 130,129,128
            root=sqrt(root)
            theta=atan(qr(-0.5d0*r+root)+qr(-0.5d0*r-root)-p/3.0d0)
            go to 180
         else if (root.eq.0.0d0) then
c evaluate 2 solutions; find theta (0<theta<alpha)
            cc=qr(-r/2.0d0)
            alpha=atan2(sa,ca)
            theta=atan(2.0d0*cc-p/3.0d0)
            alt=atan(-cc-p/3.0d0)
            if(alt*(alpha-alt).gt.theta*(alpha-theta)) theta=alt
            go to 180
         else
c evaluate 3 solutions; find theta (0<theta<alpha)
            phi=cos(-r/2.0d0/sqrt(-q3))/3.0d0
            cc=2.0d0*sqrt(-q/3.0d0)
            alpha=atan2(sa,ca)
            save=-1.0d0
            do 131 l=1,3
               test=atan(cc*cos(phi+(l-1)*2.0943951d0) - p/3.0d0)
               alt=test*(alpha-test)
               if(alt.lt.save) go to 131
               theta=test
               save=alt
 131        continue
            go to 180
         endif
 150     theta=sa
         if(ca.lt.0.0d0) theta=3.1415926536d0-sa
         theta=theta/2.0d0
 180     w=theta-atan2(x1,y1)
c save cos and sin of w, plus a and b
         cw=cos(w)
         sw=sin(w)
         xp=x1*cw+y1*sw
         yp=y1*cw-x1*sw
         cw2=cw*cw
         xp2=xp*xp
         c=4.0d0*sw*yp
c note c and xp2 can never be simultaneously zero, nor xp2 and yp, 
c nor yp and cw2, nor c and cw2 (unless x1=0, which is not legal 
c in parrot)
         if(abs(c).gt.1.0d6*xp2*cw2) go to 190
         a=yp/xp2
         b=c/xp2
         iline=0
         go to 200
 190     a=xp2/yp
         b=xp2/c
         iline=1
         go to 200
 195     cw=1.0d0
         sw=0.0d0
         a=0.0d0
         b=0.0d0
         iline=0
c calculate resulting y for this parabola
 200     if(x0.ne.x(j-1)) go to 201
         if(i.eq.1) go to 300
         y0=y(j-1)
         return
 201     if(x0.ne.x(j)) go to 202
         if(i.eq.1) go to 300
         y0=y(j)
         return
 202     cw2=cw*cw
         if(iline.eq.1) go to 250
         c=b*xo
         aa=a
         if(abs(c).gt.1.0d8*cw2) go to 240
 210     if(abs(c).lt.0.0001d0*cw2) go to 230
         if(abs(sw).lt.0.0001d0) go to 220
c find intersection of parabola and x0 line in rotated coordinates;
c calculate result back in unrotated coordinates.
         c=cw2-c
         if(c.lt.0.0d0) go to 420
         root=sqrt(c)
         c=xo*cw/sw
         d=2.0d0*aa*sw*sw
         res(i)=(cw-root)/d - c
         alt=(cw+root)/d - c
c choose root in interval [y(j-1),y(j)] (or nearest to that interval)
         if((y1-alt)*alt.gt.(y1-res(i))*res(i)) res(i)=alt
         go to 280
c special cases 220,230,240,250,260 below:
c rotation angle close to zero, but a not small
 220     res(i)=xo*(aa*xo/(1.0d0-c/2.0d0)+sw)/cw
         go to 280
c parabola degenerates towards a straight horizontal line
 230     c=xo/cw
         res(i)=c*(aa*c/cw+sw)
         go to 280
c xo very large
 240     xa=xo/a
         go to 260
c iline = 1
 250     xa=xo*a
         if(abs(xo).gt.1.0d6*abs(b)*cw2) go to 260
c xo is always nonzero, therefore, b, xp2, and a must be nonzero;
c return to normal case     
         c=xo/b
         aa=1.0d0/a
         go to 210
c parabola degenerates towards a folded vertical line
 260     root=sqrt(-xa/sw)
         c=xo*cw
         res(i)=(root-c)/sw
         alt=(-root-c)/sw
         if((y1-alt)*alt.gt.(y1-res(i))*res(i)) res(i)=alt
 280     res(i)=res(i)*fy+y(j0)
 300  enddo
c weight two parabolas according to distance from (j-1) and (j)
      dx=x(j)-x0
      dy=y(j)-res(2)
      wt1=sqrt(dx*dx+dy*dy)
      dx=x0-x(j-1)
      dy=res(1)-y(j-1)
      wt2=sqrt(dx*dx+dy*dy)
      y0=(res(1)*wt1+res(2)*wt2)/(wt1+wt2)
      return
 400  ierr=1
      y0=y(1)+(y(1)-y(2))*(x0-x(1))/(x(1)-x(2))
      if((x0-x(n))*(x(n)-x(n-1)).gt.0.0d0) 
     1     y0=y(n)+(y(n)-y(n-1))*(x0-x(n))/(x(n)-x(n-1))
      return
 410  ierr=2
      return
 420  ierr=3
      return
      end
cccccccccc
      function qr(x)
      implicit none
      double precision qr, x
      qr = 0.0d0
      if(x.ne.0.d0) qr = sign(exp(log(abs(x))/3.0d0),x)
      return
      end
c-----------------------------------------------------------------------
!      function lognormal(m,b,a)
!c     the integral of the log-normal IMF from M=a to M=b 
!c     dN/dM ~ exp( -(log(M)-log(Mc))^2 / (2*sigma^2) )
!      implicit none
!      double precision lognormal,m,a,b
!      lognormal=-erf((a*a-log(1.0d1)*m+log(b))/(sqrt(2.0d0)*a))
!      return
!      end
      function lognormal(log10M,mu,sigma)
      implicit none
      double precision lognormal, lnM, log10M, sigma, mu
c input m is log10(M), so convert to lnM
c mu and sigma are in linear, Msun units
      lnM=log(1d1)*log10M
      lognormal = 0.5d0*erfc( (log(mu)-lnM)/(sqrt(2d0)*log(sigma)))
      end
c-----------------------------------------------------------------------
      function filter(iphot,ifilt)
c     returns the name of the photometric system or the filter
c     iphot=the index of the photometric system
c     ifilt=the index of the filter
      implicit none
      integer iphot, ifilt, nphot, nfilt
      parameter(nphot=9,nfilt=8)
      character bands(nphot,nfilt)*5,filter*7,bandWFC3(77)*7
c          ---------------------------------------------------------
c          | BVRI | Strmgrn |UBVRIJHK|  ACS  | WFPC2 | SDSS | GC T |
c          ---------------------------------------------------------
      data bands/
     *     '  B  ','  u  ','  U  ','F435W','F336W','  u  ',2*'  B  ',
     *     '[3.6]',
     *     '  V  ','  v  ','  B  ','F475W','F439W','  g  ',2*'  V  ',
     *     '[4.5]',
     *     '  R  ','  b  ','  V  ','F555W','F450W','  r  ',2*'  I  ',
     *     '[5.8]',
     *     '  I  ','  y  ','  R  ','F606W','F555W','  i  ',2*'F606W',
     *     '[8.0]',
     *     '     ','     ','  I  ','F625W','F606W','  z  ',2*'F814W',
     *     '     ',
     *     '     ','     ','  J  ','F775W','F791W','     ',2*'F6062',
     *     '     ',
     *     '     ','     ','  H  ','F814W','F814W','     ',2*'F8142',
     *     '     ',
     *     '     ','     ','  Ks ','F850L','F850L','     ',3*'     '/
      data bandWFC3/
     *           'uvf200l','uvf218w','uvf225w','uvf275w','uvf280n',
     * 'uvf300x','uvf336w','uvf343n','uvf350l','uvf373n','uvf390m',
     * 'uvf390w','uvf395n','uvf410m','uvf438w','uvf467m','uvf469n',
     * 'uvf475w','uvf475x','uvf487n','uvf502n','uvf547m','uvf555w',
     * 'uvf600l','uvf606w','uvf621m','uvf625w','uvf631n','uvf645n',
     * 'uvf656n','uvf657n','uvf658n','uvf665n','uvf673n','uvf680n',
     * 'uvf689m','uvf763m','uvf775w','uvf814w','uvf845m','uvf850l',
     * 'uvf953n','uvfq232','uvfq243','uvfq378','uvfq387','uvfq422',
     * 'uvfq436','uvfq437','uvfq492','uvfq508','uvfq575','uvfq619',
     * 'uvfq634','uvfq672','uvfq674','uvfq727','uvfq750','uvfq889',
     * 'uvfq906','uvfq924','uvfq937','irf098m','irf105w','irf110w',
     * 'irf125w','irf126n','irf127m','irf128n','irf130n','irf132n',
     * 'irf139m','irf140w','irf153m','irf160w','irf164n','irf167n'/
      if(iphot.eq.10)then
         filter=bandWFC3(ifilt)
      else
         filter=bands(iphot,ifilt)
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
