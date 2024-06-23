      dimension x(360),y(360),s(10000),u(10000),v(10000),w(10000)    
      dimension a(10000),c(10000),f(2),g(2)!,tr(6)
      integer*2 map(1000,1000,1000)!,a(10000),c(10000),f(2),g(2)!,tr(6)
      character*4 string
      real mass(10000)
      data map/1000000000*0/
      data idum/4365685/,j/0/,k/0/!,tr/0.,1.,0.,0.,0.,1./
      print *,'time steps?',' degrade spatial resolution?'
      read *,many,ten
      idum=idum+many**2
  1   x0=ran1(idum)
      y0=ran1(idum)
      z0=ran1(idum)
      k=k+1
      rmean=float(j+1)/float(many)/4.
      r=rmean!*(1.+ran1(idum))/2.
      if(x0+r.gt.1..or.x0-r.lt.0.)go to 1
      if(y0+r.gt.1..or.y0-r.lt.0.)go to 1
      if(z0+r.gt.1..or.z0-r.lt.0.)go to 1
      j=j+1!
      s(j)=r
      u(j)=x0
      v(j)=y0
      w(j)=z0
c     print *,j,x0,y0,r
c see if there's space
      it=1000./ten
      do m=1,it
      em=m
      em=em/1000.*ten
      do n=1,it
      en=n
      en=en/1000.*ten
      do nz=1,it
      ez=nz
      ez=ez/1000.*ten
      if((em-x0)**2+(en-y0)**2+(ez-z0)**2.le.r*r)then
      if(map(m,n,nz).eq.1)then
              if(mod(k,1000).eq.0)then!print *,j,k,'failed',int(p),'%'
              p=float(k)/float(many)
              print *,j,k,'failed',int(p),'%'
      end if
              if(k.gt.many*100)go to 99
              j=j-1
              go to 1
      end if
      end if
      end do
      end do
      end do
c if so bag it but shrink previous circles (comoving coords)
      do i=1,it
      do m=1,it
      do n=1,it
      map(i,m,n)=0
      end do
      end do
      end do
      iflag=0
      do m=1,it
      em=m
      em=em/1000.*ten
      do n=1,it
      en=n
      en=en/1000.*ten
      do nz=1,it
      ez=nz
      ez=ez/1000.*ten
      do l=1,j
      shrink=float(l)/float(j)
      r=s(l)*shrink
      if((em-u(l))**2+(en-v(l))**2+(ez-w(l))**2.lt.r*r)then
      map(m,n,nz)=1
      iflag=1
      end if
      end do
      end do
      end do
      end do
      if(iflag.eq.1)then
      print *,j,'mapped'
      mass(j)=alog10(s(j)*2.e-6)
      end if
      if(j.le.many)go to 1
  99  call pgbegin(0,'?',1,1)
      call pgenv(0.,1.,0.,1.,1,-1)
      do l=1,j
      shrink=float(l)/float(j)
      do i=1,360
      b=i
      b=b*.0174533
      x(i)=shrink*s(l)*cos(b)+u(l)
      y(i)=shrink*s(l)*sin(b)+v(l)
      end do
      call pgline(360,x,y)
      call pgnumb(l,0,0,string,nc)
      call pgsch(.5)
      call pgtext(u(l),v(l),string)
      end do
      sx=0
      do i=1,it
      do m=1,it
      do n=1,it
      sx=sx+map(i,m,n)
      end do
      end do
      end do
      it=it**3
      print *,'packing fraction',sx/float(it)*100.,'%'
      m=sx/float(it)*100.
      call pgnumb(m,0,0,string,nc)
      call pgsch(1.4)
      call pglabel(' ',' ','3D packing fraction')
      call pgtext(.73,1.08,string)
      call pgtext(.8,1.08,'%')
      call pgslw(3)
      call pghist(j,mass,-10.,-6.5,35,0)
      call pglabel('log M PBH M\d\(2281)',' ',' ')
      f(1)=-9
      f(2)=-9
      how=many
      g(1)=how/25.*1.6
      g(2)=how/30.*1.6
      call pgline(2,f,g)
      call pgtext(f(1)-.3,how/20.,'10.4 MeV')
      f(1)=-7
      f(2)=-7
c      g(1)=20
c      g(2)=19
      call pgline(2,f,g)
      call pgtext(f(1)-.3,how/20.,'1.2 \(0638)s')
      do i=1,j
      index=(mass(i)+12.)*10.
      a(index)=a(index)+1
      end do
      p=0
      do i=1,j
      c(i)=float(i)/10.-12.
      if(a(i).eq.0.)then
              a(i)=-0.1
      else
      a(i)=alog10(a(i))
      p=amax1(a(i),p)
      end if
      end do
      print *,p
      p=p*1.05
c     call pgenv(a(1),a(j),0.,alog10(float(j))-1.,0,30)
      call pgenv(-10.,-6.5,0.,p,0,30)
      call pglabel('log M PBH M\d\(2281)','N','log N vs log M')
      call pgline(j,c,a)
      f(1)=-9.5
      f(2)=-7.5
      g(1)=0
      g(2)=2
      call pgsls(2)
      call pgline(2,f,g)
      call pgend
      write(93,*)mass
      end

      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END  




