      subroutine bezier_sing( l,N,pts,optit,otype,segnum,split,Hko,
     &                         x,y,wt,pdis,Pin,Pout,xb,yb  )
      implicit none
      integer, intent(in) :: l,N,pts,optit,otype,segnum 
      real(kind=8), intent(in) :: split 
      real(kind=8), dimension(5), intent(in) :: pdis
      real(kind=8), dimension(l),   intent(in)  :: x,y,wt
      real(kind=8), dimension(2,N),  intent(in) :: Pin
      real(kind=8), dimension(2*N,2*N), intent(in) :: Hko 
      real(kind=8), dimension(pts), intent(out) :: xb,yb
      real(kind=8), dimension(2,N), intent(out) :: Pout 
!     local variables      
      integer :: i,j,mpts,lpts,tpts 
      real(kind=8) :: pcte,pcle,lep,tep,rlep,rtep,scl,shft 
      real(kind=8), dimension(pts,N) :: ttemp,tmat
      real(kind=8), dimension(N,N)    :: a
      real(kind=8), dimension(pts,2) :: mptsu 
      real(kind=8), dimension(pts)   :: tv
      real(kind=8), dimension(2,N)    :: gradvec,Ptemp,stepd
      real(kind=8), dimension(2*N)    :: stepdir,pveck,fltgrdk
      real(kind=8), dimension(2*N,2*N) :: Hk 
      real(kind=8), dimension(:), allocatable :: tms,tx,ty,tl,tt
      real(kind=8), dimension(:,:), allocatable :: mat,mtemp,dpts
      ! Store info about spacing  
      pcle = pdis(2)
      pcte = pdis(3)
      lep  = pdis(4)
      tep  = pdis(5) 

      ! Get the polynomial basis in matrix form 
      call bernstein_matrix(N,a)
      ! This will get an array of points that represents the entire
      ! upper/lower surface as the number of points input is for the
      ! whole surface as well 
      ! tv -- entire surface 't' vector
      ! tl -- leading edge 't' vector
      ! tt -- trailing edge 't' vector 
   
      if (pdis(1) > 0.0D0) then
         call gettvec2(pts,pcle,pcte,lep,tep,tv)
      else
         call gettvec(pts,tv) 
      end if 
 
      ! Split the points on the particular surface based on the
      ! parameter 'split' 

      lpts = int( split*real(pts,kind=8) )
      tpts = pts - lpts  
      allocate( tl(lpts), tt(tpts) ) 
      tl = tv(1:lpts)
      tt = tv(lpts+1:pts)
      !
      ! Now need to rescale things a bit
      scl = 1.0D0/tl(lpts) 
      do i = 1,lpts
         tl(i) = tl(i)*scl
      end do 

      shft = tt(1)
      scl  = 1.0D0/(tt(tpts)-tt(1))
      do i = 1,tpts 
         tt(i) = (tt(i)-shft)*scl
      end do 

      ! Need to do different things based on which segment we are
      ! looking at  
      
      ! This should help with spacing
      ! *************************************************************
      mpts = 15000
      allocate( tms(mpts),tx(mpts),ty(mpts) )  
      allocate( mat(mpts,N), mtemp(mpts,N), dpts(mpts,2) ) 
      ! Get an initial equi-spaced vector with many points
      call gettvec(mpts,tms)      
      call tmatrix(tms,mpts,N,mat)
      mtemp   = matmul(mat, transpose(a)  )  
      ! Get dense collection of points from upper surface
      dpts    = matmul(mtemp,transpose(Pin))
      tx      = dpts(:,1)
      ty      = dpts(:,2)
      ! this is an operation heavy function so only call  
      if (segnum .eq. 1 .or. segnum .eq. 3) then
         call evenspace(mpts,lpts,tx,ty,tms,tl)
      else
         call evenspace(mpts,lpts,t


      Pout = Pin
      xb   = tv
      yb   = tv
      end subroutine



      subroutine bezier_opt_main( lu,ll,N,ptsu,ptsl,optit,otyp,lesc,Hk,
     &                            xu,yu,xl,yl,wtu,wtl,pdis,Pinu,Pinl,
     &                            Poutu,Poutl,xbu,ybu,xbl,ybl) 
      implicit none
      integer, intent(in) :: lu,ll,N,ptsu,ptsl,optit,otyp
      real(kind=8), intent(in) :: lesc
      real(kind=8), dimension(5),    intent(in)  :: pdis 
      real(kind=8), dimension(lu),   intent(in)  :: xu,yu,wtu
      real(kind=8), dimension(ll),   intent(in)  :: xl,yl,wtl
      real(kind=8), dimension(2,N),  intent(in)  :: Pinu,Pinl
      real(kind=8), dimension(2*N,2*N), intent(in) :: Hk 
      real(kind=8), dimension(ptsu), intent(out) :: xbu,ybu
      real(kind=8), dimension(ptsl), intent(out) :: xbl,ybl
      real(kind=8), dimension(2,N),  intent(out) :: Poutu,Poutl 
      ! variables used in subroutine
      integer :: i,j,mpts
      real(kind=8) :: pcte,pcle,lep,tep,rlep,rtep
      real(kind=8), dimension(ptsu,N) :: ttempu,tmatu
      real(kind=8), dimension(N,N)    :: a
      real(kind=8), dimension(ptsu,2) :: mptsu 
      real(kind=8), dimension(ptsu)   :: tu
      real(kind=8), dimension(ptsl,N) :: ttempl,tmatl
      real(kind=8), dimension(ptsl,2) :: mptsl 
      real(kind=8), dimension(ptsl)   :: tl
      real(kind=8), dimension(2,N)    :: gradvec,Ptemp,stepd
      real(kind=8), dimension(2*N)    :: stepdir,pveck,fltgrdk
      real(kind=8), dimension(2*N,2*N) :: Hku,Hkl  
      real(kind=8), dimension(:), allocatable :: tms,tx,ty
      real(kind=8), dimension(:,:), allocatable :: mat,mtemp,dpts
      ! Initialize both Quasi-Newton Matrices
      Hkl = Hk
      Hku = Hk 
      ! Store info about spacing  
      pcle = pdis(2)
      pcte = pdis(3)
      lep  = pdis(4)
      tep  = pdis(5) 

      ! get the polynomial basis in matrix form 
      call bernstein_matrix(N,a)
      ! 
      if (pdis(1) > 0.0D0) then
         call gettvec2(ptsu,pcle,pcte,lep,tep,tu)
         call gettvec2(ptsl,pcle,pcte,lep,tep,tl)
      else
         call gettvec(ptsu,tu)
         call gettvec(ptsl,tl)
      end if 
      ! This should help with making the leading edge spacing constant
      ! *************************************************************
      mpts = 15000
      allocate( tms(mpts),tx(mpts),ty(mpts) )  
      allocate( mat(mpts,N), mtemp(mpts,N), dpts(mpts,2) ) 
      ! Get an initial equi-spaced vector with many points
      call gettvec(mpts,tms)      
      call tmatrix(tms,mpts,N,mat)
      mtemp   = matmul(mat, transpose(a)  )  
      ! Get dense collection of points from upper surface
      dpts    = matmul(mtemp,transpose(Pinu))
      tx      = dpts(:,1)
      ty      = dpts(:,2)
      call evenspace(mpts,ptsu,tx,ty,tms,tu)

      ! Lower Surface 
      mtemp   = matmul(mat, transpose(a)  )  
      ! Get dense collection of points from upper surface
      dpts    = matmul(mtemp,transpose(Pinl))
      tx      = dpts(:,1)
      ty      = dpts(:,2)
      call evenspace(mpts,ptsl,tx,ty,tms,tl)
      deallocate( tms,tx,ty )
      deallocate( mat, mtemp, dpts )
      !***************************************************************
      
      ! generate an (npts x N) matrix [t^N,...,t,1]  
      call tmatrix(tu,ptsu,N,tmatu)
      ttempu  = matmul(tmatu, transpose(a)  )  
      mptsu   = matmul(ttempu,transpose(Pinu))            
      xbu  = mptsu(:,1)
      ybu  = mptsu(:,2)


      call tmatrix(tl,ptsl,N,tmatl)
      ttempl = matmul(tmatl, transpose(a)  )  
      mptsl   = matmul(ttempl,transpose(Pinl))            
      xbl  = mptsl(:,1)
      ybl  = mptsl(:,2)


      if (otyp == 1) then 
      Ptemp = Pinu
      do i = 1,optit
        call compgrad(lu,N,ptsu,lesc,xu,yu,wtu,ttempu,Ptemp,gradvec)
        call linesearch(lu,N,ptsu,lesc,
     &                  xu,yu,wtu,ttempu,Ptemp,gradvec,Poutu)
        Ptemp = Poutu
      end do
      mptsu   = matmul(ttempu,transpose(Ptemp))            
      xbu  = mptsu(:,1)
      ybu  = mptsu(:,2)
      else 
!     ----------------------------------------------------------      
      pveck   = 0.0D0 
      fltgrdk = 0.0D0  
      Ptemp   = Pinu
      do i = 1,optit
         call compqn(i,lu,N,ptsu,lesc,xu,yu,wtu,ttempu,Ptemp,Hku,
     &                      pveck,fltgrdk,stepdir)
         stepd(1,:)   = stepdir(1:N) 
         stepd(2,:)   = stepdir(N+1:2*N)
         call linesearch(ll,N,ptsu,lesc,
     &                  xu,yu,wtu,ttempu,Ptemp,stepd,Poutu)
         Ptemp = Poutu
      end do
       
      mptsu   = matmul(ttempu,transpose(Ptemp))            
      xbu  = mptsu(:,1)
      ybu  = mptsu(:,2)
!     -----------------------------------------------------------
      end if 

      call tmatrix(tl,ptsl,N,tmatl)
      ttempl = matmul(tmatl, transpose(a)  ) 

      if (otyp == 1) then
      Ptemp = Pinl
      do i = 1,optit
        call compgrad(ll,N,ptsl,lesc,xl,yl,wtl,ttempl,Ptemp,gradvec)
        call linesearch(ll,N,ptsl,lesc,
     &                  xl,yl,wtl,ttempl,Ptemp,gradvec,Poutl)
        Ptemp = Poutl
      end do
      mptsl   = matmul(ttempl,transpose(Ptemp))            
      xbl    = mptsl(:,1)
      ybl    = mptsl(:,2)
      else
!     ----------------------------------------------------------      
      pveck   = 0.0D0 
      fltgrdk = 0.0D0  
      Ptemp   = Pinl
      do i = 1,optit
         call compqn(i,ll,N,ptsl,lesc,xl,yl,wtl,ttempl,Ptemp,Hkl,
     &                      pveck,fltgrdk,stepdir)
         stepd(1,:)   = stepdir(1:N) 
         stepd(2,:)   = stepdir(N+1:2*N)
         call linesearch(ll,N,ptsl,lesc,
     &                  xl,yl,wtl,ttempl,Ptemp,stepd,Poutl)
         Ptemp = Poutl
      end do
       
      mptsl   = matmul(ttempl,transpose(Ptemp))            
      xbl  = mptsl(:,1)
      ybl  = mptsl(:,2)
!     -----------------------------------------------------------
      end if

      if (optit == 0) then
         Poutu = Pinu
         Poutl = Pinl
      end if

      end subroutine 


!********************************************************************
!     Subroutine that returns equi-spaced vector for points on curve
!********************************************************************
      subroutine gettvec(npts,t) 
      integer, intent(in) :: npts
      real(kind=8), dimension(npts),intent(out) :: t
      integer      :: i 
      real(kind=8) :: space
       
      space    = 1.0D0/real(npts-1)
      t(1)    = 0.0D0
      t(npts) = 1.0D0
      do i = 2,npts-1
         t(i) = space*real(i-1) 
      end do    
      end subroutine


!********************************************************************
!     Subroutine that returns vector clustered around leading/trailing
!     edge for points on curve
!********************************************************************
      subroutine gettvec2(npts,pcle,pcte,lep,tep,tu) 
      implicit none
      integer,      intent(in) :: npts
      real(kind=8), intent(in) :: pcte,pcle,lep,tep
      real(kind=8), dimension(npts),intent(out) :: tu
      integer      :: i,mpts,epts,lpts,tpts,ct 
      real(kind=8) :: space,mspace,lesp,tesp,lstr,tstr,errr,x0,temp
      real(kind=8) :: lhigh,llow,epcti,epct,mc 

      epct     = pcle + pcte 
      mpts     = int((1.0 - epct)*real(npts))
      ! Figure out how much of the chord the 'middle' section is
      ! covering 
      mc       = 1.0 - lep - tep
      mspace   = mc/real(mpts,kind=8)
      epts     = npts - mpts 
      lpts     = int( real(npts)*pcle )  
      tpts     = epts - lpts
      ! use golden section rule to find appropriate stretching ratio
      lstr  = 1.0D0
      errr  = 9999.9D0
      ct    = 1
      llow  = 0.25D0 
      lhigh = 4.0D0 
      do while (abs(errr) > 1.0e-10 .and. ct < 40)
         x0 = mspace/( (lstr)**(lpts-1) ) 
         temp = 0.0D0
         do i = 0, (lpts-2)
            temp = temp + lstr**i 
         end do 
         errr = lep - temp*x0 
         if (errr > 0.0) then
            lhigh = lstr 
         else 
            llow  = lstr
         end if
         lstr = (llow + lhigh)/2.0D0 
         ct = ct +1 
      end do
      ! Distribute the points in the leading edge 
      tu(1)    = 0.0D0 
      do i = 2,lpts-1
         tu(i) = tu(i-1) + x0*lstr**(i-2)
      end do 
      tu(lpts) = lep 
      ! Evenly distribute points along mid section 
      do i = 1,mpts - 1
         tu(lpts + i) = lep + mspace*real(i) 
      end do  
      ! Figure out trailing edge spacing 
      tstr  = 1.0D0
      errr  = 9999.9D0
      llow  = 0.25D0 
      lhigh = 4.0D0
      ct    = 1
      do while (abs(errr) > 1.0e-10 .and. ct < 40) 
         x0 = mspace/( (tstr)**(tpts-1) ) 
         temp = 0.0D0
         do i = 0,(tpts-1)
            temp = temp + tstr**i 
         end do 
         errr = tep - temp*x0 
         if (errr > 0.0) then
            lhigh = tstr 
         else 
            llow  = tstr
         end if
         tstr = (llow + lhigh)/2.0 
         ct = ct + 1
      end do
      
      tu(lpts + mpts) = 1.0D0 - tep 
      do i = 1,tpts-1 
         tu(lpts+mpts+i) = tu(lpts+mpts+i-1) + x0*((tstr)**(tpts-i))
      end do
      tu(npts) = 1.0D0 

      end subroutine


!********************************************************************
!     Trying to even out the leading edge spacing
!********************************************************************
      subroutine evenspace(dpts,npts,x,y,tlong,t)  
      implicit none
      integer :: dpts,npts 
      real(kind=8), dimension(dpts), intent(in) :: x,y,tlong
      real(kind=8), dimension(npts), intent(inout) :: t 
      integer :: i,j
      real(kind=8), dimension(npts) :: ta 
      real(kind=8) :: distot,dist,ltemp,ttemp,distm,ai,bi,ci
      logical :: lrch,trch
      
      distot = 0.0D0 
      do i = 1,dpts-1
         distot = distot + ((x(i+1)-x(i))**2 + (y(i+1)-y(i))**2)**0.5
      end do

      j     = 2 
      ta(1) = 0.0D0

      dist  = 0.0D0  
      do i = 1,dpts-1
         distm = dist 
         dist  = dist + ((x(i+1)-x(i))**2 + (y(i+1)-y(i))**2)**0.5 
         if (dist > (t(j)*distot)) then
            ai = distm
            bi = dist
            ci = t(j)*distot 
            ta(j)=tlong(i)*((bi-ci)/(bi-ai)) + 
     &            tlong(i+1)*(1.0D0-((bi-ci)/(bi-ai)))
            j = j+1 
         end if  
      end do 
   
      ta(npts) = 1.0D0 
      t = ta 

      end subroutine
         
       
           




!********************************************************************
!     Subroutine that returns matrix [t^N ... t,1] x npts
!********************************************************************
      subroutine tmatrix(tu,npts,N,tmat)
      implicit none
      integer,intent(in) :: N,npts
      real(kind=8),dimension(npts) :: tu
      real(kind=8),dimension(npts,N),intent(out) :: tmat
      
      integer :: i,j
      do j = 0,N-1
         do i = 1,npts
            tmat(i,j+1) = tu(i)**j
         end do
      end do
      end subroutine  


!********************************************************************
!     Subroutine that returns weighted 2-norm
!********************************************************************
      subroutine compnorm(l,npts,lesc,x,y,wtv,pts,norm)
      implicit none 
      integer, intent(in) :: l,npts
      real(kind=8), intent(in) :: lesc
      real(kind=8), dimension(l),      intent(in) :: x,y,wtv
      real(kind=8), dimension(npts,2), intent(in) :: pts
      real(kind=8), intent(out) :: norm
!    
      integer :: i,j 
      real(kind=8) :: dist,minv,dydx,dy,dx
      
      norm = 0.0 
      do i = 1,l
         minv = 999999.0D0
         do j = 1,npts
            dist = (x(i)-pts(j,1))**2 + (y(i)-pts(j,2))**2 
            if (dist < minv) then
               minv = dist
            end if 
         end do
         norm = norm + minv*wtv(i) 
      end do

      ! Include a component that puts a constraint on the slope at the
      ! leading edge  
      dy = pts(2,2) - pts(1,2)
      dx = pts(2,1) - pts(1,1)
      dydx = dy/dx
      norm = norm + lesc*(1/abs(dydx))      
      
      end subroutine

!********************************************************************
!     Subroutine that returns gradient vector 
!********************************************************************
      subroutine compgrad(l,N,npts,lesc,x,y,wtv,ttemp,Pin,gradvec)
      integer, intent(in) :: l,N,npts 
      real(kind=8), intent(in) :: lesc
      real(kind=8), dimension(l),      intent(in) :: x,y,wtv
      real(kind=8), dimension(2,N),    intent(in) :: Pin
      real(kind=8), dimension(npts,N), intent(in) :: ttemp
      real(kind=8), dimension(2,N),    intent(out) :: gradvec
!    
      integer :: i,j,ct  
      real(kind=8) :: diff,minv,stepsize,orignorm,norm,dy,dx
      real(kind=8) :: gradmag
      real(kind=8),dimension(2,N) :: Ptemp,gradtemp
      real(kind=8),dimension(npts,2) :: pts
      stepsize = 1e-8
  
      gradtemp = 0.0D0
      gradvec  = 0.0D0 

      pts = matmul(ttemp,transpose(Pin))   
      call compnorm(l,npts,lesc,x,y,wtv,pts,orignorm) 

      ! only need to index 2,N-1 to skip first and last point 
      do j = 2,N-1
         do i = 1,2
            Ptemp      = Pin
            Ptemp(i,j) = Pin(i,j) + Pin(i,j)*stepsize 
            pts        = matmul(ttemp,transpose(Ptemp))
            call compnorm(l,npts,lesc,x,y,wtv,pts,norm) 
            dy = norm - orignorm
            dx = Pin(i,j)*stepsize
            gradvec(i,j) = dy/dx 
         end do 
      end do

      ! normalize gradient vector 
      gradmag = 0.0 
      do j = 1,N
         do i = 1,2
            gradmag = gradmag + gradvec(i,j)**2
         end do
      end do 
      do j = 1,N
         do i = 1,2
            gradvec(i,j) = gradvec(i,j)/sqrt(gradmag) 
         end do
      end do 

      end subroutine 

!********************************************************************
!     Subroutine that returns quasi-newton step direction 
!********************************************************************
      subroutine compqn(opiter,l,N,npts,lesc,x,y,wtv,ttemp,Pin,Hk,
     &                   pveck,fltgrdk,stepdir)
      implicit none 
      integer, intent(in) :: opiter,l,N,npts 
      real(kind=8), intent(in) :: lesc
      real(kind=8), dimension(l),      intent(in) :: x,y,wtv
      real(kind=8), dimension(2,N),    intent(in) :: Pin
      real(kind=8), dimension(npts,N), intent(in) :: ttemp
      real(kind=8), dimension(2*N),    intent(out) :: stepdir
      real(kind=8), dimension(2*N,2*N),intent(inout) :: Hk
      real(kind=8), dimension(2*N),    intent(inout) :: pveck,fltgrdk
!   
      integer :: i,j,ct,ncon  
      real(kind=8) :: diff,stepsize,orignorm,dy,dx,mag,rhok
      real(kind=8), dimension(2,N) :: Ptemp,gradvec
      real(kind=8), dimension(npts,2) :: pts
      real(kind=8), dimension(2*N) :: tpmvec,pk,sk,yk,fltgrd,Pvec 
      real(kind=8), dimension(2*N,2*N) :: Hpone,tmpl,tmpr,tmpa 

      ncon     = 2*N
      tmpl     = 0.0D0 
      tmpr     = 0.0D0
      stepdir  = 0.0D0 
      do i = 1,ncon
         tmpl(i,i) = 1.0D0
         tmpr(i,i) = 1.0D0 
      end do 

      call compgrad(l,N,npts,lesc,x,y,wtv,ttemp,Pin,gradvec)
      !pts = matmul(ttemp,transpose(Pin))      
      !call compnorm(l,npts,lesc,x,y,wtv,pts,orignorm)

      ! store the gradient vector in a single dimensional array
      fltgrd(1:N)   = gradvec(1,:)    
      fltgrd(N+1:2*N) = gradvec(2,:)    
      Pvec(1:N)     = Pin(1,:) 
      Pvec(N+1:2*N)   = Pin(2,:)

      if (opiter == 1) then
         call dotmv(Hk,fltgrd,ncon,pk)
      else
         do i = 1,ncon
            sk(i) = Pvec(i)   - pveck(i) 
            yk(i) = fltgrd(i) - fltgrdk(i)
         end do
         call dotvv(sk,yk,ncon,rhok) 
         if (rhok > 0.0) then
            rhok = 1.0D0/rhok
         end if
         do i = 1,ncon
            do j = 1,ncon
               tmpl(i,j) = tmpl(i,j) - rhok*sk(i)*yk(j)
               tmpr(i,j) = tmpr(i,j) - rhok*sk(j)*yk(i) 
               tmpa(i,j) = rhok*sk(i)*sk(j) 
            end do 
         end do 
         Hpone = matmul(matmul(tmpl,Hk),tmpr) 
         do i = 1,ncon
            do j = 1,ncon
               Hpone(i,j) = Hpone(i,j) + tmpa(i,j) 
            end do 
         end do 
         pk = matmul(Hpone,fltgrd) 
         Hk = Hpone 
      end if
      stepdir = pk
!     normalize       
      mag = 0.0 
      do i = 1,ncon 
         mag = mag + stepdir(i)**2
      end do 
      do i = 1,ncon
         stepdir(i) = stepdir(i)/sqrt(mag) 
      end do
      fltgrdk = fltgrd
      pveck   = Pvec      

      end subroutine 
!********************************************************************
!     Subroutine that computes linesearch given gradient vector 
!********************************************************************
      subroutine linesearch(l,N,npts,lesc,x,y,
     &                      wtv,ttemp,Pin,gradvec,Pout)
      integer, intent(in) :: l,N,npts 
      real(kind=8), intent(in) :: lesc
      real(kind=8), dimension(l),      intent(in) :: x,y,wtv
      real(kind=8), dimension(2,N),    intent(in) :: Pin,gradvec
      real(kind=8), dimension(npts,N), intent(in) :: ttemp
      real(kind=8), dimension(2,N),    intent(out) :: Pout
!    

      integer :: i,j,ct  
      real(kind=8) :: diff,minv,steplen,stepfrac,orignorm,norm
      real(kind=8),dimension(2,N) :: Ptemp
      real(kind=8),dimension(npts,2) :: pts

      steplen = 1.0D0
      stepfrac = 0.85D0 
  
      pts = matmul(ttemp,transpose(Pin))   
      call compnorm(l,npts,lesc,x,y,wtv,pts,orignorm) 
!
      do j = 1,N
         do i = 1,2
            Ptemp(i,j) = Pin(i,j) - steplen*gradvec(i,j)  
         end do
      end do 
      pts = matmul(ttemp,transpose(Ptemp))   
      call compnorm(l,npts,lesc,x,y,wtv,pts,norm) 
      
      
      do while ( norm > orignorm )
         steplen = steplen*stepfrac      
         do j = 1,N
            do i = 1,2
               Ptemp(i,j) = Pin(i,j) - steplen*gradvec(i,j)  
            end do
         end do 
         pts = matmul(ttemp,transpose(Ptemp))   
         call compnorm(l,npts,lesc,x,y,wtv,pts,norm) 
      end do 
      Pout = Ptemp
      !write(*,*) norm    
      end subroutine 


!********************************************************************
!     Subroutine that compute dot products 
!********************************************************************
      subroutine dotvv(x,y,length,dotp)
      implicit none
      integer, intent(in) :: length
      real(kind=8), dimension(length), intent(in) :: x,y
      real(kind=8), intent(out) :: dotp
      integer :: i
      dotp = 0.0D0
      do i = 1,length
         dotp = dotp + x(i)*y(i) 
      end do 
      end subroutine


      subroutine dotmv(x,y,length,dotp)
      implicit none
      integer, intent(in) :: length
      real(kind=8), dimension(length), intent(in) :: y
      real(kind=8), dimension(length,length), intent(in) :: x
      real(kind=8), dimension(length), intent(out) :: dotp
      integer :: i,j
      dotp = 0.0D0
      do i = 1,length
         do j = 1,length
            dotp(i) = dotp(i) + x(i,j)*y(i)
         end do 
      end do 
      end subroutine

!********************************************************************
!     Subroutine that returns Bernstein polynomial matrix 
!********************************************************************
      subroutine bernstein_matrix(n,a)
      implicit none
      integer, intent(in) :: n
      real(kind=8), dimension(n,n), intent(out) :: a
      integer :: i0,j0,n0
      real(kind=8) :: r8_choose,r8_mop
       
      a  = 0.0D0
      n0 = n - 1
      do j0 = 0,n0
         do i0 = 0,j0
            a(i0+1,j0+1) = r8_mop(j0-i0)*r8_choose(n0-i0,j0-i0)* 
     &                     r8_choose(n0,i0)  
         end do
      end do 

      end subroutine 


!********************************************************************        
      function r8_choose(n,k)
      implicit none
      integer, intent(in) :: n,k
      integer :: i,mn,mx
      real(kind=8) :: value,r8_choose  
!
      mn = min(k,n-k)
      if (mn < 0) then
         value = 0.0D0 
      else if (mn == 0) then
         value = 1.0D0
      else 
         mx = max(k,n-k)
         value = real(mx+1,kind=8)
         do i = 2,mn 
            value = value*real(mx+i,kind=8)/real(i,kind=8) 
         end do 
      end if 
      r8_choose = value
      end function 
!********************************************************************       
      
      function r8_mop(i) 
      implicit none 
      integer, intent(in) :: i 
      real(kind=8) :: r8_mop 
      
      if (mod(i,2) == 0) then
         r8_mop = 1.0D0
      else
         r8_mop = -1.0D0
      end if 
      end function 












 
       




       



