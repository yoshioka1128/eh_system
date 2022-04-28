      Program Tmatrix1d
      implicit none
! test: kmax=10, dw=25
! present: kmax=20, dw=75

    ! me and mh="electron and hole mass", mr="reduced mass"
      real(8),parameter :: me=0.0665d0,mh=0.11d0,mr=me*mh/(me+mh)
    ! mass
      real(8),parameter :: mass(2)=(/me/mr,mh/mr/)
    ! E1d="1d exciton binding energy (E1d)",a1d="1d exciton Bohr radius (a3d)"
      real(8),parameter :: E1d=4.68569374122079d0,a1d=1.0d0/dsqrt(E1d)
    ! E1d*Esc="energy unit"
      real(8),parameter :: Esc=2.0d0

    ! initial condition iread=0: SHF result, read=1: read previdous result
      integer(4),parameter :: iread=0

    ! temp="temperature", num="carrier density"
      real(8),parameter :: temp=0.5d0*E1d/Esc,beta=1.0d0/temp,num=0.1d0/a1d
    ! const="SPP parameter"
      real(8),parameter :: const=4.0d0

    ! imaginary number and circular constant
      complex(8),parameter :: im=(0.0d0,1.0d0)
      real(8),parameter :: pai=dacos(-1.0d0)

! parameter --------------------------------------------------------
    ! alpha="mix parameter"
      real(8),parameter :: alpha=1.0d0
    ! gmm0="initial broadening", cnv="order of convergence"
      real(8),parameter :: gmm0=0.1d0*E1d/Esc,cnv=1.0d-4
    ! iend="iteration limit"
      integer(4),parameter :: iend=500
    ! initial shift of chemical potential
      real(8),parameter :: dmu0=0.05d0
! ------------------------------------------------------------------


! set momentum ------------------------------------------------------
  ! integral parameter in main loop
    ! kmax="k cutoff", divk="k step size"
      real(8),parameter :: kmax=10.0d0,divk=0.2d0
    ! nk="number of k point"
      integer(4),parameter :: nk=int(kmax/divk)
    ! xk and yk="point and weight"
      real(8) :: xk(0:nk),yk(0:nk)

  ! integral parameter in Coulomb hole
    ! kCH="k cutoff"
      real(8),parameter :: kCH=100.0d0
    ! nCH="number of k point"
      integer(4),parameter :: nCH=10000
    ! xCH and yCH="point and weight"
      real(8) :: xCH(nCH),yCH(nCH)

  ! integral parameter for 0 point correction
    ! n0ptk="number of fine k points"
      integer(4),parameter :: n0ptk=500
    ! x0ptk and y0ptk="point and weight"
      real(8) :: x0ptk(n0ptk),y0ptk(n0ptk)
! ----------------------------------------------------------------------


! set omega -------------------------------------------------------
    ! divw="step size of omega", dw="extarnal part"
      real(8),parameter :: divw=0.05d0,dw=25.0d0
    ! wemax="highest energy of free electron"
      real(8), parameter :: wemax=0.5d0*kmax**2/mass(1)
    ! whmax="highest energy of free hole energy"
      real(8), parameter :: whmax=0.5d0*kmax**2/mass(2)
  ! single particle
    ! nwn1="low energy edge point of electron or hole energy"
      integer(4),parameter :: nwn1=-int(dw/divw)
    ! nwe2="high energy edge point of electron energy"
      integer(4),parameter :: nwe2=int((dw+wemax)/divw)
    ! nwh2="high energy edge point of hole energy"
      integer(4),parameter :: nwh2=int((dw+whmax)/divw)
    ! set
      integer(4),parameter :: nwg1(2)=(/nwe2,nwh2/)

  ! two particle
    ! nwnn1="low energy edge point of pair energy"
      integer(4),parameter :: nwnn1=nwn1*2
    ! nwee2, nwhh2, and nweh2 ="high energy edge point of pair energy"
      integer(4),parameter :: nwee2=nwe2*2,nwhh2=nwh2*2,nweh2=nwe2+nwh2
    ! set
!      integer(4),parameter :: nwg2(3)=(/nwee2,nweh2,nwhh2/)
      integer(4),parameter :: nwg2(2,2)=(/nwee2,nweh2,nweh2,nwhh2/)
    ! omega points
      real(8) :: xw(nwnn1:nwee2)
! ------------------------------------------------------------------


! pair green's function -----------------------------------------------------
    ! gkk="bare pair green's function"
    ! gkk(nwnn1:nwee2,0:nk,0:nk,2,2)
      complex(8),allocatable :: gkk(:,:,:,:,:)
    ! Imbgkk="bare pair green's function without bose distribution function"
    ! Imbgkk(nwnn1:nwee2,0:nk,0:nk,2,2)
      real(8),allocatable :: Imbgkk(:,:,:,:,:)
! -------------------------------------------------------------------------


! effective interaction ---------------------------------------------------
    ! Wsc and V0="screened and bare interaction"
      real(8) :: Wsc(-2*nk:2*nk),V0(-2*nk:2*nk)
    ! WscV0="Wsc*V0"
      real(8) :: WscV0(-2*nk:2*nk)
! --------------------------------------------------------------------------


! Tmatrix and Montroll Ward ------------------------------------------
    ! nknum="number of k points", IPIV="use in ZGESV"
      integer(4) :: nknum(0:2*nk),IPIV(2*nk+1)
    ! ImT="Im part of Tmatrix+MW"
    ! ImT(nwnn1:nwee2,-nk:nk,-nk:nk,2,2)
      real(8),allocatable :: ImT(:,:,:,:,:) 
    ! ImbT="Im part of Tmatrix+MW without bose distribution function"
    ! ImbT(nwnn1:nwee2,-nk:nk,-nk:nk,2,2)
      real(8),allocatable :: ImbT(:,:,:,:,:) 
    ! Wsc2nd="Wsc*V0-Wsc**2", ImMW="Montroll-Word term", ImbMW=bfncImMW
      real(8) :: Wsc2nd(-2*nk:2*nk),ImMW(2*nk+1),ImbMW(2*nk+1)
    ! matrix and T2="use in ZGESV"
      complex(8) :: matrix(2*nk+1,2*nk+1),T2(2*nk+1,2*nk+1)
    ! dmatrix and dT2="use in DGESV"
      real(8) :: dmatrix(2*nk+1,2*nk+1),dT2(2*nk+1,2*nk+1)
    ! nkbd="boundary of momentum integration"
      integer(4) :: nkbd(0:2*nk)
! ------------------------------------------------------------------

        
! self energy ---------------------------------------------------------    
    ! "exchange term or SHF self energy"
      real(8) :: self0(0:nk,2)
    ! selfT and self="correlation term and total selfenergy"
      complex(8) :: selfT(nwn1:nwe2,0:nk,2),self(nwn1:nwe2,0:nk,2)
      complex(8) :: selfold(nwn1:nwe2,0:nk,2)
! --------------------------------------------------------------------


! single particle -------------------------------------------------------
    ! Dos and Akw="density of states and spectral weight"
      real(8) :: Dos(nwn1:nwe2,2),Akw(nwn1:nwe2,0:nk,2),Akwold(nwn1:nwe2,0:nk,2)
    ! engkf, engkSHF, and engkQP="free, SHF, and QP energy"
      real(8) :: engkf(0:nk,2),engkSHF(0:nk,2),engkQP(0:nk,2)
    ! fermi and bose="fermi and bosedistribution function"
      real(8) :: fermi0(0:nk,2),fermi(nwn1:nwe2,2),bose(nwnn1:nwee2,3)
    ! "pauli blocking factor"
      real(8) :: pblock(nwn1:nwe2,nwn1:nwe2,2,2)
    ! numeh and numehold ="carrier density", numeh0="QP carrier density"
      real(8) :: numeh(2),numehold(2),numeh0(2),mu(2)
! --------------------------------------------------------------------


! optical response -----------------------------------------------------
    ! g2="q=0 pair green's function"
      complex(8) :: g2(nwnn1:nweh2,0:nk)
    ! opt and opt0="full and bare optical susceptabilitity"
      complex(8) :: opt(nwnn1:nweh2),opt0(nwnn1:nweh2)
    ! bopt="opt times bosedistribution function"
      real(8) :: bopt(nwnn1:nweh2)
    ! vertex="work array"
      complex(8) :: vertex(2*nk+1)
! -----------------------------------------------------------------------


! Principal Integration ---------------------------------------------------
    ! log part of principal integration
      real(8) :: logrm0(nwnn1:nwee2,nwnn1:nwee2)
    ! Sw and fnc="input and output function"
      real(8) :: Sw(nwnn1:nwee2),fnc(nwnn1:nwee2)
! ---------------------------------------------------------------------


! Coulomb Potential ----------------------------------------------
    ! widx and widy="quantum well width of x and y direction (a3d)"
      real(8),parameter :: aspect=7.0d0/3.0d0
      real(8),parameter :: widx=3.0d0*0.154d0,widy=aspect*widx
    ! nxmax and nymax="number of real x and y"
      integer(4),parameter :: nxmax=200,nymax=int(aspect*dble(nxmax))
    ! divx and divy="real x and y step size"
      real(8),parameter :: divx=widx/dble(nxmax),divy=widy/dble(nymax)
    ! phi2x and phi2y="sinx**2 and siny**2 function"
      real(8) :: phi2x(0:nxmax),phi2y(0:nymax)
    ! Vtot,V0pt,VCH="Coulomb potential for main loop,0point correction,CH"
      real(8) :: Vtot(2*nk),V0pt(n0ptk),VCH(nCH)
  ! integral parameter for 0 point correction
    ! n0ptx="number of fine x points"
      integer(4),parameter :: n0ptx=100
    ! x0ptx and y0ptx="point and weight"
      real(8) :: x0ptx(n0ptx),y0ptx(n0ptx)
    ! Work alignment
      real(8) :: midint(0:nymax),bessk(0:nxmax,0:nymax)
! -------------------------------------------------------------------


! exact diagonalization -----------------------------------------
      character(1) JOBZ,UPLO
      data JOBZ/'V'/,UPLO/'U'/
      integer(4),parameter :: LWORK=3*(2*nk+1)
      real(8) :: eigen(2*nk+1),H0(2*nk+1,2*nk+1),WORK(LWORK)
! -------------------------------------------------------------------


      integer(4) :: &
        i,it,j,iw,iw1,iw2,iw3,iw4,i1,i2,i3,j1,j2,j3,l,ieh,ieh1,ieh2,&
        ik1,ik2,ik3,ik,ip1,ip2,ip,iq,iq0,nwmax,icent,icent0,ite,ite2,&
        icount,iave,iabc,ite3,INFO,ix,iy
      real(8) :: &
        k1,k2,absorb,selfCH,ktf2,kdh2,kappa,ffnc,w,ww,www,w0,sum,bfnc,&
        dense,p,q,t,q0,CH,CH2,abc,def,ghi,jkl,S1,Img2,a,a2,&
        screening,cffInt,xx,yy,Coulomb
        







      write(6,*) 'omega points',nwee2-nwnn1+1
      write(6,*) 'ImT size',&
           dble(nwee2-nwnn1+1)*dble(2*nk+1)*dble(2*nk+1)*4.0d0*8.0d-9,'G byte'

      write(6,*) 'T (E1d)',Esc/beta/E1d
      write(6,*) 'carrier density (1/a1d)',num*a1d
      write(6,*)
      write(6,*) 'setup'
      write(6,*) 'kmax=',kmax*a1d
      write(6,*) 'nk=',nk
      write(6,*) 'dw=',dw*Esc/E1d
      write(6,*) 'divw=',divw*Esc/E1d
      write(6,*) 'alpha=',alpha
      write(6,*) 'iread=',iread
      write(6,*)


    
! set momentum ----------------------------------------
  ! integral parameter in main loop
    ! xk="k points"
      do ik=0,nk ! k
         xk(ik)=divk*dble(ik)
      end do
    ! yk="weight"
      yk=divk
      yk(0)=0.5d0*divk
      
    ! integral parameter in Coulomb hole
      call gauleg(0.0d0,kCH,xCH,yCH,nCH)
    ! integral parameter in 0pt correction
      call gauleg(0.0d0,divk,x0ptk,y0ptk,n0ptk)
! -----------------------------------------------------------

      
! set omega -----------------------------------------------
    ! xw="omega points"
      do iw=nwnn1,nwee2
         xw(iw)=divw*dble(iw)
      end do
! -----------------------------------------------------------


! set real space integral parameter -----------------------
      call gauleg(0.0d0,divx,x0ptx,y0ptx,n0ptx)
! -----------------------------------------------------------


! principal integration -----------------------------------------
      logrm0=0.0d0
      do iw=nwnn1,nwee2
        w=xw(iw)
        do iw2=nwnn1,nwee2-1
          if(iw2.eq.iw-1.or.iw2.eq.iw) then
            abc=1.0d0
          else
            abc=(xw(iw2+1)-w)/(xw(iw2)-w)
          end if
          logrm0(iw2,iw)=dlog(abc)
        end do
      end do
! -----------------------------------------------------------


! set boundary of k at q -------------------------------------
    ! boundary
      do iq=0,2*nk
        nkbd(iq)=max0(-nk,iq-nk)
      end do
    ! number of k points at q
      do iq=0,2*nk
        nknum(iq)=nk-nkbd(iq)+1
      end do
! -----------------------------------------------------------      


! free particle energy ---------------------------------------------
      do ieh=1,2
        do ik1=0,nk ! k1
          engkf(ik1,ieh)=0.5d0*xk(ik1)**2/mass(ieh)
        end do
      end do
! --------------------------------------------------------------


! bare Coulomb potential ----------------------------------------
    ! set wave function
      do ix=0,nxmax,1
        xx=divx*dble(ix)
        phi2x(ix)=2.0d0/widx*(dsin(pai*xx/widx))**2
      end do
      do iy=0,nymax,1
        yy=divy*dble(iy)
        phi2y(iy)=2.0d0/widy*(dsin(pai*yy/widy))**2
      end do

    ! for main loop
!$omp parallel do private(iq,q,bessk,midint)
      do iq=1,2*nk
        q=divk*dble(iq)
        Vtot(iq)=Coulomb(q,x0ptx,y0ptx,n0ptx,divx,divy,&
          bessk,nxmax,nymax,phi2x,phi2y,midint)
      end do
!$omp end parallel do
    ! for 0point correction
!$omp parallel do private(iq,q,bessk,midint)
      do iq=1,n0ptk
        q=x0ptk(iq)
        V0pt(iq)=Coulomb(q,x0ptx,y0ptx,n0ptx,divx,divy,&
          bessk,nxmax,nymax,phi2x,phi2y,midint)
      end do
!$omp end parallel do
    ! for Coulomb hole
!$omp parallel do private(iq,q,bessk,midint)
      do iq=1,nCH
        q=xCH(iq)
        VCH(iq)=Coulomb(q,x0ptx,y0ptx,n0ptx,divx,divy,&
          bessk,nxmax,nymax,phi2x,phi2y,midint)
      end do
!$omp end parallel do
! -----------------------------------------------------------      


! set bare interaction ---------------------------------------
      kappa=0.0d0
      cffInt=0.0d0
      call Intexp(0,nk,V0,divk,kappa,&
      n0ptk,x0ptk,y0ptk,Vtot,V0pt,cffInt)
! -------------------------------------------------------


! exact diagonalization -------------------------------------------
    ! kinetic term
      H0=0.0d0
      do ik1=-nk,nk
        H0(ik1+nk+1,ik1+nk+1)=0.5d0*xk(abs(ik1))**2
      end do
    ! Interaction term
      do ik2=-nk,nk
        do ik1=-nk,nk
          H0(ik1+nk+1,ik2+nk+1)=H0(ik1+nk+1,ik2+nk+1)&
          -V0(ik1-ik2)/(2.0d0*pai)*divk
        end do
      end do

      call DSYEV(JOBZ,UPLO,2*nk+1,H0,2*nk+1,eigen,WORK,LWORK,INFO)
      write(6,*) 'binding energy (E1d)',-eigen(1)*Esc/E1d
      write(6,*)
! ----------------------------------------------------------------


! set initial condition -----------------------------------------------      
      if(iread.eq.1) then
      ! previous mu and selfenergy read write(89)
        read(5,*) mu(1),mu(2)
        write(6,*) 'previous mu Ry'
        write(6,*) Esc*mu(1)/E1d,Esc*mu(2)/E1d

        do ieh1=1,2
          do ik1=0,nk
            do iw=nwn1,nwg1(ieh1)
              read(5,*) self(iw,ik1,ieh1)
            end do
          end do
        end do
        
      else
      ! SHF mu and selfenergy
        do i=1,2 ! initial mu (classical)
          mu(i)=dlog(num*dsqrt(0.5d0*pai*beta/mass(i)))/beta
        end do

        call SHF(mu,beta,num,nk,engkSHF,engkf,fermi0,dmu0,&
          yk,Wsc,Vtot,V0pt,VCH,n0ptk,x0ptk,y0ptk,divk,xk,&
          nCH,xCH,yCH,self0,cnv,const,iend,E1d,a1d)
        do ieh1=1,2
          do ik1=0,nk
            do iw=nwn1,nwg1(ieh1)
              self(iw,ik1,ieh1)=self0(ik1,ieh1)-im*gmm0
            end do
          end do
        end do

      end if
! ----------------------------------------------------------------------


! set Akw Dos and calculate carrier density(numeh) --------------------      
      call Dosnumber(Akw,engkf,engkQP,nk,nwg1,xw,divw,&
      nwn1,nwe2,nwnn1,nwee2,self,Dos,yk,numeh,numeh0,&
      mu,beta,E1d,a1d)
      write(6,*) 'e-h particle density (1/a1d)'
      write(6,*) 0,numeh(1)*a1d,numeh(2)*a1d
! ----------------------------------------------------------------------      

      
! resutlts in progress -------------------------------------------------      
      open(29,file='T0.5n0.1_E_num_ite.txt',status='unknown')
      open(30,file='T0.5n0.1_H_num_ite.txt',status='unknown')
      open(129,file='T0.5n0.1_kappa_ite.txt',status='unknown')
      open(130,file='T0.5n0.1_mu_ite.txt',status='unknown')
      open(133,file='T0.5n0.1_alpha_ite.txt',status='unknown')
      write(29,*) 0,numeh(1)*a1d,num
      write(30,*) 0,numeh(2)*a1d,num*a1d
      write(130,*) 0,Esc*(mu(1)+mu(2))/E1d
! ----------------------------------------------------------------------      


! set initial condition ---------------------------------
      selfold=self
      numehold=numeh
      Akwold=Akw
! -------------------------------------------------------

    ! main loop 700
      icount=0
      do 700
    ! icount="number of iteration"
      icount=icount+1
      
! Mix selfenergy by the rate, alpha --------------------------
      do ieh1=1,2
        do ik1=0,nk
          do iw=nwn1,nwg1(ieh1)
            self(iw,ik1,ieh1)=alpha*self(iw,ik1,ieh1)+&
            (1.0d0-alpha)*selfold(iw,ik1,ieh1)
          end do
        end do 
      end do  
      selfold=self
! --------------------------------------------------------------


! set chemical potential ----------------------------------------
    ! set Akw Dos and calculate carrier density(numeh)
      call Dosnumber(Akw,engkf,engkQP,nk,nwg1,xw,divw,&
      nwn1,nwe2,nwnn1,nwee2,self,Dos,yk,numeh,numeh0,&
      mu,beta,E1d,a1d)

    ! adjust mu (numeh=num)
      do ieh=1,2
        call shiftmu(num,ieh,mu,dmu0,numeh,&
        beta,Dos,xw,divw,nwn1,nwe2,nwnn1,nwee2,nwg1(2))
        do iw=nwn1,nwe2
          fermi(iw,ieh)=ffnc(beta,xw(iw),mu(ieh)) ! fermi distribution
        end do
      end do
    ! pauli blocking factor
      do ieh2=1,2
        do ieh1=1,2
          do iw2=nwn1,nwg1(ieh2)
            do iw1=nwn1,nwg1(ieh1)
               pblock(iw1,iw2,ieh1,ieh2)=1.0d0-fermi(iw1,ieh1)-fermi(iw2,ieh2)
            end do
          end do
        end do
      end do
    ! QP fermi distribution
      do ieh=1,2
        do ik=0,nk ! k
          fermi0(ik,ieh)=ffnc(beta,engkQP(ik,ieh),mu(ieh))
        end do
      end do

    ! bose distribution (0point correction)
!      call bosedis(bose,mu,nwnn1,nwee2,xw,beta,divw)
! ----------------------------------------------------------------


! update interaction -------------------------------------------
    ! screening parameter
      kappa=screening(mu,beta,fermi0,nk,xk,yk,kappa)
      cffInt=const*kappa/16.0d0/num
    ! effective interaction
      call Intexp(1,nk,Wsc,divk,kappa,&
        n0ptk,x0ptk,y0ptk,Vtot,V0pt,cffInt)
    ! set WscV0exp="Wsc*V0" (for Montroll Word term)
      call Int2nd(1,nk,WscV0,Wsc,divk,kappa,&
        n0ptk,x0ptk,y0ptk,Vtot,V0pt,cffInt)
! ----------------------------------------------------------------

      
! HF selfenergy (bare interaction) ------------------------
      self0=0.0d0
      do ieh=1,2
        do ik2=0,nk ! k2
          abc=0.0d0
          do iw=nwn1,nwg1(ieh)
            abc=abc+fermi(iw,ieh)*Akw(iw,ik2,ieh)
          end do
          abc=yk(ik2)*divw*abc/(2.0d0*pai)
          do ik1=0,nk ! k1
            self0(ik1,ieh)=self0(ik1,ieh)+&
            (V0(ik1-ik2)+V0(ik1+ik2))*abc
          end do
        end do
      end do
      self0=-self0/(2.0d0*pai)
! ---------------------------------------------------------------------------


! set pair Green's fnc. ----------------------------------------------      
      allocate (gkk(nwnn1:nwee2,0:nk,0:nk,2,2))
      allocate (Imbgkk(nwnn1:nwee2,0:nk,0:nk,2,2))
    ! set Im part of gkk
      call Impair(gkk,Imbgkk,nweh2,nk,nwn1,nwe2,nwnn1,nwee2,&
        xw,divw,nwg1,nwg2,Akw,pblock,fermi,xk,E1d,a1d)

    ! Kramers-Kronig transformation
      do ieh2=1,2
        do ieh1=1,2
          do ik2=0,nk
            do ik1=0,ik2
              do iw=nwnn1,nwg2(ieh1,ieh2)
                fnc(iw)=dimag(gkk(iw,ik1,ik2,ieh1,ieh2))
              end do
              do iw2=nwnn1,nwg2(ieh1,ieh2)-1
                abc=(fnc(iw2+1)-fnc(iw2))/(xw(iw2+1)-xw(iw2))
!$omp parallel do private(iw,w)
                do iw=nwnn1,nwg2(ieh1,ieh2)
                  w=xw(iw)
                  gkk(iw,ik1,ik2,ieh1,ieh2)=gkk(iw,ik1,ik2,ieh1,ieh2)+&
                  logrm0(iw2,iw)*((w-xw(iw2))*abc+fnc(iw2))/pai
                end do
!$omp end parallel do
              end do
            end do
          end do
        end do
      end do

    ! fill lower triangle of gkk and Imbgkk
      do ieh1=1,2
        do ieh2=ieh1,2
          do ik2=0,nk
            do ik1=0,ik2
              do iw=nwnn1,nwg2(ieh1,ieh2)
                gkk(iw,ik2,ik1,ieh2,ieh1)=gkk(iw,ik1,ik2,ieh1,ieh2)
                Imbgkk(iw,ik2,ik1,ieh2,ieh1)=Imbgkk(iw,ik1,ik2,ieh1,ieh2)
              end do
            end do
          end do
        end do
      end do
! --------------------------------------------------------------


! set Tmatrix and Montroll Word -----------------------------------------
    ! set WscV0exp="Wsc*V0-Wsc**2"
      do ik1=-2*nk,2*nk
        Wsc2nd(ik1)=WscV0(ik1)-Wsc(ik1)**2
      end do

      allocate (ImT(nwnn1:nwee2,-nk:nk,-nk:nk,2,2))
      allocate (ImbT(nwnn1:nwee2,-nk:nk,-nk:nk,2,2))
    ! imaginary part of correlation term
      call TmatrixMW(ImT,ImbT,gkk,Imbgkk,nwnn1,nk,divk,nkbd,nknum,&
      Wsc,Wsc2nd,IPIV,ImMW,ImbMW,dmatrix,matrix,dT2,T2,&
      nwee2,nwg2,xk,xw,E1d,a1d)
    ! set g0(k,-k,q=0,w)
      do iw=nwnn1,nweh2
        do ik1=0,nk
          g2(iw,ik1)=-gkk(iw,ik1,ik1,2,1)
        end do
      end do
      deallocate(gkk,Imbgkk)

      do iq=0,2*nk
         do ik1=nkbd(iq),nk
            ik2=iq-ik1
            do iw=nwnn1,nweh2
             ! parity symmetry
               ImT(iw,-ik1,-ik2,2,1)=ImT(iw,ik1,ik2,2,1)
               ImT(iw,-ik1,-ik2,1,1)=ImT(iw,ik1,ik2,1,1)
               ImT(iw,-ik1,-ik2,2,2)=ImT(iw,ik1,ik2,2,2)
             ! exchange electron and hole
               ImT(iw,ik2,ik1,1,2)=ImT(iw,ik1,ik2,2,1)
               ImT(iw,-ik2,-ik1,1,2)=ImT(iw,ik1,ik2,2,1)

             ! parity symmetry
               ImbT(iw,-ik1,-ik2,2,1)=ImbT(iw,ik1,ik2,2,1)
               ImbT(iw,-ik1,-ik2,1,1)=ImbT(iw,ik1,ik2,1,1)
               ImbT(iw,-ik1,-ik2,2,2)=ImbT(iw,ik1,ik2,2,2)
             ! exchange electron and hole
               ImbT(iw,ik2,ik1,1,2)=ImbT(iw,ik1,ik2,2,1)
               ImbT(iw,-ik2,-ik1,1,2)=ImbT(iw,ik1,ik2,2,1)
            end do
         end do
      end do
! ------------------------------------------------------------


! Tmatrix selfenergy ------------------------------------------------------
    ! imaginary part of self energy
      call ImselfTmatrixMW(selfT,ImT,ImbT,nwn1,nwe2,&
      nwnn1,nwee2,nk,nwg1,xw,divw,divk,fermi,bose,Akw)
      deallocate(ImT,ImbT)
    ! check sign of Im part selfenergy
      do ieh=1,2
        do ik1=0,nk ! k1
          do iw=nwn1,nwg1(ieh)
            if(dimag(selfT(iw,ik1,ieh)).gt.0.0d0) then
              write(6,*) 'im_self is positive'
              selfT(iw,ik1,ieh)=0.0d0
            end if
          end do
        end do
      end do

    ! Kramers-Kronig transformation
      do ieh1=1,2
        do ik1=0,nk
          do iw=nwn1,nwg1(ieh1)
            fnc(iw)=dimag(selfT(iw,ik1,ieh1))
          end do
          do iw2=nwn1,nwg1(ieh1)-1
            abc=(fnc(iw2+1)-fnc(iw2))/(xw(iw2+1)-xw(iw2))
!$omp parallel do private(iw,w)
            do iw=nwn1,nwg1(ieh1)
              w=xw(iw)
              selfT(iw,ik1,ieh1)=selfT(iw,ik1,ieh1)+&
              logrm0(iw2,iw)*((w-xw(iw2))*abc+fnc(iw2))/pai
            end do
!$omp end parallel do
          end do
        end do
      end do 
! -----------------------------------------------------------------


! total selfenergy ------------------------------------------------      
      do ieh1=1,2
        do ik1=0,nk
          do iw=nwn1,nwg1(ieh1)
            self(iw,ik1,ieh1)=selfT(iw,ik1,ieh1)+self0(ik1,ieh1)
          end do
        end do
      end do
! --------------------------------------------------------------


! total carrier density and set Akw ----------------------------------      
      call Dosnumber(Akw,engkf,engkQP,nk,nwg1,xw,divw,&
      nwn1,nwe2,nwnn1,nwee2,self,Dos,yk,numeh,numeh0,&
      mu,beta,E1d,a1d)
! --------------------------------------------------------------------
      

! results in progress -------------------------------------------------
      call output(nk,nwg1,xk,xw,nwn1,nwe2,nwnn1,nwee2,&
        engkQP,kappa,Akw,numeh,numeh0,mu,self,E1d,a1d)
      write(6,*) icount,numeh(1)*a1d,numeh(2)*a1d
      write(29,*) icount,numeh(1)*a1d,num*a1d
      write(30,*) icount,numeh(2)*a1d,num*a1d
      write(129,*) icount,kappa
      write(130,*) icount,Esc*(mu(1)+mu(2))/E1d        
      write(133,*) icount,(numeh0(1)/numeh(1)+numeh0(2)/numeh(2))/2.0d0
! --------------------------------------------------------------------------


! optical response ----------------------------------------------------
      call optical(nk,nwnn1,nwee2,nweh2,xk,Wsc,beta,opt,opt0,bopt,&
        xw,yk,numeh,g2,vertex,matrix,IPIV,mu,divw,E1d,a1d,num)
! --------------------------------------------------------------------

      
! convergence check --------------------------------------------------
      if(abs(numeh(1)-numehold(1))/numeh(1).lt.cnv.and.&
      abs(numeh(2)-numehold(2))/numeh(2).lt.cnv.and.&
      abs(numeh(1)-numeh(2))/numeh(1).lt.cnv.and.&
      abs(numeh(1)-num)/num.lt.cnv) then
        write(6,*) 'converged !'
        exit
      else if(icount.ge.iend) then
        write(6,*) 'convergence error'
        exit
      end if
! ---------------------------------------------------------------------

    ! set previous density and Akw
      numehold=numeh
      Akwold=Akw      


700   end do
      close(29)
      close(30)
      close(129)
      close(130)
      close(133)

      
! check difference of Akw ----------------------------------------------
      open(70,file='T0.5n0.1_convergence_check.txt',status='unknown')
      write(70,*) 'E-H num=',numeh(1),numeh(2)
      do ieh1=1,2
        do iw=nwn1,nwg1(ieh1)
          w=xw(iw)
          do ik1=0,nk
            k1=xk(ik1)
            abc=Akw(iw,ik1,ieh1)-Akwold(iw,ik1,ieh1)
            if(abc/2.0d0*pai.gt.1.0d-3) then
              write(70,*) k1*a1d,abc
            end if
          end do
        end do
      end do
      close(70)
! --------------------------------------------------------------------------


! check difference of Akw ----------------------------------------------
      call optical(nk,nwnn1,nwee2,nweh2,xk,Wsc,beta,opt,opt0,bopt,&
        xw,yk,numeh,g2,vertex,matrix,IPIV,mu,divw,E1d,a1d,num)

      call output(nk,nwg1,xk,xw,nwn1,nwe2,nwnn1,nwee2,&
        engkQP,kappa,Akw,numeh,numeh0,mu,self,E1d,a1d)
! ----------------------------------------------------------------------
        
      
      end program Tmatrix1d
! *******************************************************************






      
      
      
! ***********************************************************************
      subroutine SHF(mu,beta,num,nk,engkSHF,engkf,fermi0,dmu0,&
        yk,Wsc,Vtot,V0pt,VCH,n0ptk,x0ptk,y0ptk,divk,xk,&
        nCH,xCH,yCH,self0,cnv,const,iend,E1d,a1d)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0
      integer(4) :: i,ieh,ik,ik1,ik2,nk,icount,n0ptk,nCH,iend
      real(8) :: mu(2),engkSHF(0:nk,2),engkf(0:nk,2),fermi0(0:nk,2),&
           yk(0:nk),Wsc(-2*nk:2*nk),Vtot(2*nk),V0pt(n0ptk),VCH(nCH),&
           x0ptk(n0ptk),y0ptk(n0ptk),xk(0:nk),&
           xCH(nCH),yCH(nCH),self0(0:nk,2),numeh(2),numehold(2)
      real(8) :: beta,num,const,dmu0,ffnc,kappa,cffInt,divk,selfCH,CH,&
           abc,k1,cnv,screening,E1d,a1d

    ! initial energy and fermi distribution
      do ieh=1,2 
        do ik=0,nk ! k
          engkSHF(ik,ieh)=engkf(ik,ieh)
          fermi0(ik,ieh)=ffnc(beta,engkSHF(ik,ieh),mu(ieh))
        end do
      end do

      selfCH=0.0d0
      write(6,*) 'SHF e-h particle density (1/a1d)'
! SHF iteration -------------------------------------------------
      icount=0
      do 600
        icount=icount+1
      ! set mueh
        do ieh=1,2 
          call shiftmuSHF(mu,ieh,dmu0,num,selfCH,beta,nk,xk,yk,engkSHF)
          do ik=0,nk ! k
            fermi0(ik,ieh)=ffnc(beta,engkSHF(ik,ieh),mu(ieh))
          end do
        end do

      ! effective interaction
        kappa=screening(mu,beta,fermi0,nk,xk,yk,kappa)
        cffInt=const*kappa/16.0d0/num
        call Intexp(1,nk,Wsc,divk,kappa,&
        n0ptk,x0ptk,y0ptk,Vtot,V0pt,cffInt)

      ! Coulomb hole enery
        selfCH=CH(kappa,nCH,xCH,yCH,VCH,cffInt)

      ! exchange energy
        self0=0.0d0
        do ieh=1,2 
          do ik2=0,nk ! k2
            abc=yk(ik2)*fermi0(ik2,ieh)
            do ik1=0,nk ! k1
              self0(ik1,ieh)=self0(ik1,ieh)+&
                (Wsc(ik1-ik2)+Wsc(ik1+ik2))*abc
            end do
          end do
        end do
        self0=-self0/(2.0d0*pai)+selfCH

      ! total energy
        do ieh=1,2 
          do ik=0,nk ! k
            engkSHF(ik,ieh)=engkf(ik,ieh)+self0(ik,ieh)
          end do
        end do

      ! total density
        numeh=0.0d0
        do ieh=1,2 
           abc=0.0d0
           do ik1=0,nk ! k1
              abc=abc+yk(ik1)*ffnc(beta,engkSHF(ik1,ieh),mu(ieh))
           end do
           numeh(ieh)=abc
        end do
        numeh=2.0d0*numeh/pai

        write(6,*) icount,numeh(1)*a1d,numeh(2)*a1d

      ! convergence check
        if(abs(numeh(1)-numehold(1))/numeh(1).lt.cnv.and.&
             abs(numeh(2)-numehold(2))/numeh(2).lt.cnv.and.&
             abs(numeh(1)-numeh(2))/numeh(1).lt.cnv.and.&
             abs(numeh(1)-num)/num.lt.cnv) exit

        
        if(icount.eq.iend) then
          write(6,*) 'convergence error'
          stop
        end if
        
        numehold=numeh

600   end do
!-------------------------------------------------------------------

! output SHF results -----------------------------------------------------
      write(6,*) 'SHF results'
      write(6,*) 'total mu (E1d)'
      write(6,*) Esc*(mu(1)+mu(2))/E1d
      write(6,*) 'BGR (E1d)'
      write(6,*) (self0(0,1)+self0(0,2))*Esc/E1d
      write(6,*) 'Coulomb Hole (E1d)'
      write(6,*) selfCH*Esc/E1d

      open(27,file='T0.5n0.1_SHF_results.txt',status='unknown')
      write(27,*) 'screening parameter',kappa
      write(27,*) 'Coulomb hole self energy',selfCH*Esc/E1d
      write(27,*) 'E-H screened exchange self energy',&
           (self0(0,1)-selfCH)*Esc/E1d,(self0(0,2)-selfCH)*Esc/E1d
      write(27,*) 'BGR (E1d)',&
           (self0(0,1)+self0(0,2))*Esc/E1d
! -----------------------------------------------------------------

      end subroutine SHF
! ***********************************************************************
            
      
      
      
      
      
      
      
      
      
! ***************************************************************************
      subroutine shiftmuSHF(mu,ieh,dmu0,num,selfCH,beta,nk,xk,yk,engk)
      implicit none
      integer(4) :: i,j,sign,nk,ik1,ieh
      real(8),parameter :: pai=dacos(-1.0d0),hdig=1.0d-10
      real(8) :: &
        mu(2),mu2,num,dmu0,beta,mumax,mumin,dense,selfCH,ffnc
      real(8) :: xk(0:nk),yk(0:nk),engk(0:nk,2)
      
! electron ------------------------------------------------------------
! check upper limit
      dense=0.0d0
      do ik1=0,nk
        dense=dense+yk(ik1)*ffnc(beta,engk(ik1,ieh),mu(ieh))
      end do
      dense=2.0d0*dense/pai

      if(dense.gt.num) then
        sign=-1
        mumax=mu(ieh)
      else
        sign=1
        mumin=mu(ieh)
      end if
! narrow down the chemical potential
      i=0
      do
        i=i+sign
        mu2=mu(ieh)+dble(i)*dmu0
        
        dense=0.0d0
        do ik1=0,nk
          dense=dense+yk(ik1)*ffnc(beta,engk(ik1,ieh),mu2)
        end do
        dense=2.0d0*dense/pai
        if(dense.lt.num) then
          mumin=mu2
          if(sign.eq.-1) exit
        else
          mumax=mu2
          if(sign.eq.1) exit
        end if
      end do
      
      i=0
      do
        i=i+1
        mu2=(mumax+mumin)/2.0d0
      
      ! SHF density
        dense=0.0d0
        do ik1=0,nk
          dense=dense+yk(ik1)*ffnc(beta,engk(ik1,ieh),mu2)
        end do
        dense=2.0d0*dense/pai

        if(abs(dense-num)/num.lt.hdig) then
          exit
        else if(dense.lt.num) then
          mumin=mu2
        else
          mumax=mu2
        end if
        if(i.gt.100) then
! chemical potential located in the energy gap of finite size mesh pt.
          write(6,*) 'shiftmuSHF error'
          stop
        end if
      end do
      mu(ieh)=mu2

      return
      end subroutine shiftmuSHF
!************************************************************************



!**************************************************************************** 
      subroutine shiftmu(num,ieh,mu,dmu0,numeh,&
           beta,Dos,xw,divw,nwn1,nwe2,nwnn1,nwee2,nwmax)
      implicit none
      integer(4) :: i,j,sign,iw,nwmax,nwn1,nwe2,nwnn1,nwee2,ieh
      real(8),parameter :: pai=dacos(-1.0d0),hdig=1.0d-7
      real(8) :: &
        mu0,mu2,num,dmu0,beta,mumax,mumin,dense
      real(8) :: Dos(nwn1:nwe2,2),xw(nwnn1:nwee2),divw,&
      ffnc,w,mu(2),numeh(2)
      
! electron ------------------------------------------------------------
! check upper limit
      mu0=mu(ieh)
      dense=0.0d0
      do iw=nwn1,nwmax
        w=xw(iw)
        dense=dense+Dos(iw,ieh)*ffnc(beta,w,mu0)
      end do
      dense=dense*divw

      if(dense.gt.num) then
        sign=-1
        mumax=mu0
      else
        sign=1
        mumin=mu0
      end if
! narrow down the chemical potential
      i=0
      do
        i=i+sign
        mu2=mu0+dble(i)*dmu0
        
        dense=0.0d0
        do iw=nwn1,nwmax
          w=xw(iw)
          dense=dense+Dos(iw,ieh)*ffnc(beta,w,mu2)
        end do
        dense=dense*divw
        
        if(dense.lt.num) then
          mumin=mu2
          if(sign.eq.-1) exit
        else
          mumax=mu2
          if(sign.eq.1) exit
        end if
        
      end do
      
      i=0
      do
        i=i+1
        mu2=(mumax+mumin)/2.0d0
        
        dense=0.0d0
        do iw=nwn1,nwmax
          w=xw(iw)
          dense=dense+Dos(iw,ieh)*ffnc(beta,w,mu2)
        end do
        dense=dense*divw
        
        if(abs(dense-num)/num.lt.hdig) then
          exit
        else if(dense.lt.num) then
          mumin=mu2
        else
          mumax=mu2
        end if
        if(i.gt.100) then
          write(6,*) 'shiftmu error'
          stop
        end if
      end do
      mu(ieh)=mu2
      numeh(ieh)=dense

      return
      end subroutine shiftmu
!*********************************************************************




!**************************************************************************
      function screening(mu,beta,fermi0,nk,xk,yk)
      implicit none
      integer(4) :: nk,ik,ieh
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0
      real(8) :: beta,screening,k,fermi0(0:nk,2),mu(2),&
        xk(0:nk),yk(0:nk),E1d
      
      screening=0.0d0
      do ieh=1,2
        do ik=0,nk
          k=xk(ik)
          screening=screening+yk(ik)*fermi0(ik,ieh)*&
          (1.0d0-fermi0(ik,ieh))
        end do
      end do
      screening=2.0d0*beta*screening/pai

      return
      end function screening
!******************************************************************************






! ***************************************************************************
      subroutine Intexp(isc,nk,Wsc,divk,kappa,&
        n0ptk,x0ptk,y0ptk,Vtot,V0pt,cffInt)
      implicit none
      integer(4) :: nk,iq,i,ieh1,ieh2,n0ptk,isc
      real(8),parameter :: pai=dacos(-1.0d0)
      real(8) :: kappa,k1,k2,S1,q,divk,Wscq0
      real(8) :: &
        Wsc(-2*nk:2*nk),Vtot(2*nk),V0pt(n0ptk),&
        sign(2,2),x0ptk(n0ptk),y0ptk(n0ptk),abc,cffInt
      Wsc=0.0d0
   ! interaction (no 0poit correction)
      do iq=1,2*nk
        q=divk*dble(iq)
        abc=1.0d0/Vtot(iq)+kappa/(1.0d0+q**2*cffInt)
        Wsc(iq)=1.0d0/abc
        Wsc(-iq)=Wsc(iq)
      end do

    ! interaction with 0point correction at iq=0
      S1=0.0d0
      do iq=1,n0ptk
        q=x0ptk(iq)
        S1=S1+y0ptk(iq)/(1.0d0/V0pt(iq)+&
        kappa/(1.0d0+q**2*cffInt))
      end do
      Wsc(0)=2.0d0*(S1-divk*Wsc(1)/2.0d0)/divk
      
      return
      end subroutine Intexp
!*************************************************************************





!***************************************************************************
      subroutine Int2nd(isc,nk,WscV0,Wsc,divk,kappa,&
        n0ptk,x0ptk,y0ptk,Vtot,V0pt,cffInt)
      implicit none
      integer(4) :: nk,iq,i,ieh1,ieh2,n0ptk,isc
      real(8),parameter :: pai=dacos(-1.0d0)
      real(8) :: kappa,k1,k2,S1,q,divk,abc,cffInt
      real(8) :: WscV0(-2*nk:2*nk),Wsc(-2*nk:2*nk),&
        sign(2,2),x0ptk(n0ptk),y0ptk(n0ptk),Vtot(2*nk),V0pt(n0ptk)
      
      WscV0=0.0d0
   ! interaction (no 0poit correction)
      do iq=1,2*nk
        q=divk*dble(iq)
        abc=1.0d0/Vtot(iq)+kappa/(1.0d0+q**2*cffInt)
        WscV0(iq)=Vtot(iq)/abc
        WscV0(-iq)=Vtot(iq)/abc
      end do
      
    ! interaction with 0point correction at iq=0      
      S1=0.0d0
      do iq=1,n0ptk
        q=x0ptk(iq)
        S1=S1+V0pt(iq)*y0ptk(iq)/(1.0d0/V0pt(iq)+&
        kappa/(1.0d0+q**2*cffInt))
      end do
      WscV0(0)=2.0d0*(S1-divk*WscV0(1)/2.0d0)/divk
      
      
      return
      end subroutine Int2nd
!*************************************************************************



!*************************************************************************
      function bfnc(beta,ek)
      implicit none
      real(8) :: bfnc,beta,ek
         bfnc=1.0d0/(dexp(beta*ek)-1.0d0)
      end function bfnc
!*************************************************************************
!*************************************************************************
      function ffnc(beta,ek,mu)
      implicit none
      real(8) :: ffnc,beta,ek,mu
         ffnc=1.0d0/(dexp(beta*(ek-mu))+1.0d0)
      end function ffnc
!*************************************************************************





!***********************************************************************
      function CH(kappa,nCH,xCH,yCH,VCH,cffInt)
      implicit none
      integer(4) :: nCH,i
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0
      real(8) :: &
      Opl2,kappa,u,CH,xCH(nCH),yCH(nCH),q,abc,&
      VCH(nCH),cffInt
      
      CH=0.0d0
      do i=1,nCH
        q=xCH(i)
        abc=VCH(i)**2*kappa/&
        (1.0d0+VCH(i)*kappa+q**2*cffInt)
        
        CH=CH+yCH(i)*abc
        if(abc/CH.le.1.0d-10/Esc) exit
      end do
      CH=-CH/2.0d0/pai
      
      return
      end function CH
!*************************************************************************



!*************************************************************************
      SUBROUTINE gauleg(x1,xk,x,w,n)
      real(8) :: pai=dacos(-1.0d0)
      INTEGER(4) :: n
      REAL(8) :: x1,xk,x(n),w(n)
      REAL(8),PARAMETER :: EPS=3.d-14
      INTEGER(4) :: i,j,m
      REAL(8) :: p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(xk+x1)
      xl=0.5d0*(xk-x1)
      do 12 i=1,m
        z=cos(pai*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END SUBROUTINE gauleg
!*************************************************************************











!******************************************************************
      subroutine Dosnumber(Akw,engkf,engkQP,nk,nwg1,xw,divw,&
      nwn1,nwe2,nwnn1,nwee2,self,Dos,yk,numeh,numeh0,&
      mu,beta,E1d,a1d)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0
      complex(8),parameter :: im=(0.0d0,1.0d0)
      integer(4) :: &
      nk,nwg1(2),nwn1,nwe2,nwnn1,nwee2,ieh1,ik1,iw,iw1,iw2,&
      iph,icent0,iwph
      real(8) :: &
      Akw(nwn1:nwe2,0:nk,2),&
      engkf(0:nk,2),engkQP(0:nk,2),xw(nwnn1:nwee2),&
      divw,Dos(nwn1:nwe2,2),xk(0:nk),&
      yk(0:nk),numeh(2),numeh0(2),mu(2),&
      beta,a,w,ww,wph,abc,def,ghi,jkl,ffnc,k1,E1d,a1d
      complex(8) :: self(nwn1:nwe2,0:nk,2)

    ! Akw
      Akw=0.0d0
      do ieh1=1,2
        do ik1=0,nk
          do iw=nwn1,nwg1(ieh1)
            w=xw(iw)
            Akw(iw,ik1,ieh1)=&
            dimag(1.0d0/(w-engkf(ik1,ieh1)-self(iw,ik1,ieh1)))
          end do 
        end do 
      end do 
      Akw=-2.0d0*Akw

    ! Dos
      Dos=0.0d0
      do ieh1=1,2
        do iw=nwn1,nwg1(ieh1)
          do ik1=0,nk
            Dos(iw,ieh1)=Dos(iw,ieh1)+yk(ik1)*Akw(iw,ik1,ieh1)
          end do
        end do
      end do
      Dos=2.0d0*Dos/(2.0d0*pai)/pai
      
    ! total density
      numeh=0.0d0
      do ieh1=1,2
        do iw=nwn1,nwg1(ieh1)
          w=xw(iw)
          numeh(ieh1)=numeh(ieh1)+Dos(iw,ieh1)*ffnc(beta,w,mu(ieh1))
        end do
      end do
      numeh=divw*numeh
      
      
      
    ! quasi-particle dispersion 
      do ieh1=1,2
        do ik1=0,nk
          do iw=nwg1(ieh1)-1,nwn1,-1
            w=xw(iw)
            abc=(w-(engkf(ik1,ieh1)+dreal(self(iw,ik1,ieh1))))
            if(abc.lt.0.0d0) then
              ww=xw(iw+1)
              def=(ww-(engkf(ik1,ieh1)+dreal(self(iw+1,ik1,ieh1))))
              ghi=-abc/(def-abc)*(ww-w)
              exit
            end if
          end do
          engkQP(ik1,ieh1)=w+ghi
        end do
      end do
      
      open(44,file='T0.5n0.1_qp_pole_check.txt',status='unknown')
      do ieh1=1,2
        do iw=nwn1,nwg1(ieh1)
          ww=xw(iw)
          def=(ww-(engkf(0,ieh1)+dreal(self(iw,0,ieh1))))
          write(44,*) ww,def
        end do
      end do

    ! quasi-particle density
      numeh0=0.0d0
      do ieh1=1,2
        do ik1=0,nk
          numeh0(ieh1)=numeh0(ieh1)+&
          yk(ik1)*ffnc(beta,engkQP(ik1,ieh1),mu(ieh1))
        end do
      end do
      numeh0=2.0d0*numeh0/pai

      close(44)


! Akw check -------------------------------------------------------
      open(20,file='T0.5n0.1_Akwshift_e.txt',status='unknown')
      open(21,file='T0.5n0.1_Akwshift_h.txt',status='unknown')

      do ik1=0,nk,nk/25
        k1=xk(ik1)
        write(20,*) 'k=',k1*a1d
        write(21,*) 'k=',k1*a1d
        ieh1=1
        do iw=nwn1,nwg1(ieh1)
          w=xw(iw)
          write(20,*) w*Esc/E1d,Akw(iw,ik1,ieh1)*E1d/Esc+k1*a1d
        end do
        ieh1=2
        do iw=nwn1,nwg1(ieh1)
          w=xw(iw)
          write(21,*) w*Esc/E1d,Akw(iw,ik1,ieh1)*E1d/Esc+k1*a1d
        end do
      end do

      close(20)
      close(21)
! ---------------------------------------------------------------

! Akw check2 (k=0) ------------------------------------------
      open(33,file='T0.5n0.1_Akwk0_e.txt',status='unknown')
      open(63,file='T0.5n0.1_Akwk0_h.txt',status='unknown')
      ieh1=1
      ik1=0
      k1=xk(0)
      write(33,*) 'k=',k1*a1d
      write(63,*) 'k=',k1*a1d
      do iw=nwn1,nwg1(ieh1)
        w=xw(iw)
        write(33,*) w*Esc/E1d,Akw(iw,ik1,ieh1)*E1d/Esc
      end do
      ieh1=1
      do iw=nwn1,nwg1(ieh1)
        w=xw(iw)
        write(63,*) w*Esc,Akw(iw,ik1,ieh1)*E1d/Esc
      end do
      close(33)
      close(63)
! --------------------------------------------------------------

! Dos check -------------------------------------------------------
      open(45,file='T0.5n0.1_Dos_e.txt',status='unknown')
      open(46,file='T0.5n0.1_Dos_h.txt',status='unknown')
      ieh1=1
      do iw=nwn1,nwg1(ieh1)
        w=xw(iw)
        write(45,*) w*Esc/E1d,Dos(iw,ieh1)/Esc*E1d
      end do
      ieh1=2
      do iw=nwn1,nwg1(ieh1)
        w=xw(iw)
        write(46,*) w*Esc/E1d,Dos(iw,ieh1)/Esc*E1d
      end do
      close(45)
      close(46)
! ------------------------------------------------------------------


! qp check ------------------------------------------------------
      open(43,file='T0.5n0.1_qp_dispersion_e.txt',status='unknown')
      open(44,file='T0.5n0.1_qp_dispersion_h.txt',status='unknown')
      ieh1=1
      do ik1=0,nk
        k1=xk(ik1)
        write(43,*) k1*a1d,engkQP(ik1,ieh1)*Esc/E1d
      end do
      ieh1=2
      do ik1=0,nk
        k1=xk(ik1)
        write(44,*) k1*a1d,engkQP(ik1,ieh1)*Esc/E1d
      end do
      close(43)
      close(44)
! ----------------------------------------------------------------

! fermi distribution check -------------------------------------------
      open(176,file='T0.5n0.1_E_distribution.txt',status='unknown')
      open(177,file='T0.5n0.1_H_distribution.txt',status='unknown')
      do iw=nwn1,nwg1(1)
        w=xw(iw)
        write(176,*) w*Esc/E1d,ffnc(beta,w,mu(1))
      end do
      do iw=nwn1,nwg1(2)
        w=xw(iw)
        write(177,*) w*Esc/E1d,ffnc(beta,w,mu(2))
      end do
      close(176)
      close(177)
! ----------------------------------------------------------------

      return
      end subroutine Dosnumber
!******************************************************************








! ************************************************************************
      subroutine output(nk,nwg1,xk,xw,nwn1,nwe2,nwnn1,nwee2,&
        engkQP,kappa,Akw,numeh,numeh0,mu,self,E1d,a1d)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0
      integer(4) :: &
      ieh1,ieh2,ik1,nk,iw,nwg1(2),nwn1,nwe2,nwnn1,nwee2,&
      nave,i,j
      real(8) :: xk(0:nk),xw(nwnn1:nwee2),&
      engkQP(0:nk,2),Akw(nwn1:nwe2,0:nk,2),&
      numeh(2),numeh0(2),mu(2),k1,w,kappa,E1d,a1d
      complex(8) :: self(nwn1:nwe2,0:nk,2)

! selfenergy output for read-in ------------------------------------------
      open(89,file='T0.5n0.1_results.txt',status='unknown')
      write(89,*) mu(1),mu(2)

      do ieh1=1,2
        do ik1=0,nk
          do iw=nwn1,nwg1(ieh1)
            write(89,*) self(iw,ik1,ieh1)
          end do
        end do
      end do
      close(89)
! ----------------------------------------------------------------------------

! output Akw (contour plot) -------------------------------------------------
      open(25,file='T0.5n0.1_Akwcontour_e.txt',status='unknown')
      open(26,file='T0.5n0.1_Akwcontour_h.txt',status='unknown')
      ieh1=1
      do iw=nwn1,nwg1(ieh1)
        w=xw(iw)
        if(abs(w).gt.10.0d0) goto 20
        do ik1=0,nk
          k1=xk(ik1)
          if(abs(k1).gt.10.0d0) exit
          write(25,'(3F20.15)') k1*a1d,w*Esc/E1d,&
          Akw(iw,ik1,ieh1)/Esc*E1d

        end do
        write(25,*)
20    end do
      ieh1=2
      do iw=nwn1,nwg1(ieh1)
        w=xw(iw)
        if(abs(w).gt.10.0d0) goto 21
        do ik1=0,nk
          k1=xk(ik1)
          if(abs(k1).gt.10.0d0) exit
          write(26,'(3F20.15)') k1*a1d,w*Esc/E1d,&
          Akw(iw,ik1,ieh1)/Esc*E1d

        end do
        write(26,*)
21    end do
      close(25)
      close(26)
! -------------------------------------------------------------------


! thermodynamic results -----------------------------------------------------
      open(47,file='T0.5n0.1_thermodynamic_results.txt',status='unknown')
      write(47,*) 'E chemical potential (E1d)'
      write(47,*) (numeh(1)+numeh(2))*a1d/2.0d0,mu(1)*Esc/E1d
      write(47,*) 'H chemical potential (E1d)'
      write(47,*) (numeh(1)+numeh(2))*a1d/2.0d0,mu(2)*Esc/E1d
      write(47,*) 'total chemical potential (E1d)'
      write(47,*) (numeh(1)+numeh(2))*a1d/2.0d0,(mu(1)+mu(2))*Esc/E1d
      write(47,*)

      write(47,*) 'E-H carrier density (1/a1d)'
      write(47,*) numeh(1)*a1d,numeh(2)*a1d,&
      (numeh(1)+numeh(2))*a1d/2.0d0
      write(47,*)

      write(47,*) 'electron ionization ratio'
      write(47,*) (numeh(1)+numeh(2))*a1d/2.0d0,&
      numeh0(1)/numeh(1)
      write(47,*) 'hole ionization ratio'
      write(47,*) (numeh(1)+numeh(2))*a1d/2.0d0,&
      numeh0(2)/numeh(2)
      write(47,*) 'mean ionization ratio'
      write(47,*) (numeh(1)+numeh(2))*a1d/2.0d0,&
      (numeh0(1)/numeh(1)+numeh0(2)/numeh(2))/2.0d0

      write(47,*)
      write(47,*) 'screening parameter'
      write(47,*) (numeh(1)+numeh(2))*a1d/2.0d0,kappa

      write(47,*) 'BGR (E1d)'
      write(47,*) (numeh(1)+numeh(2))*a1d/2.0d0,&
        (engkQP(0,1)+engkQP(0,2))*Esc/E1d

      close(47)
! --------------------------------------------------------------------------

      return
      end subroutine output
! ****************************************************************************









!*********************************************************************
      subroutine optical(nk,nwnn1,nwee2,nweh2,xk,Wsc,beta,opt,opt0,bopt,&
        xw,yk,numeh,g2,vertex,matrix,IPIV,mu,divw,E1d,a1d,num)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0
      integer(4) :: nk,IPIV(2*nk+1),nwnn1,nwee2,nweh2,i
      real(8) :: xk(0:nk),Wsc(-2*nk:2*nk),mu(2),b,beta,E1d,a1d,&
        xw(nwnn1:nwee2),yk(0:nk),numeh(2),num,bopt(nwnn1:nweh2)
      complex(8) :: g2(nwnn1:nweh2,0:nk),vertex(2*nk+1),matrix(2*nk+1,2*nk+1)
      
      integer(4) :: iq,ik1,ik2,iw,INFO,iw0
      real(8) :: w,absorb,bfnc,abc,def,mu0,divw,peakw,minw,maxw
      complex(8) :: x,zabc,zdef,opt(nwnn1:nweh2),opt0(nwnn1:nweh2)

      do iw=nwnn1,nweh2
        w=xw(iw)
        
        matrix=0.0d0
        do ik2=0,nk
          matrix(ik2+1,ik2+1)=1.0d0
          zabc=yk(ik2)*g2(iw,ik2)/(2.0d0*pai)
          do ik1=0,nk
            matrix(ik1+1,ik2+1)=matrix(ik1+1,ik2+1)&
            -(Wsc(ik1-ik2)+Wsc(ik1+ik2))*zabc
          end do
        end do
        
        vertex=1.0d0
        call ZGESV(nk+1,1,matrix,2*nk+1,IPIV,vertex,2*nk+1,INFO)
        if(INFO.ne.0) then
          write(6,*) 'ZGESV error INFO=',INFO
          stop
        end if
        
        zabc=0.0d0
        zdef=0.0d0
        do ik1=0,nk
          zabc=zabc+yk(ik1)*vertex(ik1+1)*g2(iw,ik1)
          zdef=zdef+yk(ik1)*g2(iw,ik1)
        end do
        
        opt(iw)=zabc
        opt0(iw)=zdef
      end do
      opt=2.0d0*opt/pai
      opt0=2.0d0*opt0/pai


! output optical response ----------------------------------------------
      open(222,file='T0.5n0.1_Imx.txt',status='unknown')
      open(125,file='T0.5n0.1_Rex.txt',status='unknown')
      open(123,file='T0.5n0.1_PL.txt',status='unknown')
      open(32,file='T0.5n0.1_Imx0.txt',status='unknown')
      open(132,file='T0.5n0.1_Rex0.txt',status='unknown')
      open(35,file='T0.5n0.1_PL0.txt',status='unknown')
      open(42,file='T0.5n0.1_Imx_per_Imx0.txt',status='unknown')
      open(142,file='T0.5n0.1_Rex_per_Rex0.txt',status='unknown')

      do iw=nwnn1,nweh2
        w=xw(iw)
        if(abs(w*Esc/E1d).lt.20.0d0) then

        write(222,*) w*Esc/E1d,dimag(opt(iw))*a1d*E1d/Esc
        write(125,*) w*Esc/E1d,dreal(opt(iw))*a1d*E1d/Esc
        write(123,*) &
        w*Esc/E1d,bfnc(beta,w-mu(1)-mu(2))*dimag(opt(iw))*a1d*E1d/Esc
        write(32,*) w*Esc/E1d,dimag(opt0(iw))*a1d*E1d/Esc
        write(132,*) w*Esc/E1d,dreal(opt0(iw))*a1d*E1d/Esc
        write(35,*) &
          w*Esc/E1d,&
          bfnc(beta,w-mu(1)-mu(2))*dimag(opt0(iw))*a1d*E1d/Esc
        write(42,*) w*Esc/E1d,&
          dimag(opt(iw))/dimag(opt0(iw))*a1d*E1d/Esc
        write(142,*) w*Esc/E1d,&
          dreal(opt(iw))/dreal(opt0(iw))*a1d*E1d/Esc
        end if
      end do
      close(222)
      close(125)
      close(123)
      close(32)
      close(132)
      close(35)
      close(42)
      close(142)
! ------------------------------------------------------------------


! output pair susceptibility --------------------------------------------
      open(56,file='T0.5n0.1_pair_susceptibility.txt',status='unknown') 
      mu0=(mu(1)+mu(2))
      do iw=nwnn1,nweh2
        w=xw(iw)
        if(w-mu0.gt.0) then
          iw0=iw-1
          exit
        end if
      end do
      abc=(dreal(opt(iw0+1))-dreal(opt(iw0)))/divw*(mu0-xw(iw0))+&
           dreal(opt(iw0))
      def=(dreal(opt0(iw0+1))-dreal(opt0(iw0)))/divw*(mu0-xw(iw0))+&
           dreal(opt0(iw0))

      write(56,*) "Rex at w=mu (1/(a1d*E1d))" 
      write(56,*) num*a1d,abc*a1d*E1d/Esc,1.0d0/beta*Esc/E1d
      write(56,*) "Rex0 at w=mu (1/(a1d*E1d))"
      write(56,*) num*a1d,def*a1d*E1d/Esc,1.0d0/beta*Esc/E1d
      write(56,*) "Rex/Rex0 at w=mu" 
      write(56,*) num*a1d,abc/def,1.0d0/beta*Esc/E1d
      write(56,*) "Rex-Rex0 at w=mu (1/(a1d*E1d))"
      write(56,*) num*a1d,(abc-def)*a1d*E1d/Esc,1.0d0/beta*Esc/E1d
      close(56)
! ---------------------------------------------------------------------


! peak position ------------------------------------------------------
      open(30,file='T0.5n0.1_peak_position.txt',status='unknown')
      call peakposition(dimag(opt),xw,nwnn1,nweh2,nwee2,peakw,minw,maxw)
      write(30,*) "absorption_peak_position"
      write(30,*) num*a1d,peakw*Esc/E1d,minw*Esc/E1d,maxw*Esc/E1d
      call peakposition(-dimag(opt),xw,nwnn1,nweh2,nwee2,peakw,minw,maxw)
      write(30,*) "gain_peak_position"
      write(30,*) num*a1d,peakw*Esc/E1d,minw*Esc/E1d,maxw*Esc/E1d
      call peakposition(bopt,xw,nwnn1,nweh2,nwee2,peakw,minw,maxw)
      write(30,*) "PL_peak_position"
      write(30,*) num*a1d,peakw*Esc/E1d,minw*Esc/E1d,maxw*Esc/E1d
! ---------------------------------------------------------------------


      end subroutine optical
!*********************************************************************


! **************************************************************************
      subroutine peakposition(spec,xw,nwnn1,nweh2,nwee2,peakw,minw,maxw)
      implicit none
      real(8) :: spec(nwnn1:nweh2),xw(nwnn1:nwee2),maxspec,peakw,minw,maxw,w
      integer(4) :: nwnn1,nweh2,nwee2,imin,i
      maxspec=0.0d0
      do i=nwnn1,nweh2
        if(spec(i).gt.maxspec) then
          maxspec=spec(i)
          peakw=xw(i)
        end if
      end do

      do i=nwnn1,nweh2
        if(spec(i).gt.0.5d0*maxspec) then
          minw=(xw(i)-xw(i-1))/(spec(i)-spec(i-1))*&
          (0.5d0*maxspec-spec(i-1))+xw(i-1)
          imin=i
          exit
        end if
      end do

      do i=imin,nweh2
        if(spec(i).lt.0.5d0*maxspec) then
          maxw=(xw(i)-xw(i-1))/(spec(i)-spec(i-1))*&
          (0.5d0*maxspec-spec(i-1))+xw(i-1)
          exit
        end if
      end do
      end subroutine peakposition
! **************************************************************************





!*********************************************************************
      subroutine &
      Impair(gkk,Imbgkk,nweh2,nk,nwn1,nwe2,nwnn1,nwee2,&
        xw,divw,nwg1,nwg2,Akw,pblock,fermi,xk,E1d,a1d)
      implicit none
      complex(8),parameter :: im=(0.0d0,1.0d0)
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0
      integer(4) :: &
        ik,ik1,ik2,nk,iw,iw1,iw2,nwn1,nwe2,nwg2(2,2),nwg1(2),&
        ieh,ieh1,ieh2,&
        nwnn1,nwee2,nweh2
      real(8) :: &
        xw(nwnn1:nwee2),divw,w,ffnc,S1,xk(0:nk),abc,&
        Akw(nwn1:nwe2,0:nk,2),&
        pblock(nwn1:nwe2,nwn1:nwe2,2,2),fermi(nwn1:nwe2,2),k,E1d,a1d
      complex(8) :: gkk(nwnn1:nwee2,0:nk,0:nk,2,2)
      real(8) :: Imbgkk(nwnn1:nwee2,0:nk,0:nk,2,2)

      gkk=0.0d0
      Imbgkk=0.0d0
      do ieh1=1,2
        do ieh2=1,2
          do ik2=0,nk
            do ik1=0,ik2
              do iw1=nwn1,nwg1(ieh1)
!$omp parallel do private(iw2,iw,abc)
                do iw2=nwn1,nwg1(ieh2)
                  iw=iw1+iw2
                  abc=Akw(iw1,ik1,ieh1)*Akw(iw2,ik2,ieh2)

                  gkk(iw,ik1,ik2,ieh1,ieh2)=gkk(iw,ik1,ik2,ieh1,ieh2)+&
                    pblock(iw1,iw2,ieh1,ieh2)*abc
                ! without bose distribution function
                  Imbgkk(iw,ik1,ik2,ieh1,ieh2)=Imbgkk(iw,ik1,ik2,ieh1,ieh2)+&
                    fermi(iw1,ieh1)*fermi(iw2,ieh2)*abc
                end do
!$omp end parallel do
              end do
            end do
          end do
        end do
      end do
      gkk=-im*divw*gkk/(4.0d0*pai)
      Imbgkk=-divw*Imbgkk/(4.0d0*pai)


! check pair Green's function -------------------------------------
      open(75,file='T0.5n0.1_Imgkk_ee.txt',status='unknown')
      open(76,file='T0.5n0.1_Imgkk_eh.txt',status='unknown')
      open(78,file='T0.5n0.1_Imgkk_hh.txt',status='unknown')
      do ik=0,nk,nk/25 
        k=xk(ik)
        write(75,*) 'k=',k*a1d
        write(76,*) 'k=',k*a1d
        write(78,*) 'k=',k*a1d
        do iw=nwnn1,nwg2(1,1)
          w=xw(iw)
          write(75,*) w*Esc/E1d,&
            -dimag(gkk(iw,ik,ik,1,1))*E1d/Esc+k*a1d*0.25d0
        end do
        do iw=nwnn1,nwg2(2,1)
          w=xw(iw)
          write(76,*) w*Esc/E1d,&
            -dimag(gkk(iw,ik,ik,2,1))*E1d/Esc+k*a1d*0.25d0
        end do
        do iw=nwnn1,nwg2(2,2)
          w=xw(iw)
          write(78,*) w*Esc/E1d,&
            -dimag(gkk(iw,ik,ik,2,2))*E1d/Esc+k*a1d*0.25d0
        end do
      end do
      close(75)
      close(76)
      close(78)
! --------------------------------------------------------

      return
      end subroutine Impair
!*********************************************************************
      






!*********************************************************************
      subroutine &
      TmatrixMW(ImT,ImbT,gkk,Imbgkk,nwnn1,nk,divk,nkbd,nknum,&
      Wsc,Wsc2nd,IPIV,ImMW,ImbMW,dmatrix,matrix,dT2,T2,&
      nwee2,nwg2,xk,xw,E1d,a1d)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0
      real(8),parameter :: sgn(3)=(/1.0d0,-1.0d0,1.0d0/)
      real(8),parameter :: cffTex(3)=(/1.0d0,0.0d0,1.0d0/)
      integer(4) :: &
        nwnn1,nwee2,nk,iw,iq,i1,i2,i3,ik1,ik3,ik2,&
        nkbd(0:2*nk),nwg2(2,2),&
        nknum(0:2*nk),INFO,IPIV(2*nk+1),ieh1,ieh2,ieh
      real(8) :: &
        ImT(nwnn1:nwee2,-nk:nk,-nk:nk,2,2),ImbMW(2*nk+1),&
        ImbT(nwnn1:nwee2,-nk:nk,-nk:nk,2,2),divk,ImMW(2*nk+1),&
        Wsc(-2*nk:2*nk),Wsc2nd(-2*nk:2*nk),abc,xk(0:nk),&
        xw(nwnn1:nwee2),k1,w,E1d,a1d,Imbgkk(nwnn1:nwee2,0:nk,0:nk,2,2),&
        dmatrix(2*nk+1,2*nk+1),dT2(2*nk+1,2*nk+1),def
      complex(8) :: &
        matrix(2*nk+1,2*nk+1),T2(2*nk+1,2*nk+1),&
        gkk(nwnn1:nwee2,0:nk,0:nk,2,2),zabc
        
      
      IPIV=0
      ImT=0.0d0
      ImbT=0.0d0
    ! Tmatrix contribution
      ieh=0
      do ieh2=1,2
        do ieh1=ieh2,2
          ieh=ieh+1
!$omp parallel do private(iw,iq,matrix,dmatrix,ik2,ik1,ik3,zabc,&
          T2,dT2,IPIV,i3,i2,i1,abc)
      do iw=nwnn1,nwg2(ieh1,ieh2)
        do iq=0,2*nk

! Tmatrix ---------------------------------------------------------        
        ! set right hand side
          T2=0.0d0
          do ik2=nkbd(iq),nk
            i2=ik2-nkbd(iq)+1
            do ik1=nkbd(iq),nk
              i1=ik1-nkbd(iq)+1
              T2(i1,i2)=sgn(ieh)*Wsc(ik1-ik2)
            end do
          end do

        ! set left hand side
          matrix=0.0d0
          do ik3=nkbd(iq),nk
            i3=ik3-nkbd(iq)+1
            ik2=iq-ik3 ! pair wavenumber ik3,ik2 in gkk
            zabc=-sgn(ieh)*divk*&
              gkk(iw,abs(ik3),abs(ik2),ieh1,ieh2)/(2.0d0*pai)
            do ik1=nkbd(iq),nk
              i1=ik1-nkbd(iq)+1
              matrix(i1,i3)=Wsc(ik1-ik3)*zabc
            end do
            matrix(i3,i3)=matrix(i3,i3)+1.0d0
          end do

          call ZGESV(nknum(iq),nknum(iq),matrix,2*nk+1,&
          IPIV,T2,2*nk+1,INFO)
          
          do ik1=nkbd(iq),nk
            i1=ik1-nkbd(iq)+1
            ik2=iq-ik1
            i2=ik2-nkbd(iq)+1            
            ImT(iw,ik1,ik2,ieh1,ieh2)=&
                 dimag(2.0d0*T2(i1,i1)-cffTex(ieh)*T2(i1,i2))
          end do

! Tmatrix multipled by bose distribution function ---------------------
        ! set right hand side
          dT2=0.0d0
          do ik2=nkbd(iq),nk
            i2=ik2-nkbd(iq)+1
            do ik3=nkbd(iq),nk
              i3=ik3-nkbd(iq)+1
              abc=sgn(ieh)*dreal(T2(i3,i2))*&
                   Imbgkk(iw,abs(ik3),abs(iq-ik3),ieh1,ieh2)
              do ik1=nkbd(iq),nk
                i1=ik1-nkbd(iq)+1
                dT2(i1,i2)=dT2(i1,i2)+Wsc(ik1-ik3)*abc
              end do
            end do
          end do
          dT2=divk*dT2/(2.0d0*pai)

        ! set left hand side
          dmatrix=0.0d0
          do ik3=nkbd(iq),nk
            i3=ik3-nkbd(iq)+1
            ik2=iq-ik3 ! pair wavenumber ik3,ik2 in gkk
            abc=-sgn(ieh)*divk*&
              dreal(gkk(iw,abs(ik3),abs(ik2),ieh1,ieh2))/(2.0d0*pai)
            do ik1=nkbd(iq),nk
              i1=ik1-nkbd(iq)+1
              dmatrix(i1,i3)=Wsc(ik1-ik3)*abc
            end do
            dmatrix(i3,i3)=dmatrix(i3,i3)+1.0d0
          end do

          call DGESV(nknum(iq),nknum(iq),dmatrix,2*nk+1,&
          IPIV,dT2,2*nk+1,INFO)
          
          do ik1=nkbd(iq),nk
            i1=ik1-nkbd(iq)+1
            ik2=iq-ik1
            i2=ik2-nkbd(iq)+1            
            ImbT(iw,ik1,ik2,ieh1,ieh2)=&
              2.0d0*dT2(i1,i1)-cffTex(ieh)*dT2(i1,i2)
          end do

        end do
      end do
!$omp end parallel do
        end do
      end do

    ! Montroll-Word contribution
      ieh=0
      do ieh2=1,2
        do ieh1=ieh2,2
          ieh=ieh+1
!$omp parallel do private(iw,iq,ImMW,ImbMW,ik3,ik2,abc,def,i1,ik1)
      do iw=nwnn1,nwg2(ieh1,ieh2)
        do iq=0,2*nk
        
          ImMW=0.0d0
          ImbMW=0.0d0
          do ik3=nkbd(iq),nk
            ik2=iq-ik3
            abc=divk*dimag(gkk(iw,abs(ik3),abs(ik2),ieh1,ieh2))/(2.0d0*pai)
            def=divk*Imbgkk(iw,abs(ik3),abs(ik2),ieh1,ieh2)/(2.0d0*pai)
            do ik1=nkbd(iq),nk
              i1=ik1-nkbd(iq)+1
              ImMW(i1)=ImMW(i1)+Wsc2nd(ik1-ik3)*abc
              ImbMW(i1)=ImbMW(i1)+Wsc2nd(ik1-ik3)*def
            end do
          end do
  
          do ik1=nkbd(iq),nk
            i1=ik1-nkbd(iq)+1
            ik2=iq-ik1
            ImT(iw,ik1,ik2,ieh1,ieh2)=&
                 ImT(iw,ik1,ik2,ieh1,ieh2)+2.0d0*ImMW(i1)
            ImbT(iw,ik1,ik2,ieh1,ieh2)=&
                 ImbT(iw,ik1,ik2,ieh1,ieh2)+2.0d0*ImbMW(i1)
          end do
        end do
      end do
!$omp end parallel do
        end do
      end do

! check imaginary part of Tmatrix -------------------------------------
      open(65,file='T0.5n0.1_ImTpq_ee.txt',status='unknown')
      open(66,file='T0.5n0.1_ImTpq_eh.txt',status='unknown')
      open(68,file='T0.5n0.1_ImTpq_hh.txt',status='unknown')
      do ik1=0,nk,nk/25 ! p
        k1=xk(ik1)
        write(65,*) 'k=',k1*a1d
        write(66,*) 'k=',k1*a1d
        write(68,*) 'k=',k1*a1d
        do iw=nwnn1,nwg2(1,1)
          w=xw(iw)
          write(65,*) w*Esc/E1d,&
            -ImT(iw,ik1,ik1,1,1)*Esc/E1d*a1d**2+k1*a1d*15.0d0
        end do
        do iw=nwnn1,nwg2(2,1)
          w=xw(iw)
          write(66,*) w*Esc/E1d,&
            -ImT(iw,ik1,ik1,2,1)*Esc/E1d*a1d**2+k1*a1d*15.0d0
        end do
        do iw=nwnn1,nwg2(2,2)
          w=xw(iw)
          write(68,*) w*Esc/E1d,&
            -ImT(iw,ik1,ik1,2,2)*Esc/E1d*a1d+k1*a1d*15.0d0
        end do
      end do
      close(65)
      close(66)
      close(68)
! ---------------------------------------------------------------------

      return
      end subroutine TmatrixMW
!*********************************************************************








!*********************************************************************
      subroutine &
      ImselfTmatrixMW(selfT,ImT,ImbT,nwn1,nwe2,&
      nwnn1,nwee2,nk,nwg1,xw,divw,divk,fermi,bose,Akw)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0)
      complex(8),parameter :: im=(0.0d0,1.0d0)
      integer(4) :: &
        ik1,ik2,nk,nwg1(2),nwn1,nwe2,nwnn1,nwee2,&
        iw1,iw2,iw,ieh1,ieh2,ieh,ik1e,ik2h
      real(8) :: &
        xw(nwnn1:nwee2),divw,&
        w,ww,ImT(nwnn1:nwee2,-nk:nk,-nk:nk,2,2),divk,abc,def,&
        Akw(nwn1:nwe2,0:nk,2),fermi(nwn1:nwe2,2),bose(nwnn1:nwee2,3),&
        ImbT(nwnn1:nwee2,-nk:nk,-nk:nk,2,2)
      complex(8) :: selfT(nwn1:nwe2,0:nk,2)

      selfT=0.0d0
      do ieh2=1,2
        do ieh1=1,2
          if(ieh1.eq.ieh2) then
             if(ieh1.eq.1) ieh=1 ! e-e pair
             if(ieh1.eq.2) ieh=3 ! h-h pair
          else
             ieh=2 ! e-h pair
          end if

          do ik1=0,nk
!$omp parallel do private(iw1,iw2,ik2,iw,abc,def)
            do iw1=nwn1,nwg1(ieh1)
              do iw2=nwn1,nwg1(ieh2)
                iw=iw1+iw2
                abc=0.0d0
                def=0.0d0
                do ik2=-nk,nk
                  abc=abc+ImT(iw,ik1,ik2,ieh1,ieh2)*Akw(iw2,abs(ik2),ieh2)
                  def=def+ImbT(iw,ik1,ik2,ieh1,ieh2)*Akw(iw2,abs(ik2),ieh2)
                end do
                selfT(iw1,ik1,ieh1)=selfT(iw1,ik1,ieh1)&
                     +abc*fermi(iw2,ieh2)+def
              end do
            end do
!$omp end parallel do
          end do
        end do
      end do
      selfT=im*divw*divk*selfT/(2.0d0*pai)**2
      
      return
      end subroutine ImselfTmatrixMW
!*********************************************************************








! *********************************************************************
      subroutine &
      bosedis(bose,mu,nwnn1,nwee2,xw,beta,divw)
      implicit none
      integer(4) :: nwnn1,nwee2,ieh,ieh1,ieh2,iw,iw2,iw0
      real(8) :: &
        xw(nwnn1:nwee2),mu(2),bfnc,bose(nwnn1:nwee2,3),beta,mu0,&
        abc,S1,S2,S3,divw,w,a1,a2,a3,a4

    ! bose distribution function (no 0pt correction)
      ieh=0
      do ieh2=1,2
        do ieh1=ieh2,2
          ieh=ieh+1
          do iw=nwnn1,nwee2 ! w
            bose(iw,ieh)=bfnc(beta,xw(iw)-mu(ieh1)-mu(ieh2))
          end do
        end do
      end do


    ! bose distribution function with 0pt correction
      ieh=0
      do ieh2=1,2
        do ieh1=ieh2,2
          ieh=ieh+1
          mu0=mu(ieh1)+mu(ieh2) ! sum of mu
        ! set sign inversion point
          iw0=nwnn1
          do iw=nwnn1,nwee2 ! w
            if(xw(iw)-mu0.gt.0.0d0) then
              iw0=iw
              exit
            end if
          end do

          if(iw0.ne.nwnn1) then
          ! set close to mu point
            if(xw(iw0)-mu0.gt.divw*0.5d0) iw0=iw0-1

            a1=dexp(-beta*(xw(iw0-1)-mu0))-1.0d0
            a4=dexp(-beta*(xw(iw0+1)-mu0))-1.0d0
            S2=dlog(abs(a4/a1))/beta

          ! remove adjacent contribution
            bose(iw0,ieh)=&
            (S2-divw*(bose(iw0-1,ieh)+bose(iw0+1,ieh))/2.0d0)/divw
          end if

        end do
      end do

      end subroutine bosedis
! ********************************************************************





! ******************************************************************
      function Coulomb(q,x0ptx,y0ptx,n0ptx,divx,divy,&
        bessk,nxmax,nymax,phi2x,phi2y,midint)
      implicit none
      integer(4) :: n0ptx,nxmax,nymax,ix,ix1,ix2,iy,iy1,iy2,i
      real(8) :: &
        q,x0ptx(n0ptx),y0ptx(n0ptx),divx,divy,bessk(0:nxmax,0:nymax),&
        phi2x(0:nxmax),phi2y(0:nymax),midint(0:nymax),&
        bessk0,S1,bess0pt,def,Coulomb,xx,yy
        
      ! bessel function (no 0point correction)
        do iy=0,nymax
          yy=divy*dble(iy)
          do ix=0,nxmax
            xx=divx*dble(ix)
            bessk(ix,iy)=bessk0(dsqrt(xx**2+yy**2)*q)
          end do
        end do
      ! bessel function with 0point correction for x-integration
        S1=0.0d0
        do i=1,n0ptx
          S1=S1+bessk0(x0ptx(i)*q)*y0ptx(i)
        end do
        bessk(0,0)=2.0d0*(S1-divx*bessk(1,0)/2.0d0)/divx

      ! x-integration
        do iy=0,nymax
          def=0.0d0
          do ix2=0,nxmax
            do ix1=0,nxmax
              def=def+bessk(abs(ix2-ix1),iy)*&
              phi2x(ix1)*phi2x(ix2)
            end do
          end do
          midint(iy)=def*divx**2
        end do
      ! y-integration
        def=0.0d0
        do iy2=0,nymax
          do iy1=0,nymax
            def=def+midint(abs(iy2-iy1))*phi2y(iy1)*phi2y(iy2)
          end do
        end do
        
        Coulomb=2.0d0*def*divy**2
      
      return
      end function Coulomb
! **************************************************************************






! ***********************************************************************
      FUNCTION bessk0(x)
      REAL(8) :: bessk0,x
!CU    USES bessi0
      REAL(8) :: bessi0
      REAL(8) :: p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,&
      0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,&
      -0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
      if (x.le.2.0d0) then
        y=x*x/4.0d0
        bessk0=(-log(x/2.0d0)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*&
      (p6+y*p7))))))
      else
        y=(2.0d0/x)
        bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
      q7))))))
      endif
      return
      END FUNCTION bessk0
! ************************************************************************

! ***********************************************************************
      FUNCTION bessi0(x)
      REAL(8) :: bessi0,x
      REAL(8) :: ax
      REAL(8) ::  p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,&
      1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,&
      0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,&
      -0.1647633d-1,0.392377d-2/
      if (abs(x).lt.3.75d0) then
        y=(x/3.75d0)**2
        bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
      else
        ax=abs(x)
        y=3.75d0/ax
        bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*&
      (q7+y*(q8+y*q9))))))))
      endif
      return
      END FUNCTION bessi0
! *******************************************************************      
      
      


