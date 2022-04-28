      Program SBE1d
      implicit none
    ! me and mh="electron and hole mass", mr="reduced mass"
      real(8),parameter :: me=0.0665d0,mh=0.11d0,mr=me*mh/(me+mh)
    ! mass
      real(8),parameter :: mass(2)=(/me/mr,mh/mr/)
    ! E1d="1d exciton binding energy (E1d)",a1d="1d exciton Bohr radius (a3d)"
      real(8),parameter :: E1d=4.68569374122079d0,a1d=1.0d0/dsqrt(E1d)
    ! E1d*Esc="energy unit"
      real(8),parameter :: Esc=2.0d0

    ! initial condition read=1: read previdous result
      integer(4),parameter :: iread=0

    ! temp="temperature", num="carrier density"
      real(8),parameter :: temp=0.5d0*E1d/Esc,beta=1.0d0/temp,num=1.0d0/a1d
    ! const="SPP parameter"
      real(8),parameter :: const=4.0d0

    ! imaginary number and circular constant
      complex(8),parameter :: im=(0.0d0,1.0d0)
      real(8),parameter :: pai=dacos(-1.0d0)

! parameter --------------------------------------------------------
    ! gamma="\chi broadening", cnv="order of convergence"
      real(8),parameter :: gamma=1.0d0/14.0d0*E1d/Esc,cnv=1.0d-4
    ! iend="iteration limit"
      integer(4),parameter :: iend=500
    ! initial shift of chemical potential
      real(8),parameter :: dmu0=0.05d0
! ------------------------------------------------------------------


! set momentum ------------------------------------------------------
  ! integral parameter in main loop
    ! kmax="k cutoff", divk="k step size"
      real(8),parameter :: kmax=60.0d0,divk=0.1d0
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
      real(8),parameter :: divw=0.05d0,dw=75.0d0
    ! wemax="highest energy of free electron"
      real(8), parameter :: wemax=0.5d0*kmax**2/mass(1)
    ! whmax="highest energy of free hole energy"
      real(8), parameter :: whmax=0.5d0*kmax**2/mass(2)
  ! single particle
    ! nwe2="high energy edge point of electron energy"
      integer(4),parameter :: nwe2=int((dw+wemax)/divw)
    ! nwh2="high energy edge point of hole energy"
      integer(4),parameter :: nwh2=int((dw+whmax)/divw)

  ! two particle
    ! nwnn1="low energy edge point of pair energy"
      integer(4),parameter :: nwnn1=-int(dw/divw)*2
    ! nweh2 ="high energy edge point of pair energy"
      integer(4),parameter :: nweh2=nwe2+nwh2
    ! omega points
      real(8) :: xw(nwnn1:nweh2)
! ------------------------------------------------------------------


! effective interaction ---------------------------------------------------
    ! Wsc="screened interaction"
      real(8) :: Wsc(-2*nk:2*nk)
! --------------------------------------------------------------------------

        
! self energy ---------------------------------------------------------    
    ! "exchange term or SHF self energy"
      real(8) :: self0(0:nk,2)
! --------------------------------------------------------------------


! single particle -------------------------------------------------------
    ! engkf and engkSHF="free and SHF energy",engkSHFeh="sum of e-h engkSHF"
      real(8) :: engkf(0:nk,2),engkSHF(0:nk,2),engkSHFeh(1:2*nk+1)
    ! fermi0="fermi and distribution function"
      real(8) :: fermi0(0:nk,2)
    ! numeh and numehold ="carrier density", numeh0="QP carrier density"
      real(8) :: numeh(2),numehold(2),numeh0(2),mu(2)
! --------------------------------------------------------------------


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


! generalized diagonalization -----------------------------------------
    ! work arrays for DGEEV
      character(1) JOBVL,JOBVR
        data JOBVL/'V'/,JOBVR/'V'/
      integer(4),parameter :: LWORK=4*(2*nk+1)
      real(8) :: H0(2*nk+1,2*nk+1),WORK(LWORK),&
        pblock0(2*nk+1),VL(2*nk+1,2*nk+1),VR(2*nk+1,2*nk+1),&
        WR(2*nk+1),WI(2*nk+1),VL2(2*nk+1),VR2(2*nk+1)
! ---------------------------------------------------------------------

! optical response ---------------------------------------------------
    ! weight of pole
      real(8) :: weight(2*nk+1),weight0(2*nk+1)
! -------------------------------------------------------------------


      integer(4) :: &
        i,it,j,iw,iw1,iw2,iw3,iw4,i1,i2,i3,j1,j2,j3,l,ieh,ieh1,ieh2,&
        ik1,ik2,ik3,ik,ip1,ip2,ip,iq,iq0,nwmax,icent,icent0,ite,ite2,&
        icount,iave,iabc,ite3,INFO,ix,iy,ieigen
      real(8) :: &
        k1,k2,absorb,selfCH,ktf2,kdh2,kappa,ffnc,w,ww,www,w0,sum,bfnc,&
        dense,p,q,t,q0,CH,CH2,abc,def,ghi,jkl,S1,Img2,a,a2,&
        screening,cffInt,xx,yy,Coulomb
      complex(8) :: zabc,zdef




      write(6,*) 'T (E1d)',Esc/beta/E1d
      write(6,*) 'carrier density (1/a1d)',num*a1d
      write(6,*)
      write(6,*) 'setup'
      write(6,*) 'kmax=',kmax*a1d
      write(6,*) 'nk=',nk
      write(6,*) 'dw=',dw*Esc/E1d
      write(6,*) 'divw=',divw*Esc/E1d
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
      do iw=nwnn1,nweh2
         xw(iw)=divw*dble(iw)
      end do
! -----------------------------------------------------------


! read previdous results ----------------------------------------
      if(iread.eq.1) then
        read(5,*) mu(1),mu(2)
        do ieigen=1,2*nk+1
          read(5,*) WR(ieigen),weight(ieigen),weight0(ieigen),&
               engkSHFeh(ieigen)
        end do
        goto 460
      end if 
! ------------------------------------------------------------------


! set real space integral parameter -----------------------
      call gauleg(0.0d0,divx,x0ptx,y0ptx,n0ptx)
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


! Screened HF --------------------------------------------------------
    ! SHF mu and selfenergy
      do i=1,2 ! initial mu (classical)
        mu(i)=dlog(num*dsqrt(0.5d0*pai*beta/mass(i)))/beta
      end do

      call SHF(mu,beta,num,nk,engkSHF,engkf,fermi0,dmu0,&
        yk,Wsc,Vtot,V0pt,VCH,n0ptk,x0ptk,y0ptk,divk,xk,&
        nCH,xCH,yCH,self0,cnv,const,iend,E1d,a1d,kappa,selfCH)
! ----------------------------------------------------------------------


! output SHF thermodynamic results ------------------------------
      write(6,*) 'SHF results'
      write(6,*) 'total mu (E1d)'
      write(6,*) Esc*(mu(1)+mu(2))/E1d
      write(6,*) 'BGR (E1d)'
      write(6,*) (self0(0,1)+self0(0,2))*Esc/E1d
      write(6,*) 'Coulomb Hole (E1d)'
      write(6,*) selfCH*Esc/E1d

      open(47,file='T0.5n1_SBE_thermodynamic_results.txt',status='unknown')
      write(47,*) 'E chemical potential (E1d)'
      write(47,*) num*a1d,mu(1)*Esc/E1d
      write(47,*) 'H chemical potential (E1d)'
      write(47,*) num*a1d,mu(2)*Esc/E1d
      write(47,*) 'total chemical potential (E1d)'
      write(47,*) num*a1d,(mu(1)+mu(2))*Esc/E1d
      write(47,*)

      write(47,*) 'screening parameter'
      write(47,*) num*a1d,kappa
      write(47,*) 'BGR (E1d)'
      write(47,*) num*a1d,(self0(0,1)+self0(0,2))*Esc/E1d
      write(47,*) 

      write(47,*) 'Coulomb hole self energy'
      write(47,*) num*a1d,selfCH*Esc/E1d
      write(47,*) 'E screened exchange self energy'
      write(47,*) num*a1d,(self0(0,1)-selfCH)*Esc/E1d
      write(47,*) 'H screened exchange self energy'
      write(47,*) num*a1d,(self0(0,2)-selfCH)*Esc/E1d
      write(47,*) 'total screened exchange self energy'
      write(47,*) num*a1d,(self0(0,1)+self0(0,2)-2.0d0*selfCH)*Esc/E1d
      write(47,*)
! -----------------------------------------------------------------


! exact diagonalization -------------------------------------------
    ! sum of e-h pair energy
      do ik1=-nk,nk
        engkSHFeh(ik1+nk+1)=engkSHF(abs(ik1),1)+engkSHF(abs(ik1),2)
      end do

    ! pauli blocking factor
      do ik1=-nk,nk
        pblock0(ik1+nk+1)=1.0d0-fermi0(abs(ik1),1)-fermi0(abs(ik1),2)
      end do

    ! kinetic term
      H0=0.0d0
      do ieigen=1,2*nk+1
        H0(ieigen,ieigen)=engkSHFeh(ieigen)
      end do

    ! Coulomb term
      do ik2=-nk,nk
        do ik1=-nk,nk
          H0(ik1+nk+1,ik2+nk+1)=H0(ik1+nk+1,ik2+nk+1)&
          -dsign(1.0d0,pblock0(ik1+nk+1))*dsqrt(abs(pblock0(ik1+nk+1)))*&
          Wsc(ik1-ik2)*dsqrt(abs(pblock0(ik2+nk+1)))/(2.0d0*pai)*divk
        end do
      end do

    ! exact diagonalization
      call DGEEV(JOBVL,JOBVR,2*nk+1,H0,2*nk+1,WR,WI,VL,2*nk+1,VR,&
      2*nk+1,WORK,LWORK,INFO)
      if(INFO.ne.0) then
        write(6,*) 'DGEEV error, INFO=',INFO
      end if
      
    ! check complex eigen value
      do ieigen=1,2*nk+1
        if(abs(WI(ieigen)).gt.1.0d-10) then
          write(6,*) 'complex eigen value',WI(ieigen)*Esc/E1d
        end if
      end do

    ! weight of pole
      weight0=0.0d0
      weight=0.0d0
      do ieigen=1,2*nk+1
        do ik2=1,2*nk+1
          do ik1=1,2*nk+1
            weight(ieigen)=weight(ieigen)+dsqrt(abs(pblock0(ik2)))*&
            VR(ik2,ieigen)*VL(ik1,ieigen)*dsqrt(abs(pblock0(ik1)))*&
            dsign(1.0d0,pblock0(ik1))
          end do
        end do
        weight0(ieigen)=pblock0(ieigen)
      end do
      weight=weight*divk/(2.0d0*pai)
      weight0=weight0*divk/(2.0d0*pai)

    ! sort by enery level
      do ik2=1,2*nk
        do ik1=2*nk,ik2,-1
        if(WR(ik1).gt.WR(ik1+1)) then
          abc=WR(ik1+1)
          def=weight(ik1+1)
          do i=1,2*nk+1
            VL2(i)=VL(i,ik1+1)
            VR2(i)=VR(i,ik1+1)
          end do
          WR(ik1+1)=WR(ik1)
          WR(ik1)=abc
          weight(ik1+1)=weight(ik1)
          weight(ik1)=def
          do i=1,2*nk+1
            VL(i,ik1+1)=VL(i,ik1)
            VR(i,ik1+1)=VR(i,ik1)
          end do
          do i=1,2*nk+1
            VL(i,ik1)=VL2(i)
            VR(i,ik1)=VR2(i)
          end do
        end if
        end do
      end do

! output SBE results -------------------------------------------------
      open(89,file='T0.5n1_SBE_results.txt',status='unknown')
      write(89,*) mu(1),mu(2)
      do ieigen=1,2*nk+1
        write(89,*) WR(ieigen),weight(ieigen),weight0(ieigen),engkSHFeh(ieigen)
      end do
      write(6,*) 'SBE lowest energy level'
      write(6,*) WR(1)*Esc/E1d
      write(47,*) 'SBE lowest energy level'
      write(47,*) num*a1d,WR(1)*Esc/E1d
! -------------------------------------------------------------------


! optical response --------------------------------------------------
460   open(222,file='T0.5n1_SBE_Imx.txt',status='unknown')
      open(125,file='T0.5n1_SBE_Rex.txt',status='unknown')
      open(123,file='T0.5n1_SBE_PL.txt',status='unknown')
      open(32,file='T0.5n1_SBE_Imx0.txt',status='unknown')
      open(132,file='T0.5n1_SBE_Rex0.txt',status='unknown')
      open(35,file='T0.5n1_SBE_PL0.txt',status='unknown')
      open(42,file='T0.5n1_SBE_Imx_per_Imx0.txt',status='unknown')
      open(142,file='T0.5n1_SBE_Rex_per_Rex0.txt',status='unknown')

      do iw=nwnn1,nweh2
        if(abs(xw(iw))*Esc/E1d.lt.20.0d0) then
        zabc=0.0d0
        zdef=0.0d0
        abc=0.0d0
        def=0.0d0
        w=xw(iw)
        do ieigen=1,2*nk+1
          zabc=zabc+weight(ieigen)/((w-WR(ieigen))+im*gamma)
          zdef=zdef+weight0(ieigen)/(w-engkSHFeh(ieigen)+im*gamma)

          abc=abc+bfnc(beta,WR(ieigen)-mu(1)-mu(2))*weight(ieigen)*&
               dimag(1.0d0/(w-WR(ieigen)+im*gamma))
          def=def+bfnc(beta,engkSHFeh(ieigen)-mu(1)-mu(2))*weight0(ieigen)*&
               dimag(1.0d0/(w-engkSHFeh(ieigen)+im*gamma))
        end do
        zabc=-2.0d0*zabc
        zdef=-2.0d0*zdef
        abc=-2.0d0*abc
        def=-2.0d0*def
        write(222,*) w*Esc/E1d,dimag(zabc)*a1d*E1d/Esc
        write(125,*) w*Esc/E1d,dreal(zabc)*a1d*E1d/Esc
        write(123,*) w*Esc/E1d,abc*a1d*E1d/Esc
        write(32,*) w*Esc/E1d,dimag(zdef)*a1d*E1d/Esc
        write(132,*) w*Esc/E1d,dreal(zdef)*a1d*E1d/Esc
        write(35,*) w*Esc/E1d,def*a1d*E1d/Esc
        write(42,*) w*Esc/E1d,dimag(zabc)/dimag(zdef)
        write(142,*) w*Esc/E1d,dreal(zabc)/dreal(zdef)
        end if
      end do
! ------------------------------------------------------------

      
! pair susceptibility (gamma=0) ----------------------------------
      open(56,file='T0.5n1_SBE_pair_susceptibility.txt',status='unknown') 
      abc=0.0d0
      def=0.0d0
      do ieigen=1,2*nk+1
        abc=abc+weight(ieigen)/(mu(1)+mu(2)-WR(ieigen))
        def=def+weight0(ieigen)/(mu(1)+mu(2)-engkSHFeh(ieigen))
      end do
      abc=-2.0d0*abc
      def=-2.0d0*def
      write(56,*) "Rex at w=mu (1/(a1d*E1d))" 
      write(56,*) num*a1d,abc*a1d*E1d/Esc
      write(56,*) "Rex0 at w=mu (1/(a1d*E1d))"
      write(56,*) num*a1d,def*a1d*E1d/Esc
      write(56,*) "Rex/Rex0 at w=mu" 
      write(56,*) num*a1d,abc/def
      write(56,*) "Rex-Rex0 at w=mu (1/(a1d*E1d))"
      write(56,*) num*a1d,(abc-def)*a1d*E1d/Esc
      close(56)
! -------------------------------------------------------------

      stop
      end program SBE1d
! *******************************************************************






      
      
      
! ***********************************************************************
      subroutine SHF(mu,beta,num,nk,engkSHF,engkf,fermi0,dmu0,&
        yk,Wsc,Vtot,V0pt,VCH,n0ptk,x0ptk,y0ptk,divk,xk,&
        nCH,xCH,yCH,self0,cnv,const,iend,E1d,a1d,kappa,selfCH)
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

    ! fermi distribution function
      do ieh=1,2 
        do ik=0,nk ! k
          fermi0(ik,ieh)=ffnc(beta,engkSHF(ik,ieh),mu(ieh))
        end do
      end do

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
      real(8) :: x1,xk,x(n),w(n)
      real(8),PARAMETER :: EPS=3.d-14
      INTEGER(4) :: i,j,m
      real(8) :: p1,p2,p3,pp,xl,xm,z,z1
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
      real(8) :: bessk0,x
!CU    USES bessi0
      real(8) :: bessi0
      real(8) :: p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y
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
      real(8) :: bessi0,x
      real(8) :: ax
      real(8) ::  p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
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
      
      


