      Program Tmatrix2d
      implicit none
    ! me and mh="electron and hole mass", mr="reduced mass"        
      real(8),parameter :: me=0.0665d0,mh=0.11d0,mr=me*mh/(me+mh)
    ! mass
      real(8),parameter :: mass(2)=(/me/mr,mh/mr/)
    ! sum of mass
      real(8),parameter :: M(3)=(/2.0d0*me/mr,(me+mh)/mr,2.0d0*mh/mr/)
    ! a2d="2d exciton Bohr radius per a3d", E2d="2d exciton binding energy per E3d"
      real(8),parameter :: a2d=1.0d0/2.0d0,E2d=4.0d0

    ! me0,mh0,mr="3d GaAs e,h,reduced mass"
      real(8),parameter :: me0=0.0665d0,mh0=0.457d0,mr0=me0*mh0/(me0+mh0)

    ! E3d0="3d GaAs binding energy (eV)"
      real(8),parameter :: E3d0=4.2d-3

    ! E2d*Esc="energy unit",Eg="band-gap(1.5 eV)"
      real(8),parameter :: Esc=2.0d0,Eg=1.5d0/(E3d0*mr/mr0)/Esc

    ! initial condition read=1: read previdous result
      integer(4),parameter :: iread=0

    ! temp="temperature", num="carrier density"
      real(8),parameter :: temp=0.5d0*E2d/Esc,beta=1.0d0/temp,num=5.0d-4/a2d**2
    ! const="SPP parameter"
      real(8),parameter :: const=4.0d0

    ! imaginary number and circular constant
      complex(8),parameter :: im=(0.0d0,1.0d0)
      real(8),parameter :: pai=dacos(-1.0d0)

! parameter --------------------------------------------------------
    ! gamma="\chi broadening", cnv="order of convergence"
      real(8),parameter :: gamma=0.05d0*E2d/Esc,cnv=1.0d-4
    ! iend="iteration limit"
      integer(4),parameter :: iend=500
    ! initial shift of chemical potential
      real(8),parameter :: dmu0=0.1d0/Esc
! ------------------------------------------------------------------


! set momentum ------------------------------------------------------
  ! integral parameter in main loop
    ! kmax="k cutoff", divk="k step size"
      real(8),parameter :: kmax=30.0d0/a2d,divk=0.05d0/a2d
    ! nk="number of k point"
      integer(4),parameter :: nk=int(kmax/divk)
    ! xk and xyk="point and weight"
      real(8) :: xk(0:nk),xyk(0:nk),xyksqrt(0:nk)

  ! integral parameter in Coulomb hole
    ! kCH="k cutoff", divCH="k step size"
      real(8),parameter :: kCH=500.0d0/a2d,divCH=1.0d-4/a2d
    ! nCH="number of k point"
      integer(4),parameter :: nCH=int(kCH/divCH)
    ! xCH and yCH="point and weight"
      real(8) :: xCH(0:nCH),yCH(0:nCH)

  ! integral parameter for 0 point correction
    ! n0ptk="number of fine k points"
      integer(4),parameter :: n0ptk=500
    ! x0ptk and y0ptk="point and weight"
      real(8) :: x0ptk(n0ptk),y0ptk(n0ptk)
! ----------------------------------------------------------------------


! set omega -------------------------------------------------------
    ! divw="step size of omega", dw="extarnal part"
      real(8),parameter :: divw=1.0d0/40.0d0*E2d/Esc,dw=37.5d0*E2d/Esc
    ! wemax and whmax="highest energy of free electron and hole energy"
      real(8), parameter :: wemax=0.5d0*kmax**2/mass(1),whmax=0.5d0*kmax**2/mass(2)

  ! single particle
    ! nwe2 and nwh2="high energy edge point of electron and hole energy"
      integer(4),parameter :: nwe2=int((dw+wemax)/divw),nwh2=int((dw+whmax)/divw)

  ! two particle
    ! nwnn1="low energy edge point of pair energy"
      integer(4),parameter :: nwnn1=-int(dw/divw)*2
    ! nweh2 ="high energy edge point of pair energy"
      integer(4),parameter :: nweh2=nwe2+nwh2
    ! omega points
      real(8) :: xw(nwnn1:nweh2)
! ------------------------------------------------------------------


! set angle average -------------------------------------------
    ! nt0="number of angle point"
      integer(4),parameter :: nt0=60
    ! divt="step size of angle"
      real(8),parameter :: divt=pai/dble(nt0)
    ! xt0 and yt0="point and weight"
      real(8) :: xt0(0:nt0),yt0(0:nt0)
! ------------------------------------------------------------------


! effective interaction --------------------------------------------------
    ! Wscexp="screened interaction"
      real(8) :: Wscexp(0:nk,0:nk)
! --------------------------------------------------------------------------


! self energy --------------------------------------------------------
    ! "exchange term or SHF self energy"
      real(8) :: self0(0:nk,2)
! --------------------------------------------------------------------


! single particle -------------------------------------------------------
    ! engkf and engkSHF="free and SHF energy",engkSHFeh="sum of e-h engkSHF"
      real(8) :: engkf(0:nk,2),engkSHF(0:nk,2),engkSHFeh(nk+1)
    ! fermi="fermi distribution function"
      real(8) :: fermi0(0:nk,2)
    ! numeh and numehold ="carrier density", numeh0="QP carrier density"
      real(8) :: numeh(2),numehold(2),numeh0(2),mu(2)
! -----------------------------------------------------------------------


! generalized diagonalization ---------------------------------------
    ! work array for DGEEV
      character(1) JOBVL,JOBVR
      data JOBVL/'V'/,JOBVR/'V'/
      integer(4),parameter :: LWORK=4*(nk+1)
      real(8) :: H0(nk+1,nk+1),WORK(LWORK),&
        VL(nk+1,nk+1),VR(nk+1,nk+1),pblock0(nk+1),&
        VL2(nk+1),VR2(nk+1),WR(nk+1),WI(nk+1)
! ----------------------------------------------------------------------


! optical response -------------------------------------------------------
    ! weight of pole
      real(8) :: weight0(nk+1),weight(nk+1)
    ! dipole matrix element
      real(8) :: dk(0:nk)
! -----------------------------------------------------------------------


      integer(4) :: i,j,iw,iw1,iw2,l,ieh,ieh1,ieh2,&
           ik1,ik2,ik,ip,iq,it,it0,nwmax,icount,INFO,ieigen,icheck
      real(8) :: k1,k2,selfCH,ksc,ffnc,bfnc,w,ww,p,q,t0,&
           CH,abc,def,ghi,jkl,aion,screening,cffInt
      complex(8) :: zabc,zdef



      write(6,*) 'T (E2d)',Esc/beta/E2d
      write(6,*) 'carrier density (1/a2d^2)',num*a2d**2
      write(6,*)
      write(6,*) 'setup'
      write(6,*) 'kmax=',kmax*a2d
      write(6,*) 'nk=',nk
      write(6,*) 'dw=',dw*Esc/E2d
      write(6,*) 'divw=',divw*Esc/E2d
      write(6,*) 'iread=',iread
      write(6,*)








! set momentum ----------------------------------------
  ! integral parameter in main loop
    ! xk="k points"
      do ik=0,nk ! k
         xk(ik)=divk*dble(ik)
      end do
    ! xyk="weight of radial integral of k"
      xyk(0)=0.0d0
      do ik=1,nk
         xyk(ik)=divk**2*dble(ik)
      end do
    ! xyksqrt=dsqrt(xyk)
      xyksqrt(0)=0.0d0
      do ik=1,nk
         xyksqrt(ik)=dsqrt(xyk(ik))
      end do

  ! integral parameter in Coulomb hole
    ! xk="k points"
      do ik=0,nCH
         xCH(ik)=divCH*dble(ik)
      end do
    ! xyk="weight of radial integral of k"
      yCH(0)=0.5d0*divCH
      do ik=1,nCH
         yCH(ik)=divCH
      end do

  ! integral parameter in 0pt correction
      call gauleg(0.0d0,divk,x0ptk,y0ptk,n0ptk)
! -----------------------------------------------------------


! set omega -----------------------------------------------
    ! xw="omega points"
      do iw=nwnn1,nweh2
         xw(iw)=divw*dble(iw)
      end do
! -----------------------------------------------------------


! set angle average paramter (interaction) ------------
    ! Cos expansion
      call gauleg(0.0d0,pai,xt0,yt0,nt0) ! numerical recipes
! -----------------------------------------------------------------


! free particle energy ---------------------------------------------
      do ieh=1,2 
        do ik1=0,nk ! k1
          engkf(ik1,ieh)=0.5d0*xk(ik1)**2/mass(ieh)
        end do
      end do
! --------------------------------------------------------------


! read previdous results ----------------------------------------
      if(iread.eq.1) then
        read(5,*) mu(1),mu(2)
        do ieigen=1,nk+1
          read(5,*) WR(ieigen),weight(ieigen),weight0(ieigen),&
               engkSHFeh(ieigen)
        end do
        goto 460
      end if 
! --------------------------------------------------------------


! Screened HF ---------------------------------------------------------
      ! SHF mu and selfenergy 
      ! set initial mu
        do i=1,2
          mu(i)=dlog(dexp(beta*pai*num/mass(i))-1.0d0)/beta ! free mu
        end do
      ! SHF
        call SHF(mu,beta,num,nk,engkSHF,engkf,fermi0,dmu0,&
           xyk,nt0,yt0,Wscexp,n0ptk,x0ptk,y0ptk,divk,xk,xt0,&
           nCH,xCH,yCH,self0,cnv,const,iend,ksc,selfCH)
! --------------------------------------------------------------------------


! output SHF results -----------------------------------------------------
      write(6,*) 'SHF results'
      write(6,*) 'total mu (E2d)'
      write(6,*) Esc*(mu(1)+mu(2))/E2d
      write(6,*) 'BGR (E2d)'
      write(6,*) (self0(0,1)+self0(0,2))*Esc/E2d
      write(6,*) 'Coulomb Hole (E2d)'
      write(6,*) selfCH*Esc/E2d

      open(47,file='T0.5n5d-4_SBE_thermodynamic_results.txt',status='unknown')
      write(47,*) 'E chemical potential (E2d)'
      write(47,*) num*a2d**2,mu(1)*Esc/E2d
      write(47,*) 'H chemical potential (E2d)'
      write(47,*) num*a2d**2,mu(2)*Esc/E2d
      write(47,*) 'total chemical potential (E2d)'
      write(47,*) num*a2d**2,(mu(1)+mu(2))*Esc/E2d
      write(47,*)
      
      write(47,*) 'screening wave number (1/a2d)'
      write(47,*) num*a2d**2,ksc*a2d
      write(47,*) 'BGR (E2d)'
      write(47,*) num*a2d**2,(self0(0,1)+self0(0,2))*Esc/E2d
      
      write(47,*) 'Coulomb hole self energy'
      write(47,*) num*a2d**2,selfCH*Esc/E2d
      write(47,*) 'E screened exchange self energy'
      write(47,*) num*a2d**2,(self0(0,1)-selfCH)*Esc/E2d
      write(47,*) 'H screened exchange self energy'
      write(47,*) num*a2d**2,(self0(0,2)-selfCH)*Esc/E2d
      write(47,*) 'total screened exchange self energy'
      write(47,*) num*a2d**2,(self0(0,1)+self0(0,2)-2.0d0*selfCH)*Esc/E2d
! -----------------------------------------------------------------


! exact diagonalization -------------------------------------------
    ! sum of e-h pair energy
      do ik1=0,nk
         engkSHFeh(ik1+1)=engkSHF(ik1,1)+engkSHF(ik1,2)
      end do

    ! pauli blocking factor
      do ik1=0,nk
        pblock0(ik1+1)=1.0d0-fermi0(ik1,1)-fermi0(ik1,2)
      end do

    ! dipole matrix element
      do ik1=0,nk
        dk(ik1)=Eg/(Eg+0.5d0*xk(ik1)**2)
      end do
 
   ! kinetic term
      H0=0.0d0
      do ieigen=1,nk+1
        H0(ieigen,ieigen)=engkSHFeh(ieigen)
      end do

    ! Coulomb term
      do ik2=0,nk
        do ik1=0,nk
          H0(ik1+1,ik2+1)=H0(ik1+1,ik2+1)&
          -dsign(1.0d0,pblock0(ik1+1))*dsqrt(abs(pblock0(ik1+1)))*xyksqrt(ik1)*&
          Wscexp(ik1,ik2)*dsqrt(abs(pblock0(ik2+1)))*xyksqrt(ik2)/(2.0d0*pai)
        end do
      end do
      
    ! exact diagonalization
      call DGEEV(JOBVL,JOBVR,nk+1,H0,nk+1,WR,WI,VL,nk+1,VR,&
      nk+1,WORK,LWORK,INFO)
      if(INFO.ne.0) then
        write(6,*) 'DGEEV error, INFO=',INFO
      end if

    ! check cmplex eigen value
      do ieigen=1,nk+1
        if(abs(WI(ieigen)).ge.1.0d-10) then
          write(6,*) 'complex eigen value',WI(ieigen)*Esc/E2d
        end if
      end do
      
    ! weight of pole
      weight=0.0d0
      weight0=0.0d0
      do ieigen=1,nk+1
        do ik2=0,nk
          do ik1=0,nk
            weight(ieigen)=weight(ieigen)+&
            dk(ik2)*xyksqrt(ik2)*VR(ik2+1,ieigen)*dsqrt(abs(pblock0(ik2+1)))*&
            dk(ik1)*xyksqrt(ik1)*VL(ik1+1,ieigen)*dsqrt(abs(pblock0(ik1+1)))*&
            dsign(1.0d0,pblock0(ik1+1))
          end do
        end do
        weight0(ieigen)=weight0(ieigen)+&
        dk(ieigen-1)**2*xyk(ieigen-1)*pblock0(ieigen)
      end do
      weight=weight/(2.0d0*pai)
      weight0=weight0/(2.0d0*pai)

    ! sort by energy level
      do ik2=1,nk
        do ik1=nk,ik2,-1
        if(WR(ik1).gt.WR(ik1+1)) then
          abc=WR(ik1+1)
          def=weight(ik1+1)
          do i=1,nk+1
            VL2(i)=VL(i,ik1+1)
            VR2(i)=VR(i,ik1+1)
          end do
          WR(ik1+1)=WR(ik1)
          WR(ik1)=abc
          weight(ik1+1)=weight(ik1)
          weight(ik1)=def
          do i=1,nk+1
            VL(i,ik1+1)=VL(i,ik1)
            VR(i,ik1+1)=VR(i,ik1)
          end do
          do i=1,nk+1
            VL(i,ik1)=VL2(i)
            VR(i,ik1)=VR2(i)
          end do
          
        end if
        end do
      end do

! output SBE results -----------------------------------------------
      open(89,file='T0.5n5d-4_SBE_results.txt',status='unknown')
      write(89,*) mu(1),mu(2)
      do ieigen=1,nk+1
        write(89,*) WR(ieigen),weight(ieigen),weight0(ieigen),engkSHFeh(ieigen)
      end do
      close(89)
      write(6,*) 'SBE lowest energy level'
      write(6,*) WR(1)*Esc/E2d
      write(47,*) 'SBE lowest energy level'
      write(47,*) num*a2d**2,WR(1)*Esc/E2d
! -----------------------------------------------------------------      
      

! optical response -------------------------------------------------
460   open(222,file='T0.5n5d-4_SBE_Imx.txt',status='unknown')
      open(125,file='T0.5n5d-4_SBE_Rex.txt',status='unknown')
      open(123,file='T0.5n5d-4_SBE_PL.txt',status='unknown')
      open(32,file='T0.5n5d-4_SBE_Imx0.txt',status='unknown')
      open(132,file='T0.5n5d-4_SBE_Rex0.txt',status='unknown')
      open(35,file='T0.5n5d-4_SBE_PL0.txt',status='unknown')
      open(42,file='T0.5n5d-4_SBE_Imx_per_Imx0.txt',status='unknown')
      open(142,file='T0.5n5d-4_SBE_Rex_per_Rex0.txt',status='unknown')

      do iw=nwnn1,nweh2
        zabc=0.0d0
        zdef=0.0d0
        abc=0.0d0
        def=0.0d0
        w=xw(iw)
        if(abs(w*Esc/E2d).lt.20.0d0) then
        do ieigen=1,nk+1
          zabc=zabc+weight(ieigen)/((w-WR(ieigen))+im*gamma)
          zdef=zdef+weight0(ieigen)/&
          (w-engkSHFeh(ieigen-1)+im*gamma)
          abc=abc+bfnc(beta,WR(ieigen)-mu(1)-mu(2))*weight(ieigen)*&
               dimag(1.0d0/((w-WR(ieigen))+im*gamma))
          def=def+bfnc(beta,engkSHFeh(ieigen)-mu(1)-mu(2))*weight0(ieigen)*&
               dimag(1.0d0/((w-WR(ieigen))+im*gamma))
        end do
        zabc=-2.0d0*zabc
        zdef=-2.0d0*zdef
        abc=-2.0d0*abc
        def=-2.0d0*def
        write(222,*) w*Esc/E2d,dimag(zabc)*a2d**2*E2d/Esc
        write(125,*) w*Esc/E2d,dreal(zabc)*a2d**2*E2d/Esc
        write(123,*) w*Esc/E2d,abc*a2d**2*E2d/Esc
        write(32,*) w*Esc/E2d,dimag(zdef)*a2d**2*E2d/Esc
        write(132,*) w*Esc/E2d,dreal(zdef)*a2d**2*E2d/Esc
        write(35,*) w*Esc/E2d,def*a2d**2*E2d/Esc
        write(42,*) w*Esc/E2d,dimag(zabc)/dimag(zdef)
        write(142,*) w*Esc/E2d,dreal(zabc)/dreal(zdef)
        end if
      end do
! -----------------------------------------------------------------

! pair susceptibility (gamma=0) -----------------------------------
      open(56,file='T0.5n5d-4_SBE_pair_susceptibility.txt',status='unknown')
      abc=0.0d0
      def=0.0d0
      do ieigen=1,nk+1
        abc=abc+weight(ieigen)/((mu(1)+mu(2)-WR(ieigen)))
        def=def+weight0(ieigen)/((mu(1)+mu(2)-engkSHFeh(ieigen)))
      end do
      abc=-2.0d0*abc
      def=-2.0d0*def
      write(56,*) "Rex at w=mu (1/(a2d**2*E2d))" 
      write(56,*) num*a2d**2,abc*a2d**2*E2d/Esc
      write(56,*) "Rex0 at w=mu (1/(a2d**2*E2d))"
      write(56,*) num*a2d**2,def*a2d**2*E2d/Esc
      write(56,*) "Rex/Rex0 at w=mu" 
      write(56,*) num*a2d**2,abc/def
      write(56,*) "Rex-Rex0 at w=mu (1/(a2d**2*E2d))"
      write(56,*) num*a2d**2,(abc-def)*a2d**2*E2d/Esc
      close(56)
! -------------------------------------------------------------

      stop

      
      end program Tmatrix2d
! ************************************************************************













      
      
! ***********************************************************************
      subroutine SHF(mu,beta,num,nk,engkSHF,engkf,fermi0,dmu0,&
           xyk,nt0,yt0,Wscexp,n0ptk,x0ptk,y0ptk,divk,xk,xt0,&
           nCH,xCH,yCH,self0,cnv,const,iend,ksc,selfCH)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0,E2d=4.0d0,a2d=0.5d0
      integer(4) :: i,ieh,ik,ik1,ik2,nk,icount,nt0,n0ptk,nCH,iend
      real(8) :: mu(2),engkSHF(0:nk,2),engkf(0:nk,2),fermi0(0:nk,2),&
           xyk(0:nk),Wscexp(0:nk,0:nk),yt0(nt0),&
           x0ptk(n0ptk),y0ptk(n0ptk),xk(0:nk),xt0(nt0),&
           xCH(0:nCH),yCH(0:nCH),self0(0:nk,2),numeh(2),numehold(2)
      real(8) :: beta,num,const,dmu0,f,ksc,cffInt,divk,selfCH,CH,&
           abc,k1,cnv,screening,ffnc

    ! initial energy and fermi distribution
      do ieh=1,2 
        do ik=0,nk ! k
          engkSHF(ik,ieh)=engkf(ik,ieh)
          fermi0(ik,ieh)=ffnc(beta,engkSHF(ik,ieh),mu(ieh))
        end do
      end do

      write(6,*) 'SHF e-h particle density (1/a2d^2)'
      icount=0
! SHF iteration -------------------------------------------------
      do 600
        icount=icount+1
      ! set mueh
        do ieh=1,2 
          call shiftmuSHF(mu,ieh,dmu0,num,beta,nk,xyk,engkSHF)
          do ik=0,nk ! k
            fermi0(ik,ieh)=ffnc(beta,engkSHF(ik,ieh),mu(ieh))
          end do
        end do

      ! effective interaction
        ksc=screening(beta,engkSHF,nk,xyk,fermi0)
        cffInt=const*ksc/32.0d0/num/pai
        call Intexp0(nk,Wscexp,yt0,ksc,cffInt,nt0,&
        n0ptk,x0ptk,y0ptk,divk,xk,xt0)        

      ! Coulomb hole enery
        selfCH=CH(ksc,nCH,xCH,yCH,cffInt)

      ! exchange energy
        self0=0.0d0
        do ieh=1,2 
           do ik2=0,nk ! k2
              abc=xyk(ik2)*fermi0(ik2,ieh)
              do ik1=0,nk ! k1
                 self0(ik1,ieh)=self0(ik1,ieh)+Wscexp(ik1,ik2)*abc
              end do
           end do
        end do
        self0=-self0/(2.0d0*pai)+selfCH

      ! total energy
        do ieh=1,2 
          do ik=0,nk ! k
            k1=xk(ik)
            engkSHF(ik,ieh)=engkf(ik,ieh)+self0(ik,ieh)
          end do
        end do

      ! total density
        numeh=0.0d0
        do ieh=1,2 
           abc=0.0d0
           do ik1=0,nk ! k1
              abc=abc+xyk(ik1)*ffnc(beta,engkSHF(ik1,ieh),mu(ieh))
           end do
           numeh(ieh)=abc
        end do
        numeh=2.0d0*numeh/(2.0d0*pai)

        write(6,*) icount,numeh(1)*a2d**2,numeh(2)*a2d**2

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
      
      
      
      
      
      
      
      
      
      
      
! ***********************************************************************
      subroutine shiftmuSHF(mu,ieh,dmu0,num,beta,nk,xyk,engk)
      implicit none
      integer(4) :: i,j,sign,nk,ik1,ieh
      real(8),parameter :: pai=dacos(-1.0d0),hdig=1.0d-7
      real(8) :: mu(2),mu2,num,dmu0,beta,mumax,mumin,dense,ffnc
      real(8) :: xyk(0:nk),engk(0:nk,2)
      
! electron ------------------------------------------------------------
! check upper limit
    ! SHF density
      dense=0.0d0
      do ik1=0,nk ! k1
        dense=dense+xyk(ik1)*ffnc(beta,engk(ik1,ieh),mu(ieh))
      end do
      dense=2.0d0*dense/(2.0d0*pai)

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
        
      ! SHF density
        dense=0.0d0
        do ik1=0,nk ! k1
          dense=dense+xyk(ik1)*ffnc(beta,engk(ik1,ieh),mu2)
        end do
        dense=2.0d0*dense/(2.0d0*pai)
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
        do ik1=0,nk ! k1
          dense=dense+xyk(ik1)*ffnc(beta,engk(ik1,ieh),mu2)
        end do
        dense=2.0d0*dense/(2.0d0*pai)

        if(abs(dense-num)/num.lt.hdig) then
!          write(6,*) (mumax+mumin)/2.0d0
          exit
        else if(dense.lt.num) then
          mumin=mu2
        else
          mumax=mu2
        end if
        if(i.gt.100) then
! chemical potential located in the energy gap of finite size mesh pt.
          write(6,*) 'step size of wavenumber (1/nk) is too large'
          stop
        end if
      end do
      mu(ieh)=mu2
!      write(6,*) 'finish',mu,dense

      return
      end subroutine shiftmuSHF
! *********************************************************************




! **********************************************************************
      function screening(beta,engk,nk,xyk,fermi0)
      implicit none
      integer(4) :: nk,ik,ieh
      real(8) :: beta,screening,engk(0:nk,2),xyk(0:nk),fermi0(0:nk,2)
      
      screening=0.0d0
      do ieh=1,2 
        do ik=0,nk ! k
          screening=screening+&
          xyk(ik)*fermi0(ik,ieh)*(1.0d0-fermi0(ik,ieh))
        end do
      end do
      close(24)
      screening=2.0d0*beta*screening

      return
      end function screening
! **********************************************************************



! **********************************************************************
      subroutine Intexp0(nk,Wscexp,yt0,ksc,cffInt,nt0,&
      n0ptk,x0ptk,y0ptk,divk,xk,xt0)
      implicit none
      integer(4) :: nk,ik1,ik2,it,mlim,ieh1,ieh2,nt0
      integer(4) :: n0ptk
      real(8),parameter :: pai=dacos(-1.0d0)
      real(8) :: &
        k1,k2,t0,q,const,Opl2,ksc,Wscexp(0:nk,0:nk),&
        yt0(nt0),S1,Stot,cffInt,x0ptk(n0ptk),y0ptk(n0ptk),&
        St(500,0:500),xint,divk,xk(0:nk),xt0(nt0)
    
    ! interaction (no 0poit correction)
      Wscexp=0.0d0
        do ik2=0,nk ! k2
          k2=xk(ik2)
          do ik1=0,nk ! k1
            k1=xk(ik1) 
            S1=0.0d0
            do it=1,nt0 ! theta_k1k2
              t0=xt0(it) 
              q=dsqrt(k1**2+k2**2-2.0d0*k1*k2*dcos(t0)) ! |q|
              S1=S1+yt0(it)/(q+ksc/(1.0d0+q**2*cffInt))
            end do
            Wscexp(ik1,ik2)=S1/pai
          end do
        end do
 
    ! interaction with 0point correction on 0<ik1=ik2<nk
        do ik1=1,nk-1 ! k1
          k1=xk(ik1)
          Stot=0.0d0
          do ik2=1,n0ptk ! dk
            S1=0.0d0
            do it=1,nt0 ! theta_k1k2
              t0=xt0(it)
              q=dsqrt(k1**2+(k1-x0ptk(ik2))**2-&
              2.0d0*k1*(k1-x0ptk(ik2))*dcos(t0)) ! |q|
              S1=S1+yt0(it)/(q+ksc/(1.0d0+q**2*cffInt))
            end do
            Stot=Stot+y0ptk(ik2)*S1/pai
          end do
          do ik2=1,n0ptk ! dk
            S1=0.0d0
            do it=1,nt0 ! theta_k1k2
              t0=xt0(it)
              q=dsqrt(k1**2+(k1+x0ptk(ik2))**2-&
              2.0d0*k1*(k1+x0ptk(ik2))*dcos(t0)) ! |q|
              S1=S1+yt0(it)/(q+ksc/(1.0d0+q**2*cffInt))
            end do
            Stot=Stot+y0ptk(ik2)*S1/pai
          end do
          
        ! remove adjacent contribution
          Wscexp(ik1,ik1)=&
          (Stot-divk*(Wscexp(ik1,ik1-1)+Wscexp(ik1,ik1+1))/2.0d0)/divk
        end do

    ! unused interaction at ik1=ik2=0
      Wscexp(0,0)=0.0d0
      Wscexp=2.0d0*pai*Wscexp

      return
      end subroutine Intexp0
! **********************************************************************






! ***********************************************************************
      function bfnc(beta,ek)
      implicit none
      real(8) :: bfnc,beta,ek
         bfnc=1.0d0/(dexp(beta*ek)-1.0d0)
      end function bfnc
! ***********************************************************************
! ***********************************************************************
      function ffnc(beta,ek,mu)
      implicit none
      real(8) :: ffnc,beta,ek,mu
         ffnc=1.0d0/(dexp(beta*(ek-mu))+1.0d0)
      end function ffnc
! ***********************************************************************





! ***********************************************************************
      function CH(ksc,nCH,xCH,yCH,cffInt)
      implicit none
      integer(4) :: nCH,iq
      real(8) :: cffInt,ksc,CH,xCH(0:nCH),yCH(0:nCH),q
      real(8),parameter :: pai=dacos(-1.0d0)
      
      CH=0.0d0
      do iq=0,nCH
        q=xCH(iq)
        
        CH=CH+yCH(iq)*ksc/(ksc+q+cffInt*q**3)

        if(ksc/(ksc+q+cffInt*q**3)/CH.le.1.0d-10) exit
      end do
      CH=-pai*CH/(2.0d0*pai)
      
      return
      end function CH
! ***********************************************************************



! ***********************************************************************
      SUBROUTINE gauleg(x1,x2,x,w,n)
      real(8) :: pai=dacos(-1.0d0)
      INTEGER(4) :: n
      REAL(8) :: x1,x2,x(n),w(n)
      REAL(8),PARAMETER :: EPS=3.d-14
      INTEGER(4) :: i,j,m
      REAL(8) :: p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
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
! ***********************************************************************




