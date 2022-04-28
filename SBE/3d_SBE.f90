      Program Tmatrix2d
      implicit none
    ! me and mh="electron and hole mass", mr="reduced mass"        
      real(8),parameter :: me=0.0665d0,mh=0.11d0,mr=me*mh/(me+mh)
    ! mass
      real(8),parameter :: mass(2)=(/me/mr,mh/mr/)
    ! sum of mass
      real(8),parameter :: M(3)=(/2.0d0*me/mr,(me+mh)/mr,2.0d0*mh/mr/)

    ! me0,mh0,mr="3d GaAs e,h,reduced mass"
      real(8),parameter :: me0=0.0665d0,mh0=0.457d0,mr0=me0*mh0/(me0+mh0)

    ! E3d0="3d GaAs binding energy (eV)"
      real(8),parameter :: E3d0=4.2d-3

    ! Esc="energy unit",Eg="band-gap(1.5 eV)"
      real(8),parameter :: Esc=2.0d0,Eg=1.5d0/(E3d0*mr/mr0)/Esc

    ! initial condition read=1: read previdous result
      integer(4),parameter :: iread=0

    ! temp="temperature", num="carrier density"
      real(8),parameter :: temp=0.5d0/Esc,beta=1.0d0/temp,num=1.0d-4
    ! const="SPP parameter"
      real(8),parameter :: const=4.0d0

    ! imaginary number and circular constant
      complex(8),parameter :: im=(0.0d0,1.0d0)
      real(8),parameter :: pai=dacos(-1.0d0)

! parameter --------------------------------------------------------
    ! gamma="\chi broadening", cnv="order of convergence"
      real(8),parameter :: gamma=0.05d0/Esc,cnv=1.0d-4
    ! iend="iteration limit"
      integer(4),parameter :: iend=500
    ! initial shift of chemical potential
      real(8),parameter :: dmu0=0.1d0/Esc
! ------------------------------------------------------------------


! set momentum ------------------------------------------------------
  ! integral parameter in main loop
    ! kmax="k cutoff", divk="k step size"
      real(8),parameter :: kmax=30.0d0,divk=0.05d0
    ! nk="number of k point"
      integer(4),parameter :: nk=int(kmax/divk)
    ! xk and xyk2="point and weight"
      real(8) :: xk(0:nk),xyk2(0:nk),xyk2sqrt(0:nk)

  ! integral parameter for 0 point correction
    ! n0ptk="number of fine k points"
      integer(4),parameter :: n0ptk=500
    ! x0ptk and y0ptk="point and weight"
      real(8) :: x0ptk(n0ptk),y0ptk(n0ptk)
! ----------------------------------------------------------------------


! set omega -------------------------------------------------------
    ! divw="step size of omega", dw="extarnal part"
      real(8),parameter :: divw=1.0d0/40.0d0/Esc,dw=37.5d0/Esc
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
      real(8) :: k1,k2,selfCH,ksc2,ffnc,bfnc,w,ww,p,q,t0,&
           abc,def,ghi,jkl,aion,screening,cffInt
      complex(8) :: zabc,zdef



      write(6,*) 'T (E3d)',Esc/beta
      write(6,*) 'carrier density (1/a3d^3)',num
      write(6,*)
      write(6,*) 'setup'
      write(6,*) 'kmax=',kmax
      write(6,*) 'nk=',nk
      write(6,*) 'dw=',dw*Esc
      write(6,*) 'divw=',divw*Esc
      write(6,*) 'iread=',iread
      write(6,*)








! set momentum ----------------------------------------
  ! integral parameter in main loop
    ! xk="k points"
      do ik=0,nk ! k
         xk(ik)=divk*dble(ik)
      end do
    ! xyk2="weight of radial integral of k"
      xyk2(0)=0.0d0
      do ik=1,nk
         xyk2(ik)=divk**3*dble(ik)**2
      end do
    ! xyk2sqrt=dsqrt(xyk2)
      xyk2sqrt(0)=0.0d0
      do ik=1,nk
         xyk2sqrt(ik)=dsqrt(xyk2(ik))
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
      call gauleg(-1.0d0,1.0d0,xt0,yt0,nt0) ! numerical recipes
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
    ! set initial mu (classical)
      do ieh=1,2
        mu(ieh)=&
        dlog(dsqrt(2.0d0*pai*beta/mass(ieh))**3*num/2.0d0)/beta
      end do
    ! SHF
      call SHF(mu,beta,num,nk,engkSHF,engkf,fermi0,dmu0,&
         xyk2,nt0,yt0,Wscexp,n0ptk,x0ptk,y0ptk,divk,xk,xt0,&
         self0,cnv,const,iend,ksc2,selfCH)
! --------------------------------------------------------------------------


! output SHF results -----------------------------------------------------
      write(6,*) 'SHF results'
      write(6,*) 'total mu (E3d)'
      write(6,*) Esc*(mu(1)+mu(2))
      write(6,*) 'BGR (E3d)'
      write(6,*) (self0(0,1)+self0(0,2))*Esc
      write(6,*) 'Coulomb Hole (E3d)'
      write(6,*) selfCH*Esc

      open(47,file='T0.5n1d-4_SBE_thermodynamic_results.txt',status='unknown')
      write(47,*) 'E chemical potential (E3d)'
      write(47,*) num,mu(1)*Esc
      write(47,*) 'H chemical potential (E3d)'
      write(47,*) num,mu(2)*Esc
      write(47,*) 'total chemical potential (E3d)'
      write(47,*) num,(mu(1)+mu(2))*Esc
      write(47,*)
      
      write(47,*) 'screening wave number (1/a3d)'
      write(47,*) num,ksc2
      write(47,*) 'BGR (E3d)'
      write(47,*) num,(self0(0,1)+self0(0,2))*Esc
      
      write(47,*) 'Coulomb hole self energy'
      write(47,*) num,selfCH*Esc
      write(47,*) 'E screened exchange self energy'
      write(47,*) num,(self0(0,1)-selfCH)*Esc
      write(47,*) 'H screened exchange self energy'
      write(47,*) num,(self0(0,2)-selfCH)*Esc
      write(47,*) 'total screened exchange self energy'
      write(47,*) num,(self0(0,1)+self0(0,2)-2.0d0*selfCH)*Esc
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

    ! no dipole matrix elemetn
      dk=1.0d0
 
   ! kinetic term
      H0=0.0d0
      do ieigen=1,nk+1
        H0(ieigen,ieigen)=engkSHFeh(ieigen)
      end do

    ! Coulomb term
      do ik2=0,nk
        do ik1=0,nk
          H0(ik1+1,ik2+1)=H0(ik1+1,ik2+1)&
          -dsign(1.0d0,pblock0(ik1+1))*dsqrt(abs(pblock0(ik1+1)))*&
          xyk2sqrt(ik1)*&
          Wscexp(ik1,ik2)*dsqrt(abs(pblock0(ik2+1)))*&
          xyk2sqrt(ik2)/(2.0d0*pai**2)
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
          write(6,*) 'complex eigen value',WI(ieigen)*Esc
        end if
      end do
      
    ! weight of pole
      weight=0.0d0
      weight0=0.0d0
      do ieigen=1,nk+1
        do ik2=0,nk
          do ik1=0,nk
            weight(ieigen)=weight(ieigen)+&
            dk(ik2)*xyk2sqrt(ik2)*VR(ik2+1,ieigen)*dsqrt(abs(pblock0(ik2+1)))*&
            dk(ik1)*xyk2sqrt(ik1)*VL(ik1+1,ieigen)*dsqrt(abs(pblock0(ik1+1)))*&
            dsign(1.0d0,pblock0(ik1+1))
          end do
        end do
        weight0(ieigen)=weight0(ieigen)+&
        dk(ieigen-1)**2*xyk2(ieigen-1)*pblock0(ieigen)
      end do
      weight=weight/(2.0d0*pai**2)
      weight0=weight0/(2.0d0*pai**2)

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
      open(89,file='T0.5n1d-4_SBE_results.txt',status='unknown')
      write(89,*) mu(1),mu(2)
      do ieigen=1,nk+1
        write(89,*) WR(ieigen),weight(ieigen),weight0(ieigen),engkSHFeh(ieigen)
      end do
      close(89)
      write(6,*) 'SBE lowest energy level'
      write(6,*) WR(1)*Esc
      write(47,*) 'SBE lowest energy level'
      write(47,*) num,WR(1)*Esc
! -----------------------------------------------------------------      
      

! optical response -------------------------------------------------
460   open(222,file='T0.5n1d-4_SBE_Imx.txt',status='unknown')
      open(125,file='T0.5n1d-4_SBE_Rex.txt',status='unknown')
      open(123,file='T0.5n1d-4_SBE_PL.txt',status='unknown')
      open(32,file='T0.5n1d-4_SBE_Imx0.txt',status='unknown')
      open(132,file='T0.5n1d-4_SBE_Rex0.txt',status='unknown')
      open(35,file='T0.5n1d-4_SBE_PL0.txt',status='unknown')
      open(42,file='T0.5n1d-4_SBE_Imx_per_Imx0.txt',status='unknown')
      open(142,file='T0.5n1d-4_SBE_Rex_per_Rex0.txt',status='unknown')

      do iw=nwnn1,nweh2
        zabc=0.0d0
        zdef=0.0d0
        abc=0.0d0
        def=0.0d0
        w=xw(iw)
        if(abs(w*Esc).lt.20.0d0) then
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
        write(222,*) w*Esc,dimag(zabc)/Esc
        write(125,*) w*Esc,dreal(zabc)/Esc
        write(123,*) w*Esc,abc/Esc
        write(32,*) w*Esc,dimag(zdef)/Esc
        write(132,*) w*Esc,dreal(zdef)/Esc
        write(35,*) w*Esc,def/Esc
        write(42,*) w*Esc,dimag(zabc)/dimag(zdef)
        write(142,*) w*Esc,dreal(zabc)/dreal(zdef)
        end if
      end do
! -----------------------------------------------------------------


! pair susceptibility (gamma=0) -----------------------------------
      open(56,file='T0.5n1d-4_SBE_pair_susceptibility.txt',status='unknown')
      abc=0.0d0
      def=0.0d0
      do ieigen=1,nk+1
        abc=abc+weight(ieigen)/((mu(1)+mu(2)-WR(ieigen)))
        def=def+weight0(ieigen)/((mu(1)+mu(2)-engkSHFeh(ieigen)))
      end do
      abc=-2.0d0*abc
      def=-2.0d0*def
      write(56,*) "Rex at w=mu (1/(a3d**3*E3d))" 
      write(56,*) num,abc/Esc
      write(56,*) "Rex0 at w=mu (1/(a3d**3*E3d))"
      write(56,*) num,def/Esc
      write(56,*) "Rex/Rex0 at w=mu" 
      write(56,*) num,abc/def
      write(56,*) "Rex-Rex0 at w=mu (1/(a3d**3*E3d))"
      write(56,*) num,(abc-def)/Esc
      close(56)
! -------------------------------------------------------------

      stop

      
      end program Tmatrix2d
! ************************************************************************













      
      
! ***********************************************************************
      subroutine SHF(mu,beta,num,nk,engkSHF,engkf,fermi0,dmu0,&
           xyk2,nt0,yt0,Wscexp,n0ptk,x0ptk,y0ptk,divk,xk,xt0,&
           self0,cnv,const,iend,ksc2,selfCH)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0,E2d=4.0d0,a2d=0.5d0
      integer(4) :: i,ieh,ik,ik1,ik2,nk,icount,nt0,n0ptk,iend
      real(8) :: mu(2),engkSHF(0:nk,2),engkf(0:nk,2),fermi0(0:nk,2),&
           xyk2(0:nk),Wscexp(0:nk,0:nk),yt0(nt0),&
           x0ptk(n0ptk),y0ptk(n0ptk),xk(0:nk),xt0(nt0),&
           self0(0:nk,2),numeh(2),numehold(2)
      real(8) :: beta,num,const,dmu0,f,ksc2,cffInt,divk,selfCH,&
           abc,k1,cnv,screening,ffnc,CH

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
          call shiftmuSHF(mu,ieh,dmu0,num,beta,nk,xyk2,engkSHF)
          do ik=0,nk ! k
            fermi0(ik,ieh)=ffnc(beta,engkSHF(ik,ieh),mu(ieh))
          end do
        end do

      ! effective interaction
        ksc2=screening(beta,engkSHF,nk,xyk2,fermi0)
        cffInt=const*ksc2/64.0d0/num/pai
        
        call Intexp0(nk,Wscexp,yt0,ksc2,cffInt,nt0,&
        n0ptk,x0ptk,y0ptk,divk,xk,xt0)        

      ! Coulomb hole enery
        selfCH=CH(num,ksc2,const)

      ! exchange energy
        self0=0.0d0
        do ieh=1,2 
           do ik2=0,nk ! k2
              abc=xyk2(ik2)*fermi0(ik2,ieh)
              do ik1=0,nk ! k1
                 self0(ik1,ieh)=self0(ik1,ieh)+Wscexp(ik1,ik2)*abc
              end do
           end do
        end do
        self0=-self0/(2.0d0*pai**2)+selfCH

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
              abc=abc+xyk2(ik1)*ffnc(beta,engkSHF(ik1,ieh),mu(ieh))
           end do
           numeh(ieh)=abc
        end do
        numeh=2.0d0*numeh/(2.0d0*pai**2)

        write(6,*) icount,numeh(1),numeh(2)

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
      subroutine shiftmuSHF(mu,ieh,dmu0,num,beta,nk,xyk2,engk)
      implicit none
      integer(4) :: i,j,sign,nk,ik1,ieh
      real(8),parameter :: pai=dacos(-1.0d0),hdig=1.0d-7
      real(8) :: mu(2),mu2,num,dmu0,beta,mumax,mumin,dense,ffnc
      real(8) :: xyk2(0:nk),engk(0:nk,2)
      
! electron ------------------------------------------------------------
! check upper limit
    ! SHF density
      dense=0.0d0
      do ik1=0,nk ! k1
        dense=dense+xyk2(ik1)*ffnc(beta,engk(ik1,ieh),mu(ieh))
      end do
      dense=2.0d0*dense/(2.0d0*pai**2)

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
          dense=dense+xyk2(ik1)*ffnc(beta,engk(ik1,ieh),mu2)
        end do
        dense=2.0d0*dense/(2.0d0*pai**2)
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
          dense=dense+xyk2(ik1)*ffnc(beta,engk(ik1,ieh),mu2)
        end do
        dense=2.0d0*dense/(2.0d0*pai**2)

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
      function screening(beta,engk,nk,xyk2,fermi0)
      implicit none
      integer(4) :: nk,ik,ieh
      real(8) :: beta,screening,engk(0:nk,2),xyk2(0:nk),fermi0(0:nk,2)
      real(8),parameter :: pai=dacos(-1.0d0)
      
      screening=0.0d0
      do ieh=1,2 
        do ik=0,nk ! k
          screening=screening+&
          xyk2(ik)*fermi0(ik,ieh)*(1.0d0-fermi0(ik,ieh))
        end do
      end do
      close(24)
      screening=4.0d0*beta*screening/pai

      return
      end function screening
! **********************************************************************



! **********************************************************************
      subroutine Intexp0(nk,Wscexp,yt0,ksc2,cffInt,nt0,&
      n0ptk,x0ptk,y0ptk,divk,xk,xt0)
      implicit none
      integer(4) :: nk,ik1,ik2,it,mlim,ieh1,ieh2,nt0
      integer(4) :: n0ptk
      real(8),parameter :: pai=dacos(-1.0d0)
      real(8) :: &
        k1,k2,t0,q,const,Opl2,ksc2,Wscexp(0:nk,0:nk),&
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
              q=dsqrt(k1**2+k2**2-2.0d0*k1*k2*t0) ! |q|
              S1=S1+yt0(it)/(q**2+ksc2/(1.0d0+q**2*cffInt))
            end do
            Wscexp(ik1,ik2)=S1/2.0d0
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
              2.0d0*k1*(k1-x0ptk(ik2))*t0) ! |q|
              S1=S1+yt0(it)/(q**2+ksc2/(1.0d0+q**2*cffInt))
            end do
            Stot=Stot+y0ptk(ik2)*S1/2.0d0
          end do
          do ik2=1,n0ptk ! dk
            S1=0.0d0
            do it=1,nt0 ! theta_k1k2
              t0=xt0(it)
              q=dsqrt(k1**2+(k1+x0ptk(ik2))**2-&
              2.0d0*k1*(k1+x0ptk(ik2))*t0) ! |q|
              S1=S1+yt0(it)/(q**2+ksc2/(1.0d0+q**2*cffInt))
            end do
            Stot=Stot+y0ptk(ik2)*S1/2.0d0
          end do
          
        ! remove adjacent contribution
          Wscexp(ik1,ik1)=&
          (Stot-divk*(Wscexp(ik1,ik1-1)+Wscexp(ik1,ik1+1))/2.0d0)/divk
        end do

    ! unused interaction at ik1=ik2=0
      Wscexp(0,0)=0.0d0
      Wscexp=4.0d0*pai*Wscexp

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
      function CH(num,ksc2,const)
      implicit none
      real(8) :: const,rs,ksc2,u,CH,num,bp,bm,cp,cm
      real(8),parameter :: pai=dacos(-1.0d0)
      integer(4) :: i
      
      u=dsqrt(4.0d0*pai*num)*4.0d0/(ksc2*dsqrt(const))
      bp=(1.0d0+1.0d0/dsqrt(1.0d0-4.0d0/u**2))
      bm=(1.0d0-1.0d0/dsqrt(1.0d0-4.0d0/u**2))
      cp=dsqrt(ksc2)*u/dsqrt(2.0d0)*&
        dsqrt(1.0d0+dsqrt(1.0d0-4.0d0/u**2))
      cm=dsqrt(ksc2)*u/dsqrt(2.0d0)*&
        dsqrt(1.0d0-dsqrt(1.0d0-4.0d0/u**2))

      if(u.le.2.0d0) then
        CH=-0.25d0*dsqrt(ksc2)*u*&
        (dsqrt(1.0d0+2.0d0/u)-1.0d0/dsqrt(1.0d0+2.0d0/u))
      else
        CH=-0.25d0*dsqrt(ksc2)*u*(&
        dsqrt(1.0d0+dsqrt(1.0d0-4.0d0/u**2))*&
        (1.0d0-1.0d0/dsqrt(1.0d0-4.0d0/u**2))/dsqrt(2.0d0)+&
        dsqrt(1.0d0-dsqrt(1.0d0-4.0d0/u**2))*&
        (1.0d0+1.0d0/dsqrt(1.0d0-4.0d0/u**2))/dsqrt(2.0d0))
      end if
      
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




