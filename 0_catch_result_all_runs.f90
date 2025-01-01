program main
  implicit none
  integer i,j,ios,ifile0,ifile
  integer flag
  integer iter0,iter,iorig0
  integer numIter0
  integer cnt
  integer,parameter::numFile=200
  integer,parameter::numIter=30000
  integer,parameter::numResoShell=9
  character(len=80) fileName0, fileName
  real meanPhaseErrorFile(numFile,0:numIter-1),meanPhaseError(0:numResoShell-1),meanPhaseError0
  real meanPhaseErrorElite(0:numIter-1)
  real corrCoefFile(numFile,0:numIter-1),corrCoef(0:numResoShell-1)
  real protMaskMatchFile(numFile,0:numIter-1),protMaskMatch
  real RfreeFile(numFile,0:numIter-1),Rfree
  real RworkFile(numFile,0:numIter-1),Rwork
  real F000File(numFile,0:numIter-1),F000
  real scaleFactFile(numFile,0:numIter-1),scaleFact
  real centMassDeviDistFile(numFile,0:numIter-1),centMassDeviDist
  real ncsAxisDeviAnglFile(numFile,0:numIter-1),ncsAxisDeviAngl
  character(len=20) :: myfmt
  character(len=80) :: string
  integer::deltaIter = 0

  ifile0=0
  cnt=0
  do ifile=0,numFile,1

    if(ifile<10)then
      write(fileName0,'(A20,I1)')'directPhasing/run_1/',ifile
    else if(ifile<100)then
      write(fileName0,'(A20,I2)')'directPhasing/run_1/',ifile
    else if(ifile<1000)then
      write(fileName0,'(A20,I3)')'directPhasing/run_1/',ifile
    end if

	!---------------------------------------------------------------------------

	fileName = trim(filename0) // '_mean_phase_error_reso_sphere.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle

	ifile0=ifile0+1

    do 
      read(10,'(I8,9(1X,F8.1))',iostat=ios) iter,(meanPhaseError(i),i=0,8)
      if(ios/=0)exit
      meanPhaseErrorFile(ifile0,iter)=meanPhaseError(0)
    end do
    close(10)
    if(meanPhaseError(0)<72) cnt=cnt+1
    close(10)
	numIter0 = iter

	!---------------------------------------------------------------------------

    fileName = 'directPhasing/run_1/0_mean_phase_error_reso_sphere_of_conv_dens.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle

    do
      read(10,'(I8,9(1X,F8.1))',iostat=ios) iter,(meanPhaseError(i),i=0,8)
      if(ios/=0)exit
      meanPhaseErrorElite(iter)=meanPhaseError(0)
    end do
    close(10)

    !---------------------------------------------------------------------------

	fileName = trim(filename0) // '_corr_coef.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle

    do 
      read(10,'(I8,9(1X,F8.3))',iostat=ios) iter,(corrCoef(i),i=0,8)
      if(ios/=0)exit
      corrCoefFile(ifile0,iter)=corrCoef(0)
    end do
    close(10)

	!---------------------------------------------------------------------------

	fileName = trim(filename0) // '_prot_mask_match.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle

    do 
      read(10,'(I8,1X,F8.3)',iostat=ios) iter,protMaskMatch
      if(ios/=0)exit
      protMaskMatchFile(ifile0,iter)=protMaskMatch
    end do
    close(10)

	!---------------------------------------------------------------------------

	fileName = trim(filename0) // '_R_free.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle
    do 
      read(10,'(I8,1X,F8.3)',iostat=ios) iter,Rfree
      if(ios/=0)exit
      RfreeFile(ifile0,iter)=Rfree
    end do
    close(10)

	!---------------------------------------------------------------------------

	fileName = trim(filename0) // '_R_work.txt'
    open(unit=10,file=fileName,status='old',iostat=ios)
    if(ios/=0)cycle
    do 
      read(10,'(I8,1X,F8.3)',iostat=ios) iter,Rwork
      if(ios/=0)exit
      RworkFile(ifile0,iter)=Rwork
    end do
    close(10)

  end do
  write(*,*)cnt

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  open(unit=11,file='0_log_of_deltaPhi_all_runs.txt',status='replace')
  open(unit=12,file='0_log_of_corrCoef_all_runs.txt',status='replace')
  open(unit=13,file='0_log_of_protMaskMatch_all_runs.txt',status='replace')
  open(unit=14,file='0_log_of_Rfree_all_runs.txt',status='replace')
  open(unit=15,file='0_log_of_Rwork_all_runs.txt',status='replace')
  open(unit=16,file='0_log_of_deltaPhi_Elite.txt',status='replace')

  do iter=0,numIter0-1,1

    if(meanPhaseErrorFile(1,iter)==0 .or. RworkFile(1,iter)==0) then
        deltaIter = deltaIter + 1
        cycle
    end if
    if(mod(iter,10)/=0)cycle

    write(myfmt, '("(I8,",I0,"(A1,F6.1))")') ifile0
    write(11, fmt=myfmt ) iter-deltaIter, (',',meanPhaseErrorFile(i,iter), i=1,ifile0)

    write(myfmt, '("(I8,",I0,"(A1,F6.3))")') ifile0
    write(12, fmt=myfmt ) iter-deltaIter, (',',corrCoefFile(i,iter), i=1,ifile0)

    write(13, fmt=myfmt ) iter-deltaIter, (',',protMaskMatchFile(i,iter), i=1,ifile0)

    if(iter/=0) write(14, fmt=myfmt ) iter-deltaIter, (',',RfreeFile(i,iter), i=1,ifile0)

    if(iter/=0) write(15, fmt=myfmt ) iter-deltaIter, (',',RworkFile(i,iter), i=1,ifile0)

    write(16, '(I8,A1,F6.1)') iter-deltaIter,',',meanPhaseErrorElite(iter)

  end do
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)
  stop
end 
