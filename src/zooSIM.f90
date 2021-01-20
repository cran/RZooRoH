subroutine zoosim(nlayer,nchr,npos,freq,posi,chr_bounds,as,Fs,gerr,nsim,genosim,Frealized,classes)
implicit none
INTEGER, PARAMETER ::DP = SELECTED_REAL_KIND(14)
integer,parameter :: k1 = selected_int_kind(1)
integer ::i,j,k,l,io,bp,nchr,nlayer,nind,npos,n,chr,state,id,nrec,chrom,chromold
integer ::nsim,nsnps,chr_bounds(nchr,2)
integer ::chr_limits(nchr,4),posi(npos),genosim(nsim,npos),classes(nsim,npos)
integer ::nF(nlayer+1),nFO(nlayer+1),Ncenso(nlayer+1)
!integer(kind=1) ::genome(300000000)
integer ::genome(300000000)
integer ::next_pos,current_pos,level
real(dp) ::a,F,xran,Ft,Glen,startp,Lsim,val,gerr
real(dp) ::freq(npos),as(nlayer),Fs(nlayer),Fsnp(nlayer+1),Fsim(nlayer+1),age(nlayer+1)
real(dp) ::Frealized(nsim,nlayer)
character*10 ::all1,all2
character*20 ::mname

call init_random_seed()
!call random_seed()

!open(12,file='geno_sim.txt',status='replace')
!open(13,file='siminfo.txt',status='replace')
!open(14,file='states_sim.txt',status='replace')

do chrom=1,nchr
 chr_limits(chrom,3)=chr_bounds(chrom,1)
 chr_limits(chrom,4)=chr_bounds(chrom,2)
 chr_limits(chrom,1)=posi(chr_limits(chrom,3))
 chr_limits(chrom,2)=posi(chr_limits(chrom,4))
enddo

Glen=0.d0
do i=1,nchr
 Glen=Glen+1.d0*chr_limits(i,2)
enddo

genosim=0
do id=1,nsim
 nrec=0;nF=0;Ft=0.000;Fsnp=0.00;Ncenso=0;Fsim=0.000;Lsim=0.00
 nFO=0
 do chr=1,nchr
  genome=1 ! all background

 do level=1,nlayer
  a=as(level);F=Fs(level)

   next_pos=0;current_pos=0
   do while (next_pos < chr_limits(chr,2))
    call sample_state
    xran=random_exponential()
    next_pos=current_pos+floor(100000000.d0*xran/a)
    if(next_pos > chr_limits(chr,2))next_pos=chr_limits(chr,2)
    if(state==1)then
      genome((current_pos+1):next_pos)=(level+1)
      genome((current_pos+1))=-(level+1)
    else
      if(level==1)genome((current_pos+1))=-1
    endif
    current_pos=next_pos
   enddo
 enddo ! layers

  do k=1,chr_limits(chr,2)
    state=abs(genome(k))
    Fsim(state)=Fsim(state)+1.d0
    if(state==1)Lsim=Lsim+1.d0
    if(k==1)nF(state)=nF(state)+1
    if(k>1)then
      if(state/=abs(genome(k-1)))nF(state)=nF(state)+1
      if(state==abs(genome(k-1)) .and. genome(k)<0)nF(state)=nF(state)+1
    endif
    if(k==chr_limits(chr,2))Ncenso(state)=Ncenso(state)+1
  enddo

  do k=chr_limits(chr,3),chr_limits(chr,4)
   state=abs(genome(posi(k)))
   Fsnp(state)=Fsnp(state)+1.d0
   if(k==chr_limits(chr,3))nFO(state)=nFO(state)+1.d0
   if(k>chr_limits(chr,3))then
    if(state/=abs(genome(posi(k-1))))nFO(state)=nFO(state)+1.d0
   endif
   classes(id,k)=state
   if(state==1)then
     state=2
   else
     state=1
   endif
   call sample_geno(state,k,id)
  enddo

  enddo ! chr
 do k=2,(nlayer+1)
  age(k)=0.d0
  if(Fsim(k)/=0.00)age(k)=1.d0*(nF(k)-Ncenso(k))/(Fsim(k)/100000000.d0)
  Fsnp(k)=Fsnp(k)/(1.d0*npos)
  Fsim(k)=Fsim(k)/Glen
 enddo
 nrec=sum(nF)

!! write(13,'(i3,<2*nlayer+1>(1x,i6),<2*nlayer>(1x,f9.6),<nlayer>(1x,i6),<nlayer>(1x,f9.2))')id,&
! write(13,*)id,&
! (nF(k),k=2,nlayer+1),nrec,(nFO(k),k=2,nlayer+1),&
! (Fsim(k),Fsnp(k),k=2,nlayer+1),(Ncenso(k),k=2,nlayer+1),(age(k),k=2,nlayer+1)
 do k=2,(nlayer+1)
  Frealized(id,k)=Fsnp(k)
 enddo

enddo ! id

all1='all1';all2='all2'

!k=0
!do chr=1,nchr
! do k=chr_limits(chr,3),chr_limits(chr,4)
!  bp=posi(k)
!    write(12,'(i2,1x,i10,2(1x,a4))',advance='no')chr,bp,all1,all2
!    write(14,'(i2,1x,i10,2(1x,a4))',advance='no')chr,bp,all1,all2
!    do id=1,nsim-1
!      write(12,'(1x,i1)',advance='no')genosim(id,k)
!      write(14,'(1x,i2)',advance='no')classes(id,k)
!    enddo
!    write(12,'(1x,i1)')genosim(nsim,k)
!    write(14,'(1x,i2)')classes(nsim,k)
! enddo
!enddo

! close(12)
! close(13)
! close(14)

contains

subroutine sample_state
implicit none
real(dp) ::ran1

call random_number(ran1)
if(ran1 < F)then
  state=1
else
  state=2
endif

end subroutine

subroutine sample_geno(hs,marker,indiv)
implicit none
integer ::hs,marker,indiv
real(dp) ::ran1

call random_number(ran1)
if(hs==1)then
 if(ran1 < freq(marker))then
  genosim(indiv,marker)=2
 else
  genosim(indiv,marker)=0
 endif
else if(hs==2)then
 if(ran1 < freq(marker)**2)then
  genosim(indiv,marker)=2
 else if(ran1 > (1.00-(1.00-freq(marker))**2))then
  genosim(indiv,marker)=0
 else
  genosim(indiv,marker)=1
 endif
endif

call random_number(ran1)
if(ran1 <= 0.00000)then
 if(ran1 <= 0.00000)then
  if(genosim(indiv,marker)==0)then
    genosim(indiv,marker)=1
  else if(genosim(indiv,marker)==1)then
    genosim(indiv,marker)=0
  else if(genosim(indiv,marker)==2)then
    genosim(indiv,marker)=0
  endif
 else
  if(genosim(indiv,marker)==0)then
    genosim(indiv,marker)=2
  else if(genosim(indiv,marker)==1)then
    genosim(indiv,marker)=2
  else if(genosim(indiv,marker)==2)then
    genosim(indiv,marker)=1
  endif
 endif
endif

end subroutine

   SUBROUTINE init_random_seed()
     INTEGER :: i, n, clock
     INTEGER, DIMENSION(:), ALLOCATABLE :: seed

     CALL RANDOM_SEED(size = n)
     ALLOCATE(seed(n))

     CALL SYSTEM_CLOCK(COUNT=clock)

     seed = clock + 37 * (/ (i - 1, i = 1, n) /)
     CALL RANDOM_SEED(PUT = seed)

     DEALLOCATE(seed)
   END SUBROUTINE init_random_seed


FUNCTION random_exponential() RESULT(fn_val)

! copied from the random.f90 module
! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
! TO EXP(-random_exponential), USING INVERSION.

REAL  :: fn_val

!     Local variable
REAL  :: r

DO
  CALL RANDOM_NUMBER(r)
!  IF (r > zero) EXIT
   IF (r > 0.d0) EXIT
END DO

fn_val = -LOG(r)
RETURN

END FUNCTION random_exponential


end subroutine
