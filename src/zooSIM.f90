subroutine zoosim(nlayer,nchr,npos,freq,posi,chr_bounds,as,Fs,gerr,nsim,genosim,Frealized,classes)
implicit none
INTEGER, PARAMETER ::DP = SELECTED_REAL_KIND(14)
integer,parameter :: k1 = selected_int_kind(1)
integer ::i,j,k,l,io,bp,nchr,nlayer,nind,npos,n,chr,state,id,nrec,chrom,chromold
integer ::nsim,nsnps,chr_bounds(nchr,2)
integer ::chr_limits(nchr,4),posi(npos),genosim(nsim,npos),classes(nsim,npos)
integer ::nF(nlayer+1),nFO(nlayer+1),Ncenso(nlayer+1)
integer(k1) ::genome(300000000),level
integer ::next_pos,current_pos
real(dp) ::a,F,xran,Ft,Glen,startp,Lsim,val,gerr,expornd
real(dp) ::freq(npos),as(nlayer),Fs(nlayer),Fsnp(nlayer+1),Fsim(nlayer+1),age(nlayer+1)
real(dp) ::Frealized(nsim,nlayer)

call rndstart()

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
    xran=expornd()
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
   if(k==chr_limits(chr,3))nFO(state)=nFO(state)+1
   if(k>chr_limits(chr,3))then
    if(state/=abs(genome(posi(k-1))))nFO(state)=nFO(state)+1
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

 do k=2,(nlayer+1)
  Frealized(id,k-1)=Fsnp(k)
 enddo

enddo ! id

call rndend()

contains

subroutine sample_state
implicit none
real(dp) ::ran1,unifrnd

ran1 = unifrnd()
if(ran1 < F)then
  state=1
else
  state=2
endif

end subroutine

subroutine sample_geno(hs,marker,indiv)
implicit none
integer ::hs,marker,indiv
real(dp) ::ran1,unifrnd

ran1 = unifrnd()
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

ran1 = unifrnd()
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

end subroutine
