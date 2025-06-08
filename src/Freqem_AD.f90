subroutine freqem3(ads,nind,freqest)
implicit none
integer, parameter :: dp = selected_real_kind(14)
integer ::nind
integer ::ads(2*nind)
real(dp) ::f1,freqest

f1=freq0(ads,nind)
freqest=freqem(ads,f1,nind)

contains

!##### initial guess
function freq0(myads,n)
implicit none
integer ::i,k
integer ::n,adt,myads(2*n)
real(dp) ::freq0

k=0;freq0=0.d0
do i=1,n
 adt=myads(2*i-1)+myads(2*i)
 if(adt == 0)cycle !### missing
 k=k+1
 freq0=freq0+(1.d0*myads(2*i-1))/(1.d0*adt)
enddo
if(k>0)freq0=freq0/(1.d0*k)
return

end function

function freqem(myads,f0,n)
implicit none
integer ::i,k,l,n,keep(n),ad1,ad2,adt
integer ::myads(2*n)
real(dp) ::eps,diff,f0,freqem,seqerr
real(dp) ::genolik(n,3),genofreq(3),genoprob(3)

keep=0;eps=1e-4;diff=1.d0;freqem=f0
genolik=0.d0;seqerr=1e-3

do i=1,n
 ad1=myads(2*i-1);ad2=myads(2*i)
 adt=ad1+ad2
 if(adt == 0)cycle !### missing
 keep(i)=1
 genolik(i,1)=log(1.d0-seqerr)*ad1+log(seqerr)*ad2
 genolik(i,2)=log(0.5)*adt   !#### (0.5**ad1)*(0.5**ad2)
 genolik(i,3)=log(1.d0-seqerr)*ad2+log(seqerr)*ad1
 genolik(i,:)=genolik(i,:)-maxval(genolik(i,:))
enddo
k=sum(keep)
if(k==0)then
 freqem=0.d0
 return
endif

l=0
do while (diff>eps)
  genofreq(1)=freqem**2; genofreq(2)=2.d0*freqem*(1.d0-freqem); genofreq(3)=(1.d0-freqem)**2
  diff=freqem; freqem=0.d0
  l=l+1
  do i=1,n
   if(keep(i)==0)cycle
   genoprob=genolik(i,:)+log(genofreq)
   genoprob=exp(genoprob)/sum(exp(genoprob))
   freqem=freqem + 2.d0*genoprob(1) + 1.d0*genoprob(2)
  end do
  freqem=freqem/(2.d0*k)
  diff=abs(diff-freqem)
enddo

return

end function

end subroutine
