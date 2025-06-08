subroutine freqem1(probin,nind,freqest)
implicit none
integer ::nind
integer, parameter :: dp = selected_real_kind(14)
real(dp) ::f1,freqest,probin(3*nind)

f1=freq0(probin,nind)
freqest=freqem(probin,f1,nind)

contains

!##### initial guess
function freq0(myprobs,n)
implicit none
integer ::i,k
integer ::n
real(dp) ::myprobs(3*n),freq0,gt

k=0;freq0=0.d0
do i=1,n
 gt=myprobs(3*i-2)+myprobs(3*i-1)+myprobs(3*i)
 if(gt < 0.05)cycle !### considered missing
 k=k+1
 freq0=freq0+2.0*myprobs(3*i-2)+myprobs(3*i-1)
enddo
if(k>0)freq0=freq0/(2.d0*k)
return

end function

function freqem(myprobs,f0,n)
implicit none
integer ::i,k,l,n,keep(n)
real(dp) ::myprobs(3*n),eps,diff,f0,freqem,gt
real(dp) ::genolik(n,3),genofreq(3),genoprob(3)

keep=0;eps=1e-4;diff=1.d0;freqem=f0
genolik=0.d0

do i=1,n
 gt=myprobs(3*i-2)+myprobs(3*i-1)+myprobs(3*i)
 if(gt < 0.05)cycle !### considered missing
 keep(i)=1
 genolik(i,1)=myprobs(3*i-2)
 genolik(i,2)=myprobs(3*i-1)
 genolik(i,3)=myprobs(3*i)
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
   genoprob=genolik(i,:)*genofreq
   genoprob=genoprob/sum(genoprob)
   freqem=freqem + 2.d0*genoprob(1) + 1.d0*genoprob(2)
  end do
  freqem=freqem/(2.d0*k)
  diff=abs(diff-freqem)
enddo

return

end function

end subroutine
