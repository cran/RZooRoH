subroutine freqem2(phredin,nind,freqest)
implicit none
integer, parameter :: dp = selected_real_kind(14)
integer ::i,nind
real(dp) ::f1,freqest
real(dp) ::phredin(3*nind)

!### rescale phreds
do i=1,nind
 phredin((3*i-2):(3*i))=phredin((3*i-2):(3*i))-minval(phredin((3*i-2):(3*i)))
enddo

f1=freq0(phredin,nind)
freqest=freqem(phredin,f1,nind)

contains

!##### initial guess
function freq0(myphreds,n)
implicit none
integer ::i,k
integer ::n
real(dp) ::myphreds(3*n),freq0,gt,g11,g12,g22

k=0;freq0=0.d0
do i=1,n
  gt=myphreds(3*i-2)+myphreds(3*i-1)+myphreds(3*i)
  if(gt<0.01)cycle
  k=k+1
  g11=10**(-myphreds(3*i-2)/10)
  g12=10**(-myphreds(3*i-1)/10)
  g22=10**(-myphreds(3*i)/10)
  gt=g11+g12+g22
  g11=g11/gt;g12=g12/gt;g22=g22/gt
  freq0=freq0+2.0*g11+g12
enddo
if(k>0)freq0=freq0/(2.d0*k)
return

end function

function freqem(myphreds,f0,n)
implicit none
integer ::i,k,l,n,keep(n)
real(dp) ::myphreds(3*n),eps,diff,f0,freqem,gt
real(dp) ::genolik(n,3),genofreq(3),genoprob(3)

keep=0;eps=1e-4;diff=1.d0;freqem=f0
genolik=0.d0

do i=1,n
 gt=myphreds(3*i-2)+myphreds(3*i-1)+myphreds(3*i)
 if(gt < 0.01)cycle !### considered missing
 keep(i)=1
 genolik(i,1)=10**(-myphreds(3*i-2)/10)
 genolik(i,2)=10**(-myphreds(3*i-1)/10)
 genolik(i,3)=10**(-myphreds(3*i)/10)
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
