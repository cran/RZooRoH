subroutine zooViterbi(nclust,nchr,npos,pemission,chr_limits,as,Fs,posi,states)
implicit none
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
integer ::i,j,l,k,nclust,nchr,npos,chr,fpos,lpos,hs,hs2,pos(1),phi(nclust,npos)
integer ::isF(nclust),chr_limits(nchr,2),IT(nclust,nclust),posi(npos),states(npos)
real(dp) ::pemission(npos,2),as(nclust),Fs(nclust),delta(nclust,npos)
real(dp) ::trans(nclust,nclust),pinit(nclust)
real(dp) ::sumF(2),F,a,r,gr,loglik,pmax,val(nclust),mind
real(dp), parameter ::Morgan=100000000.d0

pinit=Fs
!delta=-1000000000.0;phi=0;states=0
phi=0;states=0
IT=1;isF=1;isF(nclust)=0

do chr=1,nchr

fpos=chr_limits(chr,1);lpos=chr_limits(chr,2)

! initialisation

 do i=1,nclust
  delta(i,fpos)=log(pinit(i))+log(pemission(fpos,isF(i)+1))
  phi(i,fpos)=0
 enddo

! recursion

 do k=fpos+1,lpos
 trans=TM(nclust,IT,as,Fs,k,npos,posi)
 do i=1,nclust
  do j=1,nclust
    val(j)=delta(j,k-1)+log(trans(j,i))
  enddo
  pos=maxloc(val)
  phi(i,k)=pos(1)
  delta(i,k)=val(pos(1))+log(pemission(k,isF(i)+1))
 enddo
 enddo

! termination

 pmax=maxval(delta(:,lpos))
! print*,"Max prob. ::",pmax ! check that probabilities don't underflow
 pos=maxloc(delta(:,lpos))
 states(lpos)=pos(1)

! path

 do k=lpos-1,fpos,-1
  states(k)=phi(states(k+1),k+1)
 enddo

enddo ! end chr

contains

function TM(K,IT,as,Fs,snp,nsnps,mpos)
implicit none
integer ::K,hs,hs2,snp,nsnps,IT(K,K),mpos(nsnps)
real(dp) ::a,r,as(K),Fs(K),TM(K,K),sumF(2)
real(dp), parameter ::Morgan=100000000.d0


TM=0.d0
do hs=1,K
  sumF(1)=dot_product(Fs,IT(hs,:))
  a=as(hs)
  r=dexp(-a*(mpos(snp)-mpos(snp-1))/Morgan)
  do hs2=1,K
   if(IT(hs,hs2)==0)cycle
   TM(hs,hs2)=(1-r)*Fs(Hs2)/sumF(1)
  enddo
  TM(hs,hs)=TM(hs,hs)+r
enddo

end function

function pexit(K,IT,as,Fs,snp,nsnps,mpos)
implicit none
integer ::K,hs,hs2,snp,nsnps,IT(K,K),mpos(nsnps)
real(dp) ::a,r,as(K),Fs(K),pexit(K,K),sumF(2)
real(dp), parameter ::Morgan=100000000.d0

pexit=0.d0
do hs=1,K
  sumF(1)=dot_product(Fs,IT(hs,:))
  a=as(hs)
  r=dexp(-a*(mpos(snp)-mpos(snp-1))/Morgan)
  do hs2=1,K
   if(IT(hs,hs2)==0)cycle
   pexit(hs,hs2)=(1-r)*Fs(Hs2)/sumF(1)
  enddo
enddo

end function


end subroutine


