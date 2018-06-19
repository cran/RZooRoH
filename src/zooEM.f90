subroutine zooEM(nclust,nchr,npos,pemission,chr_limits,as,Fs,posi,loglik,gamma,estimateR,&
           onerate,niter,minrate,maxrate,convthr,iter)
implicit none
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
integer ::i,j,l,k,nclust,nchr,npos,chr,fpos,lpos,hs,hs2
integer ::isF(nclust),chr_limits(nchr,2),IT(nclust,nclust),posi(npos)
integer ::estimateR,onerate,niter,iter
real(dp) ::pemission(npos,2),as(nclust),Fs(nclust)
!real*8, allocatable ::alpha(nclust,npos),scaling(npos),trans(nclust,nclust),pinit(nclust)
real(dp) ::alpha(nclust,npos),scaling(npos),trans(nclust,nclust),T2(nclust,nclust),pinit(nclust)
real(dp) ::beta(nclust,npos),gamma(nclust,npos),num_pi(nclust)
real(dp) ::sumF(2),F,a,r,gr,loglik,loglik0
real(dp), parameter ::Morgan=100000000.d0
real(dp) ::lseg,fromitoj,minrate,maxrate,convf,convthr
real(dp) ::global(nclust),fromi(nclust),ini(nclust),mux(nclust),EL(nclust)

do iter=1,niter

alpha=0.d0;beta=0.d0;gamma=0.d0;scaling=0.d0;pinit=Fs
loglik=0.d0;num_pi=0.d0;global=0.d0;EL=0.d0
IT=1;isF=1;isF(nclust)=0

!############ FORWARD ALGORITHM ####################

do chr=1,nchr

fpos=chr_limits(chr,1);lpos=chr_limits(chr,2)

! initialisation

do i=1,nclust
  alpha(i,fpos)=pinit(i)*pemission(fpos,isF(i)+1)
  scaling(fpos)=scaling(fpos)+alpha(i,fpos)
enddo

 scaling(fpos)=1.0/scaling(fpos)
 alpha(:,fpos)=alpha(:,fpos)*scaling(fpos)

! induction: to get to i, two ways:
! 1) transition + jump into cluster i
! 2) no transition and in i at previous position
 
 do k=fpos+1,lpos
  trans=TM(nclust,IT,as,Fs,k,npos,posi)
  do i=1,nclust
   do j=1,nclust
!    alpha(i,k)=alpha(i,k)+alpha(j,k-1)*trans(j,i)*pemission(k,isF(i)+1)
    alpha(i,k)=alpha(i,k)+alpha(j,k-1)*trans(j,i)
   enddo
   alpha(i,k)=alpha(i,k)*pemission(k,isF(i)+1)
   scaling(k)=scaling(k)+alpha(i,k)
  enddo
 scaling(k)=1.0/scaling(k)
 alpha(:,k)=alpha(:,k)*scaling(k)
 enddo  


!############ BACKWARD ALGORITHM ####################


! initialisation

do i=1,nclust
 gamma(i,lpos)=alpha(i,lpos)*1.0  ! beta(i,j,npos)=1.0
 beta(i,lpos)=1.0*scaling(lpos)
enddo

! induction
! to arrive in k: with or without transition

do k=lpos-1,fpos,-1
 trans=TM(nclust,IT,as,Fs,k+1,npos,posi)
 do i=1,nclust
  do j=1,nclust
    beta(i,k)=beta(i,k)+trans(i,j)*pemission(k+1,isF(j)+1)*beta(j,k+1)
  enddo
  beta(i,k)=beta(i,k)*scaling(k)
  gamma(i,k)=alpha(i,k)*beta(i,k)/scaling(k)  
 enddo
enddo

do k=fpos,lpos
 gamma(:,k)=gamma(:,k)/sum(gamma(:,k))
enddo

! termination

 loglik=loglik-sum(log(scaling(fpos:lpos)))
 
!########## estimation transitions #############

do i=1,nclust
 num_pi(i)=num_pi(i)+gamma(i,fpos)
enddo

!######### mixing coefficients only ########
if(estimateR==0)then
 do k=fpos,lpos-1
  T2=pexit(nclust,IT,as,Fs,k+1,npos,posi)
  do i=1,nclust
   do j=1,nclust
    num_pi(j)=num_pi(j)+alpha(i,k)*T2(i,j)*pemission(k+1,isF(j)+1)*beta(j,k+1)
   enddo
  enddo
 enddo
endif

!######### mixing coefficients and rates ########
if(estimateR==1)then
 do k=fpos,lpos-1
  T2=pexit(nclust,IT,as,Fs,k+1,npos,posi)
  trans=TM(nclust,IT,as,Fs,k+1,npos,posi)
  fromi=0.d0;ini=0.d0
  lseg=(posi(k+1)-posi(k))/Morgan
  do i=1,nclust
   if(lseg>0.d0)mux(i)=1/as(i)-lseg/(dexp(as(i)*lseg)-1.d0)
   do j=1,nclust
    fromitoj=alpha(i,k)*T2(i,j)*pemission(k+1,isF(j)+1)*beta(j,k+1)
    num_pi(j)=num_pi(j)+fromitoj
    ini(i)=ini(i)+alpha(i,k)*trans(i,j)*pemission(k+1,isF(j)+1)*beta(j,k+1)
    fromi(i)=fromi(i)+fromitoj
   enddo
  enddo
  do i=1,nclust
   global(i)=global(i)+fromi(i)
   if(lseg>0.d0)EL(i)=EL(i)+(ini(i)-fromi(i))*(1.d0*(posi(k+1)-posi(k))/Morgan)+fromi(i)*mux(i)
!   if(lseg>0.d0)EL(i)=EL(i)+(ini(i)-fromi(i))*lseg+fromi(i)*mux(i)
  enddo
 enddo
endif

enddo ! end loop on chromosomes

if(estimateR==1)then
 if(nclust==2 .and. onerate==1)then ! joint parameter
   as(1)=sum(global)/sum(EL)
    if(as(1)<minrate)as(1)=minrate
   if(as(1)>maxrate)as(1)=maxrate
   as(2)=as(1)
 else
  do i=1,nclust
    as(i)=global(i)/EL(i)
   if(as(i)<minrate)as(i)=minrate
   if(as(i)>maxrate)as(i)=maxrate
  enddo
 endif
endif

convf=1.d0
if(iter > 1)convf=abs(loglik-loglik0)/abs(loglik)
loglik0=loglik

Fs=num_pi/sum(num_pi)

if(convf < convthr)exit

enddo ! end iterations EM

!print*,"Estimation of parameters with EM algorithm ::",iter,loglik,convf

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


