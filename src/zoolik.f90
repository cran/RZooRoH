subroutine zoolik(nclust,nchr,npos,pemission,chr_limits,as,Fs,posi,loglik)
implicit none
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
integer ::i,j,l,k,nclust,nchr,npos,chr,fpos,lpos,hs,hs2
integer ::isF(nclust),chr_limits(nchr,2),IT(nclust,nclust),posi(npos)
real(dp) ::pemission(npos,2),as(nclust),Fs(nclust)
!real(dp), allocatable ::alpha(nclust,npos),scaling(npos),trans(nclust,nclust),pinit(nclust)
real(dp) ::alpha(nclust,npos),scaling(npos),trans(nclust,nclust),pinit(nclust)
real(dp) ::sumF(2),F,a,r,gr,loglik
real(dp), parameter ::Morgan=100000000.d0


!allocate(alpha(nclust,npos),scaling(npos),pinit(nclust),trans(nclust,nclust))

alpha=0.d0;scaling=0.d0;pinit=Fs
loglik=0.d0
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


! compute transition (normally function) 
 
  trans=0.d0
  do hs=1,nclust
    sumF(1)=dot_product(Fs,IT(hs,:))
    a=as(hs)
    r=dexp(-a*(posi(k)-posi(k-1))/Morgan)
    do hs2=1,nclust
     if(IT(hs,hs2)==0)cycle
     trans(hs,hs2)=(1-r)*Fs(Hs2)/sumF(1)
    enddo
    trans(hs,hs)=trans(hs,hs)+r
  enddo



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

! termination

 loglik=loglik-sum(log(scaling(fpos:lpos)))
! if(chr==nchr)print*,loglik,as,Fs

enddo ! end chr

end subroutine


