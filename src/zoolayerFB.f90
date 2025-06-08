subroutine zoolayerFB(nclust,nchr,npos,pemission,chr_limits,zrates,zmix,posi,loglik,gamma)
implicit none
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
integer ::i,l,k,nclust,nchr,npos,chr,fpos,lpos
integer ::isF(nclust),chr_limits(nchr,2),posi(npos)
real(dp) ::pemission(npos,2),as(nclust),Fs(nclust),zmix(nclust-1),zrates(nclust-1)
!real(dp), allocatable ::alpha(nclust,npos),scaling(npos),trans(nclust,nclust),pinit(nclust)
real(dp) ::alpha(nclust,npos),scaling(npos),pinit(nclust)
real(dp) ::beta(nclust,npos),gamma(nclust,npos),alphar(nclust,npos),scalingr(npos)
real(dp) ::alpha_up(nclust),alpha_float(0:nclust)
real(dp) ::loglik,d,cst,loglik2
real(dp), parameter ::Morgan=100000000.d0
real(dp) ::ptok(nclust),pnhbd(nclust),pnorec(0:(nclust-1))

Fs(1:(nclust-1))=zmix;Fs(nclust)=1-zmix(nclust-1)
where(Fs< 1e-16)Fs=1e-16
as(1:(nclust-1))=zrates;as(nclust)=zrates(nclust-1)
pinit=Fs

alpha=0.d0;scaling=0.d0;alphar=0.d0;scalingr=0.d0;beta=0.d0;gamma=0.d0
loglik=0.d0;isF=1;isF(nclust)=0;loglik2=0.d0

pnhbd=0.d0 !#### probability to reach the non-hbd class at level k
pnhbd(1)=1.d0-Fs(1)
do k=2,(nclust-1) !#### there are nclust-1 levels
 pnhbd(k)=pnhbd(k-1)*(1.d0-Fs(k))
enddo

ptok(1)=Fs(1) !#### probabilitiy to reach state k
do k=2,(nclust-1) !### only if nclust > 2, else old 1R
 ptok(k)=pnhbd(k-1)*Fs(k)
enddo
ptok(nclust)=pnhbd(nclust-1)

!##### these probabilities also correspond to initial prob ####
 pinit=ptok

!############ FORWARD ALGORITHM ####################


do chr=1,nchr

fpos=chr_limits(chr,1);lpos=chr_limits(chr,2)

! #### initialisation

do i=1,nclust
  alpha(i,fpos)=pinit(i)*pemission(fpos,isF(i)+1)
  scaling(fpos)=scaling(fpos)+alpha(i,fpos)
enddo

scaling(fpos)=1.0/scaling(fpos)
alpha(:,fpos)=alpha(:,fpos)*scaling(fpos)

do k=fpos+1,lpos

!####### probability no recombination until layer l #######
!####### probability recombination in layer l, is then Norec(l+1)-Norec(l)

d=(posi(k)-posi(k-1))/Morgan
pnorec=1.d0
do i=1,nclust-1
 pnorec(i)=dexp(-as(i)*d)
enddo

! #### linear algorithm need to compute alpha, alpha_up and alpha_float
alpha_up=0.d0;alpha_float=0.d0

!#### we can compute recursively alpha_up (corresponds to k-1 !)
if(nclust>2)alpha_up(nclust-2)=alpha(nclust-1,k-1)+alpha(nclust,k-1) !#### probability to be > k-2 (number of layers = k-1)
if(nclust > 3)then
 do l=(nclust-3),1,-1
  alpha_up(l)=alpha_up(l+1)+alpha(l+1,k-1)
 enddo
endif

!#### we can compute alpha_float (corresponds to k-1 !)
alpha_float(1)=(alpha_up(1)+alpha(1,k-1))*(pnorec(0)-pnorec(1))*(1.d0-Fs(1))
do l=2,(nclust-2)  !##### no need for layer K (can not be floating in K+1)
 alpha_float(l)=(alpha_float(l-1)+(alpha(l,k-1)+alpha_up(l))*(pnorec(l-1)-pnorec(l)))*(1.d0-Fs(l)) !##### OK for last layer too?
enddo

!##### define rec vector, with rec0=1
!##### define float(0) = 0.d0
!##### also compute for last layer non-HBD
!#### forward variables (rescale after)
do l=1,nclust-2
 alpha(l,k)=alpha_float(l-1)*Fs(l)+(alpha_up(l)+alpha(l,k-1))*(pnorec(l-1)-pnorec(l))*Fs(l)
 alpha(l,k)=alpha(l,k)+alpha(l,k-1)*pnorec(l)
 alpha(l,k)=alpha(l,k)*pemission(k,isF(l)+1)
 scaling(k)=scaling(k)+alpha(l,k)
enddo

!### for last layer do manually
l=nclust-1
cst=alpha_float(l-1)+(alpha_up(l)+alpha(l,k-1)+alpha(l+1,k-1))*(pnorec(l-1)-pnorec(l)) !### prob in layer l is alpha(l)+alpha(l+1)
alpha(l,k)=(cst*Fs(l)+alpha(l,k-1)*pnorec(l))*pemission(k,isF(l)+1) !#### HBD in layer l
alpha(l+1,k)=(cst*Fs(l+1)+alpha(l+1,k-1)*pnorec(l))*pemission(k,isF(l+1)+1) !#### Non-HBD in layer l (state l+1)
scaling(k)=scaling(k)+alpha(l,k)+alpha(l+1,k)

scaling(k)=1.0/scaling(k)
alpha(:,k)=alpha(:,k)*scaling(k)

enddo

loglik=loglik-sum(log(scaling(fpos:lpos)))

!############ BACKWARD ALGORITHM ####################
!############ Flip the FORWARD ######################

fpos=chr_limits(chr,1);lpos=chr_limits(chr,2)

! #### initialisation

do i=1,nclust
  alphar(i,lpos)=pinit(i)*pemission(lpos,isF(i)+1)
  scalingr(lpos)=scalingr(lpos)+alphar(i,lpos)
enddo

scalingr(lpos)=1.0/scalingr(lpos)
alphar(:,lpos)=alphar(:,lpos)*scalingr(lpos)

do k=lpos-1,fpos,-1

!####### probability no recombination until layer l #######
!####### probability recombination in layer l, is then Norec(l+1)-Norec(l)

d=(posi(k+1)-posi(k))/Morgan
 pnorec=1.d0
 do i=1,nclust-1
  pnorec(i)=dexp(-as(i)*d)
 enddo

! #### linear algorithm need to compute alpha, alpha_up and alpha_float
alpha_up=0.d0;alpha_float=0.d0

!#### we can compute recursively alpha_up (corresponds to k-1 !)
if(nclust>2)alpha_up(nclust-2)=alphar(nclust-1,k+1)+alphar(nclust,k+1) !#### probability to be > k-2 (number of layers = k-1)
if(nclust > 3)then
 do l=(nclust-3),1,-1
  alpha_up(l)=alpha_up(l+1)+alphar(l+1,k+1)
 enddo
endif

!#### we can compute alpha_float (corresponds to k-1 !)
alpha_float(1)=(alpha_up(1)+alphar(1,k+1))*(pnorec(0)-pnorec(1))*(1.d0-Fs(1))
do l=2,(nclust-2)  !##### no need for layer K (can not be floating in K+1)
 alpha_float(l)=(alpha_float(l-1)+(alphar(l,k+1)+alpha_up(l))*(pnorec(l-1)-pnorec(l)))*(1.d0-Fs(l)) !##### OK for last layer too?
enddo

!##### define rec vector, with rec0=1
!##### define float(0) = 0.d0
!##### also compute for last layer non-HBD

!#### forward variables (rescale after)
do l=1,nclust-2
 alphar(l,k)=alpha_float(l-1)*Fs(l)+(alpha_up(l)+alphar(l,k+1))*(pnorec(l-1)-pnorec(l))*Fs(l)
 alphar(l,k)=alphar(l,k)+alphar(l,k+1)*pnorec(l)
 alphar(l,k)=alphar(l,k)*pemission(k,isF(l)+1)
 scalingr(k)=scalingr(k)+alphar(l,k)
enddo

!### for last layer do manually
l=nclust-1
cst=alpha_float(l-1)+(alpha_up(l)+alphar(l,k+1)+alphar(l+1,k+1))*(pnorec(l-1)-pnorec(l)) !### prob in layer l is alpha(l)+alpha(l+1)
alphar(l,k)=(cst*Fs(l)+alphar(l,k+1)*pnorec(l))*pemission(k,isF(l)+1) !#### HBD in layer l
alphar(l+1,k)=(cst*Fs(l+1)+alphar(l+1,k+1)*pnorec(l))*pemission(k,isF(l+1)+1) !#### Non-HBD in layer l (state l+1)
scalingr(k)=scalingr(k)+alphar(l,k)+alphar(l+1,k)

scalingr(k)=1.0/scalingr(k)
alphar(:,k)=alphar(:,k)*scalingr(k)

enddo

loglik2=loglik2-sum(log(scalingr(fpos:lpos)))


do k=fpos,lpos
 do l=1,nclust
  if(pinit(l)>1e-16)then
    beta(l,k)=alphar(l,k)/(pemission(k,isF(l)+1)*pinit(l))
  else
    beta(l,k)=1e-16
  endif
  gamma(l,k)=alpha(l,k)*beta(l,k)/scaling(k)
 enddo
 gamma(:,k)=gamma(:,k)/sum(gamma(:,k))
enddo

enddo ! end chr

contains

end subroutine


