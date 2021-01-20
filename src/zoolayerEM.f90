subroutine zoolayerEM(nclust,nchr,npos,pemission,chr_limits,as,Fs,posi,loglik,gamma,estimateR,&
           onerate,niter,minrate,maxrate,convthr,iter)
implicit none
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
integer ::i,j,l,k,nclust,nchr,npos,chr,fpos,lpos,hs,hs2
integer ::isF(nclust),chr_limits(nchr,2),IT(nclust,nclust),posi(npos)
integer ::estimateR,onerate,niter,iter
real(dp) ::pemission(npos,2),as(nclust),Fs(nclust)
!real*8, allocatable ::alpha(nclust,npos),scaling(npos),trans(nclust,nclust),pinit(nclust)
real(dp) ::alpha(nclust,npos),scaling(npos),trans(nclust,nclust),T2(nclust,nclust),pinit(nclust)
real(dp) ::beta(nclust,npos),gamma(nclust,npos),num_pi(nclust),denom_pi(nclust)
real(dp) ::sumF(2),F,a,r,gr,loglik,loglik0,d
real(dp), parameter ::Morgan=100000000.d0
real(dp) ::lseg,fromitoj,minrate,maxrate,convf,convthr
real(dp) ::global(nclust),fromi(nclust),ini(nclust),mux(nclust),EL(nclust)
real(dp) ::ptok(nclust),pnhbd(nclust),pchange(nclust),MTOC(nclust,nclust),val

do iter=1,niter

alpha=0.d0;beta=0.d0;gamma=0.d0;scaling=0.d0;pinit=Fs
loglik=0.d0;num_pi=0.d0;denom_pi=0.d0;global=0.d0;EL=0.d0
isF=1;isF(nclust)=0

!############ PREPARE PARAMETERS ###################
!############ FOR TRANSITION PROB ##################

!## At level 1, state 1 is HBD, all other states non-HBD
!## All other states should be considered as one block
!## With same transition probabilities to state 1
!## At level l, state l is HBD, all other states non-HBD

!#### prob to reach state k ########################
!#### after change of ancestry at first level ######

pnhbd=0.d0 !#### probability to reach the non-hbd class at level k
pnhbd(1)=1.d0-Fs(1)
do k=2,(nclust-1) !#### there are nclut-1 levels
 pnhbd(k)=pnhbd(k-1)*(1.d0-Fs(k))
enddo

ptok(1)=Fs(1) !#### probabilitiy to reach state k
do k=2,(nclust-1) !### only if nclust > 2, else old 1R
 ptok(k)=pnhbd(k-1)*Fs(k)
enddo
ptok(nclust)=pnhbd(nclust-1)

!##### these probabilities also correspond to initial prob ####
pinit=ptok

!##### MTOC = matrix of conditional probabilities #####
!##### cond. probabilities to reach state k ###########
!##### conditional on an ancestry change at level l ###
!##### after an ancestry change at level l, we can ####
!##### only reach states from level l or higher #######

MTOC=0.d00
do k=1,nclust
 MTOC(k,k:nclust)=ptok(k:nclust)/sum(ptok(k:nclust))
enddo


!############ FORWARD ALGORITHM ####################

do chr=1,nchr

fpos=chr_limits(chr,1);lpos=chr_limits(chr,2)

! ##### initialisation

do i=1,nclust
  alpha(i,fpos)=pinit(i)*pemission(fpos,isF(i)+1)
  scaling(fpos)=scaling(fpos)+alpha(i,fpos)
enddo

 scaling(fpos)=1.0/scaling(fpos)
 alpha(:,fpos)=alpha(:,fpos)*scaling(fpos)

! ##### induction: to get to i, two ways:
! ##### 1) transition + jump into cluster i
! ##### 2) no transition and in i at previous position

 do k=fpos+1,lpos

! #### compute transition (normally function)

 !#### probability that ancestry is changing at level k
 !#### but not at upper level

  pchange=0.d00
  d=(posi(k)-posi(k-1))/Morgan
  pchange(1)=1.d0-dexp(-as(1)*d)
  do i=2,nclust-1
   pchange(i)=exp(-as(i-1)*d)*(1.d0-exp(-(as(i)-as(i-1))*d))
  enddo

 !###### now loop over layers (nclust-1) #######
 !###### multiply probability of coancestry change in level k
 !###### by conditional probability column(MTOC(,k)) #######

  trans=0.d0
  do l=1,(nclust-1)
   do j=l,nclust
    trans(j,:)=trans(j,:)+MTOC(l,:)*pchange(l)
   enddo
  enddo

 !###### end by adding diagonal elements #####
 !###### no coancestry change ########

  do j=1,nclust
   trans(j,j)=trans(j,j)+exp(-as(j)*d)
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


!############ BACKWARD ALGORITHM ####################


! #### initialisation

do i=1,nclust
 gamma(i,lpos)=alpha(i,lpos)*1.0  ! beta(i,j,npos)=1.0
 beta(i,lpos)=1.0*scaling(lpos)
enddo

! #### induction
! #### to arrive in k: with or without transition

do k=lpos-1,fpos,-1
! trans=transition(k)

 ! #### compute transition (normally function)

 !#### probability that ancestry is changing at level k
 !#### but not at upper level

  pchange=0.d00
  d=(posi(k+1)-posi(k))/Morgan
  pchange(1)=1.d0-dexp(-as(1)*d)
  do i=2,nclust-1
   pchange(i)=exp(-as(i-1)*d)*(1.d0-exp(-(as(i)-as(i-1))*d))
  enddo

 !###### now loop over layers (nclust-1) #######
 !###### multiply probability of coancestry change in level k
 !###### by conditional probability column(MTOC(,k)) #######

  trans=0.d0
  do l=1,(nclust-1)
   do j=l,nclust
    trans(j,:)=trans(j,:)+MTOC(l,:)*pchange(l)
   enddo
  enddo

 !###### end by adding diagonal elements #####
 !###### no coancestry change ########

  do j=1,nclust
   trans(j,j)=trans(j,j)+exp(-as(j)*d)
  enddo


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
! if(chr==nchr)print*,loglik,as,Fs

!########## estimation transitions #############

do i=1,nclust
 num_pi(i)=num_pi(i)+gamma(i,fpos)
 denom_pi(i)=denom_pi(i)+sum(gamma(i:nclust,fpos))
enddo

!######### mixing coefficients only ########
if(estimateR==0)then
 do k=fpos,lpos-1

! #### compute transition (normally function)

 !#### probability that ancestry is changing at level k
 !#### but not at upper level

  pchange=0.d00
  d=(posi(k+1)-posi(k))/Morgan
  pchange(1)=1.d0-dexp(-as(1)*d)
  do i=2,nclust-1
   pchange(i)=exp(-as(i-1)*d)*(1.d0-exp(-(as(i)-as(i-1))*d))
  enddo

  do l=1,nclust-1 !#### number of layers is nclust-1
  do i=l,nclust   !#### after change in layer l, only transitions with states >= l (both in and out)
   do j=l,nclust  !#### after change in layer l, only transitions with states >= l (both in and out)
      val=alpha(i,k)*MTOC(l,j)*pchange(l)*pemission(k+1,isF(j)+1)*beta(j,k+1)
      num_pi(j)=num_pi(j)+val
      denom_pi(l:j)=denom_pi(l:j)+val
!#    num_pi(j)=num_pi(j)+alpha(i,k)*T2(i,j)*pemission(k+1,isF(j)+1)*beta(j,k+1)
   enddo
  enddo
  enddo ! add one additional loop !!!
 enddo
endif

!######### mixing coefficients and rates ########
if(estimateR==1)then
 ! not implemented, too slow
endif

enddo ! end loop on chromosomes

convf=1.d0
if(iter > 1)convf=abs(loglik-loglik0)/abs(loglik)
loglik0=loglik

!#Fs=num_pi/sum(num_pi)
do i=1,nclust-1
 Fs(i)=num_pi(i)/denom_pi(i)
enddo
Fs(nclust)=1-Fs(nclust-1)

if(convf < convthr)exit

enddo ! end iterations EM

!print*,"Estimation of parameters with EM algorithm ::",iter,loglik,convf

end subroutine


