subroutine zoolayerlik(nclust,nchr,npos,pemission,chr_limits,as,Fs,posi,loglik)
implicit none
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
integer ::i,j,l,k,nclust,nchr,npos,chr,fpos,lpos,hs,hs2
integer ::isF(nclust),chr_limits(nchr,2),posi(npos)
real(dp) ::pemission(npos,2),as(nclust),Fs(nclust)
real(dp) ::alpha(nclust,npos),scaling(npos),trans(nclust,nclust),pinit(nclust)
real(dp) ::sumF(2),F,a,r,gr,loglik,d
real(dp), parameter ::Morgan=100000000.d0
real(dp) ::ptok(nclust),pnhbd(nclust),pchange(nclust),MTOC(nclust,nclust)

alpha=0.d0;scaling=0.d0;pinit=Fs
loglik=0.d0;isF=1;isF(nclust)=0

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
!##### after on ancestry change at level l, we can ####
!##### only reach states from level l or higher #######

MTOC=0.d00
do k=1,nclust
 MTOC(k,k:nclust)=ptok(k:nclust)/sum(ptok(k:nclust))
enddo

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

! #### induction: to get to i, two ways:
! #### 1) transition + jump into cluster i
! #### 2) no transition and in i at previous position

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

! #### termination

 loglik=loglik-sum(log(scaling(fpos:lpos)))

enddo ! end chr

end subroutine



