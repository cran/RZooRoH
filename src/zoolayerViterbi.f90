subroutine zoolayerViterbi(nclust,nchr,npos,pemission,chr_limits,as,Fs,posi,states)
implicit none
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
integer ::i,j,l,k,nclust,nchr,npos,chr,fpos,lpos,hs,hs2,pos(1),phi(nclust,npos)
integer ::isF(nclust),chr_limits(nchr,2),IT(nclust,nclust),posi(npos),states(npos)
real(dp) ::pemission(npos,2),as(nclust),Fs(nclust),delta(nclust,npos)
real(dp) ::trans(nclust,nclust),pinit(nclust)
real(dp) ::sumF(2),F,a,r,gr,loglik,pmax,val(nclust),mind,d
real(dp), parameter ::Morgan=100000000.d0
real(dp) ::ptok(nclust),pnhbd(nclust),pchange(nclust),MTOC(nclust,nclust)


pinit=Fs
!delta=-1000000000.0;phi=0;states=0
phi=0;states=0
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
!##### after on ancestry change at level l, we can ####
!##### only reach states from level l or higher #######

MTOC=0.d00
do k=1,nclust
 MTOC(k,k:nclust)=ptok(k:nclust)/sum(ptok(k:nclust))
enddo

! ########
do chr=1,nchr

fpos=chr_limits(chr,1);lpos=chr_limits(chr,2)

! initialisation

 do i=1,nclust
  delta(i,fpos)=log(pinit(i))+log(pemission(fpos,isF(i)+1))
  phi(i,fpos)=0
 enddo

! recursion

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


