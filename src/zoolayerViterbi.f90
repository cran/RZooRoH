subroutine zoolayerViterbi(nclust,nchr,npos,pemission,chr_limits,zrates,zmix,posi,states)
implicit none
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
integer ::i,j,l,k,nclust,nchr,npos,chr,fpos,lpos,pos(1),phi(nclust,npos)
integer ::isF(nclust),chr_limits(nchr,2),posi(npos),states(npos)
real(dp) ::pemission(npos,2),as(nclust),Fs(nclust),delta(nclust,npos)
real(dp) ::trans(nclust,nclust),pinit(nclust),zmix(nclust-1),zrates(nclust-1)
real(dp) ::pmax,val(nclust),d
real(dp), parameter ::Morgan=100000000.d0
real(dp) ::ptok(nclust),pnhbd(nclust),MTOC((nclust-1),nclust)

Fs(1:(nclust-1))=zmix;Fs(nclust)=1-zmix(nclust-1)
as(1:(nclust-1))=zrates;as(nclust)=zrates(nclust-1)
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
do k=1,(nclust-1)
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

trans=TM(nclust,as,MTOC,k,npos,posi)

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

function TM(K,as,MTOC,snp,nsnps,mpos)
implicit none
integer ::K,snp,nsnps,mpos(nsnps)
real(dp) ::as(K),TM(K,K),pchange(K),MTOC((K-1),K)
real(dp), parameter ::Morgan=100000000.d0

!#### probability that ancestry is changing at level k
!#### but not at upper level

pchange=0.d00
d=(mpos(snp)-mpos(snp-1))/Morgan
 pchange(1)=1.d0-dexp(-as(1)*d)
 do i=2,K-1
!  pchange(i)=exp(-as(i-1)*d)*(1.d0-exp(-(as(i)-as(i-1))*d))
 pchange(i)=exp(-as(i-1)*d)-exp(-as(i)*d)
enddo

!###### now loop over layers (nclust-1) #######
!###### multiply probability of coancestry change in level k
!###### by conditional probability column(MTOC(,k)) #######

TM=0.d0
do l=1,(K-1)
 do j=l,K
  TM(j,:)=TM(j,:)+MTOC(l,:)*pchange(l)
 enddo
enddo

!###### end by adding diagonal elements #####
!###### no coancestry change ########

do j=1,K
  TM(j,j)=TM(j,j)+exp(-as(j)*d)
enddo

end function

end subroutine


