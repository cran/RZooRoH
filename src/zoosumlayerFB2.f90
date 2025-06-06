subroutine zoosumlayerfb2(nclust,nchr,npos,pemission,chr_limits,zrates,zmix,posi,loglik,gamma)
implicit none
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
integer ::i,j,l,k,nclust,nchr,npos,chr,fpos,lpos,hs,hs2
integer ::nG,g1,g2,g3,t,tr,t0,recfun,nd,d,d_i,lev1,lev2,fval,lval,lstep
integer ::isF(nclust),chr_limits(nchr,2),posi(npos)
real(dp) ::pemission(npos,2),as(nclust),Fs(nclust),zmix(nclust-1),zrates(nclust-1)
real(dp) ::alpha(nclust,npos),scaling(npos),pinit(nclust)
real(dp) ::beta(nclust,npos),gamma(nclust,npos),alphar(nclust,npos),scalingr(npos)
real(dp) ::sumF(2),F,a,r,gr,loglik,d2,cst,cst2,prec,loglik2
real(dp), parameter ::Morgan=100000000.d0
real(dp),allocatable ::ptok(:),pnhbd(:),pnorec(:),ctok(:)
real(dp),allocatable ::ctok2(:),precing(:),cumupk(:),upcoal(:),upnocoal(:)
real(dp) ::alpha_up(nclust),alpha_float(nclust)
real(dp), allocatable ::P1(:),P2(:),P3(:,:),P4(:,:),P5(:,:),P6(:,:),P7(:,:)
integer, allocatable ::refdist(:)

recfun=1;nd=550

!######################### make table of rounded recombination value #############

g2=nint(as(nclust-1))
allocate(refdist(nd))

!#### for recombination rates
l=0
do lev1=1,6
 if(lev1==1)fval=100
 lstep=10**(lev1+1)
 if(lev1>1)fval=lval+lstep
 lval=10**(lev1+3)
 do lev2=fval,lval,lstep
  l=l+1
  refdist(l)=lev2
  enddo
enddo

!####################################################################################

!##### as(:) = rates from the clusters (will be the end) - G will go from 1 to X_1, X_1 to X_2 (we still use R and not G)
!##### nclust = number of clusters + 1 for non-HBD, nG = number of generations
!##### Fs(:) = rates in the clusters ...

Fs(1:(nclust-1))=zmix;Fs(nclust)=1-zmix(nclust-1)
as(1:(nclust-1))=zrates;as(nclust)=zrates(nclust-1)

nG=as(nclust-1)
allocate(ptok(nG),pnhbd(nG),ctok(nG),ctok2(nG),pnorec(0:nG))
allocate(precing(nG),cumupk(nG),upcoal(nG),upnocoal(nG))
allocate(P1(nclust),P2(nclust),P3(nclust,nd),P4(nclust,nd),P5(nclust,nd))
allocate(P6(nclust,nd),P7(nclust,nd))

alpha=0.d0;scaling=0.d0;pinit=Fs
alphar=0.d0;scalingr=0.d0;beta=0.d0;gamma=0.d0
loglik=0.d0;isF=1;isF(nclust)=0;loglik2=0.d0

!#### prob to reach state k ########################
!#### after change of ancestry at first level ######

pnhbd=0.d0;ptok=0.d0;ctok=0.d0 !#### probability to reach the non-hbd class at level k
upcoal=0.d0;upnocoal=0.d0;ctok2=0.d0;cumupk=0.d0
pnhbd(1)=1.d0-Fs(1);ptok(1)=Fs(1)
do l=1,(nclust-1) !#### there are nclut-1 levels
 g1=1;if(l>1)g1=nint(as(l-1))+1;g2=nint(as(l))
 do k=g1,g2
  if(k>1)pnhbd(k)=pnhbd(k-1)*(1.d0-Fs(l))
  if(k>1)ptok(k)=pnhbd(k-1)*Fs(l)
  ctok2(k)=Fs(l)*(1.d0-Fs(l))**(k-g1)
 enddo
enddo

do l=1,nclust-1
 g1=1;if(l>1)g1=nint(as(l-1))+1;g2=nint(as(l))
!# conditional probabilities
 P1(l)=(1-Fs(l))**(g2-g1+1)
 P2(l)=1.d0-P1(l) !### could be used to replace the constant (or the second one)
 ctok(g1:g2)=ctok2(g1:g2)/P2(l)
 do k=g1,g2
  cumupk(k)=sum(ctok(k:g2))
  g3=g1+(g2-k) !#### coordinate based on number of generations after recombination in k
  upnocoal(k)=(1-Fs(l))**(g2-k+1)
  upcoal(k)=sum(ctok2(g1:g3)) !### prob. to coalesce in layer if recombination in k
 enddo
enddo
P1(nclust)=0.d0;P2(nclust)=1.d0

!##### initial state prob - stationnary distribution ####
pinit=0
do l=1,nclust-1
 g1=1;if(l>1)g1=nint(as(l-1))+1;g2=nint(as(l))
 pinit(l)=sum(ptok(g1:g2))
enddo
pinit(nclust)=pnhbd(nint(as(nclust-1))) !#### the non-hbd in the last generation from the last layer

!############ COMPUTE VALUES OF P3-P7 in a table ####

P3=0.d0 ! no recombination when in state i before
P4=0.d0 ! recombination and coalescence in state i, when in state i before
P5=0.d0 ! recombination and coalescence in state i, when in state > i before
P6=0.d0 ! recombination and floating in state i, when in state i before
P7=0.d0 ! recombination and floating in state i, when in state > i before [recombination in state i when i state i > before - rec & coalescen in i when i > before = P5]

do d=1,nd
  !####### probability no recombination until layer l #######
  !####### probability recombination in layer l, is then Norec(l+1)-Norec(l)

  pnorec=1.d0

  d2=1.d0*refdist(d)/Morgan
  if(recfun==0)then
   prec=d2**2
   precing(1)=prec
   do i=1,nG
     pnorec(i)=(1.d0-d2)**(2*i)
     if(i>1)precing(i)=prec*pnorec(i-1)
   enddo
  endif
  if(recfun==1)then
     prec=1-dexp(-2*d2)
     precing(1)=prec
   do i=1,nG
!#   pnorec(i)=dexp(-as(i)*d2)
     pnorec(i)=dexp(-2*i*d2)
     if(i>1)precing(i)=prec*pnorec(i-1)
   enddo
  endif

 do l=1,nclust-1
  g1=1;if(l>1)g1=nint(as(l-1))+1;g2=nint(as(l))
  P3(l,d)=dot_product(ctok(g1:g2),pnorec(g1:g2)) ! condition prob. generation g x norec till g
 enddo
 g2=nint(as(nclust-1));P3(nclust,d)=pnorec(g2)

  do l=1,nclust-1
   g1=1;if(l>1)g1=nint(as(l-1))+1;g2=nint(as(l))
   do g3=g1,g2
    P4(l,d)=P4(l,d)+cumupk(g3)*precing(g3)*upcoal(g3)
    P5(l,d)=P5(l,d)+precing(g3)*upcoal(g3) !#### enregistrer?
    P6(l,d)=P6(l,d)+cumupk(g3)*precing(g3)*upnocoal(g3)
   enddo
  enddo
!  g2=nint(as(nclust-1));P4(nclust,d)=pnorec(g2-1)*prec*(1-Fs(nclust-1)) ## double counting
  P5(nclust,d)=0.d0;P6(nclust,d)=0.d0

  do l=1,nclust-1
   g1=1;if(l>1)g1=nint(as(l-1))+1;g2=nint(as(l))
   cst=pnorec(g1-1)-pnorec(g2) !#### recombination between g1 and g2
   P7(l,d)=cst-P5(l,d)
  enddo
  P7(nclust,d)=0.d0

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

 do k=fpos+1,lpos

  !####### probability no recombination until layer l #######
  !####### probability recombination in layer l, is then Norec(l+1)-Norec(l)
  
  d2=1.d0*(posi(k)-posi(k-1))
  d_i=dist_index(d2)

! #### linear algorithm need to compute alpha, alpha_up and alpha_float
  alpha_up=0.d0;alpha_float=0.d0

  !#### we can compute recursively alpha_up (corresponds to k-1 !)
  !alpha_up(nclust-2)=alpha(nclust-1,k-1)+alpha(nclust,k-1) !#### probability to be > k-2 (number of layers = k-1)
  !if(nclust > 3)then
   do l=(nclust-1),1,-1
    alpha_up(l)=alpha_up(l+1)+alpha(l+1,k-1)
   enddo
  !endif

  !#### we can compute alpha_float (corresponds to k-1 !)
  alpha_float(1)=alpha_up(1)*P7(1,d_i)+alpha(1,k-1)*P6(1,d_i)
  do l=2,(nclust-1)  !##### no need for layer K (can not be floating in K+1)
   alpha_float(l)=alpha_float(l-1)*P1(l)+alpha(l,k-1)*P6(l,d_i)+alpha_up(l)*P7(l,d_i) !##### OK for last layer too?
  enddo
  alpha_float(nclust)=0.d0

  !##### define rec vector, with rec0=1
  !##### define float(0) = 0.d0
  !##### also compute for last layer non-HBD

  !#### forward variables (rescale after)
  do l=1,nclust ! should work for non-hbd too
   alpha(l,k)=alpha_float(l-1)*P2(l)+alpha_up(l)*P5(l,d_i)+alpha(l,k-1)*P4(l,d_i)
   alpha(l,k)=alpha(l,k)+alpha(l,k-1)*P3(l,d_i)
   alpha(l,k)=alpha(l,k)*pemission(k,isF(l)+1)
   scaling(k)=scaling(k)+alpha(l,k)
  enddo

  scaling(k)=1.0/scaling(k)
  alpha(:,k)=alpha(:,k)*scaling(k)

 enddo


! #### termination

 loglik=loglik-sum(log(scaling(fpos:lpos)))
! print*,loglik

!############ BACKWARD ALGORITHM ####################
!######### REVERSE FORWARD ##########################


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

  d2=1.d0*(posi(k+1)-posi(k))
  d_i=dist_index(d2)

! #### linear algorithm need to compute alpha, alpha_up and alpha_float
  alpha_up=0.d0;alpha_float=0.d0

  !#### we can compute recursively alpha_up (corresponds to k-1 !)
  !alpha_up(nclust-2)=alphar(nclust-1,k+1)+alphar(nclust,k+1) !#### probability to be > k-2 (number of layers = k-1)
  !if(nclust > 3)then
   do l=(nclust-1),1,-1
    alpha_up(l)=alpha_up(l+1)+alphar(l+1,k+1)
   enddo
  !endif

  !#### we can compute alpha_float (corresponds to k-1 !)
  alpha_float(1)=alpha_up(1)*P7(1,d_i)+alphar(1,k+1)*P6(1,d_i)
  do l=2,(nclust-1)  !##### no need for layer K (can not be floating in K+1)
   alpha_float(l)=alpha_float(l-1)*P1(l)+alphar(l,k+1)*P6(l,d_i)+alpha_up(l)*P7(l,d_i) !##### OK for last layer too?
  enddo
  alpha_float(nclust)=0.d0

  !##### define rec vector, with rec0=1
  !##### define float(0) = 0.d0
  !##### also compute for last layer non-HBD

  !#### forward variables (rescale after)
  do l=1,nclust ! should work for non-hbd too
   alphar(l,k)=alpha_float(l-1)*P2(l)+alpha_up(l)*P5(l,d_i)+alphar(l,k+1)*P4(l,d_i)
   alphar(l,k)=alphar(l,k)+alphar(l,k+1)*P3(l,d_i)
   alphar(l,k)=alphar(l,k)*pemission(k,isF(l)+1)
   scalingr(k)=scalingr(k)+alphar(l,k)
  enddo

  scalingr(k)=1.0/scalingr(k)
  alphar(:,k)=alphar(:,k)*scalingr(k)

 enddo

 loglik2=loglik2-sum(log(scalingr(fpos:lpos)))

! #### termination

 do k=fpos,lpos
  do l=1,nclust
   beta(l,k)=alphar(l,k)/(pemission(k,isF(l)+1)*pinit(l))
   gamma(l,k)=alpha(l,k)*beta(l,k)/scaling(k)
  enddo
  gamma(:,k)=gamma(:,k)/sum(gamma(:,k))
 enddo

enddo ! end chr

contains

integer function dist_index(d3)
real(dp) ::d3
integer ::i3

if(d3 <= 10000)then !### from 1 to 100
 d3=d3/100.d0
 i3=nint(d3) !### goes from 0 to 100
 if(i3==0)i3=1
else if(d3 <= 100000)then
 d3=d3/1000.d0
 i3=nint(d3) !### goes from 10 to 100 => from 100 (10k) to 190
 i3=(i3-10)+100
else if(d3 <= 1000000)then
 d3=d3/10000.d0
 i3=nint(d3) !### goes from 10 to 100 => from 190 (100k) to 280
 i3=(i3-10)+190
else if(d3 <= 10000000)then
 d3=d3/100000.d0
 i3=nint(d3) !### goes from 10 to 100 => from 280 (1M) to 370
 i3=(i3-10)+280
else if(d3 <= 100000000)then
 d3=d3/1000000.d0
 i3=nint(d3) !### goes from 10 to 100 => from 370 (10M) to 460
 i3=(i3-10)+370
else if(d3 <= 1000000000)then
 d3=d3/1000000.d0
 i3=nint(d3) !### goes from 10 to 100 => from 460 (100M) to 550
 i3=(i3-10)+460
endif

dist_index=i3
return

end function 

end subroutine






