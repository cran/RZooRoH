subroutine zoosumlayerfb(nclust,nchr,npos,pemission,chr_limits,zrates,zmix,posi,loglik,gamma)
implicit none
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
integer ::i,j,l,k,nclust,nchr,npos,chr,fpos,lpos,hs,hs2
integer ::nG,g1,g2,g3,t,tr,t0,recfun
integer ::isF(nclust),chr_limits(nchr,2),posi(npos)
real(dp) ::pemission(npos,2),as(nclust),Fs(nclust),zmix(nclust-1),zrates(nclust-1)
real(dp) ::alpha(nclust,npos),scaling(npos),pinit(nclust)
real(dp) ::beta(nclust,npos),gamma(nclust,npos),alphar(nclust,npos),scalingr(npos)
real(dp) ::sumF(2),F,a,r,gr,loglik,d,cst,cst2,prec,loglik2
real(dp), parameter ::Morgan=100000000.d0
real(dp),allocatable ::ptok(:),pnhbd(:),pnorec(:),ctok(:)
real(dp),allocatable ::ctok2(:),precing(:),cumupk(:),upcoal(:),upnocoal(:)
real(dp) ::alpha_up(nclust),alpha_float(nclust)
real(dp), allocatable ::P1(:),P2(:),P3(:),P4(:),P5(:),P6(:),P7(:)

recfun=1

!##### as(:) = rates from the clusters (will be the end) - G will go from 1 to X_1, X_1 to X_2 (we still use R and not G)
!##### nclust = number of clusters + 1 for non-HBD, nG = number of generations
!##### Fs(:) = rates in the clusters ...

Fs(1:(nclust-1))=zmix;Fs(nclust)=1-zmix(nclust-1)
as(1:(nclust-1))=zrates;as(nclust)=zrates(nclust-1)

nG=as(nclust-1)
allocate(ptok(nG),pnhbd(nG),ctok(nG),ctok2(nG),pnorec(0:nG))
allocate(precing(nG),cumupk(nG),upcoal(nG),upnocoal(nG))
allocate(P1(nclust),P2(nclust),P3(nclust),P4(nclust),P5(nclust))
allocate(P6(nclust),P7(nclust))

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
  if(recfun==0)then
   prec=d**2
   precing(1)=prec
   do i=1,nG
     pnorec(i)=(1.d0-d)**(2.d0*i)
     if(i>1)precing(i)=prec*pnorec(i-1)
   enddo
  endif
  if(recfun==1)then
     prec=1-dexp(-2.d0*d)
     precing(1)=prec
   do i=1,nG
!#   pnorec(i)=dexp(-as(i)*d)
     pnorec(i)=dexp(-i*2.d0*d)
     if(i>1)precing(i)=prec*pnorec(i-1)
   enddo

  endif

  !##### compute P3, P4, P5 and P6, based on recombination rates ######

  P3=0.d0 ! no recombination when in state i before
  do l=1,nclust-1
   g1=1;if(l>1)g1=nint(as(l-1))+1;g2=nint(as(l))
   P3(l)=dot_product(ctok(g1:g2),pnorec(g1:g2)) ! condition prob. generation g x norec till g
  enddo
  g2=nint(as(nclust-1));P3(nclust)=pnorec(g2)

  P4=0.d0 ! recombination and coalescence in state i, when in state i before
  P5=0.d0 ! recombination and coalescence in state i, when in state > i before
  P6=0.d0 ! recombination and floating in state i, when in state i before
  do l=1,nclust-1
   g1=1;if(l>1)g1=nint(as(l-1))+1;g2=nint(as(l))
   do g3=g1,g2
    P4(l)=P4(l)+cumupk(g3)*precing(g3)*upcoal(g3)
    P5(l)=P5(l)+precing(g3)*upcoal(g3) !#### enregistrer?
    P6(l)=P6(l)+cumupk(g3)*precing(g3)*upnocoal(g3)
   enddo
  enddo
!  g2=nint(as(nclust-1));P4(nclust)=pnorec(g2-1)*prec*(1-Fs(nclust-1)) ### double counting
  P5(nclust)=0.d0;P6(nclust)=0.d0

  P7=0.d0 ! recombination and floating in state i, when in state > i before [recombination in state i when i state i > before - rec & coalescen in i when i > before = P5]
  do l=1,nclust-1
   g1=1;if(l>1)g1=nint(as(l-1))+1;g2=nint(as(l))
   cst=pnorec(g1-1)-pnorec(g2) !#### recombination between g1 and g2
   P7(l)=cst-P5(l)
  enddo
  P7(nclust)=0.d0

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
  alpha_float(1)=alpha_up(1)*P7(1)+alpha(1,k-1)*P6(1)
  do l=2,(nclust-1)  !##### no need for layer K (can not be floating in K+1)
   alpha_float(l)=alpha_float(l-1)*P1(l)+alpha(l,k-1)*P6(l)+alpha_up(l)*P7(l) !##### OK for last layer too?
  enddo
  alpha_float(nclust)=0.d0

  !##### define rec vector, with rec0=1
  !##### define float(0) = 0.d0
  !##### also compute for last layer non-HBD

  !#### forward variables (rescale after)
  do l=1,nclust ! should work for non-hbd too
   alpha(l,k)=alpha_float(l-1)*P2(l)+alpha_up(l)*P5(l)+alpha(l,k-1)*P4(l)
   alpha(l,k)=alpha(l,k)+alpha(l,k-1)*P3(l)
   alpha(l,k)=alpha(l,k)*pemission(k,isF(l)+1)
   scaling(k)=scaling(k)+alpha(l,k)
  enddo

  scaling(k)=1.0/scaling(k)
  alpha(:,k)=alpha(:,k)*scaling(k)

 enddo

! #### termination

 loglik=loglik-sum(log(scaling(fpos:lpos)))

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

  d=(posi(k+1)-posi(k))/Morgan
  pnorec=1.d0
  if(recfun==0)then
   prec=d**2
   precing(1)=prec
   do i=1,nG
     pnorec(i)=(1.d0-d)**(2.d0*i)
     if(i>1)precing(i)=prec*precing(i-1)
   enddo
  endif
  if(recfun==1)then
     prec=1-dexp(-2.d0*d)
     precing(1)=prec
   do i=1,nG
!#   pnorec(i)=dexp(-as(i)*d)
     pnorec(i)=dexp(-i*2.d0*d)
     if(i>1)precing(i)=prec*pnorec(i-1)
   enddo

  endif

  !##### compute P3, P4, P5 and P6, based on recombination rates ######

  P3=0.d0 ! no recombination when in state i before
  do l=1,nclust-1
   g1=1;if(l>1)g1=nint(as(l-1))+1;g2=nint(as(l))
   P3(l)=dot_product(ctok(g1:g2),pnorec(g1:g2)) ! condition prob. generation g x norec till g
  enddo
  g2=nint(as(nclust-1));P3(nclust)=pnorec(g2)

  P4=0.d0 ! recombination and coalescence in state i, when in state i before
  P5=0.d0 ! recombination and coalescence in state i, when in state > i before
  P6=0.d0 ! recombination and floating in state i, when in state i before
  do l=1,nclust-1
   g1=1;if(l>1)g1=nint(as(l-1))+1;g2=nint(as(l))
   do g3=g1,g2
    P4(l)=P4(l)+cumupk(g3)*precing(g3)*upcoal(g3)
    P5(l)=P5(l)+precing(g3)*upcoal(g3) !#### enregistrer?
    P6(l)=P6(l)+cumupk(g3)*precing(g3)*upnocoal(g3)
   enddo
  enddo
!  g2=nint(as(nclust-1));P4(nclust)=pnorec(g2-1)*prec*(1-Fs(nclust-1)) ## double counting
  P5(nclust)=0.d0;P6(nclust)=0.d0

  P7=0.d0 ! recombination and floating in state i, when in state > i before [recombination in state i when i state i > before - rec & coalescen in i when i > before = P5]
  do l=1,nclust-1
   g1=1;if(l>1)g1=nint(as(l-1))+1;g2=nint(as(l))
   cst=pnorec(g1-1)-pnorec(g2) !#### recombination between g1 and g2
   P7(l)=cst-P5(l)
  enddo
  P7(nclust)=0.d0

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
  alpha_float(1)=alpha_up(1)*P7(1)+alphar(1,k+1)*P6(1)
  do l=2,(nclust-1)  !##### no need for layer K (can not be floating in K+1)
   alpha_float(l)=alpha_float(l-1)*P1(l)+alphar(l,k+1)*P6(l)+alpha_up(l)*P7(l) !##### OK for last layer too?
  enddo
  alpha_float(nclust)=0.d0

  !##### define rec vector, with rec0=1
  !##### define float(0) = 0.d0
  !##### also compute for last layer non-HBD

  !#### forward variables (rescale after)
  do l=1,nclust ! should work for non-hbd too
   alphar(l,k)=alpha_float(l-1)*P2(l)+alpha_up(l)*P5(l)+alphar(l,k+1)*P4(l)
   alphar(l,k)=alphar(l,k)+alphar(l,k+1)*P3(l)
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

end subroutine






