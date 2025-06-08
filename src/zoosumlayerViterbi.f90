subroutine zoosumlayerViterbi(nclust,nchr,npos,pemission,chr_limits,zrates,zmix,posi,states)
implicit none
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
integer ::i,j,l,k,nclust,nchr,npos,chr,fpos,lpos,pos(1),phi(nclust,npos)
integer ::isF(nclust),chr_limits(nchr,2),posi(npos),states(npos)
integer ::l2,l3,nG,g1,g2,g3,recfun
real(dp) ::pemission(npos,2),as(nclust),Fs(nclust),delta(nclust,npos)
real(dp) ::trans(nclust,nclust),pinit(nclust),zmix(nclust-1),zrates(nclust-1)
real(dp) ::pmax,val(nclust),cst,prec,d
real(dp), parameter ::Morgan=100000000.d0
real(dp),allocatable ::ptok(:),pnhbd(:),pnorec(:),ctok(:)
real(dp),allocatable ::ctok2(:),precing(:),cumupk(:),upcoal(:),upnocoal(:)
real(dp), allocatable ::P1(:),P2(:),P3(:),P4(:),P5(:),P6(:),P7(:)

recfun=1

Fs(1:(nclust-1))=zmix;Fs(nclust)=1-zmix(nclust-1)
as(1:(nclust-1))=zrates;as(nclust)=zrates(nclust-1)
!delta=-1000000000.0;phi=0;states=0
phi=0;states=0
isF=1;isF(nclust)=0

nG=nint(as(nclust-1))
allocate(ptok(nG),pnhbd(nG),ctok(nG),ctok2(nG),pnorec(0:nG))
allocate(precing(nG),cumupk(nG),upcoal(nG),upnocoal(nG))
allocate(P1(nclust),P2(nclust),P3(nclust),P4(nclust),P5(nclust))
allocate(P6(nclust),P7(nclust))

!############ PREPARE PARAMETERS ###################

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
!  g2=nint(as(nclust-1));P4(nclust)=pnorec(g2-1)*prec*(1-Fs(nclust-1)) ### double counting?
  P5(nclust)=0.d0;P6(nclust)=0.d0

  P7=0.d0 ! recombination and floating in state i, when in state > i before [recombination in state i when i state i > before - rec & coalescen in i when i > before = P5]
  do l=1,nclust-1
   g1=1;if(l>1)g1=nint(as(l-1))+1;g2=nint(as(l))
   cst=pnorec(g1-1)-pnorec(g2) !#### recombination between g1 and g2
   P7(l)=cst-P5(l)
  enddo
  P7(nclust)=0.d0

 !#### computation of transition matrix
  trans=0.d0
  do i=1,nclust-1 ! non-HBD is a "special" case (for transition we can cut prob to non-HBD as 1 - sum(to HBD classes)
   do j=1,nclust-1 
     l=min(i,j) !### recombination must be before layer i and j
     do l2=1,l-1 !#### loop over all the layers before i and j
       cst=P7(l2) !### no recombination till the layer x rec in the layer + floating
       do l3=(l2+1),(j-1)
        cst=cst*P1(l3) !#### keep floating in all layers till layer j-1
      enddo
      trans(i,j)=trans(i,j)+cst*P2(j) !##### multiply by coalescence in layer j
     enddo
     !##### if i < j, possibility to recombine in i when previously in i
     if(i < j)then
       cst=P6(i)
       do l3=(i+1),(j-1)
        cst=cst*P1(l3) !#### keep floating in all layers till layer j-1
       enddo
       trans(i,j)=trans(i,j)+cst*P2(j) !##### multiply by coalescence in layer j
     endif
     !##### if i > j, can recombine in layers j and coalesce in j when in higher layer
     if(i > j)trans(i,j)=trans(i,j)+P5(j)
     !##### if i == j, recombine in layer j, coalesce in j when previously in j
     if(i == j)trans(i,j)=trans(i,j)+P4(j)
     !##### if i == j, non-recombination till j when previously in j
     if(i == j)trans(i,j)=trans(i,j)+P3(j)
   enddo
 enddo
 
 !####### transition from non-hbd
 do j=1,nclust-1
  !#### recombination in a layer x keep coalesce x coalesce in j
  do l2=1,j-1
   cst=P7(l2)
   do l3=(l2+1),(j-1)
    cst=cst*P1(l3)
   enddo
   trans(nclust,j)=trans(nclust,j)+cst*P2(j)
  enddo
  !#### recombination in j and coalesce in j
  trans(nclust,j)=trans(nclust,j)+P5(j)
 enddo
 !####### transition to non-hbd
 do i=1,nclust
  trans(i,nclust)=1.d0-sum(trans(i,1:(nclust-1)))
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

end subroutine


