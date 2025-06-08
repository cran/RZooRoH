subroutine zoosumlayerViterbi2(nclust,nchr,npos,pemission,chr_limits,zrates,zmix,posi,states)
implicit none
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
integer ::i,j,l,k,nclust,nchr,npos,chr,fpos,lpos,pos(1),phi(nclust,npos)
integer ::isF(nclust),chr_limits(nchr,2),posi(npos),states(npos)
integer ::l2,l3,nG,g1,g2,g3,recfun,d,d_i,lev1,lev2,fval,lval,lstep,nd
real(dp) ::pemission(npos,2),as(nclust),Fs(nclust),delta(nclust,npos)
real(dp) ::pinit(nclust),zmix(nclust-1),zrates(nclust-1)
real(dp) ::pmax,val(nclust),cst,prec,d2
real(dp), parameter ::Morgan=100000000.d0
real(dp),allocatable ::ptok(:),pnhbd(:),pnorec(:),ctok(:),MTRANS(:,:,:)
real(dp),allocatable ::ctok2(:),precing(:),cumupk(:),upcoal(:),upnocoal(:)
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


Fs(1:(nclust-1))=zmix;Fs(nclust)=1-zmix(nclust-1)
as(1:(nclust-1))=zrates;as(nclust)=zrates(nclust-1)
!delta=-1000000000.0;phi=0;states=0
phi=0;states=0
isF=1;isF(nclust)=0

nG=nint(as(nclust-1))
allocate(ptok(nG),pnhbd(nG),ctok(nG),ctok2(nG),pnorec(0:nG))
allocate(precing(nG),cumupk(nG),upcoal(nG),upnocoal(nG))
allocate(P1(nclust),P2(nclust),P3(nclust,nd),P4(nclust,nd),P5(nclust,nd))
allocate(P6(nclust,nd),P7(nclust,nd),MTRANS(nd,nclust,nclust))

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

!############ COMPUTE VALUES OF P3-P7 in a table ####

P3=0.d0 ! no recombination when in state i before
P4=0.d0 ! recombination and coalescence in state i, when in state i before
P5=0.d0 ! recombination and coalescence in state i, when in state > i before
P6=0.d0 ! recombination and floating in state i, when in state i before
P7=0.d0 ! recombination and floating in state i, when in state > i before [recombination in state i when i state i > before - rec & coalescen in i when i > before = P5]

MTRANS=0.d0

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
!  g2=nint(as(nclust-1));P4(nclust,d)=pnorec(g2-1)*prec*(1-Fs(nclust-1)) ## double counting?
  P5(nclust,d)=0.d0;P6(nclust,d)=0.d0

  do l=1,nclust-1
   g1=1;if(l>1)g1=nint(as(l-1))+1;g2=nint(as(l))
   cst=pnorec(g1-1)-pnorec(g2) !#### recombination between g1 and g2
   P7(l,d)=cst-P5(l,d)
  enddo
  P7(nclust,d)=0.d0

!####### compute transition matrices #######

  do i=1,nclust-1 ! non-HBD is a "special" case (for transition we can cut prob to non-HBD as 1 - sum(to HBD classes)
   do j=1,nclust-1 
     l=min(i,j) !### recombination must be before layer i and j
     do l2=1,l-1 !#### loop over all the layers before i and j
       cst=P7(l2,d) !### no recombination till the layer x rec in the layer + floating
       do l3=(l2+1),(j-1)
        cst=cst*P1(l3) !#### keep floating in all layers till layer j-1
      enddo
      MTRANS(d,i,j)=MTRANS(d,i,j)+cst*P2(j) !##### multiply by coalescence in layer j
     enddo
     !##### if i < j, possibility to recombine in i when previously in i
     if(i < j)then
       cst=P6(i,d)
       do l3=(i+1),(j-1)
        cst=cst*P1(l3) !#### keep floating in all layers till layer j-1
       enddo
       MTRANS(d,i,j)=MTRANS(d,i,j)+cst*P2(j) !##### multiply by coalescence in layer j
     endif
     !##### if i > j, can recombine in layers j and coalesce in j when in higher layer
     if(i > j)MTRANS(d,i,j)=MTRANS(d,i,j)+P5(j,d)
     !##### if i == j, recombine in layer j, coalesce in j when previously in j
     if(i == j)MTRANS(d,i,j)=MTRANS(d,i,j)+P4(j,d)
     !##### if i == j, non-recombination till j when previously in j
     if(i == j)MTRANS(d,i,j)=MTRANS(d,i,j)+P3(j,d)
   enddo
 enddo
 
 !####### transition from non-hbd
 do j=1,nclust-1
  !#### recombination in a layer x keep coalesce x coalesce in j
  do l2=1,j-1
   cst=P7(l2,d)
   do l3=(l2+1),(j-1)
    cst=cst*P1(l3)
   enddo
   MTRANS(d,nclust,j)=MTRANS(d,nclust,j)+cst*P2(j)
  enddo
  !#### recombination in j and coalesce in j
  MTRANS(d,nclust,j)=MTRANS(d,nclust,j)+P5(j,d)
 enddo
 !####### transition to non-hbd
 do i=1,nclust
  MTRANS(d,i,nclust)=1.d0-sum(MTRANS(d,i,1:(nclust-1)))     
 enddo    

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
 
  d2=1.d0*(posi(k)-posi(k-1))
  d_i=dist_index(d2)

 do i=1,nclust
  do j=1,nclust
    val(j)=delta(j,k-1)+log(MTRANS(d_i,j,i))
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


