
subroutine get_BiMaxwellian(ppa0,ppe0,apa,ape,f0)
implicit none
real(8) :: ppa0,ppe0,f0
real(8) :: am,apa,ape,pi
pi=acos(-1d0)

f0=1d0/(pi**(3d0/2d0)*ape**2*apa)*exp(-1d0*(ppa0**2/apa**2+ppe0**2/ape**2))

return
end subroutine

subroutine get_kappa(ppa0,ppe0,spa,spe,ka,f0)
implicit none
real(8) :: ppa0,ppe0,f0
real(8) :: ka,pi,spa,spe
pi=acos(-1d0)

f0=GAMMA(ka+1)/(pi**(3d0/2d0)*spe**2*spa*ka**(3d0/2d0)*GAMMA(ka-0.5))*(1+(ppa0**2/ka/spa**2)+(ppe0**2/ka/spe**2))**(-ka-1)
!print*,f0
return
end subroutine

subroutine direct_kappappe(ppa1,ppe1,pi,spa,spe,ka,pped)
implicit none
integer::i,j,nx,ny
real(8)::ppe1,ppa1
real(8)::pped1,pped
real(8)::spa,spe,pi,ka
pi=acos(-1d0)

pped1=GAMMA(ka+1)/(pi**(3d0/2d0)*spe**2*spa*ka**(3d0/2d0)*GAMMA(ka-0.5))*(1+(ppa1**2/ka/spa**2)+(ppe1**2/ka/spe**2))**(-ka-2)
pped=pped1*(-ka-1)*2*ppe1/ka/spe**2

end subroutine

subroutine direct_kappappa(ppa1,ppe1,pi,spa,spe,ka,ppad)
implicit none
integer::i,j,nx,ny
real(8)::ppe1,ppa1
real(8)::ppad1,ppad
real(8)::spa,spe,pi,ka
pi=acos(-1d0)

ppad1=GAMMA(ka+1)/(pi**(3d0/2d0)*spe**2*spa*ka**(3d0/2d0)*GAMMA(ka-0.5))*(1+(ppa1**2/ka/spa**2)+(ppe1**2/ka/spe**2))**(-ka-2)
ppad=ppad1*(-ka-1)*2*ppa1/ka/spa**2

end subroutine

subroutine directppe(ppe1,ppa1,pe,apa,ape,pped)
implicit none
integer::i,j,nx,ny
real(8)::ppe1,ppa1
real(8)::pped
real(8)::apa,ape,pe

pped=1d0/(pe**(3d0/2d0)*ape**2*apa)*exp(-1d0*(ppa1**2/apa**2+ppe1**2/ape**2))*(-2d0*ppe1/ape**2)


end subroutine

subroutine directppa(ppe1,ppa1,pe,apa,ape,ppad)
implicit none
integer::i,j,nx,ny
real(8)::ppe1,ppa1
real(8)::ppad
real(8)::apa,ape,pe


ppad=1d0/(pe**(3d0/2d0)*ape**2*apa)*exp(-1d0*(ppa1**2/apa**2+ppe1**2/ape**2))*(-2d0*ppa1/apa**2)

end subroutine

subroutine sympthon(fs,xs,nx,sum)
implicit none
integer::nx,i
real(8)::fs(nx),xs(nx)
real(8)::sum1,sum2,sum,h
sum1=0
sum2=0
sum=0
h=xs(2)-xs(1)

do i=1,nx-1
   if (mod(i,2)==1)then
      sum1=sum1+fs(i)
   else
      sum2=sum2+fs(i)
   end if
end do

sum=h/3*(4*sum1+2*sum2+fs(nx))

end subroutine
   


program fm1
implicit none
integer::i,j,g1,g2,l,o,ierr,p,n
integer,parameter ::nx=10,ny=10
real(8)::ppa(nx),ppamax,ppamin
real(8)::ppe(ny),ppemax,ppemin
real(8)::f(nx,ny)
real(8)::a,b,f9,pr(ny),dx,dr(ny),gr(ny),k,or,oe,ope,x,e,B0,ac,oi
real(8)::pmax,m,eaveeV,eave,v
real(8) :: ppa0,ppe0,f0
real(8) :: am,apa,ape,pi,tpa,ppeave,ppaave,ka,spe,spa
real(8) ::bif,smf,df(nx),sumdeno,sumnume,deff,deno(nx),nume2(nx),nume1(nx),nume(nx),arel
real(8) ::def(nx),sample(nx),pe,c
real(8) ::ppedx(ny),ppadx(nx),ppe1,pped,ppa1,ppad,sumdeno1,sumnume1,numed(nx),denod(nx),arel1,vh,H,H1,oi1
real(8)::fs(nx),xs(nx),sum,sumnume2,sumdeno2,H2,arel2,oi2
real(8)::fp(nx,ny)


m=9.1d-31!9.1d-31
c=3d8
eaveeV=100d3*2.5
eave=eaveeV*1.6d-19
v=0.99*c
pmax=v/sqrt(1-v**2/c**2)
ppamax=pmax
ppamin=-pmax
ppemax=pmax
ppemin=sqrt(2d0/m*100*1.6d-19)!momentum corresponding to 100 eV
dx=ppemin/5d0
sumdeno=0
sumnume=0
sumdeno1=0
sumnume1=0
tpa=25d3*1.6d-19
am=2
pe=acos(-1d0)
apa=(2*tpa/m)**(0.5)
ape=(am+1)*apa
ppaave=tpa/m
arel1=0
vh=0.01
H=0
H1=0
!e=1.6d-19
!B0=1.2d-6
!oe=e*B0/m
!or=x*oe
!ope=3d0*oe
ac=0
B0=1.2d-6
oi=0
oi1=0
ka=2
print*,v/c
!stop
do i=1,nx
   ppa(i)=(ppamax-ppamin)/dble(nx-1)*dble(i-1)+ppamin
   !print*,"ppa=",ppa(i)
end do
do j=1,ny
   ppe(j)=(ppemax-ppemin)/dble(ny-1)*dble(j-1)+ppemin
   !print*,"ppe=",ppe(j)
end do


!print*,ppaave

do i=1,nx
   ppa0=ppa(i)
   do j=1,ny
      ppe0=ppe(j)
      apa=(2*tpa/m)**(0.5)
      ape=((am+1)*apa**2d0)**(0.5)
      spa=((2*ka-3)/ka)**(0.5)*ppaave**(0.5)
      spe=((am+1)*spa**2)**(0.5)
      
      
      !call get_BiMaxwellian(ppa0,ppe0,apa,ape,f0)
      call get_kappa(ppa0,ppe0,spa,spe,ka,f0)
      f(i,j)=f0
   !print*,i,j,f(i,j)
   end do   
end do

!print*,spa,spe
do n=2,9
   x=0.1*n
   print*,n
call prr(ny,ppe,pr,dr,gr,k,or,oe,ope,x)


ac=or/(abs(oe)-or)


do o=1,ny!（解析値）
   ppe1=ppe(o)
   
   ppa1=pr(o)
   !print*,pr(o)
   !call directppe(ppe1,ppa1,pe,apa,ape,pped)
   call direct_kappappe(ppa1,ppe1,pe,spa,spe,ka,pped)
   ppedx(o)=pped
   !print*,x,o,pped
   
   !call directppa(ppe1,ppa1,pe,apa,ape,ppad)
   call direct_kappappa(ppa1,ppe1,pe,spa,spe,ka,ppad)
   ppadx(o)=ppad
   !print*,"ppadx",ppadx(o)
   write(10,*)o,ppe(o)*ppadx(o),pr(o)*ppedx(o)
   
   denod(o)=ppe(o)**2*ppedx(o)/dr(o)
   numed(o)=ppe(o)**2*(ppe(o)*ppadx(o)-pr(o)*ppedx(o))/dr(o)/gr(o)
  
end do
!print*,(ppe(5)*ppadx(5)-pr(5)*ppedx(5)),numed(5)
close(10)
!print*,ppedx(1)
sumdeno1=0
sumnume1=0
do l=2,nx-2
   sumdeno1=sumdeno1+(denod(l+1)+denod(l))*(ppe(l+1)-ppe(l))/2d0
   !print*,deno(l+1),deno(l),ppe(l+1),ppe(l)
   sumnume1=sumnume1+(numed(l+1)+numed(l))*(ppe(l+1)-ppe(l))/2d0

end do

do o=1,nx
fs(o)=denod(o)
xs(o)=ppe(o)
end do


call sympthon(fs,xs,nx,sum)
sumdeno2=sum

do o=1,nx
fs(o)=numed(o)
xs(o)=ppe(o)
end do

call sympthon(fs,xs,nx,sum)
sumnume2=sum

!print*,sumnume1,sumnume2
!write(11,*)sumnume1
arel1=(k/(or-abs(oe)))*sumnume1/sumdeno1
H1=pe*vh*((or-abs(oe))/k)*sumdeno1
oi1=(pe*ope**2*H1*(arel1-ac))/(2*or+ope**2*abs(oe)/(or-abs(oe))**2)

arel2=(k/(or-abs(oe)))*sumnume2/sumdeno2
H2=pe*vh*((or-abs(oe))/k)*sumdeno2
oi2=(pe*ope**2*H2*(arel2-ac))/(2*or+ope**2*abs(oe)/(or-abs(oe))**2)
!write(11,*)x,arel1
!write(12,*)x,H1
!write(13,*)x,oi1/oe
!print*,"arel=",arel1,sumnume1,(k/(or-abs(oe)))/sumdeno1

!print*,pe,ope**2,H1*(arel1-ac)
!write(13,*)x,H1
!end do
!write(13,*)x,H1!oi1/abs(oe)

!do n=1,99
 !  x=0.01*n

call prr(ny,ppe,pr,dr,gr,k,or,oe,ope,x)
!print*,n,pr
ac=or/(abs(oe)-or)
do o=2,nx-1!(線形補間)
   a=pr(o)
   b=ppe(o)+dx
   call fm2(nx,ny,ppa,ppe,f,a,b,f9,ierr)
   bif=f9
   !ppa0=a
   !ppe0=b
   !call get_BiMaxwellian(ppa0,ppe0,apa,ape,f0)
   !print*,f9

   
   b=ppe(o)-dx
   call fm2(nx,ny,ppa,ppe,f,a,b,f9,ierr)
   smf=f9
   !ppa0=a
   !ppe0=b
   !call get_BiMaxwellian(ppa0,ppe0,apa,ape,f0)
   !print*,f9
  

   call defferential(bif,smf,dx,deff)
   !sample(o)=1d0/(pe**(3d0/2d0)*ape**2*apa)*exp(-1d0*(pr(o)**2/apa**2+ppe(o)**2/ape**2))*(-2d0*ppe(o)/ape**2)
   !print*,n,(deff-ppedx(o))/ppedx(o)
   !print*,x,o,(deff-sample(o))/sample(o),bif,smf,sample(o)
   !write(13,*)x,o,(deff-sample(o))/sample(o),bif,smf,sample(o)
   deno(o)=ppe(o)**2*deff/dr(o)
   nume2(o)=pr(o)*deff
   !print*,"deno(o)=",deno(o)
  
   a=pr(o)+dx
   b=ppe(o)
   call fm2(nx,ny,ppa,ppe,f,a,b,f9,ierr)
   bif=f9
   ppa0=a
   ppe0=b
   call get_kappa(ppa0,ppe0,spa,spe,ka,f0)
   !print*,n,(f9-f0)/f0


   a=pr(o)-dx
   call fm2(nx,ny,ppa,ppe,f,a,b,f9,ierr)
   smf=f9
   ppa0=a
   ppe0=b
   call get_kappa(ppa0,ppe0,spa,spe,ka,f0)
   !print*,n,(f9-f0)/f0
   !print*,(smf-f0)/f0
   call defferential(bif,smf,dx,deff)
   !write(11,*)n,deff,ppadx(o),bif,smf
   nume1(o)=ppe(o)*deff
   !print*,deff
   nume(o)=ppe(o)**2*(nume1(o)-nume2(o))/gr(o)/dr(o)
   !print*,nume1(o),nume2(o),ppe(o)
   
end do


sumdeno=0
sumnume=0
do l=2,nx-2
   sumdeno=sumdeno+(deno(l+1)+deno(l))*(ppe(l+1)-ppe(l))/2
   sumnume=sumnume+(nume(l+1)+nume(l))*(ppe(l+1)-ppe(l))/2

end do
arel=(k/(or-abs(oe)))*sumnume/sumdeno
H=pe*vh*((or-abs(oe))/k)*sumdeno
oi=(pe*ope**2*H*(arel-ac))/(2*or+ope**2*abs(oe)/(or-abs(oe))**2)
oi1=(pe*ope**2*H1*(arel1-ac))/(2*or+ope**2*abs(oe)/(or-abs(oe))**2)
write(12,*)x,H1,H
!write(12,*)x,(oi/oe),(oi1/oe)

end do

end program fm1

subroutine defferential(bif,smf,dx,deff)
implicit none
real(8)::bif,smf,deff,dx
deff=(bif-smf)/(2*dx)

end subroutine defferential


subroutine fm2(nx,ny,ppa,ppe,f,a,b,f9,ierr)!a,b与えて線形補間
implicit none
integer::i,j,nx,ny,g1,g2,k,l,ierr

real(8)::a,f1,f2,ppa(nx)
real(8)::b,f3,f4,ppe(ny)
real(8)::f(nx,ny),pr(ny)
real(8)::f5,f6,f7,f8,f9
ierr=0
!a=1d5
!b=1d5
g1=1000
g2=1000
if (a<ppa(1).or.a>ppa(nx).or.b<ppe(1).or.b>ppe(ny)) then
   ierr=1
   !print*,a,ppa(1),ppa(nx),b,ppe(1),ppe(ny)
   !print*,"out"
   

end if

do i=1,nx-1
   if (ppa(i)<a)then
      f1=ppa(i)
      f2=ppa(i+1)
      g1=i
      !print*,i,ppa(i),a
end if
end do
!print*,a,ppa(1),ppa(10)

if(g1==1000)then
!print*,"g1 is out"
!stop
end if


do j=1,ny-1
   if (ppe(j)<b)then
      f3=ppe(j)
      f4=ppe(j+1)
      g2=j
      !print*,b
end if
end do

if(g2==1000)then
print*,"g2 is out",b,ppe(1),ppe(ny)
stop
end if
!print*,g1,g2

f5=f(g1,g2)*(f2-a)/(f2-f1)*(f4-b)/(f4-f3)
f6=f(g1,g2+1)*(f2-a)/(f2-f1)*(b-f3)/(f4-f3)
f7=f(g1+1,g2+1)*(a-f1)/(f2-f1)*(b-f3)/(f4-f3)
f8=f(g1+1,g2)*(a-f1)/(f2-f1)*(f4-b)/(f4-f3)
f9=f5+f6+f7+f8
print*,f5,f6,f7,f8
!print*,f(g1,g2)!,(f2-a),(f2-f1),(f4-b),(f4-f3)

do i=1,nx-1
   !print*,i,g1,ppa(i),a,f(i,g2)
end do


end subroutine fm2

subroutine prr (ny,ppe,pr,dr,gr,k,or,oe,ope,x)
implicit none 
integer::  i,ny,j
integer,parameter ::nor=1000
real(8)::  gr(ny),c,k,or,oe,ope,pr(ny),dr(ny)
real(8)::oor(nor),kk(nor)
real(8)::  e,B0,me,n,x
real(8)::  ppe1(ny),ppe(ny)
n=1d2
c=3d8
e=1.6d-19
B0=1.2d-6
me=9.1d-31
oe=e*B0/me
or=x*oe
ope=3d0*oe
k=((or**2d0/c**2d0)-(or*ope**2d0/c**2d0/(or-abs(oe))))**(0.5d0)
!print*,"a",k
do i=1,ny
   
gr(i)=(-1d0+(c*k/or)*(((c*k/or)**2d0-1d0)*(1d0+ppe(i)**2d0/c**2d0)*(or/abs(oe))**2d0+1d0)**(0.5d0))/(((c*k/or)**2d0-1d0)*(or/oe))
pr(i)=(gr(i)*or-abs(oe))/k
dr(i)=1-or*pr(i)/c**2/k/gr(i)
!print*,x,gr(i),dr(i)
!write(11,*)x,gr(i),dr(i)
end do
!print*,x,gr(5),dr(5),pr(5)
!print*,k
end subroutine prr







