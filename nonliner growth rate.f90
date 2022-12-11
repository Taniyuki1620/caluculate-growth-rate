program nonliner
implicit none
integer::i,j
integer,parameter ::nu=10
real(8)::Q,oph,vg,gam,gam2,Utpa,xi,omega,ow,del,Upeh,pi,c,vr,gamman
real(8)::Nh,e,m,ip
real(8)::Wh,rho,beta,betah,N0,Vpe0
real(8)::a,h,A0,Utpe
real(8)::vp,Bw,oe,B0,gn,ope
real(8)::a1,a2,a3
real(8)::vpe,omega1,oph1,Upe0,Utpa1,vr1,s2,oth,oth1,oth2,gnoth
real(8)::gparam(4)
real(8)::tau,vp1,vg1,oop,oop1,oop2,gnoop,Ramda,v,vpe01
real(8)::dw,w,s0,s1,s3,wc,ar,lv,re,fce,vh
real(8)::iparam(5),oparam(4),g,gop,gth
real(8)::u(nu),p(nu),nume(nu),deno(nu),sumnume,sumdeno,Upe01,Upe02




Utpa=0.25
beta=0.3
Utpe=0.3

pi=acos(-1d0)
fce=7.6d3
re=6380d3
lv=4.5
c=3d8
wc=2*pi*fce
ar=4.5/(lv*re)**2
a=ar*c**2/wc**2
vh=2d-3
dw=0.1

iparam(1)=Utpa
iparam(2)=Utpe
iparam(3)=beta
iparam(4)=vh
iparam(5)=a

sumnume=0
sumdeno=0





do j=2,8
if (j==2.or.j==4.or.j==8) then
    
    ope=j
   
do i=1,10
    w=dw*i

    !print*,w
    call nl(w,ope,iparam,oparam)
    

    write(13,*)i,oparam(3)
    write(14,*)i,oparam(4)
end do
else

ope=0

end if

end do

do i=1,nu

u(i)=i/10
p(i)=1./(1-beta)*(exp(-1d0*u(i)**2/(2*Utpe**2))-exp(-1d0*u(i)**2/(2*beta*Utpe**2)))
deno(i)=u(i)*p(i)*2*pi*u(i)
nume(i)=p(i)*2*pi*u(i)

end do

do i=1,nu-1
sumdeno=(u(i+1)-u(i))*(deno(i+1)+deno(i))/2
sumnume=(u(i+1)-u(i))*(nume(i+1)+nume(i))/2
end do

Upe01=sumdeno/sumnume
Upe02=(pi/2)**(0.5)*(1-beta**(1.5))/(1-beta)*Utpe
print*,u(10),p(10),deno(10),nume(10)
print*,Upe01,Upe02

end program




subroutine nl(w,ope,iparam,oparam)
implicit none
integer::i,j
real(8)::Q,oph,vg,gam,gam2,Utpa,xi,omega,ow,del,Upeh,pi,c,vr,gamman
real(8)::Nh,e,m,ip
real(8)::Wh,rho,beta,betah,N0,Vpe0
real(8)::a,h,A0,Utpe
real(8)::vp,Bw,oe,B0,gn,ope
real(8)::a1,a2,a3
real(8)::vpe,omega1,oph1,Upe0,Utpa1,vr1,s2,oth,oth1,oth2,gnoth
real(8)::gparam(4)
real(8)::tau,vp1,vg1,oop,oop1,oop2,gnoop,Ramda,v,vpe01
real(8)::dw,w,s0,s1,s3,wc,ar,lv,re,fce,vh
real(8)::iparam(5),oparam(4),g,gop,gth,gam3

Q=0.5
tau=0.5
c=3d8
dw=0.01
pi=acos(-1d0)



Utpa=iparam(1)
Utpe=iparam(2)
beta=iparam(3)
vh=iparam(4)
a=iparam(5)


oph=ope*vh**(0.5)
vpe=Utpe*(pi/2)**(0.5)*(1-beta**(1.5))/(1-beta)
xi=(w*(1d0-w)/ope**2d0)**(0.5)
del=(1d0/(1d0+xi**2d0))**(0.5)
vp=del*xi
vr=(w**2d0-(w**4d0+(w**2d0+vp**2d0)*(1d0-w**2d0-vpe**2d0))**(0.5))*vp/(vp**2d0+w**2d0)

call dichotomy(w,ope,iparam,gam3)

gam=1d0/(1d0-(vr**2d0+vpe**2d0))**(0.5)
gam2=1/(1-(vp**2+vpe**2))**(0.5)!
vg=(xi/del)/(xi**2+1d0/(2*(1d0-w)))
s0=del*vpe/xi 
s1=gam*(1d0-vr/vg)**2
Upe0=vpe*gam
s3=gam*w*(vpe**2)-(2+del**2*(1d0-gam*w)/(1d0-w))*vr*vp
s2=s3/(2*xi*del)
oop1=0.81*pi**(-2.5)*Q*s1*vg/(s0*w*tau*Utpa)
oop2=(oph*Upe0*del/gam)**2*exp(-1d0*0.5*(gam*vr/Utpa)**2)
oop=oop1*oop2
oth1=100*pi**3*gam**4*xi/(w*oph**4*Upe0**5*del**5)
oth2=(a*s2*Utpa/Q)**2*exp((gam*vr/Utpa)**2)
oth=oth1*oth2
g=Q*oph**2*vg/(2*gam*Utpa)*(xi/w)**(0.5)*(del*Upe0/pi)**(1.5)*exp(-0.5*(gam*vr/Utpa)**2)
gop=g*oop**(-1d0*0.5)
gth=g*oth**(-1d0*0.5)

!print*,gam3

oparam(1)=oop
oparam(2)=oth
oparam(3)=gop
oparam(4)=gth


end subroutine

subroutine func_gamma(w,ope,gam,iparam,f)

real(8)::gam,iparam(5),Utpa,vr,vp,c,chi,xi,omega,oe,ope,vpa,vpe,f,del
real(8)::pi,row,b,Utpe,A0,Wh,bh,Upah,beta,w
c=3d8
pi=acos(-1d0)

Utpa=iparam(1)
Utpe=iparam(2)
beta=iparam(3)



xi=((w*(1-w))/ope**2)**(0.5)!O2021(eq4)
del=(1./(1+xi**2))**(0.5)!O2021(eq6)
vp=del*xi
vr=(1-1/gam/w)*vp!O2021(eq22)
vpa=vr
A0=Utpe**2/Utpa**2-1
Wh=1d0
Upah=(pi/2)**(0.5)*Utpe*(1d0-beta**(1.5))/(1d0-beta)
vpe=Upah/gam
f=gam-1./(1-(vpa**2+vpe**2))**(0.5)

!
!print*,gam
return 

end subroutine



subroutine dichotomy(w,ope,iparam,gam3)
implicit none
integer::loopmax=100,i
real(8)::Utpa,gam,c,oe,ope,omega,iparam(5)
real(8)::f,a,b,eps,res,fa,fc,gam3,w


a=1
b=2
c=0
eps=1d-6

!print*,gparam
do i=1,loopmax
!print*,a,b,c
    c=0.5*(a+b)
    gam=a
    !print*,"a=",a
    call func_gamma(w,ope,gam,iparam,f)
    
    fa=f
    gam=c
    !print*,"c=",c
    call func_gamma(w,ope,gam,iparam,f)
    
    fc=f
    !print*,fa,fc
    if(fa*fc<0)then
        b=c
        else
        a=c
    end if
    res=abs(a-b)
    !print*,a,b,c
    if(res<eps)exit
end do
!print*,gam
gam3=gam

!print*,a,b,c

end subroutine


subroutine f_K_a
implicit none
integer::i,j
integer,parameter ::nK=5,na=90
real(8)::f(nK,na),K(nK),a(na)

do j=1,na
    a(j)=10*na
    do i=1,nK
        K(i)=10**i
    end do
end do

end subroutine
