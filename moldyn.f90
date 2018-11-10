real function ranf(seed)
integer::a=16807,c=2147483647,q=127773,r=2386,h,t,l,seed
real::cd
cd=c
h=seed/q
l=mod(seed,q)
t=a*l-r*h
if(t>0)then
seed=t
else
seed=c+t
endif
ranf=seed/cd
end function

!*****************************************************//Force// ********************************************************
subroutine force(f,r,box)
real(kind=8)::f(3,300),d(3),r(3,300),rij,box(3),sig=3.2874,rc,eps=500
integer::n=300
rc=2.0**(1/6)*sig

do i=1,n,1 
	do k=1,3,1
	f(k,i)=0.0
	enddo
enddo

do i=1,n,1
	do j=i+1,n,1
	rij=0.0
		do k=1,3,1
		d(k)=r(k,i)-r(k,j)
		d(k)=d(k)-box(k)*nint(d(k)/box(k))
		rij=rij+d(k)**2
		enddo
	rij=sqrt(rij)
		if(rij.lt.rc)then
		do k=1,3,1
		ff=-48*eps/(rij)**2*((sig/rij)**12-0.5*(sig/rij)**6)
		f(k,i)=f(k,i)+ff*d(k)
		f(k,j)=f(k,j)-ff*d(k)
		enddo
		endif
	enddo
enddo

end subroutine

!*****************************Calculating bond force using two subroutines****************************

!first subroutine
subroutine fbcal(fi,fj,ri,rj,box)
real(kind=8)::fi(3),fj(3),ri(3),rj(3),kf,rij,ff,rc,r0,d(3),box(3)
kf=1.18E+5;r0=2.34;rc=2.0**(1/6)*sig
rij=0.0
do k=1,3,1
	d(k)=ri(k)-rj(k)
	d(k)=d(k)-box(k)*nint(d(k)/box(k))
	rij=rij+d(k)**2
enddo
rij=sqrt(rij)
	if(rij.lt.rc)then
	ff=-kf*(rij-r0)/rij
	do k=1,3,1
	fi(k)=fi(k)+ff*d(k)
	fj(k)=fj(k)-ff*d(k)
	enddo
	endif
end subroutine
	
!final subroutine that calls fbcal
subroutine fbond(f,r,box)

real(kind=8)::f(3,300),r(3,300),d(3),rij,box(3),sig=3.2874,rc,eps=500
integer::n=300
rc=2.0**(1/6)*sig

do i=1,n,30
	do j=i,i+29,1
	if(mod(j,30).eq.1)then
	call fbcal(f(:,j),f(:,j+1),r(:,j),r(:,j+1),box)
	elseif(mod(j,30).eq.0)then
	call fbcal(f(:,j-1),f(:,j),r(:,j-1),r(:,j),box)
	else
	call fbcal(f(:,j-1),f(:,j),r(:,j-1),r(:,j),box)
	call fbcal(f(:,j),f(:,j+1),r(:,j),r(:,j+1),box)
	endif
	enddo
enddo
end subroutine

!**********************************************Radial distribution function*****************************************
subroutine raddf(r,box,gr)
real(kind=8)::r(3,300),box(3),dr,dv,rho,nid,d(3),rij,m=1.0
real(kind=8),dimension(0:99)::gr
integer::nb=100,n=300,bnum
dr=15.0/(nb-1)
rho=n/(box(1)*box(2)*box(3))

do i=0,nb-1,1
	gr(i)=0.0
enddo

do i=1,n-1,1
	do j=i+1,n,1
	rij=0.0
	do k=1,3,1
	d(k)=r(k,i)-r(k,j)
	d(k)=d(k)-box(k)*nint(d(k)/box(k))
	rij=rij+d(k)**2
	enddo
	rij=sqrt(rij)
	if(rij.lt.(15.0))then
	bnum=nint(rij/dr)
	gr(bnum)=gr(bnum)+2
	endif
	enddo
enddo

do i=0,nb-1,1
	dv=4/3.0*22/7.0*((i+1)**3-i**3)*dr**3
	nid=rho*dv
	gr(i)=gr(i)/(nid*n)
enddo

end subroutine

!***************************************************Bond angle*****************************************************
!Calculating bond angle for three atoms
real function bangle(ri,rj,rk,box)
real(kind=8)::ri(3),rj(3),rk(3),a(3),b(3),box(3),x,pi=22.0/7.0
real(kind=8)::a_mag=0.0,b_mag=0.0,adotb=0.0
do l=1,3,1
	a(l)=ri(l)-rj(l)
	b(l)=rk(l)-rj(l)
	a(l)=a(l)-box(l)*nint(a(l)/box(l))
	b(l)=b(l)-box(l)*nint(b(l)/box(l))
	a_mag=a_mag+a(l)**2
	b_mag=b_mag+b(l)**2
	adotb=adotb+a(l)*b(l)
enddo
a_mag=sqrt(a_mag)
b_mag=sqrt(b_mag)
x=acos(adotb/(a_mag*bmag))
x=abs(x*180.0/pi)
if(x.gt.180.0)then
bangle=360.0-x
else
bangle=x
endif
end function
!Using the previous function to generate the bond angle distribution
subroutine bacalc(r,box,Nth)
real(kind=8)::r(3,300),box(3),dth
real(kind=8),dimension(0:99)::Nth
integer::c=1,nbth=100,bnum,n=300
do i=0,nbth-1,1
Nth(i)=0.0
enddo
dth=180.0/(nbth-1)
do i=1,n,30
	do j=i,i+27,1
	bnum=nint(bangle(r(:,j),r(:,j+1),r(:,j+2),box)/dth)
	Nth(0:bnum)=Nth(0:bnum)+1
	enddo
enddo


end subroutine
!************************************************Dihedral angle calculation*****************************************

!Cross product
subroutine crossprod(a,b,c)
real(kind=8)::a(3),b(3),c(3)
c(1)=a(2)*b(3)-b(2)*a(3)
c(2)=b(1)*a(3)-a(1)*b(3)
c(3)=a(1)*b(2)-b(1)*a(2)
end subroutine


!Angle between two planes
real function dangle(r1,r2,r3,r4)
real(kind=8)::r1(3),r2(3),r3(3),r4(3),a(3),b(3),c(3),axb(3),bxc(3),axb_mag=0.0,bxc_mag=0.0,axbdotbxc=0.0,x


do i=1,3,1
a(i)=r1(i)-r2(i)
b(i)=r2(i)-r3(i)
c(i)=r3(i)-r4(i)
enddo


call crossprod(a,b,axb)
call crossprod(b,c,bxc)


do i=1,3,1
axbdotbxc=axbdotbxc+axb(i)*bxc(i)
axb_mag=axb_mag+axb(i)**2
bxc_mag=bxc_mag+bxc(i)**2
enddo


axb_mag=sqrt(axb_mag)
bxc_mag=sqrt(bxc_mag)


x=acos(axbdotbxc/(axb_mag*bxc_mag))
pi=22.0/7.0
x=abs(180*x/pi)
if(x.gt.180)then
dangle=360-x
else
dangle=x
endif


end function


!Dihedral angle dist fn
subroutine dacalc(r,box,Nd)
real(kind=8)::r(3,300),box(3),dd
real(kind=8),dimension(0:99)::Nd
integer::nbd=100,n=300,bnum
dd=180.0/(nbd-1)


do i=0,nbd-1,1
Nd(i)=0.0
enddo

do i=1,n,30
  do j=i,i+26,1
  bnum=nint(dangle(r(:,j),r(:,j+1),r(:,j+2),r(:,j+3))/dd)
 Nd(0:bnum)=Nd(0:bnum)+1
  enddo
enddo


end subroutine



!****************************************************//PBC//********************************************************
subroutine pbc(d,box)
real(kind=8)::d(3),box(3)
do k=1,3,1
if(d(k)>box(k))d(k)=d(k)-box(k)
if(d(k)<0.0)d(k)=d(k)+box(k)
enddo
end subroutine

!*************************************************//Main~~~~~~~~~~~~~~~~~~Program//***************************************************
!******************************************************~~~~~~~Starts~~~~~~~~**********************************************************
!*************************************************************************************************************************************
program md
real(kind=8)::r(3,300),d(3),rij,box(3),sig=3.2874,ac=30.0,v(3,300),sv(3),sv2,kb=1.0,m=1.0,fs,T,rp(3,300),f(3,300),temp(3)
real(kind=8)::dr,rho,dth,dd
real(kind=8),dimension(0:99)::Nth,Ntht,Nd,Ndt
real(kind=8),dimension(0:99)::gr,grt
integer::s(3),c=0,n=300,nb=100,nbth=100,nbd=100,cg=0,cba=0,cda=0
s(1)=1234; s(2)=4567;s(3)=7890
box(1)=ac;box(2)=ac;box(3)=ac
open(1,file="pairs.txt")
write(1,*)"Pair#             			Rij"
do i=1,n,1
11	do k=1,3,1
	r(k,i)=(box(k)-0.0)*ranf(s(k))+0.0
	s(k)=s(k)+100
	enddo

	do j=i-1,1,-1
	rij=0.0
		do k=1,3,1
		d(k)=(r(k,i)-r(k,j))
		d(k)=d(k)-box(k)*nint(d(k)/box(k))
		rij=rij+d(k)**2
		enddo
	rij=sqrt(rij)
		if(rij.le.sig)then
		goto 11
		else
		c=c+1
		write(1,*)c,rij
	write(*,*)c
		endif
	enddo
enddo
close(1)

open(2,file="pairs_check.txt")
write(1,*)"Pair#             			Rij"
c=0
do i=1,n,1
	do j=i+1,n,1
	rij=0.0
		do k=1,3,1
		d(k)=(r(k,i)-r(k,j))
		d(k)=d(k)-box(k)*nint(d(k)/box(k))
		rij=rij+d(k)**2
		enddo
	rij=sqrt(rij)
	c=c+1
	write(2,*)c,rij
	enddo
enddo
close(2)

!Define velocities
a=-0.5
b=0.5
sv2=0.0
open(3,file="velin.txt")
write(3,*)"#          vx(#)               vy(#)             vz(#)"
do i=1,n,1
	do k=1,3,1	
	v(k,i)=(b-a)*ranf(s(k))+a
	s(k)=s(k)+100
	sv2=sv2+v(k,i)**2
	enddo
write(3,*)i,v(1,i),v(2,i),v(3,i)
enddo
close(3)

!Scaling the velocities
T=200.0
fs=sqrt(3*n*kb*T/(m*sv2))
do k=1,3,1
sv(k)=sum(v(k,:))
sv2=0.0
open(4,file="scvelin.txt")
enddo
do i=1,n,1
	do k=1,3,1
	v(k,i)=(v(k,i)-sv(k)/n)*fs
	sv2=sv2+v(k,i)**2
	enddo
write(4,*)i,v(1,i),v(2,i),v(3,i)
enddo
close(4)

T=(m*sv2/(3*n*kb))
print*,"temperature is:",T

!*****************************************************************Time Loop************************************************************

!time values
ti=0.0
dt=1E-15
tm=100.0

!previous positions
do i=1,n,1
	do k=1,3,1
	rp(k,i)=r(k,i)-v(k,i)*dt
	enddo
call pbc(rp(:,i),box)
enddo

!initialise radial dist fn
do i=0,nb-1,1
	gr(i)=0.0
enddo

!initialise bond angle dist fn
do i=0,nbth-1,1
	Nth(i)=0.0
enddo

!initialise dihedral angle dist fn
do i=0,nbd-1,1
	Ndt(i)=0.0
enddo

!time looping
c=0
cg=0
cba=0
cda=0
open(5,file="tloop.txt")
write(5,*)"TLoop#                  Time                        Temperature"
do while(ti<tm)   
sv2=0.0
call force(f,r,box)
call fbond(f,r,box)
	do i=1,n,1
		do k=1,3,1
		temp(k)=r(k,i)
		r(k,i)=2*r(k,i)-rp(k,i)-f(k,i)*dt**2/m
		v(k,i)=r(k,i)-rp(k,i)
		v(k,i)=v(k,i)-box(k)*nint(v(k,i)/box(k))
		v(k,i)=v(k,i)/(2*dt)
	sv2=sv2+v(k,i)**2
		rp(k,i)=temp(k)
		enddo
	call pbc(r(:,i),box)
	enddo
if(c.gt.482)then

	cg=cg+1
	call raddf(r,box,grt)
	do i=0,nb-1,1
	gr(i)=gr(i)+grt(i)
	enddo

	cba=cba+1
	
    call bacalc(r,box,Ntht)
	do i=0,nbth-1,1
	Nth(i)=Nth(i)+Ntht(i)
	enddo
	
	cda=cda+1
	call dacalc(r,box,Ndt)
	do i=0,nbd-1
	Nd(i)=Nd(i)+Ndt(i)
	enddo
	
endif

T=(m*sv2/(3*n*kb))
c=c+1
write(5,*)c,ti,T
ti=ti+dt
write(*,*)c
if(c.eq.10000)exit
enddo

close(5)


!storing radial dist fn
open(6,file="raddf.txt")

do i=0,nb-1,1
	gr(i)=gr(i)/cg
enddo

dr=15.0/(nb-1)
rho=n*m/(box(1)*box(2)*box(3))
write(6,*)"Value of r                         g(r)"

do i=0,nb-1,1
	write(6,*)i*dr,4*22/7.0*(i*dr)**2*rho*(gr(i)-1.0)
enddo
close(6)



!storing bond angle dist fn
open(7,file='bangle.txt')

do i=0,nbth-1,1
	Nth(i)=Nth(i)/cba
enddo

write(7,*)" Value of theta                  N(theta)"
dth=180.0/(nbth-1)

do i=0,nbth-1,1
	write(7,*)i*dth,Nth(i)
enddo
close(7)

!Storing dihedral angle dist fn
open(8,file="dangle.txt")

do i=0,nbd-1,1
	Nd(i)=Nd(i)/cda
enddo

write(8,*)"Value of theta_d                N(theta_d)"

dd=180.0/(nbd-1)

do i=0,nbd-1,1
	write(8,*)i*dd,Nd(i)
enddo
close(8)




end program
