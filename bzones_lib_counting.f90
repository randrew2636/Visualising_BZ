subroutine gen_kp(m1,m2,m3,n1,n2,n3,bz,GS,G,l)

implicit none

integer, intent(in):: n1,n2,n3
Integer, intent(in):: m1,m2,m3,bz
double precision, intent(in):: GS(3,3)
integer, intent(in):: l
double precision, intent(in):: G(l,3)

double precision  :: kp(3)
double precision :: kr, kt

integer :: i,j,k,count,s
logical :: keep

open(1, file="BZ_ZONE.dat", status="new", action="write")


!$omp parallel private(kp,kr,count,kt,keep) 
!$omp do
do i=0,n1-1
	do j=0,n2-1
		do k=0,n3-1
			count=0
			keep = .true.
			kp(1)=m1*(-1.0+2.0*i/(n1-1))*Gs(1,1)+m2*(-1.0+2.0*j/(n2-1))*Gs(2,1)+m3*(-1.0+2.0*k/(n3-1))*Gs(3,1)
			kp(2)=m1*(-1.0+2.0*i/(n1-1))*Gs(1,2)+m2*(-1.0+2.0*j/(n2-1))*Gs(2,2)+m3*(-1.0+2.0*k/(n3-1))*Gs(3,2)
			kp(3)=m1*(-1.0+2.0*i/(n1-1))*Gs(1,3)+m2*(-1.0+2.0*j/(n2-1))*Gs(2,3)+m3*(-1.0+2.0*k/(n3-1))*Gs(3,3)
			kt=(kp(1)**2+kp(2)**2+kp(3)**2)
			do s=1,l
				kr=(G(s,1)-kp(1))**2+(G(s,2)-kp(2))**2+(G(s,3)-kp(3))**2
				if(kt.gt.kr) then
					count=count+1
				end if
				if(count .ge. bz) then
					keep = .false.
					exit 
				endif				
			enddo
			if (((count+1) .eq. bz) .AND. keep) write(1,666) kp(1), kp(2),kp(3), SQRT(kt), count+1
		enddo
	enddo
enddo
!$omp end do
!$omp end parallel
 close(1)
666 format(f20.3,1x,f20.3,1x,f20.3,1x,f20.3,1x,i3)

end subroutine gen_kp