!=- fortran-libraries
!=- © Stanislav Shirokov, 2014-2020

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
	module math
		use global

      integer,parameter :: N_medians = 1d2

		integer :: N_grid_LT = 2d2, N_grid_LT_a = 2d2, N_grid_LT_b = 2d1, N_points

		real(8) :: RightLimit = 1d150, LeftLimit = -1d150, pi = 4*datan(1d0), &
         a_error_lower, b_error_lower, a_error_upper, b_error_upper, a_start = 1d-4, b_start = 1d-4, a_max = 1d1, b_max = 1d0, &
         XYdXdY(4,N_GRB), mx, my, a_step, b_step, X2=1d150, va, vb, weights, Sx, Sy, &
         va_error_lower,vb_error_lower , va_error_upper, vb_error_upper, medar(4,N_medians), &
         add_median_left_border, add_median_right_border, x_step,x_step_start,x_step_final, &
         add_median_x,add_median_y,add_median_dx,add_median_dy, log_medians_parameter = 0, &
         add_median2_x,add_median2_y,add_median2_dx,add_median2_dy

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
		contains

         real function distance(x,y,z)
            real(8) x,y,z
               distance = ( x**2 + y**2 + z**2 )**0.5d0
            end function

			integer function sign_real(x) ; real(8) x;
					sign_real=1 ; if (x<0) sign_real=-1
				end function

			real(8) function mydexp(x) 	; real(8) x;
					n=1 ; mydexp = ( 1d0 + x/1d5 )**1d5
				end function

			real(8) function factorial(n)	; integer n;
					factorial=1d0; do while (n>0); factorial=factorial*n;	n=n-1; enddo
				end function

			real(8) function mydsin(x) 	; real(8) x;
					n=1; mydsin=x; do; mydsin=mydsin+(-1d0)**n * x**(2*n+1) / factorial (2*n+1)
						n=n+1 ; if (x**(2*n+1) / factorial (2*n+1) < 1d-19) exit ; end do
				end function

			subroutine logtrend(XY,a,b,Lb,Rb)
				real(8) XY(2,N_GRB),G(2,N_GRB),a,b,sumx,sumx2,sumy,sumxy,Lb,Rb;
					G(:,:)=0d0; forall (i=1:N_GRB,j=1:2, XY(j,i)>0) G(j,i)=log10(XY(j,i))
					j=0;k=0;sumx=0d0;sumx2=0d0;sumy=0d0;sumxy=0d0; N=0
					do j=1,N_GRB; if (G(2,j).ne.0d0 .and. G(1,j).le.log10(Rb) .and. G(1,j).ge.log10(Lb)) then!		write(*,*) N,G(1,j),G(2,j),XY(1,j),XY(2,j)
						sumx=sumx+G(1,j) ; sumx2=sumx2+G(1,j)**2d0 ; sumy=sumy+G(2,j) ; sumxy=sumxy+G(1,j)*G(2,j) ; k=k+1 ; endif ; enddo
					if (k.ne.0) then ; a=(k*sumxy-sumx*sumy)/(k*sumx2-sumx**2) ; b=(sumy-a*sumx)/k ; endif
				end subroutine

			subroutine logtrendX2(Array)
				real(8) Array(:,:)
               !=- XY(1,i) - x, XY(3,i) - dx
               !=- XY(2,i) - y, XY(4,i) - dy
					XYdXdY(:,:)=0d0 ; N = size(Array(1,:))
					N_points = 0
					do i=1,N
                     if ( Array(2,i)>0 .and. Array(1,i).le.RightLimit .and. Array(1,i).ge.LeftLimit ) N_points = N_points + 1
                  enddo

               forall (i=1:N,j=1:2, Array(j,i)>0) XYdXdY(j,i)=log10(Array(j,i))
					forall (i=1:N,j=3:4, Array(j,i)>0) XYdXdY(j,i)= Array(j,i)

               mx = sum(XYdXdY(1,:))/N_points
               my = sum(XYdXdY(2,:))/N_points
if ( abs(mx) < 1 ) write(*,*) 'bad solution'
               forall (i=1:N, Array(1,i)>0) XYdXdY(1,i) = XYdXdY(1,i) - mx
               forall (i=1:N, Array(2,i)>0) XYdXdY(2,i) = XYdXdY(2,i) - my

               a_step   =  dsqrt((a_max-a_start)/N_grid_LT_a**2)
               b_step   =  dsqrt((b_max-b_start)/N_grid_LT_b**2)
               X2=1d150
					do j=-N_grid_LT_b,N_grid_LT_b
                  do i=-N_grid_LT_a,N_grid_LT_a

                     if (i.ne.0) va = i/abs(i)*( a_start + (i*a_step)**2 )
                     if (j.ne.0) vb = j/abs(j)*( b_start + (j*b_step)**2 )

                     if (i==0) va = a_start
                     if (j==0) vb = b_start

                     weights = 0d0
                     do k=1,N
                        if (XYdXdY(2,k).ne.0d0) then

                           Sx = XYdXdY(3,k)	;	if (Sx==0) Sx=1
                           Sy = XYdXdY(4,k)	;	if (Sy==0) Sy=1

                           weights = weights +  ( ( XYdXdY(2,k)-va*XYdXdY(1,k)-vb )  )**2 / Sy**2 / Sx**2

                           endif
                        enddo

                     if ( weights < X2 ) then
                        X2 = weights

								a = va
								b = vb +  my - mx*a

                        va_error_lower =  ( (i*a_step)**2 - ((i-1)*a_step)**2 )
                        vb_error_lower =  ( (j*b_step)**2 - ((j-1)*b_step)**2 )

                        va_error_upper =  ( ((i+1)*a_step)**2 - (i*a_step)**2 )
                        vb_error_upper =  ( ((j+1)*b_step)**2 - (j*b_step)**2 )

                        a_error_lower  = - va_error_lower
                        b_error_lower  = - vb_error_lower - mx*va_error_lower

                        a_error_upper = va_error_upper
                        b_error_upper = vb_error_upper + mx*va_error_upper

                        endif
                     enddo
                  enddo

               forall (i=1:N, Array(1,i)>0) XYdXdY(1,i) = XYdXdY(1,i) + mx
               forall (i=1:N, Array(2,i)>0) XYdXdY(2,i) = XYdXdY(2,i) + my

				end subroutine

			subroutine SORT2(Array)
				real(8) Array(:,:),buf(2)
					do i=1,size(Array(1,:))
						do j=i,size(Array(1,:))
							if (Array(1,j).gt.Array(1,i) .and. Array(1,j).ne.0) then;
								Buf(:)=Array(:,j);Array(:,j)=Array(:,i); Array(:,i)=Buf(:); endif
							enddo
						enddo
				end subroutine

			subroutine SORT(Array)
				real(8) Array(:),abuf
					do i=1,size(Array(:))
						do j=i,size(Array(:))
							if ( Array(j) .gt. Array(i) .and. Array(j).ne.0 ) then;
								aBuf=Array(i);Array(i)=Array(j); Array(j)=aBuf; endif
							enddo
						enddo
				end subroutine



			subroutine compute_medians(input_array,x_start_in,x_final_in,dx)
				real(8), intent(in) :: input_array(:,:),dx,x_start_in,x_final_in
				real(8) x_start, x_final
            real(8),allocatable,dimension(:,:) :: XY

					allocate(XY(size(input_array(:,1)),size(input_array(1,:))))

					medar(:,:)=0d0 ; XY(:,:)=0d0
					x_start=x_start_in
               x_final=x_final_in
					XY(:,:)=input_array(:,:)
               do i=1,size(input_array(1,:)) ; if (XY(1,i)==0) XY(1,i)=-99 ; enddo

					if (log_medians_parameter/=0) then
                     do i=1,size(input_array(1,:))
                        if (XY(1,i)/=0) XY(1,i)=log10(XY(1,i))
                        if (XY(1,i)==0) XY(1,i)=-99
                        enddo
                     x_start=log10(x_start_in)
                     x_final=log10(x_final_in)
                  endif

					call sort2(XY)
					k=0; sm=0 ; n=int((x_final-x_start)/dx)+1
					do ii=0,n
						m=0; mm=0; x_step_final = x_final-ii*dx; x_step_start = x_step_final-dx
						do i=1,size(input_array(1,:))
							if (XY(1,i)==-99) exit
							if ( XY(1,i) <= x_step_final .and. XY(1,i) >= x_step_start ) then  !   .and. XY(1,i)/=0
								if (m==0) m=i
								mm=mm+1
								endif
							end do
						mm=mm+m-1

						if (mm-m>2) then
							k=k+1;

							medar(1,k) = x_step_start + 0.5d0*(x_step_final-x_step_start)
							medar(2,k) = median(XY(2,m:mm))
							medar(3,k) = ( sum( (XY(2,m:mm)-medar(2,k))**2 )/(mm-m) )**0.5d0

							if ( log_medians_parameter /= 0 ) then
                        medar(1,k) = 1d1**medar(1,k)
                        medar(4,k) = ( sum( (1d1**XY(1,m:mm)-medar(1,k))**2 )/(mm-m) )**0.5d0
                        else
                           medar(4,k) = ( sum( (XY(1,m:mm)-medar(1,k))**2 )/(mm-m) )**0.5d0
                        endif

							endif
						enddo

!=-

						m=1 ; mm=0
                  if (log_medians_parameter==0) then
                     do while (XY(1,mm+1)>4.2) ; mm=mm+1 ; enddo
                     else
                        do while (1d1**XY(1,mm+1)>4.2) ; mm=mm+1 ; enddo !  ; write(*,*) m,1d1**XY(1,mm) ; pause
                  end if

						if (log_medians_parameter==0) then
                     add_median_x  = median(XY(1,m:mm))
                     add_median_dx = ( sum( (XY(1,m:mm)-add_median_x)**2 )/(mm-m) )**0.5d0
                     else
                        add_median_x  = median(1d1**XY(1,m:mm))
                        add_median_dx = ( sum( (1d1**XY(1,m:mm)-add_median_x)**2 )/(mm-m) )**0.5d0
                     endif

						add_median_y  = median(XY(2,m:mm))
						add_median_dy = ( sum( (XY(2,m:mm)-add_median_y)**2 )/(mm-m) )**0.5d0

!=-

						m=1 ; mm=1
                  if (log_medians_parameter==0) then
                     do i=1,size(input_array(1,:))
                        if ( XY(1,i)>0 .and. XY(1,i)<0.4 ) then
                           mm=i
                           if (m==1) m=i
                           endif
                           enddo

                        !do while

                     else
                        do i=1,size(input_array(1,:))
                        if ( 1d1**XY(1,i)>0 .and. 1d1**XY(1,i)<0.35 ) then
                           mm=i
                           if (m==1) m=i
                           endif
                           enddo !  ; write(*,*) m,1d1**XY(1,mm) ; pause
                     end if

						if (log_medians_parameter==0) then
                     add_median2_x  = median(XY(1,m:mm))
                     add_median2_dx = ( sum( (XY(1,m:mm)-add_median2_x)**2 )/(mm-m) )**0.5d0
                     else
                        add_median2_x  = median(1d1**XY(1,m:mm))
                        add_median2_dx = ( sum( (1d1**XY(1,m:mm)-add_median2_x)**2 )/(mm-m) )**0.5d0
                     endif

						add_median2_y  = median(XY(2,m:mm))
						add_median2_dy = ( sum( (XY(2,m:mm)-add_median2_y)**2 )/(mm-m) )**0.5d0

!=-

               deallocate(XY)
				end subroutine



			real(8) function median(Array)
				real(8) Array(:)
				real(8),allocatable,dimension(:) :: Barray
					allocate(Barray(size(Array)))
					Barray(:)=Array(:)
					call SORT(Barray)

					do i=1,size(Barray)
						if (Barray(i+1)==0 .or. i==size(Barray)) then
							if (mod(i,2)==1) then
								median=Barray(int(i/2)+1)
								else
									median=0.5d0*(Barray(i/2)+Barray(i/2+1))
								endif
							exit
							end if
               end do

               deallocate(Barray)
         end function



         subroutine EclCoord(RA,Dec,l2,b2) !=- from galactic to ecliptical -=1
            real(8) RA,RA2,Dec,l,l2,b,b2,RAl,Decl,l0,pi,r,d
               pi   = 4*datan(1d0)
               l0   = 32.93192 + 90d0
               RAl  = 192.858333
               Decl = 27.128333/180d0*pi
               b    = b2/180d0*pi
               l    = l2/180d0*pi
               Dec = dasin( dsin(Decl) * dsin(b) + dcos(Decl) * dcos(b) * dcos( l0/180d0*pi - l ))
               RA  = dasin( dcos(b)    * dsin( l0/180d0*pi - l )/dcos(Dec) )
               RA2 = ( dcos(Decl) * dsin(b) - dsin(Decl) * dcos(b) * dcos( l0/180d0*pi - l ) )/dcos(Dec)
               Dec = Dec*180d0/pi
               if (RA2<=0) RA = pi - RA
               RA  = RA*180d0/pi - RAl + 25.71667480468742
               if (RA <=0) RA = RA + 360d0
            end subroutine EclCoord



         subroutine GalCoord(RA,Dec,l,b) !=- from ecliptical to galactic -=1
            real(8) RA,Dec,l,l2,b,RAl,Decl,l0,pi,r,d
               pi   = 4*datan(1d0)
               l0   = 32.93192 + 90d0
               RAl  = 192.858333 / 180d0*pi
               Decl = 27.128333  / 180d0*pi
            r=RA ; d=Dec
               RA  = RA/180d0*pi
               Dec = Dec/180d0*pi
               b   = dasin( dsin(Dec) * dsin(Decl) + dcos(Dec) * dcos(Decl) * dcos( RA - RAl ) )
               l   = dasin( dcos(Dec) * dsin( RA - RAl ) / dcos(b)                             )
               l2  = ( dcos(Decl) * dsin(Dec) - dsin(Decl) * dcos(Dec) * dcos( RA - RAl ) )/dcos(b)

               if (l2<=0) l = pi - l
               l   = l0 - l*180d0/pi
               if (l <=0) l = l+360d0
               b   = b*180d0/pi
            RA=r ; Dec=d
            end subroutine GalCoord



         subroutine SphCoord(x,y,z,b,l,r) !=- from x,y,z to a,b,r -=1
            real(8) x,y,z,l,b,r  !=- b -> l , l -> b, valid -=1
               r = dsqrt( x**2 + y**2 + z**2 )
               l=0 ; b=0
               if (r.ne.0) then
                  l = datan( z/dsqrt( x**2 + y**2 ) )
                  if ( x>0  .and. y>=0 ) b = datan( y/x )
                  if ( x<0             ) b = pi + datan( y/x )
                  if ( x>0  .and. y<0  ) b = 2*pi+datan( y/x )
                  if ( x==0 .and. y>0  ) b = pi/2
                  if ( x==0 .and. y<0  ) b = 3d0*pi/2
                  endif
               l = l * 180d0/pi
               b = b * 180d0/pi
            end subroutine



         subroutine SphCoord_radians(x,y,z,l,b,r) !=- from x,y,z to a,b,r -=1
            real(8) x,y,z,l,b,r
               r = dsqrt( x**2 + y**2 + z**2 )
               l=0 ; b=0
               if (r.ne.0) then
                  l = datan( z/dsqrt( x**2 + y**2 ) )
                  if ( x>0  .and. y>=0 ) b = datan(y/x)
                  if ( x<0             ) b = pi+datan(y/x)
                  if ( x>0  .and. y<0  ) b = 2*pi+datan(y/x)
                  if ( x==0 .and. y>0  ) b = pi/2
                  if ( x==0 .and. y<0  ) b = 3d0*pi/2
                  endif
            end subroutine



         subroutine DescCoord(a,b,r,x,y,z) !=- from a,b,r to x,y,z -=1
            real(8) x,y,z,a,b,r
               x = r * dcos(a) * dcos(b)
               y = r * dcos(a) * dsin(b)
               z = r * dsin(a)
            end subroutine



   end module
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
