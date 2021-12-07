!=- fortran-libraries
!=- © Stanislav Shirokov, 2014-2020

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
	module math
		use global
		use GNUplot

      logical           :: approx_save_datafiles = .false.

      character(len)    :: columns_for_approx   = 'NaN'                       , & !=- x y | x y dy | x y dx dy
                           polylog_a_data       = 'polylog_a.dat'             , &
                           approx_datafile      = 'data_approx.dat'           , &
                           approx_coefficients  = 'approx_coefficients.dat'

      integer,parameter :: N_medians = 1d2 , approx_max_order = 10 , N_temp_approx = 1d3

      real(8),parameter :: pi = 4*datan(1d0)

		integer ::  N_grid_LT = 2d2, N_grid_LT_a = 2d2, N_grid_LT_b = 2d1, target_point_number , &
                  data_column_x , data_column_y , data_column_dx , data_column_dy, &

                  approximating_step , step_N_temp , approx_order = 2 , &
                  approx_minimum_order = 99

		real(8) :: RightLimit = 1d150, LeftLimit = -1d150,  &
         a_error_lower, b_error_lower, a_error_upper, b_error_upper, a_start = 1d-4, b_start = 1d-4, a_max = 1d1, b_max = 1d0, &
         mx, my, a_step, b_step, X2=1d150, va, vb, weights, Sx, Sy, &
         va_error_lower,vb_error_lower , va_error_upper, vb_error_upper, medar(4,N_medians), &
         add_median_left_border, add_median_right_border, x_step,x_step_start,x_step_final, &
         add_median_x,add_median_y,add_median_dx,add_median_dy, log_medians_parameter = 0, &
         add_median2_x,add_median2_y,add_median2_dx,add_median2_dy, &

         polylog_a(approx_max_order+1) = 0d0 , polylog_a_borders(2,approx_max_order+1) = 1d0 , approx_scaling_x = 1d2 , &
         extremums_polylog_a( approx_max_order+1 , approx_max_order+1 ) , &
         extremums_chi2(approx_max_order+1) = NaN, linear_approx(2) = 0d0 , minimum_chi2 = NaN, &
         temp_a(approx_max_order+1) = 0d0 , approx_minimum_chi2 = NaN , approx_minimum_polylog_a(approx_max_order+1) = NaN

      real(8),allocatable,dimension(:,:) :: XYdXdY

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
		contains

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine math_default

            columns_for_approx   = 'NaN'
            polylog_a_data       = 'polylog_a.dat'
            approx_datafile      = 'data_approx.dat'
            approx_coefficients  = 'approx_coefficients.dat'
            approx_minimum_order = 99
            approx_minimum_chi2  = NaN
            approx_minimum_polylog_a(:) = NaN

            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine approx_default

            polylog_a_data       = 'polylog_a.dat'
            approx_datafile      = 'data_approx.dat'
            approx_coefficients  = 'approx_coefficients.dat'

            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine approx_columns_checking
            integer :: x = 0 , y = 0 , dx = 0 , dy = 0

               if (columns_for_approx=='NaN') then
                  x=1 ; y=2
                  else

                     read(columns_for_approx,*,end=11,err=11) x , y
                     read(columns_for_approx,*,end=11) x , y , dy
                     read(columns_for_approx,*,end=11) x , y , dy , dx
                     11 continue ; if ( x==0 .or. y==0 ) write(*,*) '   scripting_error: approx_columns_checking: x=0 or y=0'

                  endif
               data_column_x  = x
               data_column_y  = y
               data_column_dx = dx
               data_column_dy = dy
            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine read_approx_catalog(data_file)
            character(len) :: data_file , line , catalog(200) = '0'
            integer i

            unit_1 = random_unit()
            open(unit_1,file=data_file,status='old')

               i=1
               read(unit_1,'(A4096)',end=22) line   !=- 1-st line reading (head of data table)
               do
                  read(unit_1,'(A4096)',end=22) line
                  read(line,*,end=11,err=11) catalog
                  11 continue

                  read(catalog( data_column_x ),*) XYdXdY( 1 , i )
                  read(catalog( data_column_y ),*) XYdXdY( 2 , i )
                  if (data_column_dx/=0) read(catalog( data_column_dx ),*) XYdXdY( 3 , i )
                  if (data_column_dy/=0) read(catalog( data_column_dy ),*) XYdXdY( 4 , i )
                  i=1+i

                  enddo ; 22 continue
               close(unit_1)

            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         real(8) function polylog( a , x , order )
            integer i,order
            real(8) x,a(order+1)

            polylog = 0d0
            do i=0,order
               polylog = polylog + a(1+i) * log10 (x) ** i
               !write(*,*) temp_a(1:2) , temp_a(2) * log10 (x)
               enddo

            polylog = polylog + 5d0*log10 (x) !+ 25d0 !- 1d1

            end function

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         real(8) function chi2_polylog(order)
            integer i,order,N
            real(8) polylog_model

            N = size(XYdXdY(1,:))

            ! XYdXdY( 4 , : ) = 1d0 !=- test

            chi2_polylog = 0d0
            do i=1,N
                  if (XYdXdY( 4 , i )==0) XYdXdY( 4 , i ) =1d0
                  polylog_model = polylog( temp_a , XYdXdY( 1 , i ) , order )
               chi2_polylog = chi2_polylog + ( XYdXdY( 2 , i ) - polylog_model ) ** 2 / &
                  XYdXdY( 4 , i )**2 / dabs(polylog_model) / N
               enddo

            end function

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine write_polylog_a( order )
            integer order ; real(8) chi2

               chi2 = chi2_polylog( order )
               chi2 = log10(chi2)

               unit_2 = random_unit()
            open(unit_2,file=polylog_a_data,status='old',position='append')

            theformat = '(' // trim(adjustl(inttostr(approx_max_order))) // '(E20.8))'
            write(unit_2,theformat) chi2 , temp_a(1:order+1)

               if ( chi2 < minimum_chi2 ) then
                  minimum_chi2   = chi2
                  polylog_a(:)   = temp_a(:)
                  endif

               close(unit_2)
            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine create_data_polylog_a_temp
               unit_2 = random_unit()
            open(unit_2,file=polylog_a_data,status='replace')
               close(unit_2)
            end subroutine

         subroutine write_data_approx( order )
            integer i,order
            real(8) :: temp_a(7) = (/  33.158168604381189         , &
                                       6.2594070261452195E-002    , &
                                       2.7994852249554239E-002    , &
                                       2.7962159075643332E-002    , &
                                       1.8916688159884384E-002    , &
                                       3.0028150691020179E-002    , &
                                       -1.2242181029368278E-002   /)

               !=- polylog_a(1:4) = (/ 44.1 , 6.28 , 0.57 , 0.08 /)

               if (approx_save_datafiles) then
                  approx_datafile = trim(file_path(approx_datafile)) // trim(adjustl( inttostr(order) )) // '_' // &
                     trim(name_ext(approx_datafile))
                  approx_coefficients = trim(file_path(approx_coefficients)) // trim(adjustl( inttostr(order) )) // '_' // &
                     trim(name_ext(approx_coefficients))
               end if

               unit_3 = random_unit()
            open(unit_3,file=approx_datafile,status='replace')
               do i=1,size(XYdXdY(1,:))
                  !write(unit_3,*) XYdXdY(1,i)/approx_scaling_x, XYdXdY(2:4,i) , polylog( polylog_a , XYdXdY( 1 , i ) , order ) , &
                  !   linear_approx(1) + log10(XYdXdY(1,i))*(linear_approx(2)) , dm_LCDM( XYdXdY(1,i)/approx_scaling_x , 7d1 ), &
                  !   polylog( temp_a , XYdXdY( 1 , i ) , size(temp_a)-1 ) , 43d0 + 5d0*log10( XYdXdY(1,i)/approx_scaling_x )
                  enddo
               close(unit_3)

               unit_4 = random_unit()
            open(unit_4,file=approx_coefficients,status='replace')
               write(unit_4,*) order , 1d1**minimum_chi2 , polylog_a(1:order+1)
               close(unit_4)

            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine approx_polylog(data_file)
            character(len) data_file
            real(8) :: shift(approx_max_order+1) = NaN , center , target_shift , accuracy = 1d-4
            integer i,N,order , condition

               N = file_volume(data_file) - 1
            allocate ( XYdXdY( 4 , N ) ) ; XYdXdY(:,:)=1d0

                  call approx_columns_checking
               call read_approx_catalog(data_file)

            XYdXdY(1,:) = approx_scaling_x * XYdXdY(1,:)

               call logtrend_x(XYdXdY(1:2,:),linear_approx(2),linear_approx(1),0d0*approx_scaling_x,0.1*approx_scaling_x)
                  !=- write(*,*) '   LA: ' , linear_approx(1) + 25 , linear_approx(2) - 5

                  if (approx_order<0) approx_order=0
                  if (approx_order>approx_max_order) approx_order=approx_max_order
               do j=1,approx_order

                  call approx_default

                  order                = j
                  minimum_chi2         = NaN
                  approximating_step   = 1
                  condition            = 0

                  polylog_a_borders(1,1) = 0.9 * ( linear_approx(1) ) !=- + 25
                  polylog_a_borders(2,1) = 1.1 * ( linear_approx(1) ) !=- + 25
                  polylog_a_borders(1,2:approx_max_order) = -5d0
                  polylog_a_borders(2,2:approx_max_order) = 5d0

                  do while ( condition < 3 )

                     call create_data_polylog_a_temp

                     do i = 1,N_temp_approx
                           step_N_temp = i
                           call random_polylog_a(N_temp_approx)
                        call write_polylog_a( order )
                        enddo

                     do i=1,order+1
                        center         = polylog_a_borders(1,i) + (polylog_a_borders(2,i) - polylog_a_borders(1,i))*0.5d0
                        target_shift   = ( polylog_a_borders(2,i) - center ) * 0.9d0
                        shift(i)       = polylog_a(i) - center

                        if ( shift(i) < target_shift ) then
                           polylog_a_borders(1,i) = polylog_a(i) - target_shift
                           polylog_a_borders(2,i) = polylog_a(i) + target_shift
                           else
                              polylog_a_borders(:,i) = polylog_a_borders(:,i) + shift(i)
                           end if
                        enddo

                        call plot_polylog_a( order )

                           extremums_polylog_a(j,:)   = polylog_a(:)
                           extremums_chi2(j)          = minimum_chi2
                           approximating_step         = 1 + approximating_step

                           if ( minimum_chi2 < approx_minimum_chi2 ) then
                              approx_minimum_chi2         = minimum_chi2
                              approx_minimum_polylog_a(:) = polylog_a(:)
                              approx_minimum_order        = j
                              endif

                        if (sum(dabs(shift(1:order+1))) < accuracy*(order+1)) condition = 1 + condition

                     enddo
                  call write_data_approx( order )
                  call plot_data_approx( data_file , order )
                  end do
                  write(*,*) '   math_approx: the best order:', approx_minimum_order
                  write(*,*) '   math_approx: the best coefficients:', &
                     approx_minimum_polylog_a(1) , approx_minimum_polylog_a(2:approx_minimum_order+1)
                  write(*,*) '   math_approx: the best chi2:', 1d1**approx_minimum_chi2 , approx_minimum_chi2
               deallocate(XYdXdY)
            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine plot_polylog_a( order )
            integer order
               call clear_plot_fields

               graph_name = trim(adjustl(inttostr(order))) // '_' // trim(adjustl(inttostr(approximating_step)))

               pngterm  = 'set term pngcairo enhanced font "Verdana,20" size 1024, 1024'
               epsterm  = 'set term postscript portrait enhanced color font "Verdana,10" size 6, 6'

               GNUfields (legend)   =  'off'

               GNUfields (add1)     =  'set palette rgb 7,5,15'
               GNUfields (add2)     =  'set view map'
               GNUfields (add3)     =  'set cblabel "log {/Symbol C}^2"'
               GNUfields (add4)     =  'set format cb "%2.0t{/Symbol \327}10^{%L}"'

               GNUfields (title)    =  '{/Symbol C}^2-phase plane a_0--a_1, approximating step ' // &
                  trim(adjustl(inttostr(approximating_step))) // ' (CosMod v1.0)'
               GNUfields (xlabel)   =  'a_0'
               GNUfields (ylabel)   =  'a_1'

               GNUfields (plot1)    =  'splot "' // trim((polylog_a_data)) // &
                                          '" u 2:3:1 w p lw 2 pt 7 ps 0.5 dt 1 palette'

               GNUfields(extention_out_figure)='eps' ; call plot(polylog_a_data)

            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine plot_data_approx( source_data , order )
            integer order
            character(len) source_data
            character(2) :: delta_column = '7'

               call clear_plot_fields

               graph_name = trim(adjustl(inttostr(order))) // '_' // trim(name(source_data))

               pngterm  = 'set term pngcairo enhanced font "Verdana,20" size 1400, 1050'
               epsterm  = 'set term postscript enhanced color font "Verdana,14" size 8.5, 6.3'

               GNUfields (legend)   =  'right bottom'

               GNUfields (title)    =  ' approximation of the ' // trim(name(source_data)) // &
                  ' data (CosMod v1.0)'
               GNUfields (xlabel)   =  'z'
               GNUfields (ylabel)   =  '{/Symbol m}'
               GNUfields (logscale) =  'set logscale x'

               GNUfields (ls1)      =  'set linestyle 1 lw 1 pt 7 ps 0.7 lt rgb "blue"'
               GNUfields (ls2)      =  'set linestyle 2 lw 2 pt 7 ps 0.7 lt rgb "black"'
               GNUfields (ls3)      =  'set linestyle 3 lw 4 pt 7 ps 0.7 lt rgb "yellow"'
               GNUfields (ls4)      =  'set linestyle 4 lw 3 pt 7 ps 0.7 lt rgb "red"'
               GNUfields (ls5)      =  'set linestyle 5 lw 3 pt 7 ps 0.7 lt rgb "dark-cyan"'

               GNUfields (plot1)    =  '"' // trim((approx_datafile)) // '" u 1:2:4 w yerrorb ls 1'
               GNUfields (plot2)    =  '"' // trim((approx_datafile)) // '" u 1:5   w l ls 3'
               GNUfields (plot3)    =  '"' // trim((approx_datafile)) // '" u 1:6   w l ls 2'
               GNUfields (plot4)    =  '"' // trim((approx_datafile)) // '" u 1:7   w l ls 4'

               GNUfields (title1)   =  'the ' // trim(name(source_data)) // ' data'
               GNUfields (title2)   =  'polylog approximation of power ' // trim( adjustl( inttostr(order) ) )
               GNUfields (title3)   =  'first linear approximation (FLA)'
               GNUfields (title4)   =  'LCDM (70,0.7)'


               GNUfields(extention_out_figure)='png' ; call plot(approx_datafile)

               graph_name = 'delta_LCDM_' // trim(adjustl(inttostr(order))) // '_' // trim(name(source_data))

               GNUfields (ylabel)   =  '{/Symbol D}{/Symbol m} = {/Symbol m} - {/Symbol m}_{LCDM(70,0.7)}'

                  delta_column = '7'
               GNUfields (plot1)    =  '"' // trim((approx_datafile)) // '" u 1:($2-$' // delta_column // '):4 w yerrorb ls 1'
               GNUfields (plot2)    =  '"' // trim((approx_datafile)) // '" u 1:($5-$' // delta_column // ')   w l ls 3'
               GNUfields (plot3)    =  '"' // trim((approx_datafile)) // '" u 1:($6-$' // delta_column // ')   w l ls 2'
               GNUfields (plot4)    =  '"' // trim((approx_datafile)) // '" u 1:($7-$' // delta_column // ')   w l ls 4'
               GNUfields (plot5)    =  '"' // trim((approx_datafile)) // '" u 1:($8-$' // delta_column // ')   w l ls 5'

               GNUfields (title5)   =  'testing model'

               GNUfields(extention_out_figure)='png' ; call plot(approx_datafile)

               graph_name = '5log_x_' // trim(adjustl(inttostr(order))) // '_' // trim(name(source_data))

               GNUfields (ylabel)   =  '{/Symbol D}{/Symbol m} = {/Symbol m} - {/Symbol m}_{55+5log x}'

                  delta_column = '8'
               GNUfields (plot1)    =  '"' // trim((approx_datafile)) // '" u 1:($2-$' // delta_column // '):4 w yerrorb ls 1'
               GNUfields (plot2)    =  '"' // trim((approx_datafile)) // '" u 1:($5-$' // delta_column // ')   w l ls 3'
               GNUfields (plot3)    =  '"' // trim((approx_datafile)) // '" u 1:($6-$' // delta_column // ')   w l ls 2'
               GNUfields (plot4)    =  '"' // trim((approx_datafile)) // '" u 1:($7-$' // delta_column // ')   w l ls 4'
               GNUfields (plot5)    =  '"' // trim((approx_datafile)) // '" u 1:($8-$' // delta_column // ')   w l ls 5'

               GNUfields(extention_out_figure)='png' ; call plot(approx_datafile)

               graph_name = 'delta_linear_' // trim(adjustl(inttostr(order))) // '_' // trim(name(source_data))

               GNUfields (ylabel)   =  '{/Symbol D}{/Symbol m} = {/Symbol m} - {/Symbol m}_{FLA}'
               GNUfields (legend)   =  'left top'

                  delta_column = '6'
               GNUfields (plot1)    =  '"' // trim((approx_datafile)) // '" u 1:($2-$' // delta_column // '):4 w yerrorb ls 1'
               GNUfields (plot2)    =  '"' // trim((approx_datafile)) // '" u 1:($5-$' // delta_column // ')   w l ls 3'
               GNUfields (plot3)    =  '"' // trim((approx_datafile)) // '" u 1:($6-$' // delta_column // ')   w l ls 2'
               GNUfields (plot4)    =  '"' // trim((approx_datafile)) // '" u 1:($7-$' // delta_column // ')   w l ls 4'
               GNUfields (plot5)    =  '"' // trim((approx_datafile)) // '" u 1:($8-$' // delta_column // ')   w l ls 5'

               GNUfields(extention_out_figure)='png' ; call plot(approx_datafile)

            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine random_polylog_a(N_temp)
            real(8) x , stretching , center
            integer i , N_temp

            do i=1,approx_max_order

               call random_number(x)
               x=x-0.5d0

            center     = polylog_a_borders(1,i) + (polylog_a_borders(2,i) - polylog_a_borders(1,i))*0.5d0
            stretching = x*step_N_temp**2 * ( polylog_a_borders(2,i) - polylog_a_borders(1,i) ) / N_temp**2


               temp_a(i) = center + stretching
               !temp_a(i) = polylog_a_borders(1,i) + a * ( polylog_a_borders(2,i) - polylog_a_borders(1,i) )

               end do
            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         real function distance(x,y,z)
            real(8) x,y,z
               distance = ( x**2 + y**2 + z**2 )**0.5d0
            end function

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

			integer function sign_real(x) ; real(8) x;
					sign_real=1 ; if (x<0) sign_real=-1
				end function

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

			real(8) function mydexp(x) 	; real(8) x;
					n=1 ; mydexp = ( 1d0 + x/1d5 )**1d5
				end function

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

			real(8) function factorial(n)	; integer n;
					factorial=1d0; do while (n>0); factorial=factorial*n;	n=n-1; enddo
				end function

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

			real(8) function mydsin(x) 	; real(8) x;
					n=1; mydsin=x; do; mydsin=mydsin+(-1d0)**n * x**(2*n+1) / factorial (2*n+1)
						n=n+1 ; if (x**(2*n+1) / factorial (2*n+1) < 1d-19) exit ; end do
				end function

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

			subroutine logtrend(XY,a,b,Lb,Rb)
				real(8) XY(2,N_GRB),G(2,N_GRB),a,b,sumx,sumx2,sumy,sumxy,Lb,Rb;
					G(:,:)=0d0; forall (i=1:N_GRB,j=1:2, XY(j,i)>0) G(j,i)=log10(XY(j,i))
					j=0;k=0;sumx=0d0;sumx2=0d0;sumy=0d0;sumxy=0d0; N=0
					do j=1,N_GRB; if (G(2,j).ne.0d0 .and. G(1,j).le.log10(Rb) .and. G(1,j).ge.log10(Lb)) then!		write(*,*) N,G(1,j),G(2,j),XY(1,j),XY(2,j)
						sumx=sumx+G(1,j) ; sumx2=sumx2+G(1,j)**2d0 ; sumy=sumy+G(2,j) ; sumxy=sumxy+G(1,j)*G(2,j) ; k=k+1 ; endif ; enddo
					if (k.ne.0) then ; a=(k*sumxy-sumx*sumy)/(k*sumx2-sumx**2) ; b=(sumy-a*sumx)/k ; endif
				end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

			subroutine logtrend_x(XY,a,b,Lb,Rb)
            integer N,i,j,k
				real(8) XY(:,:) ,a,b,sumx,sumx2,sumy,sumxy,Lb,Rb
				real(8),allocatable,dimension(:,:) :: G

               N = size(XY(2,:))
            allocate(G(2,N))
					G(:,:)=XY(:,:); forall (i=1:N, XY(1,i)>0) G(1,i)=log10(XY(1,i))

					k=0;sumx=0d0;sumx2=0d0;sumy=0d0;sumxy=0d0
					do j=1,N; if (G(2,j).ne.NaN .and. G(1,j).le.log10(Rb) .and. G(1,j).ge.log10(Lb)) then!		write(*,*) N,G(1,j),G(2,j),XY(1,j),XY(2,j)
						sumx=sumx+G(1,j) ; sumx2=sumx2+G(1,j)**2d0 ; sumy=sumy+G(2,j) ; sumxy=sumxy+G(1,j)*G(2,j) ; k=k+1 ; endif ; enddo
					if (k.ne.0) then ; a=(k*sumxy-sumx*sumy)/(k*sumx2-sumx**2) ; b=(sumy-a*sumx)/k ; endif

               deallocate(G)
				end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

			subroutine logtrendX2(Array)
				real(8) Array(:,:)
               !=- XY(1,i) - x, XY(3,i) - dx
               !=- XY(2,i) - y, XY(4,i) - dy
					XYdXdY(:,:)=0d0 ; N = size(Array(1,:))
					target_point_number = 0
					do i=1,N
                     if ( Array(2,i)>0 .and. Array(1,i).le.RightLimit .and. Array(1,i).ge.LeftLimit ) &
                        target_point_number = target_point_number + 1
                  enddo

               forall (i=1:N,j=1:2, Array(j,i)>0) XYdXdY(j,i)=log10(Array(j,i))
					forall (i=1:N,j=3:4, Array(j,i)>0) XYdXdY(j,i)= Array(j,i)

               mx = sum(XYdXdY(1,:))/target_point_number
               my = sum(XYdXdY(2,:))/target_point_number
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

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

			subroutine SORT2(Array)
				real(8) Array(:,:),buf(2)
					do i=1,size(Array(1,:))
						do j=i,size(Array(1,:))
							if (Array(1,j).gt.Array(1,i) .and. Array(1,j).ne.0) then;
								Buf(:)=Array(:,j);Array(:,j)=Array(:,i); Array(:,i)=Buf(:); endif
							enddo
						enddo
				end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

			subroutine SORT(Array)
				real(8) Array(:),abuf
					do i=1,size(Array(:))
						do j=i,size(Array(:))
							if ( Array(j) .gt. Array(i) .and. Array(j).ne.0 ) then;
								aBuf=Array(i);Array(i)=Array(j); Array(j)=aBuf; endif
							enddo
						enddo
				end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

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

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

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

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine EclCoord(RA,Dec,l2,b2) !=- from galactic to 2-Equatorial ones -=1
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

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine GalCoord(RA,Dec,l,b) !=- from 2-Equatorial to galactic ones -=1
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

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

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

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

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

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine DescCoord(a,b,r,x,y,z) !=- from a,b,r to x,y,z -=1
            real(8) x,y,z,a,b,r
               x = r * dcos(a) * dcos(b)
               y = r * dcos(a) * dsin(b)
               z = r * dsin(a)
            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine create_logscale_grid(left,right,N_segments,Array)
            integer i , N_segments
            real(8) left,right,Array(N_segments)
               if (N_segments==1) then
                  Array(1) = left
                  else
                     forall(i=1:N_segments) Array(i) = 1d1 ** (log10(left * 1.0) + (i - 1) &
                        * (log10(right * 1.0) - log10(left * 1.0)) / (N_segments - 1)) + 0.5
                  endif
            end

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine create_logscale_grid_int(left,right,N_segments,Array)
            integer i , left , right , N_segments , Array(N_segments)
               if (N_segments==1) then
                  Array(1) = left
                  else
                     forall(i=1:N_segments) Array(i) = floor(1d1 ** (log10(left * 1.0) + (i - 1) &
                        * (log10(right * 1.0) - log10(left * 1.0)) / (N_segments - 1)) + 0.5)
                  endif
            end

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine create_linear_grid(left,right,N_segments,Array)
            integer i , N_segments
            real(8) left,right,Array(N_segments)
               if (N_segments==1) then
                  Array(1) = left
                  else
                     forall(i=1:N_segments) Array(i) = left + (i - 1) * (right - left) / (N_segments - 1)
                  endif
            end

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine create_linear_grid_int(left,right,N_segments,Array)
            integer i , N_segments,left,right,Array(N_segments)
               if (N_segments==1) then
                  Array(1) = left
                  else
                     forall(i=1:N_segments) Array(i) = floor( 1d0*D_left + (i - 1) * (right - left) / (N_segments - 1))
                  endif
            end

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

   end module
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
