!=- fortran-libraries
!=- © Stanislav Shirokov, 2014-2020

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
	module cosmology
		use global
		use GNUplot
		use math

      logical ::  RDR_calculating = .true.   , &   !=- .false.
                  RDR_error       = .true.

		integer,parameter :: N_models = 9   , &
                  N_RDR_grid        = 1d3 , &   !=- graph grid
                  MS_grid           = 2d2 , &
                  MS_model_count    = 1d2

		integer ::  N_cosmogrid       = 2d2 , &
                  N_RDR_grid_mode   = 3   , &   !=- progression power
                  RDR_type          = 1   , &   !=- 0 is LCDM(0.3,70) , 1 is LCDM(\O_m,H_0) , 2 is wCDM(w,\O_m,H_0,O_k)
                  MS_real_model_count

		real(8) ::	c   = 2.99792458d10	, &	!=- cm/s
                  c_km_s = 2.99792458d5, &   !=- km/s
						H_0 = 70d0           , &	!=- km/s/Mpc
						H70 = 7d1            , &
						H100 = 1d2           , &
						dl  = 4421.71767d0   , &	!=- c/H_0/1d5 in GCS

					 L_sun = 3.827d33       , &   !=- erg/s
					 M_sun = 4.83d0         , &

						q_0 = 0.5d0				, &
						O_m = 0.308d0 			, &
						O_w = 1d0				, &
						O_k = 0d0            , &
					 z_max = 2d1				, &
             MS_z_max = 2d1            , &
            RDR_z_max = 2d1            , &

					 O_3	 = 0.3d0				, &
					 O_7	 = 0.7d0          , &

					 MS_H_0 , MS_O_v , MS_O_k , MS_w , MS_O_m_bug

      character(length) :: model_set (MS_model_count,2) = ' model ' , &
                        model_set_path = 'model_set.dat' , model_name , MS_columns(3)='' , &
                        MS_fields(30) , MS_parameters(30)

		real(8) ::  LumDist(N_models) = 1d0, MetrDist(N_models)= 1d0, DistMod(N_models)= 1d0, & !=- z for FCM > 0.0025
                  dz, x_max, I1, RDR(2,N_RDR_grid), &

                  MS_table ( 1+3*MS_model_count , MS_grid ) = 0d0 , &

                  vr_H0 , vr_distance , vr_velosity , vr_redshift

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
		contains



         real(8) function visible_redshift(x,y,z,vx,vy,vz)  !=- v = cz = Hr
            real(8) x,y,z,vx,vy,vz  !=- in units of kpc

            vr_distance = distance(x,y,z)          !=- in units of Mpc (because function `fun_from_R_to_z` is normed in Mpc scale)
            vr_velosity = (x+y+z)       !=- in units of km/s

            vr_redshift = fun_from_R_to_z( vr_distance )
            vr_redshift = vr_redshift + vr_velosity / c_km_s

            visible_redshift = vr_redshift

            end function visible_redshift



         subroutine make_model_set

            !=-   write(*,*) ' the model_set format:'
            !=-   write(*,*) ' wCDM H_0  w  Omega_v  Omega_k'
            !=-   write(*,*) ' FF H_0'
            !=-   write(*,*) ' TL H_0'

               model_set ( 1 , 1 ) = ' wCDM  70    -1 0.7   0.0 '
               model_set ( 1 , 2 ) = ' the standard LCDM model '

               model_set ( 2 , 1 ) = ' wCDM  38.4  -1 0.15  -0.15 '
               model_set ( 2 , 2 ) = ' the dVMS model '

               model_set ( 3 , 1 ) = ' wCDM  70    -1 1.0   0.0 '
               model_set ( 3 , 2 ) = ' the PV model '

               !model_set ( 7 , 1 ) = ' TL    70     '
               !model_set ( 7 , 2 ) = ' the TL model '

               !model_set ( 6 , 1 ) = ' wCDM  70    -1 0.8   0.0 '
               !model_set ( 6 , 2 ) = ' the LCDM {/Symbol W}_{/Symbol L} = 0.8 model '

               model_set ( 4 , 1 ) = ' wCDM  70    -1 0.9   0.0 '
               model_set ( 4 , 2 ) = ' the LCDM {/Symbol W}_{/Symbol L} = 0.9 model '

               model_set ( 5 , 1 ) = ' FF  70    -1 0.9   0.0 '
               model_set ( 5 , 2 ) = ' the FF model '

               call compute_model_set
               call write_model_set
               call plot_model_set

               end subroutine make_model_set



         subroutine compute_model_set
            integer i,j ; MS_table(:,:)=0d0

               inquire( file = model_set_path , exist = file_exists )

         if (.not. file_exists ) then

               do i=1,MS_model_count

                     read(model_set(i,1),*) model_name
                  if ( model_name /= 'model' ) then

                        do j=1,MS_grid

                           if (i==1) MS_table ( 1 , j ) = MS_z_max / MS_grid**2 * j**2 !(i+0.01d0)

                              select case (model_name)

                                 case('wCDM')

                                     MS_H_0 = 0d0 ; MS_w=0d0 ; MS_O_v=0d0 ; MS_O_k=0d0
                                    read(model_set(i,1),*) model_name , MS_H_0 , MS_w , MS_O_v , MS_O_k

                                    MS_O_m_bug = 1 + MS_O_k - MS_O_v !=-   the Friedmann equation
                           !write(*,*) 'MS_compute: Omega_m = ', MS_O_m_bug
                                    MS_table ( 2+3*(i-1) , j ) =  &
                                        r_wCDM( MS_table ( 1 , j ) , MS_H_0 , MS_w , MS_O_m_bug , MS_O_k )
                                    MS_table ( 3+3*(i-1) , j ) =  &
                                        d_wCDM( MS_table ( 1 , j ) , MS_H_0 , MS_w , MS_O_m_bug , MS_O_k )
                                    MS_table ( 4+3*(i-1) , j ) =  &
                                       dm_wCDM( MS_table ( 1 , j ) , MS_H_0 , MS_w , MS_O_m_bug , MS_O_k )

                                 case('FF')
                                    read(model_set(i,1),*) model_name , MS_H_0

                                    MS_table ( 2+3*(i-1) , j ) =  &
                                        r_FCM( MS_table ( 1 , j ) , MS_H_0 )
                                    MS_table ( 3+3*(i-1) , j ) =  &
                                        d_FCM( MS_table ( 1 , j ) , MS_H_0 )
                                    MS_table ( 4+3*(i-1) , j ) =  &
                                       dm_FCM( MS_table ( 1 , j ) , MS_H_0 )

                                 case('TL')
                                    read(model_set(i,1),*) model_name , MS_H_0

                                    MS_table ( 2+3*(i-1) , j ) =  &
                                        r_TL( MS_table ( 1 , j ) , MS_H_0 )
                                    MS_table ( 3+3*(i-1) , j ) =  &
                                        d_TL( MS_table ( 1 , j ) , MS_H_0 )
                                    MS_table ( 4+3*(i-1) , j ) =  &
                                       dm_TL( MS_table ( 1 , j ) , MS_H_0 )

                                 case default
                                    write(*,*) 'MS_compute: unknown model mask: ' // trim(model_name)
                                 end select

                           enddo
                     else
                        MS_real_model_count = i-1
                        write(*,*) 'MS_compute: model_count = ', MS_real_model_count
                        exit
                     end if

                  end do
            endif
            end subroutine compute_model_set



         subroutine write_model_set
            integer i
            character(length) title(MS_real_model_count)

               inquire( file = model_set_path , exist = file_exists )

         if (MS_model_count>0 .and. .not. file_exists ) then

            MS_columns(1) = ':metric_r'
            MS_columns(2) = ':luminosity_d_L'
            MS_columns(3) = ':distance_modulus'

            open(unit_1,file=model_set_path,status='replace')

               line=' '
               do i=1,MS_real_model_count
                  if (i>9) line=''
                  title(i) = 'model-' // trim(adjustl(inttostr(i))) // ': ' // trim(line) // trim(model_set(i,2))
                  write(unit_1,'(A2,A200)') '# ' , title(i)
                  end do

               theformat = '' ; theformat = '(A2,A23,' // trim(adjustl(inttostr(1+3*MS_real_model_count))) // &
                  '(A10,A10))'
               write(unit_1,theformat) '# ',' ', &
                  ( title(i),' ',title(i),' ',title(i),' ',i=1,MS_real_model_count )

               theformat = '' ; theformat = '(A2,A15,A8,' // trim(adjustl(inttostr(1+3*MS_real_model_count))) // '(i3,A17))'
               write(unit_1,theformat) '# ','1:redshift',' ',( 1+i,MS_columns(fix_index(i,3)),i=1,3*MS_real_model_count )

               theformat = '' ; theformat = '(' // trim(adjustl(inttostr(1+3*MS_real_model_count))) // '(E20.8))'
               write(unit_1,theformat) MS_table(1:1+3*MS_real_model_count,:)

               close(unit_1)
               else
                  if ( MS_model_count==0 ) write(*,*) 'MS: critical error'
            endif
            end subroutine write_model_set



         integer function fix_index( in_index , max_border )
            integer in_index , out_index , max_border
            out_index = in_index
            do while ( out_index > max_border )
               out_index = out_index - max_border
               enddo
               fix_index = out_index
            end function fix_index



         subroutine plot_model_set
            integer i,j,k
               call clear_plot_fields

               GNUfields(logscale)  =  'set logscale'
               GNUfields(format_y)  =  'set format y "10^{%L}"'
               GNUfields(xrange)    =  'set xrange [0.01:*]'

               GNUfields(legend)    =  'set key right bottom'
               GNUfields(xlabel)    =  'set xlabel "redshift z"'

               MS_fields(1) = 'r(z), [Mpc]'
               MS_fields(2) = 'd_L(z), [Mpc]'
               MS_fields(3) = '{/Symbol m}(z)'

               MS_fields(4) = 'log10 r(z)/r_{model-1}(z)'
               MS_fields(5) = 'log10 d_L(z)/d_{L,model-1}(z)'
               MS_fields(6) = 'log10 {/Symbol m}(z) - mu_{model-1}(z)'

               MS_fields(7) = 'log10 {/Symbol D}r(z)'
               MS_fields(8) = 'log10 {/Symbol D}d_L(z)'
               MS_fields(9) = 'log10 {/Symbol D}{/Symbol m}(z)'

               MS_fields(10) = 'r(z)'
               MS_fields(11) = 'd_L(z)'
               MS_fields(12) = 'mu(z)'

               MS_fields(13) = 'delta_r(z)'
               MS_fields(14) = 'delta_d_L(z)'
               MS_fields(15) = 'delta_mu(z)'

               do k=1,2
                  do j=1,3

                     GNUfields(title)  = 'set title "'// trim(MS_fields(j+6*(k-1))) //' for different models"'

                     if (j==3 .and. k/=2) then
                        GNUfields(logscale)  = '#set logscale'
                        GNUfields(format_y)  = '#set format y "10^{%L}"'
                        endif

                     do i=1,MS_real_model_count
                        select case (k)
                           case (1)
                              GNUfields(ylabel) = 'set ylabel "' // trim(MS_fields(j+6*(k-1))) // '"'

                              GNUfields (plot1+(i-1)*2) = ', "' // trim(slashfix(model_set_path)) // &
                                 '" u 1:' // trim(inttostr(1+j+3*(i-1))) // ' w l ls ' // trim(inttostr(i))
                              GNUfields (title1+(i-1)*2) = ' title "' // trim( model_set (i,2) ) // '"'
                           case(2)
                              GNUfields(yrange) =  'set yrange [-0.5:1]'
                              if (j==3) GNUfields(yrange) =  'set yrange [-0.025:0.05]'

                              GNUfields(ylabel) = 'set ylabel "' // trim(MS_fields(j+6*(k-1))) // ' = ' // &
                                 trim(MS_fields(j+3*(k-1))) // '"'

                              GNUfields (plot1+(i-1)*2) = ', "' // trim(slashfix(model_set_path)) // &
                                 '" u 1:(log10($' // trim(inttostr(1+j+3*(i-1))) // '/$' // &
                                 trim(inttostr(1+j)) // ')) w l ls ' // trim(inttostr(i))
                              GNUfields (title1+(i-1)*2) = ' title "' // trim( model_set (i,2) ) // '"'
                           end select
                        end do
                     GNUfields (plot1) = GNUfields (plot1)(2:length)

                     !graph_name = MS_fields(j+9+3*(k-1)) // '-logscale'
                     !   GNUfields(extention_out_figure)='#eps' ; call plot(model_set_path)
                     graph_name = MS_fields(j+9+3*(k-1)) // '-logscale'
                        GNUfields(extention_out_figure)='#png' ; call plot(model_set_path)
                     enddo
                  enddo
            end subroutine plot_model_set



      subroutine redshift_distance_relation  !=- RDR( redshift , luminosity distance )
         if (RDR_calculating) then
            do i=1,N_RDR_grid
               z = RDR_z_max**(1d0/N_RDR_grid_mode)/N_RDR_grid * i ; z = z**N_RDR_grid_mode
               RDR(1,i) = z
               RDR(2,i) = fun_from_z_to_R(z,H_0)
               end do
               RDR_calculating = .false.
            endif
         end subroutine

		real(8) function fun_from_R_to_z(luminosity_distance)
			real(8) luminosity_distance ; fun_from_R_to_z = 0d0
			   call redshift_distance_relation
            do i=1,N_RDR_grid
               if ( RDR(2,i) <= luminosity_distance ) fun_from_R_to_z = RDR(1,i)
               end do
			end function

		real(8) function fun_from_z_to_R(redshift,H0)
			real(8) redshift,H0 ; fun_from_z_to_R = 0d0
            select case (RDR_type)
                  case(0)
                     fun_from_z_to_R = D_wCDM(redshift,H70,-1d0,O_3,0d0)    !=- LCDM(0.3,0.7,H_0)
                  case(1)
                     fun_from_z_to_R = D_wCDM(redshift,H0,-1d0,O_m,0d0)    !=- LCDM(\O_m,H_0)
                  case(2)
                     fun_from_z_to_R = D_wCDM(redshift,H0,w,O_m,O_k)       !=- wCDM()
                  case(3)
                     fun_from_z_to_R = R_wCDM(redshift,H70,-1d0,O_3,0d0)    !=- LCDM(0.3,0.7,H_0)
                  case(4)
                     O_m = 1 + O_k - O_v !=-   the Friedmann equation   (bug)
                     fun_from_z_to_R = R_wCDM(redshift,H0,w,O_m,O_k)       !=- wCDM()
                  case default
                     if (RDR_error) write(*,*) ' do need choose the redshift-distance type '
               end select
			end function

      real(8) function mag_to_Lum( visible_magnitude , luminosity_distance )
         real(8) visible_magnitude,luminosity_distance ; Mag_to_Lum = 0d0
            if (visible_magnitude .ne. 0) Mag_to_Lum = L_sun * 1d1**( 0.4d0*( M_sun - &
               vis_to_abs_mag( visible_magnitude , luminosity_distance ) ) )
         end function

      real(8) function Lum_to_Mag(Luminosity)
         real(8) Luminosity ; Lum_to_Mag = 0d0
            Lum_to_Mag = M_sun - 2.5d0 * log10 ( Luminosity / L_sun )
         end function

      real(8) function abs_to_vis_mag( absolute_magnitude , luminosity_distance )
         real(8) absolute_magnitude,luminosity_distance ; abs_to_vis_mag = 0d0
            if (absolute_magnitude .ne. 0) abs_to_vis_mag = absolute_magnitude - 5d0 + 5d0 * log10 ( luminosity_distance * 1d5 ) !=- 10 pc
         end function

      real(8) function vis_to_abs_mag( visible_magnitude , luminosity_distance )
         real(8) visible_magnitude,luminosity_distance ; vis_to_abs_mag = 0d0
            if (visible_magnitude .ne. 0) vis_to_abs_mag = visible_magnitude + 5d0 - 5d0 * log10 ( luminosity_distance * 1d5 ) !=- 10 pc
         end function



			!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
			!real(8) function R_LCDM(z)			; real(8) z 		; R_LCDM  = dl * IntHwCDM(z,-1d0,O_m,0d0)	;	end function
   real(8) function R_LCDM(z,H0)			      ; real(8) z,H0 			; R_LCDM  = c/H0/1d5 * IntHwCDM(z,-1d0,O_3,0d0)	;	end function
   real(8) function R_PV(z,H0)					; real(8) z,H0 			; R_PV    = c/H0/1d5 * IntHwCDM(z,-1d0,0d0,0d0)	;	end function
   real(8) function R_EdeS(z,H0)				   ; real(8) z,H0 			; R_EdeS  = c/H0/1d5 * IntHwCDM(z,-1d0,1d0,0d0) ;	end function
   real(8) function R_wCDM(z,H0,w,Om,Ok)  	; real(8) z,H0,w,Om,Ok  ; R_wCDM  = c/H0/1d5 * RzwCDM(z,w,Om,Ok)	      ;	end function
   real(8) function R_CSS(z,H0)					; real(8) z,H0 			; R_CSS   = c/H0/1d5 * z								;	end function
   real(8) function R_TL(z,H0)					; real(8) z,H0 			; R_TL    = c/H0/1d5 * log(1d0+z)		 			;	end function
   real(8) function R_FCM(z,H0)					; real(8) z,H0 			; R_FCM   = c/H0/1d5 * YW(z) 						   ;	end function
   real(8) function R_MM(z,H0)					; real(8) z,H0 			; R_MM    = c/H0/1d5 * log(1d0+z) 					;	end function
			!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
   real(8) function D_LCDM(z,H0)					; real(8) z,H0		 		; D_LCDM  = (1d0+z)		* R_LCDM(z,H0) 			;	end function
   real(8) function D_PV(z,H0)					; real(8) z,H0		 		; D_PV    = (1d0+z)		* R_PV(z,H0) 				;	end function
   real(8) function D_EdeS(z,H0)					; real(8) z,H0		 		; D_EdeS	 = (1d0+z) 		* R_EdeS(z,H0) 			;	end function
   real(8) function D_wCDM(z,H0,w,Om,Ok)  	; real(8) z,H0,w,Om,Ok  ; D_wCDM	 = (1d0+z) 		* R_wCDM(z,H0,w,Om,Ok)  ;	end function
   real(8) function D_CSS(z,H0)					; real(8) z,H0	 		   ; D_CSS   = (1d0+z)		* R_CSS(z,H0) 				;	end function
   real(8) function D_TL(z,H0)					; real(8) z,H0				; D_TL  	 = dsqrt(1d0+z)* R_TL(z,H0) 				;	end function
   real(8) function D_FCM(z,H0)					; real(8) z,H0				; D_FCM 	 = (1d0+z) 		* R_FCM(z,H0) 				;	end function
   real(8) function D_MM(z,H0)					; real(8) z,H0				; D_MM  	 = (1d0+z) 		* R_MM(z,H0) 				;	end function
			!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
   real(8) function mu(dL)						   ; real(8) dL 		   	; mu		 = 25d0+5d0*log10(dL)					   ;	end function
			!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
   real(8) function dm_LCDM(z,H0)				; real(8) z,H0	 		   ; dm_LCDM = mu( D_LCDM(z,H0)			)		      ;	end function
   real(8) function dm_PV(z,H0)					; real(8) z,H0		 		; dm_PV	 = mu( D_PV(z,H0)				)		      ;	end function
   real(8) function dm_EdeS(z,H0)			   ; real(8) z,H0 			; dm_EdeS = mu( D_EdeS(z,H0)			)		      ;	end function
   real(8) function dm_wCDM(z,H0,w,Om,Ok)	   ; real(8) z,H0,w,Om,Ok  ; dm_wCDM = mu( D_wCDM(z,H0,w,Om,Ok))  	      ;	end function
   real(8) function dm_CSS(z,H0)					; real(8) z,H0		 		; dm_CSS	 = mu( D_CSS(z,H0)			)		      ;	end function
   real(8) function dm_TL(z,H0)					; real(8) z,H0		 		; dm_TL	 = mu( D_TL(z,H0)				)		      ;	end function
   real(8) function dm_FCM(z,H0)					; real(8) z,H0		 		; dm_FCM	 = mu( D_FCM(z,H0)			)		      ;	end function
   real(8) function dm_MM(z,H0)					; real(8) z,H0		 		; dm_MM	 = mu( D_MM(z,H0)				)		      ;	end function
			!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

			real(8) function YW(z)
				real(8) z
					YW=0d0 ; N=1d4 ; x_max=6d0

					if (z.ne.0) then

						if (z<5) x_max=3d0		!=- x_max / N =
						if (z<2) x_max=1.5d0		!=- x_max / N =
						if (z<0.5) x_max=0.5d0	!=- x_max / N =
						if (z<0.1) x_max=0.1d0	!=- x_max / N =

						do iii=N,1,-1
							x=x_max*iii/N	!	x[0:6] <=> z[0:20]
							y=4d0*x**0.5d0
							t=y/3.75d0

							if (-1d0 <= t .and. t <= 1d0) then
								I1 = y* (0.5d0 + 0.87890594*t**2d0 + 0.51498869*t**4d0 + 0.15084934*t**6d0 + 0.02658733*t**8d0 + &
									0.00301532*t**10d0 + 0.00032411*t**12d0)
								else
									I1 = dexp(y)/y**0.5d0*(0.39894228d0 - 0.03988024/t - 0.00362018/t**2d0 + 0.00163801/t**3d0 - &
										0.01031555/t**4d0 + 0.02282967/t**5d0 - 0.02895312/t**6d0 + 0.01787654/t**7d0 - 0.00420059/t**8d0)
								endif

							W=(2d0/y*I1)**0.5d0

							if ( dabs(z+1d0-W) <= (1d0+z)*x_max/N ) exit	!if (dabs(z+1d0-W)<1d-3) YW=W-1d0

							enddo
						YW=x ; !	write(*,*) z,x
						endif
				end function

			real(8) function IntHwCDM(z,w,Om,Ok)
            integer i,N
				real(8) z,w,Om,Ok,Ow
					E=1d0; N=1d4 ; IntHwCDM=0d0 ;  	!= N=1d4 +- 1Mpc =1

					if (z.ne.0) then
						dz = z/N ; Ow = 1d0 - Om + Ok;
						do i=1,N ; t = dz*(i-0.5d0) + 1d0
							!if (Om.ne.1d0)
							E=( Ow*t**( 3d0 + w*3d0 ) + Om*t*t*t - Ok*t*t )**0.5d0	!
							IntHwCDM = IntHwCDM + dz/E
							enddo
						endif
					end function IntHwCDM

			real(8) function RzwCDM(z,w,Om,Ok)
				real(8) z,w,Om,Ok

					RzwCDM=0d0 ; N=1d4	!=- N=1d4 +- 1Mpc

					if ( Ok == 0d0 ) then
						RzwCDM = IntHwCDM(z,w,Om,Ok);
						else
							if (Ok > 0d0) then
								RzwCDM = dsin( dabs(Ok)**0.5d0 * IntHwCDM(z,w,Om,Ok)  ) / dabs(Ok)**0.5d0
								else
									RzwCDM = dsinh( dabs(Ok)**0.5d0 * IntHwCDM(z,w,Om,Ok)  ) / dabs(Ok)**0.5d0
								end if
						end if
				end function RzwCDM
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
			subroutine models_head(titles)
				character(length) titles(N_titles)
	            titles(1)	='z	'										;
	            titles(2)	='z+1'										;
					titles(3)	='log(z+1)'									;

					titles(4)	='dLCDM/rLCDM'								;	titles(22)	='dLCDM'
					titles(5)	='d_CSS(z)/rLCDM'							;	titles(23)	='d_CSS(z)'
					titles(6)	='d_TL(z)/rLCDM'							;	titles(24)	='d_TL(z)'
					titles(7)	='d_FCM(z)/rLCDM'							;	titles(25)	='d_FCM(z)'
					titles(8)	='d_MM(z)/rLCDM'							;	titles(26)	='d_MM(z)'
					titles(9)	='d_wCDM(z,-1.0,0.3,-0.1)/rLCDM'		;	titles(27)	='d_wCDM(z,-1.0,0.3,-0.1)'
					titles(10)	='d_wCDM(z,-1.0,0.3, 0.1)/rLCDM'		;	titles(28)	='d_wCDM(z,-1.0,0.3, 0.1)'
					titles(11)	='d_wCDM(z,-1.0,0.1, 0.0)/rLCDM'		;	titles(29)	='d_wCDM(z,-1.0,0.1, 0.0)'
					titles(12)	='d_wCDM(z,-0.5,0.5, 0.2)/rLCDM'		;	titles(30)	='d_wCDM(z,-0.5,0.5, 0.2)'

					titles(13)	='rLCDM /rLCDM'							;	titles(31)	='rLCDM'
					titles(14)	='r_CSS(z)/rLCDM'							;	titles(32)	='r_CSS(z)'
					titles(15)	='r_TL(z)/rLCDM'							;	titles(33)	='r_TL(z)'
					titles(16)	='r_FCM(z)/rLCDM'							;	titles(34)	='r_FCM(z)'
					titles(17)	='r_MM(z)/rLCDM'							;	titles(35)	='r_MM(z)'
					titles(18)	='r_wCDM(z,-1.0,0.3,-0.1)/rLCDM'		;	titles(36)	='r_wCDM(z,-1.0,0.3,-0.1)'
					titles(19)	='r_wCDM(z,-1.0,0.3, 0.1)/rLCDM'		;	titles(37)	='r_wCDM(z,-1.0,0.3, 0.1)'
					titles(20)	='r_wCDM(z,-1.0,0.1, 0.0)/rLCDM'		;	titles(38)	='r_wCDM(z,-1.0,0.1, 0.0)'
					titles(21)	='r_wCDM(z,-0.5,0.5, 0.2)/rLCDM'		;	titles(39)	='r_wCDM(z,-0.5,0.5, 0.2)'

					titles(40)	='dm_LCDM-dm_LCDM'						;	titles(49)	='dm_LCDM'
					titles(41)	='dm_CSS(z)-dm_LCDM'						;	titles(50)	='dm_CSS(z)'
					titles(42)	='dm_TL(z)-dm_LCDM'						;	titles(51)	='dm_TL(z)'
					titles(43)	='dm_FCM(z)-dm_LCDM'						;	titles(52)	='dm_FCM(z)'
					titles(44)	='dm_MM(z)-dm_LCDM'						;	titles(53)	='dm_MM(z)'
					titles(45)	='dm_wCDM(z,-1.0,0.3,-0.1)-dm_LCDM'	;	titles(54)	='dm_wCDM(z,-1.0,0.3,-0.1)'
					titles(46)	='dm_wCDM(z,-1.0,0.3, 0.1)-dm_LCDM'	;	titles(55)	='dm_wCDM(z,-1.0,0.3, 0.1)'
					titles(47)	='dm_wCDM(z,-1.0,0.1, 0.0)-dm_LCDM'	;	titles(56)	='dm_wCDM(z,-1.0,0.1, 0.0)'
					titles(48)	='dm_wCDM(z,-0.5,0.5, 0.2)-dm_LCDM'	;	titles(57)	='dm_wCDM(z,-0.5,0.5, 0.2)'

				end subroutine

			subroutine models_calculating(datafile)
				character(length) datafile

					open(1,file=datafile, status='replace');	call models_head(titles)
						theformat='(A2,'//trim(inttostr(N_titles))//'(i2,1x,A22,1x))'
						do j=1,N_titles
                     write(1,'(A2,i2,1x,A40)') '# ',j,titles(j)
                     enddo
						write(1,theformat) '# ',(j,titles(j),j=1,N_titles)

						do i=1,N_cosmogrid
							z = dsqrt(z_max)/N_cosmogrid * i!(i+0.01d0)
							z = z*z

							LumDist(1)	= d_LCDM (z,H_0)					         ; MetrDist(1)	= r_LCDM (z,H_0)
							LumDist(2)	= d_CSS  (z,H_0)					         ; MetrDist(2)	= r_CSS  (z,H_0)
							LumDist(3)	= d_TL   (z,H_0)					         ; MetrDist(3)	= r_TL   (z,H_0)
							LumDist(4)	= d_FCM  (z,H_0)					         ; MetrDist(4)	= r_FCM  (z,H_0)
							LumDist(5)	= d_MM   (z,H_0)					         ; MetrDist(5)	= r_MM   (z,H_0)
							LumDist(6)	= d_wCDM (z,H_0,-1.0d0,0.85d0, 0.0d0)	; MetrDist(6)	= r_wCDM (z,H_0,-1.0d0,0.85d0, 0.0d0)
							LumDist(7)	= d_wCDM (z,H_0,-1.0d0,0.20d0, 0.0d0)	; MetrDist(7)	= r_wCDM (z,H_0,-1.0d0,0.20d0, 0.0d0)
							LumDist(8)	= d_wCDM (z,H_0,-1.0d0,0.1d0, 0.0d0)	; MetrDist(8)	= r_wCDM (z,H_0,-1.0d0,0.1d0, 0.0d0)
							LumDist(9)	= d_wCDM (z,H_0,-0.5d0,0.5d0, 0.2d0)	; MetrDist(9)	= r_wCDM (z,H_0,-0.5d0,0.5d0, 0.2d0)

							DistMod(1) = dm_LCDM (z,H_0)
							DistMod(2) = dm_CSS  (z,H_0)
							DistMod(3) = dm_TL   (z,H_0)
							DistMod(4) = dm_FCM  (z,H_0)
							DistMod(5) = dm_MM   (z,H_0)
							DistMod(6) = dm_wCDM (z,H_0,-1.0d0,0.85d0, 0.0d0)
							DistMod(7) = dm_wCDM (z,H_0,-1.0d0,0.20d0, 0.0d0)
							DistMod(8) = dm_wCDM (z,H_0,-1.0d0,0.1d0, 0.0d0)
							DistMod(9) = dm_wCDM (z,H_0,-0.5d0,0.5d0, 0.2d0)

							do j=1,N_models;	if ( LumDist(j)==0)  LumDist(j)=1d0
								if (MetrDist(j)==0) MetrDist(j)=1d0	;	end do

							theformat2='('//trim(inttostr(N_titles))//'(F20.8,6x))'
							write(1,theformat2) z,z+1,log(z+1)	, &
								( LumDist(j)/ LumDist(1),j=1,N_models)		, &
								(MetrDist(j)/MetrDist(1),j=1,N_models)		, &
								( LumDist(j)				,j=1,N_models)		, &
								(MetrDist(j)				,j=1,N_models)		, &
								( DistMod(j)- DistMod(1),j=1,N_models)		, &
								( DistMod(j)            ,j=1,N_models)
							end do

						close(1)

				end subroutine

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
			subroutine models(datafile)
				character(length) datafile

               call models_calculating(datafile)

 					!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
					GNUfields(1)  = 'set linestyle 1 lw 3 pt 7 ps 0.9 lt rgb "red"'
					GNUfields(2)  = 'set linestyle 2 lw 3 pt 7 ps 0.9 lt rgb "grey"'
					GNUfields(3)  = 'set linestyle 3 lw 3 pt 7 ps 0.9 lt rgb "black"'
					GNUfields(4)  = 'set linestyle 4 lw 3 pt 7 ps 0.9 lt rgb "green"'

					GNUfields(5)  = 'set linestyle 5 lw 3 pt 7 ps 0.9 lt rgb "brown"'
					GNUfields(22) = 'set linestyle 6 lw 3 pt 7 ps 0.9 lt rgb "blue"'
					GNUfields(23) = 'set linestyle 9 lw 3 pt 7 ps 0.9 lt rgb "brown"'

					GNUfields(10) = 'set yrange[-0.6:1.0]'
					GNUfields(9)  = 'set xlabel "z"' ; j = 1
					GNUfields(14) = 'set ylabel "{/Symbol D}log d_L(z) = log d_L(z)/d_L_{,{/Symbol L}-CDM}(z)"'
					GNUfields(7)  = 'set title "{/Symbol D}log d_L(z) for different models"'
					!GNUfields(24) = 'set key opaque'

					!GNUfields(31) = 'set key left bottom'
					GNUfields(31) = 'set nokey'

					do i=1,N_models
						GNUfields(31+i+i) = ', "" u '//trim(inttostr(j))//':(log10($'//trim(inttostr(i+3+N_models*0))// &
						'))  w l ls '//trim(inttostr(i))
						enddo

					GNUfields(34) = ' title "{{/Symbol L}-CDM}({/Symbol W}_{/Symbol L}=0.7)"'
					GNUfields(36) = ' title "{{/Symbol L}-CDM}({/Symbol W}_{/Symbol L}=1.0),CSS"'
					GNUfields(38) = ' title "{TL}"'
					GNUfields(40) = ' title "{FF}"'
					GNUfields(42) = ' title "{{/Symbol L}-CDM}({/Symbol W}_{/Symbol L}=0.9)"'
					GNUfields(44) = ' title "   w-CDM(w=-0.5,{/Symbol W}_{/Symbol L}=0.7,{/Symbol W}_k=0.2)"'

					GNUfields(20) = '#set logscale y'
					GNUfields(12) = '#'
					GNUfields(17) = '#'
					GNUfields(11) = 'set mxtics 2'
					GNUfields(16) = 'set mytics 2'

					graph_name = 'Shirokov_fig6' ; GNUfields(extention_out_figure)='#eps' ; call plot(datafile)
					graph_name = '_d_L(z)' ; GNUfields(extention_out_figure)='#png' ; call plot(datafile)
					!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
					GNUfields(14) = 'set ylabel "{/Symbol D}log r(z) = log r(z)/r_{{/Symbol L}-CDM}(z)"'
					GNUfields(7)  = 'set title "{/Symbol D}log r(z) for different models"'

					do i=1,N_models
						GNUfields(31+i+i) = ', "" u '//trim(inttostr(j))//':(log10($'//trim(inttostr(i+3+N_models*1))// &
						'))  w l ls '//trim(inttostr(i))
						enddo

					graph_name = 'Shirokov_fig4' ; GNUfields(extention_out_figure)='#eps' ; call plot(datafile)
					graph_name = '_d_r(z)' ; GNUfields(extention_out_figure)='#png' ; call plot(datafile)
					!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
					GNUfields(11) = 'set mxtics 5'
					GNUfields(16) = 'set mytics 5'

					GNUfields(9)  = 'set xlabel "z"' ; ii = 1
					!GNUfields(10) = 'set xrange[0.1:2]'
					!GNUfields(15) = 'set yrange[100:10000]'
					GNUfields(10) = 'set xrange[1:20]'
					GNUfields(15) = 'set yrange[500:1000000]'
					GNUfields(20) = 'set logscale'

					GNUfields(17) = 'set format y "10^{%L}"'
					GNUfields(14) = 'set ylabel "d_L(z), [Mpc]"'
					GNUfields(7)  = 'set title "d_L(z) for different models in logscales"'

               GNUfields(31) = 'set key right bottom'

					do i=1,N_models
						GNUfields(31+i+i) = ', "" u '//trim(inttostr(j))//':(($'//trim(inttostr(i+3+N_models*2))// &
						'))  w l ls '//trim(inttostr(i))
						enddo

					graph_name = 'Shirokov_fig5' ; GNUfields(extention_out_figure)='#eps' ; call plot(datafile)
					graph_name = '_d_L(z)-logscale' ; GNUfields(extention_out_figure)='#png' ; call plot(datafile)
					!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
					GNUfields(14) = 'set ylabel "r(z), [Mpc]"'
					GNUfields(7)  = 'set title "r(z) for different models in logscales"'
					GNUfields(15) = 'set yrange[500:100000]'

					do i=1,N_models
						GNUfields(31+i+i) = ', "" u '//trim(inttostr(j))//':(($'//trim(inttostr(i+3+N_models*3))// &
						'))  w l ls '//trim(inttostr(i))
						enddo

					graph_name = 'Shirokov_fig3' ; GNUfields(extention_out_figure)='#eps' ; call plot(datafile)
					graph_name = '_r(z)-logscale' ; GNUfields(extention_out_figure)='#png' ; call plot(datafile)
					!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
					GNUfields(14) = 'set ylabel "{/Symbol D}{/Symbol m}(z) = {/Symbol m}(z) - {/Symbol m}_{{/Symbol L}-CDM}(z)"'
					GNUfields(7)  = 'set title "{/Symbol D}{/Symbol m}(z) for different models"'
					GNUfields(15) = 'set yrange[-3:5]'
					GNUfields(10) = 'set xrange[0.001:20]'
					GNUfields(20) = 'unset logscale'
					GNUfields(17) = 'set format y "%g"'
					GNUfields(11) = 'set mxtics 2'
					GNUfields(16) = 'set mytics 2'

               GNUfields(31) = 'set nokey'

					do i=1,N_models
						GNUfields(31+i+i) = ', "" u '//trim(inttostr(j))//':(($'//trim(inttostr(i+3+N_models*4))// &
						'))  w l ls '//trim(inttostr(i))
						enddo

					GNUfields(34) = ' title "{/Symbol m}_{{/Symbol L}-CDM}({/Symbol W}_{/Symbol L}=0.7)"'
					GNUfields(36) = ' title "{/Symbol m}_{{/Symbol L}-CDM}({/Symbol W}_{/Symbol L}=1.0),CSS"'
					GNUfields(38) = ' title "{/Symbol m}_{TL}"'
					GNUfields(40) = ' title "{/Symbol m}_{FF}"'
					GNUfields(42) = ' title "{/Symbol m}_{{/Symbol L}-CDM}({/Symbol W}_{/Symbol L}=0.9)"'
					GNUfields(44) = ' title "   {/Symbol m}_{w-CDM}(w=-0.5,{/Symbol W}_{/Symbol L}=0.7,{/Symbol W}_k=0.2)"'

					graph_name = 'Shirokov_fig8' ; GNUfields(extention_out_figure)='#eps' ; call plot(datafile)
					graph_name = '_u(z)' ; GNUfields(extention_out_figure)='#png' ; call plot(datafile)
					!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
					GNUfields(14) = 'set ylabel "{/Symbol m}(z)"'
					GNUfields(7)  = 'set title "{/Symbol m}(z) for different models in z-logscale"'
					GNUfields(15) = 'set yrange[32:57]'
					GNUfields(10) = 'set xrange[0.1:20]'
					GNUfields(20) = 'set logscale x'
					GNUfields(11) = 'set mxtics 5'
					GNUfields(16) = 'set mytics 5'

					GNUfields(31) = 'set key right bottom'

					do i=1,N_models
						GNUfields(31+i+i) = ', "" u '//trim(inttostr(j))//':(($'//trim(inttostr(i+3+N_models*5))// &
						'))  w l ls '//trim(inttostr(i))
						enddo

					graph_name = 'Shirokov_fig7'  ; GNUfields(extention_out_figure)='#eps' ; call plot(datafile)
					graph_name = '_u(z)-logscale' ; GNUfields(extention_out_figure)='#png' ; call plot(datafile)
					!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


         end subroutine


		end module
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
