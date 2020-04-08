!=- fortran-libraries
!=- © Stanislav Shirokov, 2014-2020

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
	module cosmology
		use global
		use GNUplot

      logical ::  RDR_calculating = .true.   , &   !=- .false.
                  RDR_error       = .true.

		integer,parameter :: N_models = 9   , &
                  N_RDR_grid        = 2d3 , &   !=- graph grid
                  N_RDR_grid_mode   = 3         !=- progression power

		integer ::  N_grid            = 2d2 , &

                  RDR_z_max         = 2d1 , &
                  RDR_type          = 1      !=- 0 is LCDM(0.3,70) , 1 is LCDM(\O_m,H_0) , 2 is wCDM(w,\O_m,H_0,O_k)


		real(8) ::	c   = 2.99792458d10	, &	!=- cm/s
						H_0 = 70d0         , &	!=- km/s/Mpc
						H70 = 7d1            , &
						dl  = 4421.71767d0   , &	!=- c/H_0/1d5 in GCS

					 L_sun = 3.827d33       , &   !=- erg/s
					 M_sun = 4.83d0         , &

						q_0 = 0.5d0				, &
						O_m = 0.308d0 			, &
						O_w = 1d0				, &
						O_k = 0d0            , &
					 z_max = 2d1				, &
					 O_3	 = 0.3d0				, &
					 O_7	 = 0.7d0


		real(8) :: LumDist(N_models) = 1d0, MetrDist(N_models)= 1d0, DistMod(N_models)= 1d0, & !=- z for FCM > 0.0025
			dz, x_max, I1, RDR(2,N_RDR_grid)

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
		contains

      subroutine redshift_distance_relation  !=- RDR( redshift , luminosity distance )
         if (RDR_calculating) then
            do i=1,N_RDR_grid
               z = z_max**(1d0/N_RDR_grid_mode)/N_grid * i ; z = z**N_RDR_grid_mode
               RDR(1,i) = z
               RDR(2,i) = fun_from_z_to_R(z)
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

		real(8) function fun_from_z_to_R(redshift)
			real(8) redshift ; fun_from_z_to_R = 0d0
            select case (RDR_type)
                  case(0)
                     fun_from_z_to_R = D_wCDM(redshift,-1d0,O_3,0d0)    !=- LCDM(0.3,0.7,H_0)
                  case(1)
                     fun_from_z_to_R = D_wCDM(redshift,-1d0,O_m,0d0)    !=- LCDM(\O_m,H_0)
                  case(2)
                     fun_from_z_to_R = D_wCDM(redshift,w,O_m,O_k)       !=- wCDM()
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
			real(8) function R_LCDM(z)					; real(8) z 			; R_LCDM  = c/H_0/1d5 * IntHwCDM(z,-1d0,O_3,0d0)		;	end function
			real(8) function R_PV(z)					; real(8) z 			; R_PV    = c/H_0/1d5 * IntHwCDM(z,-1d0,0d0,0d0)		;	end function
			real(8) function R_EdeS(z)					; real(8) z 			; R_EdeS  = c/H_0/1d5 * IntHwCDM(z,-1d0,1d0,0d0) 	;	end function
			real(8) function R_wCDM(z,w,O_m,O_k)	; real(8) z,w,O_m,O_k; R_wCDM  = c/H_0/1d5 * RzwCDM(z,w,O_m,O_k)			;	end function
			real(8) function R_CSS(z)					; real(8) z 			; R_CSS   = c/H_0/1d5 * z									;	end function
			real(8) function R_TL(z)					; real(8) z 			; R_TL    = c/H_0/1d5 * log(1d0+z)		 				;	end function
			real(8) function R_FCM(z)					; real(8) z 			; R_FCM   = c/H_0/1d5 * YW(z) 								;	end function
			real(8) function R_MM(z)					; real(8) z 			; R_MM    = c/H_0/1d5 * log(1d0+z) 						;	end function
			!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
			real(8) function D_LCDM(z)					; real(8) z		 		; D_LCDM  = (1d0+z)		* R_LCDM(z) 			;	end function
			real(8) function D_PV(z)					; real(8) z		 		; D_PV    = (1d0+z)		* R_PV(z) 				;	end function
			real(8) function D_EdeS(z)					; real(8) z		 		; D_EdeS	 = (1d0+z) 		* R_EdeS(z) 			;	end function
			real(8) function D_wCDM(z,w,O_m,O_k)	; real(8) z,w,O_m,O_k; D_wCDM	 = (1d0+z) 		* R_wCDM(z,w,O_m,O_k);	end function
			real(8) function D_CSS(z)					; real(8) z		 		; D_CSS   = (1d0+z)		* R_CSS(z) 				;	end function
			real(8) function D_TL(z)					; real(8) z				; D_TL  	 = dsqrt(1d0+z)* R_TL(z) 				;	end function
			real(8) function D_FCM(z)					; real(8) z				; D_FCM 	 = (1d0+z) 		* R_FCM(z) 				;	end function
			real(8) function D_MM(z)					; real(8) z				; D_MM  	 = (1d0+z) 		* R_MM(z) 				;	end function
			!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
			real(8) function mu(dL)						; real(8) dL 			; mu		 = 25d0+5d0*log10(dL)					;	end function
			!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
			real(8) function dm_LCDM(z)				; real(8) z		 		; dm_LCDM = mu( D_LCDM(z)			)				;	end function
			real(8) function dm_PV(z)					; real(8) z		 		; dm_PV	 = mu( D_PV(z)				)				;	end function
			real(8) function dm_EdeS(z)				; real(8) z 			; dm_EdeS = mu( D_EdeS(z)			)				;	end function
			real(8) function dm_wCDM(z,w,O_m,O_k)	; real(8) z,w,O_m,O_k; dm_wCDM = mu( D_wCDM(z,w,O_m,O_k)	)			;	end function
			real(8) function dm_CSS(z)					; real(8) z		 		; dm_CSS	 = mu( D_CSS(z)			)				;	end function
			real(8) function dm_TL(z)					; real(8) z		 		; dm_TL	 = mu( D_TL(z)				)				;	end function
			real(8) function dm_FCM(z)					; real(8) z		 		; dm_FCM	 = mu( D_FCM(z)			)				;	end function
			real(8) function dm_MM(z)					; real(8) z		 		; dm_MM	 = mu( D_MM(z)				)				;	end function
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

			real(8) function IntHwCDM(z,w,O_m,O_k)
				real(8) z,w,O_m,O_k
					E=1d0; N=1d4 ; IntHwCDM=0d0 ;  	!= N=1d4 +- 1Mpc =1

					if (z.ne.0) then
						dz = z/N ; O_w = 1d0 - O_m + O_k;
						do iii=1,N ; t = dz*(iii-0.5d0) + 1d0
							if (O_m.ne.1d0) E=( O_w*t**( 3d0 + w*3d0 ) + O_m*t*t*t - O_k*t*t )**0.5d0	!
							IntHwCDM = IntHwCDM + dz/E
							enddo
						endif
					end function IntHwCDM

			real(8) function RzwCDM(z,w,O_m,O_k)
				real(8) z,w,O_m,O_k
					RzwCDM=0d0 ; N=1d4	!=- N=1d4 +- 1Mpc

					if ( O_k == 0d0 ) then
						RzwCDM = IntHwCDM(z,w,O_m,O_k);
						else
							if (O_k > 0d0) then
								RzwCDM = dsin( dabs(O_k)**0.5d0 * IntHwCDM(z,w,O_m,O_k)  ) / dabs(O_k)**0.5d0
								else
									RzwCDM = dsinh( dabs(O_k)**0.5d0 * IntHwCDM(z,w,O_m,O_k)  ) / dabs(O_k)**0.5d0
								end if
						end if
				end function RzwCDM
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
			subroutine models_head(titles)
				character(len) titles(N_titles)
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
				character(len) datafile

					open(1,file=datafile, status='replace');	call models_head(titles)
						theformat='(A2,'//trim(inttostr(N_titles))//'(i2,1x,A22,1x))'
						do j=1,N_titles
                     write(1,'(A2,i2,1x,A40)') '# ',j,titles(j)
                     enddo
						write(1,theformat) '# ',(j,titles(j),j=1,N_titles)

						do i=1,N_grid
							z = dsqrt(z_max)/N_grid * i!(i+0.01d0)
							z = z*z

							LumDist(1)	= d_LCDM(z)					         ; MetrDist(1)	= r_LCDM(z)
							LumDist(2)	= d_CSS(z)					         ; MetrDist(2)	= r_CSS(z)
							LumDist(3)	= d_TL(z)					         ; MetrDist(3)	= r_TL(z)
							LumDist(4)	= d_FCM(z)					         ; MetrDist(4)	= r_FCM(z)
							LumDist(5)	= d_MM(z)					         ; MetrDist(5)	= r_MM(z)
							LumDist(6)	= d_wCDM(z,-1.0d0,0.85d0, 0.0d0)	; MetrDist(6)	= r_wCDM(z,-1.0d0,0.85d0, 0.0d0)
							LumDist(7)	= d_wCDM(z,-1.0d0,0.20d0, 0.0d0)	; MetrDist(7)	= r_wCDM(z,-1.0d0,0.20d0, 0.0d0)
							LumDist(8)	= d_wCDM(z,-1.0d0,0.1d0, 0.0d0)	; MetrDist(8)	= r_wCDM(z,-1.0d0,0.1d0, 0.0d0)
							LumDist(9)	= d_wCDM(z,-0.5d0,0.5d0, 0.2d0)	; MetrDist(9)	= r_wCDM(z,-0.5d0,0.5d0, 0.2d0)

							DistMod(1) = dm_LCDM(z)
							DistMod(2) = dm_CSS(z)
							DistMod(3) = dm_TL(z)
							DistMod(4) = dm_FCM(z)
							DistMod(5) = dm_MM(z)
							DistMod(6) = dm_wCDM(z,-1.0d0,0.85d0, 0.0d0)
							DistMod(7) = dm_wCDM(z,-1.0d0,0.20d0, 0.0d0)
							DistMod(8) = dm_wCDM(z,-1.0d0,0.1d0, 0.0d0)
							DistMod(9) = dm_wCDM(z,-0.5d0,0.5d0, 0.2d0)

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
				character(len) datafile

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
