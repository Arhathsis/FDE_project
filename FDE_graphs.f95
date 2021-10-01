!=- Fractal Dimension Estimation
!=- © Stanislav Shirokov, 2014-2020

module FDE_graphs
	use GNUplot
   use global

	use FDE_config
	use FDE_paths
	use FDE_generators
	use FDE_methods

      logical :: plot_catalog_eps_flag = .false.

		integer, parameter :: N_graphs = 5

      integer :: method

		character(length) :: graths_names(N_graphs)='' , FDE_method_mask , method_name , FDE_geometry_mask , &
                        borders , y_axis_range

	contains



         subroutine plot_scaling(logscaling)
            logical logscaling

               call clear_plot_fields

					GNUfields(ls1)			=	'set linestyle 1 lw 3 pt 7 ps 0.2 lt rgb "blue"'
					GNUfields(ls2)			=	'set linestyle 2 lw 3 pt 7 ps 0.2 lt rgb "red"'
					GNUfields(ls3)			=	'set linestyle 3 lw 3 pt 7 ps 0.2 lt rgb "purple"'
					GNUfields(ls4)			=	'set linestyle 4 lw 3 pt 7 ps 0.2 lt rgb "orange"'

               GNUfields(legend)    =  'set key bottom right'
               GNUfields(grid)	   =	'set grid xtics ytics mxtics mytics'

               if (logscaling) then
                  GNUfields(format_y)	=	'set format y "10^{%L}"'
                  GNUfields(logscale)	=	'set logscale'
                  GNUfields(mxtics)		=	'10'
                  GNUfields(mytics)		=	'10'
                  else
                     GNUfields(logscale)	=	'unset logscale'
                     GNUfields(format_y)	=	'#set format y'
                     GNUfields(mxtics)		=	'5'
                     GNUfields(mytics)		=	'5'
                  end if

					GNUfields(format_x)	=	'#'
					GNUfields(yrange)    =  '[50:450]'

               GNUfields(title)     =  'The {/Symbol L}CDM model scaling'
					GNUfields(xlabel)		=	'z'
					GNUfields(ylabel)		=	'Mpc'

					GNUfields(plot1)		=	'"'//trim(slashfix(scaling_filepath))//'" u 1:2 w l ls 1'
               GNUfields(title1)    =  'r^{70}(z)'
					GNUfields(plot2)		=	'"'//trim(slashfix(scaling_filepath))//'" u 1:3 w l ls 2'
               GNUfields(title2)    =  'd_L^{70}(z)'
					GNUfields(plot3)		=	'"'//trim(slashfix(scaling_filepath))//'" u 1:4 w l ls 3'
               GNUfields(title3)    =  'r^{100}(z)'
					GNUfields(plot4)		=	'"'//trim(slashfix(scaling_filepath))//'" u 1:5 w l ls 4'
               GNUfields(title4)    =  'd_L^{100}(z)'

					graph_name = ''
                  GNUfields(extention_out_figure)='#png' ; call plot(scaling_filepath)

            end subroutine plot_scaling



         subroutine FDE_GNUplots

            do method=methods_count,1,-1
               FDE_method = method
               if (method/=1) call FDE_trend_log_xy( FDE_data_files_paths(method) )
               call FDE_plot(FDE_data_files_paths(method))
               enddo

            end subroutine FDE_GNUplots



         subroutine plot_means
            character(length) text
            integer i

            text = 'MEAN'
            do i=1,4
               FDE_method_names(i) = trim(text) // ' ' // trim(FDE_method_names(i))
               enddo

!                  ' \n the real parameters: N = ' // trim(adjustl(inttostr(FDE_points_amount))) // &
!                     ', D = '//trim(realtostrf(FDE_D(geometry),5,2)) // &
            FDE_points_amount = mean_FDE_points_amount

            call FDE_GNUplots
            call default_means

            end subroutine



      character(length) function FDE_method_mask_searching( data_file_path )
         character(length) data_file_path ; integer i

            l1=1 ; l2=1
            do i=1,length
               if ( data_file_path(i:i+3) == '.dat' ) l2=i-1
               if ( data_file_path(i:i)   == '_'    ) l1=1+i
               end do

            FDE_method_mask_searching = '' ; FDE_method_mask_searching = data_file_path(l1:l2)
         end function FDE_method_mask_searching



      subroutine FDE_plot(data_file_path)
         character(length) data_file_path

            FDE_method_mask = FDE_method_mask_searching( data_file_path )
            select case ( FDE_method_mask )

               case ('NP')
                  GNUfields(ylabel)		=	'N_{ops}, ops -- one-point-spheres'
                  method_name = FDE_method_names(1) !=- 'K-NEAREST NEIGHBORS ALGORITHM'
                  y_axis_range = '[ 0.1*' // trim(realtostrE(minval(FDE_out_NP(2,:) , FDE_out_NP(2,:) .ne. 0))) // &
                     ' :10*' // trim(realtostrE(maxval(FDE_out_NP(2,:) , FDE_out_NP(2,:) .ne. 0))) // ' ]'
               case ('MD')
                  GNUfields(ylabel)		=	'N_{segments}'
                  method_name = FDE_method_names(2) !=- 'MUTUAL DISTANCES ALGORITHM'
                  y_axis_range = '[ 0.1*' // trim(realtostrE(minval(FDE_out_MD(2,:) , FDE_out_MD(2,:) .ne. 0))) // &
                     ' :10*' // trim(realtostrE(maxval(FDE_out_MD(2,:) , FDE_out_MD(2,:) .ne. 0))) // ' ]'
               case ('IntCD')
                  GNUfields(ylabel)		=	'CD_{integral}'
                  method_name = FDE_method_names(3) !=- 'INTEGRAL CONDITIONAL DENSITY ALGORITHM'
                  y_axis_range = '[ 0.1*' // trim(realtostrE(minval(FDE_out_intCD(2,:) , FDE_out_intCD(2,:) .ne. 0))) // &
                     ' :10*' // trim(realtostrE(maxval(FDE_out_intCD(2,:) , FDE_out_intCD(2,:) .ne. 0))) // ' ]'
               case ('DiffCD')
                  GNUfields(ylabel)		=	'CD_{differential}'
                  method_name = FDE_method_names(4) !=- 'DIFFERENCIAL CONDITIONAL DENSITY ALGORITHM'
                  y_axis_range = '[ 0.1*' // trim(realtostrE(minval(FDE_out_diffCD(2,:) , FDE_out_diffCD(2,:) .ne. 0))) // &
                     ' :10*' // trim(realtostrE(maxval(FDE_out_diffCD(2,:) , FDE_out_diffCD(2,:) .ne. 0))) // ' ]'
               case default
                  write(*,*) 'FDE_plot: unknown method mask'
                  write(*,*) trim(FDE_method_mask)
                  FDE_method_mask = 'FDE plot: unknown method mask'
               end select

               do geometry=1,geometries_count

                  select case (geometry)
                     case(1)
                        FDE_geometry_mask = 'all sphere'
                     case(2)
                        FDE_geometry_mask = 'north hemisphere'
                     case(3)
                        FDE_geometry_mask = 'south hemisphere'
                     case(4)
                        FDE_geometry_mask = 'both north and south hemispheres'
                     case default
                        write(*,*) 'FDE_plot: unknown geometry mask'
                        FDE_geometry_mask = 'FDE plot: unknown geometry mask'
                     end select

                  call FDE_plot_script(data_file_path)
                  end do

         end subroutine



      subroutine FDE_plot_script(data_file_path)
         character(length) data_file_path

               call clear_plot_fields

					GNUfields(ls1)			=	'set linestyle 1 lw 1 pt 7 ps 0.2 lt rgb "blue"'
					GNUfields(ls2)			=	'set linestyle 2 lw 1 pt 7 ps 0.2 lt rgb "red"'
					GNUfields(ls3)			=	'set linestyle 3 lw 1 pt 7 ps 0.2 lt rgb "black"'

					GNUfields(logscale)	=	'set logscale'
               GNUfields(legend)    =  'set key bottom right'
               GNUfields(grid)	   =	'set grid xtics ytics mxtics mytics'

					GNUfields(format_y)	=	'set format y "10^{%L}"'
					GNUfields(format_x)	=	'#'
					GNUfields(mxtics)		=	'10'
					GNUfields(mytics)		=	'10'

               GNUfields(add3) = 'set trange ' // trim(y_axis_range)
               GNUfields(add4) = 'set yrange ' // trim(y_axis_range)
               GNUfields(add5) = 'set parametric'

            if ( FDE_method==3 .or. FDE_method==4 ) then
               FDE_D(:) = 3d0+FDE_LS_coefficients(:,1)
               else
                  FDE_D(:) = FDE_LS_coefficients(:,1)
               endif

            if ( FDE_method<3 ) then !=- case
               lb=FDE_left_border_MD
               rb=FDE_right_border_MD
               else
                  lb=FDE_left_border
                  rb=FDE_right_border
            end if
               theformat = '(F6.1)'
                  borders = '[' // trim(realtostrFF(lb,theformat)) // &
                     ' :' // trim(realtostrff(rb,theformat)) // ' ]'
					GNUfields(title) 		=	'set title "' // trim(method_name) // &
                  ' \n the geometry: ' // trim(FDE_geometry_mask) // &
                  ' \n the effective parameters: N = ' // trim(adjustl(inttostr(N_points))) // &
                     ', D = '//trim(realtostrf(Fractal_Dimension,5,2)) // &
                     ', R_{max} = ' // trim(inttostr(int(radius_limit))) // ' Mpc' // &
                  ' \n the real parameters: N = ' // trim(adjustl(inttostr(FDE_points_amount))) // &
                     ', D = '//trim(realtostrf(FDE_D(geometry),5,2)) // &
                     ', R_{nn} = ' // trim(realtostrFF(R_nn(geometry),theformat)) // ' Mpc \' !// &
               GNUfields(1+title)	=	'\n the approximation borders: ' // trim(borders) // '"'  !=-  noenhanced
					GNUfields(xlabel)		=	'R, [Mpc]'

					GNUfields(plot1)		=	'"'//trim(slashfix(data_file_path))// &
                  '" u 1:' // trim(inttostr(2*geometry)) // ':' // trim(inttostr(1+2*geometry)) // ' w yerrorb ls 1'
                  GNUfields(title1) = 'N_{ops}'
					GNUfields(plot2)		=	', ' // trim(realtostr(R_nn(geometry))) // ',t ls 2'
                  GNUfields(title2) = 'R_{nn}'
               if (method/=1) GNUfields(plot3)		=	', "" u 1:($1**' // &   !=-    ' // trim(borders) // '
                  trim(realtostr(FDE_LS_coefficients(geometry,1))) // &
                  ' * 10**' // trim(realtostr(FDE_LS_coefficients(geometry,2))) // ') w l ls 3'
               if (method/=1) GNUfields(title3)    = 'trend approximation'
					if (method/=1) GNUfields(plot4)		=	', ' // trim(realtostr(lb))  // ',t ls 3'
               if (method/=1) GNUfields(title4) = ' notitle'
               if (method/=1) GNUfields(plot5)		=	', ' // trim(realtostr(rb)) // ',t ls 3'
               if (method/=1) GNUfields(title5) = ' notitle'

					graph_name = FDE_geometry_mask(1:1)
                  GNUfields(extention_out_figure)='#png' ; call plot(data_file_path)

         end subroutine



		subroutine plot_catalog(catalog_path)
			character(length) catalog_path

                  write(*,*) '   plot_catalog: plotting'

               call clear_plot_fields
               call Generate_Cyrle(radius_limit)

					GNUfields(ls1)			=	'set linestyle 1 lw 1 pt 7 ps 0.2 lt rgb "blue"'
					GNUfields(ls2)			=	'set linestyle 2 lw 1 pt 7 ps 0.2 lt rgb "black"'

					GNUfields(logscale)	=	'#set logscale'
					GNUfields(extention_out_figure) = '#png'
               GNUfields(grid)	   =	'set grid xtics ytics mxtics mytics'

					GNUfields(format_y)	=	'#'
					GNUfields(format_x)	=	'#'
					GNUfields(mxtics)		=	'1'
					GNUfields(mytics)		=	'1'

               GNUfields(title)   =  'set title "catalog: ' // trim(name(catalog_path)) // ' \'
					GNUfields(title+1) 		=	'\n X vs Y plan \'

					GNUfields(title+2) 	=	'\\n N = '//trim(inttostrf(final_counter,10))// &
                  ', D = '//trim(realtostrf(Fractal_Dimension,3,1))//' " noenhanced'
					GNUfields(xlabel)		=	'Y, [Mpc]'
					GNUfields(ylabel)		=	'X, [Mpc]'
					GNUfields(xrange)		=	'['//trim(realtostr(-radius_limit))//':'// &
                  trim(realtostr(radius_limit))//']'
					GNUfields(yrange)		=	'['//trim(realtostr(-radius_limit))//':'// &
                  trim(realtostr(radius_limit))//']'

					GNUfields(plot1)		=	'"'//trim(slashfix(catalog_path))//'" u '// &
                  trim(inttostr( N_col_cat_y ))//':'//trim(inttostr( N_col_cat_x ))//' w p ls 1'
               GNUfields(title1)		=	'set points'
					GNUfields(plot2)		=	', "'//trim(slashfix(cyrle_path))//'" u 1:2 w l ls 2'
               GNUfields(title2)		=	'set border'

					if (plot_catalog_eps_flag) then
                  graph_name = 'XY_plan' ; GNUfields(extention_out_figure)='#eps' ; call plot(catalog_path)
                  endif
					graph_name = 'XY_plan' ; GNUfields(extention_out_figure)='#png' ; call plot(catalog_path)



					GNUfields(title+1) 		=	'\\n X vs Z plan \'
					GNUfields(xlabel)		=	'Z, [Mpc]'
					GNUfields(ylabel)		=	'X, [Mpc]'


					GNUfields(plot1)		=	'"'//trim(slashfix(catalog_path))//'" u '// &
                  trim(inttostr( N_col_cat_z	))//':'//trim(inttostr( N_col_cat_x ))//' w p ls 1'

					if (plot_catalog_eps_flag) then
                  graph_name = 'XZ_plan' ; GNUfields(extention_out_figure)='#eps' ; call plot(catalog_path)
                  endif
					graph_name = 'XZ_plan' ; GNUfields(extention_out_figure)='#png' ; call plot(catalog_path)



					GNUfields(title+1)		=	'\\n Y vs Z plan \'
					GNUfields(xlabel)		=	'Z, [Mpc]'
					GNUfields(ylabel)		=	'Y, [Mpc]'

					GNUfields(plot1)		=	'"'//trim(slashfix(catalog_path))//'" u '// &
                  trim(inttostr( N_col_cat_z	))//':'//trim(inttostr( N_col_cat_y ))//' w p ls 1'

					if (plot_catalog_eps_flag) then
                  graph_name = 'YZ_plan' ; GNUfields(extention_out_figure)='#eps' ; call plot(catalog_path)
                  endif
					graph_name = 'YZ_plan' ; GNUfields(extention_out_figure)='#png' ; call plot(catalog_path)

               GNUfields(plot2:title2) = ''


					GNUfields(xrange)		=	'[*:*]'
					GNUfields(yrange)		=	'[*:*]'
					GNUfields(title+1)		=	'\\n m vs R plan \'
					GNUfields(xlabel)		=	'R, [Mpc]'
					GNUfields(ylabel)		=	'm'

					GNUfields(plot1)		=	'"'//trim(slashfix(catalog_path))//'" u '// &
                  trim(inttostr( N_col_cat_dl	))//':'//trim(inttostr( N_col_cat_mag ))//' w p ls 1'

					if (plot_catalog_eps_flag) then
                  graph_name = 'mag_vs_R' ; GNUfields(extention_out_figure)='#eps' ; call plot(catalog_path)
                  endif
					graph_name = 'mag_vs_R' ; GNUfields(extention_out_figure)='#png' ; call plot(catalog_path)



					GNUfields(title+1)		=	'\\n M vs R plan \'
					GNUfields(xlabel)		=	'R, [Mpc]'
					GNUfields(ylabel)		=	'M'

					GNUfields(plot1)		=	'"'//trim(slashfix(catalog_path))//'" u '// &
                  trim(inttostr( N_col_cat_dl	))//':'//trim(inttostr( N_col_cat_M ))//' w p ls 1'

					if (plot_catalog_eps_flag) then
                  graph_name = 'M_vs_R' ; GNUfields(extention_out_figure)='#eps' ; call plot(catalog_path)
                  endif
               graph_name = 'M_vs_R' ; GNUfields(extention_out_figure)='#png' ; call plot(catalog_path)



					GNUfields(title+1)		=	'\\n R vs z plan \'
					GNUfields(xlabel)		=	'z'
					GNUfields(ylabel)		=	'R, [Mpc]'

					GNUfields(plot1)		=	'"'//trim(slashfix(catalog_path))//'" u '// &
                  trim(inttostr( N_col_cat_rs	))//':'//trim(inttostr( N_col_cat_dl))//' w p ls 1'

					if (plot_catalog_eps_flag) then
                  graph_name = 'R_vs_z' ; GNUfields(extention_out_figure)='#eps' ; call plot(catalog_path)
                  endif
					graph_name = 'R_vs_z' ; GNUfields(extention_out_figure)='#png' ; call plot(catalog_path)



					GNUfields(title+1)		=	'\\n L vs R plan \'
					GNUfields(xlabel)		=	'R, [Mpc]'
					GNUfields(ylabel)		=	'L, [erg]'
					GNUfields(logscale)  =  '#set logscale y'
					GNUfields(format_y)  =  'set format y "10^{%L}"'

               GNUfields(yrange)    =  '[*:*]'
               !GNUfields(yrange)		=	'[10**'//trim(realtostr( log10( minimal_luminosity )	))//':10**'// &
               !   trim(realtostr( log10( maximal_luminosity )	))//']'

					GNUfields(plot1)		=	'"'//trim(slashfix(catalog_path))//'" u '// &
                  trim(inttostr( N_col_cat_dl	))//':'//trim(inttostr( N_col_cat_Lum ))//' w p ls 1'

					if (plot_catalog_eps_flag) then
                  graph_name = 'L_vs_R' ; GNUfields(extention_out_figure)='#eps' ; call plot(catalog_path)
                  endif
					graph_name = 'L_vs_R' ; GNUfields(extention_out_figure)='#png' ; call plot(catalog_path)

         end subroutine

end module
