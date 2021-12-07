!=- Fractal Dimension Estimation
!=- © Stanislav Shirokov, 2014-2022

module FDE_scripts
	use math
	use GNUplot
	use cosmology

	use FDE_config
	use FDE_catalogs
	use FDE_paths
   use FDE_methods
	use FDE_graphs
   use FDE_generators
   use FDE_tables

	contains

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine gen_FDE_matrix_script
            integer D_amount, N_left, N_right, N_amount, ls_left, ls_right, ls_amount
            real(8) D_left, D_right

      !      D_left   = 2d0 ; D_right    = 2d0   ; D_amount   = 1
      !      N_left   = 5d3 ; N_right    = 15d3  ; N_amount   = 5
      !      ls_left  = 5   ; ls_right   = 20    ; ls_amount  = 4

            D_left   = 1.9d0  ; D_right    = 2.1d0    ; D_amount   = 2
            N_left   = 5d3    ; N_right    = 15d3     ; N_amount   = 1
            ls_left  = 5      ; ls_right   = 20       ; ls_amount  = 1

            call create_table_of_means(D_left, D_right, D_amount, N_left, N_right, N_amount, ls_left, ls_right, ls_amount)

            end

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine FDE_matrix_preset

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
               !=- all dimensions depend on the multiness level

               MC_catalog_prefix = 'Matrix_MCantor' !=- the prefix text for new files

               MC_uniform_approach_parameter( 1 ) = 1 !=- uniformity on V_0 level (on the multiness level) ; default = 0

               MC_generations( 1 )  = 3   !=- number of hierarchical levels (on the multiness level) ; 1 < g < 10
               MC_generations( 2 )  = 6   !=- number of hierarchical levels (on the multiness level) ; 1 < g < 10

               FDE_sequence_length  = 0   !=- length of series of sets , default 0 (changes automatically on 1) , old set_number

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
               !=- uniform selection of set points
               !=- there is two methods of uniform selection of points, changed by FDE_selected_generators in FDE_config
               !=- default .true. is uniform selection at each generation on each multiness level
               !=- .false.  is uniform selection only at final generation on each multiness level

               FDE_probability      = 1 !=- selection probability ; for selection disabling, set to 1

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
               !=- default parameters changing

                  MC_target_radius         = 1d2   !=- default is 1d2
                  MC_target_point_number   = 0     !=- default is 1d4

                  MC_Fractal_Dimension          ( : ) = 0d0 !=- default is 2.0

                  MC_extensionality             ( : ) = 0   !=- default is 2

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

               FDE_approximation_ranges( 1, FDE_method_MD ) = 1d0   !=- 1d1 * GCC_radius_limit / 2**(GCC_generations+GCC_Super_generations-2)
               FDE_approximation_ranges( 2, FDE_method_MD ) = 0.8d0 * MC_target_radius   !=- 1d2 * GCC_radius_limit / 2**(GCC_generations+GCC_Super_generations-2)

               FDE_approximation_ranges( 1, FDE_method_intCD ) = 6d0
               FDE_approximation_ranges( 2, FDE_method_intCD ) = 3d1

               FDE_approximation_ranges( 1, FDE_method_diffCD ) = 5d0
               FDE_approximation_ranges( 2, FDE_method_diffCD ) = 4d1

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

            end

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine create_table_of_means(D_left, D_right, D_amount, N_left, N_right, N_amount, ls_left, ls_right, ls_amount)

            integer geometry,iloop,jloop,kloop

            integer, intent(in) :: D_amount, N_left, N_right, N_amount, ls_left, ls_right, ls_amount
            real(8), intent(in) :: D_left, D_right
            integer, dimension(N_amount)  :: N_grid
            real(8), dimension(D_amount)  :: D_grid
            integer, dimension(ls_amount) :: ls_grid
            real(8), dimension(D_amount,N_amount,ls_amount,geometries_count) :: table

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

               call create_linear_grid       (  D_left   ,  D_right  ,  D_amount ,  D_grid   )
               call create_logscale_grid_int (  N_left   ,  N_right  ,  N_amount ,  N_grid   )
               call create_logscale_grid_int (  ls_left  ,  ls_right ,  ls_amount,  ls_grid  )

               do iloop = 1,D_amount
                  do jloop = 1,N_amount
                     do kloop = 1,ls_amount
                        call FDE_testing_means_script(D_grid(iloop), N_grid(jloop), ls_grid(kloop))
                        table(iloop,jloop,kloop,:) = mean_FDE_D(2,:)
                     end do
                  end do
               end do

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

               do geometry=1,geometries_count

                  do iloop = 1,D_amount

                     unit_1 = random_unit()

                     open(unit_1, file=trim(folders( folder_Statistics_Report_tables )) // 'D_table_' // &
                        adjustl(trim(inttostr(N_right))) // '_' // &
                        adjustl(trim(inttostr(ls_right))) // '_' // &
                        adjustl(trim(realtostrf(D_grid(iloop),4,2))) // '_' // &
                        trim(FDE_geometry_mask(geometry)(1:1)) // '.dat',status='replace')

                     write(unit_1,*) 'ls|N', N_grid
                     do jloop = 1,ls_amount
                        write(unit_1,*) ls_grid(jloop), table(iloop,:,jloop,geometry)
                        end do

                        close(unit_1)
                     end do
                  enddo

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

         end subroutine create_table_of_means

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine FDE_testing_means_script(D_arg, N_arg, ls_arg)
            real(8), intent(in) :: D_arg
            integer, intent(in) :: N_arg, ls_arg

            call FDE_matrix_preset

            FDE_sequence_length     = ls_arg
            MC_Fractal_Dimension(:) = D_arg
            MC_target_point_number  = N_arg

            call creating_MC_means
            call default_means

         end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine FDE_default_preset

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
               !=- all dimensions depend on the multiness level

               MC_catalog_prefix = 'MCantor' !=- the prefix text for new files

               MC_uniform_approach_parameter( 1 ) = 1 !=- uniformity on V_0 level (on the multiness level) ; default = 0

               MC_generations( 1 )  = 3   !=- number of hierarchical levels (on the multiness level) ; 1 < g < 10
               MC_generations( 2 )  = 5   !=- number of hierarchical levels (on the multiness level) ; 1 < g < 10

               FDE_sequence_length  = 1   !=- length of series of sets

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
               !=- uniform selection of set points
               !=- there is two methods of uniform selection of points, changed by FDE_selected_generators in FDE_config
               !=- default .true. is uniform selection at each generation on each multiness level
               !=- .false.  is uniform selection only at final generation on each multiness level

               FDE_probability      = 0.5 !=- selection probability ; for selection disabling, set to 1

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
               !=- default parameters changing

                  MC_target_radius         = 1d2   !=- default is 1d2
                  MC_target_point_number   = 0     !=- default is 1d4

                  MC_Fractal_Dimension          ( : ) = 0d0 !=- default is 2.0

                  MC_extensionality             ( : ) = 0   !=- default is 2

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

               FDE_approximation_ranges( 1, FDE_method_MD ) = 1d0   !=- 1d1 * GCC_radius_limit / 2**(GCC_generations+GCC_Super_generations-2)
               FDE_approximation_ranges( 2, FDE_method_MD ) = 0.8d0 * MC_target_radius   !=- 1d2 * GCC_radius_limit / 2**(GCC_generations+GCC_Super_generations-2)

               FDE_approximation_ranges( 1, FDE_method_intCD ) = 6d0
               FDE_approximation_ranges( 2, FDE_method_intCD ) = 3d1

               FDE_approximation_ranges( 1, FDE_method_diffCD ) = 5d0
               FDE_approximation_ranges( 2, FDE_method_diffCD ) = 4d1

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

            end

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine FDE_testing_preset

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
               !=- all dimensions depend on the multiness level

               MC_catalog_prefix = 'MCantor' !=- the prefix text for new files

               MC_uniform_approach_parameter( 1 ) = 1 !=- uniformity on V_0 level (on the multiness level) ; default = 0

               MC_generations( 1 )  = 3   !=- number of hierarchical levels (on the multiness level) ; 1 < g < 10
               MC_generations( 2 )  = 5   !=- number of hierarchical levels (on the multiness level) ; 1 < g < 10

               FDE_sequence_length  = 2   !=- length of series of sets

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
               !=- uniform selection of set points
               !=- there is two methods of uniform selection of points, changed by FDE_selected_generators in FDE_config
               !=- default .true. is uniform selection at each generation on each multiness level
               !=- .false.  is uniform selection only at final generation on each multiness level

               FDE_probability      = 0.5 !=- selection probability ; for selection disabling, set to 1

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
               !=- default parameters changing

                  MC_target_radius         = 1d2   !=- default is 1d2
                  MC_target_point_number   = 0     !=- default is 1d4

                  MC_Fractal_Dimension          ( : ) = 0d0 !=- default is 2.0

                  MC_extensionality             ( : ) = 0   !=- default is 2

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

               FDE_approximation_ranges( 1, FDE_method_MD ) = 1d0   !=- 1d1 * GCC_radius_limit / 2**(GCC_generations+GCC_Super_generations-2)
               FDE_approximation_ranges( 2, FDE_method_MD ) = 0.8d0 * MC_target_radius   !=- 1d2 * GCC_radius_limit / 2**(GCC_generations+GCC_Super_generations-2)

               FDE_approximation_ranges( 1, FDE_method_intCD ) = 6d0
               FDE_approximation_ranges( 2, FDE_method_intCD ) = 3d1

               FDE_approximation_ranges( 1, FDE_method_diffCD ) = 5d0
               FDE_approximation_ranges( 2, FDE_method_diffCD ) = 4d1

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

            end

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine creating_MC_means
            integer i ; real(8) :: time = 0

            call default_MultiCantor

            do i=1,FDE_sequence_length ; write(*,'(/,A,i5)') 'FDE_means: set number =', i
                  catalog_prefix = MC_catalog_prefix
                  call make_FDE_seed( MC_FDE_seed_order )
               MC_catalog_name = trim( make_catalog_name() )

               call Generate_Multi_Cantor_Catalog(MC_catalog_name)

               call plot_catalog ( MC_catalog_name )

               call FDE_complex  ( MC_catalog_name )
               call FDE_GNUplots

               call FDE_means(FDE_sequence_length)

               time = time + FDE_time
            enddo

            FDE_D(:,:) = mean_FDE_D(:,:)   !=- for the graphs

            call write_means
            call plot_means

            write(*,*) '  MC_means: average working time on set = ' , int(time/FDE_sequence_length) , 's'
            write(*,*) '  MC_means: total working time = ' , int(time) , 's'

            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine FDE_testing_MultiCantor_script

               !call FDE_default_preset
               call FDE_testing_preset

               call creating_MC_means
               call default_means

            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

		subroutine input_catalogs_analysis

				call catalogs_reading

				call plot_catalog(	path_catalog_CF2	)
				call plot_catalog(	path_catalog_2MRS	)

			end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

end module
