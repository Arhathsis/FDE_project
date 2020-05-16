!=- Fractal Dimension Estimation
!=- © Stanislav Shirokov, 2014-2020

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

      character(len) :: catalog_prefix = 'catalog' , catalog_name

	contains



      subroutine mean_dimension

                  !=- R_nn(geometries_count)
						!=- FDE_out_NP      ( 1 + 2*geometries_count , FDE_out_table_size ) , &
                  !=- FDE_out_MD      ( 1 + 2*geometries_count , FDE_out_table_size ) , &
                  !=- FDE_out_intCD   ( 1 + 2*geometries_count , FDE_out_table_size ) , &
                  !=- FDE_out_diffCD  ( 1 + 2*geometries_count , FDE_out_table_size ) , &

         end subroutine



		subroutine input_catalogs_analysis

				call catalogs_reading

				call plot_catalog(	path_catalog_CF2	)
				call plot_catalog(	path_catalog_2MRS	)

			end subroutine


      character(len) function make_catalog_name()
         make_catalog_name=''
         make_catalog_name = trim( folders( folder_Samples ) ) // trim(catalog_prefix) // '-' // &
            trim(inttostrf(FDE_seed,FDE_seed_order)) // '_' // trim(adjustl(inttostr(generations))) // '-' // &
            trim(adjustl(realtostrf(Fractal_Dimension,3,1))) // '-' // trim(adjustl(inttostr(N_points))) // '.dat'
         end function make_catalog_name



      subroutine cantor_catalog_analysis  !=- test
         integer i

               catalog_prefix = 'cantor'

               FDE_seed_order = 1

               generations = 8
               Fractal_Dimension = 2.0
            !read(*,*) Fractal_Dimension
               N_points = 5d3

               radius_limit = 2d2
               FDE_recalculating = .true.

               FDE_left_border   = 10
               FDE_right_border  = 100

               do i=1,1
                   call make_FDE_seed
               catalog_name = trim( make_catalog_name() )
                  call Generate_Cantor_Catalog(catalog_name)   !=- input: exact the pathname of datafile

                  call plot_catalog ( catalog_name )
                  call FDE_complex  ( catalog_name )
                  call FDE_GNUplots

                  enddo

         end subroutine



      subroutine uniform_catalog_analysis !=- test
         integer i

               catalog_prefix = 'uniform'

               FDE_seed_order = 1

               Fractal_Dimension = 3.0
               N_points = 5d3

               radius_limit = 2d2
               FDE_recalculating = .true.

               FDE_left_border   = 10
               FDE_right_border  = 100

               do i=1,2
                   call make_FDE_seed
               catalog_name = make_catalog_name()
                  call Generate_Uniform_Catalog(catalog_name)   !=- input: exact the pathname of datafile

                  call plot_catalog ( catalog_name )
                  call FDE_complex  ( catalog_name )
                  call FDE_GNUplots

                  enddo

         end subroutine




end module
