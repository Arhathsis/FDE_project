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

      character(length) :: catalog_name , model_kind , sample_path , GCC_catalog_name , catalog_path
      character(length),allocatable,dimension(:) :: string_values

      integer ::        columns_numbers(7) = 1 , set_number = 1


	contains

      subroutine means_matrix

         do iloop=0,5
            do jloop=1,4
                call means(2d0, 5000 + iloop*2000, jloop*5)
            end do
         end do

      end subroutine

      subroutine means(D_arg, N_arg, ls_arg)
         real(8), intent(in) :: D_arg
         integer, intent(in) :: N_arg, ls_arg

                  !=- R_nn(geometries_count)
						!=- FDE_out_NP      ( 1 + 2*geometries_count , FDE_out_table_size ) , &
                  !=- FDE_out_MD      ( 1 + 2*geometries_count , FDE_out_table_size ) , &
                  !=- FDE_out_intCD   ( 1 + 2*geometries_count , FDE_out_table_size ) , &
                  !=- FDE_out_diffCD  ( 1 + 2*geometries_count , FDE_out_table_size ) , &

                                  FDE_recalculating = .true.
                              GCC_sets_regenerating = .true.
                                        GSCS_replot = .true.
                               GCC_sphere_selection = .true.

                                 GCC_catalog_prefix = 'Scantor'
                                 GCC_FDE_seed_order = 1

                              GCC_Fractal_Dimension = D_arg      !=- Milli-Millennium: 2.2     !=- it is may make as an multi-fractal array
                        GCC_Super_Fractal_Dimension = D_arg     !=- Milli-Millennium: 1.35

                                    GCC_generations = 6        !=- Milli-Millennium: 5
                              GCC_Super_generations = 5        !=- Milli-Millennium: 6

                     GCC_uniform_approach_parameter = 1        !=- Milli-Millennium: 1
               GCC_Super_uniform_approach_parameter = 0        !=- Milli-Millennium: 0

                                   GCC_radius_limit = 160d0     !=- Milli-Millennium: 31
                                       GCC_N_points = N_arg    !=- Milli-Millennium: 28k

                                    FDE_left_border = 20
                                   FDE_right_border = 50

                                 FDE_left_border_MD = GCC_radius_limit / 2**(GCC_generations+GCC_Super_generations-2)
                                FDE_right_border_MD = GCC_radius_limit

                                          RDR_z_max = 0.1  !=- to need automatize
                                           FDE_grid = 100
                                         set_number = ls_arg

            call creating_means(set_number)
            call default_means

         end subroutine



         subroutine creating_means(set_number)
            integer i, set_number

            do i=1,set_number ; write(*,*) 'FDE_means: ', i
                  call GCC_default ; call make_FDE_seed
               GCC_catalog_name = trim( make_catalog_name() )

               call Generate_Super_Cantor_Catalog(GCC_catalog_name)

                  if (GSCS_replot) call plot_catalog ( skeleton )
                  if (GSCS_replot) call plot_catalog ( GCC_catalog_name )

               call FDE_complex  ( GCC_catalog_name )
               call FDE_GNUplots

               call FDE_means(set_number)

            enddo

               call write_means(set_number)
               call plot_means

            end subroutine


		subroutine input_catalogs_analysis

				call catalogs_reading

				call plot_catalog(	path_catalog_CF2	)
				call plot_catalog(	path_catalog_2MRS	)

			end subroutine




      subroutine cantor_catalog_analysis  !=- test
         integer i

               catalog_prefix = 'cantor'

               FDE_seed_order = 2

               generations = 8
               Fractal_Dimension = 2.2
            !read(*,*) Fractal_Dimension
               N_points = 26d3

               radius_limit = 31d0
               FDE_recalculating = .true.

               FDE_left_border   = 1.6
               FDE_right_border  = 11

               FDE_left_border_MD   = 0.6
               FDE_right_border_MD  = 9

               RDR_z_max = 0.01

               FDE_grid=100

               do i=1,3
                  call make_FDE_seed
               catalog_name = trim( make_catalog_name() )
                  call Generate_Cantor_Catalog(catalog_name)   !=- input: exact the pathname of datafile

                  !call plot_catalog ( catalog_name )
                  call FDE_complex  ( catalog_name )
                  !call FDE_GNUplots

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



      subroutine scaling(start_scale, final_scale, scale_step, key)
         logical :: logscaling = .false.
         character(length) key
         real(8) start_scale, final_scale, scale_step
         integer i,N

         unit_1 = random_unit()

         open(unit_1,file=scaling_filepath,status='replace')
            write(unit_1,*) '# z & r^70(z) & d_L^70(z) & r^100(z) & d_L^100(z)'

            if (key=='log') logscaling = .true.
            if (logscaling) then

               N = 100   !=- the default log-scaling accuracy for graphs (uniform log grid)
               scale_step = ( log10(final_scale) - log10(start_scale) ) / N
               do i=0,N
                  z = 10**(log10(start_scale) + i*scale_step)
                  write(unit_1,'(F10.5,4(1x,F10.1))') z , R_LCDM(z,H70) , D_LCDM(z,H70) , R_LCDM(z,H100) , D_LCDM(z,H100)
                  end do

               else

               x = (final_scale-start_scale)/scale_step  !=- int() fix
               N = x !=- write(*,*) N, (final_scale-start_scale)/scale_step
               do i=0,N
                  z = start_scale + i*scale_step
                  write(unit_1,'(F10.5,4(1x,F10.1))') z , R_LCDM(z,H70) , D_LCDM(z,H70) , R_LCDM(z,H100) , D_LCDM(z,H100)
                  end do
!df
               end if

            close(unit_1)

         call plot_scaling(logscaling)

         end subroutine



         subroutine parameter_definition_interface

            write(*,*) 'enter file path (press 1 for default)'
            read(*,'(A256)') sample_path
               select case (sample_path)

                  case('1')

                     filepath='C:\Users\Arhath\YandexDisk\Science\DATA\test\data.dat'
                     model_kind='1'

                  case('m')

                     filepath=Millennium   !=- redshift x y z velX velY velZ mag_b mag_k
                     model_kind='m'

                  case('g')

                     filepath=Galacticus
                     model_kind='g'

                  case('gal')

                     filepath = Galacticus_160_dir
                     filename = 'mdpl2_galacticus_160.csv'
                     model_kind='g'

                  case default

                     filepath=sample_path
                        write(*,*) 'please, choose a cosmological model:'
                        write(*,*) 'press 1 for default: LCDM = 0, H0 = 70, w = -1, O_v = 0.7, O_k = 0.0'
                        write(*,*) 'press 2 to enter wCDM model parameters: H0 , w, O_v, O_k'
                        write(*,*) 'press 3 for gravitational cosmological redshift (not yet)'
                     read(*,*) model_kind

               end select
            inquire( file = filepath , exist = file_exists )

            if (file_exists) then ; 55 continue

                  select case (model_kind)
                     case('1')

                        RDR_type = 3   !=- 3 is LCDM model comoving r(z) | H0 = 70, O_v = 0.7

                     case('2')

                        RDR_type = 4   !=- 4 is wCDM model comoving r(z) | H0, w, O_v, O_k
                        write(*,*) 'please, specify wCDM parameters in the format: H0 w O_v O_k'
                        do
                           read(*,'(A256)') model_kind
                           read(model_kind,*,err=33,end=33) H0,w,O_v,O_k
                              goto 44 ; 33 continue
                                 write(*,*) 'add_visible_z error: 4xfloat: incorrect specification,' // &

                                    'try use the format 70 -1 0.7 0.0'
                           enddo ; 44 continue
                     case('m')
                        RDR_type = 4
                        H0=73d0 ; w=-1d0 ; O_v=0.75 ; O_k=0.0
                     case('g')
                        RDR_type = 4
                        H0=67.8 ; w=-1 ; O_v=0.693 ; O_k=0.0
                     case default

                        write(*,*) 'add_visible_z error: incorrect model kind' ; goto 55

                     end select
                  endif

            end subroutine



         subroutine add_visible_z
            integer  n, i_Size , i,j, i1,i2
            real(8) x,y,z,vx,vy,vz,Mag , scale_factor , shift_factor

               call parameter_definition_interface

            if (file_exists) then

               write(*,*) '   add_visible_z start'

                  scale_factor = 1d0 ; shift_factor = 0d0
                  do
                     select case(sample_path)
                        case('1')
                           columns = '2 3 4 5 6 7 7'
                        case('m')
                           columns = '2 3 4 5 6 7 8'
                           shift_factor = -31d0
                           radius_limit = 31d0
                           RDR_z_max = 0.01
                           FDE_grid=40
                        case('g')
                           columns = '3 4 5 6 7 8 10'
                           shift_factor = -80d0
                           radius_limit = 80
                        case default
                           write(*,*) 'please, specify columns (x,y,z) & (vx,vy,vz) & L, e.g., 1 2 3 4 5 6 7 (press 1 for default)'
                           read(*,'(A256)') columns
                     end select

                     read(columns,*,err=22,end=22) (columns_numbers(i),i=1,7)
                        goto 33 ; 22 continue
                           write(*,*) 'add_visible_z error: incorrect specification, try use the format 1 2 3 4 5 6 7'
                     enddo ; 33 continue

                  n = maxval(columns_numbers)
               allocate(string_values(n))

                  N_points = file_volume(filepath)

               add_visible_z_out = trim(folders( folder_Samples_add_visible_z )) // trim(name(filepath)) // '.dat'
                  unit_1 = random_unit() ; unit_2 = random_unit() ; unit_3 = random_unit()
               open( unit_1 , file = filepath          , status = 'old'    )
               open( unit_2 , file = add_visible_z_out , status = 'replace')
                  FDE_catalog = trim(folders( folder_Samples_add_visible_z )) // trim(name(filepath)) // '_FDE-catalog.dat'
               open( unit_3 , file = FDE_catalog , status = 'replace')

                  call prepare_percent(N_points)
                  index_info=.false.
                  do
                     line='#'
                     do while (symbol_search(line,'#'))
                        read(unit_1,'(A256)',end=11) line
                        enddo

                        string_values(:)=''
                        if (index_info) then
                           write(*,*) 'the first string of data:'
                           write(*,*) trim(line)
                              write(*,*) 'please, enter a separating symbol of the datafile'
                              read(*,*) columns
                           index_info=.false.
                           endif

                     if (columns==' ') then

                        read(line,*,err=55,end=55) (string_values(i),i=1,n) ; 55 continue

                        else

                           j=1 ; i1=1 ; i2=1
                           do i=1,length
                              if (line(i:i)==',') then
                                 i2=i-1
                                 string_values(j)=line(i1:i2)
                                 if (j==n-1) string_values(j+1)=line(i+1:length)
                                 j=1+j ; i1=i+1
                                 if (j>n) exit
                                 end if
                              end do

                     end if

                        x=0d0 ; y=0d0 ; z=0d0 ; vx = 0d0 ; vy = 0d0 ; vz = 0d0 ; L = 0d0
                     read(string_values(columns_numbers(1)),*) x
                     read(string_values(columns_numbers(2)),*) y
                     read(string_values(columns_numbers(3)),*) z
                     read(string_values(columns_numbers(4)),*) vx
                     read(string_values(columns_numbers(5)),*) vy
                     read(string_values(columns_numbers(6)),*) vz
                     read(string_values(columns_numbers(7)),*) Mag

                     x=x  * scale_factor + shift_factor   !=- kpc or Mpc
                     y=y  * scale_factor + shift_factor   !=- kpc or Mpc
                     z=z  * scale_factor + shift_factor   !=- kpc or Mpc

                     write(unit_2,*) (trim(string_values(i)),' ',i=1,n) , visible_redshift(x,y,z,vx,vy,vz)

                        !=- creating of a FDE catalog

                        if ( vr_distance < radius_limit ) then

                           catalog_line(N_col_cat_x) = x
                           catalog_line(N_col_cat_y) = y
                           catalog_line(N_col_cat_z) = z

                           call SphCoord( catalog_line(N_col_cat_x) , catalog_line(N_col_cat_y) , catalog_line(N_col_cat_z) , &
                              catalog_line(N_col_cat_l) , catalog_line(N_col_cat_b) , catalog_line(N_col_cat_dl) )

                              catalog_line( N_col_cat_Lum ) = Mag
                              catalog_line( N_col_cat_M   ) = Mag
                              catalog_line( N_col_cat_mag ) = &
                                 abs_to_vis_mag  ( catalog_line( N_col_cat_M ) , catalog_line( N_col_cat_dl ) )
                              catalog_line( N_col_cat_rs ) = vr_redshift

                           write( unit_3 , catalog_format ) ( catalog_line(j) , j=1,N_col_std )

                           endif

                        call write_percent
                     enddo ; 11 continue

                  close(unit_1) ; close(unit_2) ; close(unit_3)

                  deallocate(string_values)
                  N_points = file_volume(FDE_catalog)
                  write(*,*) '   add_visible_z succsessfully ended'
               else

                  write(*,*) 'add_visible_z error: file ' // trim(filepath) // ' do not exist'

               end if

            end subroutine



         subroutine add_visible_z_universe
            integer  n, i_Size , i,j, i1,i2 , catalog_index
            real(8) x,y,z,vx,vy,vz,Mag , scale_factor , shift_factor

               call parameter_definition_interface

            if (file_exists) then

               write(*,*) '   add_visible_z_universe started'

                  scale_factor = 1d0 ; shift_factor = 0d0
                  do
                     select case(sample_path)
                        case('1')
                           columns = '2 3 4 5 6 7 7'
                        case('m')
                           columns = '2 3 4 5 6 7 8'
                           shift_factor = -31d0
                           radius_limit = 31d0
                           RDR_z_max = 0.01
                           FDE_grid=40
                        case('gal')
                           columns = '2 3 4 5 6 7 8'
                           shift_factor = -80d0
                           radius_limit = 80
                        case default
                           write(*,*) 'please, specify columns (x,y,z) & (vx,vy,vz) & L, e.g., 1 2 3 4 5 6 7 (press 1 for default)'
                           read(*,'(A256)') columns
                     end select

                     read(columns,*,err=22,end=22) (columns_numbers(i),i=1,7)
                        goto 33 ; 22 continue
                           write(*,*) 'add_visible_z_universe error: incorrect specification, try use the format 1 2 3 4 5 6 7'
                     enddo ; 33 continue

                  n = maxval(columns_numbers)
               allocate(string_values(n))

            catalog_index = 0
            do
               catalog_index = 1+catalog_index

               catalog_path         = trim(filepath) // trim(adjustl(inttostr(catalog_index))) // '_' // trim(filename)
                  add_visible_z_out = trim(folders( folder_Samples_add_visible_z )) // trim(name(catalog_path)) // '.dat'
                  FDE_catalog       = trim(folders( folder_Samples_add_visible_z )) // trim(name(catalog_path)) &
                     // '_FDE-catalog.dat'

                  N_points = file_volume(catalog_path)

         inquire( file = FDE_catalog , exist = file_exists )

            if ( .not. file_exists ) then

                  unit_1 = random_unit() ; unit_2 = random_unit() ; unit_3 = random_unit()
               open( unit_1 , file = catalog_path        , status = 'old'    ,err=121)
               open( unit_2 , file = add_visible_z_out   , status = 'replace')
               open( unit_3 , file = FDE_catalog         , status = 'replace')

                  call prepare_percent(N_points)
                  index_info=.false.
                  do
                     line='#'
                     do while (symbol_search(line,'#'))
                        read(unit_1,'(A256)',end=11) line
                        enddo

                        string_values(:)=''
                        if (index_info) then
                           write(*,*) 'the first string of data:'
                           write(*,*) trim(line)
                              write(*,*) 'please, enter a separating symbol of the datafile'
                              read(*,*) columns
                           index_info=.false.
                           endif

                     if (columns==' ') then

                        read(line,*,err=55,end=55) (string_values(i),i=1,n) ; 55 continue

                        else

                           j=1 ; i1=1 ; i2=1
                           do i=1,length
                              if (line(i:i)==',') then
                                 i2=i-1
                                 string_values(j)=line(i1:i2)
                                 if (j==n-1) string_values(j+1)=line(i+1:length)
                                 j=1+j ; i1=i+1
                                 if (j>n) exit
                                 end if
                              end do

                     end if

                        x=0d0 ; y=0d0 ; z=0d0 ; vx = 0d0 ; vy = 0d0 ; vz = 0d0 ; L = 0d0
                     read(string_values(columns_numbers(1)),*) x
                     read(string_values(columns_numbers(2)),*) y
                     read(string_values(columns_numbers(3)),*) z
                     read(string_values(columns_numbers(4)),*) vx
                     read(string_values(columns_numbers(5)),*) vy
                     read(string_values(columns_numbers(6)),*) vz
                     read(string_values(columns_numbers(7)),*) Mag

                     x=x  * scale_factor + shift_factor   !=- kpc or Mpc
                     y=y  * scale_factor + shift_factor   !=- kpc or Mpc
                     z=z  * scale_factor + shift_factor   !=- kpc or Mpc

                     write(unit_2,*) (trim(string_values(i)),' ',i=1,n) , visible_redshift(x,y,z,vx,vy,vz)

                        !=- creating of a FDE catalog

                        if ( vr_distance < radius_limit ) then

                           catalog_line(N_col_cat_x) = x
                           catalog_line(N_col_cat_y) = y
                           catalog_line(N_col_cat_z) = z

                           call SphCoord( catalog_line(N_col_cat_x) , catalog_line(N_col_cat_y) , catalog_line(N_col_cat_z) , &
                              catalog_line(N_col_cat_l) , catalog_line(N_col_cat_b) , catalog_line(N_col_cat_dl) )

                              catalog_line( N_col_cat_Lum ) = Mag
                              catalog_line( N_col_cat_M   ) = Mag
                              catalog_line( N_col_cat_mag ) = &
                                 abs_to_vis_mag  ( catalog_line( N_col_cat_M ) , catalog_line( N_col_cat_dl ) )
                              catalog_line( N_col_cat_rs ) = vr_redshift

                           write( unit_3 , catalog_format ) ( catalog_line(j) , j=1,N_col_std )

                           endif

                        call write_percent
                     enddo ; 11 continue

                  close(unit_1) ; close(unit_2) ; close(unit_3)

                  N_points = file_volume(FDE_catalog)
                  final_counter = N_points
                  call plot_catalog ( FDE_catalog )

               endif

                           FDE_left_border   = 10
                           FDE_right_border  = 30
                           FDE_left_border_MD = 10
                           FDE_right_border_MD = 40
                           radius_limit = 80d0
                           FDE_grid = 100

                     FDE_recalculating = .false.
                     call FDE_complex  ( FDE_catalog )
                     call FDE_GNUplots

               enddo ; 121 continue

                  deallocate(string_values)
                  write(*,*) '   add_visible_z succsessfully ended'
               else

                  write(*,*) 'add_visible_z_universe error: file ' // trim(filepath) // ' do not exist'

               end if

            end subroutine



         subroutine analysis_external_FDE_catalog

               FDE_recalculating = .true.

                  call parameter_definition_interface

               if (file_exists) then ; 11 continue

                     select case(sample_path)

                        case('m')
                           FDE_left_border   = 3
                           FDE_right_border  = 15
                           FDE_left_border_MD = 2
                           FDE_right_border_MD = 15
                           radius_limit = 31d0
                           FDE_grid = 100

                        case('gal')
                           FDE_left_border   = 3
                           FDE_right_border  = 15
                           FDE_left_border_MD = 2
                           FDE_right_border_MD = 15
                           radius_limit = 80d0
                           FDE_grid = 100

                        case default
                              write(*,*) 'please, enter maximum sample radius'
                           read(*,*,err=11) radius_limit
                              write(*,*) 'please enter FDE_LSM_border through space'
                           read(*,*,err=11,end=11) FDE_left_border, FDE_right_border

                     end select

                     FDE_catalog = trim(folders( folder_Samples_add_visible_z )) // trim(name(filepath)) // '_FDE-catalog.dat'

                     call FDE_complex  ( FDE_catalog )
                     call FDE_GNUplots

                  else

                     write(*,*) 'add_visible_z error: file ' // trim(filepath) // ' do not exist'

                  end if
            end subroutine analysis_external_FDE_catalog



         subroutine Super_Cantor_analysis(D, N, catalog_name_out)
            character(length), allocatable, dimension(:), intent(out) :: catalog_name_out
            integer :: peremennaya_cycla
            real(8) D

               allocate(catalog_name_out(N))

                                  FDE_recalculating = .true.
                              GCC_sets_regenerating = .true.
                                        GSCS_replot = .true.
                               GCC_sphere_selection = .true.

                                 GCC_catalog_prefix = 'Scantor_task1'
                                 GCC_FDE_seed_order = 4

                              GCC_Fractal_Dimension = D      !=- Milli-Millennium: 2.2     !=- it is may make array of multi-fractal
                        GCC_Super_Fractal_Dimension =  GCC_Fractal_Dimension   !=- Milli-Millennium: 1.35

                                    GCC_generations = 5        !=- Milli-Millennium: 5
                              GCC_Super_generations = 6        !=- Milli-Millennium: 6

                     GCC_uniform_approach_parameter = 1        !=- Milli-Millennium: 1
               GCC_Super_uniform_approach_parameter = 0        !=- Milli-Millennium: 0

                                   GCC_radius_limit = 160    !=- Milli-Millennium: 31
                                       GCC_N_points = 1d4    !=- Milli-Millennium: 28k

                                    FDE_left_border = 3
                                   FDE_right_border = GCC_radius_limit*0.5d0

                                 FDE_left_border_MD = GCC_radius_limit*0.01d0
                                FDE_right_border_MD = GCC_radius_limit*0.5d0

                                          RDR_z_max = 0.1  !=- to need automatize
                                           FDE_grid = 100
                                         set_number = N

         do peremennaya_cycla=1,set_number
                           write(*,*) 'test1 ', peremennaya_cycla
                  call GCC_default ; call make_FDE_seed
               GCC_catalog_name = trim( make_catalog_name() )
               catalog_name_out(peremennaya_cycla) = GCC_catalog_name

               call Generate_Super_Cantor_Catalog(GCC_catalog_name)
                  if (GSCS_replot) call plot_catalog ( skeleton )
                  if (GSCS_replot) call plot_catalog ( GCC_catalog_name )

               call FDE_complex  ( GCC_catalog_name )
               call FDE_GNUplots
            enddo

            end subroutine



         subroutine example1

            integer :: i

            !do i=1,10
                !call Super_Cantor_analysis(1.0 + i*0.15, 3)
                !call Super_Cantor_analysis(1.0 + i*0.15, 8)
                !call Super_Cantor_analysis(1.0 + i*0.15, 21)
                !call Super_Cantor_analysis(1.0 + i*0.15, 55)
            !enddo
            !call Super_Cantor_analysis(2.0, 5, catalog_name_out)

         end subroutine example1



         subroutine mean_MD(D, N)
            real(8), allocatable, dimension(:) :: Rmin, Rmax, R_setka, A_mean, B_mean, North_mean, &
                                                  South_mean, Adel_mean, Bdel_mean, Ndel_mean, Sdel_mean
            real(8), allocatable, dimension(:,:) :: A, B, North, South, A_delta, B_delta, N_delta, S_delta
            real(8), dimension(9) :: line_from_file
            character(length), allocatable, dimension(:) :: catalog_name_out1
            character(length) :: catalog_name_out2
            real(8) :: R_right_border, R_left_border , D !=- the input parameters must be declared!
            integer :: i,j,N,nstrings    !=- yes, the inner cycle indexes must be declared

            allocate( catalog_name_out1(N) , Rmin(N) , Rmax(N) )

            call Super_Cantor_analysis(D, N, catalog_name_out1)

            do i=1,N
               write(*,*) catalog_name_out1(i)
            enddo

            do i=1,N
               !catalog_name_out2 = catalog_name_out1(i)
               !catalog_name_out2 = catalog_name_out2(23:53)
               catalog_name_out2 = 'Main_Workspace/MD/add-files'//trim(catalog_name_out2)//'_100_MD.dat'
               catalog_name_out1(i) = catalog_name_out2
               write(*,*)catalog_name_out2

               open(1,file=catalog_name_out2,status='replace') !=- just in case, the open's status should be typed
               read(1,*) line_from_file
               Rmin(i) = line_from_file(1)
               do
                  read(1,*, end=99) line_from_file
               end do
            99 Rmax(i) = line_from_file(1)
               close(1)
            enddo

            R_left_border = maxval(Rmin)
            R_right_border = minval(Rmax)

            open(1,file=catalog_name_out2)
            nstrings = 0
            do
               read(1,*, end=97) line_from_file
               if ((line_from_file(1) >= R_left_border).and.(line_from_file(1) <= R_right_border)) then
                   nstrings = nstrings+1
               end if
            end do
         97 close(1)

            allocate(R_setka(nstrings))
            allocate(A(N,nstrings))
            allocate(B(N,nstrings))
            allocate(North(N,nstrings))
            allocate(South(N,nstrings))
            allocate(A_delta(N,nstrings))
            allocate(B_delta(N,nstrings))
            allocate(N_delta(N,nstrings))
            allocate(S_delta(N,nstrings))
            allocate(A_mean(nstrings))
            allocate(B_mean(nstrings))
            allocate(North_mean(nstrings))
            allocate(South_mean(nstrings))
            allocate(Adel_mean(nstrings))
            allocate(Bdel_mean(nstrings))
            allocate(Ndel_mean(nstrings))
            allocate(Sdel_mean(nstrings))

            open(1,file=catalog_name_out2)
            i = 1
            do
               read(1,*, end=98) line_from_file
               if ((line_from_file(1) >= R_left_border).and.(line_from_file(1) <= R_right_border)) then
                   R_setka(i) = line_from_file(1)
                   i = i+1
               end if
            end do
         98 close(1)
            write(*,*) R_setka

            do i=1,N
                open(1,file=catalog_name_out1(i))
                j=0
                do
                    read(1,*, end=105) line_from_file
                    if ((line_from_file(1) >= R_left_border).and.(line_from_file(1) <= R_right_border)) then
                        j = j+1
                        A(i,j) = line_from_file(2)
                        A_delta(i,j) = line_from_file(3)
                        North(i,j) = line_from_file(4)
                        N_delta(i,j) = line_from_file(5)
                        South(i,j) = line_from_file(6)
                        S_delta(i,j) = line_from_file(7)
                        B(i,j) = line_from_file(8)
                        B_delta(i,j) = line_from_file(9)
                     end if
                end do
            105 close(1)
            end do

            write(*,*) 'testtest'   !=- here, the pause command can be used, e.g. pause 1, pause 23, etc.

            write(*,*) A(:,1)

            do j=1,nstrings
                A_mean(j) = sum(A(:,j)) / N
                B_mean(j) = sum(B(:,j)) / N
                North_mean(j) = sum(North(:,j)) / N
                South_mean(j) = sum(South(:,j)) / N
                Adel_mean(j) = sum(A_delta(:,j)) / N
                Bdel_mean(j) = sum(B_delta(:,j)) / N
                Ndel_mean(j) = sum(N_delta(:,j)) / N
                Sdel_mean(j) = sum(S_delta(:,j)) / N
            end do

            write(*,*) 'testtest2'

            open(2, file='Main_Workspace/MD/add-files/Mean_values.dat',status='replace')
            do j=1,nstrings
                write(2,FDE_data_format) R_setka(j), A_mean(j), Adel_mean(j), North_mean(j), &
                                         Ndel_mean(j), South_mean(j), Sdel_mean(j), B_mean(j), Bdel_mean(j)
            end do
            close(2)
         end subroutine mean_MD

         subroutine create_table_of_means(D_left, D_right, D_amount, N_left, N_right, N_amount, ls_left, ls_right, ls_amount)

            integer, intent(in) :: D_amount, N_left, N_right, N_amount, ls_left, ls_right, ls_amount
            real(8), intent(in) :: D_left, D_right
            integer, dimension(N_amount) :: N_mass
            real(8), dimension(D_amount) :: D_mass
            integer, dimension(ls_amount) :: ls_mass
            real(8), dimension(D_amount,N_amount,ls_amount,geometries_count) :: table

            if (D_amount == 1) then
                D_mass(1) = D_left
            else
                do iloop = 1,D_amount
                    D_mass(iloop) = D_left + (iloop - 1) * (D_right - D_left) / (D_amount - 1)
                end do
            end if

            do jloop = 1,N_amount
                N_mass(jloop) = floor(10 ** (log10(N_left * 1.0) + (jloop - 1)&
                 * (log10(N_right * 1.0) - log10(N_left * 1.0)) / (N_amount - 1)) + 0.5)
            end do

            do kloop = 1,ls_amount
                ls_mass(kloop) = floor(10 ** (log10(ls_left * 1.0) + (kloop - 1)&
                 * (log10(ls_right * 1.0) - log10(ls_left * 1.0)) / (ls_amount - 1)) + 0.5)
            end do

            do iloop = 1,D_amount
                do jloop = 1,N_amount
                    do kloop = 1,ls_amount
                        call means(D_mass(iloop), N_mass(jloop), ls_mass(kloop))
                        table(iloop,jloop,kloop,1) = FDE_D(1)
                        table(iloop,jloop,kloop,2) = FDE_D(2)
                        table(iloop,jloop,kloop,3) = FDE_D(3)
                        table(iloop,jloop,kloop,4) = FDE_D(4)
                    end do
                end do
            end do

            do iloop = 1,D_amount
                open(3, file='Main_Workspace/Statistics/add-files/D_table_' //&
                 trim(realtostr(D_mass(iloop))) // '_A.dat',status='replace')
                write(3,*) 'ls|N', N_mass
                do jloop = 1,ls_amount
                    write(3,*) ls_mass(jloop), table(iloop,:,jloop,1)
                end do
                close(3)
                open(4, file='Main_Workspace/Statistics/add-files/D_table_' //&
                 trim(realtostr(D_mass(iloop))) // '_N.dat',status='replace')
                write(4,*) 'ls|N', N_mass
                do jloop = 1,ls_amount
                    write(4,*) ls_mass(jloop), table(iloop,:,jloop,2)
                end do
                close(4)
                open(5, file='Main_Workspace/Statistics/add-files/D_table_' //&
                 trim(realtostr(D_mass(iloop))) // '_S.dat',status='replace')
                write(5,*) 'ls|N', N_mass
                do jloop = 1,ls_amount
                    write(5,*) ls_mass(jloop), table(iloop,:,jloop,3)
                end do
                close(5)
                open(6, file='Main_Workspace/Statistics/add-files/D_table_' //&
                 trim(realtostr(D_mass(iloop))) // '_B.dat',status='replace')
                write(6,*) 'ls|N', N_mass
                do jloop = 1,ls_amount
                    write(6,*) ls_mass(jloop), table(iloop,:,jloop,4)
                end do
                close(6)
            end do

end subroutine create_table_of_means

end module
