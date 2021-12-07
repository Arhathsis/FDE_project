!=- Fractal Dimension Estimation
!=- © Stanislav Shirokov, 2014-2020

module FDE_generators
   use FDE_paths
   use FDE_config

   use global
   use cosmology
   use math

   character(length) :: catalog_prefix = 'catalog' , secondary_catalog , skeleton , GCC_catalog_prefix , &
                        MC_catalog_prefix , MC_catalog_name

   logical ::  FDE_generator_report_flag = .true. , file_end = .false. , sphere_selection = .true. , &
               GCC_datafile = .true. , GSCS_replot = .true. , GCC_sets_regenerating , GCC_FDE_generator_report_flag, &
               GCC_sphere_selection

   real(8) ::  Lum, mag , catalog_line(N_col_std)=0d0 , chance , noise , npoint_noise , Fractal_Dimension , &
               nearest_point_to_origins(3) , min_dist , x_shift=0d0 , y_shift=0d0 , z_shift=0d0 , Super_radius_limit, &
               GCC_Super_Fractal_Dimension , GCC_Fractal_Dimension , GCC_z_shift , GCC_y_shift , GCC_x_shift , &
               GCC_radius_limit , MC_Fractal_Dimension( 100 ) , MC_target_radius , FDE_approximation_ranges(2,100) , &
               FDE_Cantor_catalog(N_col_std,max_set_size) , FDE_MultiCantor_catalog(N_col_std,max_set_size) , radius_scale , &
               FDE_set_power , FDE_probability , FDE_time

   integer ::  fractionality , extensionality , gen_start , iteration , generations , uniform_approach_parameter , &
               step , imax , imin , jmin , jmax , kmin , kmax, GCC_unit=0, skeleton_generations, &
               unit_cyrle , N_cyrle , FDE_seed = 1 , FDE_seed_order = 1, GSCS_N, GCC_Super_generations, &
               GCC_Super_uniform_approach_parameter , GCC_Super_N_points , Super_FDE_seed_order , GCC_generations , &
               GCC_N_points , GCC_uniform_approach_parameter , GCC_FDE_seed_order , Super_counter , MC_FDE_seed_order , &
               FDE_multiness , MC_generations( 100 ) , MC_extensionality( 100 ) , MC_uniform_approach_parameter( 100 ) , &
               MC_target_point_number , FDE_sequence_length , MC_counter

   contains

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine default_MultiCantor
            integer i

               catalog_prefix = MC_catalog_prefix // '_' // trim(inttostr(FDE_multiness))
               if ( MC_FDE_seed_order<1 .or. MC_FDE_seed_order>5 ) MC_FDE_seed_order=3
               FDE_seed_order = MC_FDE_seed_order

               do i=1,100
                  if (MC_Fractal_Dimension(i)<1 .or. MC_Fractal_Dimension(i)  >3 ) MC_Fractal_Dimension(i)  = 2d0
                  if (MC_generations(i)      <1 .or. MC_generations(i)        >10) MC_generations(i)        = 3
                  if (MC_extensionality(i)   <2 .or. MC_extensionality(i)     >10) MC_extensionality(i)     = 2
                  if (MC_uniform_approach_parameter(i)<0 .and. MC_uniform_approach_parameter(i)>10) &
                      MC_uniform_approach_parameter(i)=0  !=- may be 0,1,2, ...
                  enddo

               if ( FDE_multiness < 1 .or. FDE_multiness > 100 ) FDE_multiness = 2

               if (MC_target_radius       <=0) MC_target_radius         = 1d2
               if (MC_target_point_number <=0) MC_target_point_number   = 1d4

               Fractal_Dimension          = MC_Fractal_Dimension          ( 1 )
               generations                = MC_generations                ( 1 )
               extensionality             = MC_extensionality             ( 1 )
               uniform_approach_parameter = MC_uniform_approach_parameter ( 1 )

               target_point_number = MC_target_point_number
               if (FDE_probability<=0 .or. FDE_probability>1d0) FDE_probability = 1d0
               if (FDE_sequence_length<1) FDE_sequence_length=1

            end

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

      subroutine Generate_Multi_Cantor_Catalog(filepath)
         character(length) filepath
         integer i,j,k,l , multiness_power_N , fragmentation
         real(8) XYZ(N_col_std,max_set_size) , probability , x , mean_count , N_max , start_time , final_time

            inquire(file=slashfix(filepath), exist=file_exists)
            if ( file_exists .and. .not. sets_regenerating ) then
               final_counter = file_volume(filepath)
               else

               call cpu_time(start_time)

               !=- computation of set power N_max
               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
               !=- example,  2x3 * 3x4 * 2x3 | f(100) = 0.02 | f=100/(2**3 * 3**4 * 2**3) | N_max = 125kk

               fragmentation=1
               do i=1,FDE_multiness
                  fragmentation = fragmentation * MC_extensionality(i)**MC_generations(i) !=- accounting of extensionality and generations
                  enddo
               N_max = (1d0*fragmentation)**3 !=- accounting of 3D space | N_max - count of V_N
               write(*,'(/,A,i10,A,i20)') '   MC: total fractionality = ' , fragmentation , ', MC: N_max = ' , int(N_max)
                  if (N_max>8000000) write(*,*) '  MC:    attention: too large N_max'

               !=- the first skeleton creating
               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

               call set_MC_parameters(1)  !=- parameter setting to the first level of multiness

               call Generate_Cantor_Catalog(skeleton) !=- without a datafile
                  FDE_generator_report_flag = .false. !=- to hide the set parameters in terminal

               XYZ(:,:)=0d0
               multiness_power_N = 0
               x=0
               do k = 1, final_counter !=- , FDE_Cantor_catalog(1,k)>0

                     if ( .not. FDE_selected_generators ) call random_number(x)
                  if (x<FDE_probability) then
                     multiness_power_N=1+multiness_power_N
                     XYZ(:,multiness_power_N) = FDE_Cantor_catalog(:,k) !=- buffering of a fractal
                     endif
                  end do

               write(*,'("   MC: multiness = ",i3,", number of points: ",i10)') 1 , final_counter
               write(*,*) '  MC:    estimation of the number of points at next SAME multiness ' , &
                  int(1d0*final_counter**2*FDE_probability)

               !=- the MultiCantor creating
               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

               do i=2,FDE_multiness  !=- write(*,*) 'MC: multiness = ', i

                  call set_MC_parameters(i)  !=- parameter setting to the i level of multiness , including radius_scale

                  mean_count=0d0
                  MC_counter = 0
                  do j=1,multiness_power_N   !=- number of the multifractal points at i-1

                     if (MC_counter>max_set_size) then
                        write(*,*) '  MC: MC_counter > max_set_size (2^23)' ; exit
                        endif

                        x_shift=XYZ(1,j)
                        y_shift=XYZ(2,j)
                        z_shift=XYZ(3,j)

                     call Generate_Cantor_Catalog(skeleton)

                        x=0
                     do k = 1, final_counter !=- , FDE_Cantor_catalog(1,k)>0

                           if ( .not. FDE_selected_generators ) call random_number(x)
                        if ( x<FDE_probability .and. FDE_Cantor_catalog( N_col_cat_dl , k ) < MC_target_radius ) then !=- it is required for more completeness of fractal ball
                           MC_counter=1+MC_counter
                           FDE_MultiCantor_catalog(:,MC_counter) = FDE_Cantor_catalog(:,k) !=- buffering of a fractal
                           mean_count = mean_count + 1d0*final_counter
                           endif
                        end do

                     enddo

               write(*,'("   MC: multiness = ",i3,", number of points: ",i10)') i , MC_counter   !=- final MC_counter should be about MC_target_point_number
               write(*,'("   MC: selection probability = ",F5.3,", secondary fractal mean number of points = ",i10)') &
                  FDE_probability , int(1d0*mean_count/MC_counter)
               if (i<FDE_multiness) write(*,*) '  MC:    estimation of the number of points at next SAME multiness ' , &
                  int(1d0*mean_count*FDE_probability)

                  XYZ(:,:)=0d0   !=- , FDE_MultiCantor_catalog(1,k)>0
                  forall (k=1:MC_counter ) XYZ(:,k)=FDE_MultiCantor_catalog(:,k) !=- buffering of a multifractal
                  multiness_power_N = MC_counter

                  enddo

               !=- saving in FDE folder to the permanent memory
               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

                  probability = 1d0 * MC_target_point_number / MC_counter

                  unit_1 = random_unit()
                  open(unit_1,file=filepath,status='replace')
                     final_counter=0
                     do i=1,MC_counter
                        call random_number(x)
                        if ( x < probability ) then   !=- this is the main probability
                           final_counter = 1+final_counter
                           write(unit_1,*) FDE_MultiCantor_catalog(:,i)
                           endif
                        enddo
                     close(unit_1)

               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
                  call set_MC_parameters(1)  !=- for final parameters in the graphs
                  write(*,'(A,i10,A,F10.4)') '   MC: final number of points = ' , final_counter , &
                     ', power_index = ' , FDE_set_power
                  call cpu_time(final_time) ; FDE_time = final_time-start_time
                  write(*,*) '  MC: working time: ', int(FDE_time) , 's'
               end if
         end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

      subroutine set_MC_parameters(multiness_index)
         integer i, multiness_index , fragmentation

            sphere_selection = .true.
            if (multiness_index>1) sphere_selection = .false.

            Fractal_Dimension          = MC_Fractal_Dimension          ( multiness_index )
            generations                = MC_generations                ( multiness_index )
            extensionality             = MC_extensionality             ( multiness_index )
            uniform_approach_parameter = MC_uniform_approach_parameter ( multiness_index )

            target_point_number = MC_target_point_number

            GCC_datafile = .false.  !=- do not create datafiles
               catalog_prefix = 'skeleton'
            skeleton = trim( make_catalog_name() ) !=- it is just required to correct working `Generate_Cantor_Catalog`

            x_shift=0
            y_shift=0
            z_shift=0
            radius_limit = MC_target_radius

            if (multiness_index>1) then
               fragmentation = 1
               do i=1,multiness_index-1
                  fragmentation = fragmentation * MC_extensionality(i)**MC_generations(i) !=- accounting of extensionality and generations
                  enddo
               radius_scale = MC_target_radius / fragmentation
               else
            radius_scale = MC_target_radius
            end if

            npoint_noise = 0d0      !=- shifting of cantor set points
            if (multiness_index == FDE_multiness) npoint_noise = 1d0

            FDE_set_power = 1d0 * final_counter / (1+MC_counter)

         end

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

      subroutine Generate_Cantor_Catalog(path) !=- Cantor Sets Generator
         character(length) path ; integer i,j,k , ii,jj,kk , lvl
         real(8) x,y,z
         integer(1),allocatable,dimension(:,:,:) :: cube

               inquire(file=slashfix(path), exist=file_exists)
            if ( file_exists .and. .not. sets_regenerating ) then
               if (FDE_generator_report_flag) write(*,*) '   cantor set ',trim(path),' already exists'
               final_counter = file_volume(path)
               else

                  call cpu_time(start_time)
                     11 continue
                  call FDE_seed_fix

                     if (GCC_datafile)  then
                        unit_3 = random_unit()
                        open(unit_3,file=slashfix(path),status='replace')
                           write(unit_3,catalog_heads_format  ) (catalog_titles(j),j=1,N_col_std)
                           write(unit_3,catalog_heading_format) (catalog_titles(j),j=1,N_col_std)
                        endif

                        if (extensionality<2) extensionality=2

                        fractionality=extensionality**generations
                     allocate( cube(fractionality,fractionality,fractionality) )
                        forall (i=1:fractionality,j=1:fractionality,k=1:fractionality) cube(i,j,k)=1

                     do iteration=1,generations
                        lvl=extensionality**iteration
                        if (iteration.gt.uniform_approach_parameter) then ; step=fractionality/lvl
                           do ii=1,lvl ; do jj=1,lvl ; do kk=1,lvl

                              imin=(ii-1)*step+1 ; imax=ii*step
                              jmin=(jj-1)*step+1 ; jmax=jj*step
                              kmin=(kk-1)*step+1 ; kmax=kk*step
                              call random_number(chance)

                              if ( chance.gt.2**(Fractal_Dimension-3d0) .or. &
                                 FDE_selected_generators .and. chance.gt.FDE_probability ) &
                                 forall (i=imin:imax,j=jmin:jmax,k=kmin:kmax) cube(i,j,k)=0

                              enddo ; enddo ; enddo

                           endif
                        enddo

                     counter=0
                     do ii=1,fractionality ; do jj=1,fractionality ; do kk=1,fractionality
                        if ( cube(ii,jj,kk)==1 ) then
                           x = ((ii-0.5)*1d0/fractionality-0.5d0)*2
                           y = ((jj-0.5)*1d0/fractionality-0.5d0)*2
                           z = ((kk-0.5)*1d0/fractionality-0.5d0)*2

                           if ( distance(x,y,z) > 1d0+1.5d0/fractionality .and. sphere_selection .or. &
                           dabs(x)>1 .or. dabs(y)>1 .or. dabs(z)>1) then
                              cube(ii,jj,kk) = 0
                              else

                                 if ( dabs(x)<1 .and. dabs(y)<1 .and. dabs(z)<1 .and. .not. sphere_selection &
                                    .or. (sphere_selection .and. distance(x,y,z)<1 ) ) counter=counter+1
                              endif
                           endif
                        enddo ; enddo ; enddo
                        if (GCC_datafile) write(unit_3,'("# an initial count of points ",i7)') counter

                     FDE_Cantor_catalog(:,:)=0d0
                     final_counter=0 ; Super_counter = 0
                     do ii=1,fractionality ; do jj=1,fractionality ; do kk=1,fractionality
                        if ( cube(ii,jj,kk)==1 ) then

                           call npoint( ii, jj, kk, &
                              catalog_line(N_col_cat_x), catalog_line(N_col_cat_y), catalog_line(N_col_cat_z) )

                              catalog_line(N_col_cat_x) = catalog_line(N_col_cat_x)*radius_scale + x_shift
                              catalog_line(N_col_cat_y) = catalog_line(N_col_cat_y)*radius_scale + y_shift
                              catalog_line(N_col_cat_z) = catalog_line(N_col_cat_z)*radius_scale + z_shift

                              call SphCoord( catalog_line(N_col_cat_x) , catalog_line(N_col_cat_y) , &
                                 catalog_line(N_col_cat_z) , &
                                 catalog_line(N_col_cat_l) , catalog_line(N_col_cat_b) , &
                                 catalog_line(N_col_cat_dl) )

                           if (  dabs(catalog_line(N_col_cat_x))<radius_limit &
                              .and. dabs(catalog_line(N_col_cat_y))<radius_limit &
                              .and. dabs(catalog_line(N_col_cat_z))<radius_limit &
                              .and. .not. sphere_selection  &
                              .or.  sphere_selection .and. catalog_line(N_col_cat_dl) < radius_limit  ) then

                                 call make_luminosity( catalog_line(N_col_cat_Lum) )

                                 catalog_line( N_col_cat_M   ) = Lum_to_Mag      ( catalog_line( N_col_cat_Lum ) )
                                 catalog_line( N_col_cat_mag ) = &
                                    abs_to_vis_mag  ( catalog_line( N_col_cat_M ) , catalog_line( N_col_cat_dl ) )
                                 catalog_line( N_col_cat_rs ) = fun_from_R_to_z ( catalog_line( N_col_cat_dl ) )

                                 final_counter=final_counter+1
                                 FDE_Cantor_catalog(:,final_counter) = catalog_line(:)

                                 if (GCC_datafile) write( unit_3 , catalog_format ) ( catalog_line(j) , j=1,N_col_std )

                              endif
                           end if
                        enddo ; enddo ; enddo !=- write(*,*) final_counter , Super_counter , counter , target_point_number ; pause

                     if (GCC_datafile) close(unit_3)
                     deallocate(cube)

                     if (final_counter==0) goto 11    !=- a fix of null set

               endif

               call generation_report ; if (FDE_generator_report_flag) write(*,*) '  GCS: datafile: ' , trim(name(path))
               call cpu_time(final_time)

         end subroutine Generate_Cantor_Catalog

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

      subroutine make_FDE_seed(order)
         real(8) x ; integer order
            do
                  call random_number(x)
               FDE_seed = int( x*1d1**order )
               if (FDE_seed < 1d1**order .and. FDE_seed > 1d1**(order-1)-1 ) exit
               enddo
      end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

      subroutine random_number_centered(x)
         real(8) x ; call random_number(x) ; x=(x-0.5d0)*2
         end subroutine

      subroutine random_number_centered_normed(x,norm)
         real(8) x,norm ; call random_number(x) ; x=(x-0.5d0)*2*norm
         end subroutine

      real(8) function log_uniform_luminosity()
            call random_number(log_uniform_luminosity) ; log_uniform_luminosity = 1d1 ** ( log10(minimal_luminosity) + &
               log_uniform_luminosity * ( log10(maximal_luminosity) - log10(minimal_luminosity) ) )
         end function

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

      subroutine make_luminosity( Luminosity )
         real(8) Luminosity

            Luminosity = log_uniform_luminosity()

            select case ( luminosity_model )

               case default

            end select
         end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

      subroutine FDE_seed_fix
         real(8) x ; integer i
            do i=1,FDE_seed
               call random_number(x)
               end do
         end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

      subroutine generation_report
            if ( FDE_generator_report_flag ) then
               write(*,'(3x,A,5x,F3.1 )' ) 'GCS: D_eff:' , Fractal_Dimension
               write(*,'(3x,A,1x,i7     )' ) 'GCS: N_eff:' , target_point_number
               write(*,'(3x,A,2x,i6,A   )' ) 'GCS: R_max:' , int(radius_limit) , ' Mpc'  !=- !write(*,'(3x,A,1x,E10.3 )' ) 'Radius:'       , radius_limit
               !write(*,'(3x,A,1x,i7    )' ) 'N_p:  ' , final_counter
               endif
         end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

   	subroutine Generate_Uniform_Catalog(path)
         character(length) path ; integer j

            write(*,'(/,A)') 'GUC is started..'

               inquire(file=slashfix(path), exist=file_exists)
            if ( file_exists .and. .not. sets_regenerating ) then
               !=- write(*,*) '   uniform set ',trim(path),' already exists'
               else

                  call FDE_seed_fix

                  unit_4 = random_unit()

                  if (target_point_number<1) target_point_number=100
                  open(unit_4,file=slashfix(path),status='replace')
                     write(unit_4,catalog_heads_format  ) (catalog_titles(j),j=1,N_col_std)
                     write(unit_4,catalog_heading_format) (catalog_titles(j),j=1,N_col_std)

                     counter=0
                     do while ( counter < target_point_number )
                        call random_number_centered_normed( catalog_line(N_col_cat_x) , radius_limit )
                        call random_number_centered_normed( catalog_line(N_col_cat_y) , radius_limit )
                        call random_number_centered_normed( catalog_line(N_col_cat_z) , radius_limit )
                        call SphCoord( catalog_line(N_col_cat_x) , catalog_line(N_col_cat_y) , catalog_line(N_col_cat_z) , &
                           catalog_line(N_col_cat_l) , catalog_line(N_col_cat_b) , catalog_line(N_col_cat_dl) )

                        if ( catalog_line(N_col_cat_dl) < radius_limit ) then
                           call make_luminosity( catalog_line(N_col_cat_Lum) )

                           catalog_line( N_col_cat_M   ) = Lum_to_Mag      ( catalog_line( N_col_cat_Lum ) )
                           catalog_line( N_col_cat_mag ) = &
                              abs_to_vis_mag  ( catalog_line( N_col_cat_M ) , catalog_line( N_col_cat_dl ) )
                           catalog_line( N_col_cat_rs ) = fun_from_R_to_z ( catalog_line( N_col_cat_dl ) )

                           write( unit_4 , catalog_format ) ( catalog_line(j) , j=1,N_col_std )
                           counter=counter+1

                           endif
                        enddo

                     close(unit_4)
               endif

               call generation_report ; write(*,*) '  datafile: ' , trim(name(path))

            write(*,'(/,A)') 'GUC is complited'
         end subroutine Generate_Uniform_Catalog

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

      subroutine Generate_Cyrle(norma)
         real(8) norma ; integer i

            unit_cyrle = random_unit()
            open(unit_cyrle,file=slashfix(cyrle_path),status='replace')

               N_cyrle = 1d2
               do i=1,N_cyrle
                  alpha = pi*2/N_cyrle*i
                  x = dcos(alpha)*norma
                  y = dsin(alpha)*norma
                  write(unit_cyrle,*) x,y    !=-  write(*,*) pi*2,N_cyrle*i,alpha,x,y,norma    ;pause
               end do

               close(unit_cyrle)
      end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

			subroutine catalog_loading( data_file_path , Sample )
				character(length) data_file_path ; integer ii
				real(8) Sample(:,:) ; Sample(:,:)=0d0

               unit_1 = random_unit()

               theformat='' ; theformat='(A'//trim(inttostr(length))//')'
					open(unit_1,file=data_file_path,status='old',err=11)  !=-

					ii=0
               do  													!=- size()
                  read(unit_1,theformat,end=22) line
						if (.not. symbol_search(line,'#')) then

                     ii=ii+1
                     read(line,*) Sample(ii,:)

                     end if
						enddo ; goto 22

					11 continue ; write(*,*) '   catalog_loading: file ',trim(data_file_path),' does not exist'
					22 continue ; close(unit_1)
				end subroutine catalog_loading

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

			subroutine read_catalog_string( file_unit , Sample_string )
				integer file_unit
				real(8) Sample_string(:) ; Sample_string(:)=0d0
               theformat='' ; theformat='(A'//trim(inttostr(length))//')'

                  if (file_end) file_end=.false.
               do  													!=- size()
                  read(file_unit,theformat,end=22,err=11) line

						if (.not. symbol_search(line,'#')) then
                     read(line,*,end=11,err=11) Sample_string(:)
                     exit ; endif

						enddo ; goto 11 ; 22 continue ; file_end=.true.
						11 continue
				end subroutine read_catalog_string

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

      character(length) function make_catalog_name()
         make_catalog_name=''
         make_catalog_name = trim( folders( folder_Samples ) ) // trim(catalog_prefix) // '-' // &
            trim(adjustl(inttostr(FDE_seed))) // '_' // trim(adjustl(inttostr(generations))) // '-' // &
            trim(adjustl(realtostrf(Fractal_Dimension,4,2))) // '-' // trim(adjustl(inttostr(target_point_number))) // '.dat'
         end function make_catalog_name

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

      logical function count_selection( target_point_number , counter )
         integer target_point_number,counter
         real(8) p
            count_selection = .false.
            if (target_point_number<counter) then
               call random_number(p)
               if ( p > 1d0*target_point_number/counter ) count_selection = .true. !=- write(*,*) p , 1d0*target_point_number/counter , count_selection ; pause
               endif
         end function

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

      subroutine npoint(i,j,k,x,y,z)
         integer i,j,k
         real(8) x,y,z, nx,ny,nz
         intent (in) i,j,k
         intent (out) x,y,z

            if (npoint_noise<0) npoint_noise=1d0
            call random_number_centered_normed( nx , npoint_noise/fractionality )
            call random_number_centered_normed( ny , npoint_noise/fractionality )
            call random_number_centered_normed( nz , npoint_noise/fractionality )
            x = ((i-0.5)*1d0/fractionality-0.5d0)*2 + nx
            y = ((j-0.5)*1d0/fractionality-0.5d0)*2 + ny
            z = ((k-0.5)*1d0/fractionality-0.5d0)*2 + nz
         end subroutine npoint

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

end module



