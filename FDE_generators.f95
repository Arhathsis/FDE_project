!=- Fractal Dimension Estimation
!=- © Stanislav Shirokov, 2014-2020

module FDE_generators
   use FDE_paths
   use FDE_config

   use global
   use cosmology
   use math

   character(length) :: catalog_prefix = 'catalog' , secondary_catalog , skeleton , GCC_catalog_prefix

   logical ::  FDE_generator_report_flag = .true. , file_end = .false. , sphere_selection = .true. , &
               GCC_datafile = .true. , GSCS_replot = .true. , GCC_sets_regenerating , GCC_FDE_generator_report_flag, &
               GCC_sphere_selection

   real(8) ::  Lum, mag , catalog_line(N_col_std)=0d0 , chance , noise , npoint_noise , Fractal_Dimension , &
               nearest_point_to_origins(3) , min_dist , x_shift=0d0 , y_shift=0d0 , z_shift=0d0 , Super_radius_limit, &
               GCC_Super_Fractal_Dimension , GCC_Fractal_Dimension , GCC_z_shift , GCC_y_shift , GCC_x_shift , GCC_radius_limit

   integer ::  fractionality , gen_start , iteration , generations , uniform_approach_parameter , &
               step , imax , imin , jmin , jmax , kmin , kmax, GCC_unit=0, skeleton_generations, &
               unit_cyrle , N_cyrle , FDE_seed = 1 , FDE_seed_order = 1, GSCS_N, GCC_Super_generations, &
               GCC_Super_uniform_approach_parameter , GCC_Super_N_points , Super_FDE_seed_order , GCC_generations , &
               GCC_N_points , GCC_uniform_approach_parameter , GCC_FDE_seed_order , Super_counter

   contains



      subroutine make_FDE_seed
         real(8) x
            do
                  call random_number(x)
               FDE_seed = int( x*1d1**FDE_seed_order )
               if (FDE_seed < 1d1**FDE_seed_order .and. FDE_seed > 1d1**(FDE_seed_order-1)-1 ) exit
               enddo
      end subroutine



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



      subroutine make_luminosity( Luminosity )
         real(8) Luminosity

            Luminosity = log_uniform_luminosity()

            select case ( luminosity_model )

               case default

            end select
         end subroutine


      subroutine FDE_seed_fix
         real(8) x ; integer i
            do i=1,FDE_seed
               call random_number(x)
               end do
         end subroutine



      subroutine generation_report
            if ( FDE_generator_report_flag ) then
               write(*,'(/,3x,A,5x,F3.1)' ) 'D_eff:' , Fractal_Dimension
               write(*,'(3x,A,1x,i7    )' ) 'N_eff:' , N_points
               write(*,'(3x,A,2x,i6,A  )' ) 'R_max:' , int(radius_limit) , ' Mpc'  !=- !write(*,'(3x,A,1x,E10.3 )' ) 'Radius:'       , radius_limit
               !write(*,'(3x,A,1x,i7    )' ) 'N_p:  ' , final_counter
               endif
         end subroutine



   	subroutine Generate_Uniform_Catalog(path)
         character(length) path ; integer j

            write(*,'(/,A)') 'GUC is started..'

               inquire(file=slashfix(path), exist=file_exists)
            if ( file_exists .and. .not. sets_regenerating ) then
               !=- write(*,*) '   uniform set ',trim(path),' already exists'
               else

                  call FDE_seed_fix

                  unit_4 = random_unit()

                  if (N_points<1) N_points=100
                  open(unit_4,file=slashfix(path),status='replace')
                     write(unit_4,catalog_heads_format  ) (catalog_titles(j),j=1,N_col_std)
                     write(unit_4,catalog_heading_format) (catalog_titles(j),j=1,N_col_std)

                     counter=0
                     do while ( counter < N_points )
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

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1

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

						!	forall (i=1:N,j=1:3) Sample(i,j)=Sample(i,j) ; forall (i=1:N) Sample(i,6)=Sample(i,6)	- ???
					11 continue ; write(*,*) '   catalog_loading: file ',trim(data_file_path),' does not exist'
					22 continue ; close(unit_1)
				end subroutine catalog_loading



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



      character(length) function make_catalog_name()
         make_catalog_name=''
         make_catalog_name = trim( folders( folder_Samples ) ) // trim(catalog_prefix) // '-' // &
            trim(adjustl(inttostr(FDE_seed))) // '_' // trim(adjustl(inttostr(generations))) // '-' // &
            trim(adjustl(realtostrf(Fractal_Dimension,3,1))) // '-' // trim(adjustl(inttostr(N_points))) // '.dat'
         end function make_catalog_name

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1

      subroutine Generate_Super_Cantor_Catalog(filepath)
         character(length) filepath ; integer i , N_test
         real(8) XYZ(3), FD

            inquire(file=slashfix(filepath), exist=file_exists)

            if ( file_exists .and. .not. sets_regenerating ) then
               final_counter = file_volume(filepath)
               else

                     catalog_prefix = 'skeleton'
                     FDE_seed_order = 1
                  skeleton = trim( make_catalog_name() )

                     npoint_noise=0d0
                     GCC_datafile = .true.
                     FDE_generator_report_flag = .true.
                  call Generate_Cantor_Catalog(skeleton)
                     FDE_generator_report_flag = .false.
                     sets_regenerating = .true.
                     GCC_datafile = .false.
                  N_points = GCC_N_points/file_volume(skeleton) * 1.25  !=- test, why the factor is 1.25 ???

                     unit_4 = 1212
                     unit_6 = 1313
                  open(unit_4,file=skeleton,status='old')
                  open(unit_6,file=filepath,status='replace')
                     write(unit_6,catalog_heads_format  ) (catalog_titles(j),j=1,N_col_std)
                     write(unit_6,catalog_heading_format) (catalog_titles(j),j=1,N_col_std)
                     write(unit_6,'("# an initial count of points ",i7)') final_counter

                  radius_limit = ( GCC_radius_limit / 2**GCC_generations ) / ( 1d0 - min_dist )   !=- to need optimize
                  FDE_seed_order = 1
                  generations=GCC_Super_generations

                  GCC_unit = unit_6
                  npoint_noise = 0.98d0
                  sphere_selection = .false.
                  uniform_approach_parameter=GCC_Super_uniform_approach_parameter
                  Fractal_Dimension=GCC_Super_Fractal_Dimension

                  i=0 ; GSCS_N=0 ; N_test=0
                  do
                     call read_catalog_string( unit_4 , XYZ )
                        if (file_end) exit

                     i=1+i

                     x_shift=XYZ(1)
                     y_shift=XYZ(2)
                     z_shift=XYZ(3)

                        catalog_prefix = 'GSCC_' // trim(adjustl(inttostr(i)))
                        call make_FDE_seed
                     secondary_catalog = trim( make_catalog_name() )

                     call Generate_Cantor_Catalog(secondary_catalog)
     ! write(*,*) final_counter , N_points , 3*GCC_N_points/N_points ; pause

                           N_test=N_test+final_counter
                        if ( i == 20 ) then
                           !write(*,*) i,N_points , N_test/20
                           !N_points = N_points * N_points/N_test * 20
                           N_test=0
                        end if

                     enddo ; 11 continue
                     close(unit_4) ; close(unit_6)

                  final_counter = GSCS_N

                  call GCC_default

              end if
         end subroutine



      subroutine GCC_default

                       catalog_prefix = GCC_catalog_prefix              ;  if (GCC_catalog_prefix=='') catalog_prefix='catalog'

                    sets_regenerating = GCC_sets_regenerating
            FDE_generator_report_flag = .true.
                     sphere_selection = GCC_sphere_selection
                         GCC_datafile = .true.

                    Fractal_Dimension = GCC_Fractal_Dimension           ;  if ( GCC_Fractal_Dimension < 1 .or. &
                                                                              Fractal_Dimension>3 ) Fractal_Dimension = 3
                          generations = GCC_generations                 ;  if ( GCC_generations<1) generations=5
                             N_points = GCC_N_points                    ;  if ( GCC_N_points < 1 ) N_points=1d9     !=- corrected

                             GCC_unit = 0                               ;  if (GCC_unit<1) GCC_unit=0
                              z_shift = 0
                              y_shift = 0
                              x_shift = 0
                         radius_limit = GCC_radius_limit                ;  if (GCC_radius_limit<=0) GCC_radius_limit=1
           uniform_approach_parameter = GCC_uniform_approach_parameter  ;  if (GCC_uniform_approach_parameter<0) &
                                                                              uniform_approach_parameter=0
                       FDE_seed_order = GCC_FDE_seed_order              ;  if (GCC_FDE_seed_order<1) GCC_FDE_seed_order=1

         end subroutine



      subroutine Generate_Cantor_Catalog(path) !=- Cantor Sets Generator
         character(length) path ; integer i,j,k , ii,jj,kk
         real(8) x,y,z
         integer(1),allocatable,dimension(:,:,:) :: cube

            if (FDE_generator_report_flag) write(*,'(/,A)') 'GCC is started..'

               inquire(file=slashfix(path), exist=file_exists)
            if ( file_exists .and. .not. sets_regenerating ) then
               if (FDE_generator_report_flag) write(*,*) '   cantor set ',trim(path),' already exists'
               final_counter = file_volume(path)
               else

                  call FDE_seed_fix

                     if (GCC_datafile)  then
                        unit_3 = random_unit()
                        open(unit_3,file=slashfix(path),status='replace')
                           write(unit_3,catalog_heads_format  ) (catalog_titles(j),j=1,N_col_std)
                           write(unit_3,catalog_heading_format) (catalog_titles(j),j=1,N_col_std)
                        endif

                        fractionality=2**generations
                     allocate( cube(fractionality,fractionality,fractionality) )
                        forall (i=1:fractionality,j=1:fractionality,k=1:fractionality) cube(i,j,k)=1

                           if ( gen_start<1) gen_start=1
                        do i=1,gen_start ; call random_number(chance) ; enddo

                     do iteration=1,generations
                        if (iteration.gt.uniform_approach_parameter) then ; step=fractionality/2**iteration
                           do ii=1,2**iteration ; do jj=1,2**iteration ; do kk=1,2**iteration

                              imin=(ii-1)*step+1 ; imax=ii*step
                              jmin=(jj-1)*step+1 ; jmax=jj*step
                              kmin=(kk-1)*step+1 ; kmax=kk*step
                              call random_number(chance)

                              if (chance.gt.2**(Fractal_Dimension-3d0)) &
                                 forall (i=imin:imax,j=jmin:jmax,k=kmin:kmax) cube(i,j,k)=0

                              enddo ; enddo ; enddo

                           endif
                        enddo

                     min_dist=1d99
                     do ii=1,fractionality ; do jj=1,fractionality ; do kk=1,fractionality
                        if ( cube(ii,jj,kk)==1 ) then
                           x = ((ii-0.5)*1d0/fractionality-0.5d0)*2
                           y = ((jj-0.5)*1d0/fractionality-0.5d0)*2
                           z = ((kk-0.5)*1d0/fractionality-0.5d0)*2

                           if ( distance(x,y,z) < min_dist ) then
                              min_dist = distance(x,y,z)
                              nearest_point_to_origins(1)=x
                              nearest_point_to_origins(2)=y
                              nearest_point_to_origins(3)=z
                              end if

                           endif
                        enddo ; enddo ; enddo

                     counter=0
                     do ii=1,fractionality ; do jj=1,fractionality ; do kk=1,fractionality
                        !#write(*,*)'test'
                        if ( cube(ii,jj,kk)==1 ) then
                           !#write(*,*)'test2'
                           x = ((ii-0.5)*1d0/fractionality-0.5d0)*2
                           y = ((jj-0.5)*1d0/fractionality-0.5d0)*2
                           z = ((kk-0.5)*1d0/fractionality-0.5d0)*2

                           x = ( x - nearest_point_to_origins(1) + 1d-20 ) / ( 1d0 - min_dist )
                           y = ( y - nearest_point_to_origins(2) + 1d-20 ) / ( 1d0 - min_dist )
                           z = ( z - nearest_point_to_origins(3) + 1d-20 ) / ( 1d0 - min_dist )

                           if ( distance(x,y,z) > 1d0+1.5d0/fractionality .and. sphere_selection .or. &
                           dabs(x)>1 .or. dabs(y)>1 .or. dabs(z)>1) then
                              !#write(*,*)'test3'
                              cube(ii,jj,kk) = 0
                              else
                                 if ( dabs(x)<1 .and. dabs(y)<1 .and. dabs(z)<1 .and. .not. sphere_selection &
                                    .or. (GCC_sphere_selection .and. distance(x,y,z)<1 ) ) then
                                    !#write(*,*)'test4'
                                    counter=counter+1
                                 endif
                              endif
                           endif
                        enddo ; enddo ; enddo
                        if (GCC_datafile) write(unit_3,'("# an initial count of points ",i7)') counter

                     final_counter=0 ; Super_counter = 0
                     do ii=1,fractionality ; do jj=1,fractionality ; do kk=1,fractionality
                        if ( cube(ii,jj,kk)==1 ) then

                           call npoint( ii, jj, kk, &
                              catalog_line(N_col_cat_x), catalog_line(N_col_cat_y), catalog_line(N_col_cat_z) )

                              catalog_line(N_col_cat_x) = catalog_line(N_col_cat_x)*radius_limit + x_shift
                              catalog_line(N_col_cat_y) = catalog_line(N_col_cat_y)*radius_limit + y_shift
                              catalog_line(N_col_cat_z) = catalog_line(N_col_cat_z)*radius_limit + z_shift

                              call SphCoord( catalog_line(N_col_cat_x) , catalog_line(N_col_cat_y) , &
                                 catalog_line(N_col_cat_z) , &
                                 catalog_line(N_col_cat_l) , catalog_line(N_col_cat_b) , &
                                 catalog_line(N_col_cat_dl) )

                           if (  dabs(catalog_line(N_col_cat_x))<GCC_radius_limit &
                              .and. dabs(catalog_line(N_col_cat_y))<GCC_radius_limit &
                              .and. dabs(catalog_line(N_col_cat_z))<GCC_radius_limit &
                              .and. .not. sphere_selection  &
                              .or.  GCC_sphere_selection .and. catalog_line(N_col_cat_dl) < GCC_radius_limit  ) then

                              if ( .not. count_selection( N_points , counter ) ) then  !=-  .or. 5>3
                                 call make_luminosity( catalog_line(N_col_cat_Lum) )

                                 catalog_line( N_col_cat_M   ) = Lum_to_Mag      ( catalog_line( N_col_cat_Lum ) )
                                 catalog_line( N_col_cat_mag ) = &
                                    abs_to_vis_mag  ( catalog_line( N_col_cat_M ) , catalog_line( N_col_cat_dl ) )
                                 catalog_line( N_col_cat_rs ) = fun_from_R_to_z ( catalog_line( N_col_cat_dl ) )

                                 final_counter=final_counter+1
                                 if (GCC_datafile) write( unit_3 , catalog_format ) ( catalog_line(j) , j=1,N_col_std )

                                 if ( GCC_unit/=0 ) then
                                    write( GCC_unit , catalog_format ) ( catalog_line(j) , j=1,N_col_std )
                                    GSCS_N=1+GSCS_N
                                    Super_counter = 1 + Super_counter
                                    endif
                                 endif
                              endif
                           end if
                        enddo ; enddo ; enddo !=- write(*,*) final_counter , Super_counter , counter , N_points ; pause

                     if (final_counter==0 .and. FDE_generator_report_flag) write(*,*) 'GCC: 0 points' !=- need error fix

                     if (GCC_datafile) close(unit_3)
                     deallocate(cube)
               endif

               call generation_report ; if (FDE_generator_report_flag) write(*,*) '  datafile: ' , trim(name(path))

            if (FDE_generator_report_flag) write(*,'(/,A)') 'GCC is complited'
         end subroutine Generate_Cantor_Catalog


      logical function count_selection( N_points , counter )
         integer N_points,counter
         real(8) p
            count_selection = .false.
            if (N_points<counter) then
               call random_number(p)
               if ( p > 1d0*N_points/counter ) count_selection = .true. !=- write(*,*) p , 1d0*N_points/counter , count_selection ; pause
               endif
         end function

      subroutine npoint(i,j,k,x,y,z) !===================================================================1
         integer i,j,k
         real(8) x,y,z, nx,ny,nz
         intent (in) i,j,k
         intent (out) x,y,z

            if (npoint_noise<0) npoint_noise=0.98
            call random_number_centered_normed( nx , npoint_noise/fractionality )
            call random_number_centered_normed( ny , npoint_noise/fractionality )
            call random_number_centered_normed( nz , npoint_noise/fractionality )
            x=(((i-0.5)*1d0/fractionality-0.5d0)*2 + nx - nearest_point_to_origins(1) + 1d-20) / ( 1d0 - min_dist )
            y=(((j-0.5)*1d0/fractionality-0.5d0)*2 + ny - nearest_point_to_origins(2) + 1d-20) / ( 1d0 - min_dist )
            z=(((k-0.5)*1d0/fractionality-0.5d0)*2 + nz - nearest_point_to_origins(3) + 1d-20) / ( 1d0 - min_dist )
         end subroutine npoint !============================================================================1



end module



