module FDE_generators
   use FDE_paths
   use FDE_config
   use global

   use cosmology
   use math

   real(8) ::  Lum, mag , catalog_line(N_col_std)=0d0 , chance , noise , Fractal_Dimension , &
               nearest_point_to_origins(3) , min_dist

   integer ::  fractionality , gen_start , iteration , generations , uniform_approach_parameter , &
               step , imax , imin , jmin , jmax , kmin , kmax, &
               N_cyrle , FDE_seed = 1 , FDE_seed_order = 1

   contains



      subroutine make_FDE_seed
         real(8) x
            do
                  call random_number(x)
               FDE_seed = int( x*1d1**FDE_seed_order )
               if (FDE_seed /= 1d1**FDE_seed_order .and. FDE_seed /= 0 ) exit
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



   	subroutine Generate_Uniform_Catalog(path)
         character(len) path

            write(*,*) 'GUC is started..'
               inquire(file=path, exist=file_exists)

            if ( file_exists .and. .not. sets_regenerating ) then
               write(*,*) '   uniform set ',trim(path),' already exists'
               else

                  call FDE_seed_fix

                  if (N_points<1) N_points=100
                  open(N_points,file=path,status='replace')
                     write(N_points,catalog_heads_format  ) (catalog_titles(j),j=1,N_col_std)
                     write(N_points,catalog_heading_format) (catalog_titles(j),j=1,N_col_std)

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

                           write( N_points , catalog_format ) ( catalog_line(j) , j=1,N_col_std )
                           counter=counter+1

                           endif
                        enddo

                     close(N_points)
               endif

            write(*,*) 'GUC is complited'
         end subroutine Generate_Uniform_Catalog



      subroutine Generate_Cyrle(norma)
         real(8) norma
            N_cyrle = 101
            open(N_cyrle,file=cyrle_path,status='replace')

               do i=1,N_cyrle
                  alpha = pi*2/N_cyrle*i
                  x = dcos(alpha)*norma
                  y = dsin(alpha)*norma
                  write(N_cyrle,*) x,y    !=-  write(*,*) pi*2,N_cyrle*i,alpha,x,y,norma    ;pause
               end do

               close(N_cyrle)
      end subroutine



      subroutine Generate_Cantor_Catalog(path) !=- Cantor Sets Generator
         character(len) path
         integer(1),allocatable,dimension(:,:,:) :: cube

            write(*,*) 'GCC is started..'
               inquire(file=path, exist=file_exists)

            if ( file_exists .and. .not. sets_regenerating ) then
               write(*,*) '   cantor set ',trim(path),' already exists'
               else

                  call FDE_seed_fix

                     if (N_points<1) N_points=100
                  open(N_points,file=path,status='replace')
                     write(N_points,catalog_heads_format  ) (catalog_titles(j),j=1,N_col_std)
                     write(N_points,catalog_heading_format) (catalog_titles(j),j=1,N_col_std)

                        if ( Fractal_Dimension < 1 ) Fractal_Dimension = 3
                        if ( uniform_approach_parameter<=0) uniform_approach_parameter=1
                        if ( generations<1) generations=5
                        if ( gen_start<1) gen_start=1

                     fractionality=2**generations
                     allocate( cube(fractionality,fractionality,fractionality) )
                     forall (i=1:fractionality,j=1:fractionality,k=1:fractionality) cube(i,j,k)=1

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
                        if ( cube(ii,jj,kk)==1 ) then
                           x = ((ii-0.5)*1d0/fractionality-0.5d0)*2
                           y = ((jj-0.5)*1d0/fractionality-0.5d0)*2
                           z = ((kk-0.5)*1d0/fractionality-0.5d0)*2

                           x = ( x - nearest_point_to_origins(1) + 1d-5 ) / ( 1d0 - min_dist )
                           y = ( y - nearest_point_to_origins(2) + 1d-5 ) / ( 1d0 - min_dist )
                           z = ( z - nearest_point_to_origins(3) + 1d-5 ) / ( 1d0 - min_dist )

                           if ( distance(x,y,z) > 1 ) then
                              cube(ii,jj,kk) = 0
                              else
                                 counter=counter+1
                              endif
                           endif
                        enddo ; enddo ; enddo
                        write(N_points,'("# an initial count of points ",i7)') counter

                     final_counter=0
                     do ii=1,fractionality ; do jj=1,fractionality ; do kk=1,fractionality
                        if ( cube(ii,jj,kk)==1 ) then

                           call npoint( ii, jj, kk, &
                              catalog_line(N_col_cat_x),catalog_line(N_col_cat_y),catalog_line(N_col_cat_z) )
                              catalog_line(N_col_cat_x) = catalog_line(N_col_cat_x)*radius_limit
                              catalog_line(N_col_cat_y) = catalog_line(N_col_cat_y)*radius_limit
                              catalog_line(N_col_cat_z) = catalog_line(N_col_cat_z)*radius_limit

                           call SphCoord( catalog_line(N_col_cat_x) , catalog_line(N_col_cat_y) , &
                              catalog_line(N_col_cat_z) , &
                              catalog_line(N_col_cat_l) , catalog_line(N_col_cat_b) , &
                              catalog_line(N_col_cat_dl) )

                           if ( .not. count_selection( N_points , counter ) ) then
                              call make_luminosity( catalog_line(N_col_cat_Lum) )

                              catalog_line( N_col_cat_M   ) = Lum_to_Mag      ( catalog_line( N_col_cat_Lum ) )
                              catalog_line( N_col_cat_mag ) = &
                                 abs_to_vis_mag  ( catalog_line( N_col_cat_M ) , catalog_line( N_col_cat_dl ) )
                              catalog_line( N_col_cat_rs ) = fun_from_R_to_z ( catalog_line( N_col_cat_dl ) )

                              final_counter=final_counter+1
                              write( N_points , catalog_format ) ( catalog_line(j) , j=1,N_col_std )

                              endif
                           end if
                        enddo ; enddo ; enddo

                     if (final_counter==0) write(*,*) 'GCC: 0 points' !=- need error fix

                     close(N_points)
                     deallocate(cube)
               endif

            write(*,*) 'GCC is complited'
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

            if (noise<=0) noise=0.98
            call random_number_centered_normed( nx , noise/fractionality )
            call random_number_centered_normed( ny , noise/fractionality )
            call random_number_centered_normed( nz , noise/fractionality )
            x=(((i-0.5)*1d0/fractionality-0.5d0)*2 + nx - nearest_point_to_origins(1) + 1d-5) / ( 1d0 - min_dist )
            y=(((j-0.5)*1d0/fractionality-0.5d0)*2 + ny - nearest_point_to_origins(2) + 1d-5) / ( 1d0 - min_dist )
            z=(((k-0.5)*1d0/fractionality-0.5d0)*2 + nz - nearest_point_to_origins(3) + 1d-5) / ( 1d0 - min_dist )
         end subroutine npoint !============================================================================1



end module



