!=- Fractal Dimension Estimation
!=- © Stanislav Shirokov, 2014-2020

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1

	module FDE_methods
		use global
		use math

		use FDE_config
		use FDE_paths

         logical :: FDE_data_files_existing = .false.

         integer,parameter :: methods_count = 4 , geometries_count = 4 , FDE_out_table_size = 100

         character(len) :: NP_name = '_NP.dat' , MD_name = '_MD.dat' , diffCD_name = '_DiffCD.dat' , &
                           intCD_name = '_IntCD.dat' , FDE_data_file_path , FDE_data_format , &
                           FDE_NP_data_file_path , FDE_MD_data_file_path , FDE_IntCD_data_file_path , &
                           FDE_DiffCD_data_file_path , FDE_data_files_paths(methods_count)

			integer ::	FDE_points_amount , FDE_grid = 20 , k_NP(geometries_count) , geometry

			real(8) ::	FDE_R_min , FDE_R_max , FDE_l_min , FDE_l_max , FDE_b_min , FDE_b_max , FDE_D_min , FDE_D_max , &
							logS , FDE_distance , R_nn(geometries_count) , FDE_LS_coefficients(geometries_count,2) , &
							FDE_left_border , FDE_right_border , temp_XY(2,FDE_out_table_size) , &

							FDE_out_NP      ( 1 + 2*geometries_count , FDE_out_table_size ) , &
							FDE_out_MD      ( 1 + 2*geometries_count , FDE_out_table_size ) , &
							FDE_out_intCD   ( 1 + 2*geometries_count , FDE_out_table_size ) , &
							FDE_out_diffCD  ( 1 + 2*geometries_count , FDE_out_table_size ) , &

							XYY             ( 1 + 2*geometries_count,FDE_out_table_size )

         real(8),allocatable,dimension(:,:) ::  Sample , Integral_Density , Differential_Density , Mutual_Distances , &
                                                Dispersion_NP , Dispersion_MD , Dispersion_intCD , Dispersion_diffCD, &
                                                Nearest_Neighbores

         real(8),allocatable,dimension(:) :: Radii , Volumes

         integer(8),allocatable,dimension(:,:) :: working_volumes , FDE_diffCD , FDE_MD , FDE_NP



			contains

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1

			subroutine catalog_loading( data_file_path , Sample )
				character(len) data_file_path
				real(8) Sample(:,:) ; Sample(:,:)=0d0

               theformat='' ; theformat='(A'//trim(inttostr(len))//')'
					open(unit_1,file=trim(data_file_path),status='old',err=11)  !=-

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

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1

			subroutine geometries_making
            integer i

					if (FDE_R_min<=0) FDE_R_min = radius_limit/1d3
					if (FDE_R_max<=0) FDE_R_max = radius_limit
					FDE_D_min = 2*FDE_R_min
					FDE_D_max = 2*FDE_R_max
					if (FDE_b_min<=0) FDE_b_min = latitude_limit    !=- "+" -> 1 < 0 < 1 , "-" -> 0 < 1 < 0
               if (FDE_b_max<=0) FDE_b_max = -latitude_limit   !=- "-" -> 1 < 0 < 1 , "+" -> 0 < 1 < 0

					logS=(log10(FDE_R_max)-log10(FDE_R_min))/FDE_grid

					forall (i=1:FDE_grid) Radii(i)=1d1**(log10(FDE_R_min)+(i-0.5)*logS)
					forall (i=1:FDE_grid) Volumes(i)=4d0/3d0*pi*Radii(i)**3 ; Volumes(0)=0d0

				end subroutine geometries_making



         real(8) function vector_distance(A,B)
            real(8) A(:),B(:)
               vector_distance = ( (A(1)-B(1))**2 + (A(2)-B(2))**2 + (A(3)-B(3))**2 )**0.5d0
            end function



         real(8) function compute_alpha( point_b )
            real(8) point_b
               if (point_b>=0) compute_alpha = ( point_b   - FDE_b_max ) / 180d0*pi
               if (point_b< 0) compute_alpha = ( FDE_b_min - point_b   ) / 180d0*pi
            end function compute_alpha

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1

         subroutine check_FDE_data_files_existing
            FDE_data_files_existing = .true.
            do i=1,methods_count
               inquire( file = FDE_data_files_paths(i) , exist = file_exists )
               FDE_data_files_existing = FDE_data_files_existing .and. file_exists
               enddo
         end subroutine check_FDE_data_files_existing



			subroutine FDE_complex( data_file_path )
				character(len) data_file_path ; integer i , j , k , ii , jj , kk

					write(*,*) 'FDE-method is started..'

               FDE_data_file_path = name(data_file_path)
               FDE_data_format = '(' // trim(inttostr(2*geometries_count+1)) // '(E16.8))'
               if (FDE_grid<2) FDE_grid=20
                  call FDE_create_names

					FDE_points_amount = file_volume( data_file_path )
                  if (FDE_points_amount<2) write(*,*) 'FDE_complex: no data'

               call check_FDE_data_files_existing
            if ( FDE_recalculating .or. .not. FDE_data_files_existing ) then

					allocate( Sample(FDE_points_amount,N_col_std) ) !=- 7
					call catalog_loading( data_file_path , Sample )

					allocate( 	Dispersion_NP        ( geometries_count , FDE_grid ) , &
                           Dispersion_MD        ( geometries_count , FDE_grid ) , &
                           Dispersion_intCD     ( geometries_count , FDE_grid ) , &
                           Dispersion_diffCD    ( geometries_count , FDE_grid ) , &

                           FDE_NP               ( geometries_count , FDE_grid ) , &
                           FDE_MD               ( geometries_count , FDE_grid ) , &
                           FDE_diffCD	         ( geometries_count , FDE_grid ) , &
									working_volumes      ( geometries_count , FDE_grid ) , &

                           Nearest_Neighbores   ( geometries_count , FDE_grid ) , &
									Mutual_Distances     ( geometries_count , FDE_grid ) , &
									Integral_Density     ( geometries_count , FDE_grid ) , &
									Differential_Density ( geometries_count , FDE_grid ) , &

									Radii			(  FDE_grid) , &
									Volumes		(0:FDE_grid)  )

					call geometries_making

               Dispersion_NP        (:,:) = 0d0
               Dispersion_MD        (:,:) = 0d0
               Dispersion_intCD     (:,:) = 0d0
               Dispersion_diffCD    (:,:) = 0d0

               FDE_MD	            (:,:) = 0
               FDE_NP               (:,:) = 0
					FDE_diffCD	         (:,:) = 0
					working_volumes      (:,:) = 0

               Nearest_Neighbores   (:,:) = 0d0
               Mutual_Distances     (:,:) = 0d0
					Integral_Density     (:,:) = 0d0
					Differential_Density (:,:) = 0d0

					do i=1,FDE_points_amount
                  k_NP(:) = 1 + FDE_grid
                  FDE_diffCD(:,:)=0d0

						do j=1,FDE_points_amount
                     FDE_distance = vector_distance( Sample(i,:) , Sample(j,:) )

							if ( i.ne.j .and. Sample(i,6).le.FDE_R_max .and. Sample(j,6).le.FDE_R_max ) then

                           if (FDE_distance.gt.FDE_D_min) &
                              k = 1 + dint( FDE_grid*log10(FDE_distance/FDE_D_min)/(log10(FDE_D_max/FDE_D_min)) )

                        FDE_MD(1,k) = 1 + FDE_MD(1,k)
									if (Sample(i,5).le.FDE_b_max .and. Sample(j,5).le.FDE_b_max) &
                              FDE_MD(2,k) = 1 + FDE_MD(2,k)
									if (Sample(i,5).le.FDE_b_min .and. Sample(j,5).le.FDE_b_min) &
                              FDE_MD(3,k) = 1 + FDE_MD(3,k)

                        if ( FDE_distance.ge.FDE_R_min .and. FDE_distance.le.FDE_R_max ) then

                           k = 1 + dint( FDE_grid*log10(FDE_distance/FDE_R_min)/(log10(FDE_R_max/FDE_R_min)) )
                              if ( k.lt.k_NP(1) ) k_NP(1) = k  !=- a Nearest Point

                           if ( Radii(k).le.FDE_R_max-Sample(i,6) ) then

                              FDE_diffCD(1,k) = 1 + FDE_diffCD(1,k)

                              if ( Radii(k).le.dsin( compute_alpha( Sample(i,5) ))*Sample(i,6) ) then  !=- the point within geometry

                                 if (Sample(i,5)<=0) then
                                    FDE_diffCD(2,k) = 1 + FDE_diffCD(2,k)
                                       if ( k.lt.k_NP(2) ) k_NP(2) = k  !=- a Nearest Point
                                    endif
                                 if (Sample(i,5)> 0) then
                                    FDE_diffCD(3,k) = 1 + FDE_diffCD(3,k)
                                       if ( k.lt.k_NP(3) ) k_NP(3) = k  !=- a Nearest Point
                                    endif

                                 endif
                              endif
                           endif
                        endif
							enddo

                  call compute_working_volumes
                  call compute_concentration

						call countNP

						enddo

               Mutual_Distances  (:,:) = FDE_MD(:,:)
               Nearest_Neighbores(:,:) = FDE_NP(:,:)

               working_volumes(4,:) = working_volumes(2,:) + working_volumes(3,:)

            !=- the conditional density calculations
               forall ( jj=1:3 , ii=1:FDE_grid , working_volumes(jj,ii)>0 )
                  Integral_Density (jj,ii) = &
                     Integral_Density (jj,ii) / working_volumes(jj,ii) !/ Volumes(ii)
                  Differential_Density (jj,ii) = &
                     Differential_Density (jj,ii) / working_volumes(jj,ii) !/ ( Volumes(ii) - Volumes(ii-1) )
                  end forall
            !=- end of the conditional density calculations

            !=- the forth geometry calculations
               forall (ii=1:FDE_grid) Mutual_Distances (4,ii) = FDE_MD(2,ii) + FDE_MD(3,ii)
               forall (ii=1:FDE_grid) Integral_Density (4,ii) = &
                  ( Integral_Density  (2,ii) + Integral_Density  (3,ii) ) * 0.5d0
					forall (ii=1:FDE_grid) Differential_Density (4,ii) = &
                  ( Differential_Density (2,ii) + Differential_Density (3,ii) ) * 0.5d0
            !=- end of the forth geometry calculations

               call FDE_compute_dispersion

               call writeNP ; call writeMD ; call writeCD

				deallocate(  Sample , Radii , Volumes , &
                Nearest_Neighbores , Mutual_Distances , Differential_Density , Integral_Density , &
                Dispersion_NP , Dispersion_MD , Dispersion_intCD , Dispersion_diffCD )
            deallocate( FDE_NP , FDE_MD , FDE_diffCD , working_volumes ) !=- Exception: Access Violation
               else
                  write(*,*) '   FDE_complex: all data files already exist'
            endif

            call FDE_export

				write(*,*) 'FDE-method is complited'

				end subroutine FDE_complex

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1

         subroutine FDE_compute_dispersion
            integer i , j
            !=- the Poisson errors
            forall (j=1:geometries_count , i=1:FDE_grid , Nearest_Neighbores(j,i)>0 )
               Dispersion_NP(j,i)     = Nearest_Neighbores(j,i)   **2 * &
                  ( 1d0 / FDE_points_amount + 1d0 / Nearest_Neighbores(j,i) )
               end forall

            forall (j=1:geometries_count , i=1:FDE_grid , Mutual_Distances(j,i)>0 )
               Dispersion_MD(j,i)     = Mutual_Distances(j,i)     **2 * &
                  ( 1d0 / FDE_points_amount + 1d0 / Mutual_Distances(j,i) )
               end forall

            forall (j=1:geometries_count , i=1:FDE_grid , working_volumes(j,i)>0 )
               Dispersion_intCD(j,i)  = Integral_Density(j,i)     **2 * &
                  ( 1d0 / FDE_points_amount + 1d0 / sum(working_volumes(j,1:i)) )
               Dispersion_diffCD(j,i) = Differential_Density(j,i) **2 * &
                  ( 1d0 / FDE_points_amount + 1d0 / working_volumes(j,i) )
               end forall
            !=- end of the Poisson errors
            end subroutine FDE_compute_dispersion


         subroutine FDE_export

               call FDE_read_R_nn

               FDE_out_NP(:,:)=0d0 ; FDE_out_MD(:,:)=0d0 ; FDE_out_intCD(:,:)=0d0 ; FDE_out_diffCD(:,:)=0d0
               call read_FDE_data_file(FDE_data_files_paths(1))
                  FDE_out_NP      (:,:)=XYY(:,:)
               call read_FDE_data_file(FDE_data_files_paths(2))
                  FDE_out_MD      (:,:)=XYY(:,:)
               call read_FDE_data_file(FDE_data_files_paths(3))
                  FDE_out_intCD   (:,:)=XYY(:,:)
               call read_FDE_data_file(FDE_data_files_paths(4))
                  FDE_out_diffCD  (:,:)=XYY(:,:)
            end subroutine FDE_export



         subroutine FDE_create_names
            FDE_NP_data_file_path     = trim(folders(folder_NP_add_files)) // trim(FDE_data_file_path) // &
               '_' // trim(adjustl(inttostr(FDE_grid))) // trim(NP_name)
            FDE_MD_data_file_path     = trim(folders(folder_MD_add_files)) // trim(FDE_data_file_path) // &
               '_' // trim(adjustl(inttostr(FDE_grid))) // trim(MD_name)
            FDE_IntCD_data_file_path  = trim(folders(folder_CD_add_files)) // trim(FDE_data_file_path) // &
               '_' // trim(adjustl(inttostr(FDE_grid))) // trim(intCD_name)
            FDE_DiffCD_data_file_path = trim(folders(folder_CD_add_files)) // trim(FDE_data_file_path) // &
               '_' // trim(adjustl(inttostr(FDE_grid))) // trim(diffCD_name)

               FDE_data_files_paths(1) = FDE_NP_data_file_path
               FDE_data_files_paths(2) = FDE_MD_data_file_path
               FDE_data_files_paths(3) = FDE_IntCD_data_file_path
               FDE_data_files_paths(4) = FDE_DiffCD_data_file_path
            end subroutine FDE_create_names



         subroutine FDE_read_R_nn
            integer i , j
            R_nn(:)=0d0
               call read_FDE_data_file(FDE_data_files_paths(1))
            do j=1,geometries_count !=- the maximum searching (FDE_v.2.8 contains the half-integral maximum searching)
               maximum = maxval(XYY(2*j,:))

               do i=1,FDE_grid

                  if ( XYY(2*j,i) == maximum ) R_nn(j) = XYY(1,i)   !=- R_nn ~ (V / N)^{1/3}

                  end do
               end do             !=-  write(*,*) (pi*4/3*radius_limit**3/N_points)**0.33, R_nn
               FDE_out_NP(:,:) = XYY(:,:)
            end subroutine FDE_read_R_nn



         subroutine countNP
            integer i

            k_NP(4) = min(k_NP(2),k_NP(3))

            forall ( i=1:geometries_count , k_NP(i)<=FDE_grid ) FDE_NP(i,k_NP(i)) = 1 + FDE_NP(i,k_NP(i))

            end subroutine countNP



         subroutine compute_working_volumes
            integer i , k
               do i=1,3
                  do k=1,size(FDE_diffCD(i,:))
                        if ( FDE_diffCD(i,k)>0 ) working_volumes(i,k) = 1 + working_volumes(i,k)
                     enddo
                  enddo
            end subroutine compute_working_volumes



         subroutine compute_concentration
            integer i , k
               do i=1,3 ; do k=1,size(FDE_diffCD(i,:)) ; if ( FDE_diffCD(i,k)>0 ) then

                  Integral_Density(i,k)     = Integral_Density(i,k)      + sum(FDE_diffCD(1:i,k))  !+mN_ind(k) ???
                  Differential_Density(i,k) = Differential_Density(i,k)  + FDE_diffCD(i,k)         !+mN_ind(k) ???

                  endif ; enddo ; enddo
            end subroutine compute_concentration



         subroutine writeNP
            integer i

            open(unit_1,file=FDE_NP_data_file_path,status='replace')

               do i=1,FDE_grid
                 if ( Nearest_Neighbores(1,i)>0) &
                     write(unit_1,FDE_data_format) Radii(i), &
                        Nearest_Neighbores(1,i), Dispersion_NP(1,i)**0.5, &
                        Nearest_Neighbores(2,i), Dispersion_NP(2,i)**0.5, &
                        Nearest_Neighbores(3,i), Dispersion_NP(3,i)**0.5, &
                        Nearest_Neighbores(4,i), Dispersion_NP(4,i)**0.5
                  enddo

               close(unit_1)
            end subroutine writeNP



         subroutine writeMD
            open(unit_2,file=FDE_MD_data_file_path,status='replace')

               do i=1,FDE_grid
                  if (Mutual_Distances(1,i)>0)  &
                     write(unit_2,FDE_data_format) 2*Radii(i),   &
                        Mutual_Distances(1,i), Dispersion_MD(1,i)**0.5, &
                        Mutual_Distances(2,i), Dispersion_MD(2,i)**0.5, &
                        Mutual_Distances(3,i), Dispersion_MD(3,i)**0.5, &
                        Mutual_Distances(4,i), Dispersion_MD(4,i)**0.5
                 enddo

               close(unit_2)
            end subroutine writeMD



         subroutine writeCD
            open(unit_3,file=FDE_IntCD_data_file_path,status='replace')

               do i=1,FDE_grid
               if (Integral_Density(1,i)>0)  &
                  write(unit_3,FDE_data_format) Radii(i),   &
                     Integral_Density(1,i), Dispersion_intCD(1,i)**0.5, &
                     Integral_Density(2,i), Dispersion_intCD(2,i)**0.5, &
                     Integral_Density(3,i), Dispersion_intCD(3,i)**0.5, &
                     Integral_Density(4,i), Dispersion_intCD(4,i)**0.5
                  enddo

               close(unit_3)

            open(unit_4,file=FDE_DiffCD_data_file_path,status='replace')

               do i=1,FDE_grid
               if (Differential_Density(1,i)>0)  &
                  write(unit_4,FDE_data_format) Radii(i),   &
                     Differential_Density(1,i), Dispersion_diffCD(1,i)**0.5, &
                     Differential_Density(2,i), Dispersion_diffCD(2,i)**0.5, &
                     Differential_Density(3,i), Dispersion_diffCD(3,i)**0.5, &
                     Differential_Density(4,i), Dispersion_diffCD(4,i)**0.5
                  enddo

               close(unit_4)
            end subroutine writeCD



         subroutine read_FDE_data_file(data_file_path)
            character(len) data_file_path
               open(unit_7,file=data_file_path,status='old',err=22)

                  XYY(:,:) = 0d0
                  read(unit_7,*,end=11) XYY

               22 continue ; write(*,*) 'no file'
               11 continue ; close(unit_7)
            end subroutine read_FDE_data_file



			subroutine FDE_trend_log_xy( data_file_path )
            character(len) data_file_path ; integer i , j , k

            call FDE_read_R_nn

            call read_FDE_data_file(data_file_path)

				temp_XY(:,:)=0d0 ; FDE_LS_coefficients(:,:)=0d0
				forall (j=1:FDE_out_table_size , XYY(1,j)>0 ) temp_XY(1,j)=log10(XYY(1,j))
				do i=1,geometries_count
               forall (j=1:FDE_out_table_size , XYY(2*i,j)>0 ) temp_XY(2,j)=log10( XYY(2*i,j) )

					k=0;sumx=0d0;sumx2=0d0;sumy=0d0;sumxy=0d0

					do j=1,FDE_out_table_size
                  if ( temp_XY(2,j) .ne. 0d0 .and. &
                     temp_XY(1,j).ge.log10(FDE_left_border) .and. temp_XY(1,j).le.log10(FDE_right_border) ) then   !		write(*,*) N,G(1,j),G(2,j),XY(1,j),XY(2,j)

                     sumx=sumx+temp_XY(1,j) ; sumx2=sumx2+temp_XY(1,j)**2d0
                     sumy=sumy+temp_XY(2,j) ; sumxy=sumxy+temp_XY(1,j)*temp_XY(2,j)
                     k=k+1

                     endif
                  enddo

					if (k.ne.0) then
                  FDE_LS_coefficients(i,1)=( k*sumxy-sumx*sumy ) / ( k*sumx2-sumx**2 )
                  FDE_LS_coefficients(i,2)=( sumy - FDE_LS_coefficients(i,1)*sumx ) / k
                  endif

					enddo
				end subroutine



			subroutine FDE_variation_log_xy( data_file_path )   !=- have not been written
            character(len) data_file_path

				end subroutine

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1

   end module
