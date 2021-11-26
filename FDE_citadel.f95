!=- Fractal Dimension Estimation
!=- © Stanislav Shirokov, 2014-2020

module FDE_citadel
   use FDE_config
   use FDE_paths
   use FDE_scripts
   use FDE_generators

   use global

   character(length) input_command , logscaling
   real(8) :: start_scale = -1, final_scale = -1, scale_step = -1, mean_MD_D = 2.0

   integer :: mean_MD_N = 10

   contains

		subroutine preparation

				write(*,*) 'Fractal Dimension Estimation is runned' ; write(*,*)

            call define_system
            call gen_paths
            call shell_MD_Tree
            call the_catalog_heading

			end subroutine

      subroutine opening !===============================================================================1

            22 continue

				call preparation

				do ; 33 continue ; write(*,*) '--'

					read(*,'(A256)',err=33,end=33) input_command
                  read(input_command,*,err=33,end=33) command

					select case (command)

                  case ('reset')
                     goto 22

                  case ('scaling')
                     start_scale=-1 ; final_scale=-1 ; scale_step=-1 ; logscaling=''
                     read(input_command,*,iostat=iostat_value) command, start_scale, final_scale, scale_step , logscaling
                        if ( start_scale==-1 .or. final_scale==-1 .or. scale_step==-1 ) then
                           write(*,*) 'scaling input error: scaling 10 100 10 [log]'
                           else
                              call scaling(start_scale, final_scale, scale_step , logscaling)
                           endif

						case ('m')
							call means(2d0, 5000, 3)
                        case ('mm')
							call means_matrix
                        case ('mmt')
                            call create_table_of_means(2d0, 2d0, 1, 5000, 15000, 5, 5, 20, 4)
						case ('exit')
							goto 11
						case ('help')
							call help
						case ('uca')
							call uniform_catalog_analysis
						case ('cca')
							call cantor_catalog_analysis

                  case ('ex1')
                     call example1
                  case ('meanmd0')
                     call mean_MD(2.0d0, 3)  !=- subroutine's input parameters should be as variables (see below)
                  case ('meanmd')
                     read(input_command,*,iostat=iostat_value) command,mean_MD_D,mean_MD_N
                     call mean_MD( mean_MD_D, mean_MD_N )  !=- the input command's format: 'meanmd 2 3'

						case ('c')
							call input_catalogs_analysis

                  case ('os')
                     read(input_command,*,iostat=iostat_value) command,operating_system
                     text1='Windows' ; if (operating_system==2) text1='Linux'
                     write(*,*) 'The operating system has been changed to ', trim(text1)

                  case ('avz')
                     read(input_command,*,iostat=iostat_value) command,command1,command2
                     if ( command1=='universe' .or. command1=='u' ) call add_visible_z_universe

                  case ('add')
                     read(input_command,*,iostat=iostat_value) command,command1,command2
                     if ( command1=='visible' .and. command2=='z') call add_visible_z

                     call add_visible_z   !=- fix
                     call plot_catalog ( FDE_catalog ) !=- fix

                  case ('e')
                     call analysis_external_FDE_catalog

						case (char(225))	!=- rus c
							call input_catalogs_analysis

						case ('-r c')
							do
								!call catalogue ; pause
								enddo
						case ('a')
							!call analytical
						case ('md')
							call shell_MD_Tree
						case ('s')
							!call samples
						case ('rs')
							!call redshift_distance
						case ('n')
							read(*,*) n
							do i=1,n
								!call samples
								enddo
						case ('z')
							!call redshift_test
						case ('d')
							call system('>> log.txt (DEL /S/Q CD MD NP)')
						case ('dd')
							call system('>> log.txt (DEL /S/Q CD MD NP Samples Statistics)')
						case ('ds')
							call system('>> log.txt (RD /S/Q Samples && MD Samples)')
						case default
							write(*,*) '   enter the command..'

					end select ; enddo ; 11 continue ; write(*,*)
			end subroutine opening !===========================================================================1



		subroutine help !==================================================================================2
			implicit none ; character(10) text ; write(*,*)

            write(*,*) '   analysis     - main FDE sqript'
            write(*,*) '   samples      - VL samples sqript'
            write(*,*) '   readMDNP(md) - reading middle distance between points from last data'
            write(*,*)
            write(*,*) '   CSG(rdm,fractal)    - Cantor Sets Generator'
            write(*,*) '   SM(set,list(SA,SS,SN,SD)) - Samples Maker'
            write(*,*)
            write(*,*) '   MDNPM(set)          - Minimal Distance to Nearest Point'
            write(*,*) '   MDM(set,Dmd)        - Mutual Distance'
            write(*,*) '   CDM(set,Dcdd,Dcdi)  - Conditional Density'
            write(*,*) '   CM                  - Cylinders Method'
            write(*,*) ; write(*,*) '   extend (y/n)?'

            read(*,*) text ; if (trim(text)=='y') then

               write(*,*) '   EV(A,dim(A))           - Expected value'
               write(*,*) '   maxeq(A,dim(A))        - maximal value of A'
               write(*,*) '   mineq(A,dim(A))        - minimal value of A'
               write(*,*) '   DescCoord(a,b,r,x,y,z) - conversion from spherical to Cartesian coordinates'
               write(*,*) '   SphCoord(x,y,z,a,b,r)  - conversion from Cartesian to spherical coordinates'
               write(*,*)
               write(*,*) '   fztd(z)              - conversion from redshift to Mpc'
               write(*,*) '   GalCoord(RA,Dec,l,b) - conversion from galactic to ecliptical coordinates'
               write(*,*) '   EclCoord(RA,Dec,l,b) - conversion from ecliptical to galactic coordinates'
               write(*,*)
               write(*,*) '   filevolume(in,N)   - strings amount of file'
               write(*,*) '   readfile(in,N,set) - reading from file'
               write(*,*) '   exist(in,i)        - testing of exist of file'
               write(*,*) '   strreal(a)         - Convert an integer to string'
               write(*,*) '   strint(i)          - Convert an integer to string'
               write(*,*) '   help               - this description'
               write(*,*) '   opening            - first sqript'
               write(*,*) '   ending             - end sqript'
               write(*,*) '   last update 16.03.2017'

               write(*,*) ; write(*,*) 'end of HELP' ; write(*,*) ; endif

         end subroutine help !==============================================================================2

      subroutine ending !================================================================================1
         write(*,*) ; write(*,*) 'All algorithms are implemented successfully' ; write(*,*) '---=1---'
         end subroutine ending !============================================================================1

end module
