!=- Fractal Dimension Estimation
!=- © Stanislav Shirokov, 2014-2021

module FDE_citadel

      use global

      use FDE_config
      use FDE_paths
      use FDE_scripts
      use FDE_generators
      use FDE_overleaf

      contains

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine commands

               11 continue
            call preparation

				do
               call input_command(commandN)
					select case (commandN(1))

               !=- The FDE custom command list
               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

                  case('matrix')
                     call gen_FDE_matrix_script

               !=- The testing custom command list
               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

						case ('m')
							call FDE_testing_MultiCantor_script

               !=- removed from the release version
               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

							!call FDE_testing_script
                     !call FDE_analysis_external_catalog
                     !call add_visible_z_script
                     !call scaling_script

               !=- The FDE-system command list
               !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

						case ('help')
							call help

                  case ('reset')
                     goto 11

                  case('exit')
                     exit

						case ('delete')
                     call FDE_delete

                  case default
                     write(*,*) 'unknown command'

                  end select

               call set_default
               call ending

               enddo
            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine FDE_delete
            integer i

               select case (commandN(2))
                  case default   !=- deletes methods & samples (only files)
                     do i=1,5
                        call system('>> log.txt (DEL /S/Q '// trim(folders(i)) //')')
                        enddo
                  case('-m')     !=- deletes all methods (only files)
                     do i=2,5
                        call system('>> log.txt (DEL /S/Q '// trim(folders(i)) //')')
                        enddo
                  case('-s')     !=- deletes samples (files and folders)
                        call system('>> log.txt (RD /S/Q '// trim(folders(1)) //')')   !=- BUG: requires execution of reset command
                  case('-so')
                     do i=1,5    !=- deletes methods & samples (files and folders)
                        call system('>> log.txt (RD /S/Q '// trim(folders(i)) //')')   !=- BUG: requires execution of reset command
                        enddo
                  end select

            end

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine preparation
               call define_system
               call gen_paths
               call shell_MD_Tree
               call the_catalog_heading
               call geometry_texts
                  write(*,*) 'The FDE Citadel has been runned'
            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine input_command(commandN)
            character(len) command,commandN(:)
               33 continue ; write(*,*) '--' ; commandN(:)='NaN'

					read(*,'(A256)',err=33,end=44) command
                  read(command,*,iostat=iostat_value,err=33,end=44) commandN ; 44 continue

            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine set_default  !=- ?

            !call overleaf_default
            !call math_default

            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

		subroutine help
			implicit none ; character(10) text ; write(*,*)

            write(*,*) ; write(*,*) '   extend (y/n)?'

            read(*,*) text ; if (trim(text)=='y') then

               write(*,*) ; write(*,*) 'end of HELP' ; write(*,*) ; endif

         end

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

      subroutine ending
         write(*,*) ; write(*,*) 'All algorithms are implemented successfully' ; write(*,*) '---=1---'
         end

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

end module
