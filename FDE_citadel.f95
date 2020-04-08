module FDE_citadel
   use FDE_config
   use FDE_paths
   use FDE_sqripts
   use FDE_debug
   use FDE_generators

   contains

		subroutine preparation

				write(*,*) 'Fractalg Dimension Estimation is runned' ; write(*,*)

            call gen_paths
            call the_catalog_heading

			end subroutine

      subroutine opening !===============================================================================1

				call preparation

call debug

				do ; write(*,*) '--'

					read(*,*) command

					select case (command)
						case ('cg')
							!call compaire_graphics
						case ('exit')
							goto 11
						case ('help')
							call help
						case ('uca')
							call uniform_catalog_analysis
						case ('cca')
							call cantor_catalog_analysis
						case ('c')
							call input_catalogs_analysis

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
