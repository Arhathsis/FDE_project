!=- fortran-libraries
!=- � Stanislav Shirokov, 2014-2020

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
	module global
		implicit none

      integer ::  operating_system  =  0   !=- 1 Windows, 2 Linux

      logical :: file_exists , index_info , write_percent_fix

		integer,parameter :: length=512,  len=512, N_col_GRB = 90, N_GRB=7d3, & !=- constants
			N_titles = 57, N_comparisons = 4, N_files = 20, N_folders = 10 , N_commands = 10 !=- 3+6*N_models

		integer ::  i,ii,iii,j,jj,jjj,k,kk,kkk,n,nn,nnn,m,mm,mmm, &
                  unit_1 = 1, unit_2 = 2, unit_3 = 3 , unit_4 = 4 , unit_5 = 5 , unit_6 = 6  , unit_7 = 7 , &
                  counter , final_counter , iostat_value , i_percent , N_percent , mod_percent = 10 , m_percent

		character(length) :: theformat='', theformat2='', datafile='', titles(N_titles)='', &
                        input='', output='', line='',preformat='',path='',head_format='', &
                        text1='',text2='',text3='',text4='',text5='', command='',figure_number='', &
                        files(N_files)='',new_files(N_files)='', folders(N_folders), &
                        system_commands(10), command1, command2, filepath='', columns, filename='', commandN(N_commands)

      real(8),parameter :: NaN = 1d150

		real(8) a,b,bbb,d,e,g,l,t,tx,ty,nx,ny,aa,bb,vva,vvb,q,p,sa,sb,ssa,ssb,sm, lb,rb, &
			x,xx,xxx, y,yy,yyy, z,zz,zzz, w,ww,www, v,vv,vvv,r, maximum, &
			alpha,beta,p1,p2,p3,p4,p5,sumx,sumx2,sumy,sumxy , &
			z_max_var, log_z, GRB_shift , &

			final_time , start_time

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
		contains



         subroutine define_system
            character(length) :: line = '' , system_log_file = 'system.log'

            call system(' echo %OS% > ' // trim(system_log_file))
               line = read_last_string(system_log_file)
               if ( line == 'Windows_NT' ) operating_system = 1

            call system(' echo $SHELL  >> ' // trim(system_log_file))
               line = trim( read_last_string(system_log_file) )
               if ( line == '/bin/bash'  ) operating_system = 2

            select case (operating_system)
               case(1)  !=-   Windows
                  write(*,*) 'OS: Windows_NT'
                  call sleep (1)
                  call system('cls')
               case(2)  !=-   Linux
                  write(*,*) 'OS: Linux'
                  call sleep (1)
                  !call system('clear')
               case default   !=-   Other
                  write(*,*) 'define_system: wrong mask'
                  call sleep (10) ; operating_system = 1
               end select
            end subroutine define_system



         integer function random_unit()
            real(8) x
               do
                  call random_number(x)
                  random_unit = int(1d3*x)
                  if (random_unit>1d2) exit
                  enddo
            end function



         character(length) function read_last_string(log_pathfile)
            character(length) log_pathfile
               read_last_string = ''
                        unit_1 = random_unit()
            open(       unit_1,file=log_pathfile,status='old',err=22)
               do ; read(    unit_1,*,end=11) read_last_string ; enddo
               22 continue ; write(*,*) 'read_last_string critical error: read_last_string: file ', &
                  trim(log_pathfile),' does not exist'
               11 continue ; close(   unit_1)
            end function



         subroutine prepare_percent( N_max )
            integer N_max
            i_percent=0 ; N_percent = N_max ; call CPU_TIME(start_time) ; m_percent=0
            end subroutine prepare_percent



         subroutine write_percent   !=- before a cycle set `call prepare_percent( N_max )`
               i_percent = 1+i_percent
               a = 1d2*i_percent/N_percent
               if ( int(a)>m_percent .and. int(a)>0 .and. (mod(int(a),mod_percent) == 0 .or. int(a) == 1 .or. int(a) == 99) ) then
                  m_percent=int(a)
                     call CPU_TIME(final_time)
                  b = final_time - start_time
         write(*,'(6x,"complited",i4,"%, work time is",F10.2,"s, time left is",F10.2,"s, processed points are",i10," from",i10)') &
                     int(a),b,b/int(a)*1d2-b, i_percent,N_percent
                  endif
            end subroutine write_percent



         character(length) function make_workdir()
            character(length) :: line = '' , log_file = 'log.log'

            select case (operating_system)
               case(1)  !=-   Windows

                  call system( 'echo %cd% > ' // trim(log_file) )
                  line = read_last_string(log_file)
                  make_workdir = file_path(line)

               case(2)  !=-   Linux

                  call system( 'echo $PWD > ' // trim(log_file) )
                  line = read_last_string(log_file)
                  make_workdir = file_path(line)

               case default   !=-   Other
                  write(*,*) 'fatal error: make_workdir: unknown mask'
               end select

            end function make_workdir



         character(length) function name(path)
            character(length) path ; integer ii,jj
               call file_name(path,ii,jj); name=''; name=path(ii:jj)
            end function

         character(len) function name_ext(path)
            character(len) path; integer ii,jj
               call file_name(path,ii,jj); name_ext=''; name_ext=path(ii:len)
            end function

			subroutine file_name(path,l1,l2) !-=- file path excision
				character(length) path ; integer i,l1,l2
					do i=length,1,-1 ; if (path(i:i)=='.') then ; l2=i-1 ; exit ; endif ; enddo
					do i=length,1,-1 ; if (path(i:i)=='\'.or.path(i:i)=='/') then ; l1=i+1 ; exit ; endif ; enddo
				end subroutine file_name

			character(length) function file_path(path) !-=- file name excision
				character(length) path ; integer l1 , l2 , i
					do i=length,1,-1 ; if (path(i:i)=='\'.or.path(i:i)=='/') then ; l1=i ; exit ; endif ; enddo
               file_path='' ; file_path = path(1:l1)
				end function

			character(length) function inttostr(integer_number) !=- Convert an integer to string
				integer integer_number ; write(inttostr,*) integer_number
				end function inttostr

			character(length) function inttostrf(integer_number,order) !=- Convert an integer to string
				integer integer_number,order ; character(length) theformat
               theformat='(i'// trim(inttostr(order)) //')'
					write(inttostrf,theformat) integer_number
				end function inttostrf

			character(length) function realtostr(real_number) !=- Convert an real to string
				real(8) real_number ; write(realtostr,'(F10.4)') real_number
				end function realtostr

			character(length) function realtostrE(real_number) !=- Convert an real to string
				real(8) real_number ; write(realtostrE,'(E20.12)') real_number
				end function realtostrE

			character(length) function realtostrf(real_number,order,mantissa) !=- Convert an format real to string
				real(8) real_number ; integer order,mantissa ; character(length) theformat
					theformat='(F'// trim(inttostr(order)) //'.'// trim(inttostr(mantissa)) //')'
					write(realtostrf,theformat) real_number
				end function realtostrf

			character(length) function realtostrff(real_number,theformat) !=- INVALID
				real(8) real_number ; character(length) theformat
					write(realtostrff,'('//trim(theformat)//')') real_number
				end function realtostrff

			integer function file_volume(file_path)
				character(length) file_path
				integer unit_4

               unit_4 = random_unit() ; file_volume=0

               inquire( file = file_path , exist = file_exists )
					if ( file_exists ) then
						open(unit_4,file=file_path,status='old')
							do
								read(unit_4,*,end=11) line
									if ( .not. symbol_search(line,'#') ) file_volume = 1 + file_volume
								enddo
							11 continue ; close(unit_4)
						else
						write(*,*) '   file ',trim(file_path),' does not exist'
						endif
				end function



			logical function symbol_search(line,symbol)
            integer i
				character(length) line ; character(1) symbol ; symbol_search = .false.
					do i=1,length
						if (line(i:i)==symbol) then
                     symbol_search = .true.
                     exit
                     endif
						end do
				end function

			logical function word_search(line,goal)   !=- of only the first word - ?
				character(length) word,line ; character(*) goal ; word_search = .false.
               read(line,*) word
               if (word==goal) word_search = .true.
				end function



         subroutine shell_MD_Tree
            integer i
               do i=1,size(folders(:))
                  if (folders(i).ne.'NaN') call create_folder(folders(i))
                  enddo
            end subroutine shell_MD_Tree



         subroutine move_file(whence,here)
            character(len) whence , here
               select case (operating_system)
                  case(1)  !=-   Windows
                     whence = backslashfix(whence)
                       here = backslashfix(here)
                     if ( whence.ne.'' .and. here.ne.'' .and. whence.ne.'NaN' .and. here.ne.'NaN' ) &
                        call system( 'move ' // trim(whence) // ' ' // trim(here)  // ' >> log.log ' )
                  case(2)  !=-   Linux
                     whence = slashfix(whence)
                       here = slashfix(here)
                     if ( whence.ne.'' .and. here.ne.'' .and. whence.ne.'NaN' .and. here.ne.'NaN' ) &
                        call system( 'mv ' // trim(whence) // ' ' // trim(here)  // ' >> log.log ' )
                  end select
            end subroutine


         subroutine copy_file(whence,here)
            character(len) whence , here
               select case (operating_system)
                  case(1)  !=-   Windows
                     whence = backslashfix(whence)
                       here = backslashfix(here)
                     if ( whence.ne.'' .and. here.ne.'' .and. whence.ne.'NaN' .and. here.ne.'NaN' ) &
                        call system( 'copy ' // trim(whence) // ' ' // trim(here)  // ' >> log.log ' )
                  case(2)  !=-   Linux
                     whence = slashfix(whence)
                       here = slashfix(here)
                     if ( whence.ne.'' .and. here.ne.'' .and. whence.ne.'NaN' .and. here.ne.'NaN' ) &
                     call system( 'cp ' // trim(whence) // ' ' // trim(here)  // ' >> log.log ' )
                  end select
            end subroutine


         subroutine create_folder(path_folder)
            character(length) path_folder
               select case (operating_system)
                  case(1)  !=-   Windows
                     path_folder = backslashfix(path_folder)
                     if (path_folder.ne.'') call system( 'If Not Exist ' // trim(path_folder) // &
                        ' MD ' // trim(path_folder)  // ' >> log.log ' )
                  case(2)  !=-   Linux
                     path_folder = slashfix(path_folder)
                     if (path_folder.ne.'') call system( 'mkdir -p ' // trim(path_folder)  // ' >> log.log ' )
                  end select
            end subroutine



			character(length) function slashfix(path)
            integer i
				character(length) path; slashfix = ''; slashfix = path
					do i=length,1,-1 ; if (path(i:i)=='\') slashfix(i:i)='/' ; enddo
				end function

			character(length) function backslashfix(path)
            integer i
				character(length) path; backslashfix = ''; backslashfix = path
					do i=length,1,-1 ; if (path(i:i)=='/') backslashfix(i:i)='\' ; enddo
				end function



		end module
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1

