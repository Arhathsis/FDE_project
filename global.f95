!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
	module global
		implicit none

      logical :: file_exists

		integer,parameter :: len=256,  N_col_GRB = 90, N_GRB=7d3, & !=- constants
			N_titles = 57, N_comparisons = 4, N_files = 20, N_folders = 10 !=- 3+6*N_models

		integer ::  i,ii,iii,j,jj,jjj,k,kk,kkk,n,nn,nnn,m,mm,mmm, &
                  unit_1 = 1, unit_2 = 2, unit_3 = 3 , unit_4 = 4 , unit_5 = 5 , unit_6 = 6  , unit_7 = 7 , &
                  counter , final_counter

		character(len) :: theformat='', theformat2='', datafile='', titles(N_titles)='', &
                        input='', output='', line='',preformat='',path='',head_format='', &
                        text1='',text2='',text3='',text4='',text5='', command='',figure_number='', &
                        files(N_files)='',new_files(N_files)='', folders(N_folders) , &

                     !dir = 'D:\Arhath\YandexDisk\'          !=- D:\Arhath or C:\Users\
                     dir = 'C:\Users\Arhath\YandexDisk\'   !=- Dell notebook

		real(8) a,b,bbb,d,e,g,l,t,tx,ty,nx,ny,aa,bb,vva,vvb,q,p,sa,sb,ssa,ssb,sm, &
			x,xx,xxx, y,yy,yyy, z,zz,zzz, w,ww,www, v,vv,vvv,r, maximum, &
			alpha,beta,p1,p2,p3,p4,p5,sumx,sumx2,sumy,sumxy , &
			z_max_var, log_z, GRB_shift , &

			final_time , start_time

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
		contains



         character(len) function name(path)
            character(len) path; call file_name(path,ii,jj); name=''; name=path(ii:jj)
            end function

			subroutine file_name(path,l1,l2) !-=- file path excision
				character(len) path ; integer i,l1,l2
					do i=len,1,-1 ; if (path(i:i)=='.') then ; l2=i-1 ; exit ; endif ; enddo
					do i=len,1,-1 ; if (path(i:i)=='\'.or.path(i:i)=='/') then ; l1=i+1 ; exit ; endif ; enddo
				end subroutine file_name

			character(len) function file_path(path) !-=- file name excision
				character(len) path ; integer l1,l2
					do i=len,1,-1 ; if (path(i:i)=='\'.or.path(i:i)=='/') then ; l1=i ; exit ; endif ; enddo
               file_path='' ; file_path = path(1:l1)
				end function

			character(len) function inttostr(integer_number) !=- Convert an integer to string
				integer integer_number ; write(inttostr,*) integer_number
				end function inttostr

			character(len) function inttostrf(integer_number,order) !=- Convert an integer to string
				integer integer_number,order ; theformat='(i'// trim(inttostr(order)) //')'
					write(inttostrf,theformat) integer_number
				end function inttostrf

			character(len) function realtostr(real_number) !=- Convert an real to string
				real(8) real_number ; write(realtostr,'(F10.4)') real_number
				end function realtostr

			character(len) function realtostrE(real_number) !=- Convert an real to string
				real(8) real_number ; write(realtostrE,'(E20.12)') real_number
				end function realtostrE

			character(len) function realtostrf(real_number,order,mantissa) !=- Convert an format real to string
				real(8) real_number ; integer order,mantissa
					theformat='(F'// trim(inttostr(order)) //'.'// trim(inttostr(mantissa)) //')'
					write(realtostrf,theformat) real_number
				end function realtostrf

			character(len) function realtostrff(real_number,theformat) !=- INVALID
				real(8) real_number ; character(len) theformat
					write(realtostrff,trim(theformat)) real_number
				end function realtostrff

			integer function file_volume(file_path)
				character(len) file_path ; file_volume=0

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
				character(len) line ; character(1) symbol ; symbol_search = .false.
					do i=1,len
						if (line(i:i)==symbol) symbol_search = .true.
						end do
				end function


         subroutine shell_MD_Tree
               do i=1,size(folders(:))
						if (folders(i).ne.'') call system( 'MD '// trim(folders(i)) )
                  enddo
            end subroutine shell_MD_Tree


		end module
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1

