module tests
	use global
	use GRBurst_paths
	use math
	use GNUplot

		integer,parameter :: N_col_test = 1d1, N_raw_test = 2d3

		integer 				:: j_test = 0, N_test = 2d2, N_sample = 2d2, N_grid_test = 1d2

		real(8) 				:: sample(N_col_test,N_raw_test) = 0d0, a_test = 1d0

	contains

		subroutine regression_test
call GNUplot_verRegression; pause
				open(21,file=VerRegress,status='replace')
				sa=0d0; sb=0d0 ; ssa=0d0 ; ssb=0d0
					do j_test=1,N_grid_test

						N_grid_LT_a	=	4d2
						N_sample		=	6d1
						a_test		=	2.5d0

						a_max			=	3d0

						aa = j_test*(a_test*2/N_grid_test) - a_test	;	bb = 52d0

						call create_sample(sample,aa,bb,N_sample)
						call logtrendX2(sample(1:4,:))	;	write(*,*) j_test

						write(21,'(10(1x,F6.2)2(1x,i5))') aa-a, bb-b, &
							aa, a, a_error_lower, a_error_upper, bb, b, b_error_lower, b_error_upper, &
							N_sample,N_grid_LT
sa=sa+aa-a
sb=sb+bb-b
ssa=ssa+abs(aa-a)
ssb=ssb+abs(bb-b)
						write(*,'(10(1x,F6.2)2(1x,i5))') aa-a, bb-b, &
							aa, a, a_error_lower, a_error_upper, bb, b, b_error_lower, b_error_upper, &
							N_sample,N_grid_LT
						 !=-call GNUplot_sample(Sample,aa,bb,N_sample) ; pause
						enddo
					close(21)
					sa =sa/N_grid_test
					sb =sb/N_grid_test
					ssa=ssa/N_grid_test
					ssb=ssb/N_grid_test
write(*,*) 'sa=,sb=',sa,sb,ssa,ssb
call GNUplot_verRegression

			end subroutine

      subroutine GNUplot_verRegression
				GNUfields(1)  = 'set linestyle 1 lw 2 pt 7 ps 0.5 lt rgb "gray"'
				GNUfields(2)  = 'set linestyle 2 lw 2 pt 7 ps 0.9 lt rgb "black"'
				GNUfields(3)  = 'set linestyle 3 lw 2 pt 7 ps 0.9 lt rgb "blue"'
				GNUfields(4)  = 'set linestyle 4 lw 1 pt 7 ps 0.9 lt rgb "red"'

		GNUfields(9)  = 'set xlabel "a"'
		GNUfields(14) = 'set ylabel "a_e"'
		GNUfields(7)  = 'set title "testing ( y + y_n ) = a ( x + x_n ) + b, N = '//trim(inttostr(N_sample))//'"'

				GNUfields(20) = '#set logscale xy'
				GNUfields(17) = '#set format y "10^{%L}"'
				GNUfields(12) = '#set format x "10^{%L}"'
				GNUfields(31) = 'set key bottom right'
				GNUfields(29) = '#set key off'
				GNUfields(10) = '#set xrange[]'
				!GNUfields(15) = 'set yrange [0.00001:1]'

				GNUnorm = log10( 1d-100 )

				GNUfields(33) = ', "" u 3:($1*0) w l ls 1'
				GNUfields(34) = ' notitle "a - a_e"'
				GNUfields(35) = ', "" u 3:3 w l ls 1'
				GNUfields(36) = ' title "x*x"'

				GNUfields(37) = ', "" u 3:($1*0+'//trim(realtostr(ssa))//') w l ls 4'
				GNUfields(38) = ' notitle "1{/Symbol s}(a_e)"'
				GNUfields(39) = ', "" u 3:($1*0-'//trim(realtostr(ssa))//') w l ls 4'
				GNUfields(40) = ' notitle "a - a_e"'
				GNUfields(41) = ', "" u 3:($3+'//trim(realtostr(ssa))//') w l ls 4'
				GNUfields(42) = ' title "1{/Symbol s}(a_e)='//trim(realtostr(ssa))//'"'
				GNUfields(43) = ', "" u 3:($3-'//trim(realtostr(ssa))//') w l ls 4'
				GNUfields(44) = ' notitle "x*x"'

				GNUfields(45) = ', "" u 3:(-$1) w l ls 2'
				GNUfields(46) = ' title "a - a_e"'
				GNUfields(47) = ', "" u 3:4 w l ls 3'
				GNUfields(48) = ' title "a_e = a_e(a)"'

				call plot(VerRegress);pause 2
GNUfields(33:48)=''
				GNUfields(15) = 'set yrange [47:57]'
            GNUfields(9)  = 'set xlabel "b"'
            GNUfields(14) = 'set ylabel "b_e"'

				GNUfields(33) = ', "" u 3:7 w l ls 1'
				GNUfields(34) = ' title "b"'

GNUnorm = log10( 52d0 )
				GNUfields(37) = ', "" u 3:($7+'//trim(realtostr(ssb))//') w l ls 4'
				GNUfields(38) = ' title "1{/Symbol s}(b_e)='//trim(realtostr(ssb))//'"'
				GNUfields(39) = ', "" u 3:($7-'//trim(realtostr(ssb))//') w l ls 4'
				GNUfields(40) = ' notitle "b - b_e"'

				GNUfields(47) = ', "" u 3:8 w l ls 3'
				GNUfields(48) = ' title "b_e"'

				call plot(VerRegress);pause 2

      end subroutine


		subroutine create_sample(Sample,aa,bb,N_sample)
			real(8) aa,bb,Sample(:,:) ; integer N_sample
				Sample(:,:)=0d0	;	tx =	3d0

				do i=1,N_sample
				call random_number(Sample(1,i))
				call random_number(Sample(2,i))
				call random_number(Sample(3,i))
				call random_number(Sample(4,i))

				nx = 2*(Sample(1,i)-0.5) !+ 1.5d0
				ny = 2*(Sample(2,i)-0.5)

				!Sample(3:4,i)=0d0

				Sample(1,i) = 1d1**( nx + tx )
				Sample(2,i) = 1d1**( ny + log10(Sample(1,i))*aa + bb )
				Sample(3,i) = 0.2d0+Sample(3,i)*2d-1
				Sample(4,i) = 0.2d0+Sample(4,i)*2d-1
				!write(*,*) ny + log10(Sample(1,i))*aa + bb,ny,log10(Sample(1,i)),aa ,bb
				!write(*,*) Sample(1,i),Sample(2,i),Sample(3,i),Sample(4,i); pause
				!write(*,*) ( Sample(j,i) ,j=1,4) ; pause
				!write(*,*) ( log10(Sample(j,i)) ,j=1,4) ; pause
				enddo
			end subroutine

		subroutine write_sample(Sample,aa,bb,N_sample)
			real(8) aa,bb,Sample(:,:) ; integer N_sample
				open(12,file=Regression,status='replace')
					do i=1,N_sample
						Sample(:,i)=log10(Sample(:,i))
						write(12,*) Sample(1,i),Sample(2,i),Sample(3,i),Sample(4,i), &
							1d1**(aa*log10(Sample(1,i))+bb),1d1**(a*log10(Sample(1,i))+b)
						write(*,'(4(F6.2),6(E12.4))') a, aa, b, bb, &
							Sample(1,i),Sample(2,i),Sample(3,i),Sample(4,i), aa*Sample(1,i)+bb, a*Sample(1,i)+b
						end do
					close(12)
			end subroutine

		subroutine GNUplot_sample(Sample,aa,bb,N_sample)
			real(8) aa,bb,Sample(:,:) ; integer N_sample

				call write_sample(Sample,aa,bb,N_sample)

				GNUfields(1)  = 'set linestyle 1 lw 1 pt 7 ps 0.5 lt rgb "gray"'
				GNUfields(2)  = 'set linestyle 2 lw 2 pt 7 ps 0.9 lt rgb "black"'
				GNUfields(3)  = 'set linestyle 3 lw 2 pt 7 ps 0.9 lt rgb "red"'

				GNUfields(20) = 'set logscale xy'
				GNUfields(17) = '#set format y "10^{%L}"'
				GNUfields(12) = '#set format x "10^{%L}"'
				GNUfields(31) = '#set key off'
				GNUfields(29) = '#set key off'
				GNUfields(10) = '#set xrange[]'
				!GNUfields(15) = 'set yrange [0.00001:1]'

				GNUnorm = log10( aa*log10(Sample(1,1))+bb )
				GNUfields(33) = ', "" u 1:2:3:4 w xyerrorb ls 1'
				GNUfields(34) = ' notitle "{/Symbol m}_{SN}(Pantheon)"'
				GNUfields(35) = ', "" u 1:5 w l ls 2'
				GNUfields(36) = ' title "original"'
				GNUfields(37) = ', "" u 1:6 w l ls 3'
				GNUfields(38) = ' title "approx"'
				!GNUfields(39) = ', "" u 1:7 w l ls 2'
				!GNUfields(40) = ' notitle "approx"'
				!GNUfields(41) = ', "" u 1:8 w l ls 2'
				!GNUfields(42) = ' notitle "approx"'

				call plot(Regression); pause 2

				GNUfields(33) = ', "" u 6:7:8:9 w xyerrorb ls 1'
				GNUfields(34) = ' notitle "{/Symbol m}_{SN}(Pantheon)"'
				GNUfields(35) = ', "" u 6:10 w p ls 2'
				GNUfields(36) = ' notitle "{/Symbol m}_{SN}(Pantheon)"'

				call plot(Regression); pause 3

			end subroutine

end module
