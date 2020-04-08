!=- fortran-libraries
!=- © Stanislav Shirokov, 2014-2020

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
	module GNUplot
		use global
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
		integer, parameter	::	N_GNUfields=60, &
			ls1=1, ls2=2, ls3=3, ls4=4, ls5=5, title=7,xlabel=9,xrange=10,mxtics=11,format_x=12, &
			ylabel=14,yrange=15,mytics=16,format_y=17,grid=19,logscale=20,add1=23,add2=24,add3=25, add4=26, add5=27, add6=28, &
			extention_out_figure=29,legend=31,plot1=33,title1=34,plot2=35,title2=36,plot3=37,title3=38,plot4=39,title4=40, &
			plot5=41,title5=42,plot6=43,title6=44,plot7=45,title7=46,plot8=47,title8=48,plot9=49,title9=50,plot10=51,title10=52, &
			plot11=53,title11=54,plot12=55,title12=56,plot13=57,title13=58,plot14=59,title14=60,set_parametric=22

		character(len) ::	pngterm  = 'set term pngcairo enhanced font "Verdana,10" size 850, 630'				, &
								epsterm  = 'set term postscript enhanced color font "Verdana,14" size 8.5, 6.3'	, &
								scriptpath  = '', figurepath  = '',	figoutput  = '', GNUdatafile = '', graph_name='', &
								fig_dir = '', &

								GNUfields(N_GNUfields)=(/ &
		' set linestyle 1 lw 1 pt 7 ps 0.7 lt rgb "blue"  # 1  line style 1       ', &
		' set linestyle 2 lw 1 pt 7 ps 0.7 lt rgb "gray"  # 2  line style 2       ', &
		' set linestyle 3 lw 3 pt 2 ps 0.5 lt rgb "red"   # 3  line style 3       ', &
		' set linestyle 4 lw 3 pt 2 ps 0.5 lt rgb "green" # 4  line style 4       ', &
		' set linestyle 5 lw 3 pt 2 ps 0.5 lt rgb "gray"  # 5  line style 5       ', &
		'                                                                         ', &
		' set title "title"                               # 7  figure title       ', &
		'                                                                         ', &
		' set xlabel "x"                                  # 9  x-axis label       ', &
		' #set xrange [0:1]                                # 10                   ', &
		' set mxtics 5                                    # 11                    ', &
		' #set format x "10^{%L}"                          # 12                   ', &
		'                                                                         ', &
		' set ylabel "y"                                  # 14 y-axis label       ', &
		' #set yrange [0:1]                                # 15                   ', &
		' set mytics 5                                    # 16                    ', &
		' #set format y "10^{%L}"                          # 17                   ', &
		'                                                                         ', &
		' set grid xtics ytics mxtics mytics              # 19                    ', &
		' #set logscale                                    # 20                   ', &
		'                                                                         ', &
		' #set parametric                                  # 22 constants plots   ', &
		' #approx(x) = a*x + b                            # 23 LSM function       ', &
		' #fit approx(x) "file.dat" u 1:2 via a,b         # 24 LSM solution       ', &
		' #additional field 1                             # 25                    ', &
		' #additional field 2                             # 26                    ', &
		' #additional field 3                             # 27                    ', &
		' #additional field 4                             # 28                    ', &
		' #eps                                            # 29 output expansion   ', &
		'                                                                         ', &
		' set key off                                     # 31 legend option      ', &
		'                                                                         ', &
		', "" u 1:2 w l                                                           ', &		!-=- 33 plot options 1
		' notitle                                                                 ', &		!-=- 34 legend title 1

		'                                                                         ', &		!-=- 35 plot options 2
		'                                                                         ', &		!-=- 36 legend title 2

		'                                                                         ', &		!-=- 37 plot options 3
		'                                                                         ', &		!-=- 38 legend title 3
		'                                                                         ', &		!-=- 39 plot options 4
		'                                                                         ', &		!-=- 40 legend title 4
		'                                                                         ', &		!-=- 41 plot options 5
		'                                                                         ', &		!-=- 42 legend title 5
		'                                                                         ', &		!-=- 43 plot options 6
		'                                                                         ', &		!-=- 44 legend title 6
		'                                                                         ', &		!-=- 45 plot options 7
		'                                                                         ', &		!-=- 46 legend title 7
		'                                                                         ', &		!-=- 47 plot options 8
		'                                                                         ', &		!-=- 48 legend title 8
		'                                                                         ', &		!-=- 49 plot options 9
		'                                                                         ', &		!-=- 50 legend title 9
		'                                                                         ', &		!-=- 51 plot options 10
		'                                                                         ', &		!-=- 52 legend title 10
		'                                                                         ', &		!-=- 53 plot options 11
		'                                                                         ', &		!-=- 54 legend title 11
		'                                                                         ', &		!-=- 55 plot options 12
		'                                                                         ', &		!-=- 56 legend title 12
		'                                                                         ', &		!-=- 57 plot options 13
		'                                                                         ', &		!-=- 58 legend title 13
		'                                                                         ', &		!-=- 59 plot options 14
		'                                                                         '/)			!-=- 60 plot options 14

		real(8) :: GNUnorm=0d0

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
		contains

			character(len) function slashfix(path)
				character(len) path; slashfix = ''; slashfix = path;
					do i=len,1,-1 ; if (path(i:i)=='\') slashfix(i:i)='/' ; enddo
				end function

			character(len) function backslashfix(path)
				character(len) path; backslashfix = ''; backslashfix = path;
					do i=len,1,-1 ; if (path(i:i)=='/') backslashfix(i:i)='\' ; enddo
				end function

			integer function exptension(key) !=- 0 = default, 1 = eps, ...
				character(len) key; exptension=0; do i=1,len-2 ; if (key(i:i+2)=='eps') exptension=1 ; enddo
				end function

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
			subroutine plot(datafile) !-=- plot autoscripting
				character(len) datafile
					GNUdatafile=slashfix(datafile); call file_name(datafile,ii,jj)
					scriptpath = datafile(1:ii-1)//'plots\GNUscripts\'//datafile(ii:jj)//'_'//trim(graph_name)//'.plt'

					call system('MD ' // datafile(1:ii-1)//'plots\GNUscripts\' // ' >> log.log ' )
					call system('MD ' // datafile(1:ii-1)//'plots\'//datafile(ii:jj)//'\'//trim(fig_dir) // ' >> log.log ' )

					open(9,file=scriptpath,status='replace')
						if (exptension(GNUfields(29))==1) then
							write(9,*) epsterm
							figoutput = GNUdatafile(1:ii-1)//'plots/GNUscripts/' &
                        //datafile(ii:jj)//'_'//trim(graph_name)//'.eps'
							else
								write(9,*) pngterm
								figoutput = GNUdatafile(1:ii-1)//'plots/GNUscripts/' &
                           //datafile(ii:jj)//'_'//trim(graph_name)//'.png'
								endif;

                  write(9,*) 'set output "'//trim(figoutput)//'"'
                     do i=1,32
                        write(9,'(A)') trim(GNUfields(i))
                        enddo
						!write(9,*) 'plot "<echo -1 10**'//trim(realtostr(GNUnorm))//'" w p ls 10 notitle \'
						write(9,*) 'plot ',(trim(GNUfields(i)),i=33,N_GNUfields) ; close(9);

					figoutput=backslashfix(figoutput)
					call system( 'gnuplot ' // trim(scriptpath) // ' && move ' // trim(figoutput) // ' ' // &
						datafile(1:ii-1)//'plots\'//datafile(ii:jj)//'\'//trim(fig_dir) )

					open(8,file=datafile(1:ii-1)//'plots\GNUscripts\clear.bat',status='replace')
						write(8,*) 'taskkill /im wgnuplot*' ; close(8)



					graph_name='' ; GNUnorm=0d0

				end subroutine



      subroutine clear_plot_fields

         GNUfields( : ) = ''  !=- default
         fig_dir=''
         pngterm  = 'set term pngcairo enhanced font "Verdana,10" size 850, 630' !=- default
         epsterm  = 'set term postscript enhanced color font "Verdana,14" size 8.5, 6.3'

         end subroutine



		end module
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
