!=- fortran-libraries
!=- © Stanislav Shirokov, 2014-2020

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
	module GNUplot
		use global
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
		integer, parameter	::	N_GNUfields = 200

      integer :: ls1=1, ls2=2, ls3=3, ls4=4, ls5=5, title=6,xlabel=9,xrange=10,mxtics=11,format_x=12, &
			ylabel=14,yrange=15,mytics=16,format_y=17,grid=19,logscale=20,add1=23,add2=24,add3=25, add4=26, add5=27, add6=28, &
			extention_out_figure=29,legend=31,plot1=33,title1=34,plot2=35,title2=36,plot3=37,title3=38,plot4=39,title4=40, &
			plot5=41,title5=42,plot6=43,title6=44,plot7=45,title7=46,plot8=47,title8=48,plot9=49,title9=50,plot10=51,title10=52, &
			plot11=53,title11=54,plot12=55,title12=56,plot13=57,title13=58,plot14=59,title14=60,set_parametric=22, &
			thickness , comma_fix , dash_type

      real(8) :: external_parameter = 0

		character(length) ::	pngterm  = 'set term pngcairo enhanced font "Verdana,12" size 1400, 1050'				, &
								epsterm  = 'set term postscript enhanced color font "Verdana,14" size 8.5, 6.3'	, &
								scriptpath  = '', figurepath  = '',	figoutput  = '', GNUdatafile = '', graph_name='', &
								fig_dir = '', GNUplots = '', GNUscripts = '', script_file_path='', plot_dir='', &
                        plot_extension='', GNUdatafile_win='', &

                        GNUfields(N_GNUfields), GNUplot_colors( N_GNUfields/2 )

!		' set linestyle 1 lw 1 pt 7 ps 0.7 lt rgb "blue"  # 1  line style 1       ', &
!		'                                                                         ', &
!		' set title "title"                               # 7  figure title       ', &
!		'                                                                         ', &
!		' set xlabel "x"                                  # 9  x-axis label       ', &
!		' set xrange [0:1]                                # 10                    ', &
!		' set mxtics 5                                    # 11                    ', &
!		' set format x "10^{%L}"                          # 12                    ', &
!		'                                                                         ', &
!		' set ylabel "y"                                  # 14 y-axis label       ', &
!		' set yrange [0:1]                                # 15                    ', &
!		' set mytics 5                                    # 16                    ', &
!		' set format y "10^{%L}"                          # 17                    ', &
!		'                                                                         ', &
!		' set grid xtics ytics mxtics mytics              # 19                    ', &
!		' set logscale                                    # 20                    ', &
!		'                                                                         ', &
!		' set parametric                                  # 22 constants plots    ', &
!		' #approx(x) = a*x + b                            # 23 LSM function       ', &
!		' #fit approx(x) "file.dat" u 1:2 via a,b         # 24 LSM solution       ', &
!		' #additional field 1                             # 25                    ', &
!		' #eps                                            # 29 output expansion   ', &
!		'                                                                         ', &
!		' set key off                                     # 31 legend option      ', &
!		'                                                                         ', &
!		', "" u 1:2 w l                                                           ', &		!-=- 33 plot options 1
!		' notitle                                                                 ', &		!-=- 34 legend title 1

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
		contains

			integer function extension(key) !=- 0 = default, 1 = eps, ...
				character(length) key; extension=0; do i=1,length-2 ; if (key(i:i+2)=='eps') extension=1 ; enddo
				end function

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1

			subroutine plot(datafile) !-=- plot autoscripting
				character(length) datafile

                  unit_2 = random_unit()

                  select case (operating_system)
                     case(1)  !=- Windows
                        datafile    = backslashfix(datafile)
                        GNUplots    = 'plots\'
                        GNUscripts  = 'plots\GNUscripts\'
                        plot_dir    = trim(file_path(datafile)) // trim(GNUplots) // trim(name(datafile)) &
                           // '\' // trim(fig_dir)
                     case(2)  !=- Linux
                        datafile    = slashfix(datafile)
                        GNUplots    = 'plots/'
                        GNUscripts  = 'plots/GNUscripts/'
                        plot_dir    = trim(file_path(datafile)) // trim(GNUplots) // trim(name(datafile)) &
                           // '/' // trim(fig_dir)
                     end select

                  output = trim(file_path(datafile)) // trim(GNUscripts)
                  call create_folder( output )
                  call create_folder( plot_dir )

                  figoutput = trim(file_path(datafile)) // trim(GNUscripts) // trim(name(datafile)) // &
                     '_' // trim(graph_name)
                  script_file_path  = trim(figoutput) // '.plt'
					open(unit_2,file = script_file_path ,status='replace')

						select case ( extension( GNUfields(extention_out_figure) ) )
                     case (1)
                        write(unit_2,*) epsterm
                        plot_extension = '.eps'
                     case default
								write(unit_2,*) pngterm
                        plot_extension = '.png'
                     end select

                  figoutput = trim(figoutput) // trim(plot_extension)

                     call GNUfields_fix
                  write(unit_2,*) 'set output "'//trim(slashfix(figoutput))//'"'
                     do i=1,plot1-1
                        write(unit_2,'(A)') trim(GNUfields(i))
                        enddo
						!write(9,*) 'plot "<echo -1 10**'//trim(realtostr(GNUnorm))//'" w p ls 10 notitle \'
						write(unit_2,*) 'plot ',( trim(GNUfields(i)) , i = plot1 , N_GNUfields )

						close(unit_2)

					call system( 'gnuplot ' // trim(script_file_path) // ' >> log.log ' )
               call move_file( figoutput , plot_dir )

					graph_name=''

				end subroutine



      logical function GNUcorrect(field)
         character(length) field
            GNUcorrect = .false.
            if ( field /= '' ) then
               read(field,*) line
               if ( line == 'set' .or. line == '#set' .or. line == '#' ) GNUcorrect = .true.
               else
                  GNUcorrect = .true.
               endif
         end function GNUcorrect



      subroutine GNUfields_fix

            if ( .not. GNUcorrect( GNUfields( title ) ) ) &
               GNUfields( title ) = 'set title "' // trim(GNUfields( title )) // '"'

            if ( .not. GNUcorrect( GNUfields( xlabel ) ) ) &
               GNUfields( xlabel ) = 'set xlabel "' // trim(GNUfields( xlabel )) // '"'
            if ( .not. GNUcorrect( GNUfields( xrange ) ) ) &
               GNUfields( xrange ) = 'set xrange '  // trim(GNUfields( xrange )) // ''
            if ( .not. GNUcorrect( GNUfields( mxtics ) ) ) &
               GNUfields( mxtics ) = 'set mxtics '  // trim(GNUfields( mxtics )) // ''

            if ( .not. GNUcorrect( GNUfields( ylabel ) ) ) &
               GNUfields( ylabel ) = 'set ylabel "' // trim(GNUfields( ylabel )) // '"'
            if ( .not. GNUcorrect( GNUfields( yrange ) ) ) &
               GNUfields( yrange ) = 'set yrange '  // trim(GNUfields( yrange )) // ''
            if ( .not. GNUcorrect( GNUfields( mytics ) ) ) &
               GNUfields( mytics ) = 'set mytics '  // trim(GNUfields( mytics )) // ''

            if ( .not. GNUcorrect( GNUfields( legend ) ) ) &
               GNUfields( legend ) = 'set key '  // trim(GNUfields( legend )) // ''

            if ( .not. GNUcorrect( GNUfields( grid ) ) ) &
               GNUfields( grid ) = 'set grid '  // trim(GNUfields( grid )) // ''

            do i=title1,N_GNUfields,2

               if ( i > title1 .and. GNUfields(i-1)(1:1) /= ',' .and. GNUfields(i-1) /= '' ) &
                  GNUfields(i-1) = ', ' // trim(GNUfields(i-1))

               if ( GNUfields(i) /= '' .and. .not. word_search(GNUfields(i),'title') &
                  .and. .not. word_search(GNUfields(i),'notitle') ) &
                  GNUfields(i) = ' title "' // trim(GNUfields(i)) // '"'
               end do

            if ( .not. GNUcorrect( GNUfields( format_y ) ) ) &
               GNUfields( format_y ) = 'set format y "' // trim(GNUfields( format_y )) // '"'
            if ( .not. GNUcorrect( GNUfields( format_x ) ) ) &
               GNUfields( format_x ) = 'set format x "' // trim(GNUfields( format_x )) // '"'

            if ( GNUfields(plot1)(1:1) == ',' ) GNUfields(plot1) = GNUfields(plot1)(2:length)

         end subroutine GNUfields_fix


      subroutine clear_plot_fields

         GNUfields( : ) = ''  !=- default
         fig_dir  = ''
         pngterm  = 'set term pngcairo enhanced font "Verdana,16" size 1400, 1050'
         epsterm  = 'set term postscript enhanced color font "Verdana,14" size 8.5, 6.3'

         !=- Verdana , Palatino-Roman

         end subroutine



		end module
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
