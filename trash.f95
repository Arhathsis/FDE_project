
			GNUfields_ls1			=' set linestyle 1 lw 1 pt 7 ps 0.7 lt rgb "blue"  # 1  line style 1       ',	&

			GNUfields_ls2			='',	GNUfields_ls3			='',	GNUfields_ls4		='',	GNUfields_ls5	='', &

			GNUfields_title		=' set title "title"                               # 7  figure title       ',	&
			GNUfields_xlabel		=' set xlabel "x"                                  # 9  x-axis label       ',	&
			GNUfields_xrange		=' #set xrange [0:1]                                # 10                   ',	&
			GNUfields_mxtics		=' set mxtics 5                                    # 11                    ',	&
			GNUfields_format_x	=' #set format x "10^{%L}"                          # 12                   ',	&
			GNUfields_ylabel		=' set ylabel "y"                                  # 14 y-axis label       ',	&
			GNUfields_yrange		=' #set yrange [0:1]                                # 15                   ',	&
			GNUfields_mytics		=' set mytics 5                                    # 16                    ',	&
			GNUfields_format_y	=' #set format y "10^{%L}"                          # 17                   ',	&
			GNUfields_grid			=' set grid xtics ytics mxtics mytics              # 19                    ',	&
			GNUfields_logscale	=' #set logscale                                    # 20                   ',	&

			GNUfields_add1			='',	GNUfields_add2		='', &
			GNUfields_add3			='',	GNUfields_add4		='',	GNUfields_add5		='',	GNUfields_add6		='', &

			GNUfields_extention	=' #eps                                            # 29 output expansion   ',	&
			GNUfields_legend		=' set key off                                     # 31 legend option      ',	&
			GNUfields_plot1		=', "" u 1:2 w l                                                           ',	&
			GNUfields_title1		=' notitle                                                                 ',	&

			GNUfields_plot2		='',	GNUfields_title2	='',	GNUfields_plot3	='',	GNUfields_title3	='', &
			GNUfields_plot4		='',	GNUfields_title4	='',	GNUfields_plot5	='',	GNUfields_title5	='', &
			GNUfields_plot6		='',	GNUfields_title6	='',	GNUfields_plot7	='',	GNUfields_title7	='', &
			GNUfields_plot8		='',	GNUfields_title8	='',	GNUfields_plot9	='',	GNUfields_title9	='', &
			GNUfields_plot10		='',	GNUfields_title10	='',	GNUfields_plot11	='',	GNUfields_title11	='', &
			GNUfields_plot12		='',	GNUfields_title12	='',	GNUfields_plot13	='',	GNUfields_title13	='', &
			GNUfields_plot14		='',	GNUfields_title14	=''

			subroutine new_interface

					GNUfields(1)  = GNUfields_ls1			;	GNUfields(16) = GNUfields_mytics
					GNUfields(2)  = GNUfields_ls2			;	GNUfields(17) = GNUfields_format_y
					GNUfields(3)  = GNUfields_ls3			;	GNUfields(19) = GNUfields_grid
					GNUfields(4)  = GNUfields_ls4			;	GNUfields(20) = GNUfields_logscale
					GNUfields(5)  = GNUfields_ls5			;	GNUfields(23) = GNUfields_add1
					GNUfields(7)  = GNUfields_title		;	GNUfields(24) = GNUfields_add2
					GNUfields(9)  = GNUfields_xlabel		;	GNUfields(25) = GNUfields_add3
					GNUfields(10) = GNUfields_xrange		;	GNUfields(26) = GNUfields_add4
					GNUfields(11) = GNUfields_mxtics		;	GNUfields(27) = GNUfields_add5
					GNUfields(12) = GNUfields_format_x	;	GNUfields(28) = GNUfields_add6
					GNUfields(14) = GNUfields_ylabel		;	GNUfields(29) = GNUfields_extention
					GNUfields(15) = GNUfields_yrange		;	GNUfields(31) = GNUfields_legend

					GNUfields(33) = GNUfields_plot1		;	GNUfields(47) = GNUfields_plot8
					GNUfields(34) = GNUfields_title1		;	GNUfields(48) = GNUfields_title8
					GNUfields(35) = GNUfields_plot2		;	GNUfields(49) = GNUfields_plot9
					GNUfields(36) = GNUfields_title2		;	GNUfields(50) = GNUfields_title9
					GNUfields(37) = GNUfields_plot3		;	GNUfields(51) = GNUfields_plot10
					GNUfields(38) = GNUfields_title3		;	GNUfields(52) = GNUfields_title10
					GNUfields(39) = GNUfields_plot4		;	GNUfields(53) = GNUfields_plot11
					GNUfields(40) = GNUfields_title4		;	GNUfields(54) = GNUfields_title11
					GNUfields(41) = GNUfields_plot5		;	GNUfields(55) = GNUfields_plot12
					GNUfields(42) = GNUfields_title5		;	GNUfields(56) = GNUfields_title12
					GNUfields(43) = GNUfields_plot6		;	GNUfields(57) = GNUfields_plot13
					GNUfields(44) = GNUfields_title6		;	GNUfields(58) = GNUfields_title13
					GNUfields(45) = GNUfields_plot7		;	GNUfields(59) = GNUfields_plot14
					GNUfields(46) = GNUfields_title7		;	GNUfields(60) = GNUfields_title14

				end subroutine

			subroutine plot_new(datafile) !-=- plot autoscripting
				character(len) datafile

					call new_interface

					GNUdatafile=slashfix(datafile); call file_name(datafile,ii,jj)
					scriptpath = datafile(1:ii-1)//'plots\GNUscripts\'//datafile(ii:jj)//trim(graph_name)//'.plt'

					call system('MD ' // datafile(1:ii-1)//'plots\GNUscripts\' // ' >> log.log ' )

					open(9,file=scriptpath,status='replace')
						if (exptension(GNUfields(29))==1) then
							write(9,*) epsterm
							figoutput = GNUdatafile(1:ii-1)//'plots/GNUscripts/'//GNUdatafile(ii:jj)//trim(graph_name)//'.eps'
							else
								write(9,*) pngterm
								figoutput = GNUdatafile(1:ii-1)//'plots/GNUscripts/'//GNUdatafile(ii:jj)//trim(graph_name)//'.png'
								endif; write(9,*) 'set output "'//trim(figoutput)//'"'; write(9,'(A256)') GNUfields(1:32)

						write(9,*)	'plot "',trim(GNUdatafile),'" u ($1*0):($1*0+10**'//trim(realtostr(GNUnorm))//') w l notitle \'
						write(9,*) (trim(GNUfields(i)),i=33,N_GNUfields) ; close(9);

					figoutput=backslashfix(figoutput)
					call system( 'gnuplot ' // trim(scriptpath) // ' && move ' // trim(figoutput) // ' ' // datafile(1:ii-1)//'plots\' )

					open(8,file=datafile(1:ii-1)//'plots\GNUscripts\clear.bat',status='replace')
						write(8,*) 'taskkill /im wgnuplot*' ; close(8)

					graph_name='' ; GNUnorm=0d0

				end subroutine


