!=- Fractal Dimension Estimation
!=- © Stanislav Shirokov, 2014-2020

module FDE_config
   use global

   logical ::  sets_regenerating = .false. , &  !=- .true. | .false.
               FDE_recalculating = .false.  , &
               FDE_replot        = .false.      !=- have not been written yet

   integer ::  luminosity_model = 0  !=- 0 is uniform geometry (default), 1 is functions, 2 is ...

   real(8) ::  minimal_luminosity   =  1d36  , &   !=- ergs
               maximal_luminosity   =  1d44  , &   !=- ergs
               radius_limit         =  2d2   , &   !=- Mpc
               latitude_limit       =  5     , &   !=- galactic belt in latitude degrees
               longitude_limit      =  -1          !=- the covering in longitude degrees, without = -1

   integer,parameter :: N_col_std = 10, N_col_cat_x = 1, N_col_cat_y = 2, N_col_cat_z = 3, N_col_cat_l = 4, &
		N_col_cat_b = 5, N_col_cat_dl = 6, N_col_cat_Lum = 7, N_col_cat_mag = 8, N_col_cat_M = 9, N_col_cat_rs = 10

   character(len) :: catalog_titles(N_col_std)='',catalog_heading_format='',catalog_format='',catalog_heads_format=''

   contains

      subroutine the_catalog_heading

         catalog_titles	(	N_col_cat_x		) 	= 'x, Mpc'
         catalog_titles	(	N_col_cat_y		) 	= 'y, Mpc'
         catalog_titles	(	N_col_cat_z		) 	= 'z, Mpc'
         catalog_titles	(	N_col_cat_l		) 	= 'l, deg'
         catalog_titles	(	N_col_cat_b		) 	= 'b, deg'
         catalog_titles	(	N_col_cat_dl	) 	= 'dl, Mpc'
         catalog_titles	(	N_col_cat_Lum	) 	= 'L, erg'
         catalog_titles	(	N_col_cat_mag	)	= 'm, mag'
         catalog_titles	(	N_col_cat_M		) 	= 'M, mag'
         catalog_titles	(	N_col_cat_rs	)	= 'redshift'

         catalog_format          = '('//trim(inttostr(N_col_std))//'(E16.8))'
         do i=1,N_col_std
            catalog_titles(i) 	= ' # '//trim(inttostrf(i,2))//' '//adjustl(catalog_titles(i))
            enddo
         catalog_heads_format    = '(A16)'
         catalog_heading_format  = '('//trim(inttostr(N_col_std))//'(A16))'

         end subroutine

end module
