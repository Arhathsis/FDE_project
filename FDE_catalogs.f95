module FDE_catalogs
   use FDE_paths
   use FDE_config

   use math
   use Cosmology

   integer,parameter :: N_CF2 = 8188, N_col_CF2 = 40, N_2MRS = 43533, N_col_2MRS = 22

   integer           :: N_col_l, N_col_b, N_col_dl, N_col_mag

   real(8)           :: data_CF2(N_col_CF2,N_CF2)=0d0, data_2MRS(N_col_2MRS,N_2MRS)=0d0, &
                           catalog_CF2(N_col_std,N_CF2)=0d0, catalog_2MRS(N_col_std,N_2MRS)=0d0

   contains

      subroutine catalogs_reading



			!=- The Cosmic Flows 2 catalog

				call reading( N_CF2  , 41 , path_CF2  , data_CF2  ) !=- 8188

					N_col_l = 8 ; N_col_b = 9 ; N_col_dl = 3 ; N_col_mag = 15

				call make_the_catalog( path_catalog_CF2  , data_CF2 , catalog_CF2)



			!=- The 2MASS Redshift Survey catalog

				call reading( N_2MRS , 41 , path_2MRS , data_2MRS ) !=- 43533

					N_col_l = 3 ; N_col_b = 4 ; N_col_dl = 21 ; N_col_mag = 8

					data_2MRS(21,:) = dabs(data_2MRS(21,:)) / H_0
				call make_the_catalog( path_catalog_2MRS , data_2MRS , catalog_2MRS )



         end subroutine



      subroutine reading(N,shift,path,data_catalog)
         character(len) path ; integer N,shift ; real(8) data_catalog(:,:)

         open(N,file=path,status='old')
            do i=1,shift ; read(N,*,end=11) line ; enddo

            do i=1,N
               read(N,*,end=11) data_catalog(:,i)
               enddo
            11 continue ; close(N)

         end subroutine



      subroutine make_the_catalog( catalog_path , data_array , catalog_array )
         character(len) catalog_path
         real(8) data_array(:,:),catalog_array(:,:)

         N_mc = size( data_array(1,:) )

         open( N_mc , file = catalog_path , status = 'replace' )
            write(N_mc,catalog_heads_format  ) (catalog_titles(j),j=1,N_col_std)
            write(N_mc,catalog_heading_format) (catalog_titles(j),j=1,N_col_std)

            do i=1,N_mc

               catalog_array(	N_col_cat_l		,i	) 	= data_array( N_col_l	,i )
               catalog_array(	N_col_cat_b		,i	) 	= data_array( N_col_b	,i )
               catalog_array(	N_col_cat_dl	,i	) 	= data_array( N_col_dl	,i )
               catalog_array(	N_col_cat_mag	,i	)	= data_array( N_col_mag	,i )

               catalog_array(	N_col_cat_Lum	,i	)	= &
                  Mag_to_Lum      ( catalog_array(N_col_cat_mag,i) , catalog_array(N_col_cat_dl,i) ) !=- extinction corrected , 10 pc
               catalog_array(	N_col_cat_M		,i	) 	= &
                  vis_to_abs_mag  ( catalog_array(N_col_cat_mag,i) , catalog_array(N_col_cat_dl,i) ) !=- extinction corrected , 10 pc

               call DescCoord(	&
						catalog_array( N_col_cat_l ,i )	,	catalog_array( N_col_cat_b ,i )	,	catalog_array( N_col_cat_dl,i )	, &   !=- Mpc
                  catalog_array( N_col_cat_x ,i )	,	catalog_array( N_col_cat_y ,i )	,	catalog_array( N_col_cat_z ,i )	)

					if (data_array( N_col_rs	,i )==0) &
						catalog_array( N_col_cat_rs	,i ) = fun_from_R_to_z(catalog_array(	N_col_cat_dl	,i	))	!=- default = 0
					if (data_array( N_col_R		,i )==0) &
						catalog_array( N_col_cat_dl	,i ) = fun_from_z_to_R(catalog_array(	N_col_cat_rs	,i	))	!=- default = 0

               write(N_mc,catalog_format) (catalog_array(j,i),j=1,N_col_std)

               enddo
            close(N_mc)

         end subroutine

end module
