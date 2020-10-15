!=- Fractal Dimension Estimation
!=- © Stanislav Shirokov, 2014-2020

module FDE_paths
   use global

		character(len) WorkDir , path_CF2 , path_catalog_CF2 , path_2MRS , path_catalog_2MRS , &
                        uniform_test , fractal_test , cyrle_path , scaling_filepath , add_visible_z_out, &
                        FDE_catalog , Millennium , Galacticus , Galacticus_160_dir

      integer ::     folder_Samples                   = 1   , &
                     folder_NP_add_files              = 2   , &
                     folder_CD_add_files              = 3   , &
                     folder_MD_add_files              = 4   , &
                     folder_Statistics_add_files      = 5   , &
                     folder_Statistics_Report_tables  = 6   , &
                     folder_Catalogue_add_files       = 7   , &
                     folder_Samples_add_visible_z     = 8

		contains

         subroutine gen_paths

            !=- external paths

            Millennium = 'C:\Users\Arhath\YandexDisk\Science\DATA\Milli-Millennium\Millennium_z_0_48k.dat'
            Galacticus = 'C:\Users\Arhath\YandexDisk\Science\DATA\CosmoSim\MDPL2_Galacticus\MDPL2_Galacticus_160_10kk.csv'
            Galacticus_160_dir = 'C:\Users\Arhath\YandexDisk\Science\DATA\Ann\cats-20200809T160400Z-001\cats\'

            !=- internal paths

            WorkDir           =  trim( make_workdir() ) // 'Main_Workspace/'   !=- ./

            path_CF2          =  trim(WorkDir)//'Catalogue/CF2.dat'
            path_2MRS         =  trim(WorkDir)//'Catalogue/2MRS-grouped.dat'

            path_catalog_CF2  =  trim(WorkDir)//'Catalogue/catalog_CF2.dat'
            path_catalog_2MRS =  trim(WorkDir)//'Catalogue/catalog_2MRS.dat'

            cyrle_path        =  trim(WorkDir)//'Samples/cyrle.dat'
            scaling_filepath  =  trim(WorkDir)//'Samples/scaling.dat'

            folders(:)=''
            folders( folder_Samples                   )  = trim(WorkDir) //'Samples/'
            folders( folder_NP_add_files              )  = trim(WorkDir) //'NP/add-files/'
            folders( folder_CD_add_files              )  = trim(WorkDir) //'CD/add-files/'
            folders( folder_MD_add_files              )  = trim(WorkDir) //'MD/add-files/'

            folders( folder_Statistics_add_files      )  = trim(WorkDir) //'Statistics/add-files/'
            folders( folder_Statistics_Report_tables  )  = trim(WorkDir) //'Statistics/Report_tables/'
            folders( folder_Catalogue_add_files       )  = trim(WorkDir) //'Catalogue/add-files/'
            folders( folder_Samples_add_visible_z     )  = trim(WorkDir) //'Samples/add_visible_z/'

            end subroutine gen_paths



end module
