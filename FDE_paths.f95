module FDE_paths
   use global

		character(len) WorkDir , path_CF2 , path_catalog_CF2 , path_2MRS , path_catalog_2MRS , &
                        uniform_test , fractal_test , cyrle_path

      integer ::     folder_Samples                   = 1   , &
                     folder_NP_add_files              = 2   , &
                     folder_CD_add_files              = 3   , &
                     folder_MD_add_files              = 4   , &
                     folder_Statistics_add_files      = 5   , &
                     folder_Statistics_Report_tables  = 6   , &
                     folder_Catalogue_add_files       = 7   , &
                     folder_Temp_Generators_tests     = 8

		contains

         subroutine gen_paths

            WorkDir           =  'Science\FDE\FDE_Workspace\Main_Workspace\'

            path_CF2          =  trim(dir)//trim(WorkDir)//'Catalogue\CF2.dat'
            path_2MRS         =  trim(dir)//trim(WorkDir)//'Catalogue\2MRS-grouped.dat'

            path_catalog_CF2  =  trim(dir)//trim(WorkDir)//'Catalogue\catalog_CF2.dat'
            path_catalog_2MRS =  trim(dir)//trim(WorkDir)//'Catalogue\catalog_2MRS.dat'

            cyrle_path        =  trim(dir)//trim(WorkDir)//'Samples\cyrle.dat'


            uniform_test      =  trim(dir)//trim(WorkDir)//'Temp\Generators_tests\uni_test.dat'
            fractal_test      =  trim(dir)//trim(WorkDir)//'Temp\Generators_tests\fra_test.dat'

            folders(:)=''
            folders( folder_Samples                   )  = trim(dir)//trim(WorkDir) //'Samples\'
            folders( folder_NP_add_files              )  = trim(dir)//trim(WorkDir) //'NP\add-files\'
            folders( folder_CD_add_files              )  = trim(dir)//trim(WorkDir) //'CD\add-files\'
            folders( folder_MD_add_files              )  = trim(dir)//trim(WorkDir) //'MD\add-files\'

            folders( folder_Statistics_add_files      )  = trim(dir)//trim(WorkDir) //'Statistics\add-files\'
            folders( folder_Statistics_Report_tables  )  = trim(dir)//trim(WorkDir) //'Statistics\Report_tables\'
            folders( folder_Catalogue_add_files       )  = trim(dir)//trim(WorkDir) //'Catalogue\add-files\'
            folders( folder_Temp_Generators_tests     )  = trim(dir)//trim(WorkDir) //'Temp\Generators_tests\'

            end subroutine gen_paths



end module
