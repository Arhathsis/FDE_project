!=- fortran-libraries
!=- © Stanislav Shirokov, 2014-2020

module FDE_overleaf
   use global
   use FDE_paths

      character(len) :: overleaf_info , overleaf_CanEdit , overleaf_ReadOnly

   contains

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine overleaf_making

            character(len) filename, new_filename

            overleaf_info = trim(overleaf_dir) // trim(paper_name) // '\'
               call create_folder(overleaf_info)

            forall ( i=1:N_files , files(i)/='NaN' ) files(i) = trim(WorkDir) // trim(files(i))

            do i=1,N_files

               if (files(i)/='NaN') call copy_file( files(i) , overleaf_info )

                      filename = trim(overleaf_info) // trim(name(files(i))) //'.eps'
                  new_filename = trim(overleaf_info) // trim( new_files(i))
               if (files(i)/='NaN') call move_file( filename , new_filename )

               enddo

            overleaf_info = trim(overleaf_info) // trim(paper_name) // '.info'

               unit_1 = random_unit()
            open(unit_1,file=overleaf_info,status='replace')
               write(unit_1,*) 'Anyone with this link can edit this project'
               write(unit_1,*) trim(overleaf_CanEdit)
               write(unit_1,*) 'Anyone with this link can view this project'
               write(unit_1,*) trim(overleaf_ReadOnly)
               write(unit_1,*) 'Local filepath: ', trim(overleaf_info)
               close(unit_1)

            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

         subroutine overleaf_default

            paper_name     = 'NaN'
            new_files(:)   = 'NaN'
            files    (:)   = 'NaN'

            end subroutine

         !=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=%=!

   end module
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- truncated=136-=1
