#ifdef MPAS_DEBUG
#define COMMA ,
#define PACKAGE_LOGIC_PRINT(M) call mpas_log_write(M)
#else
#define PACKAGE_LOGIC_PRINT(M) ! M
#endif

   !
   ! WARNING: This function is automatically generated at compile time.
   !          Any modifications to this code will be lost when MPAS is recompiled.
   !
   function init_atm_setup_packages_when(configPool, packagePool) result(ierr)

      use mpas_derived_types, only : mpas_pool_type
      use mpas_log, only : mpas_log_write

      implicit none

      integer :: ierr

      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool


      ierr = 0

      call mpas_log_write('')
      call mpas_log_write('Configuring registry-specified packages...')


      call mpas_log_write('----- done configuring registry-specified packages -----')
      call mpas_log_write('')

   end function init_atm_setup_packages_when

