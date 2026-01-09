   function atm_setup_derived_dimensions(readDimensions, dimensionPool, configPool) result(iErr)

      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_log, only : mpas_log_write

      implicit none

      type (mpas_pool_type), intent(inout) :: readDimensions !< Input: Pool to pull read dimensions from
      type (mpas_pool_type), intent(inout) :: configPool !< Input: Pool containing namelist options with configs
      type (mpas_pool_type), intent(inout) :: dimensionPool !< Input/Output: Pool to add dimensions into

      integer :: iErr, errLevel

      integer, pointer :: nCells
      integer, pointer :: nEdges
      integer, pointer :: maxEdges
      integer, pointer :: maxEdges2
      integer, pointer :: nVertices
      integer, pointer :: TWO
      integer, pointer :: THREE
      integer, pointer :: vertexDegree
      integer, pointer :: FIFTEEN
      integer, pointer :: TWENTYONE
      integer, pointer :: R3
      integer, pointer :: nVertLevels
      integer, pointer :: nVertLevelsP1
      integer, pointer :: lon_for
      integer, pointer :: lat_leg
      integer, pointer :: nMonths
      integer, pointer :: nSoilLevels
      integer, pointer :: nLags
      integer, pointer :: nOznLevels
      integer, pointer :: nAerLevels
      integer, pointer :: cam_dim1
      integer, pointer :: nVertLevelsP2
      integer, pointer :: lat
      integer, pointer :: days
      integer, pointer :: nIsoLevelsT
      integer, pointer :: nIsoLevelsZ
      integer, pointer :: t_iso_levels
      integer, pointer :: FIVE
      integer, pointer :: SIXTY
      integer, pointer :: nSoilComps
      integer, pointer :: nSnowLevels
      integer, pointer :: nzSnowLevels
      integer, pointer :: months
      integer, pointer :: num_soil_layers
      integer, pointer :: input_soil_temperature_lag
      integer, pointer :: noznlev
      integer, pointer :: naerlev
      integer, pointer :: camdim1

      iErr = 0
      errLevel = mpas_pool_get_error_level()
      call mpas_pool_set_error_level(MPAS_POOL_SILENT)

      nullify(months)
      call mpas_pool_get_config(configPool, 'months', months)
      nullify(num_soil_layers)
      call mpas_pool_get_config(configPool, 'num_soil_layers', num_soil_layers)
      nullify(input_soil_temperature_lag)
      call mpas_pool_get_config(configPool, 'input_soil_temperature_lag', input_soil_temperature_lag)
      nullify(noznlev)
      call mpas_pool_get_config(configPool, 'noznlev', noznlev)
      nullify(naerlev)
      call mpas_pool_get_config(configPool, 'naerlev', naerlev)
      nullify(camdim1)
      call mpas_pool_get_config(configPool, 'camdim1', camdim1)

      nullify(nCells)
      call mpas_pool_get_dimension(dimensionPool, 'nCells', nCells)
      nullify(nEdges)
      call mpas_pool_get_dimension(dimensionPool, 'nEdges', nEdges)
      nullify(maxEdges)
      call mpas_pool_get_dimension(dimensionPool, 'maxEdges', maxEdges)
      nullify(maxEdges2)
      call mpas_pool_get_dimension(dimensionPool, 'maxEdges2', maxEdges2)
      nullify(nVertices)
      call mpas_pool_get_dimension(dimensionPool, 'nVertices', nVertices)
      nullify(TWO)
      call mpas_pool_get_dimension(dimensionPool, 'TWO', TWO)
      nullify(THREE)
      call mpas_pool_get_dimension(dimensionPool, 'THREE', THREE)
      nullify(vertexDegree)
      call mpas_pool_get_dimension(dimensionPool, 'vertexDegree', vertexDegree)
      nullify(FIFTEEN)
      call mpas_pool_get_dimension(dimensionPool, 'FIFTEEN', FIFTEEN)
      nullify(TWENTYONE)
      call mpas_pool_get_dimension(dimensionPool, 'TWENTYONE', TWENTYONE)
      nullify(R3)
      call mpas_pool_get_dimension(dimensionPool, 'R3', R3)
      nullify(nVertLevels)
      call mpas_pool_get_dimension(dimensionPool, 'nVertLevels', nVertLevels)
      nullify(nVertLevelsP1)
      call mpas_pool_get_dimension(dimensionPool, 'nVertLevelsP1', nVertLevelsP1)
      nullify(lon_for)
      call mpas_pool_get_dimension(dimensionPool, 'lon_for', lon_for)
      nullify(lat_leg)
      call mpas_pool_get_dimension(dimensionPool, 'lat_leg', lat_leg)
      nullify(nMonths)
      call mpas_pool_get_dimension(dimensionPool, 'nMonths', nMonths)
      nullify(nSoilLevels)
      call mpas_pool_get_dimension(dimensionPool, 'nSoilLevels', nSoilLevels)
      nullify(nLags)
      call mpas_pool_get_dimension(dimensionPool, 'nLags', nLags)
      nullify(nOznLevels)
      call mpas_pool_get_dimension(dimensionPool, 'nOznLevels', nOznLevels)
      nullify(nAerLevels)
      call mpas_pool_get_dimension(dimensionPool, 'nAerLevels', nAerLevels)
      nullify(cam_dim1)
      call mpas_pool_get_dimension(dimensionPool, 'cam_dim1', cam_dim1)
      nullify(nVertLevelsP2)
      call mpas_pool_get_dimension(dimensionPool, 'nVertLevelsP2', nVertLevelsP2)
      nullify(lat)
      call mpas_pool_get_dimension(dimensionPool, 'lat', lat)
      nullify(days)
      call mpas_pool_get_dimension(dimensionPool, 'days', days)
      nullify(nIsoLevelsT)
      call mpas_pool_get_dimension(dimensionPool, 'nIsoLevelsT', nIsoLevelsT)
      nullify(nIsoLevelsZ)
      call mpas_pool_get_dimension(dimensionPool, 'nIsoLevelsZ', nIsoLevelsZ)
      nullify(t_iso_levels)
      call mpas_pool_get_dimension(dimensionPool, 't_iso_levels', t_iso_levels)
      nullify(FIVE)
      call mpas_pool_get_dimension(dimensionPool, 'FIVE', FIVE)
      nullify(SIXTY)
      call mpas_pool_get_dimension(dimensionPool, 'SIXTY', SIXTY)
      nullify(nSoilComps)
      call mpas_pool_get_dimension(dimensionPool, 'nSoilComps', nSoilComps)
      nullify(nSnowLevels)
      call mpas_pool_get_dimension(dimensionPool, 'nSnowLevels', nSnowLevels)
      nullify(nzSnowLevels)
      call mpas_pool_get_dimension(dimensionPool, 'nzSnowLevels', nzSnowLevels)

call mpas_log_write('Assigning remaining dimensions from definitions in Registry.xml ...')
      if ( .not. associated(nCells) ) then
         allocate(nCells)
         nCells = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nCells', nCells)
      end if

      if ( .not. associated(nEdges) ) then
         allocate(nEdges)
         nEdges = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nEdges', nEdges)
      end if

      if ( .not. associated(maxEdges) ) then
         allocate(maxEdges)
         maxEdges = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'maxEdges', maxEdges)
      end if

      if ( .not. associated(maxEdges2) ) then
         allocate(maxEdges2)
         maxEdges2 = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'maxEdges2', maxEdges2)
      end if

      if ( .not. associated(nVertices) ) then
         allocate(nVertices)
         nVertices = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nVertices', nVertices)
      end if

      call mpas_pool_get_dimension(dimensionPool, 'TWO', TWO)
      if ( .not. associated(TWO) ) then
         allocate(TWO)
         TWO = 2
call mpas_log_write('       TWO = $i', intArgs=(/2/))
         call mpas_pool_add_dimension(dimensionPool, 'TWO', TWO)
          else if ( TWO == MPAS_MISSING_DIM ) then
         TWO = 2
          end if

      call mpas_pool_get_dimension(dimensionPool, 'THREE', THREE)
      if ( .not. associated(THREE) ) then
         allocate(THREE)
         THREE = 3
call mpas_log_write('       THREE = $i', intArgs=(/3/))
         call mpas_pool_add_dimension(dimensionPool, 'THREE', THREE)
          else if ( THREE == MPAS_MISSING_DIM ) then
         THREE = 3
          end if

      if ( .not. associated(vertexDegree) ) then
         allocate(vertexDegree)
         vertexDegree = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'vertexDegree', vertexDegree)
      end if

      call mpas_pool_get_dimension(dimensionPool, 'FIFTEEN', FIFTEEN)
      if ( .not. associated(FIFTEEN) ) then
         allocate(FIFTEEN)
         FIFTEEN = 15
call mpas_log_write('       FIFTEEN = $i', intArgs=(/15/))
         call mpas_pool_add_dimension(dimensionPool, 'FIFTEEN', FIFTEEN)
          else if ( FIFTEEN == MPAS_MISSING_DIM ) then
         FIFTEEN = 15
          end if

      call mpas_pool_get_dimension(dimensionPool, 'TWENTYONE', TWENTYONE)
      if ( .not. associated(TWENTYONE) ) then
         allocate(TWENTYONE)
         TWENTYONE = 21
call mpas_log_write('       TWENTYONE = $i', intArgs=(/21/))
         call mpas_pool_add_dimension(dimensionPool, 'TWENTYONE', TWENTYONE)
          else if ( TWENTYONE == MPAS_MISSING_DIM ) then
         TWENTYONE = 21
          end if

      call mpas_pool_get_dimension(dimensionPool, 'R3', R3)
      if ( .not. associated(R3) ) then
         allocate(R3)
         R3 = 3
call mpas_log_write('       R3 = $i', intArgs=(/3/))
         call mpas_pool_add_dimension(dimensionPool, 'R3', R3)
          else if ( R3 == MPAS_MISSING_DIM ) then
         R3 = 3
          end if

      if ( .not. associated(nVertLevels) ) then
         allocate(nVertLevels)
         nVertLevels = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nVertLevels', nVertLevels)
      end if

      call mpas_pool_get_dimension(dimensionPool, 'nVertLevelsP1', nVertLevelsP1)
      if ( .not. associated(nVertLevelsP1) ) then
         allocate(nVertLevelsP1)
         nVertLevelsP1 = nVertLevels+1
call mpas_log_write('       nVertLevelsP1 = $i', intArgs=(/nVertLevels+1/))
         call mpas_pool_add_dimension(dimensionPool, 'nVertLevelsP1', nVertLevelsP1)
          else if ( nVertLevelsP1 == MPAS_MISSING_DIM ) then
         nVertLevelsP1 = nVertLevels+1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'lon_for', lon_for)
      if ( .not. associated(lon_for) ) then
         allocate(lon_for)
         lon_for = 504
call mpas_log_write('       lon_for = $i', intArgs=(/504/))
         call mpas_pool_add_dimension(dimensionPool, 'lon_for', lon_for)
          else if ( lon_for == MPAS_MISSING_DIM ) then
         lon_for = 504
          end if

      call mpas_pool_get_dimension(dimensionPool, 'lat_leg', lat_leg)
      if ( .not. associated(lat_leg) ) then
         allocate(lat_leg)
         lat_leg = 248
call mpas_log_write('       lat_leg = $i', intArgs=(/248/))
         call mpas_pool_add_dimension(dimensionPool, 'lat_leg', lat_leg)
          else if ( lat_leg == MPAS_MISSING_DIM ) then
         lat_leg = 248
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nMonths', nMonths)
      if ( .not. associated(nMonths) ) then
         allocate(nMonths)
         nMonths = months
call mpas_log_write('       nMonths = $i (months)', intArgs=(/months/))
         call mpas_pool_add_dimension(dimensionPool, 'nMonths', nMonths)
          else if ( nMonths == MPAS_MISSING_DIM ) then
         nMonths = months
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nSoilLevels', nSoilLevels)
      if ( .not. associated(nSoilLevels) ) then
         allocate(nSoilLevels)
         nSoilLevels = num_soil_layers
call mpas_log_write('       nSoilLevels = $i (num_soil_layers)', intArgs=(/num_soil_layers/))
         call mpas_pool_add_dimension(dimensionPool, 'nSoilLevels', nSoilLevels)
          else if ( nSoilLevels == MPAS_MISSING_DIM ) then
         nSoilLevels = num_soil_layers
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nLags', nLags)
      if ( .not. associated(nLags) ) then
         allocate(nLags)
         nLags = input_soil_temperature_lag
call mpas_log_write('       nLags = $i (input_soil_temperature_lag)', intArgs=(/input_soil_temperature_lag/))
         call mpas_pool_add_dimension(dimensionPool, 'nLags', nLags)
          else if ( nLags == MPAS_MISSING_DIM ) then
         nLags = input_soil_temperature_lag
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nOznLevels', nOznLevels)
      if ( .not. associated(nOznLevels) ) then
         allocate(nOznLevels)
         nOznLevels = noznlev
call mpas_log_write('       nOznLevels = $i (noznlev)', intArgs=(/noznlev/))
         call mpas_pool_add_dimension(dimensionPool, 'nOznLevels', nOznLevels)
          else if ( nOznLevels == MPAS_MISSING_DIM ) then
         nOznLevels = noznlev
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nAerLevels', nAerLevels)
      if ( .not. associated(nAerLevels) ) then
         allocate(nAerLevels)
         nAerLevels = naerlev
call mpas_log_write('       nAerLevels = $i (naerlev)', intArgs=(/naerlev/))
         call mpas_pool_add_dimension(dimensionPool, 'nAerLevels', nAerLevels)
          else if ( nAerLevels == MPAS_MISSING_DIM ) then
         nAerLevels = naerlev
          end if

      call mpas_pool_get_dimension(dimensionPool, 'cam_dim1', cam_dim1)
      if ( .not. associated(cam_dim1) ) then
         allocate(cam_dim1)
         cam_dim1 = camdim1
call mpas_log_write('       cam_dim1 = $i (camdim1)', intArgs=(/camdim1/))
         call mpas_pool_add_dimension(dimensionPool, 'cam_dim1', cam_dim1)
          else if ( cam_dim1 == MPAS_MISSING_DIM ) then
         cam_dim1 = camdim1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nVertLevelsP2', nVertLevelsP2)
      if ( .not. associated(nVertLevelsP2) ) then
         allocate(nVertLevelsP2)
         nVertLevelsP2 = nVertLevels+2
call mpas_log_write('       nVertLevelsP2 = $i', intArgs=(/nVertLevels+2/))
         call mpas_pool_add_dimension(dimensionPool, 'nVertLevelsP2', nVertLevelsP2)
          else if ( nVertLevelsP2 == MPAS_MISSING_DIM ) then
         nVertLevelsP2 = nVertLevels+2
          end if

      if ( .not. associated(lat) ) then
         allocate(lat)
         lat = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'lat', lat)
      end if

      if ( .not. associated(days) ) then
         allocate(days)
         days = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'days', days)
      end if

      call mpas_pool_get_dimension(dimensionPool, 'nIsoLevelsT', nIsoLevelsT)
      if ( .not. associated(nIsoLevelsT) ) then
         allocate(nIsoLevelsT)
         nIsoLevelsT = 22
call mpas_log_write('       nIsoLevelsT = $i', intArgs=(/22/))
         call mpas_pool_add_dimension(dimensionPool, 'nIsoLevelsT', nIsoLevelsT)
          else if ( nIsoLevelsT == MPAS_MISSING_DIM ) then
         nIsoLevelsT = 22
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nIsoLevelsZ', nIsoLevelsZ)
      if ( .not. associated(nIsoLevelsZ) ) then
         allocate(nIsoLevelsZ)
         nIsoLevelsZ = 13
call mpas_log_write('       nIsoLevelsZ = $i', intArgs=(/13/))
         call mpas_pool_add_dimension(dimensionPool, 'nIsoLevelsZ', nIsoLevelsZ)
          else if ( nIsoLevelsZ == MPAS_MISSING_DIM ) then
         nIsoLevelsZ = 13
          end if

      call mpas_pool_get_dimension(dimensionPool, 't_iso_levels', t_iso_levels)
      if ( .not. associated(t_iso_levels) ) then
         allocate(t_iso_levels)
         t_iso_levels = 22
call mpas_log_write('       t_iso_levels = $i', intArgs=(/22/))
         call mpas_pool_add_dimension(dimensionPool, 't_iso_levels', t_iso_levels)
          else if ( t_iso_levels == MPAS_MISSING_DIM ) then
         t_iso_levels = 22
          end if

      call mpas_pool_get_dimension(dimensionPool, 'FIVE', FIVE)
      if ( .not. associated(FIVE) ) then
         allocate(FIVE)
         FIVE = 5
call mpas_log_write('       FIVE = $i', intArgs=(/5/))
         call mpas_pool_add_dimension(dimensionPool, 'FIVE', FIVE)
          else if ( FIVE == MPAS_MISSING_DIM ) then
         FIVE = 5
          end if

      call mpas_pool_get_dimension(dimensionPool, 'SIXTY', SIXTY)
      if ( .not. associated(SIXTY) ) then
         allocate(SIXTY)
         SIXTY = 60
call mpas_log_write('       SIXTY = $i', intArgs=(/60/))
         call mpas_pool_add_dimension(dimensionPool, 'SIXTY', SIXTY)
          else if ( SIXTY == MPAS_MISSING_DIM ) then
         SIXTY = 60
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nSoilComps', nSoilComps)
      if ( .not. associated(nSoilComps) ) then
         allocate(nSoilComps)
         nSoilComps = 8
call mpas_log_write('       nSoilComps = $i', intArgs=(/8/))
         call mpas_pool_add_dimension(dimensionPool, 'nSoilComps', nSoilComps)
          else if ( nSoilComps == MPAS_MISSING_DIM ) then
         nSoilComps = 8
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nSnowLevels', nSnowLevels)
      if ( .not. associated(nSnowLevels) ) then
         allocate(nSnowLevels)
         nSnowLevels = 3
call mpas_log_write('       nSnowLevels = $i', intArgs=(/3/))
         call mpas_pool_add_dimension(dimensionPool, 'nSnowLevels', nSnowLevels)
          else if ( nSnowLevels == MPAS_MISSING_DIM ) then
         nSnowLevels = 3
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nzSnowLevels', nzSnowLevels)
      if ( .not. associated(nzSnowLevels) ) then
         allocate(nzSnowLevels)
         nzSnowLevels = 7
call mpas_log_write('       nzSnowLevels = $i', intArgs=(/7/))
         call mpas_pool_add_dimension(dimensionPool, 'nzSnowLevels', nzSnowLevels)
          else if ( nzSnowLevels == MPAS_MISSING_DIM ) then
         nzSnowLevels = 7
          end if

      call mpas_log_write(' ')
      call mpas_log_write(' ----- done assigning dimensions from Registry.xml -----')
      call mpas_log_write(' ')
      call mpas_log_write(' ')
      call mpas_pool_set_error_level(errLevel)

   end function atm_setup_derived_dimensions


   function atm_setup_decomposed_dimensions(block, manager, readDimensions, dimensionPool, totalBlocks) result(iErr)

      use mpas_derived_types
      use mpas_decomp
      use mpas_pool_routines
      use mpas_io_units
      use mpas_abort, only : mpas_dmpar_global_abort
      use mpas_log, only : mpas_log_write

      implicit none

      type (block_type), intent(inout) :: block !< Input: Pointer to block
      type (mpas_streamManager_type), intent(inout) :: manager !< Input: Stream manager
      type (mpas_pool_type), intent(inout) :: readDimensions !< Input: Pool to pull read dimensions from
      type (mpas_pool_type), intent(inout) :: dimensionPool !< Input/Output: Pool to add dimensions into
      integer, intent(in) :: totalBlocks !< Input: Number of blocks

      integer :: iErr
      type (field1DInteger), pointer :: ownedIndices
      procedure (mpas_decomp_function), pointer :: decompFunc


      iErr = 0
      call mpas_log_write('Processing decomposed dimensions ...')

      call mpas_log_write(' ')
      call mpas_log_write(' ----- done processing decomposed dimensions -----')
      call mpas_log_write(' ')
      call mpas_log_write(' ')

   end function atm_setup_decomposed_dimensions
