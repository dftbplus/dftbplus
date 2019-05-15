module helpsetupgeom
  use dftbp_accuracy
  use dftbp_globalenv
  use dftbp_constants
  use dftbp_message
  use dftbp_sorting
  use dftbp_simplealgebra  
  use dftbp_wrappedIntr
  use dftbp_linkedlist
  use dftbp_typeGeometry
  implicit none

  public :: setupGeometry

  contains

  subroutine setupGeometry(geom, iAtInRegion, ContVec, plCutoff, nPLs, translVec, tfold, &
              & printDebug)
    type(TGeometry), intent(inout) :: geom
    type(wrappedInt1), intent(inout) :: iAtInRegion(:)
    real(dp), intent(inout) :: contVec(:,:)
    real(dp), intent(in) :: plCutoff
    integer, intent(in) :: nPLs(:)
    real(dp), intent(in) :: translVec(:)
    logical, intent(in) :: tfold
    logical, intent(in) :: printDebug
   

    type(listIntR1) :: PLlist
    integer, allocatable :: contdir(:)
    integer :: icont, ncont

    if (geom%tfracCoord) then
      call error("setup geometry does not work for fractional coordinates yet")    
    end if

    ! get contact directions (positive or negative)
    ncont = size(iAtInRegion)-1    
    allocate(contDir(ncont))
    do icont=1,ncont
      contDir(icont) = get_contdir(contVec(:,icont))
    end do

    ! 1. Set the indeces of atoms in the device region
    call assignDeviceAtoms(geom, iAtInRegion)

    ! 2. Sort each contact along the contact direction
    call sortContacts(geom, iAtInRegion, contDir)

    ! Print atom lists
    if (printDebug) then
      call print_debug(geom, iAtInRegion)
    end if

    ! 3. resort the second contact PL to be a shifted copy of the first
    call arrangeContactPLs(geom, iAtInRegion, contVec, contDir, nPLs, plcutoff)

    ! 4. Define device PLs
    call definePLs(geom, iAtInRegion, plcutoff, PLlist)
  
    ! 5. Translate or fold cell
    call TranslateAndFold(geom, translVec, tfold)

    ! 6. write ordered geometry
    call print_gen(geom, iAtInRegion, PLlist, contDir, plcutoff)
 
    write(stdOut,*)
    write(stdOut,*) "Written processed geometry file 'processed.gen'"
    write(stdOut,*) "Written input lines in 'transport.hsd'"

  end subroutine setupGeometry

  ! -----------------------------------------------------------------------------------------------| 
  subroutine assignDeviceAtoms(geom, iAtInRegion)
    type(TGeometry), intent(in) :: geom
    type(wrappedInt1), intent(inout) :: iAtInRegion(:)

    integer :: ii, jj, icont, ncont
    logical, allocatable :: mask(:)
    character(15) :: sindx

    ncont = size(iAtInRegion)-1

    allocate(mask(geom%nAtom))
    mask = .true.
    
    do icont = 1, ncont
       do ii = 1, size(iAtInRegion(icont)%data)
         jj = iAtInRegion(icont)%data(ii)
         if (.not.mask(jj)) then 
           write(sindx,'(I10)') jj    
           call error("atom "//adjustl(sindx)//" is found in more than one region")  
         else 
           mask(jj) = .false.    
         end if
       end do
    end do
  
    jj = count(mask)
    if (allocated(iAtInRegion(ncont+1)%data)) then
      deallocate(iAtInRegion(ncont+1)%data)
    end if 
    allocate(iAtInRegion(ncont+1)%data(jj))
    
    jj = 0 
    do ii = 1, geom%nAtom
      if (mask(ii)) then
        jj = jj + 1
        iAtInRegion(ncont+1)%data(jj) = ii
      end if  
    end do
          
  end subroutine assignDeviceAtoms
  
  ! -----------------------------------------------------------------------------------------------| 
  subroutine sortContacts(geom, iAtInRegion, contDir)
    type(TGeometry), intent(in) :: geom
    type(wrappedInt1), intent(inout) :: iAtInRegion(:)
    integer, intent(inout) :: contDir(:)

    integer :: ii, icont, ncont
    real(dp), allocatable :: subarray(:)
    integer, allocatable :: indxs(:), buffer(:)
    real(dp) :: mean

    ncont = size(iAtInRegion)-1
      
    do icont = 1, ncont
      ! check position of contact w.r.t. device atoms 
      ! and decide sign of contact direction
      mean = 0.0_dp
      associate(data=>iAtInRegion(icont)%data, dir=>contDir(icont))
      do ii = 1, size(data)
        mean = mean + geom%coords(dir, data(ii))
      end do
      if (mean/size(data) < minval(geom%coords(dir,iAtInRegion(ncont+1)%data ))) then
        dir = -dir
      end if   
      end associate

      allocate(subarray(size(iAtInRegion(icont)%data)))
      allocate(indxs(size(subarray)))
      if (contDir(icont)>0) then
        subarray(:)=geom%coords(abs(contDir(icont)), iAtInRegion(icont)%data(:))
      else
        subarray(:)=-geom%coords(abs(contDir(icont)), iAtInRegion(icont)%data(:))
      end if  
      call index_heap_sort(indxs, subarray)
      buffer = iAtInRegion(icont)%data(indxs)
      iAtInRegion(icont)%data = buffer
      deallocate(subarray,indxs,buffer)
    end do
  end subroutine sortContacts
  
  ! -----------------------------------------------------------------------------------------------| 
  subroutine arrangeContactPLs(geom, iAtInRegion, contVec, contDir, nPLs, plcutoff)
    type(TGeometry), intent(inout) :: geom
    type(wrappedInt1), intent(inout) :: iAtInRegion(:)
    real(dp), intent(inout) :: contVec(:,:)
    integer, intent(in) :: contDir(:)
    integer, intent(in) :: nPLs(:)
    real(dp), intent(in) :: plcutoff

    integer :: ii, jj, kk, icont, ncont, PLsize, nAddPLs, iPL
    real(dp) :: vec(3), uu(3), vv(3), tol, bestcross, bestdiff
    character(10) :: sindx


    ncont = size(iAtInRegion)-1

    do icont = 1, ncont
      if (nPLs(icont)==2) then
        associate(data=>iAtInRegion(icont)%data)       
        uu = 0.0_dp
        uu(abs(contDir(icont)))=1.0_dp
        vv = contvec(1:3,icont)
        tol = norm2(vv)*contvec(4,icont)
        PLsize = size(data)/2
        write(stdOut, *) "Contact",icont
        write(stdOut, *) "PL size",PLsize
        write(stdOut, *) "Number of PLs",nPLs(icont)
        write(stdOut, *) "contact vector",contvec(1:3,icont)
        write(stdOut, *) "tolerance",tol
        ! check PL size    
        if (norm2(contvec(1:3,icont))<plcutoff) then
          write(stdOut, *)    
          write(stdOut, *) "The size of the contact PL is shorter than SK cutoff" 
          write(stdOut, *) "Check your input geometry or force SKTruncation"
          write(stdOut, *) "Alternatively, specify 1 PL and let the code adjust contact size"
          call error("")  
        end if
        ! 
        do ii = 1, PLsize
          bestcross = 1e10
          bestdiff = 1e10
          do jj = PLsize+1, 2*PLsize 
            vec(:) = geom%coords(:,data(jj))-geom%coords(:,data(ii))
            if (norm2(cross3(vec,uu))<bestcross) then
               bestcross = norm2(cross3(vec,uu))  
            end if
            if (norm2(vec-vv)<bestdiff) then
               bestdiff = norm2(vec-vv)   
            end if
            if (norm2(vec+vv)<bestdiff) then
               bestdiff = norm2(vec+vv)   
            end if
            if (norm2(cross3(vec,uu))< tol .and. &
                  & (norm2(vec-vv)< tol .or. norm2(vec+vv)<tol) ) then
              exit
            end if  
          end do
          if (bestcross>tol .or. bestdiff>tol) then
             write(stdOut, *) "Atom ",data(ii)   
             write(stdOut, *) "Best cross vector:", bestcross   
             write(stdOut, *) "Best difference:", bestdiff
             call error("Atom not found")
          end if  
          call swap(data(PLsize+ii),data(jj))
        end do  
        end associate
      else if (nPLs(icont)==1) then
        PLsize = size(iAtInRegion(icont)%data)
        write(stdOut, *) "PL size",PLsize
        contVec(1:3,icont) = contVec(1:3,icont)*real(sign(1,contDir(icont)),dp)
        write(stdOut, *) "contact vector",contvec(1:3,icont)
        ! counting number of added PLs. Factor of 2 is needed to get 
        ! always an even total number
        nAddPLs = floor(plcutoff/norm2(contvec(1:3,icont)))*2 + 1
        write(stdOut, *) "Adding",nAddPLs,"PL(s)"
        if (nAddPLs>=3) then
          call warning("More than 1 PL needs to be added") 
        end if
        call reallocate_int(iAtInRegion(icont)%data, nAddPLs*PLsize)
        call reallocate_coords(geom%coords, nAddPLs*PLsize)
        call reallocate_int(geom%species, nAddPLs*PLsize)
        associate(data=>iAtInRegion(icont)%data)
        do iPL = 1, nAddPLs
          do ii = 1, PLsize
            jj = iPL*PLsize+ii
            data(jj) = geom%nAtom + ii !new atoms to the end of structure
            geom%coords(:,data(jj)) = geom%coords(:,data(ii)) + iPL*contVec(1:3,icont)
            geom%species(data(jj)) = geom%species(data(ii))
          end do
        end do
        end associate
        geom%nAtom = geom%nAtom + PLsize*nAddPLs
      else
        call error("Number of PLs in contact must be either 1 or 2")
      end if     
    end do

    contains 
    subroutine reallocate_int(array, addsize)
      integer, intent(inout), allocatable :: array(:)
      integer, intent(in) :: addsize

      integer :: siz
      integer, allocatable :: tmpdata(:) 

      siz = size(array)
      allocate(tmpdata(siz))
      tmpdata=array
      deallocate(array)
      allocate(array(siz+addsize))
      array(1:siz)=tmpdata
      deallocate(tmpdata)
    end subroutine reallocate_int

    subroutine reallocate_coords(array, addsize)
      real(dp), intent(inout), allocatable :: array(:,:)
      integer, intent(in) :: addsize

      integer :: siz
      real(dp), allocatable :: tmpcoords(:,:) 

      siz = size(array,2)
      allocate(tmpcoords(3,siz))
      tmpcoords=array
      deallocate(array)
      allocate(array(3,siz+addsize))
      array(:,1:siz)=tmpcoords
      deallocate(tmpcoords)
    end subroutine reallocate_coords

  end subroutine arrangeContactPLs
  
  ! -----------------------------------------------------------------------------------------------| 
  subroutine definePLs(geom, iAtInRegion, plcutoff, PLlist)
    type(TGeometry), intent(in) :: geom
    type(wrappedInt1), intent(inout) :: iAtInRegion(:)
    real(dp), intent(in) :: plCutoff
    type(listIntR1), intent(out) :: PLlist

    integer :: ii, jj, iPL, sizeL, sizeD, ncont
    logical, allocatable :: mask(:)
    type(listInt) :: atomsInPL
    integer, allocatable :: buffer(:)
    real(dp) :: vec(3)
    
    ! March from contact 1 in the device. Put atoms within cutoff
    ! in a list then restart from that list to the next 
    call init(PLlist)
    buffer = iAtInRegion(1)%data
    ncont = size(iAtInRegion)-1
    sizeL = size(buffer)    
    sizeD = size(iAtInRegion(ncont+1)%data)
    allocate(mask(sizeD))
    mask=.true.
    iPL = 0

    ! put array of contact atoms in the first node of PLlist
    associate(dataD=>iAtInRegion(ncont+1)%data)
    do while (count(mask)>0) 
      call init(atomsInPL)
      do ii = 1, sizeL
        do jj = 1, sizeD
          vec(:) = geom%coords(:,buffer(ii))-geom%coords(:,dataD(jj))
          if (norm2(vec)<plcutoff .and. mask(jj)) then
            call append(atomsInPL, dataD(jj))
            mask(jj) = .false.   
          end if
        end do 
      end do 
      ! A. transform list to array and add the array to PLlist
      deallocate(buffer)
      allocate(buffer(len(atomsInPL)))
      call asArray(atomsInPL, buffer)
      call destruct(atomsInPL)
      call append(PLlist, buffer) 
      sizeL = size(buffer)
      iPL = iPL + 1
    end do   
    end associate

  end subroutine definePLs
  ! -----------------------------------------------------------------------------------------------| 

  function get_contdir(contvec) result(dir)
     real(dp), intent(in) :: contvec(:)
     integer :: dir

     real(dp) :: uu(3), vv(3), tol

     vv = contvec(1:3)
     tol = norm2(vv)*contvec(4)
     dir = 0
     uu=0.0_dp
     uu(1) = 1.0_dp
     if (norm2(cross3(vv,uu))<tol) then
       dir = 1
     end if
     uu(1) = 0.0_dp
     uu(2) = 1.0_dp
     if (norm2(cross3(vv,uu))<tol) then
       dir = 2 
     end if
     uu(2) = 0.0_dp
     uu(3) = 1.0_dp
     if (norm2(cross3(vv,uu))<tol) then
       dir = 3 
     end if
     if (dir == 0) then
        call error("Cannot find contact direction")
     end if   

  end function get_contdir

  ! -----------------------------------------------------------------------------------------------| 
  subroutine swap(a,b)
    integer :: a,b
    integer :: tmp
    tmp=a
    a=b
    b=tmp
  end subroutine swap

  ! -----------------------------------------------------------------------------------------------| 
  ! debug subroutine
  subroutine print_debug(geom, iAtInRegion)
    type(TGeometry), intent(in) :: geom
    type(wrappedInt1), intent(in) :: iAtInRegion(:)
    
    integer :: icont, ncont, PLsize

    ncont = size(iAtInRegion) - 1

    do icont = 1, ncont
      if (allocated(iAtInRegion(icont)%data)) then
        PLsize = size(iAtInRegion(icont)%data)/2  
        write(stdOut,*) 'Atoms in contact',icont,':'
        write(stdOut,*) iAtInRegion(icont)%data(1:PLsize)
        write(stdOut,*) iAtInRegion(icont)%data(PLsize+1:2*PLsize)
        write(stdOut,*) 
      else
        call error("Atom list not allocated")
      end if    
    end do
    write(stdOut,*) 'Atoms in device:'
    write(stdOut,*) iAtInRegion(ncont+1)%data
    write(stdOut,*) 

  end subroutine print_debug

  ! -----------------------------------------------------------------------------------------------| 
  subroutine translateAndFold(geom, translVec, tfold)
    type(TGeometry), intent(inout) :: geom
    real(dp), intent(in) :: translVec(:)
    logical, intent(in) :: tfold

    integer :: ii

    if (geom%tPeriodic .and. tfold) then
      geom%coords(1,:) = geom%coords(1,:) - minval(geom%coords(1,:))
      geom%coords(2,:) = geom%coords(2,:) - minval(geom%coords(2,:))
      geom%coords(3,:) = geom%coords(3,:) - minval(geom%coords(3,:))
      call foldCoordToUnitCell(geom%coords, geom%latVecs, geom%recVecs2p)
    end if

    do ii = 1, geom%nAtom
      geom%coords(:,ii) = geom%coords(:,ii) + translVec 
    end do

  end subroutine translateAndFold

        
  ! -----------------------------------------------------------------------------------------------| 
  subroutine print_gen(geom, iAtInRegion, PLlist, contDir, plCutoff)
    type(TGeometry), intent(in) :: geom
    type(wrappedInt1), intent(in) :: iAtInRegion(:)
    type(listIntR1), intent(inout) :: PLlist
    integer, intent(in) :: contDir(:)
    real(dp), intent(in) :: plCutoff


    integer, allocatable :: atomsInPL(:)
    integer :: ii, jj, kk, icont, ncont, fd1, fd2, fd3, PLsize
    character(10) :: sindx
    
    ncont = size(iAtInRegion)-1    

    open(newunit=fd1, file='processed.gen')
    open(newunit=fd2, file='transport.hsd')
    write(fd2,'(A)') 'Transport{'

    if (geom%tPeriodic) then
      write(fd1,*) geom%natom, 'S'
    else  
      write(fd1,*) geom%natom, 'C'
    endif
    do ii = 1, geom%nSpecies
      write(fd1,'(1x,A,1x)',advance='No') trim(geom%SpeciesNames(ii))
    end do
    write(fd1,*)

    ! Write Device Atoms
    write(fd2,'(2x,A)') 'Device{'
    write(fd2,'(4x,A)', advance='no') 'FirstLayerAtoms={ '
    kk = 0
    do jj = 1, len(PLlist)
      write(sindx,'(I10)') kk+1
      write(fd2,'(A)', advance='no') ' '//trim(adjustl(sindx))
      call get(PLlist, atomsInPL, jj)
      do ii = 1, size(atomsInPL)
        kk = kk + 1
        write(fd1,*) kk, geom%species(atomsInPL(ii)), &
              geom%coords(:,atomsInPL(ii))*Bohr__AA
      end do  
      deallocate(atomsInPL)
    end do  
    write(fd2,*) '}' !close FirstLayerAtoms
    write(sindx,'(I10)') kk
    write(fd2,'(4x,A)') 'AtomRange= 1 '//trim(adjustl(sindx))
    write(fd2,'(2x,A)') '}' !close Device
    ! Write Contact Atoms
    do icont = 1, ncont
      write(fd2,'(2x,A)') 'Contact{'
      write(sindx,'(I10)') kk+1
      write(fd2,'(4x,A)',advance='no') 'AtomRange= '//trim(adjustl(sindx))
      do ii = 1, size(iAtInRegion(icont)%data)
         kk = kk + 1
         write(fd1,*) kk, geom%species(iAtInRegion(icont)%data(ii)), &
            geom%coords(:,iAtInRegion(icont)%data(ii))*Bohr__AA
      end do
      write(sindx,'(I10)') kk
      write(fd2,'(A)') ' '//trim(adjustl(sindx))
      write(fd2,'(2x,A)') '}' !close Contact
    end do
    write(fd2,'(A)') '}' !close Transport
    write(fd2,'(A)') '+Hamiltonian = DFTB{'
    write(fd2,'(2x,A)') '*TruncateSKRange = {'
    write(fd2,'(4x,A,F8.4)') '*SKMaxDistance = ', plCutoff
    write(fd2,'(4x,A)') '*HardCutoff = Yes'
    write(fd2,'(2x,A)') '}'
    write(fd2,'(A)') '}'

    if (geom%tPeriodic) then
      write(fd1,*) geom%origin*Bohr__AA
      write(fd1,*) geom%latVecs(:,1)*Bohr__AA
      write(fd1,*) geom%latVecs(:,2)*Bohr__AA
      write(fd1,*) geom%latVecs(:,3)*Bohr__AA
    end if
    close(fd1)
    close(fd2)
  end subroutine print_gen

  !> Fold coordinates back in the central cell.
  !>
  !> Throw away the integer part of the relative coordinates of every atom. If the resulting
  !> coordinate is very near to 1.0 (closer than 1e-12 in absolute length), fold it to 0.0 to make
  !> the algorithm more predictable and independent of numerical noise.
  subroutine foldCoordToUnitCell(coord, latVec, recVec2p, invShift)

    !> Contains the original coordinates on call and the folded ones on return.
    real(dp), intent(inout) :: coord(:,:)

    !> Lattice vectors (column format).
    real(dp), intent(in) :: latVec(:,:)

    !> Reciprocal vectors in units of 2pi (column format).
    real(dp), intent(in) :: recVec2p(:,:)

    !> Contains difference vectors old_coords - new_coords.
    real(dp), intent(out), optional :: invShift(:,:)


    integer :: nAtom
    integer :: ii, jj
    real(dp) :: frac(3), frac2(3), tmp3(3), vecLen(3)

    nAtom = size(coord, dim=2)

    vecLen(:) = sqrt(sum(latVec(:,:)**2, dim=1))
    do ii = 1, nAtom
      do jj = 1, 3
        frac(jj) = dot_product(recVec2p(:,jj), coord(:,ii))
      end do
      tmp3(:) = coord(:,ii)
      frac2(:) = frac(:) - real(floor(frac(:)), dp)
      where (abs(vecLen*(1.0_dp - frac2)) < 1e-12_dp) frac2 = 0.0_dp
      coord(:, ii) = matmul(latVec, frac2)
      if (present(invShift)) then
        invShift(:,ii) = tmp3(:) - coord(:,ii)
      end if
    end do

  end subroutine foldCoordToUnitCell

end module helpsetupgeom


