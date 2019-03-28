module setupgeom
  use accuracy
  use constants
  use message
  use sorting
  use simplealgebra  
  use wrappedIntrinsics
  use linkedlist
  use typeGeometry
  implicit none

  public :: setupGeometry

  contains

  subroutine setupGeometry(geom, iAtInRegion, ContVec, plCutoff, nPLs)
    type(TGeometry), intent(inout) :: geom
    type(wrappedInt1), intent(inout) :: iAtInRegion(:)
    real(dp), intent(inout) :: contVec(:,:)
    real(dp), intent(in) :: plCutoff
    integer, intent(in) :: nPLs(:)

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

    ! 3. resort the second contact PL to be a shifted copy of the first
    call arrangeContactPLs(geom, iAtInRegion, contVec, contDir, nPLs)

    ! 4. Define device PLs
    call definePLs(geom, iAtInRegion, plcutoff, PLlist)
   
    !5. write ordered geometry
    call print_gen(geom, iAtInRegion, PLlist, contDir)


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
  subroutine arrangeContactPLs(geom, iAtInRegion, contVec, contDir, nPLs)
    type(TGeometry), intent(inout) :: geom
    type(wrappedInt1), intent(inout) :: iAtInRegion(:)
    real(dp), intent(inout) :: contVec(:,:)
    integer, intent(in) :: contDir(:)
    integer, intent(in) :: nPLs(:)

    integer :: ii, jj, kk, icont, ncont, PLsize
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
            if (norm2(cross3(vec,uu))< tol .and. norm2(vec-vv)< tol ) then
              exit
            end if  
          end do
          if (bestcross>tol .or. bestdiff>tol) then
             write(sindx,*) data(ii)
             call error("Atom not found for atom "//trim(adjustl(sindx)))
          end if  
          call swap(data(PLsize+ii),data(jj))
        end do  
        end associate
      else if (nPLs(icont)==1) then
        PLsize = size(iAtInRegion(icont)%data)
        contVec(1:3,icont) = contVec(1:3,icont)*real(sign(1,contDir(icont)),dp)
        call reallocate_int(iAtInRegion(icont)%data, PLsize)
        call reallocate_coords(geom%coords, PLsize)
        call reallocate_int(geom%species, PLsize)
        associate(data=>iAtInRegion(icont)%data)
        do ii = 1, PLsize
          data(PLsize+ii) = geom%nAtom + ii !new atoms to the end of structure
          geom%coords(:,data(PLsize+ii)) = geom%coords(:,data(ii)) + contVec(1:3,icont)
          geom%species(data(PLsize+ii)) = geom%species(data(ii))
        end do
        end associate
        geom%nAtom = geom%nAtom + PLsize
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
    type(wrappedInt1), intent(inout) :: iAtInRegion(:)
    
    integer :: icont, ncont, PLsize

    do icont = 1, ncont
      if (allocated(iAtInRegion(icont)%data)) then
          PLsize = size(iAtInRegion(icont)%data)/2  
          print*,'Atoms in contact',icont,':'
          print*,iAtInRegion(icont)%data(1:PLsize)
          print*,iAtInRegion(icont)%data(PLsize+1:2*PLsize)
      endif
      print*
    end do
    print*,'Atoms in device:'
    print*,iAtInRegion(ncont+1)%data
    print*

  end subroutine print_debug

  ! -----------------------------------------------------------------------------------------------| 
  subroutine print_gen(geom, iAtInRegion, PLlist, contDir)
    type(TGeometry), intent(in) :: geom
    type(wrappedInt1), intent(in) :: iAtInRegion(:)
    type(listIntR1), intent(inout) :: PLlist
    integer, intent(in) :: contDir(:)


    integer, allocatable :: atomsInPL(:)
    integer :: ii, jj, kk, icont, ncont, fd1, fd2, fd3, PLsize
    character(10) :: sindx
    
    ncont = size(iAtInRegion)-1    

    open(newunit=fd1, file='processed.gen')
    open(newunit=fd2, file='transport.hsd')
    write(fd2,*) 'Transport{'

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
    if (geom%tPeriodic) then
      write(fd1,*) geom%origin*Bohr__AA
      write(fd1,*) geom%latVecs(:,1)*Bohr__AA
      write(fd1,*) geom%latVecs(:,2)*Bohr__AA
      write(fd1,*) geom%latVecs(:,3)*Bohr__AA
    end if
    close(fd1)
    close(fd2)
  end subroutine print_gen

end module setupgeom


