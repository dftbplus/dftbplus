!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module dftbp_helpsetupgeom
  use dftbp_accuracy
  use dftbp_globalenv
  use dftbp_constants
  use dftbp_message
  use dftbp_sorting
  use dftbp_simplealgebra  
  use dftbp_wrappedintr
  use dftbp_linkedlist
  use dftbp_typegeometry
  use dftbp_f08math
  use libnegf_vars, only : contactInfo
  implicit none

  public :: setupGeometry

contains

  !> Set up transport structure with contacts
  subroutine setupGeometry(geom, iAtInRegion, contacts, plCutoff, nPLs, printDebug)

    !> Container of the system geometry    
    type(TGeometry), intent(inout) :: geom

    !> Variable size vectors containing atom indices in each region
    type(TWrappedInt1), intent(inout) :: iAtInRegion(:)
    
    !> Contact vector (contact PL repetition)
    type(contactInfo), intent(in),  allocatable :: contacts(:)
    
    !> Slater-Koster cutoff for partitioning
    real(dp), intent(in) :: plCutoff
    
    !> Number of PLs provided
    integer, intent(in) :: nPLs(:)
    
    !> Whether debug infos should be printed
    logical, intent(in) :: printDebug
   

    type(TListIntR1) :: PLlist
    integer, allocatable :: contdir(:)
    integer :: icont, ncont
    real(dp), allocatable :: contVec(:,:)

    if (geom%tfracCoord) then
      call error("setup geometry does not work for fractional coordinates yet")    
    end if

    ! number of contacts
    ncont = size(iAtInRegion)-1  

    ! this is to keep old code as much as possible
    allocate(contVec(4,ncont)) 

    ! get contact directions (positive or negative)
    allocate(contDir(ncont))
    do icont=1,ncont
      contVec(1:3,icont) = contacts(icont)%lattice(1:3)
      contVec(4,icont) = contacts(icont)%shiftAccuracy
      contDir(icont) = get_contdir(contVec(:,icont))
    end do
    
    ! 1. Translate or fold cell
    call TranslateAndFold(geom, .true.)

    ! 2. Set the indices of atoms in the device region
    call assignDeviceAtoms(geom, iAtInRegion)

    ! 3. Sort each contact along the contact direction
    call sortContacts(geom, iAtInRegion, contDir)

    ! Print atom lists
    if (printDebug) then
      call print_debug(iAtInRegion)
    end if

    ! 4. re-sort the second contact PL to be a shifted copy of the first
    call arrangeContactPLs(geom, iAtInRegion, contacts, contVec, contDir, nPLs, plcutoff)

    ! 5. Define device PLs
    call defineDevicePLs(geom, iAtInRegion, plcutoff, contVec, PLlist)
  
    ! 6. write ordered geometry
    call print_gen(geom, contacts, iAtInRegion, PLlist, plcutoff)
 
    write(stdOut,*)
    write(stdOut,*) "Written processed geometry file 'processed.gen'"
    write(stdOut,*) "Written input lines in 'transport.hsd'"

  end subroutine setupGeometry

  ! -----------------------------------------------------------------------------------------------
  subroutine assignDeviceAtoms(geom, iAtInRegion)
    type(TGeometry), intent(in) :: geom
    type(TWrappedInt1), intent(inout) :: iAtInRegion(:)

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
  
  ! -----------------------------------------------------------------------------------------------
  subroutine sortContacts(geom, iAtInRegion, contDir)
    type(TGeometry), intent(in) :: geom
    type(TWrappedInt1), intent(inout) :: iAtInRegion(:)
    integer, intent(inout) :: contDir(:)

    integer :: ii, jj, icont, ncont
    real(dp), allocatable :: subarray(:)
    integer, allocatable :: indxs(:), buffer(:)
    real(dp) :: mean, contRange(3)
    integer :: visitOrder(3)

    ncont = size(iAtInRegion)-1
      
    do icont = 1, ncont
      ! check position of contact w.r.t. device atoms 
      ! and decide sign of contact direction
      associate(data=>iAtInRegion(icont)%data, dir=>contDir(icont))
      mean = sum(geom%coords(dir, data))/real(size(data),dp)
      if (mean < minval(geom%coords(dir,iAtInRegion(ncont+1)%data ))) then
        dir = -dir
      end if   
      end associate

      allocate(subarray(size(iAtInRegion(icont)%data)))
      allocate(indxs(size(subarray)))
      ! identify contact direction wrt device
      if (contDir(icont)>0) then
        subarray(:)=geom%coords(abs(contDir(icont)), iAtInRegion(icont)%data(:))
      else
        subarray(:)=-geom%coords(abs(contDir(icont)), iAtInRegion(icont)%data(:))
      end if
      ! sort contact directions perpendicular to device by largest coord range in that direction
      contRange(:) = abs(maxval(geom%coords(:,iAtInRegion(icont)%data(:)),dim=2)&
          & - minval(geom%coords(:,iAtInRegion(icont)%data(:)),dim=2))
      contRange(abs(contDir(icont))) = huge(1.0)
      call index_heap_sort(visitOrder, contRange)
      do ii = 3, 1, -1
        jj = visitOrder(ii)
        if (jj == abs(contDir(icont))) then
          cycle
        end if
        subarray(:) = subarray * (abs(maxval(geom%coords(jj,iAtInRegion(icont)%data(:))) -&
            & minval(geom%coords(jj,iAtInRegion(icont)%data(:)))) + epsilon(1.0))
        subarray(:) = subarray + geom%coords(jj, iAtInRegion(icont)%data(:))
      end do
      call index_heap_sort(indxs, subarray)
      buffer = iAtInRegion(icont)%data(indxs)
      iAtInRegion(icont)%data = buffer
      deallocate(subarray,indxs,buffer)
    end do
  end subroutine sortContacts
  
  ! -----------------------------------------------------------------------------------------------!
  subroutine arrangeContactPLs(geom, iAtInRegion, contacts, contVec, contDir, nPLs, plcutoff)
    type(TGeometry), intent(inout) :: geom
    type(TWrappedInt1), intent(inout) :: iAtInRegion(:)
    type(contactInfo), intent(in) :: contacts(:)
    real(dp), intent(inout) :: contVec(:,:)
    integer, intent(in) :: contDir(:)
    integer, intent(in) :: nPLs(:)
    real(dp), intent(in) :: plcutoff

    integer :: ii, jj, icont, ncont, PLsize, nAddPLs, iPL
    real(dp) :: vec(3), uu(3), vv(3), tol, bestcross, bestdiff, mindist
    character(lc) :: errmess(4)

    ncont = size(iAtInRegion)-1

    do icont = 1, ncont
      write(stdOut, *) "Contact",icont, '"'//trim(contacts(icont)%name)//'"'
      if (nPLs(icont)==2) then
        associate(data=>iAtInRegion(icont)%data)       
        uu = 0.0_dp
        uu(abs(contDir(icont)))=1.0_dp
        vv = contvec(1:3,icont)
        tol = norm2(vv)*contvec(4,icont)
        PLsize = size(data)/2
        write(stdOut, *) "PL size:",PLsize
        write(stdOut, *) "Number of PLs:",nPLs(icont)
        write(stdOut, *) "contact vector:",contvec(1:3,icont)*Bohr__AA,'(A)'
        write(stdOut, *) "tolerance:",tol
        ! check PL size   
        mindist=minDist2ndPL(geom%coords,data,PLsize,contvec(1:3,icont)) 
        write(stdOut, *) "minimum distance 2nd neighbour PL:", mindist*Bohr__AA,'(A)'
        if (mindist<plcutoff) then
          errmess(1) = "The size of the contact PL is shorter than SK cutoff" 
          errmess(2) = "Check your input geometry or force SKTruncation"
          errmess(3) = "Alternatively, specify 1 PL and let the code adjust contact size"
          errmess(4) = ""
          call error(errmess)  
        end if
        ! 
        write(stdOut, *) "check and reorder contact PLs"
        do ii = 1, PLsize
          bestcross = huge(bestcross)
          bestdiff = huge(bestdiff)
          do jj = PLsize+1, 2*PLsize 
            vec(:) = geom%coords(:,data(jj))-geom%coords(:,data(ii))
            bestcross = min(bestcross, norm2(cross3(vec,uu)))
            bestdiff = min(bestdiff, norm2(vec-vv))
            bestdiff = min(bestdiff, norm2(vec+vv))
            if (norm2(cross3(vec,uu))< tol .and. &
                  & (norm2(vec-vv)< tol .or. norm2(vec+vv)<tol) ) then
              exit
            end if  
          end do
          if (bestcross>tol .or. bestdiff>tol) then
             write(stdOut, *) "Atom ",data(ii)
             write(stdOut, *) "Best cross vector:", bestcross*Bohr__AA, '(A)'
             write(stdOut, *) "Best difference:", bestdiff*Bohr__AA, '(A)'
             call error("Atom not found")
          end if  
          call swap(data(PLsize+ii),data(jj))
        end do  
        end associate
        write(stdOut, *) "contact done"
      else if (nPLs(icont)==1) then
        PLsize = size(iAtInRegion(icont)%data)
        write(stdOut, *) "PL size:",PLsize
        contVec(1:3,icont) = contVec(1:3,icont)*real(sign(1,contDir(icont)),dp)
        write(stdOut, *) "contact vector:",contvec(1:3,icont)*Bohr__AA, '(A)'
        ! counting number of added PLs. Factor of 2 is needed to get 
        ! always an even total number
        mindist=minDist2ndPL(geom%coords,iAtInRegion(icont)%data,PLsize,contvec(1:3,icont))
        nAddPLs = floor(plcutoff/mindist)*2 + 1
        if (nAddPLs==1) then
          write(stdOut, *) "Adding",nAddPLs,"PL"
        else
          write(stdOut, *) "Adding",nAddPLs,"PLs"
          if (nAddPLs>=3) then
            call warning("More than 1 copy of the supplied PL will be added to this contact")
          end if
        end if
        call reallocateInt(iAtInRegion(icont)%data, nAddPLs*PLsize)
        call reallocateCoords(geom%coords, nAddPLs*PLsize)
        call reallocateInt(geom%species, nAddPLs*PLsize)

        ! add atoms to the end of the structure
        do iPL = 1, nAddPLs
          do ii = 1, PLsize
            jj = iPL*PLsize+ii
            iAtInRegion(icont)%data(jj) = geom%nAtom + jj - PLsize
            geom%coords(:,iAtInRegion(icont)%data(jj)) = geom%coords(:,iAtInRegion(icont)%data(ii)) &
                  &+ iPL*contVec(1:3,icont)
            geom%species(iAtInRegion(icont)%data(jj)) = geom%species(iAtInRegion(icont)%data(ii))
          end do
        end do
        geom%nAtom = geom%nAtom + PLsize*nAddPLs
      else
        call error("Number of PLs in contact must be either 1 or 2")
      end if
    end do

    contains
 
    subroutine reallocateInt(array, addsize)
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
    end subroutine reallocateInt

    subroutine reallocateCoords(array, addsize)
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
    end subroutine reallocateCoords

  end subroutine arrangeContactPLs
  
  ! -----------------------------------------------------------------------------------------------
  subroutine defineDevicePLs(geom, iAtInRegion, plcutoff, contVec, PLlist)
    type(TGeometry), intent(in) :: geom
    type(TWrappedInt1), intent(inout) :: iAtInRegion(:)
    real(dp), intent(in) :: plCutoff
    real(dp), intent(in) :: contVec(:,:)
    type(TListIntR1), intent(out) :: PLlist

    integer :: ii, jj, kk, sizeL,  sizeD, ncont
    integer :: icx, icy, icz, nc(3)
    logical, allocatable :: mask(:)
    type(TListInt) :: atomsInPL
    integer, allocatable :: buffer(:)
    real(dp) :: vec(3)
    logical :: addAllR
    
    ! March from contact 1 in the device. Put atoms within cutoff
    ! in a list then restart from that list to the next 
    call init(PLlist)
    buffer = iAtInRegion(1)%data
    sizeL = size(buffer)    
    ncont = size(iAtInRegion)-1
    sizeD = size(iAtInRegion(ncont+1)%data)
    allocate(mask(sizeD))
    mask=.true.
    addAllR=.false.

    nc = 0
    if (geom%tPeriodic .and. allocated(geom%latVecs)) then
      nc = 1
      ! Remove loops on periodic copies along contact direction
      do ii = 1, size(contVec,2)
        nc(maxloc(abs(contvec(:,ii)))) = 0
      end do
    end if
       

    write(stdOut,*) "Partitioning device into PLs"
    ! put array of contact atoms in the first node of PLlist
    associate(dataD=>iAtInRegion(ncont+1)%data, dataR=>iAtInRegion(2)%data)
     
    ! Loop until there are atoms to process      
    do while (count(mask)>0) 
      call init(atomsInPL)
      ! Loop on atoms of current PL
      lpL: do ii = 1, sizeL
        ! Loop on all other atoms including periodic copies 
        do icx = -nc(1), nc(1)
          do icy = -nc(2), nc(2)
            do icz = -nc(3), nc(3)
              do jj = 1, sizeD
                vec(:) = geom%coords(:,buffer(ii))-geom%coords(:,dataD(jj))
                if (allocated(geom%latVecs)) then
                  vec(:) = vec(:)+icx*geom%latVecs(:,1)+icy*geom%latVecs(:,2)+icz*geom%latVecs(:,3)
                end if 
                if (norm2(vec)<plcutoff .and. mask(jj)) then
                  call append(atomsInPL, dataD(jj))
                  mask(jj) = .false.   
                  ! check distance from contact 2. If d < cutoff all remaining atoms
                  ! will go in the last PL
                  do kk = 1, size(dataR)
                    vec(:) = geom%coords(:,dataD(jj))-geom%coords(:,dataR(kk))
                    if (norm2(vec)<plcutoff) then
                      addAllR=.true.
                      exit lpL
                    end if      
                  end do
                end if
              end do 
            end do ! cell z
          end do ! cell y
        end do ! cell x 
      end do lpL
      ! Add all remaining atoms
      if (addAllR) then
        do jj = 1, sizeD
          if (mask(jj)) then
            call append(atomsInPL, dataD(jj))
            mask(jj) = .false.   
          end if 
        end do
      end if
      ! A. transform list to array and add the array to PLlist
      deallocate(buffer)
      allocate(buffer(len(atomsInPL)))
      call asArray(atomsInPL, buffer)
      call destruct(atomsInPL)
      call append(PLlist, buffer) 
      sizeL = size(buffer)
      if (sizeL==0) then
        call error("Found layer of 0 size")
      end if  
      write(stdOut,*) "* Layer size:",sizeL
    end do
    end associate
    write(stdOut,*) "Done."

  end subroutine defineDevicePLs
  ! -----------------------------------------------------------------------------------------------

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

  ! -----------------------------------------------------------------------------------------------
  subroutine swap(a,b)
    integer :: a,b
    integer :: tmp
    tmp=a
    a=b
    b=tmp
  end subroutine swap

  ! -----------------------------------------------------------------------------------------------
  ! debug subroutine
  subroutine print_debug(iAtInRegion)
    type(TWrappedInt1), intent(in) :: iAtInRegion(:)
    
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

  ! -----------------------------------------------------------------------------------------------
  subroutine translateAndFold(geom, tfold)
    type(TGeometry), intent(inout) :: geom
    logical, intent(in) :: tfold


    if (geom%tPeriodic .and. tfold) then
      geom%coords(1,:) = geom%coords(1,:) - minval(geom%coords(1,:))
      geom%coords(2,:) = geom%coords(2,:) - minval(geom%coords(2,:))
      geom%coords(3,:) = geom%coords(3,:) - minval(geom%coords(3,:))
      call foldCoordToUnitCell(geom%coords, geom%latVecs, geom%recVecs2p)
    end if

  end subroutine translateAndFold

        
  ! -----------------------------------------------------------------------------------------------
  subroutine print_gen(geom, contacts, iAtInRegion, PLlist, plCutoff)
    type(TGeometry), intent(in) :: geom
    type(contactInfo), intent(in) :: contacts(:)
    type(TWrappedInt1), intent(in) :: iAtInRegion(:)
    type(TListIntR1), intent(inout) :: PLlist
    real(dp), intent(in) :: plCutoff


    integer, allocatable :: atomsInPL(:)
    character(10) :: sindx
    integer :: ii, jj, kk, icont, ncont, fd1, fd2

    ncont = size(iAtInRegion)-1    

    open(newunit=fd2, file='transport.hsd')
    write(fd2,'(A)') 'Transport{'

    ! Write Device Atoms
    write(fd2,'(2x,A)') 'Device{'
    write(fd2,'(4x,A)', advance='no') 'FirstLayerAtoms={ '
    kk = 0
    do jj = 1, len(PLlist)
      write(sindx,'(I10)') kk+1
      write(fd2,'(A)', advance='no') ' '//trim(adjustl(sindx))
      call get(PLlist, atomsInPL, jj)
      kk = kk + size(atomsInPL)
      deallocate(atomsInPL)
    end do
    write(fd2,*) '}' !close FirstLayerAtoms
    write(sindx,'(I10)') kk
    write(fd2,'(4x,A)') 'AtomRange= 1 '//trim(adjustl(sindx))
    write(fd2,'(2x,A)') '}' !close Device

    ! Write Contact Atoms
    do icont = 1, ncont
      write(fd2,'(2x,A)') 'Contact{'
      write(fd2,'(4x,A)') 'Id = "'//trim(contacts(icont)%name)//'"'
      write(sindx,'(I10)') kk+1
      write(fd2,'(4x,A)',advance='no') 'AtomRange= '//trim(adjustl(sindx))
      kk = kk + size(iAtInRegion(icont)%data)
      write(sindx,'(I10)') kk
      write(fd2,'(A)') ' '//trim(adjustl(sindx))
      write(fd2,'(2x,A)') '}' !close Contact
    end do
    write(fd2,'(A)') '}' !close Transport
    write(fd2,'(A)') 'Hamiltonian = DFTB{'
    write(fd2,'(2x,A)') 'TruncateSKRange = {'
    write(fd2,'(4x,A,F8.4)') 'SKMaxDistance = ', plCutoff
    write(fd2,'(4x,A)') 'HardCutoff = Yes'
    write(fd2,'(2x,A)') '}'
    write(fd2,'(A)') '}'

    close(fd2)

    ! gen output

    open(newunit=fd1, file='processed.gen')
    if (geom%tPeriodic) then
      write(fd1,*) geom%natom, 'S'
    else
      write(fd1,*) geom%natom, 'C'
    endif
    do ii = 1, geom%nSpecies
      write(fd1,'(1x,A,1x)',advance='No') trim(geom%SpeciesNames(ii))
    end do
    write(fd1,*)

102 format(I5,I2,3E20.10)

    ! Write Device Atoms
    write(fd1,"(4X,A)") '# device atoms'
    do jj = 1, len(PLlist)
      kk = 0
      call get(PLlist, atomsInPL, jj)
      do ii = 1, size(atomsInPL)
        kk = kk + 1
        write(fd1,102) kk, geom%species(atomsInPL(ii)), geom%coords(:,atomsInPL(ii))*Bohr__AA
      end do
      deallocate(atomsInPL)
    end do

    ! Write Contact Atoms
    do icont = 1, ncont
      write(fd1,"(4X,A)") '# contact "'//trim(contacts(icont)%name)//'" atoms'
      do ii = 1, size(iAtInRegion(icont)%data)
        kk = mod(ii - 1, size(iAtInRegion(icont)%data) / 2) + 1
        write(fd1,102) kk, geom%species(iAtInRegion(icont)%data(ii)), &
            geom%coords(:,iAtInRegion(icont)%data(ii))*Bohr__AA
      end do
    end do

    if (geom%tPeriodic) then
      write(fd1,*) geom%origin*Bohr__AA
      write(fd1,*) geom%latVecs(:,1)*Bohr__AA
      write(fd1,*) geom%latVecs(:,2)*Bohr__AA
      write(fd1,*) geom%latVecs(:,3)*Bohr__AA
    end if
    close(fd1)

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

  function minDist2ndPL(coord, iAt, PLsize, contvec) result(mindist)
    
    !> Contains the original coordinates on call and the folded ones on return.
    real(dp), intent(in) :: coord(:,:)

    !> Indices of the atoms in a contact
    integer, intent(in) :: iAt(:)

    !> Number of atoms in each contact PL
    integer, intent(in) :: PLsize

    !> Lattice vectors (column format).
    real(dp), intent(in) :: contVec(3)

    !> Output minimal distance between 2nd-neighbour PLs
    real(dp) :: mindist

    integer :: ii, jj
    real(dp) :: vec2(3)

    vec2=2.0_dp*contVec
    mindist = 1e10

    do ii = 1, PLsize
      do jj = 1, PLsize
        mindist = min(mindist,norm2(coord(:,iAt(ii))+vec2-coord(:,iAt(jj))))
      end do   
    end do

  end function minDist2ndPL

end module dftbp_helpsetupgeom
