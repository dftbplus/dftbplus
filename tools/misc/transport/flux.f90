program current
  implicit none

  integer, parameter :: dp=8 
  integer :: natoms    
  integer, dimension(:,:), allocatable :: neig     
  integer, dimension(:), allocatable :: nn     
  real(8), dimension(:,:), allocatable :: Inm    
  real(8), dimension(:,:), allocatable :: coord
  integer :: m, tmp, i, j, n, k, maxnn, maxnn_valid
  real(8) :: z, Iz, dr(3), e(3), Imax, frac, width, norm_dr, arr_len
  character(64) :: arg, filename
  character(4) :: id
  logical :: bondflux

  if (iargc().lt.4) then 
     write(*,*) "flux -a|b 'lcurrent.dat' natoms maxneigbors [fraction width]"
     write(*,*) " -a : atom currents; -b: bond currents"
     write(*,*) " natoms: number of atoms in flux file"
     write(*,*) " maxneighbours: neighbours considered in flux calculation"
     write(*,*) " fraction [1.0]: The arrow lengths are normalized to I_max/fraction"
     write(*,*) " width [0.2]: arrows are given a width depending on length"   
     write(*,*) "              arrows too long are made thicker and rescaled"   
     stop 
  endif

  CALL getarg(1, arg)
  if (trim(arg).eq."-b") bondflux = .true.
  if (trim(arg).eq."-a") bondflux = .false.

  CALL getarg(2, filename)

  CALL getarg(3, arg)
  read(arg,*) natoms
  n = 40
  CALL getarg(4, arg)
  read(arg,*) maxnn

  if (iargc().gt.4 ) then
    CALL getarg(5, arg)
    read(arg,*) frac
  else
    frac = 1.0
  endif

  if (iargc().gt.5 ) then
    CALL getarg(6, arg)
    read(arg,*) width
  else
    width = 0.2      
  endif  

  allocate(neig(natoms,n))
  allocate(nn(natoms))
  allocate(Inm(natoms,n))
  allocate(coord(3,natoms))

  Inm = 0.d0

  open(105,file=trim(filename))
  do m=1, natoms    
     read(105,*) tmp, coord(1:3,m), nn(m), (neig(m,i), Inm(m,i), i=1,nn(m)) 
  enddo
  close(105)

  coord=coord*0.529177_dp

  Imax=maxval(abs(Inm(1:natoms,1:maxnn)))
  print*,'# Imax=',Imax

  if(Imax.eq.0.d0) Imax=1.0_dp
  k=1

  ! Bond flux
  if (bondflux) then
    do m = 1, natoms 

       maxnn_valid = min(maxnn,nn(m))

       do j=1, maxnn_valid
       
         ! vector joining atoms
         n=neig(m,j)
         ! vector directed along the bond:
         dr(:)=coord(:,n)-coord(:,m)
         norm_dr = sqrt(dot_product(dr,dr))
         !draw the bond current only if it follows the bond
         arr_len = frac*Inm(m,j)/Imax

         if (arr_len .gt. 0.2d0) then
            write(*,*) '# ',n,m,Inm(m,j)/Imax 
            write(id,'(i4.4)') k
            arg='draw arr'//id//' arrow width'

            if(arr_len .lt. 0.4d0) then   
               e(:)=coord(:,m)+arr_len*dr(:)*2.0
               write(*,'(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)') trim(arg),0.5*width,' {', &
                     & coord(1:3,m),'}{',e(1:3),'} color yellow' 
            endif
            if(arr_len .ge. 0.4 .and. arr_len .lt. 1.d0) then   
               e(:)=coord(:,m)+arr_len*dr(:)
               write(*,'(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)') trim(arg),width,' {', &
                     & coord(1:3,m),'}{',e(1:3),'} color yellow' 
            endif
            if(arr_len .ge. 1.d0 .and. arr_len .lt. 1.5d0) then   
               e(:)=coord(:,m)+dr(:)*0.8
               write(*,'(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)') trim(arg),1.5*width,' {', &
                     & coord(1:3,m),'}{',e(1:3),'} color yellow' 
            endif
            if(arr_len .ge. 1.5d0 .and. arr_len .lt. 2.0d0) then   
               e(:)=coord(:,m)+dr(:)*0.8
               write(*,'(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)') trim(arg),2.0*width,' {', &
                     & coord(1:3,m),'}{',e(1:3),'} color yellow' 
            endif
            if(arr_len .ge. 2.0d0 .and. arr_len .lt. 2.5d0) then   
               e(:)=coord(:,m)+dr(:)*0.8
               write(*,'(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)') trim(arg),2.5*width,' {', &
                     & coord(1:3,m),'}{',e(1:3),'} color yellow' 
            endif
            if(arr_len .ge. 2.5d0 ) then   
               e(:)=coord(:,m)+dr(:)*0.8
               write(*,'(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)') trim(arg),3.0*width,' {', &
                     & coord(1:3,m),'}{',e(1:3),'} color yellow' 
            endif
            k=k+1
         endif     
 
       enddo
    enddo

  else  ! Atomic flux

    do m = 1, natoms 
       Iz =0.0_dp
       e = 0.0_dp
 
       do j=1, maxnn
       
         ! vector joining atoms
         n=neig(m,j)
         ! vector directed along the bond:
         dr(:)= coord(:,n)-coord(:,m)
         e(:) = e(:) + frac*Inm(m,j)*dr(:)/Imax
       
       enddo     
    
       write(id,'(i4.4)') k
       arg='draw arr'//id//' arrow width'
       write(*,'(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)') trim(arg),width,' {', &
           & coord(1:3,m),'}{',coord(:,m)+e(1:3),'} color yellow' 
       k=k+1
 
    enddo

  endif 

  deallocate(neig,nn,Inm,coord)

end program current      
