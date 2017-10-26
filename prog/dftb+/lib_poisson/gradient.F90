module gradient


real(8), allocatable, save :: gr(:,:)     !3,NNDIM
real(8), allocatable, save :: hgrad(:,:)  !3,3*NNDIM
integer, allocatable, save :: conat(:)    !NNDIM+1
real(8), allocatable, save :: convec(:,:) !3,NNDIM
logical, save :: constr


end module gradient
