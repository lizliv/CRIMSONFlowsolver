subroutine restar_mpiio (code,  qp, acp)
!
!----------------------------------------------------------------------
!  This routine is the restart option.
!
! input:
!  code                 : restart option on the primitive variables
!                          eq. 'in  ', read from [restar.inp]
!                          eq. 'out ', write to  [restar.out]
!
! input or output:
!  q     (nshg,ndof)   : the variables to be read/written
!
!
! Farzin Shakib, Summer 1985.
!----------------------------------------------------------------------
!
use GlobalArrays          ! used to access acold
use readarrays          ! used to access qold
use phcommonvars  
IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
include "mpif.h"
!
character*4 code
character*8 mach2
character*20 fname1,  fmt1
character*5  cname

!
dimension qp(nshg,ndof),acp(nshg,ndof)
dimension qw(nshg,ndof),acw(nshg,ndof)
!
!.... -------------------------->  'in (unused now since initialization of y and ac moved to readnblk) '  <---------------------------
!
        !if (code .eq. 'in  ') then
!
! incompressible orders velocity, pressure, temperature unlike compressible
! which is what we have our files set up for
!
        !  qp(:,1:3)=qold(:,2:4)
        !  qp(:,4)=qold(:,1)
        !  if(ndof.gt.4)  qp(:,5:ndof)=qold(:,5:ndof)

        !  acp(:,1:3)=acold(:,2:4)
        !  acp(:,4)=acold(:,1)
        !  if(ndof.gt.4)  acp(:,5:ndof)=acold(:,5:ndof)

        !  deallocate(qold)
        !  return
        !endif
!
!.... -------------------------->  'out '  <---------------------------
!
if (code .eq. 'out ') then
        qw(:,2:4) = qp(:,1:3)
        qw(:,1)   = qp(:,4)
        if(ndof.gt.4) qw(:,5:ndof)   = qp(:,5:ndof)
! 
        acw(:,2:4) = acp(:,1:3)
        acw(:,1)   = acp(:,4)
        if(ndof.gt.4) acw(:,5:ndof)   = acp(:,5:ndof)
!           
!       call write_restart(myrank, currentTimestepIndex, nshg, ndof, & 
!         qw, acw)

        ! First, we just send data to a single file from all processors
        call par_write_restart()


        ! Previous numstart update
        if (myrank.eq.master) then 
        open(unit=72,file='numstart.dat',status='old')
        write(72,*) currentTimestepIndex
        close(72)

        endif
!$$$          ttim(75) = ttim(75) + tmr()
        return
endif
!
!.... ---------------------->  Error Handling  <-----------------------
!
!.... Error handling
!
call error ('restar  ',code//'    ',0)
!
!.... file error handling
!
call error ('restar  ','opening ', irstin)
call error ('restar  ','opening ', irstou)
!
!.... end
!
end

subroutine par_write_restart()
        ! Get variables like currentTimeStepIndex, nshg, and ndof
        ! Look at common.f90 to see a list/description of available variables
        use phcommonvars 
        use, intrinsic :: iso_fortran_env
        implicit none
        include "mpif.h"
        
        ! if(rank.eq.0) write(*,"(a)") __FILE__//":"//itoa(__LINE__) ! Useful debugging trick
        
        ! character*4 code
        ! dimension qp(nshg,ndof),acp(nshg,ndof)
        
        integer :: fileno
        integer, dimension(mpi_status_size) :: wstatus
        integer :: ierr, rank, comsize
        
        integer i1,i2
        integer :: magicNumber
        character,parameter :: lf=achar(10) !linefeed aka newline

        character(len=:), allocatable :: fh, sh, message
        ! cannot use literal 0, use 0_MPI_OFFSET_KIND for zero
        integer(mpi_offset_kind) displacement 
        integer(mpi_offset_kind) :: offset, fhoffset, shoffset

        magicNumber = 362436
        ! variables that don't exist in CRIMSON but won't change
        i1=5
        i2=1

        call MPI_Comm_size(MPI_COMM_WORLD, comsize, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
        
        ! NOTE: Should be using the total number of nshg not just nshg
        fh = "# PHASTA Input File Version 2.0"//lf//"# Byte Order Magic Number : "//itoa(magicNumber)//lf
        fh = fh//"number of modes : < "//itoa(0)//" > "//itoa(nshg)//lf
        fh = fh//"number of variables : < "//itoa(0)//" > "//itoa(ndof)//lf
        fh = fh//"byteorder magic number  : < "//itoa(i1)//" > "//itoa(i2)//lf
                
        sh = lf//"solution  : < "//itoa(nshg*ndof*8+1)//" > "//itoa(nshg)//itoa(ndof)//itoa(currentTimeStepIndex)//lf
        
        message = 'Hello from proc '//itoa(rank)//lf
        
        call MPI_File_open(MPI_COMM_WORLD, "helloworld.txt",     &
        ior(MPI_MODE_CREATE,MPI_MODE_WRONLY), &
        MPI_INFO_NULL, fileno, ierr)
        
        ! Set our file view to be in bytes
        displacement = 0_MPI_OFFSET_KIND
        fhoffset = sizeof(fh) ! prints to 101, aka num of characters = num of bytes (1 char = 1 byte)
        shoffset = fhoffset + sizeof(magicNumber)

        call MPI_File_set_view(fileno, displacement, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
        
        if (rank == 0) then
                ! Write the file header
                call MPI_File_write_at(fileno, displacement, fh, len(fh), MPI_CHARACTER, &
                        wstatus, ierr)
                ! Write the byte order magic number
                call MPI_File_write_at(fileno, fhoffset, magicNumber, 1, MPI_INT, &
                        wstatus, ierr)
                ! Write the solution header
                call MPI_File_write_at(fileno, shoffset, sh, len(sh), MPI_CHARACTER, &
                        wstatus, ierr)
        endif
        
        displacement = shoffset + sizeof(sh)
        call MPI_File_set_view(fileno, displacement, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)

        offset = rank*sizeof(message)
        call MPI_File_write_at(fileno, offset, message, len(message), MPI_CHARACTER, &
                wstatus, ierr)
        ! call MPI_File_write(fileno, message, len(message), MPI_CHARACTER, &
        !         wstatus, ierr)

        call MPI_File_close(fileno, ierr)
contains
function itoa(i) result(res)
        implicit none
        character(:), allocatable :: res
        INTEGER,intent(in) :: i
        character(range(i)+2) :: tmp
        write(tmp,'(i0)') i
        res = trim(tmp)
     end function
end