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
        call par_write_restart(qp,acp)

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

! ltg(nshg)       : the local to global node mapping
! ise(2,numpe)    : global start and end for each process
! sendsize(numpe) : the number of nodes to send to each process
! sendbuf         : this is a derived type. The index is rank of the 
!                   process to send the data to. Inside the derived type,
!                   there is 2 arrays to be sent.
! recvbuf         : similar to sendbuf. Data received from other processes. 
!                   first index gives the global node number and
!                   the second index gives the value
! sendcount(numpe): 
!      the number of nodes to send to each process
!
! send_array_size(numpe):
!      the actual size of the sendbuf array. This is different
!      from the actual number of points sent. The number of points being to any process
!      is not known apriori. a simple resiable array is implemented to deal with this.
!      This variable stores the actual size of the array for each rank
!
! send_dest(nshg) : 
!      the rank of the process to which the node is to be sent
!
! shadow(nshg)    :
!      a logical array to indicate whether the node is a shadow node or not
!
! jcount(numpe)   :
!   loop counter helper
subroutine par_write_restart(qp,acp)
        ! Get variables like currentTimeStepIndex, nshg, and ndof
        ! Look at common.f90 to see a list/description of available variables
        use phcommonvars 
        use, intrinsic :: iso_fortran_env, only: wp => real64
        use globalArrays, only: ilwork
        implicit none
        include "mpif.h"
        
        ! if(rank.eq.0) write(*,"(a)") __FILE__//":"//itoa(__LINE__) ! Useful debugging trick
        

        integer, dimension(mpi_status_size) :: wstatus
        integer :: ierr, rank, comsize

        ! nshg will equal nno
        ! 
        real(wp) qp(nshg,ndof), acp(nshg,ndof)
        real(wp), allocatable :: qw(:,:),acw(:,:)
        integer,save :: gnno,nno
        
        integer avgnp,temp,sendto,extranp
        integer global_node_id, new_local_node_id

        ! arrays that need to be saved
        integer, allocatable,save :: ise(:,:),send_count(:),send_dest(:),send_array_size(:),jcount(:)
      
        integer, allocatable :: ltg(:),isize(:)
      
        logical, allocatable,save :: shadow(:)
      
        integer numtask,itkbeg,itask,iacc,is,numseg,isgbeg,lenseg,isgend
        integer itag,nrecv
      
        integer istatus(MPI_STATUS_SIZE,numpe),ireq(numpe),jstatus(numpe)
        logical comm_finished(numpe),flag
      
        type mytype
          integer,  dimension(:),   pointer  :: global_node_id
          real(wp), dimension(:,:), pointer  :: solution
        endtype mytype
      
        type(mytype), allocatable,save :: sendbuf(:), recvbuf(:)
      
        integer i,j,k
        logical,save :: first=.true.

        if(first) call setup_communications

        call do_communications
        call write_restart_serial

        call write_restart_parallel

        ! call write_pht
        
contains
        ! Utility function to convert integer to string
        function itoa(i) result(res)
                implicit none
                character(:), allocatable :: res
                INTEGER,intent(in) :: i
                character(range(i)+2) :: tmp
                write(tmp,'(i0)') i
                res = trim(tmp)
        end function

         ! getproc returns the processor rank (0-based) that should own a given node
         ! based on the node's global number for writing the restart file
        function getproc(node_number) result (proc)
                integer node_number,proc
                integer pc
                proc = -99
                do pc=1,numpe
                   if((node_number.ge.ise(1,pc)) .and. (node_number.le.ise(2,pc)) )then
                      proc = pc-1
                      return
                   endif
                enddo
                write(*,*)'error in getproc on rank '//itoa(myrank)
                stop
        endfunction getproc
          
        ! arr_rsz_int is a subroutine to resize an integer array
        subroutine arr_rsz_int(array,n)
                implicit none
                integer, intent(in)  :: n
                integer, dimension(:), pointer :: array
                integer, dimension(:), pointer :: array_temp
                integer :: n_now,n_new,i
                real(WP), parameter :: coeff_up   = 1.3_WP
                real(WP), parameter :: coeff_down = 0.7_WP
                real(WP) :: up,down
          
                ! Resize array array to size n
                if (.NOT.associated(array)) then
                   ! array is of size 0
                   if (n.eq.0) then
                      ! Nothing to do, that's what we want
                   else
                      ! Allocate directly of size n
                      allocate(array(n))
                      array(1:n)=0
                   end if
                else if (n.eq.0) then
                   ! array is associated, be we want to empty it
                   deallocate(array)
                   nullify(array)
                else
                   ! Update non zero size to another non zero size
                   n_now = size(array,dim=1)
                   up  =real(n_now,WP)*coeff_up
                   down=real(n_now,WP)*coeff_down
                   if (n.gt.n_now) then
                      ! Increase from n_now to n_new
                      n_new = max(n,int(up))
                      allocate(array_temp(n_new))
                      do i=1,n_now
                         array_temp(i) = array(i)
                      end do
                      deallocate(array)
                      nullify(array)
                      array => array_temp
                      array(n_now+1:n_new)=-99
                   else if (n.lt.int(down)) then
                      ! Decrease from n_now to n_new
                      allocate(array_temp(n))
                      do i=1,n
                         array_temp(i) = array(i)
                      end do
                      deallocate(array)
                      nullify(array)
                      array => array_temp
                   end if
                end if
          
                return
        endsubroutine arr_rsz_int
          
        ! arr_rsz_5xdble is a subroutine to resize a (ndof,:) 64bit float  array
        subroutine arr_rsz_5xdble(array,n)
                implicit none
                integer, intent(in)  :: n
                real(wp), dimension(:,:), pointer :: array
                real(wp), dimension(:,:), pointer :: array_temp
                integer :: n_now,n_new,i
                real(WP), parameter :: coeff_up   = 1.3_WP
                real(WP), parameter :: coeff_down = 0.7_WP
                real(WP) :: up,down
          
                ! Resize array array to size n
                if (.NOT.associated(array)) then
                   ! array is of size 0
                   if (n.eq.0) then
                      ! Nothing to do, that's what we want
                   else
                      ! Allocate directly of size n
                      allocate(array(ndof,n))
                      array(:,1:n)=0.0_WP
                   end if
                else if (n.eq.0) then
                   ! array is associated, be we want to empty it
                   deallocate(array)
                   nullify(array)
                else
                   ! Update non zero size to another non zero size
                   n_now = size(array,dim=1)
                   up  =real(n_now,WP)*coeff_up
                   down=real(n_now,WP)*coeff_down
                   if (n.gt.n_now) then
                      ! Increase from n_now to n_new
                      n_new = max(n,int(up))
                      allocate(array_temp(ndof,n_new))
                      do i=1,n_now
                         array_temp(:,i) = array(:,i)
                      end do
                      deallocate(array)
                      nullify(array)
                      array => array_temp
                      array(:,n_now+1:n_new)=0.0_WP
                   else if (n.lt.int(down)) then
                      ! Decrease from n_now to n_new
                      allocate(array_temp(ndof,n))
                      do i=1,n
                         array_temp(:,i) = array(:,i)
                      end do
                      deallocate(array)
                      nullify(array)
                      array => array_temp
                   end if
                end if
          
                return
        endsubroutine arr_rsz_5xdble


        ! plan is to have each process write a contiguous block of the file.
        ! for n processes, i'th process writes from (i-1)*nshg/numpe+1 to i*nshg/numpe. 
        ! note nshg here is global total (and not the per process nshg)
        ! this is done by communicating the data to each process from all other processes
        ! we use the local to global map (ltg.dat.rank) to determine which nodes to send to each process
      
        !   if(myrank .eq. 0) write(*,'(a)')'Entering write_restart_parallel'
        subroutine setup_communications()
                implicit none
                include "mpif.h"

                ! stuff for writing the file
                character(len=256) filename,str
      
                ireq = MPI_REQUEST_NULL
                comm_finished = .false.
        
                comm_finished(myrank+1) = .true.! this process does not need to communicate with itself
                
                ! this part only needs to be done once per simulation. but some of the data needs to persist!
                ! if(first)then
                first=.false.
                write(filename,*) myrank
                filename = "ltg.dat."//adjustl(trim(filename))
                open(1,file=filename)
                read(1,*) gnno ! global number of nodes in the mesh (all processes) i.e. not nshg which is nodes per process
                read(1,*) nno  ! number of nodes on this process
                if (allocated(ltg)) then
                deallocate(ltg)
                endif
                allocate(ltg(nno))
                read(1,*) ltg
                close(1)
                
                allocate(ise(2,numpe),isize(numpe))
                ise = 0
                isize = 0
                
                avgnp = int(gnno/numpe)
                extranp = gnno - avgnp*numpe
                
                do i=1,numpe
                isize(i) = avgnp
                if(i.le.extranp)isize(i) = isize(i)+1
                enddo
        
                ! ise(1,i) is the starting node number for process i
                ! ise(2,i) is the ending node number for process i
                ! inclusive of both ends
                
                ise(1,1) = 1
                ise(2,1) = ise(1,1) + isize(1) - 1
                do i=2,numpe
                ise(1,i) = ise(1,i-1) + isize(i-1)
                ise(2,i) = ise(1,i) + isize(i) - 1
                enddo
        
                !  call DPR1("is: "//itoa(ise(1,myrank+1))//", ie: "//&
                !  itoa(ise(2,myrank+1))//", isize: "//itoa(isize(myrank+1)))
                
        
                !identify shadowed nodes
                allocate(shadow(nno))
        
                shadow = .false.
                numtask = ilwork(1)
                itkbeg = 1
                
                do itask = 1, numtask
                
                iacc   = ilwork (itkbeg + 2)
                numseg = ilwork (itkbeg + 4)
                
                if (iacc .eq. 0) then
                        do is = 1,numseg
                        isgbeg = ilwork (itkbeg + 3 + 2*is)
                        lenseg = ilwork (itkbeg + 4 + 2*is)
                        isgend = isgbeg + lenseg - 1
                        shadow(isgbeg:isgend) = .TRUE.
                        enddo
                endif
                
                itkbeg = itkbeg + 4 + 2*numseg
                
                enddo
        
                ! sanity check. count the number of non-shadow nodes and 
                ! make sure the sum is equal to the global number of nodes
                j=0
                do i=1,nno
                if (.not. shadow(i)) j=j+1
                enddo
                call MPI_ALLREDUCE(j,i,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
                if (i .ne. gnno) call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
        
        
        
                allocate(send_dest(nno))
                send_dest(:) = -99
        
                allocate(jcount(numpe))
        
                ! allocate sendbuf
                allocate(sendbuf(numpe),recvbuf(numpe))
                allocate(send_array_size(numpe),send_count(numpe))
        
                send_array_size(:) = gnno/numpe
                send_count(:) = 0
        
                do i=1,numpe
                allocate(sendbuf(i)%global_node_id(send_array_size(i)))
                allocate(sendbuf(i)%solution(ndof,send_array_size(i)))
                sendbuf(i)%global_node_id(:) = -99
                sendbuf(i)%solution(:,:) = 0.0_wp
                enddo
        
                ! prepare a send array:
                ! need to know what to send and where
                do i=1,nno
                if (.not. shadow(i))then
        
                global_node_id = ltg(i)
                
                !identify which node to send to. plus 1 to address arrays correctly
                sendto = 1 + getproc(global_node_id)
                
                if(sendto .lt. 0)then
                        write(*,'(a)')'Error in write_restart_parallel: node not found'
                        stop
                endif
        
        
                send_dest(i) = sendto
                
                ! count the number of nodes to this process
                send_count(sendto) = send_count(sendto) + 1
                if(send_count(sendto) .gt. send_array_size(sendto))then
                
                        call arr_rsz_int   (sendbuf(sendto)%global_node_id, send_count(sendto))
                        call arr_rsz_5xdble(sendbuf(sendto)%solution,       send_count(sendto))
                
                        send_array_size(sendto) = size(sendbuf(sendto)%global_node_id,DIM=1) ! muy importante
                endif
                
                ! add the global node number to the send array
                sendbuf(sendto)%global_node_id(send_count(sendto)) = global_node_id
        
                endif ! shadow
                enddo   ! i
        
        
                ! do the communication now of global_node_ids now
                ! first part: non-blocking send to all processes
                do i=1,numpe
                if((i-1) .ne. myrank)then
                itag = myrank
                !global node id only needs to be communicated once per run as it does not change with time
                call MPI_ISEND(sendbuf(i)%global_node_id, send_count(i), MPI_INTEGER, (i-1), itag, MPI_COMM_WORLD, ireq(i), ierr)
                endif
                enddo
                
                !second part: probe and receive
                do
                flag = .false.
                call MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, flag, jstatus, ierr)
                if(flag)then
                i=jstatus(MPI_SOURCE)+1
                itag = jstatus(MPI_TAG)
                
                call MPI_GET_COUNT(jstatus, MPI_INTEGER, nrecv, ierr)
                allocate(recvbuf(i)%global_node_id(nrecv))
                call MPI_RECV(recvbuf(i)%global_node_id, nrecv, MPI_INTEGER, (i-1), itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                comm_finished(i) = .true.
                endif
                if(ALL(comm_finished)) exit
                enddo
                
                ! wait for all sends to finish
                call MPI_WAITALL(numpe, ireq, MPI_STATUSES_IGNORE, ierr)
        
                ! endif ! first

        endsubroutine setup_communications

        ! node number communication is done. Now we move on to the solution data. same procedure as above
        ! prepare the data arrays for the send
        subroutine do_communications()
                implicit none
                include "mpif.h"

                jcount(:) = 1
                do i=1,nno
                if (.not. shadow(i))then
                j = jcount(send_dest(i))
                sendbuf(send_dest(i))%solution(:,j) = qp(i,:)
                jcount(send_dest(i)) = jcount(send_dest(i)) + 1
                endif
                enddo

                comm_finished(:) = .false.
                comm_finished(myrank+1) = .true.

                do i=1,numpe
                if((i-1) .ne. myrank)then
                itag = myrank
                !global node id only needs to be communicated once per run as it does not change with time
                call MPI_ISEND(sendbuf(i)%solution, ndof*send_count(i), MPI_DOUBLE_PRECISION, (i-1), itag, MPI_COMM_WORLD, ireq(i), ierr)
                endif
                enddo
        
                ! second part: probe and receive
                do
                flag = .false.
                call MPI_IPROBE(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, flag, jstatus, ierr)
                if(flag)then
                i=jstatus(MPI_SOURCE)+1
                itag = jstatus(MPI_TAG)
        
                call MPI_GET_COUNT(jstatus, MPI_DOUBLE_PRECISION, nrecv, ierr)
                allocate(recvbuf(i)%solution(ndof,nrecv))
                call MPI_RECV(recvbuf(i)%solution, ndof*nrecv, MPI_DOUBLE_PRECISION, (i-1), itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                comm_finished(i) = .true.
                endif
        
                if(ALL(comm_finished)) exit
                enddo
        
                call MPI_WAITALL(numpe, ireq, MPI_STATUSES_IGNORE, ierr)
        
                ! TODO: send the acceleration data
                ! now we have all the data in the recvbuf arrays. we need to put it in the right place
        
                allocate(qw(ise(1,myrank+1):ise(2,myrank+1),ndof))
        
                do i=1,numpe
                if((i-1) .ne. myrank)then
                !check size of recvbuf(i)%global_node_id and loop over it
                k = size(recvbuf(i)%global_node_id,DIM=1)
                ! call DPR1('received '//itoa(k)//' nodes')
                if (k>0)then
                do j=1,k
                        global_node_id = recvbuf(i)%global_node_id(j)
                        qw(global_node_id,1) = recvbuf(i)%solution(4,j)
                        qw(global_node_id,2) = recvbuf(i)%solution(1,j)
                        qw(global_node_id,3) = recvbuf(i)%solution(2,j)
                        qw(global_node_id,4) = recvbuf(i)%solution(3,j)
                        qw(global_node_id,5) = recvbuf(i)%solution(5,j)
                enddo
                endif
                endif
                enddo
                
                !now deal with this process's own data
                k=send_count(myrank+1)
                !   call DPR1('self: '//itoa(k)//' nodes')
                do j=1,k
        
                global_node_id = sendbuf(myrank+1)%global_node_id(j)
        
                ! check if global_node_id is in the range of this process
                if(global_node_id.lt.ise(1,myrank+1).or.global_node_id.gt.ise(2,myrank+1))then
                write(*,'(a)')"k"//itoa(k)//', global_node_id : '//itoa(global_node_id)//' is not in the range of this process'
                call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                endif
                qw(global_node_id,1) = sendbuf(myrank+1)%solution(4,j)
                qw(global_node_id,2) = sendbuf(myrank+1)%solution(1,j)
                qw(global_node_id,3) = sendbuf(myrank+1)%solution(2,j)
                qw(global_node_id,4) = sendbuf(myrank+1)%solution(3,j)
                qw(global_node_id,5) = sendbuf(myrank+1)%solution(5,j)
                enddo
        endsubroutine do_communications

        ! Hardcoded writing to restart file in serial order
        subroutine write_restart_serial()
                implicit none
                include "mpif.h"

                integer fh
                character,parameter :: lf=achar(10) !linefeed aka newline
                ! stuff for writing the file
                character(len=256) filename,str

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !! we have the data where we want it. now we can write it to file 
                !! as first pass, I will write it in serial in rank order. 
                !! Then we can do it using MPI_IO in one shot
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              
                fh = 20
                write(filename, '(a, i0, a)') 'restart_ser.', currentTimestepIndex, '.0'
                if(myrank .eq. 0)then
                   ! root writes the header
                   ! open the file
                   open(unit=fh, file=trim(filename), status='new', action='write', &
                   form='unformatted', access='stream')
                   write(fh) "# PHASTA Input File Version 2.0"//lf
                   write(fh) "# Byte Order Magic Number : 362436"//lf
                !    write(fh) "# Output Generated by CRIMSON"//lf
                   
                   write(str,'(a,i0,a,i0,x,a)')"byteorder magic number  : < ",5," > ",1,lf
                   write(fh)trim(str)
                   write(fh) 362436
                   write(fh)lf
                   
                   str=''
                   write(str,'(a,i0,a,i0,a)')"number of modes : < ",0," > ",gnno,lf
                   write(fh)trim(str)
                   
                   str=''
                   write(str,'(a,i0,a,i0,a)')"number of variables : < ",0," > ",ndof,lf
                   write(fh)trim(str)
                   
                   str=''
                   write(str,'(a,i0,a,3(x,i0),a)')"solution  : < ",(gnno*ndof*8+1)," > ",gnno,ndof,currenttimestepindex,lf
                   write(fh)trim(str)
                   close(fh)
                 endif
              
              
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
              
              
                !now write the data variable by variable in rank order
                ! Done separately for each DOF
                !1
                if(myrank .gt. 0)then
                  call MPI_RECV(temp, 1, MPI_INTEGER, myrank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                endif
              
                open(unit=fh, file=trim(filename), status='old', action='write', &
                 form='unformatted', access='stream', position='append')
                write(fh)qw(:,1)
                close(fh)
              
                if(myrank .lt. (numpe-1))then
                  call MPI_SEND(temp, 1, MPI_INTEGER, myrank+1, 0, MPI_COMM_WORLD, ierr)
                endif
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
              
                !repeat as needed for the other variables
                !2
                if(myrank .gt. 0)then
                 call MPI_RECV(temp, 1, MPI_INTEGER, myrank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                endif
              
                open(unit=fh, file=trim(filename), status='old', action='write', &
                form='unformatted', access='stream', position='append')
                write(fh)qw(:,2)
                close(fh)
                
                if(myrank .lt. (numpe-1))then
                   call MPI_SEND(temp, 1, MPI_INTEGER, myrank+1, 0, MPI_COMM_WORLD, ierr)
                endif
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
              
                !3
                if(myrank .gt. 0)then
                  call MPI_RECV(temp, 1, MPI_INTEGER, myrank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                endif
              
                open(unit=fh, file=trim(filename), status='old', action='write', &
                 form='unformatted', access='stream', position='append')
                write(fh)qw(:,3)
                close(fh)
              
                if(myrank .lt. (numpe-1))then
                  call MPI_SEND(temp, 1, MPI_INTEGER, myrank+1, 0, MPI_COMM_WORLD, ierr)
                endif
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                
                !4
                if(myrank .gt. 0)then
                 call MPI_RECV(temp, 1, MPI_INTEGER, myrank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                endif
              
                open(unit=fh, file=trim(filename), status='old', action='write', &
                form='unformatted', access='stream', position='append')
                write(fh)qw(:,4)
                close(fh)
              
                if(myrank .lt. (numpe-1))then
                 call MPI_SEND(temp, 1, MPI_INTEGER, myrank+1, 0, MPI_COMM_WORLD, ierr)
                endif
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
              
                ! 5
                if(myrank .gt. 0)then
                 call MPI_RECV(temp, 1, MPI_INTEGER, myrank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                endif
              
                open(unit=fh, file=trim(filename), status='old', action='write', &
                form='unformatted', access='stream', position='append')
                write(fh)qw(:,5)
                close(fh)
              
                if(myrank .lt. (numpe-1))then
                 call MPI_SEND(temp, 1, MPI_INTEGER, myrank+1, 0, MPI_COMM_WORLD, ierr)
                endif
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
              
                ! if(myrank .eq. 0) write(*,'(a)')'Exiting write_restart_parallel'
        
        endsubroutine

        subroutine write_restart_parallel()
                implicit none
                include "mpif.h"

                integer :: fileno
                character(len=256) filename
                
                character*100 hostname
                integer*4 status
                integer date_time(8)
                character*10 b(3)
                character*200 date_time_str

                integer i1,i2,dofIdx
                integer :: magicNumber
                character,parameter :: lf=achar(10) !linefeed aka newline

                character(len=:), allocatable :: fh, sh, message
                ! cannot use literal 0, use 0_MPI_OFFSET_KIND for zero
                integer(mpi_offset_kind) displacement, dofoffset

                call hostnm(hostname, status)
                call date_and_time(b(1), b(2), b(3), date_time)
                write(date_time_str, '(a,i0.2,a,i0.2,a,i0.2,x,i0.2,a,i0.2,a,i0.2)') '# ',date_time(1),'/',date_time(2),'/',date_time(3),date_time(5),':',date_time(6),':',date_time(7)
                
                write(filename, '(a, i0, a)') 'restart_par.', currentTimestepIndex, '.0'
                ! variables that don't exist in the flowsolver but won't change?
                magicNumber = 362436
                i1=5
                i2=1

                call MPI_Comm_size(MPI_COMM_WORLD, comsize, ierr)
                call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
                
                
                ! NOTE: Should be using the total number of nshg not just nshg
                ! Also, this file header will be stored on every processor. A bit silly but necessary to be able to set the offsets correctly.
                fh = "# PHASTA Input File Version 2.0"//lf
                fh = fh//"# Byte Order Magic Number : "//itoa(magicNumber)//lf
                fh = fh//"# This result was produced on: "//trim(hostname)//lf
                fh = fh//trim(date_time_str)//lf//lf

                fh = fh//"number of modes : < "//itoa(0)//" > "//itoa(gnno)//lf
                fh = fh//"number of variables : < "//itoa(0)//" > "//itoa(ndof)//lf//lf
                fh = fh//"byteorder magic number  : < "//itoa(i1)//" > "//itoa(i2)//lf
                        
                sh = lf//lf//"solution  : < "//itoa(gnno*ndof*8+1)//" > "//itoa(gnno)//" "//itoa(ndof)//" "//itoa(currentTimeStepIndex)//lf
                
               
                call MPI_File_open(MPI_COMM_WORLD, filename, ior(MPI_MODE_CREATE,MPI_MODE_WRONLY), MPI_INFO_NULL, fileno, ierr)
                
                ! Set our file view to be in bytes
                displacement = 0_MPI_OFFSET_KIND

                call MPI_File_set_view(fileno, displacement, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
                
                if (rank == 0) then
                        ! Write the file header
                        call MPI_File_write_at(fileno, displacement, fh, len(fh), MPI_CHARACTER, wstatus, ierr)
                        ! Write the byte order magic number
                        displacement = displacement + sizeof(fh)
                        call MPI_File_write_at(fileno, displacement, magicNumber, 1, MPI_INT, wstatus, ierr)
                        ! Write the solution header
                        displacement = displacement + sizeof(magicNumber)
                        call MPI_File_write_at(fileno, displacement, sh, len(sh), MPI_CHARACTER, wstatus, ierr)
                endif
                ! For all procs, set the displacement to be the end of the solution header
                displacement = sizeof(fh) + sizeof(magicNumber) + sizeof(sh)
                
                ! Write the dof data
                ! 8 for size of solution data (double precision)
                do dofIdx=1,ndof
                        dofoffset = displacement + (dofIdx-1)*gnno*8 + (ise(1,rank+1)-1)*8
                        call MPI_File_set_view(fileno, dofoffset, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
                        call MPI_File_write_all(fileno, qw(:,dofIdx), isize(rank+1), MPI_DOUBLE_PRECISION, wstatus, ierr)
                enddo

                call MPI_File_close(fileno, ierr)
        endsubroutine write_restart_parallel

      !   interface
      !    integer( kind = C_INT ) function gethostname( hname, len ) bind( C, name = 'gethostname' )
      !          use iso_c_binding
      !          implicit none
      !          character( kind = C_CHAR ) :: hname( * )
      !          integer( kind = C_INT ), VALUE :: len
      !    end function gethostname
      !    end interface
     end