     subroutine e3b (ul,      yl,      acl, &
                     iBCB,    BCB, &
                     shpb,    shglb, &
                     xlb,     xdistl,  xdnvl, &
                     rl,      sgn,     dwl,     xKebe, &
                     SWB)
!
!----------------------------------------------------------------------
!
!   This routine calculates the 3D RHS residual of the fluid boundary 
!   elements.
!
! input:
!  yl     (npro,nshl,ndof)      : Y variables
!  ylotwn (npro,nshl,ndof)      : Y variables at off-the-wall nodes
!  iBCB   (npro,ndiBCB)         : boundary condition code (iBCB(:,1) is
!      a bit tested boundary integral flag i.e.
!                  if set to value of BCB      if set to floating value
!      iBCB(:,1) : convective flux * 1            0  (ditto to all below)
!                  pressure   flux * 2
!                  viscous    flux * 4
!                  heat       flux * 8
!                  turbulence wall * 16
!                  scalarI   flux  * 16*2^I 
!                  (where I is the scalar number)
!
!      iBCB(:,2) is the srfID given by the user in MGI that we will
!                collect integrated fluxes for.
!
!  BCB    (npro,nshlb,ndBCB)    : Boundary Condition values
!                                  BCB (1) : mass flux 
!                                  BCB (2) : pressure 
!                                  BCB (3) : viscous flux in x1-direc.
!                                  BCB (4) : viscous flux in x2-direc.
!                                  BCB (5) : viscous flux in x3-direc.
!                                  BCB (6) : heat flux
!  shpb   (nen,ngaussb)           : boundary element shape-functions
!  shglb  (nsd,nen,ngaussb)       : boundary element grad-shape-functions
!  xlb    (npro,nenl,nsd)       : nodal coordinates at current step
!
! output:
!  rl     (npro,nshl,nflow)      : element residual
!
! Note: Always the first side of the element is on the boundary.  
!       However, note that for higher-order elements the nodes on 
!       the boundary side are not the first nshlb nodes, see the 
!       array mnodeb.
!
!
! Zdenek Johan, Summer 1990.  (Modified from e2b.f)
! Zdenek Johan, Winter 1991.  (Fortran 90)
! Alberto Figueroa, Winter 2004.  CMM-FSI
! Irene Vignon, Spring 2004
!----------------------------------------------------------------------
!
        use incpBC
        use LagrangeMultipliers 
        use multidomain, only: hrt, multidom
        use cpp_interface
!
        use phcommonvars  
        use nnw, only : get_shear_rate 
        
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension yl(npro,nshl,ndof),          iBCB(npro,ndiBCB), &
                  BCB(npro,nshlb,ndBCB),       shpb(nshl,ngaussb), &
                  shglb(nsd,nshl,ngaussb),     xlb(npro,nenl,nsd), &
                  ul(npro,nshl,nsd),           acl(npro,nshl,ndof), &
                  rl(npro,nshl,nflow), &
                  xdistl(npro,nshl),           xdnvl(npro,nshl,nsd), &
                  SWB(npro,nProps)
!
        dimension g1yi(npro,ndof),             g2yi(npro,ndof), &
                  g3yi(npro,ndof),             WdetJb(npro), &
                  bnorm(npro,nsd)
!
        dimension u1(npro),                    u2(npro), &
                  u3(npro),                    rho(npro), &
                  unm(npro),                   pres(npro), &
                  vdot(npro,nsd),              rlKwall(npro,nshlb,nsd), &
                  usup1(npro,nsd),             usup2(npro,nsd), &
                  velsup(npro,nsd)
!
        dimension rou(npro),                   rmu(npro)
!
        dimension tau1n(npro), &
                  tau2n(npro),                 tau3n(npro)
!
        dimension lnode(27),                   sgn(npro,nshl), &
                  shapeVar(npro,nshl),         shdrv(npro,nsd,nshl), &
                  rNa(npro,4)

        integer   ienb(npro,nshl) ! added by KD Lau 04/01/13 
        real*8    gpres(npro)     ! added by KD Lau 04/01/13
        real*8    pval            ! added by KD Lau 22/01/13
        real*8    unm_copy(npro)  ! added by KD Lau 14/06/13   
        real*8    rou_copy(npro)  ! added by KD Lau 14/06/13   
        real*8    stb_pres(npro)  ! added by KD Lau 13/08/13   
        real*8    s_pres          ! added by KD Lau 13/08/13 

        real*8    xmudmi(npro,ngauss),         dwl(npro,nshl)
        integer flowIsPermitted
!
      	dimension xKebe(npro,9,nshl,nshl)

        character(len=50):: file_name, file_name_aux  ! added by MA 02/06/2016
        integer :: result                             ! added by MA 02/06/2016
        real*8 :: gamma_shear(npro)
!
!.... compute the nodes which lie on the boundary (hierarchic)
!
        call getbnodes(lnode)
!
!       ! **************************** !
!       ! *** zero gpres/stb_press *** !
!       ! **************************** !
!
        gpres(:) = real(0.0,8) 
        stb_pres(:) = real(0.0,8)               
!
!.... loop through the integration points
!
        if(lcsyst.eq.3.or.lcsyst.eq.4) then
           ngaussb = nintb(lcsyst)
        else
           ngaussb = nintb(lcsyst)
        endif
        
        do intp = 1, ngaussb
!
!.... get the hierarchic shapeVar functions at this int point
!
        call getshp(shpb,        shglb,        sgn, &
                    shapeVar,       shdrv)
!
!     NOTE I DID NOT PASS THE lnode down.  It is not needed
!     since the shapeVar functions are zero on the boundary
!
!     Note that xmudmi is not calculated at these quadrature 
!     points so you give it a zero.  This has implications.
!     the traction calculated by this approach will include 
!     molecular stresses ONLY.  This is why we will use the 
!     consistent flux method to obtain the forces when doing
!     effective viscosity wall modeling.  When doing slip velocity
!     this is not a problem since the traction is given from the
!     log law relation (not the viscosity).
!
        xmudmi=zero
!
        ! call get_shear_rate_boundary (yl, shg, nshl, ndof, nsd, npro, gamma_shear)
!.... get necessary fluid properties (including eddy viscosity)
!
        call getdiff(dwl, yl,     shapeVar,     xmudmi, xlb,   rmu, rho)
!
!.... calculate the integraton variables
!
        call e3bvar (yl,              acl,             ul, &
                     iBCB,             BCB, &
                     shapeVar, &
                     shdrv, &         
                     xlb,             xdistl,          xdnvl, &
                     lnode,           SWB, &
                     WdetJb, &
                     bnorm,           pres, &
                     u1,              u2,              u3, &
                     rmu,             unm, &
                     tau1n,           tau2n,           tau3n, &
                     vdot, &            
                     usup1,           usup2, &
                     velsup, &
                     rlKwall, &
                     xKebe)
        
        ! write(*,*) "**************************************************"
        ! write(*,*) "DEBUG_ALE = ",DEBUG_ALE
        ! write(*,*) "**************************************************"

#if DEBUG_ALE == 1 
        ! write(*,*) "**************************************************"
        ! write(*,*) "DEBUG_ALE inside = ",DEBUG_ALE
        ! write(*,*) "**************************************************"
        ! write(*,*) 'printing shpb inside e3b'
        ! write(file_name_aux,'(i5.5)') intp
        ! file_name = trim("shpb_")
        ! file_name = trim("shpb_")//trim(file_name_aux)
        ! file_name = trim(file_name)//'.dat'
        ! file_name = trim(file_name)
        ! open(795,file=file_name,status='new')
        ! do i = 1, npro
        !   write(795,'(1(i10),4(e20.10))') i,shapeVar(i,1), shapeVar(i,2), shapeVar(i,3),&
        !                           shapeVar(i,4)                      
        ! end do 
        ! close(795)

        ! write(*,*) 'printing shglb inside e3b'
        ! write(file_name_aux,'(i5.5)') intp
        ! file_name = trim("shglb_")
        ! file_name = trim("shglb_")//trim(file_name_aux)
        ! file_name = trim(file_name)//'.dat'
        ! file_name = trim(file_name)
        ! open(793,file=file_name,status='new')
        ! do i = 1, npro
        !   do j = 1, nsd
        !   write(793,'(2(i10),4(e20.10))') i,j,shdrv(i,j,1), shdrv(i,j,2), shdrv(i,j,3),&
        !                           shdrv(i,j,4) ! nenl = 4
        !   end do                          
        ! end do 
        ! close(793)

        ! write(*,*) 'printing shpb inside e3b'
        ! write(file_name_aux,'(i5.5)') intp
        ! file_name = trim("shpb_")
        ! file_name = trim("shpb_")//trim(file_name_aux)
        ! file_name = trim(file_name)//'.dat'
        ! file_name = trim(file_name)
        ! open(795,file=file_name,status='new')
        ! do i = 1, npro
        !   write(795,'(1(i10),4(e20.10))') i,shapeVar(i,1), shapeVar(i,2), shapeVar(i,3),&
        !                           shapeVar(i,4)                      
        ! end do 
        ! close(795)

#endif
!.... -----------------> boundary conditions <-------------------
!
        do iel = 1, npro
!
!  if we have a nonzero value then
!  calculate the fluxes through this surface 
!
           iface = abs(iBCB(iel,2))
           if (iface .ne. 0 .and. ires.ne.2) then
              flxID(1,iface) =  flxID(1,iface) + WdetJb(iel)! measure area too
              flxID(2,iface) =  flxID(2,iface) - WdetJb(iel) * unm(iel)
              flxID(3,iface) = flxID(3,iface) &
                         - ( tau1n(iel) - bnorm(iel,1)*pres(iel)) &
                         * WdetJb(iel) 
              flxID(4,iface) = flxID(4,iface) &
                         - ( tau2n(iel) - bnorm(iel,2)*pres(iel)) &
                         * WdetJb(iel) 
              flxID(5,iface) = flxID(5,iface) &
                         - ( tau3n(iel) - bnorm(iel,3)*pres(iel)) &
                         * WdetJb(iel) 

           endif
!          
!
! **************************************************************************** ! 
! *** start of influx stabilisation                                        *** !
! ***                                                                      *** !
! *** here a copy of the velocity flux is made: unm = \vec{U}\cdot\vec{n}  *** !
! *** the resultant stabilisation pressure is also calculated here         *** ! 
! ***                                                                      *** !
! **************************************************************************** ! 
!          
          unm_copy(iel) = real(0.0,8)
          if (btest(iBCB(iel,1),1)) then                      
            
            ! check if this element/integration point has influx
            if(unm(iel) .lt. real(0.0,8))then              
              ! store influx to add to residual later
              unm_copy(iel) = unm(iel)
              ! calculate stabilisation pressure on this integration point
              s_pres = stabflux_coeff*rho(iel)*unm(iel)*unm(iel)
              ! integrate stabilisation pressure over element for this integration point             
              stb_pres(iel) = s_pres*WdetJb(iel)             
            end if                         

          end if
!          
! **************************************************************************** ! 
! **************************************************************************** ! 
!
!.... mass flux
!

           if (btest(iBCB(iel,1),0)) then
              unm(iel)  = zero
              do n = 1, nshlb
                 nodlcl = lnode(n)
                 unm(iel) = unm(iel) &
                          + shapeVar(iel,nodlcl) * BCB(iel,n,1)
              enddo
           endif
!
!.... pressure
!
! **************************************************************************** ! 
! ***                                                                      *** !
! *** start of zero pressure / multidomain boundary pressure BC            *** !
! ***                                                                      *** !
! *** in this section we zero the elemental pressure on zero pressure and  *** !
! *** multidomain boundary faces, if using the heart model the pressure is *** !
! *** not zeroed out                                                       *** ! 
! ***                                                                      *** !
! *** note that shape is now shapeVar - KDL, NAN & UWM                     *** !  
! ***                                                                      *** !
! **************************************************************************** ! 
!
        ! check if heart model 
        if (iheart .gt. int(0)) then
          
          ! bit test for zero pressure boundary condition          
          if (btest(iBCB(iel,1),1)) then
            
            ! if heart model surface and valve is closed, do nothing
            if (hrt%hassurfid(iBCB(iel,2)) .and. hrt%isavopen() == int(0)) then                          

            ! else zero apply zero pressure boundary condition              
            else 

              ! zero element pressure
              pres(iel) = zero

              ! apply BCB value, as this is zero pressure is this required?!
              do n = 1, nshlb
                 nodlcl = lnode(n)
                 pres(iel) = pres(iel) + shapeVar(iel,nodlcl) * BCB(iel,n,2)
              end do
            
            end if 
          
          end if
  
        ! if not using heart model        
        else 

          ! bit test for zero pressure boundary 
          if (btest(iBCB(iel,1),1)) then
            ! Check for Netlist boundary which is currently in a state which stops flow
            ! across the boundary, due to closed diodes
            call callCPPDiscoverWhetherFlowPermittedAcrossSurface(iBCB(iel,2),flowIsPermitted)
            if (flowIsPermitted .eq. int(0)) then

            ! else zero apply zero pressure boundary condition              
            else
              ! zero element pressure            
              pres(iel) = zero
              
              ! apply BCB value, as this is zero pressure is this required?!
              do n = 1, nshlb
                nodlcl = lnode(n)
                pres(iel) = pres(iel) + shapeVar(iel,nodlcl) * BCB(iel,n,2)
              end do
            end if
          
          end if
        
        end if 
!        
! **************************************************************************** ! 
! **************************************************************************** !   
!
!.... viscous flux
!        
           if (btest(iBCB(iel,1),2)) then
              tau1n(iel) = zero
              tau2n(iel) = zero
              tau3n(iel) = zero
              do n = 1, nshlb
                 nodlcl = lnode(n)
                 tau1n(iel) = tau1n(iel) &
                            + shapeVar(iel,nodlcl)*BCB(iel,n,3)
                 tau2n(iel) = tau2n(iel) &
                            + shapeVar(iel,nodlcl)*BCB(iel,n,4)
                 tau3n(iel) = tau3n(iel) &
                            + shapeVar(iel,nodlcl)*BCB(iel,n,5)
              enddo
           endif
!
!.... turbulence wall (as a way of checking for deformable wall stiffness)
!
           if (btest(iBCB(iel,1),4)) then
              
!              rlKwall(iel,:,:) = rlKwall(iel,:,:) / ngaussb ! divide by number of gauss points 
              pres(iel) = zero                              ! to avoid the gauss point loop
              tau1n(iel) = zero                             ! and make the traction contribution
              tau2n(iel) = zero                             ! zero
              tau3n(iel) = zero                              
           else
              
              !if (intp.eq.ngaussb) then
              !   rlKwall(iel,:,:) = zero                    ! this is not a deformable element
              !end if
              rlKwall(iel,:,:) = zero

              vdot(iel,:) = zero                             
              usup1(iel,:) = zero
              usup2(iel,:) = zero
              velsup(iel,:) = zero
              
              ! we do this because the mass matrix contribution is already in xKebe via e3bvar
              xKebe(iel,:,:,:) = zero  

           endif
!
!..... to calculate inner products for Lagrange Multipliers
!
           if(Lagrange.gt.zero) then
              do k=1, numLagrangeSrfs
                 if (iBCB(iel,2).eq.nsrflistLagrange(k)) then
                    do n = 1, nshlb
                       nodlcl = lnode(n)
                       do m=1, nsd
                          do i=1, nshlb
                             nodlcl2 = lnode(i)
                             do l=1, nsd
                                do j=1, 3
                                   num=(m-1)*nsd+l
                                   loclhsLag(iel,num,n,i,j)= &
                                      loclhsLag(iel,num,n,i,j)+ &
                                      shapeVar(iel,nodlcl)*WdetJb(iel)* &
                                      shapeVar(iel,nodlcl2)* &
                                      LagInplaneVectors(m,j,k)* &
                                      LagInplaneVectors(l,j,k)
                                enddo 
                             enddo
                          enddo
                       enddo
                    enddo
                 endif
              enddo
           endif  ! end of if(Lagrange.gt.zero)
!
        enddo                                               ! end of bc loop
!
! **************************************
! *** start of multidomain container ***
! **************************************
!
        ! sum up stabilistion pressure from in multidomain container
        call multidom%addstb_pres(npro,stb_pres,iBCB(:,2))        
!
! **************************************
! **************************************  
!

!$$$c.... if we are computing the bdry for the consistent
!$$$c     boundary forces, we must not include the surface elements
!$$$c     in the computataion (could be done MUCH more efficiently!)--->
                                                                  !this
                                                                  !comment should read as for the consistent flux calculation rather than boundary forces
!$$$c
        if (ires .eq. 2) then
           do iel = 1, npro 
              if (nsrflist(iBCB(iel,2)) .ne. 0) then
                 unm(iel) = zero
                 tau1n(iel) = zero
                 tau2n(iel) = zero
                 tau3n(iel) = zero
!                 pres(iel) = zero
!
! whatever is zeroed here will beome part of the post-processed surface
!                 "traction force"
!
! uncomment the next lines to get all of the t vector coming from
!                 Alberto's wall motion model.
                 vdot(iel,:)=zero
                 usup1(iel,:)=zero
                 usup2(iel,:)=zero
                 velsup(iel,:)=zero
                 rlKwall(iel,:,:)=zero

!
! uncomment the next 8 lines to get only the tangential part
!

!                  vn=dot_product(vdot(iel,:),bnorm(iel,:))
!                  vdot(iel,:)=vn*bnorm(iel,:)
!                  walln1=dot_product(rlkwall(iel,1,:),bnorm(iel,:))
!                  walln2=dot_product(rlkwall(iel,2,:),bnorm(iel,:))
!                  walln3=dot_product(rlkwall(iel,3,:),bnorm(iel,:))
!                  rlkwall(iel,1,:)=walln1*bnorm(iel,:)
!                  rlkwall(iel,2,:)=walln2*bnorm(iel,:)
!                  rlkwall(iel,3,:)=walln3*bnorm(iel,:)
              endif
           enddo
        endif
!
!.... assemble the contributions
!
        rNa(:,1) = -tau1n + bnorm(:,1)*pres 
        rNa(:,2) = -tau2n + bnorm(:,2)*pres 
        rNa(:,3) = -tau3n + bnorm(:,3)*pres 
        rNa(:,4) = unm
        
        if (ideformwall.eq.1) then

           rNa(:,1) = rNa(:,1) + vdot(:,1)
           rNa(:,2) = rNa(:,2) + vdot(:,2)
           rNa(:,3) = rNa(:,3) + vdot(:,3)

           if (iwalldamp.eq.1) then  
              rNa(:,1) = rNa(:,1) + velsup(:,1)
              rNa(:,2) = rNa(:,2) + velsup(:,2)
              rNa(:,3) = rNa(:,3) + velsup(:,3)
           endif
        
           if (iwallsupp.eq.1) then  
              rNa(:,1) = rNa(:,1) + usup1(:,1)
              rNa(:,2) = rNa(:,2) + usup1(:,2)
              rNa(:,3) = rNa(:,3) + usup1(:,3)
           endif
        
           if (idistancenudge.eq.1) then
              rNa(:,1) = rNa(:,1) + usup2(:,1)
              rNa(:,2) = rNa(:,2) + usup2(:,2)
              rNa(:,3) = rNa(:,3) + usup2(:,3)
           endif
           
        endif
        
!
        if(iconvflow.eq.1) then     ! conservative form was integrated
                                    ! by parts and has a convective 
                                    ! boundary integral
!
!.... assemble the contributions
!
           rou=rho*unm
           rNa(:,1) = rNa(:,1) + rou*u1 
           rNa(:,2) = rNa(:,2) + rou*u2
           rNa(:,3) = rNa(:,3) + rou*u3
        endif
!       
! *************************************
! *** start of influx stabilisation ***
! *************************************
! 
        ! calculate rho (density) x unm
        rou_copy=rho*unm_copy
        ! minus stabilised term 
        rNa(:,1) = rNa(:,1) - stabflux_coeff*rou_copy*u1                
        rNa(:,2) = rNa(:,2) - stabflux_coeff*rou_copy*u2
        rNa(:,3) = rNa(:,3) - stabflux_coeff*rou_copy*u3        
!        
! *************************************
! *************************************
!
        rNa(:,1) = rNa(:,1)*WdetJb
        rNa(:,2) = rNa(:,2)*WdetJb
        rNa(:,3) = rNa(:,3)*WdetJb
        rNa(:,4) = rNa(:,4)*WdetJb
        
!
!.... ------------------------->  Residual  <--------------------------
!
!.... add the flux to the residual
!
        do n = 1, nshlb
           nodlcl = lnode(n)

           rl(:,nodlcl,1) = rl(:,nodlcl,1) - shapeVar(:,nodlcl) * rNa(:,1)
           rl(:,nodlcl,2) = rl(:,nodlcl,2) - shapeVar(:,nodlcl) * rNa(:,2)
           rl(:,nodlcl,3) = rl(:,nodlcl,3) - shapeVar(:,nodlcl) * rNa(:,3)
           rl(:,nodlcl,4) = rl(:,nodlcl,4) - shapeVar(:,nodlcl) * rNa(:,4)

        enddo
        
        !if(ideformwall.eq.1 .and. intp.eq.ngaussb) then
        if(ideformwall.eq.1) then
           rl(:,1,1) = rl(:,1,1) - rlKwall(:,1,1)
           rl(:,1,2) = rl(:,1,2) - rlKwall(:,1,2)
           rl(:,1,3) = rl(:,1,3) - rlKwall(:,1,3)

           rl(:,2,1) = rl(:,2,1) - rlKwall(:,2,1)
           rl(:,2,2) = rl(:,2,2) - rlKwall(:,2,2)
           rl(:,2,3) = rl(:,2,3) - rlKwall(:,2,3)

           rl(:,3,1) = rl(:,3,1) - rlKwall(:,3,1)
           rl(:,3,2) = rl(:,3,2) - rlKwall(:,3,2)
           rl(:,3,3) = rl(:,3,3) - rlKwall(:,3,3)
        endif
!
!.... -------------------->  Aerodynamic Forces  <---------------------
!
          if (((ires .ne. 2) .and. (iter .eq. nitr))) then ! \todo-binary i removed the itwmod line as it is never initialised
        ! if (((ires .ne. 2) .and. (iter .eq. nitr)) &
        !                    .and. (abs(itwmod).eq.1)) then
!
!.... compute the forces on the body
!
           where (nsrflist(iBCB(:,2)).eq.1)
              tau1n = ( tau1n - bnorm(:,1)*pres)  * WdetJb
              tau2n = ( tau2n - bnorm(:,2)*pres)  * WdetJb
              tau3n = ( tau3n - bnorm(:,3)*pres)  * WdetJb
           elsewhere
              tau1n = zero
              tau2n = zero
              tau3n = zero
           endwhere
!
! Note that the sign has changed from the compressible code to 
! make it consistent with the way the bflux calculates the forces
! Note also that Hflux has moved to e3btemp
!
          Force(1) = Force(1) - sum(tau1n)
          Force(2) = Force(2) - sum(tau2n)
          Force(3) = Force(3) - sum(tau3n)

!
        endif
!
!.... end of integration loop
!
        enddo

! #if DEBUG_ALE == 1 
!         stop 
! #endif
        
!$$$        ttim(40) = ttim(40) + tmr()
!
!.... return
!
        return
        end


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*********************************************************************
!*********************************************************************


        subroutine e3bSclr (yl,      iBCB,    BCB,     shpb,    shglb, &
                            xlb,     rl,      sgn,     dwl)
        use phcommonvars  
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension yl(npro,nshl,ndof),          iBCB(npro,ndiBCB), &
                  BCB(npro,nshlb,ndBCB),       shpb(nshl,*), &
                  shglb(nsd,nshl,*), &
                  xlb(npro,nenl,nsd), &         
                  rl(npro,nshl)
!
        real*8    WdetJb(npro),                bnorm(npro,nsd)
!
        dimension lnode(27),                   sgn(npro,nshl), &
                  shapeVar(npro,nshl),            shdrv(npro,nsd,nshl), &
                  rNa(npro),                   flux(npro)
        real*8    dwl(npro,nshl)

!
!.... compute the nodes which lie on the boundary (hierarchic)
!
        call getbnodes(lnode)
!
!.... loop through the integration points
!
        if(lcsyst.eq.3.or.lcsyst.eq.4) then
           ngaussb = nintb(lcsyst)
        else
           ngaussb = nintb(lcsyst)
        endif
        do intp = 1, ngaussb
!
!.... get the hierarchic shapeVar functions at this int point
!
        call getshp(shpb,        shglb,        sgn, &
                    shapeVar,       shdrv)
!
!.... calculate the integraton variables
!
        call e3bvarSclr (yl,          shdrv,   xlb, &
                         shapeVar,       WdetJb,  bnorm, &
                         flux,        dwl )
!        
!.... -----------------> boundary conditions <-------------------
!

!
!.... heat or scalar  flux
!     
        if(isclr.eq.0) then 
           iwalljump=0
        else
           iwalljump=1  !turb wall between heat and scalar flux..jump over
        endif
        ib=4+isclr+iwalljump
        ibb=6+isclr
        do iel=1, npro
!
!  if we have a nonzero value then
!  calculate the fluxes through this surface 
!
           if (iBCB(iel,2) .ne. 0 .and. ires.ne.2) then
              iface = abs(iBCB(iel,2))
              flxID(ibb,iface) =  flxID(ibb,iface) &
                                - WdetJb(iel) * flux(iel)
           endif

           if (btest(iBCB(iel,1),ib-1)) then
              flux(iel) = zero
              do n = 1, nshlb
                 nodlcl = lnode(n)
                 flux(iel) = flux(iel) &
                           + shapeVar(iel,nodlcl) * BCB(iel,n,ibb)
              enddo           
           endif
        enddo
!
!.... assemble the contributions
!
        rNa(:) = -WdetJb * flux
!
!.... ------------------------->  Residual  <--------------------------
!
!.... add the flux to the residual
!
        do n = 1, nshlb
           nodlcl = lnode(n)
 
           rl(:,nodlcl) = rl(:,nodlcl) - shapeVar(:,nodlcl) * rNa(:)
        enddo
!
!.... -------------------->  Aerodynamic Forces  <---------------------
!
        if ((ires .ne. 2) .and. (iter .eq. nitr)) then
!
!.... compute the forces on the body
!
           if(isclr.eq.0)   HFlux    = sum(flux)
!
        endif
!
!.... end of integration loop
!
        enddo

!
!.... return
!
        return
        end
  
