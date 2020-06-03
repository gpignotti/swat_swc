      subroutine percmicro(ly1)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine computes percolation and lateral subsurface flow
!!    from a soil layer when field capacity is exceeded

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hru_slp(:)   |m/m           |average slope steepness
!!    ihru		 |none          |HRU number
!!    iwatable     |none          |high water table code:
!!                                |0 no high water table
!!                                |1 high water table
!!    ldrain(:)    |none          |soil layer where drainage tile is located
!!    slsoil(:)    |m             |slope length for lateral subsurface flow
!!    sol_fc(:,:)  |mm H2O        |amount of water available to plants in soil
!!                                |layer at field capacity (fc - wp water)
!!    sol_hk(:,:)  |none          |beta coefficient to calculate hydraulic 
!!                                |conductivity
!!    sol_k(:,:)   |mm/hr         |saturated hydraulic conductivity of soil
!!                                |layer
!!    sol_nly(:)   |none          |number of soil layers in HRU
!!    sol_st(:,:)  |mm H2O        |amount of water stored in the soil layer
!!                                |on any given day
!!    sol_sumfc(:) |mm H2O        |amount of water held in the soil profile
!!                                |at field capacity
!!    sol_sw(:)    |mm H2O        |amount of water stored in the soil profile
!!                                |on any given day
!!    sol_tmp(:,:) |deg C         |daily average temperature of soil layer
!!    sol_ul(:,:)  |mm H2O        |amount of water held in the soil layer at
!!                                |saturation (sat - wp water)
!!    sol_z(:,:)   |mm            |depth to bottom of soil layer
!!    sw_excess    |mm H2O        |amount of water in soil that exceeds field 
!!                                |capacity (gravity drained water)
!!    tdrain(:)    |hrs           |time to drain soil to field capacity
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    latlyr       |mm H2O        |lateral subsurface flow in layer
!!    lyrtile      |mm H2O        |drainage tile flow in layer for day in HRU
!!    sepday       |mm H2O        |percolation from soil layer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    adjf         |none          |adjustment factor for lateral flow
!!    dg           |mm            |depth of soil layer
!!    ho           |none          |variable to hold intermediate calculation
!!                                |result
!!    j            |none          |HRU number
!!    ly1          |none          |soil layer number
!!    ratio        |none          |ratio of seepage to (latq + sepday)
!!    yy           |mm            |depth to top of soil layer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer, intent (in) :: ly1
      integer :: j
      real :: adjf, yy, dg, ho, ratio, sol_k_sep
      real :: m, minv, n

      j = 0
      j = ihru

      adjf = 1.

      !! if temperature of layer is 0 degrees C or below
      !! there is no water flow
      if (sol_tmp(ly1,j) <= 0.) then
        sepday = 0.
        return
      end if

        !! COMPUTE LATERAL FLOW USING HILLSLOPE STORAGE METHOD
        if (ly1 == 1) then
          yy = 0.
        else
          yy = 0.
          yy = sol_z(ly1-1,j)
        end if

        dg = 0.
        ho = 0.
        latlyr = 0.
        dg = sol_z(ly1,j) - yy
        if (sol_ul(ly1,j) - sol_fc(ly1,j)==0.) then
          ho=0.
        else
          ho = 2. * sw_excess / ((sol_ul(ly1,j) - sol_fc(ly1,j)) /  dg)
        end if
        !the factor must be changed from 0.024 to 0.001 because we changing from daily to hourly calculation
!        latlyr = adjf * ho * sol_k(ly1,j) * hru_slp(j) / slsoil(j)      
!     &                                                            * .024
        latlyr = adjf * ho * sol_k(ly1,j) * hru_slp(j) / slsoil(j)      
     &                                                            * .001

      if (latlyr < 0.) latlyr = 0. 
      if (latlyr > sw_excess) latlyr = sw_excess

      sol_hk(ly1,j) = (sol_ul(ly1,j) - sol_fc(ly1,j)) / sol_k(ly1,j)

!!  septic changes 1/28/09 
      if (ly1 == i_sep(j)) then
         if (isep_opt(j) == 1) then !active system
           sol_k_sep = sol_k(ly1,j)* 
     &              (sol_st(ly1,j) - sol_fc(ly1,j))/
     &              (sol_ul(ly1,j) - sol_fc(ly1,j))
           sol_k_sep = Max(1.e-6, sol_k_sep)
           sol_k_sep = Min(sol_k(ly1,j), sol_k_sep)
           
           sol_hk(ly1,j) = (sol_ul(ly1,j) - sol_fc(ly1,j)) 
     &      / sol_k_sep
         
         elseif (isep_opt(j) == 2) then !failing system
           sol_hk(ly1,j) = 1.e10
         endif
      endif 
!!  septic changes 1/28/09

!!SWC edits by GWP
      if (iperc == 0) then ! original calculation
          sol_hk(ly1,j) = Max(2., sol_hk(ly1,j))

      !!  compute seepage to the next layer
          sepday = 0.
          sepday = sw_excess * (1. - Exp(-1. / sol_hk(ly1,j)))
          
      else if (iperc == 1) then ! Campbell-Rawls
          !if (ly1 == 1) then
			!z = sol_z(ly1,j)
          !else
			!z = sol_z(ly1,j) - sol_z(ly1-1,j)
		!endif
          if (ly1 == 1) z = sol_z(ly1,j)
          if (ly1 > 1) z = sol_z(ly1,j) - sol_z(ly1-1,j)
          sepday = 0.
          if (sol_st(ly1,j) > 0.99 * sol_ul(ly1,j)) then
              sepday = sol_k(ly1,j) * 1.
          else
              theta = sol_st(ly1,j) / z + sol_wp(ly1,j)
              se = theta  /
     &            sol_por(ly1,j)
!              if (se > 1.) then
!                  se = 1.
!              end if
              n = 3. + 2. * ca_rb_b(ly1,j)
              sol_kun = sol_k(ly1,j) * se ** n
              sepday = sol_kun * 1.
         end if
          !sepday = min(sw_excess, sepday)
          !sepday = min(sepday, sol_st(ly1,j))

!      else if (iperc == 1) then ! Campbell-Rawls
!     
!          if (ly1 == 1) z = sol_z(ly1,j)
!          if (ly1 > 1) z = sol_z(ly1,j) - sol_z(ly1-1,j)
!          sepday = 0.
!         !if (sol_st(ly1,j) > 0.99 * sol_ul(ly1,j)) then
!          !    sepday = sol_k(ly1,j) * 1.
!          !else
!          theta = sol_st(ly1,j) / z + sol_wp(ly1,j)
!          se = theta  /
!     &        sol_por(ly1,j)
!          n = 3. + 2. * ca_rb_b(ly1,j)
!         sol_kun = sol_k(ly1,j) * se ** n
!          !sepday = sol_kun * 1.
!             
!          sol_hk(ly1,j) = (sol_ul(ly1,j) - sol_fc(ly1,j)) / 
!     &        sol_kun
!          !sepday = sw_excess * (1. - Exp(-1. / sol_hk(ly1,j)))
!          sepday = sol_st(ly1,j) * (1. - Exp(-1. / sol_hk(ly1,j)))
 !             
  !       end if


      
      else if (iperc == 2) then ! VanGenuchten-Rawls
          if (ly1 == 1) z = sol_z(ly1,j)
          if (ly1 > 1) z = sol_z(ly1,j) - sol_z(ly1-1,j)
          sepday = 0.
          if (sol_st(ly1,j) > 0.99 * sol_ul(ly1,j)) then
              sepday = sol_k(ly1,j) * 1.
          else
              theta = sol_st(ly1,j) / z + sol_wp(ly1,j)
              se = (theta - rb_thr(ly1,j)) /
     &            (sol_por(ly1,j) - rb_thr(ly1,j))
              m = vg_rb_m(ly1,j)
              if (m > 0) then
                  minv = 1/m
              else
                  m = 0
              end if
              factor_1 = 1. - se ** minv
              factor_2 = factor_1 ** m
              factor_3 = 1. - factor_2
              factor_4 = factor_3 ** 2.
              if (se > 0) then
                  sol_kun = sol_k(ly1,j) * factor_4 * se ** 0.5
              else
                  sol_kun = 0.
              endif
              sepday = sol_kun * 1.
         end if
           !sepday = min(sepday, sol_st(ly1,j))
      
      else if (iperc == 3) then ! Campbell-Cosby
          if (ly1 == 1) z = sol_z(ly1,j)
          if (ly1 > 1) z = sol_z(ly1,j) - sol_z(ly1-1,j)
          sepday = 0.
          if (sol_st(ly1,j) > 0.99 * sol_ul(ly1,j)) then
              sepday = sol_k(ly1,j) * 1.
          else
              theta = sol_st(ly1,j) / z + sol_wp(ly1,j)
              se = theta  /
     &            sol_por(ly1,j)
              n = 3. + 2. * ca_co_b(ly1,j)
              sol_kun = sol_k(ly1,j) * se ** n
              sepday = sol_kun * 1.
         end if
           !sepday = min(sepday, sol_st(ly1,j))         
  
      
      else if (iperc == 4) then ! VanGenuchten-Cosby
          if (ly1 == 1) z = sol_z(ly1,j)
          if (ly1 > 1) z = sol_z(ly1,j) - sol_z(ly1-1,j)
          sepday = 0.
          if (sol_st(ly1,j) > 0.99 * sol_ul(ly1,j)) then
              sepday = sol_k(ly1,j) * 1.
          else
              theta = sol_st(ly1,j) / z + sol_wp(ly1,j)
              se = theta  /
     &            sol_por(ly1,j)
              m = vg_co_m(ly1,j)
              if (m > 0) then
                  minv = 1/m
              else
                  m = 0
              end if
              factor_1 = 1. - se ** minv
              factor_2 = factor_1 ** m
              factor_3 = 1. - factor_2
              factor_4 = factor_3 ** 2.
              sol_kun = sol_k(ly1,j) * factor_4 * se ** 0.5
              sepday = sol_kun * 1.
         end if
           !sepday = min(sepday, sol_st(ly1,j))
          
      
      else if (iperc == 5) then ! Campbell-Saxton
          if (ly1 == 1) z = sol_z(ly1,j)
          if (ly1 > 1) z = sol_z(ly1,j) - sol_z(ly1-1,j)
          sepday = 0.
          if (sol_st(ly1,j) > 0.99 * sol_ul(ly1,j)) then
              sepday = sol_k(ly1,j) * 1.
          else
              theta = sol_st(ly1,j) / z + sol_wp(ly1,j)
              se = theta  /
     &            sol_por(ly1,j)
              n = 3. + 2. * ca_sx_b(ly1,j)
              sol_kun = sol_k(ly1,j) * se ** n
              sepday = sol_kun * 1.
         end if
           !sepday = min(sepday, sol_st(ly1,j))
          
      
      else if (iperc == 6) then ! VanGenuchten-Saxton
          if (ly1 == 1) z = sol_z(ly1,j)
          if (ly1 > 1) z = sol_z(ly1,j) - sol_z(ly1-1,j)
          sepday = 0.
          if (sol_st(ly1,j) > 0.99 * sol_ul(ly1,j)) then
              sepday = sol_k(ly1,j) * 1.
          else
              theta = sol_st(ly1,j) / z + sol_wp(ly1,j)
              se = theta /
     &            sol_por(ly1,j)
              m = vg_sx_m(ly1,j)
              if (m > 0) then
                  minv = 1/m
              else
                  m = 0.
              end if
              factor_1 = 1. - se ** minv
              factor_2 = factor_1 ** m
              factor_3 = 1. - factor_2
              factor_4 = factor_3 ** 2.
              sol_kun = sol_k(ly1,j) * factor_4 * se ** 0.5
              sepday = sol_kun * 1.
         end if
           !sepday = min(sepday, sol_st(ly1,j))
          
      
      else if (iperc == 7) then ! Campbell-Wosten
          if (ly1 == 1) z = sol_z(ly1,j)
          if (ly1 > 1) z = sol_z(ly1,j) - sol_z(ly1-1,j)
          sepday = 0.
          if (sol_st(ly1,j) > 0.99 * sol_ul(ly1,j)) then
              sepday = sol_k(ly1,j) * 1.
          else
              theta = sol_st(ly1,j) / z + sol_wp(ly1,j)
              se = theta  /
     &            sol_por(ly1,j)
              n = 3. + 2. * ca_wo_b(ly1,j)
              sol_kun = sol_k(ly1,j) * se ** n
              sepday = sol_kun * 1.
         end if
           !sepday = min(sepday, sol_st(ly1,j))
          
      
      else if (iperc == 8) then ! VanGenuchten-Wosten
          if (ly1 == 1) z = sol_z(ly1,j)
          if (ly1 > 1) z = sol_z(ly1,j) - sol_z(ly1-1,j)
          sepday = 0.
          if (sol_st(ly1,j) > 0.99 * sol_ul(ly1,j)) then
              sepday = sol_k(ly1,j) * 1.
          else
              theta = sol_st(ly1,j) / z + sol_wp(ly1,j)
              se = (theta - wo_thr(ly1,j)) /
     &            (sol_por(ly1,j) - wo_thr(ly1,j))
              m = vg_wo_m(ly1,j)
              if (m > 0) then
                  minv = 1/m
              else
                  m = 0.
              end if
              factor_1 = 1. - se ** minv
              factor_2 = factor_1 ** m
              factor_3 = 1. - factor_2
              factor_4 = factor_3 ** 2.
              if (se > 0) then
                  sol_kun = sol_k(ly1,j) * factor_4 * se ** 0.5
              else
                  sol_kun = 0.
              endif
              sepday = sol_kun * 1.
         end if
           !sepday = min(sepday, sol_st(ly1,j))
               
      end if
!!End GWP edits
      
      !! limit maximum seepage from biozone layer below potential perc amount
	if(ly1 == i_sep(j).and.isep_opt(j)==1) then
	   !sepday = min(sepday,sol_k_sep *24.)
          sepday = min(sepday,sol_k_sep *1.)
	   bz_perc(j) = sepday
	end if
      
      !! restrict seepage if next layer is saturated
      if (ly1 == sol_nly(j)) then
        xx = (dep_imp(j) - sol_z(ly1,j)) / 1000.
        if (xx < 1.e-4) then
          sepday = 0.
        else
          sepday = sepday * xx / (xx + Exp(8.833 - 2.598 * xx))
        end if
      end if


      ! this is now done in percmain    
            !! check mass balance moved to percmain
!      if (sepday + latlyr > sw_excess) then
!        ratio = 0.
!        ratio = sepday / (latlyr + sepday)
!        sepday = 0.
!        latlyr = 0.
!        sepday = sw_excess * ratio
!        latlyr = sw_excess * (1. - ratio)
!      endif
!      if (sepday + lyrtile > sw_excess) then
!        sepday = 0.
!        sepday = sw_excess - lyrtile
!      endif

      
      return
      end