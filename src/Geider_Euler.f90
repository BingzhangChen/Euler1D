SUBROUTINE BIOLOGY
USE params
USE state_variables
USE forcing,          only : Temp, PAR
USE Trait_functions,  only : TEMPBOL, PHY_C2Vol, palatability
USE grid,             only : Hz, nlev, Z_r
USE Time_setting,     only : dtdays, sec_of_day
implicit none
INTEGER :: k, i, j,j1, j2, m,kk
INTEGER :: i_alpha, i_Topt, i_ESD
real    :: NO3   = 0.
real    :: DET   = 0. 
real    :: tf_z  = 0.
real    :: tf_p  = 0.
real    :: Graz  = 0.
real    :: dC_   = 0.
real    :: dN_   = 0.
real    :: dChl_ = 0.
real    :: uptake= 0.   !Total NO3 uptake
real    :: NPPc_(nlev)  = 0.  !C-based phytoplankton production (mg C m-3 d-1)
real    ::  IPAR(nlev)  = 0.  !Daily integrated PAR at each depth
real    :: pp_DZ = 0.   
real    :: pp_ND = 0.   
real    :: Pmort = 0. !Phytoplankton natural mortality 
real    ::  ZOO(NZOO) = 0. 
real    :: PHYN(NPHY) = 0. 
real    :: PHYN1(NPHY) = 0. 
real    :: PHYC(NPHY) = 0. 
real    :: PHYC1(NPHY) = 0. 
real    ::  CHL(NPHY) = 0. 
real    ::  CHL1(NPHY) = 0. 
real    ::   dC(NPHY) = 0.   !Scratch variable indicating changes in phytoplankton carbon
real    ::   dN(NPHY) = 0.   !Scratch variable indicating changes in phytoplankton nitrogen
real    :: dCHL(NPHY) = 0.   !Scratch variable indicating changes in Chl
real    ::  muC(NPHY) = 0.   !Scratch variable indicating specific growth rate of phytoplankton carbon
real    ::  muN(NPHY) = 0.   !Scratch variable indicating specific growth rate of phytoplankton nitrogen
real    :: muCHL(NPHY)= 0.   !Scratch variable indicating specific growth rate of CHL
real    :: FZoo(NZOO) = 0.   !The total amount of palatable prey (in Nitrogen)
                             !for each zooplankton size class
real    :: dNO3dt   = 0.   
real    :: dDETdt   = 0.   
real    :: dZOOdt(NZOO)   = 0.   
real    :: RES   = 0.   
real    :: EGES  = 0.   
real    :: gbar  = 0.   
real    :: INGES(NZOO) = 0.   
real    :: Zmort = 0.   
real    :: Spec_graz = 0.   !Specific grazing rate
real    :: Gmatrix(NZOO,NZOO) = 0.d0     !Grazer biomass specific grazing rate matrix
real    :: Pmatrix(NPHY,NZOO) = 0.d0    !Phytoplankton mortality rates by each zooplankton size class for each phyto. species
real, parameter   :: eta     = -1.d0*6.6  !Prey refuge parameter for nitrogen
real, parameter   :: A_g    = 21.9   !Intercept of the allometric equation of maximal zooplankton grazing rate (Ward et al. 2012)
real, parameter   :: B_g    = -0.16  !Slope of the allometric equation of maximal zooplankton grazing rate (Ward et al. 2012)
real, parameter   :: Pm     = 0.01  !Phytoplankton natural specific mortality rate (d-1)
!End of declaration

DO k = nlev, 1, -1
   NO3 = t(iNO3, k)
   DET = t(iDET, k)
   Varout(oTEMP,k) = Temp(k)
   IPAR(k) = IPAR(k) + PAR(k)*dtdays   !Unit: W m-2

   do kk = 1, NZOO
      ZOO(kk)= t(iZOO(kk), k)
   enddo

   do kk = 1, NPHY
      PHYC(kk)= t(iPHYC(kk), k)
      PHYN(kk)= t(iPHYN(kk), k)
      CHL(kk) = t(iCHL(kk),  k)
   enddo

   if (sec_of_day == 0) then
       Varout(oNPP, k) = NPPc_(k)  !NPP of the past day; this is real NPP (mg C d-1 m-3)
       Varout(oPAR, k) = IPAR(k)   !Integrated PAR of the past day (W m-2)
       NPPc_(k)   = 0d0              !Reset NPPc_
       IPAR(k)    = 0d0              !Reset IPAR
   endif

   !The multiple zooplankton size class model follows Ward et al. L&O 2012
   ! In the NPZD model, phytoplankton cells utilize DIN and are eaten by zooplankton. 
   !The ingested food by zooplankton has three fates: 
   !1) being recycled to DIN; 2) being converted to detritus; and 3) supporting zooplankton growth. 
   !The natural mortality of zooplankton are converted to detritus which is recycled to DIN and also sinks.

   tf_z = TEMPBOL(Ez,Temp(k))
   tf_p = TEMPBOL(Ep,Temp(k))

   !Calculate the total amount of prey N biomass available to each size class of zooplankton
   IF (NZOO > 1) THEN
      DO kk = 1, NZOO
         gmax = A_g * VolZOO(kk)**B_g 

         FZoo(kk) = 0d0

         !First calculate total phyto. prey 
         do m = 1, NPHY

           !The amount of patalable prey in phytoplankton m
           Pmatrix(m,kk) = palatability(VolZOO(kk), VolPHY(m), SDZoo) * PHYN(m)

           !Calculate the palatability of each prey superindividual and add to the total amount palatable prey
           FZoo(kk) = FZoo(kk) + Pmatrix(m,kk)
         enddo

	     !Second, calculate the total zooplankton prey
	     IF (kk > 1) THEN
	       do m = 1, (kk - 1)

          !Save the palatability into Gmatrix
          Gmatrix(m,kk) = palatability(VolZOO(kk), VolZOO(m), SDZoo) 

		    !Calculate the palatability of each zoo. prey and add to the total amount palatable prey
		    FZoo(kk) = FZoo(kk) + Gmatrix(m,kk) * ZOO(m)

	       enddo
	     ENDIF

        !Save the total available prey for zooplankton
        gbar = FZoo(kk)/(FZoo(kk) + Kp)*(1.d0 - exp(eta *FZoo(kk)))

        !Total ingestion of zooplankton kk (mmol N m-3 d-1)
        INGES(kk) = ZOO(kk)*gmax*tf_z*gbar

	     IF (kk > 1) THEN
	       do m = 1, (kk - 1)

            !Calculate the total ingestion rate (mmol N m-3 d-1) of zooplankton kk on zooplankton m
            if (FZOO(kk) > 0d0) then
               Gmatrix(m,kk) = Gmatrix(m,kk)*ZOO(m)/FZoo(kk)*INGES(kk)
            else
               Gmatrix(m,kk) = 0d0
            endif
           enddo
	     ENDIF

      ENDDO !End of the zooplankton loop
   ELSE
      !Calculate the sum of total phytoplankton biomass
      FZOO(1) = sum(PHYN(:))

      gbar = FZoo(1)/(FZoo(1) + Kp)*(1.d0 - exp(eta *FZoo(1)))

      !Total ingestion of zooplankton (mmol N m-3 d-1)
      INGES(1) = ZOO(1)*gmax*tf_z*gbar

   ENDIF

   !Computing zooplankton mortality
   RES  = 0d0  !Total amount of nitrogen that is excreted by zooplankton and becomes DIN
   EGES = 0d0  !Total egestion by zooplankton to detritus

   DO kk = 1, NZOO

      !Zooplankton excretion rate (-> DIN)
      RES  = RES + INGES(kk)*(1d0-GGE-unass)

      !ZOOPLANKTON EGESTION (-> Detritus)
      EGES = EGES + INGES(kk)*unass

      !Calculate zooplankton mortality
      Zmort = ZOO(kk)*mz *tf_z     !Linear Mortality term

      !Loop through all predators
      if (kk .lt. NZOO) then
        do m = (kk + 1), NZOO
          Zmort = Zmort + Gmatrix(kk,m)    !Linear Mortality term + grazing by other ZOO
        enddo
      endif

      !Update the biomass of ZOOplankton kk
      dZOOdt(kk) = GGE*INGES(kk) - Zmort
      t(iZOO(kk),k) = max(ZOO(kk) + dtdays*dZOOdt(kk), 0d0)
      Varout(iZOO(kk), k) = t(iZOO(kk), k)

   ENDDO !End of the zooplankton loop

   ! For production/destruction matrix:
   pp_ND = RDN*DET*tf_z                   !Flux from DET to DIN
   pp_DZ = EGES + mz*tf_z*sum(ZOO(:))     !Flux from ZOO to DET 
  
   !Now calculate dynamics of phytoplankton
   !Impose the zooplankton grazing
   do j = 1, NPHY
      Graz = 0d0  !Total ingested phytoplankton by zooplankton

      !Calculate all the zooplankton ingestion for this phytoplankton species (unit: mmol N m-3 d-1)
      if (NZOO > 1) then
         do m = 1, NZOO
           if (FZOO(m) > 0d0) then
              Graz = Graz + INGES(m) * Pmatrix(j, m)/FZoo(m)
           else
              Graz = Graz 
           endif
         enddo
         !Specific grazing rate (d-1)
         Spec_graz = Graz*dtdays/PHYN(j)
      else
         Spec_graz = INGES(1)*dtdays/FZoo(1) !No discrimination among different phyto. prey
      endif

      PHYN1(j) = PHYN(j)*(1d0 - Spec_graz)   !Apply grazing
      PHYC1(j) = PHYC(j)*(1d0 - Spec_graz)   !Apply grazing
       CHL1(j) =  CHL(j)*(1d0 - Spec_graz)   !Apply grazing

   enddo

   ! Calculate total phytoplankton nitrogen uptake, mortality, and PP
   uptake  = 0d0
   Pmort   = 0d0
   do j = 1, NPHY

      SELECTCASE (Model_ID)
      CASE(GMK98_simple)
         call GMK98_simple_mod(Temp(k), PAR(k), NO3, PHYC(j), PHYN(j), CHL(j), &
                                 dC(j), dN(j), dChl(j) )

      CASE(GMK98_Topt)
         stop "To be developed..."
      CASE(GMK98_Light)
         stop "To be developed..."
      CASE(GMK98_Size)
         call GMK98_Euler_Size(Temp(k), PAR(k), NO3, PHYC(j), PHYN(j), &
                             CHL(j), VolPHY(j), dC(j), dN(j), dCHL(j))
      CASE(GMK98_SizeLight)
         stop "To be developed..."
      CASE(GMK98_ToptLight)
         stop "To be developed..."
      CASE(GMK98_ToptSize)
         stop "To be developed..."
      CASE(GMK98_ToptSizeLight)
         stop "To be developed..."
      CASE DEFAULT
         stop "Model choice is wrong!!"
      END SELECT

      uptake   =   uptake + dN(j)
      NPPc_(k) = NPPc_(k) + dC(j)*12.d0*dtdays !Unit: mgC m-3 d-1
      Pmort    =    Pmort + Pm * PHYN(j)*tf_P

      if (NPHY > 1) then
         !Calculate growth rate
         muC(j)   = dC(j)/PHYC(j)
         muN(j)   = dN(j)/PHYN(j)
         muCHL(j) = dCHL(j)/CHL(j)
      endif

   enddo

   !Trait mutation
   If (NPHY > 1) Then
      do j = 1, NPHY

        !Find the index corresponding to each trait axis
        i_alpha = mod(j, N_alpha)
        if (i_alpha .eq. 0) i_alpha = N_alpha
        i_Topt  = mod(j-i_alpha, N_Topt) + 1
        i_ESD   = (j-i_alpha-(i_Topt-1)*N_alpha)/(N_alpha*N_Topt) + 1

        if(i_ESD == 1) then
         !Find the next index along ESD
         j2 = j + N_alpha * N_Topt

         !Update cellular C, N, and Chl
         PHYC1(j) = PHYC1(j) + (dC(j) - delta*muC(j)*PHYC(j) + delta*muC(j2)*PHYC(j2)-PHYC(j)*Pm*tf_p)*dtdays
         PHYN1(j) = PHYN1(j) + (dN(j) - delta*muN(j)*PHYN(j) + delta*muN(j2)*PHYN(j2)-PHYN(j)*Pm*tf_p)*dtdays
          CHL1(j) =  CHL1(j) + (dCHL(j)-delta*muCHL(j)*CHL(j)+ delta*muCHL(j2)*CHL(j2)-CHL(j)*Pm*tf_p)*dtdays
        elseif (i_ESD == N_ESD) then
         !Find the previous index along ESD
         j1 = j - N_alpha * N_Topt
         PHYC1(j) = PHYC1(j) + (dC(j) - delta*muC(j)*PHYC(j) + delta*muC(j1)*PHYC(j1)-PHYC(j)*Pm*tf_p)*dtdays
         PHYN1(j) = PHYN1(j) + (dN(j) - delta*muN(j)*PHYN(j) + delta*muN(j1)*PHYN(j1)-PHYN(j)*Pm*tf_p)*dtdays
          CHL1(j) =  CHL1(j) + (dCHL(j)-delta*muCHL(j)*CHL(j)+ delta*muCHL(j1)*CHL(j1)-CHL(j)*Pm*tf_p)*dtdays
        else
         j1 = j - N_alpha * N_Topt
         j2 = j + N_alpha * N_Topt
         PHYC1(j) = PHYC1(j) + (dC(j) - 2d0*delta*muC(j)*PHYC(j) + &
                                          delta*muC(j1)*PHYC(j1) + &
                                          delta*muC(j2)*PHYC(j2) - PHYC(j)*Pm*tf_p)*dtdays

         PHYN1(j) = PHYN1(j) + (dN(j) - 2d0*delta*muN(j)*PHYN(j) + &
                                          delta*muN(j1)*PHYN(j1) + &
                                          delta*muN(j2)*PHYN(j2) - PHYN(j)*Pm*tf_p)*dtdays

         CHL1(j) = CHL1(j) + (dCHL(j)-2d0*delta*muCHL(j)*CHL(j) +  &
                                        delta*muCHL(j1)*CHL(j1) +  &
                                        delta*muCHL(j2)*CHL(j2) - CHL(j)*Pm*tf_p)*dtdays
        endif

        !Store back onto t
        t(iPHYC(j),k) = PHYC1(j) 
        t(iPHYN(j),k) = PHYN1(j) 
        t( iCHL(j),k) =  CHL1(j) 

        Varout(iPHYC(j), k) = t(iPHYC(j),k)
        Varout(iPHYN(j), k) = t(iPHYN(j),k)
        Varout( iCHL(j), k) = t(iCHL(j), k)
      enddo
   Else
      !No mutation
      PHYC1(1) = PHYC1(1) + (dC(  1)-tf_p*Pm*PHYC(1))*dtdays
      PHYN1(1) = PHYN1(1) + (dN(  1)-Pmort          )*dtdays
       CHL1(1) =  CHL1(1) + (dCHL(1)-tf_p*Pm*CHL1(1))*dtdays

      !Store back onto t
      t(iPHYC(1),k) = PHYC1(1) 
      t(iPHYN(1),k) = PHYN1(1) 
      t( iCHL(1),k) =  CHL1(1) 
      Varout(iPHYC(1), k) = t(iPHYC(1),k)
      Varout(iPHYN(1), k) = t(iPHYN(1),k)
      Varout( iCHL(1), k) = t( iCHL(1),k)

   Endif !End of if (NPHY > 1)

   !Now calculate NO3 and DET
   dNO3dt = pp_ND + RES - uptake
   t(iNO3,k) = NO3 + dtdays*dNO3dt
   Varout(iNO3, k) = t(iNO3, k)

   dDETdt    = Pmort + pp_DZ - pp_ND
   t(iDET,k) = DET + dtdays*dDETdt
   Varout(iDET,k) = t(iDET,k)

ENDDO
END SUBROUTINE BIOLOGY