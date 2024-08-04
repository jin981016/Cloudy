c23456 NGC6884 fitting H beta

c      H beta line
       cwvl_hb = 4861.0
       p_hb = 2.1e-11
       w_hb = 0.45

c      He II 4859
       cwvl_he = 4858.8
       p_he = 3.1e-13
       w_he = 0.4

c      Raman He II 4851
       cwvl_r = 4849.74 
       p_r = p_he*0.15
       w_r = w_he*(4851.0/972.0)

       
ccc
       nmax = 1000
       wmin = 4800
       wmax = 4900

       do i = 1, nmax
        wvl = wmin + (wmax - wmin)*real(i)/real(nmax)
        
        fit_hb = p_hb/exp((wvl-cwvl_hb)**2/(w_hb**2))
        fit_he = p_he/exp((wvl-cwvl_he)**2/(w_he**2))
        fit_r = p_r/exp((wvl-cwvl_r)**2/(w_r**2))
c       print*, fit_r
c
c
        fit = fit_hb + fit_he + fit_r 
        f_cont = 3.5e-14
        tot_fit = fit + f_cont
c    
         print*, wvl, tot_fit
        enddo 
        stop
        end
