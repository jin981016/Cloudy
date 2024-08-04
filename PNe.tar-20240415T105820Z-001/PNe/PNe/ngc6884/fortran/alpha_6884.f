c23456 NGC6884 fitting H alpha

c      H alpha line
       cwvl_ha = 6561.8
       p_ha = 9.3e-11
       w_ha = 0.65

c      He II 6560
       cwvl_he60 = 6559.2
       p_he60 = 1.4e-12
       w_he60 = 0.6

c      He II 6527
       cwvl_he27 = 6526.2
       p_he27 = 8.7e-14
       w_he27 = 0.6

c      Raman He II 6545
       del_r = 21.0
       cwvl_r = cwvl_he27 + del_r
       p_r = p_he27*0.16
       w_r = w_he27*(6545.0/1025.0)

c      H alph wing
       del_w = 1.0
       cwvl_hw = cwvl_ha + del_w
       p_hw = p_ha*1.5e-4
       w_hw = w_ha*16.0
       
ccc
       nmax = 1000
       wmin = 6500
       wmax = 6590

       do i = 1, nmax
        wvl = wmin + (wmax - wmin)*real(i)/real(nmax)
        
        tmp_ha = (wvl-cwvl_ha)**2/(w_ha**2)
        fit_ha = p_ha/exp(tmp_ha)

        tmp_he60 = (wvl-cwvl_he60)**2/(w_he60**2)
        fit_he60 = p_he60/exp(tmp_he60)

        tmp_he27 = (wvl-cwvl_he27)**2/(w_he27**2)
        fit_he27 = p_he27/exp(tmp_he27)

        tmp_r = (wvl-cwvl_r)**2/(w_r**2)
        fit_r = p_r/exp(tmp_r)

        tmp_hw = (wvl-cwvl_hw)**2/(w_hw**2)
        fit_hw = p_hw/exp(tmp_hw)
      
c
c
c        print*, p_ha, fit_ha

         fits = fit_ha + fit_he60 + fit_he27 + fit_r + fit_hw
         f_cont = 3.5e-14 
         tot_fit = fits + f_cont
c    
         print*, wvl, tot_fit
        enddo 
        stop
        end
