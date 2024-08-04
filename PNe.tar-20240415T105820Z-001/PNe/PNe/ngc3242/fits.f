c23456 Hu 2-1 fitting
c      H alpha line
       cwvl_ha = 6563.
       p_ha = 1.9e-11
       w_ha = 0.5

c      He II 6560
       cwvl_he60 = 6560.4
       p_he60 = 4.1e-13
       w_he60 = 0.4

c      He II 6527
       cwvl_he27 = 6527.52
       p_he27 = 1.7e-14
       w_he27 = 0.5

c      N II 6548
       cwvl_n48_1 = 6548.06
       cwvl_n48_2 = 6548.85
       p_n48_1 = 1.1e-13
       p_n48_2 = 2.0e-14
       w_n48_1 = 0.2
       w_n48_2 = 0.2
       
c      N II 6583
       del_w = 35.5
       cwvl_n83_1 = cwvl_n48_1 + del_w
       cwvl_n83_2 = cwvl_n48_2 + del_w 
       p_n83_1 = p_n48_1 * 3.
       p_n83_2 = p_n48_2 * 3.
       w_n83_1 = w_n48_1  
       w_n83_2 = w_n48_2 
       
c      C II 6578
       cwvl_c78 = 6578.
       p_c78 = 4.5e-14
       w_c78 = 0.3
      
c      H alpha wing
       del_w = 0.2
       cwvl_hw = cwvl_ha + del_w
       p_hw = p_ha * 2.e-3
       w_hw = w_ha * 8.

c      Raman He II 6545
       del_r = 21.
       cwvl_rm = cwvl_he27 + del_r
       p_rm = p_he27 * 0.1
       w_rm = w_he27 * 6545./1025.


c      Fitting

       nmax = 1000.
       wmin = 6500.
       wmax = 6590.

       do i = 1, nmax
        wvl = wmin + (wmax - wmin)*real(i)/real(nmax)

c      H alpha 6563
       tmp_ha = (wvl-cwvl_ha)**2/(w_ha**2)
       fit_ha = p_ha/exp(tmp_ha)

c      He II 6560
       tmp_he60 = (wvl-cwvl_he60)**2/(w_he60**2)
       fit_he60 = p_he60/exp(tmp_he60)

c      He II 6527 
       tmp_he27 = (wvl-cwvl_he27)**2/(w_he27**2)
       fit_he27 = p_he27/exp(tmp_he27)

c      N II 6548_1 
       tmp_n48_1 = (wvl-cwvl_n48_1)**2/(w_n48_1**2)
       fit_n48_1 = p_n48_1/exp(tmp_n48_1)
c      N II 6548_2 
       tmp_n48_2 = (wvl-cwvl_n48_2)**2/(w_n48_2**2)
       fit_n48_2 = p_n48_2/exp(tmp_n48_2)

c      N II 6583_1 
       tmp_n83_1 = (wvl-cwvl_n83_1)**2/(w_n83_1**2)
       fit_n83_1 = p_n83_1/exp(tmp_n83_1)
c      N II 6583_2 
       tmp_n83_2 = (wvl-cwvl_n83_2)**2/(w_n83_2**2)
       fit_n83_2 = p_n83_2/exp(tmp_n83_2)

c      C II 6578 
       tmp_c78 = (wvl-cwvl_c78)**2/(w_c78**2)
       fit_c78 = p_c78/exp(tmp_c78)

c      H alpha Gaussian wing
       tmp_hw = (wvl-cwvl_hw)**2/(w_hw**2)
       fit_hw = p_hw/exp(tmp_hw)

c      Raman He II 6527 
       tmp_rm = (wvl-cwvl_rm)**2/(w_rm**2)
       fit_rm = p_rm/exp(tmp_rm)


c      
       fit1 = fit_ha + fit_he60 + fit_he27
       fit2 = fit_n48_1 + fit_n48_2 + fit_n83_1 + fit_n83_2 +fit_c78
       tot_fit = fit1 + fit2 + fit_hw + fit_rm 
       f_cont = 7e-15
       tot_fit = tot_fit + f_cont

       print*, wvl, tot_fit

c       print*, wvl, fit_n83_1
       enddo
       stop
       end







       


       
