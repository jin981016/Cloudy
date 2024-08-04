c23456 NGC6881 fitting
c     H alpha line
      cwvl_ha = 6562.21
      p_ha = 9.2e-12
      w_ha = 0.55 
c
c     He II 6560
c
      cwvl_he60 = 6559.56
      p_he60 = 2.9e-13
      w_he60 = 0.4

c     He II 6527
      cwvl_he27 = 6526.58
      p_he27 = 1.1e-14
      w_he27 = 0.45

c     N II 6548
      cwvl_n48_1 = 6547.23
      cwvl_n48_2 = 6547.75
      p_n48_1 = 1.9e-12
      p_n48_2 = 1.6e-12
      w_n48_1 = 0.27
      w_n48_2 = 0.32

c     N II 6583
      del_w = 35.40
      cwvl_n83_1 = cwvl_n48_1 + del_w
      cwvl_n83_2 = cwvl_n48_2 + del_w
      p_n83_1 = p_n48_1 * 3.1
      p_n83_2 = p_n48_2 * 3.1
      w_n83_1 = w_n48_1
      w_n83_2 = w_n48_2

c     C II 6578
      cwvl_c78 = 6575.57
      p_c78 = 6.4e-14
      w_c78 = 0.1


c     H alpha wing
      del_w = 1.4
      cwvl_hw = cwvl_ha +del_w
      p_hw = p_ha * 1.5e-4
      w_hw = w_ha * 5.

c     Raman He II 6545
      del_r = 21.5
      cwvl_rm = cwvl_he27 + del_r
      p_rm = p_he27 * 0.5
      w_rm = w_he27 * 6545./1025.
c
c
c
      nmax = 3000
      wmin = 6500.
      wmax = 6590.

      do i=1, nmax
       wvl = wmin + (wmax - wmin)*real(i)/real(nmax)
c  H alpha 6563
       tmp_ha = (wvl-cwvl_ha)**2/(w_ha**2)
       fit_ha = p_ha /exp(tmp_ha)

c  He II 6560
       tmp_he60 = (wvl-cwvl_he60)**2/(w_he60**2)
       fit_he60 = p_he60 /exp(tmp_he60)

c  He II 6527
       tmp_he27 = (wvl-cwvl_he27)**2/(w_he27**2)
       fit_he27 = p_he27 /exp(tmp_he27)

c  N II 6548_1

       tmp_n48_1 = (wvl-cwvl_n48_1)**2/(w_n48_1**2)
       fit_n48_1 = p_n48_1 /exp(tmp_n48_1)


c  N II 6548_2

       tmp_n48_2 = (wvl-cwvl_n48_2)**2/(w_n48_2**2)
       fit_n48_2 = p_n48_2 /exp(tmp_n48_2)

c  N II 6583_1

       tmp_n83_1 = (wvl-cwvl_n83_1)**2/(w_n83_1**2)
       fit_n83_1 = p_n83_1 /exp(tmp_n83_1)

c  N II 6583_2

       tmp_n83_2 = (wvl-cwvl_n83_2)**2/(w_n83_2**2)
       fit_n83_2 = p_n83_2 /exp(tmp_n83_2)

c    CII 6578
       tmp_c78 = (wvl-cwvl_c78)**2/(w_c78**2)
       fit_c78 = p_c78 /exp(tmp_c78)

c    H alpha Gaussian wing
       tmp_hw = (wvl-cwvl_hw)**2/(w_hw**2)
       fit_hw = p_hw /exp(tmp_hw)

c    Raman He II 6527
       tmp_rm = (wvl-cwvl_rm)**2/(w_rm**2)
       fit_rm = p_rm /exp(tmp_rm)

c       print*, wvl, tmp_hw, fit_hw

       fit1 = fit_ha+fit_he60 + fit_he27
       fit2 = fit_n48_1 + fit_n48_2 + fit_c78 + fit_n83_1 + fit_n83_2
       tot_fit = fit1 + fit2 + fit_hw + fit_rm
       f_cont = 5e-15
       tot_fit =  tot_fit + f_cont
       print*, wvl, tot_fit
       end do
       stop
       end
