c23456
c     H alpha line
      cwvl_ha = 6561.92
      p_ha = 1.e-10
      w_ha = 0.5 
c
c   He II 6560
c
      cwvl_he60 = 6559.2
      p_he60 = 1.4e-12
      w_he60 = 0.6

c    He II 6527
      cwvl_he27 = 6526.2
      p_he27 = 8.e-14
      w_he27 = 0.6

c    N2 6548
      cwvl_n48_1 = 6546.72
      cwvl_n48_2 = 6547.48
      p_n48_1 = 5.4e-12
      p_n48_2 = 2.4e-12
      w_n48_1 = 0.4
      w_n48_2 = 0.25

c   N2 6583
      del_w = 35.40
      cwvl_n83_1 = cwvl_n48_1 + del_w
      cwvl_n83_2 = cwvl_n48_2 + del_w
      p_n83_1 = p_n48_1 * 3.3
      p_n83_2 = p_n48_2 * 3.3
      w_n83_1 = w_n48_1
      w_n83_2 = w_n48_2

c   CII 6578
      cwvl_c78 = 6577.25
      p_c78 = 1.2e-13
      w_c78 = 0.39


c   H alpha wing
      del_w = 1.5
      cwvl_hw = cwvl_ha +del_w
      p_hw = p_ha * 1.3e-4
      w_hw = w_ha * 20.

c   Raman He II 6545
      del_r = 21.
      cwvl_rm = cwvl_he27 + del_r
      p_rm = p_he27 * 0.16
      w_rm = w_he27 * 6545./1025.
c
c
c
      nmax = 1000
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

       fit =  fit_he60 + fit_he27 + fit_n48_1 +fit_n48_2 + fit_c78
       tot_fit = fit + fit_n83_1 + fit_n83_2 +fit_hw + fit_rm
       f_cont = 3.5e-14
       tot_fit =  tot_fit + f_cont
       print*, wvl, tot_fit
       end do
       stop
       end
