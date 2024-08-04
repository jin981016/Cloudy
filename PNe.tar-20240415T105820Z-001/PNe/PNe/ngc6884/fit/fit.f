c23456
c     center wavelengths
      cwvl_ha = 6561.92
      cwvl_he60 = 6559.2
      cwvl_he27 = 6526.2
      cwvl_n48_1 = 6546.72
      cwvl_n48_2 = 6547.48
c   n2 6583
      del_w = 35.40
      cwvl_n83_1 = cwvl_n48_1 + del_w
      cwvl_n83_2 = cwvl_n48_2 + del_w

c
c     flux peak
      p_ha = 1.e-10
      p_he60 = 1.4e-12
      p_he27 = 8.e-14
      p_n48_1 = 5.4e-12
      p_n48_2 = 2.4e-12
c   n2 6583
      p_n83_1 = p_n48_1 * 3.3
      p_n83_2 = p_n48_2 * 3.3


c     profile width 
      w_ha = 0.5 
      w_he60 = 0.6
      w_he27 = 0.6
      w_n48_1 = 0.4
      w_n48_2 = 0.25
      w_n83_1 = w_n48_1
      w_n83_2 = w_n48_2


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
c       print*, wvl, tmp_n83_2, fit_n83_2



       fit =  fit_he60 + fit_he27 + fit_n48_1 +fit_n48_2
       tot_fit = fit + fit_n83_1 + fit_n83_2
       f_cont = 3.5e-14
       tot_fit =  tot_fit + f_cont
       print*, wvl, tot_fit
       end do
       stop
       end
