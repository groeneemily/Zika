library(devtools)
devtools::install_github("DARTH-git/dampack")
library(dampack)

#recoded Markov Model in according to dampack vignette: dsa_generation
rm(list = ls())

decision_model <- function(l_params_all, verbose = FALSE) {
  with(as.list(l_params_all), {
    # compute internal paramters as a function of external parameters
    pDie_background <-scan("p_asr.csv") #Probability to die in health
    pDie_mic <-scan("p_asr_mic.csv") #probability to die with microcephaly
    pDie_gbs <-scan("p_asr_gbs.csv") #probability to die with GBS
    n_t<-100
    d_c <- 0.03
    d_e <- 0.03 
    ####### RUN DECISION TREE#######################################
    #Probability of Zika infection and outcomes under each strategy
    p_Mosquito_spraying <-((1-p_mosq_elim)*(((1-p_mosq_infect)*((1*p_mosq_infect*((p_sex_transmit*1))+((1-p_sex_transmit)*0))+(1-p_mosq_infect)*0))+(p_mosq_infect*1)))
    p_Limiting_movement <-(((1-p_mosq_infect)*(1-(p_mosq_infect*p_limiting_reduction)))*((1*p_mosq_infect*p_limiting_reduction)*((p_sex_transmit*1))+((1-p_sex_transmit)*0)+((1-p_mosq_infect)*0)))+(p_mosq_infect*p_limiting_reduction)
    p_Zika_testing <-((((((((((1-p_abort_pos_test)*(p_sensitivity)+(1-p_sensitivity))*p_test_symp)+(1-p_test_symp))*p_show_symp)+(1-p_show_symp))*p_sex_transmit)+((1-p_sex_transmit)*0))*p_mosq_infect)*(1-p_mosq_infect))+(((((((1-p_abort_pos_test)*(p_sensitivity)+(1-p_sensitivity))*p_test_symp)+(1-p_test_symp))*p_show_symp)+(1-p_show_symp))*p_mosq_infect)
    p_Condom_provision <-((((((p_sex_transmit*p_condom_fail))*(p_mosq_infect))*(1-p_mosq_infect))+p_mosq_infect)*p_condom_use)+((((1-p_mosq_infect)*(p_mosq_infect*p_sex_transmit))+p_mosq_infect)*(1-p_condom_use))
    p_Do_nothing <-(p_sex_transmit*p_mosq_infect)*(1-p_mosq_infect)+p_mosq_infect
    p_Mic_s1 <- p_Mosquito_spraying*e_mic
    p_Gbs_s1 <- p_Mosquito_spraying*pGbs
    p_OA_s1 <- p_Mosquito_spraying*pOA
    p_D_s1 <- p_Mosquito_spraying*pDie_birth
    p_ZH_s1 <- p_Mosquito_spraying*(1-(e_mic+pGbs+pOA+pDie_birth))
    p_H_s1 <- 1-(p_Mic_s1+p_Gbs_s1+p_OA_s1+p_D_s1+p_ZH_s1)
    p_Mic_s2 <- p_Limiting_movement*e_mic
    p_Gbs_s2 <- p_Limiting_movement*pGbs
    p_OA_s2 <- p_Limiting_movement*pOA
    p_D_s2 <- p_Limiting_movement*pDie_birth
    p_ZH_s2 <- p_Limiting_movement*(1-(e_mic+pGbs+pOA+pDie_birth))
    p_H_s2 <- 1-(p_Mic_s2+p_Gbs_s2+p_OA_s2+p_D_s2+p_ZH_s2)
    p_Mic_s3 <- p_Zika_testing*e_mic
    p_Gbs_s3<-p_Zika_testing*pGbs
    p_OA_s3<-p_Zika_testing*pOA
    p_D_s3<-p_Zika_testing*pDie_birth
    p_ZH_s3<-p_Zika_testing*(1-(e_mic+pGbs+pOA+pDie_birth))
    p_H_s3<-1-(p_Mic_s3+p_Gbs_s3+p_OA_s3+p_D_s3+p_ZH_s3)
    p_Mic_s4<-p_Condom_provision*e_mic
    p_Gbs_s4<-p_Condom_provision*pGbs
    p_OA_s4<-p_Condom_provision*pOA
    p_D_s4<-p_Condom_provision*pDie_birth
    p_ZH_s4<-p_Condom_provision*(1-(e_mic+pGbs+pOA+pDie_birth))
    p_H_s4<-1-(p_Mic_s4+p_Gbs_s4+p_OA_s4+p_D_s4+p_ZH_s4)
    p_Mic_s5<-p_Do_nothing*e_mic
    p_Gbs_s5<-p_Do_nothing*pGbs
    p_OA_s5<-p_Do_nothing*pOA
    p_D_s5<-p_Do_nothing*pDie_birth
    p_ZH_s5<-p_Do_nothing*(1-(e_mic+pGbs+pOA+pDie_birth))
    p_H_s5<-1-(p_Mic_s5+p_Gbs_s5+p_OA_s5+p_D_s5+p_ZH_s5)
    
    ####### INITIALIZATION ##########################################
    #
    # create the cohort trace
    m_M1 <- matrix(NA, nrow = n_t + 1 , 
                  ncol = n_s,
                  dimnames = list(0:n_t, v_n))     # create Markov trace (n_t + 1 because R doesn't understand  Cycle 0)
    # create the cohort trace
    m_M2 <- matrix(NA, nrow = n_t + 1 , 
                  ncol = n_s,
                  dimnames = list(0:n_t, v_n))     # create Markov trace (n_t + 1 because R doesn't understand  Cycle 0)
    # create the cohort trace
    m_M3 <- matrix(NA, nrow = n_t + 1 , 
                  ncol = n_s,
                  dimnames = list(0:n_t, v_n))     # create Markov trace (n_t + 1 because R doesn't understand  Cycle 0)
    # create the cohort trace
    m_M4 <- matrix(NA, nrow = n_t + 1 , 
                  ncol = n_s,
                  dimnames = list(0:n_t, v_n))     # create Markov trace (n_t + 1 because R doesn't understand  Cycle 0)
    # create the cohort trace
    m_M5<- matrix(NA, nrow = n_t + 1 , 
                  ncol = n_s,
                  dimnames = list(0:n_t, v_n))     # create Markov trace (n_t + 1 because R doesn't understand  Cycle 0)
    
    
    
    
    
    #m_M[1, ] <- c(1,0,0,0,0,0)                      # initialize Markov trace
    #Calculate initial states of the transition matrix under each strategy
    m_M1[1, ] <- c(p_Mic_s1, p_Gbs_s1, p_OA_s1, p_D_s1, p_ZH_s1, p_H_s1)
    m_M2[1, ] <- c(p_Mic_s2, p_Gbs_s2, p_OA_s2, p_D_s2, p_ZH_s2, p_H_s2)
    m_M3[1, ] <- c(p_Mic_s3, p_Gbs_s3, p_OA_s3, p_D_s3, p_ZH_s3, p_H_s3)
    m_M4[1, ] <- c(p_Mic_s4, p_Gbs_s4, p_OA_s4, p_D_s4, p_ZH_s4, p_H_s4)
    m_M5[1, ] <- c(p_Mic_s5, p_Gbs_s5, p_OA_s5, p_D_s5, p_ZH_s5, p_H_s5)
    
    # create transition probability matrix for all strategies
    m_P1 <- array(0, dim = c(n_s, n_s, 100),
                  dimnames = list(v_n, v_n, 0:99))
    # create transition probability matrix for all strategies
    m_P2 <- array(0, dim = c(n_s, n_s, 100),
                  dimnames = list(v_n, v_n, 0:99))
    # create transition probability matrix for all strategies
    m_P3 <- array(0, dim = c(n_s, n_s, 100),
                  dimnames = list(v_n, v_n, 0:99))
    # create transition probability matrix for all strategies
    m_P4 <- array(0, dim = c(n_s, n_s, 100),
                  dimnames = list(v_n, v_n, 0:99))
    # create transition probability matrix for all strategies
    m_P5 <- array(0, dim = c(n_s, n_s, 100),
                  dimnames = list(v_n, v_n, 0:99))
    # fill in the transition probability array
    ### From Microcephaly
    m_P1["M", "M",] <- 1-pDie_mic
    m_P1["M", "D",] <- pDie_mic
    m_P1["M", "OA",] <-0
    m_P1["M", "GBS",] <-0
    m_P1["M", "ZH",] <-0
    m_P1["M", "H",] <-0
    ### From Other Abnormalities
    m_P1["OA", "OA",] <-1-pDie_background
    m_P1["OA", "D",] <- pDie_background
    m_P1["OA", "M",] <-0
    m_P1["OA", "GBS",]<-0
    m_P1["OA", "ZH",]<-0
    m_P1["OA", "H",]<-0
    ### From with Zika but healthy
    m_P1["ZH", "ZH",] <- 1-pDie_background - pLater_mic
    m_P1["ZH", "D",] <- pDie_background
    m_P1["ZH", "M",] <- pLater_mic
    m_P1["ZH", "OA",] <- 0
    m_P1["ZH", "GBS",] <- 0
    m_P1["ZH", "H",]<-0
    ### From Dead
    m_P1["D", "D",] <-1
    m_P1["D", "M",] <-0 
    m_P1["D", "OA",] <-0
    m_P1["D", "ZH",] <-0
    m_P1["D", "GBS",] <-0
    m_P1["D", "H",] <-0
    ### From GBS
    m_P1["GBS", "M",] <- 0 
    m_P1["GBS", "D",] <- pDie_gbs 
    m_P1["GBS", "OA",] <- 0
    m_P1["GBS", "ZH",] <- 0
    m_P1["GBS", "GBS",] <- 1-pDie_gbs
    m_P1["GBS", "H",] <-0
    ### From Healthy
    m_P1["H", "M",] <- 0
    m_P1["H", "D",] <- pDie_background
    m_P1["H", "OA",] <- 0
    m_P1["H", "ZH",] <- 0
    m_P1["H", "GBS",] <- 0
    m_P1["H", "H",] <- 1-pDie_background

    # fill in the transition probability array
    ### From Microcephaly
    m_P2["M", "M",] <- 1-pDie_mic
    m_P2["M", "D",] <- pDie_mic
    m_P2["M", "OA",] <-0
    m_P2["M", "GBS",] <-0
    m_P2["M", "ZH",] <-0
    m_P2["M", "H",] <-0
    ### From Other Abnormalities
    m_P2["OA", "OA",] <-1-pDie_background
    m_P2["OA", "D",] <- pDie_background
    m_P2["OA", "M",] <-0
    m_P2["OA", "GBS",]<-0
    m_P2["OA", "ZH",]<-0
    m_P2["OA", "H",]<-0
    ### From with Zika but healthy
    m_P2["ZH", "ZH",] <- 1-pDie_background - pLater_mic
    m_P2["ZH", "D",] <- pDie_background
    m_P2["ZH", "M",] <- pLater_mic
    m_P2["ZH", "OA",] <- 0
    m_P2["ZH", "GBS",] <- 0
    m_P2["ZH", "H",]<-0
    ### From Dead
    m_P2["D", "D",] <-1
    m_P2["D", "M",] <-0 
    m_P2["D", "OA",] <-0
    m_P2["D", "ZH",] <-0
    m_P2["D", "GBS",] <-0
    m_P2["D", "H",] <-0
    ### From GBS
    m_P2["GBS", "M",] <- 0 
    m_P2["GBS", "D",] <- pDie_gbs 
    m_P2["GBS", "OA",] <- 0
    m_P2["GBS", "ZH",] <- 0
    m_P2["GBS", "GBS",] <- 1-pDie_gbs
    m_P2["GBS", "H",] <-0
    ### From Healthy
    m_P2["H", "M",] <- 0
    m_P2["H", "D",] <- pDie_background
    m_P2["H", "OA",] <- 0
    m_P2["H", "ZH",] <- 0
    m_P2["H", "GBS",] <- 0
    m_P2["H", "H",] <- 1-pDie_background
    
    # fill in the transition probability array
    ### From Microcephaly
    m_P3["M", "M",] <- 1-pDie_mic
    m_P3["M", "D",] <- pDie_mic
    m_P3["M", "OA",] <-0
    m_P3["M", "GBS",] <-0
    m_P3["M", "ZH",] <-0
    m_P3["M", "H",] <-0
    ### From Other Abnormalities
    m_P3["OA", "OA",] <-1-pDie_background
    m_P3["OA", "D",] <- pDie_background
    m_P3["OA", "M",] <-0
    m_P3["OA", "GBS",]<-0
    m_P3["OA", "ZH",]<-0
    m_P3["OA", "H",]<-0
    ### From with Zika but healthy
    m_P3["ZH", "ZH",] <- 1-pDie_background - pLater_mic
    m_P3["ZH", "D",] <- pDie_background
    m_P3["ZH", "M",] <- pLater_mic
    m_P3["ZH", "OA",] <- 0
    m_P3["ZH", "GBS",] <- 0
    m_P3["ZH", "H",]<-0
    ### From Dead
    m_P3["D", "D",] <-1
    m_P3["D", "M",] <-0 
    m_P3["D", "OA",] <-0
    m_P3["D", "ZH",] <-0
    m_P3["D", "GBS",] <-0
    m_P3["D", "H",] <-0
    ### From GBS
    m_P3["GBS", "M",] <- 0 
    m_P3["GBS", "D",] <- pDie_gbs 
    m_P3["GBS", "OA",] <- 0
    m_P3["GBS", "ZH",] <- 0
    m_P3["GBS", "GBS",] <- 1-pDie_gbs
    m_P3["GBS", "H",] <-0
    ### From Healthy
    m_P3["H", "M",] <- 0
    m_P3["H", "D",] <- pDie_background
    m_P3["H", "OA",] <- 0
    m_P3["H", "ZH",] <- 0
    m_P3["H", "GBS",] <- 0
    m_P3["H", "H",] <- 1-pDie_background
    
    # fill in the transition probability array
    ### From Microcephaly
    m_P4["M", "M",] <- 1-pDie_mic
    m_P4["M", "D",] <- pDie_mic
    m_P4["M", "OA",] <-0
    m_P4["M", "GBS",] <-0
    m_P4["M", "ZH",] <-0
    m_P4["M", "H",] <-0
    ### From Other Abnormalities
    m_P4["OA", "OA",] <-1-pDie_background
    m_P4["OA", "D",] <- pDie_background
    m_P4["OA", "M",] <-0
    m_P4["OA", "GBS",]<-0
    m_P4["OA", "ZH",]<-0
    m_P4["OA", "H",]<-0
    ### From with Zika but healthy
    m_P4["ZH", "ZH",] <- 1-pDie_background - pLater_mic
    m_P4["ZH", "D",] <- pDie_background
    m_P4["ZH", "M",] <- pLater_mic
    m_P4["ZH", "OA",] <- 0
    m_P4["ZH", "GBS",] <- 0
    m_P4["ZH", "H",]<-0
    ### From Dead
    m_P4["D", "D",] <-1
    m_P4["D", "M",] <-0 
    m_P4["D", "OA",] <-0
    m_P4["D", "ZH",] <-0
    m_P4["D", "GBS",] <-0
    m_P4["D", "H",] <-0
    ### From GBS
    m_P4["GBS", "M",] <- 0 
    m_P4["GBS", "D",] <- pDie_gbs 
    m_P4["GBS", "OA",] <- 0
    m_P4["GBS", "ZH",] <- 0
    m_P4["GBS", "GBS",] <- 1-pDie_gbs
    m_P4["GBS", "H",] <-0
    ### From Healthy
    m_P4["H", "M",] <- 0
    m_P4["H", "D",] <- pDie_background
    m_P4["H", "OA",] <- 0
    m_P4["H", "ZH",] <- 0
    m_P4["H", "GBS",] <- 0
    m_P4["H", "H",] <- 1-pDie_background
    
    # fill in the transition probability array
    ### From Microcephaly
    m_P5["M", "M",] <- 1-pDie_mic
    m_P5["M", "D",] <- pDie_mic
    m_P5["M", "OA",] <-0
    m_P5["M", "GBS",] <-0
    m_P5["M", "ZH",] <-0
    m_P5["M", "H",] <-0
    ### From Other Abnormalities
    m_P5["OA", "OA",] <-1-pDie_background
    m_P5["OA", "D",] <- pDie_background
    m_P5["OA", "M",] <-0
    m_P5["OA", "GBS",]<-0
    m_P5["OA", "ZH",]<-0
    m_P5["OA", "H",]<-0
    ### From with Zika but healthy
    m_P5["ZH", "ZH",] <- 1-pDie_background - pLater_mic
    m_P5["ZH", "D",] <- pDie_background
    m_P5["ZH", "M",] <- pLater_mic
    m_P5["ZH", "OA",] <- 0
    m_P5["ZH", "GBS",] <- 0
    m_P5["ZH", "H",]<-0
    ### From Dead
    m_P5["D", "D",] <-1
    m_P5["D", "M",] <-0 
    m_P5["D", "OA",] <-0
    m_P5["D", "ZH",] <-0
    m_P5["D", "GBS",] <-0
    m_P5["D", "H",] <-0
    ### From GBS
    m_P5["GBS", "M",] <- 0 
    m_P5["GBS", "D",] <- pDie_gbs 
    m_P5["GBS", "OA",] <- 0
    m_P5["GBS", "ZH",] <- 0
    m_P5["GBS", "GBS",] <- 1-pDie_gbs
    m_P5["GBS", "H",] <-0
    ### From Healthy
    m_P5["H", "M",] <- 0
    m_P5["H", "D",] <- pDie_background
    m_P5["H", "OA",] <- 0
    m_P5["H", "ZH",] <- 0
    m_P5["H", "GBS",] <- 0
    m_P5["H", "H",] <- 1-pDie_background
    
    # check rows add up to 1
    for (i in 1:100) {
      m_P1_i <- as.matrix(m_P1[,, i])
      if (!isTRUE(all.equal(as.numeric(rowSums(m_P1_i)), as.numeric(rep(1, n_s))))) {
        stop("This is not a valid transition Matrix")
      }
    }
    for (i in 1:100) {
      m_P2_i <- as.matrix(m_P2[,, i])
      if (!isTRUE(all.equal(as.numeric(rowSums(m_P2_i)), as.numeric(rep(1, n_s))))) {
        stop("This is not a valid transition Matrix")
      }
    }
    for (i in 1:100) {
      m_P3_i <- as.matrix(m_P3[,, i])
      if (!isTRUE(all.equal(as.numeric(rowSums(m_P3_i)), as.numeric(rep(1, n_s))))) {
        stop("This is not a valid transition Matrix")
      }
    }
    for (i in 1:100) {
      m_P4_i <- as.matrix(m_P4[,, i])
      if (!isTRUE(all.equal(as.numeric(rowSums(m_P4_i)), as.numeric(rep(1, n_s))))) {
        stop("This is not a valid transition Matrix")
      }
    }
    for (i in 1:100) {
      m_P5_i <- as.matrix(m_P5[,, i])
      if (!isTRUE(all.equal(as.numeric(rowSums(m_P5_i)), as.numeric(rep(1, n_s))))) {
        stop("This is not a valid transition Matrix")
      }
    }
  
      ############# PROCESS ###########################################
    
    for (t in 1:n_t){                              # throughout the number of cycles (Strategy 1)
       m_M1[t + 1, ] <- m_M1[t, ] %*% m_P1[,, t]           # estimate the Markov trace for cycle the next cycle (t + 1)
    }
    for (t in 1:n_t){                              # throughout the number of cycles (Strategy 2)
      m_M2[t + 1, ] <- m_M2[t, ] %*% m_P2[,, t]           # estimate the Markov trace for cycle the next cycle (t + 1)
    }
    for (t in 1:n_t){                              # throughout the number of cycles (Strategy 3)
      m_M3[t + 1, ] <- m_M3[t, ] %*% m_P3[,, t]           # estimate the Markov trace for cycle the next cycle (t + 1)
    }
    for (t in 1:n_t){                              # throughout the number of cycles (Strategy 4)
      m_M4[t + 1, ] <- m_M4[t, ] %*% m_P4[,, t]           # estimate the Markov trace for cycle the next cycle (t + 1)
    }
    for (t in 1:n_t){                              # throughout the number of cycles (Strategy 5)
      m_M5[t + 1, ] <- m_M5[t, ] %*% m_P5[,, t]           # estimate the Markov trace for cycle the next cycle (t + 1)
    }
    
    ####### EPIDEMIOLOGICAL OUTPUT  ###########################################
    #### Overall Survival (OS) ####
    v_os1 <- 1 - m_M1[, "D"]                # calculate the overall survival (OS) probability for Strategy 1
    v_os2 <- 1 - m_M2[, "D"]                # calculate the overall survival (OS) probability for Strategy 2
    v_os3 <- 1 - m_M3[, "D"]                # calculate the overall survival (OS) probability for Strategy 3
    v_os4 <- 1 - m_M4[, "D"]                # calculate the overall survival (OS) probability for Strategy 4
    v_os5 <- 1 - m_M5[, "D"]                # calculate the overall survival (OS) probability for Strategy 5
    
    #### Microcephaly prevalence #####
    # v_prev_M <- rowSums(m_M[, c("M")])/v_os
    v_prev_M1 <- m_M1[, c("M")]/v_os1
    v_prev_M2 <- m_M2[, c("M")]/v_os2
    v_prev_M3 <- m_M3[, c("M")]/v_os3
    v_prev_M4 <- m_M4[, c("M")]/v_os4
    v_prev_M5 <- m_M5[, c("M")]/v_os5
    # v_prev_GBS<- rowSums(m_M[, c("GBS")])/v_os
    v_prev_GBS1<- m_M1[, c("GBS")]/v_os1
    v_prev_GBS2<- m_M2[, c("GBS")]/v_os2
    v_prev_GBS3<- m_M3[, c("GBS")]/v_os3
    v_prev_GBS4<- m_M4[, c("GBS")]/v_os4
    v_prev_GBS5<- m_M5[, c("GBS")]/v_os5
    
    # v_prev_OA<- rowSums(m_M[, c("OA")])/v_os
    v_prev_OA1<- m_M1[, c("OA")]/v_os1
    v_prev_OA2<- m_M2[, c("OA")]/v_os2
    v_prev_OA3<- m_M3[, c("OA")]/v_os3
    v_prev_OA4<- m_M4[, c("OA")]/v_os4
    v_prev_OA5<- m_M5[, c("OA")]/v_os5
    
    ####### RETURN OUTPUT  ###########################################
    out1 <- list(m_M1 = m_M1,
                m_P1 = m_P1,
                Surv1 = v_os1[-1],
                Prev_M1 = v_prev_M1[-1],
                Prev_GBS1 = v_prev_GBS1[-1],
                Prev_OA1 = v_prev_OA1[-1])
                #PropSick = v_prop_S1[c(11, 21, 31)]
    
    return(out1)
    out2 <- list(m_M2 = m_M2,
                m_P2 = m_P2,
                Surv2 = v_os2[-1],
                Prev_M2 = v_prev_M2[-1],
                Prev_GBS2 = v_prev_GBS2[-1],
                Prev_OA2 = v_prev_OA2[-1])
    #PropSick = v_prop_S1[c(11, 21, 31)]
    
    return(out2)
    out3 <- list(m_M3 = m_M3,
                m_P3 = m_P3,
                Surv3 = v_os3[-1],
                Prev_M3 = v_prev_M3[-1],
                Prev_GBS3 = v_prev_GBS3[-1],
                Prev_OA3 = v_prev_OA3[-1])
    #PropSick = v_prop_S1[c(11, 21, 31)]
    
    return(out3)
    out4 <- list(m_M4 = m_M4,
                m_P4 = m_P4,
                Surv4 = v_os4[-1],
                Prev_M4 = v_prev_M4[-1],
                Prev_GBS4 = v_prev_GBS4[-1],
                Prev_OA4 = v_prev_OA4[-1])
    #PropSick = v_prop_S1[c(11, 21, 31)]
    
    return(out4)
    out5 <- list(m_M5 = m_M5,
                m_P5 = m_P5,
                Surv5 = v_os5[-1],
                Prev_M5 = v_prev_M5[-1],
                Prev_GBS5 = v_prev_GBS5[-1],
                Prev_OA5 = v_prev_OA5[-1])
    #PropSick = v_prop_S1[c(11, 21, 31)]
    
    return(out5)
    
  }
  )
}
calculate_ce_out <- function(l_params_all, n_wtp = 4250262){
  # User defined
  with(as.list(l_params_all), {
    
    #probabilities of each path) times the reward of each path
    #costs
    c_Mosquito_spraying <-c_spray+c_prevent_total
    c_Limiting_movement <-c_limit+c_prevent_total
    c_Zika_testing <-c_test_total+c_prevent_total
    c_Condom_provision <-c_condom_pop+c_prevent_total
    c_Do_nothing <-c_prevent_total
    
    #and effects
    e_Mosquito_spraying_inf <- p_Mosquito_spraying*e_inf # Mosquito spraying Zika infection
    e_Limiting_movement_inf <- p_Limiting_movement*e_inf   # Limiting movement Zika infection
    e_Zika_testing_inf <- p_Zika_testing*e_inf # Zika testing Zika infection
    e_Condom_provision_inf <-p_Condom_provision*e_inf #Condom provision Zika infection
    e_Do_nothing_inf <-p_Do_nothing*e_inf #Do nothing Zika infection
    
    e_Mosquito_spraying_mic <- p_Mosquito_spraying*e_mic # Mosquito spraying Zika infection
    e_Limiting_movement_mic <- p_Limiting_movement*e_mic   # Limiting movement Zika infection
    e_Zika_testing_mic <- p_Zika_testing*e_mic # Zika testing Zika infection
    e_Condom_provision_mic <-p_Condom_provision*e_mic #Condom provision Zika infection
    e_Do_nothing_mic <-p_Do_nothing*e_mic #Do nothing Zika infection
    
    ## Create discounting vectors
    v_dwc <- 1 / ((1 + d_c) ^ (0:(n_t))) # vector with discount weights for costs
    v_dwe <- 1 / ((1 + d_e) ^ (0:(n_t))) # vector with discount weights for QALYs
    
    ## Run STM model at a parameter set for each intervention
    l_model_out_s1 <- decision_model(l_params_all = l_params_all)
    l_model_out_s2    <- decision_model(l_params_all = l_params_all)
    l_model_out_s3    <- decision_model(l_params_all = l_params_all)
    l_model_out_s4    <- decision_model(l_params_all = l_params_all)
    l_model_out_s5    <- decision_model(l_params_all = l_params_all)
    
    ## Cohort trace by treatment 
    m_M_s1 <- l_model_out_s1$m_M1 # Mosquito spraying
    m_M_s2 <- l_model_out_s2$m_M2 # Limiting movement
    m_M_s3 <- l_model_out_s3$m_M3 # Zika Testing
    m_M_s4 <- l_model_out_s4$m_M4 # Condom Provision
    m_M_s5 <- l_model_out_s5$m_M5 # Do nothing
    
      ## Vectors with costs and utilities by treatment
    v_u_s1 <- c(u_Mic, u_Gbs, u_OA, u_D, u_ZH, u_H)
    v_u_s2 <- c(u_Mic, u_Gbs, u_OA, u_D, u_ZH, u_H)
    v_u_s3 <- c(u_Mic, u_Gbs, u_OA, u_D, u_ZH, u_H)
    v_u_s4 <- c(u_Mic, u_Gbs, u_OA, u_D, u_ZH, u_H)
    v_u_s5 <- c(u_Mic, u_Gbs, u_OA, u_D, u_ZH, u_H)
    
    v_c_s1 <- c(c_Mic, c_Gbs, c_OA, c_D, c_ZH, c_H)
    v_c_s2 <- c(c_Mic, c_Gbs, c_OA, c_D, c_ZH, c_H)
    v_c_s3 <- c(c_Mic, c_Gbs, c_OA, c_D, c_ZH, c_H)
    v_c_s4 <- c(c_Mic, c_Gbs, c_OA, c_D, c_ZH, c_H)
    v_c_s5 <- c(c_Mic, c_Gbs, c_OA, c_D, c_ZH, c_H)
    
    ## Mean Costs and QALYs for Treatment and NO Treatment
    v_tu_s1 <- m_M_s1 %*% v_u_s1
    v_tu_s2 <- m_M_s2 %*% v_u_s2
    v_tu_s3 <- m_M_s3 %*% v_u_s3
    v_tu_s4 <- m_M_s4 %*% v_u_s4
    v_tu_s5 <- m_M_s5 %*% v_u_s5
    
    v_tc_s1 <- m_M_s1 %*% v_c_s1
    v_tc_s2 <- m_M_s2 %*% v_c_s2
    v_tc_s3 <- m_M_s3 %*% v_c_s3
    v_tc_s4 <- m_M_s4 %*% v_c_s4
    v_tc_s5 <- m_M_s5 %*% v_c_s5

    ## Total discounted mean Costs and QALYs
    tu_d_s1 <- t(v_tu_s1) %*% v_dwe 
    tu_d_s2 <- t(v_tu_s2) %*% v_dwe 
    tu_d_s3 <- t(v_tu_s3) %*% v_dwe 
    tu_d_s4 <- t(v_tu_s4) %*% v_dwe 
    tu_d_s5 <- t(v_tu_s5) %*% v_dwe 
    
    tc_d_s1 <- t(v_tc_s1) %*% v_dwc
    tc_d_s2 <- t(v_tc_s2) %*% v_dwc
    tc_d_s3 <- t(v_tc_s3) %*% v_dwc
    tc_d_s4 <- t(v_tc_s4) %*% v_dwc
    tc_d_s5 <- t(v_tc_s5) %*% v_dwc
    
    
    ## Vector with total discounted mean Costs and QALYs
    v_tc_d <- c(tc_d_s1, tc_d_s2, tc_d_s3, tc_d_s4, tc_d_s5)
    v_tu_d <- c(tu_d_s1, tu_d_s2, tu_d_s3, tu_d_s4, tu_d_s5)
    
    ## Vector with discounted net monetary benefits (NMB)
    v_nmb_d <- v_tu_d * n_wtp - v_tc_d
    
    ## Dataframe with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_tc_d,
                        Effect   = v_tu_d,
                        NMB      = v_nmb_d)
    
    return(df_ce)
  }
  )
}
v_names_str <- c("s1", "s2", "s3", "s4", "s5") 
#Full names of strategies "Mosquito_spraying" - s1, "Limiting_movement" - s2 "Zika_testing" - s3, "Condom_provision" - s4, "Do_nothing" -s5

## Number of strategies
n_str <- length(v_names_str)
## Markov model parameters
n_age_init <- 0     # age at baseline
n_age_max  <- 100   # maximum age of follow up
n_t  <- n_age_max - n_age_init  # time horizon, number of cycles
v_n  <- c("M","GBS","OA","ZH","D","H")               # the 6 states of the model: Microcephaly (M), Guillain-Barre Syndrome (GBS), Other Abnormalities (OA), Zika and Healthy (ZH), Dead (D), Healthy (H)
n_s <- length(v_n)                            # number of health states 


my_params_basecase <- list(pLater_mic = 0.01,
                           pGbs = 0.05,
                           pFirstri_mic = 0.13,
                           pLatetri_mic = 0.05,
                           p_mic = 0.07666,
                           pOA = 0.07,
                           pDie_micbirth = 0.051,
                           pDie_birth = 0.059,
                           pDie_gbs = 0.05,
                           cSpray = 38956000,
                           cLimit = 500000,
                           cTest = 800,
                           cAbort = 509.29,
                           cPos_test = 1317,
                           cCondom = 205.5,
                           c_Care_defects = 179760, #cOA, cMic
                           c_Premature_death = 4250262,
                           c_Gbs = 58419.58,
                           c_Mic = 179760,
                           c_OA = 9723.13,
                           c_ZH = 0,
                           c_D = 0,
                           c_H = 0,
                           u_Mic = 0.5,
                           u_Gbs = 0.87,
                           u_OA = 0.6,
                           u_D = 0,
                           u_H = 1,
                           u_ZH = 0.975,
                           d_e = 0.03,
                           d_c = 0.03,
                           p_mosq_elim = 0.77, 
                           p_mosq_infect = 0.00844, 
                           p_sex_transmit = 0.699906249, 
                           p_limiting_reduction = 0.36, 
                           p_abort_pos_test = 0.05383, 
                           p_sensitivity = 1,
                           p_pos_test_zika = 0.608695652, 
                           p_condom_use = 0.75, 
                           p_condom_fail = 0.43, 
                           p_show_symp = 0.2, 
                           p_test_symp = 0.59, 
                           p_first_tri_mic = 0.13, 
                           p_late_tri_mic = 0.05, 
                           pDie_birth = 0.059,
                           c_spray = 38956000, 
                           c_limit = 500000, 
                           c_test = 1169.07, 
                           c_abort = 509.29, 
                           c_pos_test = 1317,
                           c_condom = 93.21, 
                           c_test_total = 19220, 
                           c_prevent_total = 32839, 
                           c_condom_pop = 16725305.93,
                           e_inf = 1,
                           e_mic = 0.08666667 
                           )
# Calculate costs and effectiveness relative to Do Nothing
#v_c_rel <- v_c - rep(v_c["Do Nothing"],length(v_names_str))
# v_c_rel <- v_c_rel
# v_inf_avert <- -(v_inf - rep(v_inf["Do Nothing"],length(v_names_str)))
# v_mic_avert <- -(v_mic - rep(v_mic["Do Nothing"],length(v_names_str)))

#### 05 Compute ICERs of Decision Tree ####
m_cea <- calculate_icers(df_ce)
m_cea

#### 06 Plot frontier of Decision Tree ####
#front_ce <- getFrontier(CEmat = m_cea, plot = F)
#plot_frontier(CEmat = m_ce, frontier = front_ce)
#ggsave("figs/Zika_DecTree-CEA-Frontier.pdf", width = 8, height = 6)

my_params_range <- data.frame(pars = c("c_limit", "p_limiting_reduction"),
                              min = c(0, 0.01),
                              max = c(20000000, 0.5))
owsa_det <- run_owsa_det(params_range = my_params_range,
                         params_basecase = my_params_basecase,
                         nsamp = 100,
                         FUN = calculate_ce_out,
                         outcomes = c("Cost", "Effect", "NMB"),
                         strategies = c("s1", "s2", "s3", "s4", "s5"))
###when multiple outcomes are specified, individual owsa objects are found within the list output of run_owsa_det
plot(owsa_det$owsa_NMB)
twsa_det <- run_twsa_det(params_range = my_params_range,
                         params_basecase = my_params_basecase,
                         nsamp = 100,
                         FUN = calculate_ce_out,
                         outcomes = c("Cost", "Effect", "NMB"),
                         strategies = c("s1", "s2", "s3", "s4", "s5"))
plot(twsa_det$twsa_NMB)

