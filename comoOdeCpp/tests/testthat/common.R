CORE_FILE <- "/inc/covidage_v19.1.2.R"

load_libraries <- function() {
  library_list <- list(
    "deSolve", # ode
    "dplyr",
    "readxl", # read_excel
    "stringr" # str_replace_all
  )
  for (ll in library_list) {
    if (!requireNamespace(ll, quietly = TRUE)) {
      testthat::skip(paste(ll, "needed but not available"))
    }
    library(ll, character.only = TRUE)
  }
  library("comoOdeCpp")
}

init <- function(e) {
  load_libraries()
  load("data/data_CoMo.RData", envir = e)
}

check_parameters_list_for_na <- function(parameters_list) {
  no_na_found <- TRUE
  for (pp_name in names(parameters_list)) {
    if (is.na(parameters_list[[pp_name]]) || is.nan(parameters_list[[pp_name]])) {
      print(paste0("parameters_list[\"", pp_name, "\"] = ", parameters_list[[pp_name]]), quote = FALSE)
      testthat::expect_equal(is.na(parameters_list[[pp_name]]), FALSE)
    }
  }
  if (!no_na_found) {
    stop()
  }
}

match_outputs <- function(
    output_a,      # output matrix #1
    output_b,      # output matrix #2
    tlr = 0.0001, # tolerance
    smp = 1000    # num samples to take
  ) {

  # for (ii in 1:smp) {
  for (ii in seq_len(smp)) {
    # rr <- sample(1:nrow(output_a), 1)
    # cc <- sample(1:ncol(output_a), 1)
    rr <- sample(seq_len(nrow(output_a)), 1)
    cc <- sample(seq_len(ncol(output_a)), 1)
    # print(paste("output_a[rr,cc] =", output_a[rr,cc]))
    # print(paste("output_b[rr,cc] =", output_b[rr,cc]))

    out_a <- output_a[rr, cc]
    out_b <- output_b[rr, cc]

    testthat::expect_true(is.numeric(out_a))
    testthat::expect_true(is.numeric(out_b))

    testthat::expect_gte(out_a, 0) # >=0
    testthat::expect_gte(out_b, 0) # >=0

    if (out_a > 0) {

      res <- expect_equal(
        out_b,
        out_a,
        tolerance = tlr,
        scale = out_a
      )

      if (abs(out_b - out_a) > out_a * tlr) {
        print(paste(
          "not equal: rr=", rr,
          ", cc=", cc,
          ", pp=", pp,
          ", output_a[rr,cc]", out_a,
          ", output_b[rr,cc]", out_b
        ))
      }

    }
  }

}

match_processed_outputs <- function(
    output_a,    # processed output matrix
    output_b,    # processed output matrix
    tlr = 0.0001 # tolerance
  ) {

  testthat::expect_true(is.numeric(output_a$total_reported_deaths_end))
  testthat::expect_true(is.numeric(output_b$total_reported_deaths_end))

  testthat::expect_equal(
      output_a$total_reported_deaths_end,
      output_b$total_reported_deaths_end,
      tolerance = tlr,
      scale = output_b$total_reported_deaths_end
  )

}

make_spline_funs <- function(
    give,
    beds_available,
    icu_beds_available,
    ventilators_available
  ) {

  f <- c(1, (1 + give) / 2, (1 - give) / 2, 0)
  
  KH <- beds_available
  KICU <- icu_beds_available + ventilators_available
  Kvent <- ventilators_available
  
  x.H    <- c(0, (1 + give) * KH / 2,    (3 - give) * KH / 2,    2 * KH)
  x.ICU  <- c(0, (1 + give) * KICU / 2,  (3 - give) * KICU / 2,  2 * KICU)
  x.Vent <- c(0, (1 + give) * Kvent / 2, (3 - give) * Kvent / 2, 2 * Kvent)
  
  fH    <- splinefun(x.H,    f, method = "hyman")
  fICU  <- splinefun(x.ICU,  f, method = "hyman")
  fVent <- splinefun(x.Vent, f, method = "hyman")

  return(list(
    fH = fH,
    fICU = fICU,
    fVent = fVent
  ))
}

get_incidents <- function(
    parameters,
    out_mat,
    ihr_col,
    times
  ) {

  spline_funs <- make_spline_funs(
      give = parameters["give"],
      beds_available = parameters["beds_available"],
      icu_beds_available = parameters["icu_beds_available"],
      ventilators_available = parameters["ventilators_available"]
  )

  incd_time_series <-
    with(as.list(parameters, spline_funs), {

      critH<-c()
      crit<-c()
      critV<-c()

      for (t in 1:length(times)){
        critH[t] <- min(1-fH((sum(out_mat[t,(Hindex+1)]))+sum(out_mat[t,(ICUCindex+1)])+sum(out_mat[t,(ICUCVindex+1)])),1)
        crit[t]  <- min(1-fICU((sum(out_mat[t,(ICUindex+1)]))+(sum(out_mat[t,(Ventindex+1)]))+(sum(out_mat[t,(VentCindex+1)]))))
        critV[t] <- min(1-fVent((sum(out_mat[t,(Ventindex+1)]))),1)
      }


      # daily incidence
      incidence <-
        report     * gamma * (1 - pclin)    * out_mat[, (Eindex+1)]   %*% (1 - ihr_col) +
        reportc    * gamma * pclin          * out_mat[, (Eindex+1)]   %*% (1 - ihr_col) +
        report     * gamma * (1 - pclin)    * out_mat[, (QEindex+1)]  %*% (1 - ihr_col) +
        reportc    * gamma * pclin          * out_mat[, (QEindex+1)]  %*% (1 - ihr_col) +
        report_v   * gamma * (1 - pclin_v)  * out_mat[, (EVindex+1)]  %*% (1 - sigmaEV * ihr_col) +
        report_cv  * gamma * pclin_v        * out_mat[, (EVindex+1)]  %*% (1 - sigmaEV * ihr_col) +
        report_vr  * gamma * (1 - pclin_vr) * out_mat[, (EVRindex+1)] %*% (1 - sigmaEVR * ihr_col) +
        report_cvr * gamma * pclin_vr       * out_mat[, (EVRindex+1)] %*% (1 - sigmaEVR * ihr_col) +
        report_r   * gamma * (1 - pclin_r)  * out_mat[, (ERindex+1)]  %*% (1 - sigmaER * ihr_col) +
        report_cr  * gamma * pclin_r        * out_mat[, (ERindex+1)]  %*% (1 - sigmaER * ihr_col)
      
        # print("incidence:")
        # print(incidence)

      incidenceh <-
        gamma * out_mat[, (Eindex+1)]  %*% ihr_col * (1 - critH) * (1 - prob_icu) * reporth +
        gamma * out_mat[, (Eindex+1)]  %*% ihr_col * (1 - critH) * (1 - prob_icu) * (1 - reporth) * reporth_g +
        gamma * out_mat[, (QEindex+1)] %*% ihr_col * (1 - critH) * (1 - prob_icu) * reporth +
        gamma * out_mat[, (QEindex+1)] %*% ihr_col * (1 - critH) * (1 - prob_icu) * (1 - reporth) * reporth_g +
        gamma * sigmaEV  * out_mat[, (EVindex+1)]  %*% ihr_col * (1 - critH) * (1 - prob_icu_v)  * reporth +
        gamma * sigmaEVR * out_mat[, (EVRindex+1)] %*% ihr_col * (1 - critH) * (1 - prob_icu_vr) * reporth +
        gamma * sigmaER  * out_mat[, (ERindex+1)]  %*% ihr_col * (1 - critH) * (1 - prob_icu_r)  * reporth +
        gamma * out_mat[, (Eindex+1)]              %*% ihr_col * critH * reporth_g * (1 - prob_icu) +
        gamma * out_mat[, (QEindex+1)]             %*% ihr_col * critH * reporth_g * (1 - prob_icu) +
        gamma * sigmaEV  * out_mat[, (EVindex+1)]  %*% ihr_col * critH * reporth_g * (1 - prob_icu_v) +
        gamma * sigmaEVR * out_mat[, (EVRindex+1)] %*% ihr_col * critH * reporth_g * (1 - prob_icu_vr) +
        gamma * sigmaER  * out_mat[, (ERindex+1)]  %*% ihr_col * critH * reporth_g * (1 - prob_icu_r) +
        gamma * out_mat[, (Eindex+1)]              %*% ihr_col * prob_icu * (1 - crit) * reporth_ICU +
        gamma * out_mat[, (QEindex+1)]             %*% ihr_col * prob_icu * (1 - crit) * reporth_ICU +
        gamma * out_mat[, (Eindex+1)]              %*% ihr_col * prob_icu * crit * reporth_ICU * reporth_g +
        gamma * out_mat[, (QEindex+1)]             %*% ihr_col * prob_icu * crit * reporth_ICU * reporth_g +
        gamma * sigmaEV  * out_mat[, (EVindex+1)]  %*% ihr_col * (1 - crit) * prob_icu_v * reporth_ICU +
        gamma * sigmaEVR * out_mat[, (EVRindex+1)] %*% ihr_col * (1 - crit) * prob_icu_vr * reporth_ICU +
        gamma * sigmaER  * out_mat[, (ERindex+1)]  %*% ihr_col * (1 - crit) * prob_icu_r * reporth_ICU +
        gamma * sigmaEV  * out_mat[, (EVindex+1)]  %*% ihr_col * crit * prob_icu_v  * reporth_ICU * reporth_g +
        gamma * sigmaEVR * out_mat[, (EVRindex+1)] %*% ihr_col * crit * prob_icu_vr * reporth_ICU * reporth_g +
        gamma * sigmaER  * out_mat[, (ERindex+1)]  %*% ihr_col * crit * prob_icu_r  * reporth_ICU * reporth_g +
        gamma * out_mat[, (Eindex+1)]              %*% ihr_col * prob_icu    * (1 - reporth_ICU) * reporth_g +
        gamma * out_mat[, (QEindex+1)]             %*% ihr_col * prob_icu    * (1 - reporth_ICU) * reporth_g +
        gamma * sigmaEV  * out_mat[, (EVindex+1)]  %*% ihr_col * prob_icu_v  * (1 - reporth_ICU) * reporth_g +
        gamma * sigmaEVR * out_mat[, (EVRindex+1)] %*% ihr_col * prob_icu_vr * (1 - reporth_ICU) * reporth_g +
        gamma * sigmaER  * out_mat[, (ERindex+1)]  %*% ihr_col * prob_icu_r  * (1 - reporth_ICU) * reporth_g

      rowSums(incidence) + rowSums(incidenceh)
   
    })

  return(incd_time_series)

    # # daily incidence
    # incidence<-parameters_dup["report"]*parameters_dup["gamma"]*(1-parameters_dup["pclin"])*mat_ode[,(Eindex+1)]%*%(1-ihr[,2])+
    #   parameters_dup["reportc"]*parameters_dup["gamma"]*parameters_dup["pclin"]*mat_ode[,(Eindex+1)]%*%(1-ihr[,2])+
    #   parameters_dup["report"]*parameters_dup["gamma"]*(1-parameters_dup["pclin"])*mat_ode[,(QEindex+1)]%*%(1-ihr[,2])+
    #   parameters_dup["reportc"]*parameters_dup["gamma"]*parameters_dup["pclin"]*mat_ode[,(QEindex+1)]%*%(1-ihr[,2])+
    #   parameters_dup["report_v"]*parameters_dup["gamma"]*(1-parameters_dup["pclin_v"])*mat_ode[,(EVindex+1)]%*%(1-parameters_dup["sigmaEV"]*ihr[,2])+
    #   parameters_dup["report_cv"]*parameters_dup["gamma"]*parameters_dup["pclin_v"]*mat_ode[,(EVindex+1)]%*%(1-parameters_dup["sigmaEV"]*ihr[,2])+
    #   parameters_dup["report_vr"]*parameters_dup["gamma"]*(1-parameters_dup["pclin_vr"])*mat_ode[,(EVRindex+1)]%*%(1-parameters_dup["sigmaEVR"]*ihr[,2])+
    #   parameters_dup["report_cvr"]*parameters_dup["gamma"]*parameters_dup["pclin_vr"]*mat_ode[,(EVRindex+1)]%*%(1-parameters_dup["sigmaEVR"]*ihr[,2])+
    #   parameters_dup["report_r"]*parameters_dup["gamma"]*(1-parameters_dup["pclin_r"])*mat_ode[,(ERindex+1)]%*%(1-parameters_dup["sigmaER"]*ihr[,2])+
    #   parameters_dup["report_cr"]*parameters_dup["gamma"]*parameters_dup["pclin_r"]*mat_ode[,(ERindex+1)]%*%(1-parameters_dup["sigmaER"]*ihr[,2])
    
    # incidenceh<- parameters_dup["gamma"]*mat_ode[,(Eindex+1)]%*%ihr[,2]*(1-critH)*(1-parameters_dup["prob_icu"])*parameters_dup["reporth"]+
    #   parameters_dup["gamma"]*mat_ode[,(Eindex+1)]%*%ihr[,2]*(1-critH)*(1-parameters_dup["prob_icu"])*(1-parameters_dup["reporth"])*parameters_dup["reporth_g"]+
    #   parameters_dup["gamma"]*mat_ode[,(QEindex+1)]%*%ihr[,2]*(1-critH)*(1-parameters_dup["prob_icu"])*parameters_dup["reporth"]+
    #   parameters_dup["gamma"]*mat_ode[,(QEindex+1)]%*%ihr[,2]*(1-critH)*(1-parameters_dup["prob_icu"])*(1-parameters_dup["reporth"])*parameters_dup["reporth_g"]+
    #   parameters_dup["gamma"]*parameters_dup["sigmaEV"]*mat_ode[,(EVindex+1)]%*%ihr[,2]*(1-critH)*(1-parameters_dup["prob_icu_v"])*parameters_dup["reporth"]+
    #   parameters_dup["gamma"]*parameters_dup["sigmaEVR"]*mat_ode[,(EVRindex+1)]%*%ihr[,2]*(1-critH)*(1-parameters_dup["prob_icu_vr"])*parameters_dup["reporth"]+
    #   parameters_dup["gamma"]*parameters_dup["sigmaER"]*mat_ode[,(ERindex+1)]%*%ihr[,2]*(1-critH)*(1-parameters_dup["prob_icu_r"])*parameters_dup["reporth"]+
    #   parameters_dup["gamma"]*mat_ode[,(Eindex+1)]%*%ihr[,2]*critH*parameters_dup["reporth_g"]*(1-parameters_dup["prob_icu"])+
    #   parameters_dup["gamma"]*mat_ode[,(QEindex+1)]%*%ihr[,2]*critH*parameters_dup["reporth_g"]*(1-parameters_dup["prob_icu"])+
    #   parameters_dup["gamma"]*parameters_dup["sigmaEV"]*mat_ode[,(EVindex+1)]%*%ihr[,2]*critH*parameters_dup["reporth_g"]*(1-parameters_dup["prob_icu_v"])+
    #   parameters_dup["gamma"]*parameters_dup["sigmaEVR"]*mat_ode[,(EVRindex+1)]%*%ihr[,2]*critH*parameters_dup["reporth_g"]*(1-parameters_dup["prob_icu_vr"])+
    #   parameters_dup["gamma"]*parameters_dup["sigmaER"]*mat_ode[,(ERindex+1)]%*%ihr[,2]*critH*parameters_dup["reporth_g"]*(1-parameters_dup["prob_icu_r"])+
    #   #ICU
    #   parameters_dup["gamma"]*mat_ode[,(Eindex+1)]%*%ihr[,2]*parameters_dup["prob_icu"]*(1-crit)*parameters_dup["reporth_ICU"]+
    #   parameters_dup["gamma"]*mat_ode[,(QEindex+1)]%*%ihr[,2]*parameters_dup["prob_icu"]*(1-crit)*parameters_dup["reporth_ICU"]+
    #   parameters_dup["gamma"]*mat_ode[,(Eindex+1)]%*%ihr[,2]*parameters_dup["prob_icu"]*crit*parameters_dup["reporth_ICU"]*parameters_dup["reporth_g"]+
    #   parameters_dup["gamma"]*mat_ode[,(QEindex+1)]%*%ihr[,2]*parameters_dup["prob_icu"]*crit*parameters_dup["reporth_ICU"]*parameters_dup["reporth_g"]+
    #   parameters_dup["gamma"]*parameters_dup["sigmaEV"]*mat_ode[,(EVindex+1)]%*%ihr[,2]*(1-crit)*parameters_dup["prob_icu_v"]*parameters_dup["reporth_ICU"]+
    #   parameters_dup["gamma"]*parameters_dup["sigmaEVR"]*mat_ode[,(EVRindex+1)]%*%ihr[,2]*(1-crit)*parameters_dup["prob_icu_vr"]*parameters_dup["reporth_ICU"]+
    #   parameters_dup["gamma"]*parameters_dup["sigmaER"]*mat_ode[,(ERindex+1)]%*%ihr[,2]*(1-crit)*parameters_dup["prob_icu_r"]*parameters_dup["reporth_ICU"]+
    #   parameters_dup["gamma"]*parameters_dup["sigmaEV"]*mat_ode[,(EVindex+1)]%*%ihr[,2]*crit*parameters_dup["prob_icu_v"]*parameters_dup["reporth_ICU"]*parameters_dup["reporth_g"]+
    #   parameters_dup["gamma"]*parameters_dup["sigmaEVR"]*mat_ode[,(EVRindex+1)]%*%ihr[,2]*crit*parameters_dup["prob_icu_vr"]*parameters_dup["reporth_ICU"]*parameters_dup["reporth_g"]+
    #   parameters_dup["gamma"]*parameters_dup["sigmaER"]*mat_ode[,(ERindex+1)]%*%ihr[,2]*crit*parameters_dup["prob_icu_r"]*parameters_dup["reporth_ICU"]*parameters_dup["reporth_g"]+
    #   parameters_dup["gamma"]*mat_ode[,(Eindex+1)]%*%ihr[,2]*parameters_dup["prob_icu"]*(1-parameters_dup["reporth_ICU"])*parameters_dup["reporth_g"]+
    #   parameters_dup["gamma"]*mat_ode[,(QEindex+1)]%*%ihr[,2]*parameters_dup["prob_icu"]*(1-parameters_dup["reporth_ICU"])*parameters_dup["reporth_g"]+
    #   parameters_dup["gamma"]*parameters_dup["sigmaEV"]*mat_ode[,(EVindex+1)]%*%ihr[,2]*parameters_dup["prob_icu_v"]*(1-parameters_dup["reporth_ICU"])*parameters_dup["reporth_g"]+
    #   parameters_dup["gamma"]*parameters_dup["sigmaEVR"]*mat_ode[,(EVRindex+1)]%*%ihr[,2]*parameters_dup["prob_icu_vr"]*(1-parameters_dup["reporth_ICU"])*parameters_dup["reporth_g"]+
    #   parameters_dup["gamma"]*parameters_dup["sigmaER"]*mat_ode[,(ERindex+1)]%*%ihr[,2]*parameters_dup["prob_icu_r"]*(1-parameters_dup["reporth_ICU"])*parameters_dup["reporth_g"]
    
}


get_incidents_v16_8 <- function(
    parameters,
    out_mat,
    ihr_col,
    times
  ) {

  spline_funs <- make_spline_funs(
      give = parameters["give"],
      beds_available = parameters["beds_available"],
      icu_beds_available = parameters["icu_beds_available"],
      ventilators_available = parameters["ventilators_available"]
  )

  incd_time_series <-
    with(as.list(parameters, spline_funs), {

      critH<-c()
      crit<-c()
      critV<-c()

      for (t in 1:length(times)){
        critH[t] <- min(1-fH((sum(out_mat[t,(Hindex+1)]))+sum(out_mat[t,(ICUCindex+1)])+sum(out_mat[t,(ICUCVindex+1)])),1)
        crit[t]  <- min(1-fICU((sum(out_mat[t,(ICUindex+1)]))+(sum(out_mat[t,(Ventindex+1)]))+(sum(out_mat[t,(VentCindex+1)]))))
        critV[t] <- min(1-fVent((sum(out_mat[t,(Ventindex+1)]))),1)
      }


      # daily incidence
      incidence <-
        report     * gamma * (1 - pclin)    * out_mat[, (Eindex+1)]   %*% (1 - ihr_col) +
        reportc    * gamma * pclin          * out_mat[, (Eindex+1)]   %*% (1 - ihr_col) +
        report     * gamma * (1 - pclin)    * out_mat[, (QEindex+1)]  %*% (1 - ihr_col) +
        reportc    * gamma * pclin          * out_mat[, (QEindex+1)]  %*% (1 - ihr_col) +
        report_v   * gamma * (1 - pclin_v)  * out_mat[, (EVindex+1)]  %*% (1 - sigmaEV * ihr_col) +
        report_cv  * gamma * pclin_v        * out_mat[, (EVindex+1)]  %*% (1 - sigmaEV * ihr_col) +
        report_vr  * gamma * (1 - pclin_vr) * out_mat[, (EVRindex+1)] %*% (1 - sigmaEVR * ihr_col) +
        report_cvr * gamma * pclin_vr       * out_mat[, (EVRindex+1)] %*% (1 - sigmaEVR * ihr_col) +
        report_r   * gamma * (1 - pclin_r)  * out_mat[, (ERindex+1)]  %*% (1 - sigmaER * ihr_col) +
        report_cr  * gamma * pclin_r        * out_mat[, (ERindex+1)]  %*% (1 - sigmaER * ihr_col)
      
        # print("incidence:")
        # print(incidence)

      incidenceh <-
        gamma * out_mat[, (Eindex+1)]  %*% ihr_col * (1 - critH) * (1 - prob_icu) * reporth +
        gamma * out_mat[, (Eindex+1)]  %*% ihr_col * (1 - critH) * (1 - prob_icu) * (1 - reporth) * reporth_g +
        gamma * out_mat[, (QEindex+1)] %*% ihr_col * (1 - critH) * (1 - prob_icu) * reporth +
        gamma * out_mat[, (QEindex+1)] %*% ihr_col * (1 - critH) * (1 - prob_icu) * (1 - reporth) * reporth_g +
        gamma * sigmaEV  * out_mat[, (EVindex+1)]  %*% ihr_col * (1 - critH) * (1 - prob_icu_v)  * reporth +
        gamma * sigmaEVR * out_mat[, (EVRindex+1)] %*% ihr_col * (1 - critH) * (1 - prob_icu_vr) * reporth +
        gamma * sigmaER  * out_mat[, (ERindex+1)]  %*% ihr_col * (1 - critH) * (1 - prob_icu_r)  * reporth +
        gamma * out_mat[, (Eindex+1)]              %*% ihr_col * critH * reporth_g * (1 - prob_icu) +
        gamma * out_mat[, (QEindex+1)]             %*% ihr_col * critH * reporth_g * (1 - prob_icu) +
        gamma * sigmaEV  * out_mat[, (EVindex+1)]  %*% ihr_col * critH * reporth_g * (1 - prob_icu_v) +
        gamma * sigmaEVR * out_mat[, (EVRindex+1)] %*% ihr_col * critH * reporth_g * (1 - prob_icu_vr) +
        gamma * sigmaER  * out_mat[, (ERindex+1)]  %*% ihr_col * critH * reporth_g * (1 - prob_icu_r) +
        gamma * out_mat[, (Eindex+1)]              %*% ihr_col * prob_icu * (1 - crit) * reporth_ICU +
        gamma * out_mat[, (QEindex+1)]             %*% ihr_col * prob_icu * (1 - crit) * reporth_ICU +
        gamma * out_mat[, (Eindex+1)]              %*% ihr_col * prob_icu * crit * reporth_ICU * reporth_g +
        gamma * out_mat[, (QEindex+1)]             %*% ihr_col * prob_icu * crit * reporth_ICU * reporth_g +
        gamma * sigmaEV  * out_mat[, (EVindex+1)]  %*% ihr_col * (1 - crit) * prob_icu_v * reporth_ICU +
        gamma * sigmaEVR * out_mat[, (EVRindex+1)] %*% ihr_col * (1 - crit) * prob_icu_vr * reporth_ICU +
        gamma * sigmaER  * out_mat[, (ERindex+1)]  %*% ihr_col * (1 - crit) * prob_icu_r * reporth_ICU +
        gamma * sigmaEV  * out_mat[, (EVindex+1)]  %*% ihr_col * crit * prob_icu_v  * reporth_ICU * reporth_g +
        gamma * sigmaEVR * out_mat[, (EVRindex+1)] %*% ihr_col * crit * prob_icu_vr * reporth_ICU * reporth_g +
        gamma * sigmaER  * out_mat[, (ERindex+1)]  %*% ihr_col * crit * prob_icu_r  * reporth_ICU * reporth_g +
        gamma * out_mat[, (Eindex+1)]              %*% ihr_col * prob_icu    * (1 - reporth_ICU) * reporth_g +
        gamma * out_mat[, (QEindex+1)]             %*% ihr_col * prob_icu    * (1 - reporth_ICU) * reporth_g +
        gamma * sigmaEV  * out_mat[, (EVindex+1)]  %*% ihr_col * prob_icu_v  * (1 - reporth_ICU) * reporth_g +
        gamma * sigmaEVR * out_mat[, (EVRindex+1)] %*% ihr_col * prob_icu_vr * (1 - reporth_ICU) * reporth_g +
        gamma * sigmaER  * out_mat[, (ERindex+1)]  %*% ihr_col * prob_icu_r  * (1 - reporth_ICU) * reporth_g

      rowSums(incidence) + rowSums(incidenceh)

      day_infections<- round(rowSums( gamma * out_mat[,(Eindex+1)] +
                                      gamma * out_mat[,(QEindex+1)] +
                                      gamma * out_mat[,(EVindex+1)] +
                                      gamma * out_mat[,(EVRindex+1)] +
                                      gamma * out_mat[,(ERindex+1)]
                        ))

      list("incd" = incidence, "incdh" = incidenceh, "day_inf" = day_infections)
   
    })

  return(incd_time_series)

}