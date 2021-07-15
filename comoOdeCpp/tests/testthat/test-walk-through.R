test_that("Walk-through (v19.1.2) of the model code", {
  # skip("temp skip")

  rm(list = ls())
  source(paste0(getwd(), "/common.R"), local = environment())
  init(e = environment())

  file_path <- paste0(getwd(), "/data/templates_v19.1.2/Template_CoMoCOVID-19App_v19_all_interventions.xlsx")

  if (!exists("inputs", mode = "function")) {
    source(paste0(getwd(), CORE_FILE), local = environment())
  #   # print("before source")
  #   # source(paste0(getwd(), "/v16.3.core.R"), local = environment())
  #   source(paste0(getwd(), "/v16.4.core.mod.16.6.R"), local = environment())
  #   # print("after source")
  }
  # print("done sourcing")

  # stop()

  check_parameters_list_for_na(parameters)

  # environment(check_mortality_count) <- environment()

  # p_value_list <- seq(0.0, 0.1, by = 0.025)
  p_value_list <- seq(0.1, 0.1)

  scenario_list <- list(
    # vectors0, # Baseline
    vectors0 # Baseline
    # vectors   # Hypothetical
  )

  for (ss in scenario_list) {
    for (pp in p_value_list) {
      parameters["p"] <- pp

      param_vector <- parameters
      # param_vector[parameters_noise] <- parameters[parameters_noise]
      #   + rnorm(
      #       length(parameters_noise),
      #       mean = 0,
      #       sd = noise * abs(parameters[parameters_noise])
      #     )

      RUN_CPP <- TRUE
      RUN_R <- TRUE

      if (RUN_CPP) {
        # start_time <- Sys.time()
        covidOdeCpp_reset()
        output_message <- capture_output(
          out_cpp <- ode(
                      y = Y, times = times, method = "euler", hini = 0.05,
                      func = covidOdeCpp, parms = param_vector,
                      input = ss, A = A,
                      contact_home = contact_home,
                      contact_school = contact_school,
                      contact_work = contact_work,
                      contact_other = contact_other,
                      popbirth_col2 = popbirth[, 2],
                      popstruc_col2 = popstruc[, 2],
                      ageing = ageing,
                      ifr_col2 = ifr[, 2],
                      ihr_col2 = ihr[, 2],
                      mort_col = mort,
                      age_group_vectors = age_group_vectors
                      )

        )
        # elapsed_time <- Sys.time() - start_time

        # print("Rcpp version time:")
        # print(elapsed_time)

        expect_equal(output_message, "covidOdeCpp: splinefuns updated")

        processed_cpp_results <- multi_runs_pp_only(
              out0 = out_cpp, times = times, parameters = param_vector, input = ss
            )

        processed_cpp_results <- process_ode_outcome(
              out = processed_cpp_results, intv_vector = ss, param_used = param_vector
            )

        # processed_cpp_results <- process_ode_outcome(out_cpp, ss, param_vector)
        print("total_reported_deaths_end (Cpp):")
        print(processed_cpp_results$total_reported_deaths_end)

      }

      if (RUN_R) {
        # start_time <- Sys.time()
        # print("age_group_vectors:")
        # print(age_group_vectors)

        expect_silent(
          out_r <- ode(
                    y = Y, times = times, method = "euler", hini = 0.05,
                    func = covid, parms = param_vector, input = ss
                    )
        )

        # print("dim(out_r):")
        # print(dim(out_r))
        # elapsed_time <- Sys.time() - start_time

        # print("R version time:")
        # print(elapsed_time)

        processed_r_results <- multi_runs_pp_only(
              out0 = out_r, times = times, parameters = param_vector, input = ss
            )

        processed_r_results <- process_ode_outcome(
              out = processed_r_results, intv_vector = ss, param_used = param_vector
            )

        # processed_r_results <- process_ode_outcome(out_r, ss, param_vector)
        print("total_reported_deaths_end (R):")
        print(processed_r_results$total_reported_deaths_end)


      }

      if (RUN_R && RUN_CPP) {
        match_processed_outputs(
            output_a = processed_cpp_results,
            output_b = processed_r_results,
            tlr = 0.0001
        )

        # # sss = 1
        # # write.csv(out_cpp, paste0("out_cpp_",sss,"_",parameters["p"],".csv"),row.names = FALSE)
        # # write.csv(out_r, paste0("out_r_",sss,"_",parameters["p"],".csv"),row.names = FALSE)

        match_outputs(
          output_a = out_r,
          output_b = out_cpp,
          tlr = 0.0001,
          smp = 10000
        )
      }


    }
  }

})