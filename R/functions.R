# Finding the optimal variable cutoff.

# Variable cutoff ------

#' Test survival differences at the variable cutoff.
#'
#' @description Tests for the differences in survival for observations
#' stratified by the given variable cutoff. The testing ia accomplished by
#' the \code{\link[survival]{survdiff}} function.
#' @return a data frame with the cutoff value, chi-squared test statistic,
#' degrees of freedom, p value and n numbers in the low and high strata.
#' @param time name of the survival time variables, needs to be numeric.
#' @param event name of the binary event variable: 1 is event, 0 no event.
#' @param variable name of the variable to be cut, needs to be numeric.
#' @param cutoff the cutoff value.
#' @param ... additional arguments passed to \code{\link[survival]{survdiff}}.
#' @export

  test_cutoff <- function(data, time, event, variable, cutoff, ...) {

    ## entry control.

    stopifnot(is.data.frame(data))
    stopifnot(all(c(time, event, variable) %in% names(data)))
    stopifnot(all(is.numeric(c(data[[time]],
                               data[[event]],
                               data[[variable]],
                               cutoff))))


    ## testing

    data <- dplyr::mutate(data[c(time, event, variable)],
                          strata = cut(.data[[variable]],
                                       c(-Inf, cutoff, Inf),
                                       c('low', 'high')))

    data <- dplyr::filter(data, complete.cases(data))

    surv_obj <- survival::Surv(time = data[[time]],
                               event = data[[event]])

    surv_test <- survival::survdiff(surv_obj ~ strata,
                                    data = data, ...)

    tibble::tibble(cutoff = cutoff,
                   n_low = surv_test$n[1],
                   n_high = surv_test$n[2],
                   chisq = surv_test$chisq,
                   df = 1,
                   p_value = pchisq(q = surv_test$chisq,
                                    df = 1,
                                    lower.tail = FALSE))

  }

# Cutoff optimization ------

#' Find optimal survival cutoff.
#'
#' @description Finds the optimal cutoff of a variable corresponding to the
#' greatest difference in survival. The difference in survival is estimated with
#' \code{\link[survival]{survdiff}}.
#' @details In case, more minuma are find, the stratified data table is
#' adjusted accordingly.
#' @return an object of the 'survcut'. Contains te optimal cutoff value
#' ('cutoff'), a tibble with the numeric variable of interest and and it's
#' strata and a tibble with the cutoff values, test statistics
#' and p values for the survival differences.
#' @param data a data frame with the survival time, a numeric event index
#' variable and the variable to be cut.
#' @param time name of the survival time variables, needs to be numeric.
#' @param event name of the binary event variable: 1 is event, 0 no event.
#' @param variable name of the variable to be cut, needs to be numeric.
#' @param min_n the minimal strata size.
#' @param .parallel logical, should the search be run in parallel?
#' @param ... additional arguments passed to \code{\link[survival]{survdiff}}.
#' @export

  find_cutoff <- function(data,
                          time,
                          event,
                          variable,
                          min_n = 2,
                          .parallel = FALSE, ...) {

    ## entry control

    if(!is.data.frame(data)) {

      stop('A data frame required as the data argument.', call. = FALSE)

    }

    if(any(!c(time, event, variable) %in% names(data))) {

      stop('One of time, event, variable is missing from the data.',
           call. = FALSE)

    }

    data <- data[c(time, event, variable)]

    classes <- purrr::map_lgl(data, is.numeric)

    if(any(!classes)) {

      stop('The variables time, event and variable need to be numeric.',
           call. = FALSE)

    }

    if(any(!data[[event]] %in% c(0, 1))) {

      stop('The event variable needs to be in a 0,1 format.', call. = FALSE)

    }

    min_n <- as.integer(min_n)

    if(min_n > nrow(data)) {

      stop('The min_n value must be less than the observation number.',
           call. = FALSE)

    }

    ## benchmarking

    start_time <- Sys.time()
    message(paste('Finding the survival cutoff, n =', nrow(data)))
    on.exit(message(paste('Elapsed: ', Sys.time() - start_time)))

    ## unique numeric values, eliminating the min/max -> one group

    uni_vals <- sort(unique(data[[variable]]))

    uni_vals <- uni_vals[2:(length(uni_vals) - 1)]

    ## serial testing

    if(.parallel) {

      future::plan('multisession')

      test_res <- furrr::future_map_dfr(uni_vals,
                                        kmOptimizer::test_cutoff,
                                        data = test_data,
                                        time = 'surv_time',
                                        event = 'event_var',
                                        variable = 'test_var', ...,
                                        .options = furrr::furrr_options(seed = TRUE,
                                                                        packages = c('dplyr',
                                                                                     'furrr',
                                                                                     'survival')))

      future::plan('sequential')

    } else {

      test_res <- purrr::map_dfr(uni_vals,
                                 kmOptimizer::test_cutoff,
                                 data = test_data,
                                 time = 'surv_time',
                                 event = 'event_var',
                                 variable = 'test_var', ...)


    }

    ## formatting

    test_res <- dplyr::filter(test_res,
                              n_low >= min_n,
                              n_high >= min_n)

    cutoff_stats <- dplyr::filter(test_res,
                                  p_value == min(p_value))

    cutoff <- cutoff_stats[['cutoff']]

    if(length(cutoff == 1)) {

      new_name <- paste0(variable, '_strata')

    } else {

      new_name <- paste0(variable, '_strata_', 1:length(cutoff))

    }

    ## stratification

    for(i in new_name) {

      data <- dplyr::mutate(data,
                            .observation = 1:nrow(data),
                            !!i := cut(.data[[variable]],
                                       c(-Inf, cutoff, Inf),
                                       c('low', 'high')))

    }

    ## output

    kmOptimizer::survcut(list(cutoff = cutoff,
                              cutoff_stats = cutoff_stats,
                              data = tibble::as_tibble(data[c('.observation',
                                                              time,
                                                              event,
                                                              variable,
                                                              new_name)]),
                              test = tibble::as_tibble(test_res)))


  }


# Survival curves ------

#' Survival curves for a survcut object.
#'
#' @description Calculates survival curves for the variable stratified by
#' the optimal cutoffs.
#' @details Technically, a wrapper around \code{\link[survival]{survfit}}.
#' If more cutoffs were identified, returns a list with the survdfit objects.
#' @return a survfit object or a list of survfit objects.
#' @param x a survcut object.
#' @param ... extra arguments passed to \code{\link[survival]{survfit}}.
#' @export

  surv_fit <- function(x, ...) {

    stopifnot(kmOptimizer::is_survcut(x))

    forms <- formula(x)

    if(!is.list(forms)) {

      return(survival::survfit(forms, data = x$data, ...))

    } else {

      purrr::map(forms,
                 ~survival::survfit(.x, data = x$data, ...))

    }

  }

#' Test for survival differences between the strata.
#'
#' @description Tests for differences in survival for the variable stratified by
#' the optimal cutoffs.
#' @details Technically, a wrapper around \code{\link[survival]{survdiff}}.
#' If more cutoffs were identified, returns a list with the survdiff objects.
#' @return a survdiff object or a list of survdiff objects.
#' @param x a survcut object.
#' @param ... extra arguments passed to \code{\link[survival]{survdiff}}.
#' @export

  surv_diff <- function(x, ...) {

    stopifnot(kmOptimizer::is_survcut(x))

    forms <- formula(x)

    if(!is.list(forms)) {

      return(survival::survdiff(forms, data = x$data, ...))

    } else {

      purrr::map(forms,
                 ~survival::survdiff(.x, data = x$data, ...))

    }

  }
