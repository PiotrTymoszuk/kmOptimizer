# S3 OOP interface.

# Constructor ------

#' Generate a survcut object.
#'
#' @description Creates a survcut object on a top of a named list with the
#' 'cutoff, 'cutoff_stats', 'data' and 'test' components.
#' @param x an object.
#' @param ... extra arguments, currently none specified.
#' @export

  survcut <- function(x, ...) {

    if(!is.list(x)) {

      stop('A list is required.', call. = FALSE)

    }

    if(any(!c('cutoff', 'cutoff_stats', 'data', 'test') %in% names(x))) {

      stop('A list with cutoff, cutoff_stats, data and test elements is required.',
           call. = FALSE)

    }

    if(any(c(!is.data.frame(x$cutoff_stats),
             !is.data.frame(x$data),
             !is.data.frame(x$test)))) {

      stop('The elements cutoff_stats, data and test need to be data frames.',
           call. = FALSE)

    }

    structure(x, class = 'survcut')


  }

# Generics ------

#' Extract the optimal cutoff.
#'
#' @description Extracts the optimal cutoff from an object.
#' @param x an object.
#' @param ... extra arguments passed to the methods.
#' @export

  cutoff <- function(x, ...) {

    UseMethod('cutoff', x)

  }

#' Create survival curves.
#'
#' @description A default method, creates survival curves from a formula or
#' an object with the survival modeling results.
#' @param formula either a formula or a previously fitted model.
#' @param ... extra arguments passed to the methods.
#' @return An object of class survfit containing one or more survival curves.
#' @export

  survfit.default <- function(formula, ...) {

    survival::survfit(formula = formula, ...)

  }

# Class checking ----

#' Checks for the survcut class.
#'
#' @description Checks if the object is an instance of the survcut class.
#' @param x an object.
#' @param ... extra arguments, currently none.
#' @return a logical value.
#' @export

  is_survcut <- function(x, ...) {

    all(class(x) == 'survcut')

  }

# OOP: appearance and extraction -----

#' Appearance of the survcut object.
#'
#' @description Prints a survcut object.
#' @param x a survcut object.
#' @param ... extra arguments, currently none.
#' @return none, called for it's side effects.
#' @export

  print.survcut <- function(x, ...) {

    stopifnot(kmOptimizer::is_survcut(x))

    cat('Survcut object, optimal cutoff at:',
        paste(signif(x$cutoff, 2), collapse = ', '))

  }

#' The optimal cutoff.
#'
#' @description Extracts the optimal cutoff from a survcut object.
#' @param x a survcut object.
#' @param ... extra arguments, currently none.
#' @return a numeric vector with the optima cutoffs.
#' @export cutoff.survcut
#' @export

  cutoff.survcut <- function(x, ...) {

    stopifnot(kmOptimizer::is_survcut(x))

    x$cutoff

  }

#' Stratified data set.
#'
#' @description Extracts the data frame with the stratified data set from
#' a survcut object.
#' @param formula a survcut object.
#' @param ... extra arguments, currently none.
#' @return a data frame with the stratified variables.
#' @export

  model.frame.survcut <- function(formula, ...) {

    stopifnot(kmOptimizer::is_survcut(formula))

    formula$data

  }

#' Testing summary of a survcut object.
#'
#' @description Returns the testing result summary for a survcut object.
#' @param object a survcut object.
#' @param type type of the returned summary: 'all' returns the entire results,
#' 'cutoff' (default) returns the stats for the optimal cutoffs.
#' @param ... extra arguments, currently none.
#' @return a data frame with the testing stats.
#' @export summary.survcut
#' @export

  summary.survcut <- function(object, type = c('cutoff', 'all'), ...) {

    stopifnot(kmOptimizer::is_survcut(object))

    type <- match.arg(type[1], c('cutoff', 'all'))

    switch(type,
           cutoff = object$cutoff_stats,
           all = object$test)

  }

#' Formula of a survcut object.
#'
#' @description Extracts the formula of a survcut object, which may be used
#' survival difference testing.
#' @param x a survcut object.
#' @param ... extra arguments, currently none.
#' @return a formula or a list of formulas, if more optimal, cutoffs were found.
#' @export

  formula.survcut <- function(x, ...) {

    stopifnot(kmOptimizer::is_survcut(x))

    av_vars <- names(x$data)

    time_var <- av_vars[2]

    event_var <- av_vars[3]

    strata_vars <- av_vars[5:length(av_vars)]

    surv_chunk <- paste0('survival::Surv(',
                         time_var,
                         ', ',
                         event_var, ')')

    forms <- paste(surv_chunk, strata_vars, sep = ' ~ ')

    forms <- purrr::map(forms, as.formula)

    if(length(forms) == 1) {

      return(forms[[1]])

    } else {

      return(forms)

    }

  }

# OOP plotting ----

#' Plot a survcut object.
#'
#' @description Plots Kaplan-Meier curves for the variable strata or a series
#' of diagnostic plots visualizing the optimization process.
#' @details For Kaplan-Meier plots, a wrapper around
#' \code{\link[survminer]{ggsurvplot}}. Cutoff and p values are presented in the
#' plot subtitle, complete observation counts are shown in the plot tag.
#' If more cutoffs were identified, returns a list with the Kaplan-Meier plots.
#' @return a ggplot or a list of ggplots.
#' @param x a survcut object.
#' @param type plot type: 'km' plots a Kaplan-Meier plot, 'diagnostic' generates
#' optimization plots.
#' @param palette a color palette used in the Kaplan-Meier plot.
#' @param ... extra arguments passed to \code{\link[survminer]{ggsurvplot}}.
#' @export plot.survcut
#' @export

  plot.survcut <- function(x,
                           type = c('km', 'diagnostic'),
                           palette = c('steelblue', 'firebrick'),
                           ggtheme = survminer::theme_survminer(), ...) {

    ## entry control

    stopifnot(kmOptimizer::is_survcut(x))

    type <- match.arg(type[1], c('km', 'diagnostic'))

    stopifnot(any(class(ggtheme) == 'theme'))

    ## plot meta information

    cutoff_stats <- dplyr::mutate(x$cutoff_stats,
                                  plot_tag = paste0('low: n = ', n_low,
                                                    ', high: n = ', n_high),
                                  plot_subtitle = paste0('cutoff = ', signif(cutoff, 2)),
                                  plot_subtitle = ifelse(p_value < 0.05,
                                                         paste0(plot_subtitle,
                                                                ', p = ',
                                                                signif(p_value, 2)),
                                                         paste0(plot_subtitle,
                                                                ', ns(p = ',
                                                                signif(p_value, 2),
                                                                ')')))

    ## plotting

    if(type == 'km') {

      forms <- formula(x)

      test_var <- names(x$data)[4]

      if(!is.list(forms)) {

        plots <- survminer::ggsurvplot(fit = survminer::surv_fit(forms,
                                                                 data = x$data),
                                       data = x$data,
                                       palette = palette,
                                       legend.title = test_var,
                                       legend.labs = c('low', 'high'), ...)

        plots$plot <- plots$plot +
          ggplot2::labs(subtitle = cutoff_stats$plot_subtitle[1],
                        tag = cutoff_stats$plot_tag[1]) +
          ggplot2::theme(plot.tag.position = 'bottom')

      } else {

        plots <- purrr::map(forms,
                            ~survminer::ggsurvplot(fit = survminer::surv_fit(forms,
                                                                             data = x$data),
                                                   data = x$data,
                                                   palette = palette,
                                                   legend.title = test_var,
                                                   legend.labs = c('low', 'high'), ...))

        plots <- purrr::transpose(plots)

        plots$plot <- purrr::pmap(list(x = plots$plot,
                                       y = cutoff_stats$plot_subtitle,
                                       z = cutoff_stats$plot_tag),
                                  function(x, y, z) x +
                                    ggplot2::labs(subtitle = y,
                                                  tag = z) +
                                    ggplot2::theme(plot.tag.position = 'bottom'))

      }

      return(plots)

    } else {

      plot_tbl <- dplyr::mutate(x$test,
                                hi_ratio = n_high/(n_high + n_low),
                                neg_log_p = -log10(p_value),
                                group = 'group1')

      plots <- purrr::pmap(list(x = c('hi_ratio', 'chisq', 'neg_log_p'),
                                y = list('fraction of high', 'Chi-squared', expression('-log'[10]*' p')),
                                z = c('High fraction', 'Testing results', 'Significance')),
                           function(x, y, z, int) ggplot2::ggplot(plot_tbl,
                                                                  ggplot2::aes(x = .data[['cutoff']],
                                                                               y = .data[[x]])) +
                             ggplot2::geom_line(ggplot2::aes(group = .data[['group']]),
                                                color = 'steelblue') +
                             ggtheme +
                             ggplot2::labs(y = y,
                                           title = z))

      for(i in cutoff_stats$cutoff) {

        plots <- purrr::map(plots,
                            ~.x +
                              ggplot2::geom_vline(xintercept = i,
                                                  linetype = 'dashed',
                                                  color = 'coral3'))

      }

      plots <- rlang::set_names(plots, c('fraction', 'testing', 'significance'))

      return(plots)

    }

  }

# END ----
