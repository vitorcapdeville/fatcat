plot_chain <- function(chain, true_value = NULL, ylab = "") {
  ggplot2::ggplot(mapping = ggplot2::aes(x = seq(length(chain)))) +
    ggplot2::xlab("") +
    ggplot2::ylab(ylab) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(ggplot2::aes(y = chain)) +
    if (!is.null(true_value) & !is.na(true_value)) {
      ggplot2::geom_line(ggplot2::aes(y = true_value), color = "red")
    } else {
      NULL
    }
}
#' Função para plotar os resultados do ajuste
#'
#' @param res lista contendo os elemntos beta, f e sigma2, assim como exportado pela
#' função `fitfatcat`
#' @param true_values argumento opcional com os valores verdadeiros dos parametros,
#' para o caso de estudos de simulação.
#' @param name qual parametro deve ser plotado, beta, f ou sigma2?
#'
#' @return nada
#'
#' @export
#'
plotfatcat <- function(res, true_values = NULL, name = c("beta", "f", "sigma2")) {
  name <- match.arg(name)
  res <- res[[name]]
  if (name == "f" & dim(res)[2] > 5) {
    res <- res[, 1:5, ]
    true_values <- true_values[, 1:5]
    message("Plotting only res$f[,1:5,].")
  } else if (name == "beta") {
    if (dim(res)[2] > dim(true_values)[2]) {
      message("Matching 'true values' with the fisrt p columns of beta, since true_p > p.")
    }
  }
  true_dim <- dim(res)[-length(dim(res))]
  indices <- seq(prod(true_dim))
  size_total <- prod(dim(res))
  size_it <- prod(true_dim)
  patchwork::wrap_plots(
    lapply(
      indices,
      function(x) {
        col_ind <- floor((x + 1) / dim(res)[1])
        row_ind <- x - dim(res)[1] * (col_ind - 1)
        plot_chain(
          res[seq(from = x, to = size_total, by = size_it)],
          true_values[x],
          ylab = glue::glue("{name}[{row_ind},{col_ind}]")
        )
      }
    ),
    nrow = true_dim[1], byrow = F
  )
}
