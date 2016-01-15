library(EMD)

generate.signal <- function(period = 1.0, num.cycles = 2.0, noise.sd = 1.0) {

  angular.frequency <- 2*pi / period

  t <- seq(0, num.cycles, 0.01)

  signal <- cos(angular.frequency*t)
  noise  <- rnorm(signal, sd = noise.sd)

  observation <- signal+noise

  list(t           = t,
       signal      = signal,
       noise       = noise,
       observation = observation)

}


separate.signal <- function(t, x) {

  result <- emd(xt = x, tt = t, boundary = "periodic")

  result

}


plot.signal <- function(signal) {

  pdf("generated_signal.pdf")

  plot(signal$t, signal$observation,
       type = "p", lwd = 1, pch = ".",
       xlab = "t", ylab = "x(t)")
  lines(signal$t, signal$signal,
        col = "red")

  dev.off()

}


plot.emd <- function(signal, emd.result) {

  pdf("emd_results.pdf")

  par(mfrow = c(ceiling(emd.result$nimf/2)+1, 2))

  plot(signal$t, signal$observation,
       type = "p", lwd = 1, pch = ".")

  for(i in 1:emd.result$nimf) {

    imf <- emd.result$imf[,i]

    plot(signal$t, imf,
         type = "l")

  }

  plot(signal$t, emd.result$residue,
       type = "l")

  dev.off()

}


plot.imf.partitions <- function(signal, emd.result) {

  pdf("extracted_determinism_guesses.pdf")

  par(mfrow = c(ceiling(emd.result$nimf/2), 2))

  for(i in 2:emd.result$nimf) {

    imf.partition <- emd.result$imf[,-(1:i-1)]

    d <- emd.result$residue +
         if (!is.null(dim(imf.partition))) {
           rowSums(imf.partition)
         } else {
           imf.partition
         }

    plot(signal$t, signal$signal,
         type = "l", col = "blue")
    lines(signal$t, d,
          col = "red")

  }

  dev.off()

}



main <- function() {

  set.seed(4)

  signal <- generate.signal(noise.sd = 0.5)

  emd.result <- separate.signal(signal$t, signal$observation)

  plot.signal(signal)

  plot.emd(signal, emd.result)

  plot.imf.partitions(signal, emd.result)

}


main()