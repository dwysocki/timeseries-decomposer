library(EMD)

#Arg <- function(z) {
#  atan(Im(z) / Re(z))
#}


sinusoid <- function(t, period = 1.0) {

  angular.frequency <- 2*pi / period

  sin(angular.frequency * t)

}


sawtooth <- function(t, period = 1.0) {

  t.over.period = t / period

  2 * (t.over.period - floor(0.5 + t.over.period))

}


triangle <- function(t, period = 1.0) {

  2 * abs(sawtooth(t, period = period)) - 1

}


square <- function(t, period = 1.0) {

  sign(sinusoid(t, period = period))

}


generate.signal <- function(fn = sinusoid,
                            period = 1.0, num.cycles = 2.0, precision = 0.01,
                            noise.sd = 1.0) {

  angular.frequency <- 2*pi / period

  t <- seq(0, num.cycles*period, precision)

  signal <- fn(t, period = period)
  noise  <- rnorm(signal, sd = noise.sd)

  observation <- signal+noise

  list(t           = t,
       signal      = signal,
       noise       = noise,
       observation = observation)

}


separate.signal <- function(t, x) {

  result <- emd(xt = x, tt = t, boundary = "periodic")

}


phase.spectra <- function(imfs) {

  ffts <- apply(imfs, 2, fft)
  phases <- apply(ffts, 2, Arg)

  list(fft           = ffts,
       phase.spectra = phases,
       nimf          = dim(imfs)[2])

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


plot.phase.spectra <- function(signal, phase.spectra) {

  pdf("imf_phase_spectra.pdf")

  par(mfrow = c(ceiling(phase.spectra$nimf / 2), 2))

  for(i in 1:phase.spectra$nimf) {

    phase.spectrum <- phase.spectra$phase.spectra[,i]
    print(phase.spectra$phase.spectra)
    plot(signal$t, phase.spectrum,
         type = "l")

  }

  dev.off()

}


plot.imf.partitions <- function(signal, emd.result) {

  pdf("extracted_determinism_guesses.pdf")

  par(mfrow = c(ceiling(emd.result$nimf/2), 2))

  for(i in 1:emd.result$nimf) {

    imf.partition <- emd.result$imf[,(i:emd.result$nimf)]

    d <- if (!is.null(dim(imf.partition))) {
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

  signal <- generate.signal(square,
                            noise.sd = 0.5, num.cycles = 5)

  emd.result <- separate.signal(signal$t, signal$observation)
  phase.spectra <- phase.spectra(emd.result$imf)

  plot.signal(signal)
  plot.emd(signal, emd.result)
  plot.phase.spectra(signal, phase.spectra)
  plot.imf.partitions(signal, emd.result)

}


main()