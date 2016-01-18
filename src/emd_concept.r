library(EMD)

#Arg <- function(z) {
#  atan(Im(z) / Re(z))
#}


sinusoid <- function(t, period = 1.0) {

  angular.frequency <- 2*pi / period

  sin(angular.frequency * t)

}


sawtooth <- function(t, period = 1.0) {

  t.over.period <- t / period

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


decompose.signal <- function(t, x) {

  emd(xt = x, tt = t, boundary = "periodic")

}


fourier.transform <- function(t, imfs) {

  ffts <- apply(imfs, 2, fft)
  amplitudes <- apply(ffts, 2, function(x) { Mod(x)^2 })
  phases <- apply(ffts, 2, Arg)

  N <- length(t)
  T <- max(t) - min(t)
  delta.f <- 1/T
  f.max <- (N/2) * delta.f
  freq <- seq(-f.max, f.max, delta.f)[-floor(N/2)]

  list(freq          = freq,
       fft           = ffts,
       amplitudes    = amplitudes,
       phase.spectra = phases,
       nimf          = dim(imfs)[2])

}


extract.determinism <- function(emd.result, fourier.transform,
                                criteria = 2) {

  # placeholder, in the future this will only happen if criteria is an
  # integer, and otherwise will assume it's a function which returns an
  # integer
  index <- criteria

  components <- emd.result$imf[,(index:emd.result$nimf)]

  rowSums(components)

}


fit.statistics <- function(signal, deterministic.signal) {

  N <- length(deterministic.signal)
  RSS <- sum((signal$signal - deterministic.signal)^2)

  MSE <- RSS / N

  data.frame(MSE = MSE)

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


plot.determinism <- function(signal, deterministic.signal) {

  pdf("deterministic_signal.pdf")

  plot(signal$t, signal$signal,
       type = "l", col = "blue")
  lines(signal$t, deterministic.signal,
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


plot.amplitude.spectra <- function(fourier.transform) {

  pdf("imf_amplitude_spectra.pdf")

  par(mfrow = c(ceiling(fourier.transform$nimf / 2), 2))

  N <- length(fourier.transform$freq)
  idx <- seq(floor(N/2), N)

  for(i in 1:fourier.transform$nimf) {

    amplitudes <- fourier.transform$amplitudes[idx,i]
    plot(fourier.transform$freq[idx], amplitudes,
         type = "l")

  }

  dev.off()

}


plot.phase.spectra <- function(fourier.transform) {

  pdf("imf_phase_spectra.pdf")

  par(mfrow = c(ceiling(fourier.transform$nimf / 2), 2))

  for(i in 1:fourier.transform$nimf) {

    phase.spectrum <- fourier.transform$phase.spectra[,i]
    plot(fourier.transform$freq, phase.spectrum,
         type = "l")

  }

  dev.off()

}


plot.imf.partitions <- function(signal, emd.result) {

  pdf("extracted_determinism_candidates.pdf")

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

  signal <- generate.signal(sinusoid,
                            noise.sd = 0.5, num.cycles = 5)

  emd.result <- decompose.signal(signal$t, signal$observation)
  fourier.transform <- fourier.transform(signal$t, emd.result$imf)
  deterministic.signal <- extract.determinism(emd.result, fourier.transform,
                                              criteria = 3)
  stats <- fit.statistics(signal, deterministic.signal)

  plot.signal(signal)
  plot.determinism(signal, deterministic.signal)
  plot.emd(signal, emd.result)
  plot.amplitude.spectra(fourier.transform)
  plot.phase.spectra(fourier.transform)
  plot.imf.partitions(signal, emd.result)

  print(stats, row.names = FALSE)

}


main()
