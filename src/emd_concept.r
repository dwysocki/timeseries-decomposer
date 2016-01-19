#!/usr/bin/Rscript

suppressWarnings(suppressMessages(library(argparse)))
suppressWarnings(suppressMessages(library(EMD)))


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

signal.types <- list(
    sinusoid = sinusoid,
    sawtooth = sawtooth,
    triangle = triangle,
    square   = square
)

# placeholder
noise.types <- list(
    normal   = rnorm
)


generate.signal <- function(fn = sinusoid,
                            period = 1.0, num.cycles = 2.0, precision = 0.01,
                            noise.sd = 1.0) {

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
    index <- min(criteria, emd.result$nimf)

    if(emd.result$nimf > 1) {
        components <- emd.result$imf[,(index:emd.result$nimf)]

        if(!is.null(dim(components))) {
            rowSums(components)
        } else {
            components
        }
    } else {
        emd.result$imf
    }

}


fit.statistics <- function(signal, deterministic.signal) {

    N <- length(deterministic.signal)
    RSS <- sum((signal$signal - deterministic.signal)^2)

    MSE <- RSS / N

    data.frame(MSE = MSE)

}


plot.signal <- function(signal,
                        output = ".") {

    pdf(file.path(output, "generated_signal.pdf"))

    plot(signal$t, signal$observation,
         type = "p", lwd = 1, pch = ".",
         xlab = "t", ylab = "x(t)")
    lines(signal$t, signal$signal,
          col = "red")

    dev.off()

}


plot.determinism <- function(signal, deterministic.signal,
                             output = ".") {

    pdf(file.path(output, "deterministic_signal.pdf"))

    plot(signal$t, signal$signal,
         type = "l", col = "blue")
    lines(signal$t, deterministic.signal,
          col = "red")

    dev.off()

}


plot.emd <- function(signal, emd.result,
                     output = ".") {

    pdf(file.path(output, "emd_results.pdf"))

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


plot.amplitude.spectra <- function(fourier.transform,
                                   output = ".") {

    pdf(file.path(output, "imf_amplitude_spectra.pdf"))

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


plot.phase.spectra <- function(fourier.transform,
                               output = ".") {

    pdf(file.path(output, "imf_phase_spectra.pdf"))

    par(mfrow = c(ceiling(fourier.transform$nimf / 2), 2))

    for(i in 1:fourier.transform$nimf) {

        phase.spectrum <- fourier.transform$phase.spectra[,i]
        plot(fourier.transform$freq, phase.spectrum,
             type = "l")

    }

    dev.off()

}


plot.imf.partitions <- function(signal, emd.result,
                                output = ".") {

    pdf(file.path(output, "extracted_determinism_candidates.pdf"))

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


get.args <- function() {

    p <- ArgumentParser(
        description = "Concept script for decomposing singals with EMD.")

    add_argument <- function(...) { p$add_argument(...) }

    add_argument("signal_type", choices = names(signal.types),
        help = "Type of signal to use as input.")
    add_argument("noise_type", choices = names(noise.types),
        help = "Type of noise to add to signal.")
    add_argument("criteria", type = "integer",
        help = "Criteria for separating stochastic and deterministic signals.")

    add_argument("-o", "--output",
        help = "Output directory for plots.")

    add_argument("--noise-sd", type = "double", default = 1.0,
        help = "Standard deviation of added noise.")
    add_argument("--num-cycles", type = "double", default = 2.0,
        help = "Number of times to repeat generated signal.")
    add_argument("--period", type = "double", default = 1.0,
        help = "Period of generated signal.")
    add_argument("--time-step", type = "double", default = 0.01,
        help = "Time between samples in generated signal.")
    add_argument("--seed", type = "integer", default = 0,
        help = "Random number generator seed.")

    args <- p$parse_args()

    args$signal_type <- signal.types[[args$signal_type]]
    args$noise_type <- noise.types[[args$noise_type]]

    args

}


main <- function() {

    args <- get.args()

    set.seed(args$seed)

    if(!is.null(args$output)) {
        dir.create(args$output, recursive = TRUE)
    }

    signal <- generate.signal(args$signal_type,
                              period = args$period,
                              noise.sd = args$noise_sd,
                              precision = args$time_step,
                              num.cycles = args$num_cycles)

    emd.result <- decompose.signal(signal$t, signal$observation)
    fourier.transform <- fourier.transform(signal$t, emd.result$imf)
    deterministic.signal <- extract.determinism(emd.result, fourier.transform,
                                                criteria = args$criteria)
    stats <- fit.statistics(signal, deterministic.signal)

    if(!is.null(args$output)) {
        plot.signal(signal, output = args$output)
        plot.determinism(signal, deterministic.signal, output = args$output)
        plot.emd(signal, emd.result, output = args$output)
        plot.amplitude.spectra(fourier.transform, output = args$output)
        plot.phase.spectra(fourier.transform, output = args$output)
        plot.imf.partitions(signal, emd.result, output = args$output)
    }

    print(stats, row.names = FALSE)

}


main()
