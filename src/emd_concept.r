#!/usr/bin/Rscript

suppressWarnings(suppressMessages(library(argparse)))
suppressWarnings(suppressMessages(library(EMD)))
suppressWarnings(suppressMessages(library(c3net)))
suppressWarnings(suppressMessages(library(cross.spectral.analysis)))

## Flat signal, x(t) = 0.
flatline <- function(t, ...) {

    rep(0, length(t))

}

## 2nd order Fourier series
twosine <- function(t, period = 1.0) {

    angular.frequency <- 2*pi / period

    sin(angular.frequency * t) + sin(2*angular.frequency * t + pi/6)

}


## Sinusoid with frequency changing linearly at given rate
freqchange <- function(t, period = 1.0, rate = 0.01) {

    angular.frequency.init <- 2*pi / period
    angular.frequency <- rate*t + angular.frequency.init

    sin(angular.frequency * t)

}

## Sinusoidal signal with given period.
sinusoid <- function(t, period = 1.0) {

    angular.frequency <- 2*pi / period

    sin(angular.frequency * t)

}

## Sawtooth signal with given period.
sawtooth <- function(t, period = 1.0) {

    t.over.period <- t / period

    2 * (t.over.period - floor(0.5 + t.over.period))

}

## Triangle-wave signal with given period.
triangle <- function(t, period = 1.0) {

    2 * abs(sawtooth(t, period = period)) - 1

}

## Square-wave signal with given period.
square <- function(t, period = 1.0) {

    sign(sinusoid(t, period = period))

}

## Mapping of signal type names to generating functions.
signal.types <- list(
    flatline = flatline,
    sinusoid = sinusoid,
    sawtooth = sawtooth,
    triangle = triangle,
    square   = square,
    twosine  = twosine,
    freqchange = freqchange
)

## Mapping of noise type names to generating functions.
noise.types <- list(
    normal   = rnorm,
    foobar   = rnorm
)

## Mapping of EMD type names to functions.
emd.methods <- list(
    emd      = emd,
    semd     = semd
)

## Generate asignal with the given deterministic and stochastic generating
## functions.
generate.signal <- function(fn = sinusoid, noise = rnorm,
                            period = 1.0, num.cycles = 2.0, precision = 0.01,
                            noise.sd = 1.0) {

    t <- seq(0, num.cycles*period, precision)

    signal <- fn(t, period = period)
    noise  <- noise(n = signal, mean = 0.0, sd = noise.sd)

    observation <- signal+noise

    list(t           = t,
         signal      = signal,
         noise       = noise,
         observation = observation)

}

## Decompose a signal using EMD.
decompose.signal <- function(t, x, emd.fn = emd, emd.args = list()) {

    do.call(function(...) { emd.fn(xt = x, tt = t, ...) }, emd.args)

}

## Transform a signal into frequency-space using a Fourier transform or an
## approximation thereof.
## TODO: Add support for different approximations (e.g. periodograms).
fourier.transform <- function(t, imfs,
                              method = "FFT") {

    ## TODO: Make choice of frequency sampling an option, when using LSP.
    N <- length(t)
    T <- max(t) - min(t)
    delta.f <- 1/T
    f.max <- (N/2) * delta.f
    freq <- seq(-f.max, f.max, delta.f)[-floor(N/2)]

    if (method == "FFT") {

        ffts <- apply(imfs, 2, fft)
        amplitudes <- apply(ffts, 2, function(x) { Mod(x)^2 })
        phases <- apply(ffts, 2, Arg)


    } else if (method == "LSP") {

        omega <- 2*pi*freq
        amplitudes <- apply(
            imfs,
            2,
            function (imf) {
                amplitude.spectrum(t, imf, omega)
            }
        )
        phases <- apply(
            imfs,
            2,
            function (imf) {
                phase.spectrum(t, imf, omega)
            }
        )
    } else {

        stop("Error: `method` must be one of { FFT, LSP }")

    }


    list(freq          = freq,
#         fft           = ffts,
         amplitudes    = amplitudes,
         phase.spectra = phases,
         nimf          = dim(imfs)[2])

}

## Extract the deterministic component of an observed signal.
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

## Compute various fit statistics on extracted deterministic signal and
## real signal.
fit.statistics <- function(signal, deterministic.signal) {

    N <- length(deterministic.signal)
    RSS <- sum((signal$signal - deterministic.signal)^2)

    MSE <- RSS / N

    data.frame(MSE = MSE)

}

## Plot the original signal.
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

## Plot the true and extracted deterministic components of the signal.
plot.determinism <- function(signal, deterministic.signal,
                             output = ".") {

    pdf(file.path(output, "deterministic_signal.pdf"))

    plot(signal$t, signal$signal,
         type = "l", col = "blue",
         xlab = "t", ylab = "s(t)")
    lines(signal$t, deterministic.signal,
          col = "red")

    dev.off()

}

## Plot the original signal, as well as its IMFs and residue.
plot.emd <- function(signal, emd.result,
                     output = ".") {

    pdf(file.path(output, "emd_results.pdf"))

    par(mfrow = c(ceiling(emd.result$nimf/2)+1, 2))

    plot(signal$t, signal$observation,
         type = "p", lwd = 1, pch = ".",
         xlab = "t", ylab = "x(t)")

    for(i in 1:emd.result$nimf) {

        imf <- emd.result$imf[,i]

        plot(signal$t, imf,
             type = "l",
             xlab = "t", ylab = substitute(h[i](t), list(i=i)))

    }

    plot(signal$t, emd.result$residue,
         type = "l",
         xlab = "t", ylab = "r(t)")

    dev.off()

}

## Plot the amplitude-components of the Fourier spectra of each IMF.
plot.amplitude.spectra <- function(fourier.transform,
                                   output = ".") {

    pdf(file.path(output, "imf_amplitude_spectra.pdf"))

    par(mfrow = c(ceiling(fourier.transform$nimf / 2), 2))

    N <- length(fourier.transform$freq)
    idx <- seq(floor(N/2), N)

    for(i in 1:fourier.transform$nimf) {

        amplitudes <- fourier.transform$amplitudes[idx,i]
        plot(fourier.transform$freq[idx], amplitudes,
             type = "l",
             xlab = "f", ylab = substitute(A[i](f), list(i=i)))

    }

    dev.off()

}

## Plot the phase-components of the Fourier spectra of each IMF.
plot.phase.spectra <- function(fourier.transform,
                               output = ".") {

    pdf(file.path(output, "imf_phase_spectra.pdf"))

    par(mfrow = c(ceiling(fourier.transform$nimf / 2), 2))

    for(i in 1:fourier.transform$nimf) {

        phase.spectrum <- fourier.transform$phase.spectra[,i]
        plot(fourier.transform$freq, phase.spectrum,
             type = "l",
             xlab = "f", ylab = substitute(theta[i](f), list(i=i)))

    }

    dev.off()

}

## Plot all of the candidate deterministic signals created by
## partitioning the IMFs.
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
             type = "l", col = "blue",
             xlab = "t", ylab = substitute(tilde(s[i](t)), list(i=i)))
        lines(signal$t, d,
              col = "red")

    }

    dev.off()

}

## Plot the mutual information between successive IMF phase spectra.
plot.mutual.information <- function(mutual.information,
                                    output = ".") {

    pdf(file.path(output, "mutual_information.pdf"))

    plot(1:length(mutual.information), mutual.information,
         type = "l", col = "black",
         xlab = "i", ylab = expression(nu[i]))

    dev.off()

}

## Turn a list of key-value pairs [[k1, v1], ..., [kN, vN]] into a list where
## each key maps to its respective value, [k1 = v1, ..., kN = vN].
key.val.list <- function(kv.list) {

    ret = list()

    for (key.val.pair in kv.list) {
        key <- key.val.pair[1]
        val <- key.val.pair[2]

        ret[[key]] <- val
    }

    ret

}

## Evaluate a string containing a valid R expression
parse.and.eval <- function(string) {

    eval(parse(text=string))

}



## Evaluate the (string) keys of a list as R expressions
eval.keys <- function(list) {

    lapply(list, parse.and.eval)

}




## Parse CLI arguments.
get.args <- function() {

    p <- ArgumentParser(
        description = "Concept script for decomposing signals with EMD.")

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

    add_argument("--fourier-method", choices = c("FFT", "LSP"),
        default = "FFT",
        help = "Method used for Fourier transformation.")

    add_argument("--emd-method", choices = c("emd", "semd"),
        default = "emd",
        help = "Method used for Empirical Mode Decomposition")
    add_argument("--emd-args", action = "append", nargs = 2,
        help = "Provide a key value pair to pass as an argument to EMD method.")


    args <- p$parse_args()

    args$signal_type <- signal.types[[args$signal_type]]
    args$noise_type  <- noise.types[[args$noise_type]]
    args$emd_method  <- emd.methods[[args$emd_method]]
    args$emd_args    <- eval.keys(key.val.list(args$emd_args))

    args

}

## Main program logic.
main <- function() {

    # Parse command line arguments.
    args <- get.args()

    # Fix the random seed for reproducibility.
    set.seed(args$seed)

    # Create output directory if specified.
    if(!is.null(args$output)) {
        dir.create(args$output, recursive = TRUE)
    }

    # Generate a signal to decompose.
    signal <- generate.signal(args$signal_type,
                              period = args$period,
                              noise.sd = args$noise_sd,
                              precision = args$time_step,
                              num.cycles = args$num_cycles)

    # Decompose signal into IMFs and residue with EMD.
    emd.result <- decompose.signal(signal$t, signal$observation,
                                   args$emd_method, args$emd_args)
    # Take the Fourier transform of each IMF.
    fourier.transform <- fourier.transform(signal$t, emd.result$imf,
                                           method = args$fourier_method)
    # Determine mutual information of adjacent IMFs.
    mutual.information <- apply(fourier.transform$phase.spectra, 1, rbind)
    print(mutual.information)
    print(cor(t(mutual.information))^2)
    mutual.information <- makemim(mutual.information)
    mutual.information <-
        diag(mutual.information[-length(mutual.information), -1])
    # Partition the IMFs into deterministic and non-deterministic components,
    # and combine the deterministic ones into a single signal.
    deterministic.signal <- extract.determinism(emd.result, fourier.transform,
                                                criteria = args$criteria)
    # Determine fit statistics.
    stats <- fit.statistics(signal, deterministic.signal)

    # Create plots if output directory specified.
    if(!is.null(args$output)) {
        plot.signal(signal, output = args$output)
        plot.determinism(signal, deterministic.signal, output = args$output)
        plot.emd(signal, emd.result, output = args$output)
        plot.amplitude.spectra(fourier.transform, output = args$output)
        plot.phase.spectra(fourier.transform, output = args$output)
        plot.imf.partitions(signal, emd.result, output = args$output)
        plot.mutual.information(mutual.information, output = args$output)
    }

    # Display fit statistics.
    print(stats, row.names = FALSE)

}


main()
