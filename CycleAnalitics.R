### Cycle Analytics for Traders: Advanced Technical Trading Concepts (2013, Wiley Trading)

library (dyplr)
library (quantmod)
getSymbols.extra("SPY", from = '2015-01-01', to = '2016-10-01')
x <- Cl(SPY)

## Chapter 3

# 2-pole modified Butterworth
bwMod2 <- function (x, period = 10) {
    a1 <- exp(-sqrt(2) * pi / period)
    b1 <- 2 * a1 * cos(sqrt(2) * 1.25 * pi / period)
    c2 <- b1
    c3 <- -a1 * a1
    c1 <- 1 - c2 - c3
    filt <- c1 * x
    filt <- filter(filt, c(c2, c3), method = "recursive")
    filt <- xts(filt, order.by = index(x))
    return(filt)
}

# 3-pole modified Butterworth
bwMod3 <- function (x, period = 10) {
    a1 <- exp(-pi / period)
    b1 <- 2 * a1 * cos(1.738 * pi / period)
    c1 <- a1 * a1
    d2 <- b1 + c1
    d3 <- -(c1 + b1 * c1)
    d4 <- c1 * c1
    d1 <- 1 - d2 - d3 - d4
    filt <- d1 * x
    filt <- filter(filt, c(d2, d3, d4), method = "recursive")
    filt <- xts(filt, order.by = index(x))
    return(filt)
}

# Supersmoother
superSmoother <- function (x, period = 10) {
    a1 <- exp(-sqrt(2) * pi / period)
    b1 <- 2 * a1 * cos(sqrt(2) * pi / period)
    c2 <- b1
    c3 <- -a1 * a1
    c1 <- 1 - c2 - c3
    filt <- c1 * (x + lag(x)) / 2
    leadNAs <- sum(is.na(filt))
    filt <- filt[-c(1: leadNAs)]
    filt <- filter(filt, c(c2, c3), method = "recursive")
    filt <- c(rep(NA, leadNAs), filt)
    filt <- xts(filt, order.by = index(x))
    return(filt)
}


## Chapter 4

# Decycler (code 4.1)
decycler <- function (x, cutoff = 60) {
    alpha1 <- (cos(2 * pi / cutoff) + sin(2 * pi / cutoff) - 1) / cos(2 * pi / cutoff)
    decycl <- alpha1 * (x + lag(x)) / 2
    leadNAs <- sum(is.na(decycl))
    decycl <- decycl[-c(1: leadNAs)]
    decycl <- filter(decycl, 1 - alpha1, method = "recursive")
    decycl <- c(rep(NA, leadNAs), decycl)
    decycl <- xts(decycl, order.by = index(x))
    return(filt)
}

# Decycler Oscillator (code 4.2)
decyclerOscillator <- function (x, cutoff1 = 30, cutoff2 = 60) {
    alpha1 <- (cos(sqrt(2) * pi / cutoff1) + sin(sqrt(2) * pi / cutoff1) - 1) / cos(sqrt(2) * pi / cutoff1)
    alpha2 <- (cos(sqrt(2) * pi / cutoff2) + sin(sqrt(2) * pi / cutoff2) - 1) / cos(sqrt(2) * pi / cutoff2)
    hp1 <- (1 - alpha1 / 2) ^ 2 * (x - 2 * lag(x, 1) + lag(x, 2))
    leadNAs <- sum(is.na(hp1))
    hp1 <- hp1[-c(1: leadNAs)]
    hp1 <- filter(hp1, c(2 * (1 - alpha1), -(1 - alpha1) ^ 2), method = "recursive")
    hp1 <- c(rep(NA,leadNAs), hp1)
    hp2 <- (1 - alpha2 / 2) ^ 2 * (x - 2 * lag(x, 1) + lag(x, 2))
    leadNAs <- sum(is.na(hp2))
    hp2 <- hp2[-c(1: leadNAs)]
    hp2 <- filter(hp2, c(2 * (1 - alpha2), -(1 - alpha2) ^ 2), method = "recursive")
    hp2 <- c(rep(NA, leadNAs), hp2)    
    decycl <- hp2 - hp1
    decycl <- xts(decycl, order.by = index(x))
    return(filt)
}


## Chapter 5

# Automatic Gain Control (fast attack - slow decay AGC)
AGC <- function (loCutoff = 10, hiCutoff = 48, slope = 1.5) {   
    accSlope = -slope # acceptableSlope = 1.5 dB
    ratio = 10 ^ (accSlope / 20)
    if ((hiCutoff - loCutoff) > 0)
        factor <-  ratio ^ (2 / (hiCutoff - loCutoff));
    return (factor)
}

# Band-Pass Filter (code 5.1)
bandPassFilter <- function (x, period = 20, bandwidth = 0.3) {
    alpha2 <- (cos(0.25 * bandwidth * 2 * pi / period) + sin(0.25 * bandwidth * 2 * pi / period) - 1) / cos(0.25 * bandwidth * 2 * pi / period)
    hp <- (1 + alpha2 / 2) * (x - lag(x))
    hp <- hp[-1]
    hp <- filter(hp, (1 - alpha2), method = "recursive")
    hp <- xts(c(NA, hp), order.by = index(x))
    beta1 <- cos(2 * pi * period)
    gamma1 <- 1 / cos(bandwidth * 2 * pi / period)
    alpha1 <- gamma1 - sqrt(gamma1 ^ 2 - 1)
    bp <- (1 - alpha1) / 2 * (hp - lag(hp, 2))
    bp <- bp[-c(1: 3)]
    bp <- filter(bp, c(beta1 * (1 + alpha1), -alpha1), method = "recursive")
    bp <- c(0, 0, 0, bp)
    peak <- 0; sig <- rep(0, length(bp))
    for (i in 1: length(bp)) {
        peak <- AGC(10, 48, 1.5) * peak
        if (abs(bp[i]) > peak) peak <- abs(bp[i])
        if (peak != 0) sig[i] <- bp[i] / peak
    }
    bp <- xts(bp, order.by = index(x))
    sig <- xts(sig, order.by = index(x))
    alpha2 <- (cos(1.5 * bandwidth * 2 * pi / period) + sin(1.5 * bandwidth * 2 * pi / period) - 1) / cos(1.5 * bandwidth * 2 * pi / period)
    trig <- (1 + alpha2 / 2) * (sig - lag(sig))
    trig <- trig[-1]
    trig <- filter(trig, (1 - alpha2), method = "recursive")
    trig <- c(NA, trig)
    trig <- xts(trig, order.by = index(x))
    return(data.frame(bp, sig, trig))
}

# Dominant Cycle (code 5.2)
bandPassFilter <- function (x, period = 20, bandwidth = 0.3) {
    alpha2 <- (cos(0.25 * bandwidth * 2 * pi / period) + sin(0.25 * bandwidth * 2 * pi / period) - 1) / cos(0.25 * bandwidth * 2 * pi / period)
    hp <- (1 + alpha2 / 2) * (x - lag(x))
    hp <- hp[-1]
    hp <- filter(hp, (1 - alpha2), method = "recursive")
    hp <- xts(c(NA, hp), order.by = index(x))
    beta1 <- cos(2 * pi * period)
    gamma1 <- 1 / cos(bandwidth * 2 * pi / period)
    alpha1 <- gamma1 - sqrt(gamma1 ^ 2 - 1)
    bp <- 0.5 * (1 - alpha1) * (hp - lag(hp, 2))
    bp <- bp[-c(1: 3)]
    bp <- filter(bp, c(beta1 * (1 + alpha1), -alpha1), method = "recursive")
    bp <- c(0, 0, 0, bp)
    peak <- count <- 0; sig <- dc <- rep(0, length(bp))
    for (i in 1: length(bp)) {
        if (abs(bp[i]) > peak) peak <- abs(bp[i])
        else peak <- 0.991 * peak
        if (peak != 0) sig[i] <- bp[i] / peak
        dc[i] <- max(6, dc[max(1, i - 1)])
        count <- count + 1
        if ((sig[i] > 0) & (sig[max(1, i - 1)] < 0) | ((sig[i] < 0) & (sig[max(1, i - 1)] > 0))) {
            dc[i] <- 2 * count
            if ((2 * count) > (1.25 * dc[max(1, i - 1)])) dc[i] <- 1.25 * dc[max(1, i - 1)]
            if ((2 * count) < (0.8 * dc[max(1, i - 1)])) dc[i] <- 0.8 * dc[max(1, i - 1)]
            count <- 0
        }
    }
    return(xts(dc, order.by = index(x)))
}


## Chapter 6

# Hurst Coefficient (code 6.1)
hursCoefficient <- function (x, period = 20) {
    if (period %% 2 == 1) period <- period - 1 # period must be an even number
    a1 <- exp(-sqrt(2) * pi / period)
    b1 <- 2 * a1 * cos(sqrt(2) * pi / period)
    c2 <- b1
    c3 <- -a1 * a1
    c1 <- 1 - c2 - c3
    require(TTR)
    n3 <- (runMax(x, period) - runMin(x, period)) / period
    n1 <- (runMax(x, period / 2) - runMin(x, period / 2)) / (period / 2)
    x2 <- lag(x, period / 2)
    n2 <- (runMax(x2, period / 2) - runMin(x2, period / 2)) / (period / 2)
    dim <- (log(n1 + n2) - log(n3)) / log(2)
    hurst <- 2 - dim
    hurstSmooth <- c1 * (hurst + lag(hurst)) / 2
    leadNAs <- sum(is.na(hurstSmooth))
    hurstSmooth <- hurstSmooth[-c(1: leadNAs)]
    hurstSmooth <- filter(hurstSmooth, c(c2, c3), method = "recursive")
    hurstSmooth <- c(rep(NA, leadNAs), hurstSmooth)
}


## Chapter 7

# Roofing Filter (code 7.1)
roofingFilter <- function (x, period1 = 10, period2 = 48) {
    # high pass filter
    alpha1 <- (cos(2 * pi / period2) + sin(2 * pi / period2) - 1) / cos(2 * pi / period2)
    hp <- (1 - alpha1 / 2) * (x - lag(x))
    hp <- hp[-1]
    hp <- filter(hp, (1 - alpha1), method = "recursive")
    hp <- c(NA, hp)
    hp <- xts(hp, order.by = index(x))
    # super smoother
    a1 <- exp(-sqrt(2) * pi / period1)
    b1 <- 2 * a1 * cos(sqrt(2) * pi / period1)
    c2 <- b1
    c3 <- -a1 * a1
    c1 <- 1 - c2 - c3
    filt <- c1 * (hp + lag(hp)) / 2
    leadNAs <- sum(is.na(filt))
    filt <- filt[-c(1: leadNAs)]
    filt <- filter(filt, c(c2, c3), method = "recursive")
    filt <- c(rep(NA, leadNAs), filt)
    return(xts(filt, order.by = index(x)))
}

# Zero Mean Roofing Filter (code 7.2)
ZMroofingFilter <- function (x, period1 = 10, period2 = 48) {
    # high pass filter
    alpha1 <- (cos(2 * pi / period2) + sin(2 * pi / period2) - 1) / cos(2 * pi / period2)
    hp <- (1 - alpha1 / 2) * (x - lag(x))
    hp <- hp[-1]
    hp <- filter(hp, (1 - alpha1), method = "recursive")
    hp <- c(NA, hp)
    hp <- xts(hp, order.by = index(x))
    # super smoother
    a1 <- exp(-sqrt(2) * pi / period1)
    b1 <- 2 * a1 * cos(sqrt(2) * pi / period1)
    c2 <- b1
    c3 <- -a1 * a1
    c1 <- 1 - c2 - c3
    filt <- c1 * (hp + lag(hp)) / 2
    leadNAs <- sum(is.na(filt))
    filt <- filt[-c(1: leadNAs)]
    filt <- filter(filt, c(c2, c3), method = "recursive")
    filt <- c(rep(NA, leadNAs), filt)
    filt <- xts(filt, order.by = index(x))
    filt2 <- (1 - alpha1 / 2) * (filt - lag(filt))
    leadNAs <- sum(is.na(filt2))
    filt2 <- filt2[-c(1: leadNAs)]
    filt2 <- filter(filt2, (1 - alpha1), method = "recursive")
    filt2 <- c(rep(NA, leadNAs), filt2)
    return(xts(filt2, order.by = index(x)))
}

# Roofing Filter Indicator (code 7.3)
ZMroofingFilter <- function (x, period1 = 40, period2 = 80) {
    # high pass filter
    alpha1 <- (cos(sqrt(2) * pi / period2) + sin(sqrt(2) * pi / period2) - 1) / cos(sqrt(2) * pi / period2)
    hp <- (1 - alpha1 / 2) ^ 2 * (x - 2 * lag(x) + lag(x, 2))
    hp <- hp[-c(1, 2)]
    hp <- filter(hp, (1 - alpha1), method = "recursive")
    hp <- c(NA, NA, hp)
    hp <- xts(hp, order.by = index(x))
    # super smoother
    a1 <- exp(-sqrt(2) * pi / period1)
    b1 <- 2 * a1 * cos(sqrt(2) * pi / period1)
    c2 <- b1
    c3 <- -a1 * a1
    c1 <- 1 - c2 - c3
    filt <- c1 * (hp + lag(hp)) / 2
    leadNAs <- sum(is.na(filt))
    filt <- filt[-c(1: leadNAs)]
    filt <- filter(filt, c(c2, c3), method = "recursive")
    filt <- c(rep(NA, leadNAs), filt)
    filt <- xts(filt, order.by = index(x))
    filt2 <- (1 - alpha1 / 2) * (filt - lag(filt))
    leadNAs <- sum(is.na(filt2))
    filt2 <- filt2[-c(1: leadNAs)]
    filt2 <- filter(filt2, (1 - alpha1), method = "recursive")
    filt2 <- c(rep(NA, leadNAs), filt2)
    return(xts(filt2, order.by = index(x)))
}

# Modified Stochastic (code 7.4)
modifiedStochastic <- function (x, window = 20, period1 = 10, period2 = 48) {
    # high pass filter
    alpha1 <- (cos(sqrt(2) * pi / period2) + sin(sqrt(2) * pi / period2) - 1) / cos(sqrt(2) * pi / period2)
    hp <- (1 - alpha1 / 2) ^ 2 * (x - 2 * lag(x) + lag(x, 2))
    hp <- hp[-c(1, 2)]
    hp <- filter(hp, (1 - alpha1), method = "recursive")
    hp <- c(NA, NA, hp)
    hp <- xts(hp, order.by = index(x))
    # super smoother
    a1 <- exp(-sqrt(2) * pi / period1)
    b1 <- 2 * a1 * cos(sqrt(2) * pi / period1)
    c2 <- b1
    c3 <- -a1 * a1
    c1 <- 1 - c2 - c3
    filt <- c1 * (hp + lag(hp)) / 2
    leadNAs <- sum(is.na(filt))
    filt <- filt[-c(1: leadNAs)]
    filt <- filter(filt, c(c2, c3), method = "recursive")
    filt <- c(rep(NA, leadNAs), filt)
    filt <- xts(filt, order.by = index(x))
    mstoc <- (filt - runMin(filt, window)) / (runMax(filt, window) - runMin(filt, window))
    mstoc <- c1 * (mstoc + lag(mstoc)) / 2
    leadNAs <- sum(is.na(mstoc))
    mstoc <- mstoc[-c(1: leadNAs)]
    mstoc <- filter(mstoc, c(c2, c3), method = "recursive")
    mstoc <- c(rep(NA, leadNAs), mstoc)
    return(xts(mstoc, order.by = index(x)))
}

# Modified RSI (code 7.5)
modifiedRSI <- function (x, window = 20, period1 = 10, period2 = 48) {
    # high pass filter
    alpha1 <- (cos(sqrt(2) * pi / period2) + sin(sqrt(2) * pi / period2) - 1) / cos(sqrt(2) * pi / period2)
    hp <- (1 - alpha1 / 2) ^ 2 * (x - 2 * lag(x) + lag(x, 2))
    hp <- hp[-c(1, 2)]
    hp <- filter(hp, (1 - alpha1), method = "recursive")
    hp <- c(NA, NA, hp)
    hp <- xts(hp, order.by = index(x))
    # super smoother
    a1 <- exp(-sqrt(2) * pi / period1)
    b1 <- 2 * a1 * cos(sqrt(2) * pi / period1)
    c2 <- b1
    c3 <- -a1 * a1
    c1 <- 1 - c2 - c3
    filt <- c1 * (hp + lag(hp)) / 2
    leadNAs <- sum(is.na(filt))
    filt <- filt[-c(1: leadNAs)]
    filt <- filter(filt, c(c2, c3), method = "recursive")
    filt <- c(rep(NA, leadNAs), filt)
    filt <- xts(filt, order.by = index(x))
    posDiff <- negDiff <- diff(filt)
    posDiff[posDiff < 0] <- 0
    negDiff[negDiff > 0] <- 0
    posSum <- runSum(posDiff, window)
    negSum <- runSum(negDiff, window)
    mrsi <- posSum / (posSum - negSum)
    mrsi <- c1 * (mrsi + lag(mrsi)) / 2
    leadNAs <- sum(is.na(mrsi))
    mrsi <- mrsi[-c(1: leadNAs)]
    mrsi <- filter(mrsi, c(c2, c3), method = "recursive")
    mrsi <- c(rep(NA, leadNAs), mrsi)
    return(xts(mrsi, order.by = index(x)))
}


## Chapter 8

# Autocorrelation Indicator (code 8.2)
autocorrIndicator <- function (x, period1 = 10, period2 = 48, avgLength = 3) {
    # high pass filter
    alpha1 <- (cos(sqrt(2) * pi / period2) + sin(sqrt(2) * pi / period2) - 1) / cos(sqrt(2) * pi / period2)
    hp <- (1 - alpha1 / 2) ^ 2 * (x - 2 * lag(x) + lag(x, 2))
    hp <- hp[-c(1, 2)]
    hp <- filter(hp, (1 - alpha1), method = "recursive")
    hp <- c(NA, NA, hp)
    hp <- xts(hp, order.by = index(x))
    # super smoother
    a1 <- exp(-sqrt(2) * pi / period1)
    b1 <- 2 * a1 * cos(sqrt(2) * pi / period1)
    c2 <- b1
    c3 <- -a1 * a1
    c1 <- 1 - c2 - c3
    filt <- c1 * (hp + lag(hp)) / 2
    leadNAs <- sum(is.na(filt))
    filt <- filt[-c(1: leadNAs)]
    filt <- filter(filt, c(c2, c3), method = "recursive")
    filt <- c(rep(NA, leadNAs), filt)
    filt <- xts(filt, order.by = index(x))
    # Pearson correlation for each value of lag
    autocorr <- matrix(0, period2, length(filt))
    for (lag in 2: period2) {
        # Set the average length as M
        if (avgLength == 0) M <- lag
        else M <- avgLength
        autocorr[lag, ] <- runCor(filt, lag(filt, lag), M)
    }
    # scale correlations between 0 and 1
    autocorr <- (autocorr + 1) / 2
    heatmap(autocorr, Rowv = NA, Colv = NA, na.rm = TRUE, labCol = "")
}

# Autocorrelation Periodogram (code 8.3)
autocorrPeriodogram <- function (x, period1 = 10, period2 = 48, avgLength = 3) {
    # high pass filter
    alpha1 <- (cos(sqrt(2) * pi / period2) + sin(sqrt(2) * pi / period2) - 1) / cos(sqrt(2) * pi / period2)
    hp <- (1 - alpha1 / 2) ^ 2 * (x - 2 * lag(x) + lag(x, 2))
    hp <- hp[-c(1, 2)]
    hp <- filter(hp, (1 - alpha1), method = "recursive")
    hp <- c(NA, NA, hp)
    hp <- xts(hp, order.by = index(x))
    # super smoother
    a1 <- exp(-sqrt(2) * pi / period1)
    b1 <- 2 * a1 * cos(sqrt(2) * pi / period1)
    c2 <- b1
    c3 <- -a1 * a1
    c1 <- 1 - c2 - c3
    filt <- c1 * (hp + lag(hp)) / 2
    leadNAs <- sum(is.na(filt))
    filt <- filt[-c(1: leadNAs)]
    filt <- filter(filt, c(c2, c3), method = "recursive")
    filt <- c(rep(NA, leadNAs), filt)
    filt <- xts(filt, order.by = index(x))
    # Pearson correlation for each value of lag
    autocorr <- matrix(0, period2, length(filt))
    for (lag in 2: period2) {
        # Set the average length as M
        if (avgLength == 0) M <- lag
        else M <- avgLength
        autocorr[lag, ] <- runCor(filt, lag(filt, lag), M)
    }
    autocorr[is.na(autocorr)] <- 0
    # Discrete Fourier transform
    # Correlate autocorrelation values with the cosine and sine of each period of interest
    # The sum of the squares of each value represents relative power at each period
    cosinePart <- sinePart <- sqSum <- R <- Pwr <- matrix(0, period2, length(filt))
    for (period in period1: period2) {
        for (N in 2: period2) {
            cosinePart[period, ] = cosinePart[period, ] + autocorr[N, ] * cos(2 * N * pi / period)
            sinePart[period, ] = sinePart[period, ] + autocorr[N, ] * sin(2 * N * pi / period)
        }
        sqSum[period, ] = cosinePart[period, ] ^ 2 + sinePart[period, ] ^ 2
        R[period, ] <- EMA(sqSum[period, ] ^ 2, ratio = 0.2)
    }
    R[is.na(R)] <- 0
    # Normalising Power
    K <- AGC(period1, period2, 1.5)
    maxPwr <- rep(0, length(filt))
    for(period in period1: period2) {
        for (i in 1: length(filt)) {
            if (R[period, i] >= maxPwr[i]) maxPwr[i] <- R[period, i]
            # else maxPwr[i] <- K * maxPwr[i]
        }
    }
    for(period in 2: period2) {
        Pwr[period, ] <- R[period, ] / maxPwr
    }
    # Compute the dominant cycle using the Center of Gravity of the spectrum
    Spx <- Sp <- rep(0, length(filter))
    for(period in period1: period2) {
        Spx <- Spx + period * Pwr[period, ] * (Pwr[period, ] >= 0.5)
        Sp <- Sp + Pwr[period, ] * (Pwr[period, ] >= 0.5)
    }
    dominantCycle <- Spx / Sp
    dominantCycle[is.nan(dominantCycle)] <- 0
    heatmap(Pwr, Rowv = NA, Colv = NA, na.rm = TRUE, labCol = "", add.expr = lines(dominantCycle, col = 'blue'))
}

# Autocorrelation Reversals (code 8.4)
autocorrReversals <- function (x, period1 = 10, period2 = 48, avgLength = 3) {
    # high pass filter
    alpha1 <- (cos(sqrt(2) * pi / period2) + sin(sqrt(2) * pi / period2) - 1) / cos(sqrt(2) * pi / period2)
    hp <- (1 - alpha1 / 2) ^ 2 * (x - 2 * lag(x) + lag(x, 2))
    hp <- hp[-c(1, 2)]
    hp <- filter(hp, (1 - alpha1), method = "recursive")
    hp <- c(NA, NA, hp)
    hp <- xts(hp, order.by = index(x))
    # super smoother
    a1 <- exp(-sqrt(2) * pi / period1)
    b1 <- 2 * a1 * cos(sqrt(2) * pi / period1)
    c2 <- b1
    c3 <- -a1 * a1
    c1 <- 1 - c2 - c3
    filt <- c1 * (hp + lag(hp)) / 2
    leadNAs <- sum(is.na(filt))
    filt <- filt[-c(1: leadNAs)]
    filt <- filter(filt, c(c2, c3), method = "recursive")
    filt <- c(rep(NA, leadNAs), filt)
    filt <- xts(filt, order.by = index(x))
    # Pearson correlation for each value of lag
    autocorr <- matrix(0, period2, length(filt))
    for (lag in 2: period2) {
        # Set the average length as M
        if (avgLength == 0) M <- lag
        else M <- avgLength
        autocorr[lag, ] <- runCor(filt, lag(filt, lag), M)
    }
    autocorr[is.na(autocorr)] <- 0
    # Scale
    autocorr <- 0.5 * (autocorr + 1)
    sumDeltas <- rep(0, length(x))
    for (lag in 2: period2) {
        for (i in 1: length(x)) {
            if (((autocorr[lag, i] > 0.5) & (autocorr[lag, max(1, i - 1)] < 0.5)) | ((autocorr[lag, i] < 0.5) & (autocorr[lag, max(1, i - 1)] > 0.5))) sumDeltas[i] <- sumDeltas[i] + 1
        }
    }
    reversals <- xts(sumDeltas > 24, order.by = index(x))
}
