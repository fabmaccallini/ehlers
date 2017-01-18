# Autocorrelation Periodogram - John F. Ehlers - Technical Analysis of Stocks and Commodities, Sep 2016
# ref: http://technical.traders.com/content/TTlink.asp?mo=07&yr=2016

autocorrPeriodogram <- function(x, period1 = 10, period2 = 48, avgLength = 3) {
    # high pass filter
    # NOTE: changed from previous version in "Cycle Analytics for Traders"
    alpha1 <- (1 - sin(2 * pi / period2)) / cos(2 * pi / period2)
    hp <- (1 - alpha1) / 2 * (x - lag(x))
    hp <- hp[-1]
    hp <- filter(hp, alpha1, method = "recursive")
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
    # Pearson correlation for each value of lag
    autocorr <- matrix(0, period2, length(filt))
    for (lag in 2: period2) {
        # Set the average length as M
        if (avgLength == 0) M <- lag
        else M <- avgLength
        autocorr[lag, ] <- runCor(filt, lag(filt, lag), M)
    }
    autocorr[is.na(autocorr)] <- 0
    # Discrete Fourier Transform for each autocorrelation
    # The sum of the squares of each value represents relative power at each period
    cosinePart <- sinePart <- sqSum <- R <- Pwr <- matrix(0, period2, length(filt))
    for (period in period1: period2) {
        for (N in 2: period2) {
            cosinePart[period, ] = cosinePart[period, ] + autocorr[N, ] * cos(N * 2 * pi / period)
            sinePart[period, ] = sinePart[period, ] + autocorr[N, ] * sin(N * 2 * pi / period)
        }
        sqSum[period, ] = cosinePart[period, ] ^ 2 + sinePart[period, ] ^ 2
        R[period, ] <- EMA(sqSum[period, ] ^ 2, ratio = 0.2)
    }
    R[is.na(R)] <- 0
    # Normalising Power
    maxPwr <- rep(0, length(filt))
    for(period in period1: period2) {
        for (i in 1: length(filt)) {
            if (R[period, i] >= maxPwr[i]) maxPwr[i] <- R[period, i]
            # else maxPwr[i] <- K * maxPwr[i]
        }
    }
    for(period in 2: period2) {
        Pwr[period, ] <- R[period, ] / maxPwr
        # to enhance resolution Pwr can be elevated to a power
    }
    # Compute the dominant cycle using the Center of Gravity of the spectrum
    Spx <- Sp <- rep(0, length(filter))
    for(period in period1: period2) {
        Spx <- Spx + period * Pwr[period, ] * (Pwr[period, ] >= 0.5)
        Sp <- Sp + Pwr[period, ] * (Pwr[period, ] >= 0.5)
    }
    dominantCycle <- Spx / Sp
    dominantCycle[is.nan(dominantCycle)] <- 1
    
    heatmap(Pwr, Rowv = NA, Colv = NA, na.rm = TRUE, labCol = "", add.expr = lines(dominantCycle, col = 'blue'))
}

# Autocorrelation Reversals (code 8.4)