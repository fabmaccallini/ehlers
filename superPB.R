# Super Passband Filter - John F. Ehlers - Technical Analysis of Stocks and Commodities, Jul 2016
# ref: http://technical.traders.com/content/TTlink.asp?mo=07&yr=2016

superPB <- function(x, period1 = 10, period2 = 48, plot = TRUE) {
    a1 <- 5 / period1
    a2 <- 5 / period2
    pb <- (a1 - a2) * x + (a2 * (1 - a1) - a1 * (1 - a2)) * lag(x)
    pb <- pb[-1]
    pb <- filter(pb, c((1 - a1) + (1 - a2), -(1 - a1) * (1 - a2)), method = "recursive")
    pb <- xts(c(NA, pb), order.by = index(x))
    RMS <- runSD(pb, (period1 + period2) / 2)
    if (plot) {
        par(mfrow = c(2, 1))
        plot(x, main = 'Price')
        lines(x - pb, col = 'red')
        plot(pb, main = 'Super Passband Filter')
        lines(2 * RMS, col = 'red')
        lines(-1 * RMS, col = 'red')
        abline(h = 0, col = 'blue')
    }
    return(data.frame(pb, RMS, -RMS))
}

# Strategy
#
# Buy on the filter crossing above its -RMS line
# Short on the filter crossing below its RMS line
# Exit long when the filter either crosses below its RMS or crosses below -RMS (which signifies a false entry signal)
# Cover short when the filter either crosses above its -RMS or crosses above RMS (which signifies a false entry signal)