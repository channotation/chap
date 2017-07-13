
setwd("/home/gianni/Desktop/chap-demo/data")

library(ggplot2)
library(jsonlite)



readLineN <- function(filename, n)
{
  # open connection:
  con <- file(filename, "r")
  
  # skip lines:
  count <- 0
  while( count < n )
  {
    readLines(con, n = 1)
    count <- count + 1 
  }
  
  # read line of interest
  lineN <- readLines(con, n = 1)
  
  # close connection:
  close(con)
  
  # return line of interest:
  return(lineN)
}



filename <- "traj.json"
json.frame <- fromJSON(readLineN(filename, 1), flatten = TRUE)


ggplot(melt(as.data.frame(json.frame$pathSummary)),
       aes(x = variable,
           y = value,
           fill = variable)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ variable,
             scales = "free")





ggplot(as.data.frame(json.frame$solventPositions),
       aes(x = s)) +
  geom_histogram()









filename <- "output.json"
outfile <- fromJSON(readLines(filename, n = 1), flatten = FALSE)




ggplot(as.data.frame(outfile$pathProfile),
       aes(x = s,
           y = radiusMean,
           ymin = radiusMean - radiusSd,
           ymax = radiusMean + radiusSd)) +
  geom_ribbon(alpha = 0.2,
              aes(ymin = radiusMin,
                  ymax = radiusMax)) +
  geom_ribbon(alpha = 0.4) +
  geom_line()



ggplot(as.data.frame(outfile$pathProfile),
       aes(x = s,
           y = densityMean,
           ymin = densityMean - densitySd,
           ymax = densityMean + densitySd)) +
  geom_ribbon(alpha = 0.2,
              aes(ymin = densityMin,
                  ymax = densityMax)) +
  geom_ribbon(alpha = 0.4) +
  geom_line()










