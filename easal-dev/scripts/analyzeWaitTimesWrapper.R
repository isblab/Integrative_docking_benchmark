# This script wraps scripts/extractLogs.R to process logs of multiple runs iteratively

library(tidyverse)

# Read arguments
# args[1]: jobID
# args[2]: output file stub

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Incorrect input arguments!") 
}

print(getwd())
# source script/analyzeWaitTimes.R
source('../scripts/analyzeWaitTimes.R')

# sampler_life_times
filename <- paste0(args[1], "_part1")
sampler_life_time <- samplerLifeTimes(filename)
output <- data.frame(ID=args[1],t(as.matrix(summary(sampler_life_time$sampler_life_time))),Std_Dev=sd(sampler_life_time$sampler_life_time))
outfilename <- paste0(args[2],"_part1.csv")
write.table(output, outfilename, sep = ",", append = T, row.names = F, col.names = F)
write_csv(sampler_life_time, paste0("sampler_life_time_", args[1], ".csv"))


# wait_in_sampler_kill
filename <- paste0(args[1], "_part3")
wait_in_sampler_kill <- waitInSamplerKills(filename)
output <- data.frame(ID=args[1],t(as.matrix(summary(wait_in_sampler_kill$wait_in_sampler_kill))),Std_Dev=sd(wait_in_sampler_kill$wait_in_sampler_kill))
outfilename <- paste0(args[2],"_part3.csv")
write.table(output, outfilename, sep = ",", append = T, row.names = F, col.names = F)
write_csv(wait_in_sampler_kill, paste0("wait_sampler_kill_", args[1], ".csv"))

# wait_in_sampler_kill
filename <- paste0(args[1], "_part4")
wait_in_boundary_send <- waitInBoundarySendFast(filename)
output <- data.frame(ID=args[1],t(as.matrix(summary(wait_in_boundary_send$wait_in_boundary_send))),Std_Dev=sd(wait_in_boundary_send$wait_in_boundary_send))
outfilename <- paste0(args[2],"_part4.csv")
write.table(output, outfilename, sep = ",", append = T, row.names = F, col.names = F)
write_csv(wait_in_boundary_send, paste0("wait_boundary_send_", args[1], ".csv"))
