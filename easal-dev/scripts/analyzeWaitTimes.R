library(tidyverse)

toMicroSeconds <- function(strr) {
  strr1 <- as.character(strr)
  spltime1 <- unlist(strsplit(strr1, '[.]'))
  spltime2 <- as.numeric(unlist(strsplit(spltime1[1],":")))
  tseconds <- spltime2[1] * 3600 + spltime2[2] * 60 + spltime2[3]
  micro.seconds <- tseconds * 1000000 + as.numeric(spltime1[2])
  return(micro.seconds)
}

fixSampler_R_BoundaryColumns <- function(row) {
  if (row["message"] == "Sampler_R_Boundary") {
    row["boundary"] <- row["flip"]
    row["flip"] <- row["cayley_point"]
    row["cayley_point"] <- row["graph"]
    row["graph"] <- NA
  }
  row
}

# Message coding
# Sampling_S_Result = "SamplerActor sent sampling results for node:"
# Sampler_S_Done = "SamplerActor about to send DONE message for node:"
# ABA_R_Done = "AtlasBuilderActor received DONE message from:"
# ABA_S_Kill = "AtlasBuilderActor Sent Kill message to:"
# Sampler_R_Kill = "SamplerActor received KILL message for node:"
# Sampler_R_Kill - Sampler_S_Done
#
# Takes a space delimeted file as input
waitInSamplerKills <- function(filename) {
  data3 <- read_table2(filename,col_types = cols(time_stamp=col_character()))
  data3$time_stamp <- as.character(data3$time_stamp)
  data3$time_stamp <- sapply(data3$time_stamp,toMicroSeconds)
  data3 <- data3 %>%  pivot_wider(names_from = message, values_from = time_stamp)
  wait_in_sampler_kill <- data3$Sampler_R_Kill - data3$Sampler_S_Done
  return(data.frame(time_stamp=data3$Sampler_S_Done,wait_in_sampler_kill=wait_in_sampler_kill))
}


samplerLifeTimes <- function(filename) {
  data3 <- read_table2(filename,col_types = cols(time_stamp=col_character()))
  data3$time_stamp <- as.character(data3$time_stamp)
  data3$time_stamp <- sapply(data3$time_stamp,toMicroSeconds)
  data3 <- data3 %>%  pivot_wider(names_from = message, values_from = time_stamp)
  sampler_life_time <- data3$Sampler_R_Kill - data3$Sampler_R_Start
  return(data.frame(time_stamp=data3$Sampler_R_Start,sampler_life_time=sampler_life_time))
}

# Message coding
# Sampler_S_Region="SamplerActor sent a new region to AB Actor from node:"
# ABA_R_Region="Received New Region from:"
# Sampler_R_Boundary="SamplerActor received boundary message for node:"
# ABA_S_Boundary="AtlasBuilderActor sending boundary message to:"
#
# Takes a space delimeted file as input
waitInBoundarySend <- function(filename) {
  data4 <- read_table2(filename,col_types = cols(time_stamp=col_character()))
  data4 <- as_tibble(t(apply(data4,1,fixSampler_R_BoundaryColumns)))
  data4$time_stamp <- as.character(data4$time_stamp)
  data4$time_stamp <- sapply(data4$time_stamp,toMicroSeconds)
  data4$node <- as.numeric(data4$node)
  data4$flip <- as.numeric(data4$flip)
  data4$cayley_point <- as.numeric(data4$cayley_point)
  data4$boundary <- as.numeric(data4$boundary)
  data4$conc_col <- paste0(data4$node,"_",data4$graph,"_",data4$cayley_point,"_",data4$flip)
  data4 <- data4 %>% select(-flip, -graph, -boundary, -node, -cayley_point) %>% filter(message != "Sampler_R_Boundary") %>%
                pivot_wider(names_from = message, values_from = time_stamp)
  wait_in_boundary_send <- data4$ABA_S_Boundary - data4$Sampler_S_Region
  return(data.frame(time_stamp=data4$Sampler_S_Region,wait_in_boundary_send=wait_in_boundary_send))
}

# Takes a space delimeted file as input
# Ignores "Sampler_R_Boundary" messages for speedup
waitInBoundarySendFast <- function(filename) {
  data4<-read_table2(filename,col_types = cols(time_stamp=col_character()))
  #data4<- data4 %>% filter(message != "Sampler_R_Boundary")
  data4$time_stamp <- sapply(data4$time_stamp,toMicroSeconds)
  data4$conc_col <- paste0(data4$node,"_",data4$graph,"_",data4$cayley_point,"_",data4$flip)
  data4 <- data4 %>% select(-flip, -graph, -boundary, -node, -cayley_point)  %>%
    pivot_wider(names_from = message, values_from = time_stamp)
  wait_in_boundary_send <- data4$ABA_S_Boundary - data4$Sampler_S_Region
  return(data.frame(time_stamp=data4$Sampler_S_Region,wait_in_boundary_send=wait_in_boundary_send))
}
