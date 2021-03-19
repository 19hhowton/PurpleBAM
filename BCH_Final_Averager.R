library("readr")
library("data.table")
library("RcppHMM")
library("randomForest")
library("caret")
library("tidyverse")

# Get file names
setwd("E:/Data")
file_paths = dir("E:/Data") #Include your file path where data is stored here- create directory where all and only all of your data is found

# Review Data- do not run if confident that all data will work with main analysis - DON'T RUN ANYTHING IN THIS SECTION
# Get data tables from csvs
csvs = list()
for (file_path in file_paths) {
  csvs = append(csvs, list(as.data.table(read_csv(file_path))))
}

# Review Beginning Dates
for (idx in 1:length(csvs)) {
  print(paste(file_paths[idx], ": ", csvs[[idx]][1,1]))
}

# Review End Dates
for (idx in 1:length(csvs)) {
  print(paste(file_paths[idx], ": ", csvs[[idx]][nrow(csvs[[idx]]),1]))
}

# Review column names in data sets
for (idx in 1:length(csvs)) {
  print(paste(file_paths[idx], ": "))
  print(colnames(csvs[[idx]]))
}

# Begin Data pre-analysis and modification
# Organize data tables
csvs = list()
for (file_idx in seq(1,length(file_paths), 2)) {
  csvs = append(csvs, list(list()))
  csvs[[length(csvs)]] = append(csvs[[length(csvs)]], list(as.data.table(read_csv(file_paths[file_idx]))))
  csvs[[length(csvs)]] = append(csvs[[length(csvs)]], list(as.data.table(read_csv(file_paths[file_idx + 1]))))
}

# Add datehour column (may take some time)
for (outer_idx in 1:length(csvs)) {
  for (inner_idx in 1:length(csvs[[outer_idx]])){
    csvs[[outer_idx]][[inner_idx]]$datehour <- as.character(cut(as.POSIXct(csvs[[outer_idx]][[inner_idx]]$created_at, format= "%Y-%m-%d %H:%M:%S"), breaks = "hour"))
  }
}

# Remove times from pairs not present in both A and B sensors
for (outer_idx in 1:length(csvs)) {
  print(outer_idx)
  print("")
  print(paste("Items removed from A: ", paste(setdiff(csvs[[outer_idx]][[1]]$datehour, csvs[[outer_idx]][[2]]$datehour), collapse = ", ")))
  print(paste("Items removed from B: ", paste(setdiff(csvs[[outer_idx]][[2]]$datehour, csvs[[outer_idx]][[1]]$datehour), collapse = ", ")))
  matching = intersect(csvs[[outer_idx]][[1]]$datehour, csvs[[outer_idx]][[2]]$datehour)
  csvs[[outer_idx]][[1]] = csvs[[outer_idx]][[1]][csvs[[outer_idx]][[1]]$datehour %in% matching,]
  csvs[[outer_idx]][[2]] = csvs[[outer_idx]][[2]][csvs[[outer_idx]][[2]]$datehour %in% matching,]
  print("")
  print("")
}

# Data check- only run if not sure datetime columns will match between A-B pairs- DON'T RUN
to_be_dropped = list()
for (outer_idx in 1:length(csvs)) {
  print(paste(outer_idx, ": ", identical(levels(as.factor(csvs[[outer_idx]][[1]]$datehour)), levels(as.factor(csvs[[outer_idx]][[2]]$datehour)))))
  to_be_dropped = append(to_be_dropped, identical(levels(as.factor(csvs[[outer_idx]][[1]]$datehour)), levels(as.factor(csvs[[outer_idx]][[2]]$datehour))))
}

# Drop False results from data check- Do not run if above data check was skipped- DON'T RUN
csvs = csvs[as.vector(to_be_dropped, mode = "logical")]

# Find correlations between A-B pairs
to_be_dropped = list()
for (outer_idx in 1:length(csvs)) {
  print(paste(outer_idx, ": ", cor(aggregate(`PM2.5_CF1_ug/m3` ~ datehour, csvs[[outer_idx]][[1]], mean)[2], aggregate(`PM2.5_CF1_ug/m3` ~ datehour, csvs[[outer_idx]][[2]], mean)[2])))
  to_be_dropped = append(to_be_dropped, cor(aggregate(`PM2.5_CF1_ug/m3` ~ datehour, csvs[[outer_idx]][[1]], mean)[2], aggregate(`PM2.5_CF1_ug/m3` ~ datehour, csvs[[outer_idx]][[2]], mean)[2]))
}

# Drop correlations between A-B pairs below 0.95
csvs = csvs[to_be_dropped > 0.95]

# Drop unneccessary columns in data
for (outer_idx in 1:length(csvs)) {
  csvs[[outer_idx]][[1]][, c("created_at", "UptimeMinutes", "entry_id"):=NULL]
  
  csvs[[outer_idx]][[2]][, c("created_at", "UptimeMinutes", "entry_id", "IAQ"):=NULL]
}

# Aggregate Data by datehour (may take some time)
for (outer_idx in 1:length(csvs)) {
  for (inner_idx in 1:2){
    csvs[[outer_idx]][[inner_idx]] = aggregate(csvs[[outer_idx]][[inner_idx]], by = list(csvs[[outer_idx]][[inner_idx]]$datehour), mean)
    csvs[[outer_idx]][[inner_idx]]$datehour = as.character(csvs[[outer_idx]][[inner_idx]]$Group.1)
    csvs[[outer_idx]][[inner_idx]]$Group.1 =NULL
  }
}

# Get and aggregate BAM data
setwd("E:/")
# *IMPORTANT* Delete first 5 rows of Bam file before continuing (no data lost)
bam = as.data.table(read_csv("BAM.CSV"))
bam$datehour <- as.character(cut(as.POSIXct(bam$Time, format= "%m/%d/%Y %H:%M"), breaks = "hour"))
bam$Time =NULL
bam = aggregate(bam, by = list(bam$datehour), mean)
bam$datehour = as.character(bam$Group.1)
bam$Group.1 =NULL

# Integrate BAM errors into Data
for (outer_idx in 1:length(csvs)) {
 for (inner_idx in 1:length(csvs[[outer_idx]])) {
   tempBam = bam[bam$datehour %in% csvs[[outer_idx]][[inner_idx]]$datehour,]
   csvs[[outer_idx]][[inner_idx]]$Error_raw = csvs[[outer_idx]][[inner_idx]]$`PM2.5_CF1_ug/m3` / tempBam$`ConcHR(ug/m3)`
   csvs[[outer_idx]][[inner_idx]]$Error_cat = cut(csvs[[outer_idx]][[inner_idx]]$Error_raw, c(-Inf, .05, .15, .25, .35, .45, .55, .65, .75, .85, 
                                              .95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, Inf), labels = c("<.05", ".05 - .15", 
                                              ".15 - .25", ".25 - .35", ".35 - .45", ".45 - .55", ".55 - .65", ".65 - .75", ".75 - .85",
                                              ".85 - .95", ".95 - 1.05", "1.05 - 1.15", "1.15 - 1.25", "1.25 - 1.35", "1.35 - 1.45", "1.45 - 1.55", "1.55 - 1.65",
                                              "1.65 - 1.75", "1.75 - 1.85", "1.85 - 1.95", ">1.95"))
 }
}

# Find BAM beginning date- DON'T RUN
bam[1, 12]

# Find BAM ending date- DON'T RUN
bam[nrow(bam), 12]

# Trim data to least restrictive range- DON'T RUN
for (outer_idx in 1:length(csvs)) {
  for (inner_idx in 1:2){
    csvs[[outer_idx]][[inner_idx]] = 
      csvs[[outer_idx]][[inner_idx]][match("2019-08-19 11:00:00", csvs[[outer_idx]][[inner_idx]]$datehour):match("2019-11-04 23:00:00", csvs[[outer_idx]][[inner_idx]]$datehour),]
  }
}

# Trim bam into least restrictive range
bam = bam[match("2019-08-18 00:00:00", bam$datehour):match("2019-11-04 23:00:00", bam$datehour),]


# Get all hours in least restrictive range needed for analysis- DON'T RUN
analysis_hours = as.character(seq(from= as.POSIXct("2019-08-28 00:00:00", format= "%Y-%m-%d %H:%M:%S"), to= as.POSIXct("2019-10-28 23:00:00", format= "%Y-%m-%d %H:%M:%S"), by= "hour"))

# Test to see if all time values in range are present in BAM and data- DON'T RUN
for (outer_idx in 1:length(csvs)) {
  print(paste(outer_idx, ": "))
  print(setdiff(analysis_hours, csvs[[outer_idx]][[1]]$datehour))
}
print("BAM: ")
print(setdiff(analysis_hours, bam$datehour))

# Get training time range
csvs_training = list()
for (outer_idx in 1:length(csvs)) {
  
  csvs_training = append(csvs_training, list(list()))
  
  csvs_training[[length(csvs_training)]] = append(csvs_training[[length(csvs_training)]], list(csvs[[outer_idx]][[1]][as.POSIXct(csvs[[outer_idx]][[1]]$datehour, format= "%Y-%m-%d %H:%M:%S")
  >= as.POSIXct("2019-08-28 00:00:00", format= "%Y-%m-%d %H:%M:%S") & as.POSIXct(csvs[[outer_idx]][[1]]$datehour, format= "%Y-%m-%d %H:%M:%S") <= as.POSIXct("2019-10-28 00:00:00", format= "%Y-%m-%d %H:%M:%S"),]))
  
  csvs_training[[length(csvs_training)]] = append(csvs_training[[length(csvs_training)]], list(csvs[[outer_idx]][[2]][as.POSIXct(csvs[[outer_idx]][[2]]$datehour, format= "%Y-%m-%d %H:%M:%S")
  >= as.POSIXct("2019-08-28 00:00:00", format= "%Y-%m-%d %H:%M:%S") & as.POSIXct(csvs[[outer_idx]][[2]]$datehour, format= "%Y-%m-%d %H:%M:%S") <= as.POSIXct("2019-10-28 00:00:00", format= "%Y-%m-%d %H:%M:%S"),]))
}

# Remove rows with NA values
for (outer_idx in 1:length(csvs_training)) {
  for (inner_idx in 1:length(csvs_training[[outer_idx]])) {
    csvs_training[[outer_idx]][[inner_idx]] = csvs_training[[outer_idx]][[inner_idx]][complete.cases(csvs_training[[outer_idx]][[inner_idx]]),]
  }
}

# Remove differences between training pairs
to_be_dropped = list()
for (outer_idx in 1:length(csvs_training)) {
  to_be_dropped = c(to_be_dropped, !(nrow(csvs_training[[outer_idx]][[1]]) == 0 | nrow(csvs_training[[outer_idx]][[2]]) == 0))
  if (!(nrow(csvs_training[[outer_idx]][[1]]) == 0 | nrow(csvs_training[[outer_idx]][[2]]) == 0)) {
    print(outer_idx)
    print("")
    print(paste("Items removed from A: ", paste(setdiff(csvs_training[[outer_idx]][[1]]$datehour, csvs_training[[outer_idx]][[2]]$datehour), collapse = ", ")))
    print(paste("Items removed from B: ", paste(setdiff(csvs_training[[outer_idx]][[2]]$datehour, csvs_training[[outer_idx]][[1]]$datehour), collapse = ", ")))
    matching = intersect(csvs_training[[outer_idx]][[1]]$datehour, csvs_training[[outer_idx]][[2]]$datehour)
    csvs_training[[outer_idx]][[1]] = csvs_training[[outer_idx]][[1]][csvs_training[[outer_idx]][[1]]$datehour %in% matching,]
    csvs_training[[outer_idx]][[2]] = csvs_training[[outer_idx]][[2]][csvs_training[[outer_idx]][[2]]$datehour %in% matching,]
    print("")
    print("")
  }
}

csvs_training = csvs_training[as.vector(to_be_dropped, mode = "logical")]

# Filter out training csvs < 168
csvs_training = csvs_training[sapply(csvs_training, function(x) {nrow(x[[1]])}) >= 168]



# Get uninterrupted training segments
training = list()
for (sequence_idx in 1:length(csvs_training)) {
  training = append(training, list(list()))
  training[[length(training)]] = append(training[[length(training)]], list(list()))
  training[[length(training)]] = append(training[[length(training)]], list(list()))
}

for (outer_idx in 1:length(csvs_training)) {
  for (inner_idx in 1:length(csvs_training[[outer_idx]]))  {
    start = 1
    if (nrow(csvs_training[[outer_idx]][[inner_idx]]) != 0) {
    for (date_idx in 1:(nrow(csvs_training[[outer_idx]][[inner_idx]]) - 1)) {
      if (as.POSIXct(csvs_training[[outer_idx]][[inner_idx]][date_idx + 1, "datehour"], format= "%Y-%m-%d %H:%M:%S") 
          - as.POSIXct(csvs_training[[outer_idx]][[inner_idx]][date_idx, "datehour"], format= "%Y-%m-%d %H:%M:%S") != 1){
        training[[outer_idx]][[inner_idx]] = c(training[[outer_idx]][[inner_idx]], list(csvs_training[[outer_idx]][[inner_idx]][start:date_idx,]))
        start = date_idx + 1
      }
    }
    training[[outer_idx]][[inner_idx]] = c(training[[outer_idx]][[inner_idx]], list(csvs_training[[outer_idx]][[inner_idx]][start:nrow(csvs_training[[outer_idx]][[inner_idx]]),]))
    }
  }
}

# Randomly sample from training sequencces to diversify training set
set.seed(42)
for (outer_idx in 1:length(training)) {
  for (inner_idx in 1:length(training[[outer_idx]])) {
    training[[outer_idx]][[inner_idx]] = training[[outer_idx]][[inner_idx]][sapply(training[[outer_idx]][[inner_idx]], nrow) > 168]
  }
}

training_trimmed = list()
for (sequence_idx in 1:length(training)) {
  training_trimmed = append(training_trimmed, list(list()))
  training_trimmed[[length(training_trimmed)]] = append(training_trimmed[[length(training_trimmed)]], list(list()))
  training_trimmed[[length(training_trimmed)]] = append(training_trimmed[[length(training_trimmed)]], list(list()))
}

for (outer_idx in 1:length(training)) {
  if (length(training[[outer_idx]][[1]]) > 0) {
    for (sequence_idx in 1:length(training[[outer_idx]][[1]])) {
      upper = sample(168:nrow(training[[outer_idx]][[1]][[sequence_idx]]), 1)
      lower = sample(1:(upper - 167), 1)
      training_trimmed[[outer_idx]][[1]] = c(training_trimmed[[outer_idx]][[1]], list(as.data.table(training[[outer_idx]][[1]][[sequence_idx]][lower:upper,])))
      training_trimmed[[outer_idx]][[2]] = c(training_trimmed[[outer_idx]][[2]], list(as.data.table(training[[outer_idx]][[2]][[sequence_idx]][lower:upper,])))
      
    }
  }
}

# Drop empty lists
to_be_dropped = list()
for (outer_idx in 1:length(training_trimmed)) {
  to_be_dropped <- c(to_be_dropped, !(length(training_trimmed[[outer_idx]][[1]]) == 0))
}

training_trimmed = training_trimmed[as.vector(to_be_dropped, mode = "logical")]

# Get testing time range
csvs_testing = list()
for (outer_idx in 1:length(csvs)) {
  
  csvs_testing = append(csvs_testing, list(list()))
  
  csvs_testing[[length(csvs_testing)]] = append(csvs_testing[[length(csvs_testing)]], list(csvs[[outer_idx]][[1]][as.POSIXct(csvs[[outer_idx]][[1]]$datehour, format= "%Y-%m-%d %H:%M:%S")
                                                                                                                      < as.POSIXct("2019-08-28 00:00:00", format= "%Y-%m-%d %H:%M:%S") | as.POSIXct(csvs[[outer_idx]][[1]]$datehour, format= "%Y-%m-%d %H:%M:%S") > as.POSIXct("2019-10-28 00:00:00", format= "%Y-%m-%d %H:%M:%S"),]))
  
  csvs_testing[[length(csvs_testing)]] = append(csvs_testing[[length(csvs_testing)]], list(csvs[[outer_idx]][[2]][as.POSIXct(csvs[[outer_idx]][[2]]$datehour, format= "%Y-%m-%d %H:%M:%S")
                                                                                                                      < as.POSIXct("2019-08-28 00:00:00", format= "%Y-%m-%d %H:%M:%S") | as.POSIXct(csvs[[outer_idx]][[2]]$datehour, format= "%Y-%m-%d %H:%M:%S") > as.POSIXct("2019-10-28 00:00:00", format= "%Y-%m-%d %H:%M:%S"),]))
}

# Remove rows with NA values
for (outer_idx in 1:length(csvs_testing)) {
  for (inner_idx in 1:length(csvs_testing[[outer_idx]])) {
    csvs_testing[[outer_idx]][[inner_idx]] = csvs_testing[[outer_idx]][[inner_idx]][complete.cases(csvs_testing[[outer_idx]][[inner_idx]]),]
  }
}

# Remove differences between testing pairs
to_be_dropped = list()
for (outer_idx in 1:length(csvs_testing)) {
  to_be_dropped = c(to_be_dropped, !(nrow(csvs_testing[[outer_idx]][[1]]) == 0 | nrow(csvs_testing[[outer_idx]][[2]]) == 0))
  if (!(nrow(csvs_testing[[outer_idx]][[1]]) == 0 | nrow(csvs_testing[[outer_idx]][[2]]) == 0)) {
    print(outer_idx)
    print("")
    print(paste("Items removed from A: ", paste(setdiff(csvs_testing[[outer_idx]][[1]]$datehour, csvs_testing[[outer_idx]][[2]]$datehour), collapse = ", ")))
    print(paste("Items removed from B: ", paste(setdiff(csvs_testing[[outer_idx]][[2]]$datehour, csvs_testing[[outer_idx]][[1]]$datehour), collapse = ", ")))
    matching = intersect(csvs_testing[[outer_idx]][[1]]$datehour, csvs_testing[[outer_idx]][[2]]$datehour)
    csvs_testing[[outer_idx]][[1]] = csvs_testing[[outer_idx]][[1]][csvs_testing[[outer_idx]][[1]]$datehour %in% matching,]
    csvs_testing[[outer_idx]][[2]] = csvs_testing[[outer_idx]][[2]][csvs_testing[[outer_idx]][[2]]$datehour %in% matching,]
    print("")
    print("")
  }
}

csvs_testing = csvs_testing[as.vector(to_be_dropped, mode = "logical")]

# Filter out testing csvs < 2
csvs_testing = csvs_testing[sapply(csvs_testing, function(x) {nrow(x[[1]])}) >= 2]

# Get uninterrupted testing segments
testing = list()
for (sequence_idx in 1:length(csvs_testing)) {
  testing = append(testing, list(list()))
  testing[[length(testing)]] = append(testing[[length(testing)]], list(list()))
  testing[[length(testing)]] = append(testing[[length(testing)]], list(list()))
}

for (outer_idx in 1:length(csvs_testing)) {
  for (inner_idx in 1:length(csvs_testing[[outer_idx]]))  {
    start = 1
    if (nrow(csvs_testing[[outer_idx]][[inner_idx]]) != 0) {
    for (date_idx in 1:(nrow(csvs_testing[[outer_idx]][[inner_idx]]) - 1)) {
      if (as.POSIXct(csvs_testing[[outer_idx]][[inner_idx]][date_idx + 1, "datehour"], format= "%Y-%m-%d %H:%M:%S") 
          - as.POSIXct(csvs_testing[[outer_idx]][[inner_idx]][date_idx, "datehour"], format= "%Y-%m-%d %H:%M:%S") != 1){
        testing[[outer_idx]][[inner_idx]] = c(testing[[outer_idx]][[inner_idx]], list(csvs_testing[[outer_idx]][[inner_idx]][start:date_idx,]))
        start = date_idx + 1
      }
    }
    testing[[outer_idx]][[inner_idx]] = c(testing[[outer_idx]][[inner_idx]], list(csvs_testing[[outer_idx]][[inner_idx]][start:nrow(csvs_testing[[outer_idx]][[inner_idx]]),]))
    }
  }
}

# Ensure at least 2 items in the testing sequences
for (outer_idx in 1:length(testing)) {
  for (inner_idx in 1:length(testing[[outer_idx]])) {
    if (length(testing[[outer_idx]][[inner_idx]]) > 0) {
      testing[[outer_idx]][[inner_idx]] = testing[[outer_idx]][[inner_idx]][sapply(testing[[outer_idx]][[inner_idx]], nrow) >= 2]
    }
  }
}

# Drop empty lists
to_be_dropped = list()
for (outer_idx in 1:length(testing)) {
  to_be_dropped <- c(to_be_dropped, !(length(testing[[outer_idx]][[1]]) == 0))
}

testing = testing[as.vector(to_be_dropped, mode = "logical")]

# Begin HMM analysis
# Linearize training data
markov_training = list()
for (outer_idx in 1:length(training_trimmed)) {
  for (inner_idx in 1:length(training_trimmed[[outer_idx]])) {
    markov_training = c(markov_training, training_trimmed[[outer_idx]][[inner_idx]])
  }
}

# Linearize testing data
markov_testing = list()
for (outer_idx in 1:length(testing)) {
  for (inner_idx in 1:length(testing[[outer_idx]])) {
    markov_testing = c(markov_training, testing[[outer_idx]][[inner_idx]])
  }
}

# Define parameters
states = c("<.05", ".05 - .15", 
           ".15 - .25", ".25 - .35", ".35 - .45", ".45 - .55", ".55 - .65", ".65 - .75", ".75 - .85",
           ".85 - .95", ".95 - 1.05", "1.05 - 1.15", "1.15 - 1.25", "1.25 - 1.35", "1.35 - 1.45", "1.45 - 1.55", "1.55 - 1.65",
           "1.65 - 1.75", "1.75 - 1.85", "1.85 - 1.95", ">1.95")

# Define initial probabilities
pi = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
for (sequence in markov_training){
  pi[[match(as.character(sequence$Error_cat[[1]]), states)]] = pi[[match(as.character(sequence$Error_cat[[1]]), states)]] + 1
}
pi = pi / sum(pi)

# Define Mu, a matrix of means for the Gausian distributions
mu = matrix(sample(0, length(states), T), ncol = length(states))
for (state_idx in 1:length(states)) {
  state_vals = c()
  for (sequence in markov_training) {
    unistateseq = sequence[as.character(sequence$Error_cat) == states[[state_idx]],]
    state_vals = c(state_vals, unistateseq$`PM2.5_CF1_ug/m3`)
  }
  mu[1,state_idx] = mean(state_vals)
}

# Define A, the transition probabilities matrix
A = matrix(0, nrow = length(states), ncol = length(states))
for (sequence in markov_training) {
  for (state_idx in 1:(length(sequence$Error_cat) - 1)) {
    A[match(sequence$Error_cat[[state_idx]], states), match(sequence$Error_cat[[state_idx + 1]], states)] = 
      A[match(sequence$Error_cat[[state_idx]], states), match(sequence$Error_cat[[state_idx + 1]], states)] + 1
  }
}
A = A / rowSums(A)

# Define Sigma, the covariance matrix of each state
Sigma = array(0, dim = c(1,1,length(states)))
for (state_idx in 1:length(states)) {
  Sigma[,,state_idx] = matrix(1)
}

# Define the HMM
markov_model = verifyModel(list( "Model" = "GHMM",
                                 "StateNames" = states,
                                 "A" = A, 
                                 "Mu" = mu, 
                                 "Sigma" = Sigma, 
                                 "Pi" = pi))

# Generate Predictions for testing sequences
for (sequence_idx in 1:length(markov_testing)) {
  markov_testing[[sequence_idx]]$HMM_cat <- viterbi(markov_model, matrix(markov_testing[[sequence_idx]]$`PM2.5_CF1_ug/m3`, nrow = 1))
}

# Correct Purple Air readings
correction = function(cat, PM2.5) {
  if (cat == "<.05" || cat == ">1.95") {
    return(NA)
  } else if (cat == ".05 - .15") {
    return(PM2.5 * 10)
  } else if (cat == ".15 - .25") {
    return(PM2.5 * 5)
  } else if (cat == ".25 - .35") {
    return(PM2.5 * (1/.3))
  } else if (cat == ".35 - .45") {
    return(PM2.5 * 2.5)
  } else if (cat == ".45 - .55") {
    return(PM2.5 * 2)
  } else if (cat == ".55 - .65") {
    return(PM2.5 * (1/.6))
  } else if (cat == ".65 - .75") {
    return(PM2.5 * (1/.7))
  } else if (cat == ".75 - .85") {
    return(PM2.5 * 1.25)
  } else if (cat == ".85 - .95") {
    return(PM2.5 * (1/.9))
  } else if (cat == ".95 - 1.05") {
    return(PM2.5)
  } else if (cat == "1.05 - 1.15") {
    return(PM2.5 * (1/1.1))
  } else if (cat == "1.15 - 1.25") {
    return(PM2.5 * (1/1.2))
  } else if (cat == "1.25 - 1.35") {
    return(PM2.5 * (1/1.3))
  } else if (cat == "1.35 - 1.45") {
    return(PM2.5 * (1/1.4))
  } else if (cat == "1.45 - 1.55") {
    return(PM2.5 * (1/1.5))
  } else if (cat == "1.55 - 1.65") {
    return(PM2.5 * .625)
  } else if (cat == "1.65 - 1.75") {
    return(PM2.5 * (1/1.7))
  } else if (cat == "1.75 - 1.85") {
    return(PM2.5 * (1/1.8))
  } else if (cat == "1.85 - 1.95") {
    return(PM2.5 * (1/1.9))
  } else {
    print(cat)
    return(NA)
  }
}

for (sequence_idx in 1:length(markov_testing)) {
  markov_testing[[sequence_idx]]$HMM_correction = NA
  for (reading_idx in 1:nrow(markov_testing[[sequence_idx]])) {
    markov_testing[[sequence_idx]][[reading_idx, "HMM_correction"]] = 
      correction(markov_testing[[sequence_idx]][[reading_idx, "HMM_cat"]], markov_testing[[sequence_idx]][[reading_idx, "PM2.5_CF1_ug/m3"]])
  }
}

# Find errors in predictions
for (sequence_idx in 1:length(markov_testing)) {
  tempBam = bam[bam$datehour %in% markov_testing[[sequence_idx]]$datehour,]
  markov_testing[[sequence_idx]]$Absolute_Error = abs(markov_testing[[sequence_idx]]$HMM_correction - tempBam$`ConcHR(ug/m3)`)
  markov_testing[[sequence_idx]]$HMM_correction_Error_raw = markov_testing[[sequence_idx]]$HMM_correction / tempBam$`ConcHR(ug/m3)`
  markov_testing[[sequence_idx]]$HMM_correction_Error_cat = cut(markov_testing[[sequence_idx]]$HMM_correction_Error_raw, c(-Inf, .05, .15, .25, .35, .45, .55, .65, .75, .85, 
                                                                                             .95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, Inf), labels = c("<.05", ".05 - .15", 
                                                                                             ".15 - .25", ".25 - .35", ".35 - .45", ".45 - .55", ".55 - .65", ".65 - .75", ".75 - .85",
                                                                                             ".85 - .95", ".95 - 1.05", "1.05 - 1.15", "1.15 - 1.25", "1.25 - 1.35", "1.35 - 1.45", "1.45 - 1.55", "1.55 - 1.65",
                                                                                              "1.65 - 1.75", "1.75 - 1.85", "1.85 - 1.95", ">1.95"))
}

# Compare original accuracy to HMM accuracy by determing how many data points for each fall into the .95 - 1.05 error category
orig_accuracy = 0
HMM_accuracy = 0
stats_error = vector("numeric")
for (sequence_idx in 1:(length(markov_testing))){
  orig_accuracy = orig_accuracy + sum((markov_testing[[sequence_idx]]$Error_cat == ".95 - 1.05"), na.rm = T)
  HMM_accuracy = HMM_accuracy + sum(markov_testing[[sequence_idx]]$HMM_correction_Error_cat == ".95 - 1.05", na.rm = T)
  stats_error = c(stats_error, markov_testing[[sequence_idx]]$Absolute_Error)
}
print(paste("Original Accuracy (# of values falling in .95 - 1.05 category)", orig_accuracy/sum(sapply(markov_testing, nrow)), sep = " "))
print(paste("HMM Accuracy (# of values falling in .95 - 1.05 category)", HMM_accuracy/sum(sapply(markov_testing, nrow)), sep = " "))
print(paste("Median Error", median(stats_error, na.rm = T), sep = " "))

print(sqrt(mean(stats_error^2, na.rm = T)))

# Begin Random Forest analysis
#Create the Random Forest
library("readr")
library('data.table')
library(randomForest)
library(caret)
library(tidyverse)
library(ranger)
library(ggplot2)



#L&H Adjustments
RF_training = list()
for (sequence_idx in 1:length(training_trimmed)) {
  RF_training = append(RF_training, list(list()))
  RF_training[[length(RF_training)]] = append(RF_training[[length(RF_training)]], list(list()))
  RF_training[[length(RF_training)]] = append(RF_training[[length(RF_training)]], list(list()))
}

for (outer_idx in 1:length(training_trimmed)) {
  
  for (inner_idx in 1:length(training_trimmed[[outer_idx]])){
    if (length(training_trimmed[[outer_idx]][[inner_idx]]) > 0) {
      if (length(training_trimmed[[outer_idx]][[inner_idx]]) == 1) {
        RF_training[[outer_idx]][[inner_idx]] = training_trimmed[[outer_idx]][[inner_idx]][[1]]
      }
      else if (length(training_trimmed[[outer_idx]][[inner_idx]])== 2) {
        RF_training[[outer_idx]][[inner_idx]] = rbind(training_trimmed[[outer_idx]][[inner_idx]][[1]],training_trimmed[[outer_idx]][[inner_idx]][[2]])
        
      }
      
      else {
        temp = rbind(training_trimmed[[outer_idx]][[inner_idx]][[1]],training_trimmed[[outer_idx]][[inner_idx]][[2]])
        
        for (index in 3:length(training_trimmed[[outer_idx]][[inner_idx]])){
          temp = rbind(temp, training_trimmed[[outer_idx]][[inner_idx]][[index]])
        }
        RF_training[[outer_idx]][[inner_idx]] = temp 
      }     
      
    } 
  }
}   

RF_testing = list()
for (sequence_idx in 1:length(testing)) {
  RF_testing = append(RF_testing, list(list()))
  RF_testing[[length(RF_testing)]] = append(RF_testing[[length(RF_testing)]], list(list()))
  RF_testing[[length(RF_testing)]] = append(RF_testing[[length(RF_testing)]], list(list()))
}

for (outer_idx in 1:length(testing)) {
  
  for (inner_idx in 1:length(testing[[outer_idx]])){
    if (length(testing[[outer_idx]][[inner_idx]]) > 0) {
      if (length(testing[[outer_idx]][[inner_idx]]) == 1) {
        RF_testing[[outer_idx]][[inner_idx]] = testing[[outer_idx]][[inner_idx]][[1]]
      }
      else if (length(testing[[outer_idx]][[inner_idx]])== 2) {
        RF_testing[[outer_idx]][[inner_idx]] = rbind(testing[[outer_idx]][[inner_idx]][[1]],testing[[outer_idx]][[inner_idx]][[2]])
        
      }
      
      else {
        temp = rbind(testing[[outer_idx]][[inner_idx]][[1]],testing[[outer_idx]][[inner_idx]][[2]])
        
        for (index in 3:length(testing[[outer_idx]][[inner_idx]])){
          temp = rbind(temp, testing[[outer_idx]][[inner_idx]][[index]])
        }
        RF_testing[[outer_idx]][[inner_idx]] = temp 
      }     
      
    } 
    
  }
}     

#Create lists with BAM and PM columns together pulling from   RF_testing_training
merged_1 <- merge(RF_testing[[1]][[1]], RF_testing[[1]][[2]],by="datehour")

clean_merge_1 <- merged_1[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_1["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_1)){
  
  bam_row = bam[which(bam$datehour == clean_merge_1$datehour[i]), ]
  
  clean_merge_1$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_1_with_BAM = clean_merge_1

for (k in length(out_of_range_rows):1){
  
  clean_merge_1_with_BAM = clean_merge_1_with_BAM[-out_of_range_rows[k],]
  
}


#### 2nd Sensor
merged_2 <- merge( RF_testing[[2]][[1]], RF_testing[[2]][[2]],by="datehour")

clean_merge_2 <- merged_2[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_2["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_2)){
  
  bam_row = bam[which(bam$datehour == clean_merge_2$datehour[i]), ]
  
  clean_merge_2$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_2_with_BAM = clean_merge_2

for (k in length(out_of_range_rows):1){
  
  clean_merge_2_with_BAM = clean_merge_2_with_BAM[-out_of_range_rows[k],]
  
}


#### 3rd Sensor
merged_3 <- merge( RF_testing[[3]][[1]], RF_testing[[3]][[2]],by="datehour")

clean_merge_3 <- merged_3[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_3["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_3)){
  
  
  bam_row = bam[which(bam$datehour == clean_merge_3$datehour[i]), ]
  
  clean_merge_3$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_3_with_BAM = clean_merge_3

for (k in length(out_of_range_rows):1){
  
  clean_merge_3_with_BAM = clean_merge_3_with_BAM[-out_of_range_rows[k],]
  
}


#### 4th Sensor
merged_4 <- merge( RF_testing[[4]][[1]], RF_testing[[4]][[2]],by="datehour")

clean_merge_4 <- merged_4[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_4["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_4)){
  
  
  bam_row = bam[which(bam$datehour == clean_merge_4$datehour[i]), ]
  
  clean_merge_4$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_4_with_BAM = clean_merge_4

for (k in length(out_of_range_rows):1){
  
  clean_merge_4_with_BAM = clean_merge_4_with_BAM[-out_of_range_rows[k],]
  
}


#### 5th Sensor
merged_5 <- merge( RF_testing[[5]][[1]], RF_testing[[5]][[2]],by="datehour")

clean_merge_5 <- merged_5[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_5["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_5)){
  
  
  bam_row = bam[which(bam$datehour == clean_merge_5$datehour[i]), ]
  
  clean_merge_5$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_5_with_BAM = clean_merge_5

for (k in length(out_of_range_rows):1){
  
  clean_merge_5_with_BAM = clean_merge_5_with_BAM[-out_of_range_rows[k],]
  
}



#### 6th Sensor
merged_6 <- merge( RF_testing[[6]][[1]], RF_testing[[6]][[2]],by="datehour")

clean_merge_6 <- merged_6[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_6["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_6)){
  
  
  bam_row = bam[which(bam$datehour == clean_merge_6$datehour[i]), ]
  
  clean_merge_6$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_6_with_BAM = clean_merge_6

for (k in length(out_of_range_rows):1){
  
  clean_merge_6_with_BAM = clean_merge_6_with_BAM[-out_of_range_rows[k],]
  
}

clean_merge_6_with_BAM = clean_merge_6_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]



###Comp 7
merged_7 <- merge( RF_testing[[7]][[1]], RF_testing[[7]][[2]],by="datehour")

clean_merge_7 <- merged_7[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_7["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_7)){
  
  bam_row = bam[which(bam$datehour == clean_merge_7$datehour[i]), ]
  
  clean_merge_7$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_7_with_BAM = clean_merge_7

for (k in length(out_of_range_rows):1){
  
  clean_merge_7_with_BAM = clean_merge_7_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 8
merged_8 <- merge( RF_testing[[8]][[1]], RF_testing[[8]][[2]],by="datehour")

clean_merge_8 <- merged_8[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_8["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_8)){
  
  bam_row = bam[which(bam$datehour == clean_merge_8$datehour[i]), ]
  
  clean_merge_8$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_8_with_BAM = clean_merge_8

for (k in length(out_of_range_rows):1){
  
  clean_merge_8_with_BAM = clean_merge_8_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 9
merged_9 <- merge( RF_testing[[9]][[1]], RF_testing[[9]][[2]],by="datehour")

clean_merge_9 <- merged_9[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_9["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_9)){
  
  bam_row = bam[which(bam$datehour == clean_merge_9$datehour[i]), ]
  
  clean_merge_9$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_9_with_BAM = clean_merge_9

for (k in length(out_of_range_rows):1){
  
  clean_merge_9_with_BAM = clean_merge_9_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 10
merged_10 <- merge( RF_testing[[10]][[1]], RF_testing[[10]][[2]],by="datehour")

clean_merge_10 <- merged_10[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_10["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_10)){
  
  bam_row = bam[which(bam$datehour == clean_merge_10$datehour[i]), ]
  
  clean_merge_10$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_10_with_BAM = clean_merge_10

for (k in length(out_of_range_rows):1){
  
  clean_merge_10_with_BAM = clean_merge_10_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 11
merged_11 <- merge( RF_testing[[11]][[1]], RF_testing[[11]][[2]],by="datehour")

clean_merge_11 <- merged_11[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_11["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_11)){
  
  bam_row = bam[which(bam$datehour == clean_merge_11$datehour[i]), ]
  
  clean_merge_11$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_11_with_BAM = clean_merge_11

for (k in length(out_of_range_rows):1){
  
  clean_merge_11_with_BAM = clean_merge_11_with_BAM[-out_of_range_rows[k],]
  
}



###Comp 12
merged_12 <- merge( RF_testing[[12]][[1]], RF_testing[[12]][[2]],by="datehour")

clean_merge_12 <- merged_12[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_12["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_12)){
  
  bam_row = bam[which(bam$datehour == clean_merge_12$datehour[i]), ]
  
  clean_merge_12$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_12_with_BAM = clean_merge_12

for (k in length(out_of_range_rows):1){
  
  clean_merge_12_with_BAM = clean_merge_12_with_BAM[-out_of_range_rows[k],]
  
}


###Comp 13
merged_13 <- merge( RF_testing[[13]][[1]], RF_testing[[13]][[2]],by="datehour")

clean_merge_13 <- merged_13[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_13["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_13)){
  
  bam_row = bam[which(bam$datehour == clean_merge_13$datehour[i]), ]
  
  clean_merge_13$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_13_with_BAM = clean_merge_13

for (k in length(out_of_range_rows):1){
  
  clean_merge_13_with_BAM = clean_merge_13_with_BAM[-out_of_range_rows[k],]
  
}


###Comp 14
merged_14 <- merge( RF_testing[[14]][[1]], RF_testing[[14]][[2]],by="datehour")

clean_merge_14 <- merged_14[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_14["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_14)){
  
  bam_row = bam[which(bam$datehour == clean_merge_14$datehour[i]), ]
  
  clean_merge_14$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_14_with_BAM = clean_merge_14

for (k in length(out_of_range_rows):1){
  
  clean_merge_14_with_BAM = clean_merge_14_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 15
merged_15 <- merge( RF_testing[[15]][[1]], RF_testing[[15]][[2]],by="datehour")

clean_merge_15 <- merged_15[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_15["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_15)){
  
  bam_row = bam[which(bam$datehour == clean_merge_15$datehour[i]), ]
  
  clean_merge_15$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_15_with_BAM = clean_merge_15

for (k in length(out_of_range_rows):1){
  
  clean_merge_15_with_BAM = clean_merge_15_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 16
merged_16 <- merge( RF_testing[[16]][[1]], RF_testing[[16]][[2]],by="datehour")

clean_merge_16 <- merged_16[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_16["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_16)){
  
  bam_row = bam[which(bam$datehour == clean_merge_16$datehour[i]), ]
  
  clean_merge_16$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_16_with_BAM = clean_merge_16

for (k in length(out_of_range_rows):1){
  
  clean_merge_16_with_BAM = clean_merge_16_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 17
merged_17 <- merge( RF_testing[[17]][[1]], RF_testing[[17]][[2]],by="datehour")

clean_merge_17 <- merged_17[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_17["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_17)){
  
  bam_row = bam[which(bam$datehour == clean_merge_17$datehour[i]), ]
  
  clean_merge_17$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_17_with_BAM = clean_merge_17

for (k in length(out_of_range_rows):1){
  
  clean_merge_17_with_BAM = clean_merge_17_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 18
merged_18 <- merge( RF_testing[[18]][[1]], RF_testing[[18]][[2]],by="datehour")

clean_merge_18 <- merged_18[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_18["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_18)){
  
  bam_row = bam[which(bam$datehour == clean_merge_18$datehour[i]), ]
  
  clean_merge_18$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_18_with_BAM = clean_merge_18

for (k in length(out_of_range_rows):1){
  
  clean_merge_18_with_BAM = clean_merge_18_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 19
merged_19 <- merge( RF_testing[[19]][[1]], RF_testing[[19]][[2]],by="datehour")

clean_merge_19 <- merged_19[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_19["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_19)){
  
  bam_row = bam[which(bam$datehour == clean_merge_19$datehour[i]), ]
  
  clean_merge_19$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_19_with_BAM = clean_merge_19

for (k in length(out_of_range_rows):1){
  
  clean_merge_19_with_BAM = clean_merge_19_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 20
merged_20 <- merge( RF_testing[[20]][[1]], RF_testing[[20]][[2]],by="datehour")

clean_merge_20 <- merged_20[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_20["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_20)){
  
  bam_row = bam[which(bam$datehour == clean_merge_20$datehour[i]), ]
  
  clean_merge_20$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_20_with_BAM = clean_merge_20

for (k in length(out_of_range_rows):1){
  
  clean_merge_20_with_BAM = clean_merge_20_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 21
merged_21 <- merge( RF_testing[[21]][[1]], RF_testing[[21]][[2]],by="datehour")

clean_merge_21 <- merged_21[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_21["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_21)){
  
  bam_row = bam[which(bam$datehour == clean_merge_21$datehour[i]), ]
  
  clean_merge_21$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_21_with_BAM = clean_merge_21

for (k in length(out_of_range_rows):1){
  
  clean_merge_21_with_BAM = clean_merge_21_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 22
merged_22 <- merge( RF_testing[[22]][[1]], RF_testing[[22]][[2]],by="datehour")

clean_merge_22 <- merged_22[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_22["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_22)){
  
  bam_row = bam[which(bam$datehour == clean_merge_22$datehour[i]), ]
  
  clean_merge_22$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_22_with_BAM = clean_merge_22

for (k in length(out_of_range_rows):1){
  
  clean_merge_22_with_BAM = clean_merge_22_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 23
merged_23 <- merge( RF_testing[[23]][[1]], RF_testing[[23]][[2]],by="datehour")

clean_merge_23 <- merged_23[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_23["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_23)){
  
  bam_row = bam[which(bam$datehour == clean_merge_23$datehour[i]), ]
  
  clean_merge_23$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_23_with_BAM = clean_merge_23

for (k in length(out_of_range_rows):1){
  
  clean_merge_23_with_BAM = clean_merge_23_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 24
merged_24 <- merge( RF_testing[[24]][[1]], RF_testing[[24]][[2]],by="datehour")

clean_merge_24 <- merged_24[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_24["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_24)){
  
  bam_row = bam[which(bam$datehour == clean_merge_24$datehour[i]), ]
  
  clean_merge_24$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_24_with_BAM = clean_merge_24

for (k in length(out_of_range_rows):1){
  
  clean_merge_24_with_BAM = clean_merge_24_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 25
merged_25 <- merge( RF_testing[[25]][[1]], RF_testing[[25]][[2]],by="datehour")

clean_merge_25 <- merged_25[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_25["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_25)){
  
  bam_row = bam[which(bam$datehour == clean_merge_25$datehour[i]), ]
  
  clean_merge_25$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_25_with_BAM = clean_merge_25

for (k in length(out_of_range_rows):1){
  
  clean_merge_25_with_BAM = clean_merge_25_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 26
merged_26 <- merge( RF_testing[[26]][[1]], RF_testing[[26]][[2]],by="datehour")

clean_merge_26 <- merged_26[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_26["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_26)){
  
  bam_row = bam[which(bam$datehour == clean_merge_26$datehour[i]), ]
  
  clean_merge_26$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_26_with_BAM = clean_merge_26

for (k in length(out_of_range_rows):1){
  
  clean_merge_26_with_BAM = clean_merge_26_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 27
merged_27 <- merge( RF_testing[[27]][[1]], RF_testing[[27]][[2]],by="datehour")

clean_merge_27 <- merged_27[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_27["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_27)){
  
  bam_row = bam[which(bam$datehour == clean_merge_27$datehour[i]), ]
  
  clean_merge_27$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_27_with_BAM = clean_merge_27

for (k in length(out_of_range_rows):1){
  
  clean_merge_27_with_BAM = clean_merge_27_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 28
merged_28 <- merge( RF_testing[[28]][[1]], RF_testing[[28]][[2]],by="datehour")

clean_merge_28 <- merged_28[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_28["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_28)){
  
  bam_row = bam[which(bam$datehour == clean_merge_28$datehour[i]), ]
  
  clean_merge_28$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_28_with_BAM = clean_merge_28

for (k in length(out_of_range_rows):1){
  
  clean_merge_28_with_BAM = clean_merge_28_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 29
merged_29 <- merge( RF_testing[[29]][[1]], RF_testing[[29]][[2]],by="datehour")

clean_merge_29 <- merged_29[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_29["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_29)){
  
  bam_row = bam[which(bam$datehour == clean_merge_29$datehour[i]), ]
  
  clean_merge_29$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_29_with_BAM = clean_merge_29

for (k in length(out_of_range_rows):1){
  
  clean_merge_29_with_BAM = clean_merge_29_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 30
merged_30 <- merge( RF_testing[[30]][[1]], RF_testing[[30]][[2]],by="datehour")

clean_merge_30 <- merged_30[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_30["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_30)){
  
  bam_row = bam[which(bam$datehour == clean_merge_30$datehour[i]), ]
  
  clean_merge_30$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_30_with_BAM = clean_merge_30

for (k in length(out_of_range_rows):1){
  
  clean_merge_30_with_BAM = clean_merge_30_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 31
merged_31 <- merge( RF_testing[[31]][[1]], RF_testing[[31]][[2]],by="datehour")

clean_merge_31 <- merged_31[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_31["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_31)){
  
  bam_row = bam[which(bam$datehour == clean_merge_31$datehour[i]), ]
  
  clean_merge_31$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_31_with_BAM = clean_merge_31

for (k in length(out_of_range_rows):1){
  
  clean_merge_31_with_BAM = clean_merge_31_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 32
merged_32 <- merge( RF_testing[[32]][[1]], RF_testing[[32]][[2]],by="datehour")

clean_merge_32 <- merged_32[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_32["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_32)){
  
  bam_row = bam[which(bam$datehour == clean_merge_32$datehour[i]), ]
  
  clean_merge_32$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_32_with_BAM = clean_merge_32

for (k in length(out_of_range_rows):1){
  
  clean_merge_32_with_BAM = clean_merge_32_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 33
merged_33 <- merge( RF_testing[[33]][[1]], RF_testing[[33]][[2]],by="datehour")

clean_merge_33 <- merged_33[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_33["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_33)){
  
  bam_row = bam[which(bam$datehour == clean_merge_33$datehour[i]), ]
  
  clean_merge_33$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_33_with_BAM = clean_merge_33

for (k in length(out_of_range_rows):1){
  
  clean_merge_33_with_BAM = clean_merge_33_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 34
merged_34 <- merge( RF_testing[[34]][[1]], RF_testing[[34]][[2]],by="datehour")

clean_merge_34 <- merged_34[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_34["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_34)){
  
  bam_row = bam[which(bam$datehour == clean_merge_34$datehour[i]), ]
  
  clean_merge_34$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_34_with_BAM = clean_merge_34

for (k in length(out_of_range_rows):1){
  
  clean_merge_34_with_BAM = clean_merge_34_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 35
merged_35 <- merge( RF_testing[[35]][[1]], RF_testing[[35]][[2]],by="datehour")

clean_merge_35 <- merged_35[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_35["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_35)){
  
  bam_row = bam[which(bam$datehour == clean_merge_35$datehour[i]), ]
  
  clean_merge_35$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_35_with_BAM = clean_merge_35

for (k in length(out_of_range_rows):1){
  
  clean_merge_35_with_BAM = clean_merge_35_with_BAM[-out_of_range_rows[k],]
  
}

###########################Merge All train data

Tc1 = clean_merge_1_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc2 = clean_merge_2_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc3 = clean_merge_3_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc4 = clean_merge_4_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc5 = clean_merge_5_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc6 = clean_merge_6_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc7 = clean_merge_7_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc8 = clean_merge_8_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc9 = clean_merge_9_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc10 = clean_merge_10_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc11 = clean_merge_11_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc12 = clean_merge_12_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc13 = clean_merge_13_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc14 = clean_merge_14_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc15 = clean_merge_15_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc16 = clean_merge_16_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc17 = clean_merge_17_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc18 = clean_merge_18_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc19 = clean_merge_19_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc20 = clean_merge_20_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc21 = clean_merge_21_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc22 = clean_merge_22_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc23 = clean_merge_23_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc24 = clean_merge_24_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc25 = clean_merge_25_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc26 = clean_merge_26_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc27 = clean_merge_27_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc28 = clean_merge_28_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc29 = clean_merge_29_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc30 = clean_merge_30_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc31 = clean_merge_31_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc32 = clean_merge_32_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc33 = clean_merge_33_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc34 = clean_merge_34_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

Tc35 = clean_merge_35_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

# test_data <- rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27, c28, c29,c30,c31,c32,c33,c34,c35)

test_data <- rbind(Tc1,Tc2,Tc3,Tc4,Tc5,Tc6,Tc7,Tc8,Tc9,Tc10,Tc11,Tc12,Tc13,Tc14,Tc15,Tc16,Tc17,Tc18,Tc19,Tc20,Tc21,Tc22,Tc23,Tc24,Tc25,Tc26,c27, Tc28, Tc29,Tc30,Tc31,Tc32,Tc33,Tc34,Tc35)

names(test_data) <- c("PM1A", "PM25A", "PM10A", "PM25ATMA", "PATemperature","PAHumidity", "PM1B", "PM25B", "PM10B", "PM25ATMB", "PAPressure","ConcHRugm3")

#Create lists with BAM and PM columns together pulling from   RF_training_training


merged_1 <- merge(RF_training[[1]][[1]], RF_training[[1]][[2]],by="datehour")

clean_merge_1 <- merged_1[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_1["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_1)){
  
  bam_row = bam[which(bam$datehour == clean_merge_1$datehour[i]), ]
  
  clean_merge_1$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_1_with_BAM = clean_merge_1

for (k in length(out_of_range_rows):1){
  
  clean_merge_1_with_BAM = clean_merge_1_with_BAM[-out_of_range_rows[k],]
  
}


#### 2nd Sensor
merged_2 <- merge( RF_training[[2]][[1]], RF_training[[2]][[2]],by="datehour")

clean_merge_2 <- merged_2[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_2["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_2)){
  
  bam_row = bam[which(bam$datehour == clean_merge_2$datehour[i]), ]
  
  clean_merge_2$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_2_with_BAM = clean_merge_2

for (k in length(out_of_range_rows):1){
  
  clean_merge_2_with_BAM = clean_merge_2_with_BAM[-out_of_range_rows[k],]
  
}


#### 3rd Sensor
merged_3 <- merge( RF_training[[3]][[1]], RF_training[[3]][[2]],by="datehour")

clean_merge_3 <- merged_3[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_3["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_3)){
  
  
  bam_row = bam[which(bam$datehour == clean_merge_3$datehour[i]), ]
  
  clean_merge_3$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_3_with_BAM = clean_merge_3

for (k in length(out_of_range_rows):1){
  
  clean_merge_3_with_BAM = clean_merge_3_with_BAM[-out_of_range_rows[k],]
  
}


#### 4th Sensor
merged_4 <- merge( RF_training[[4]][[1]], RF_training[[4]][[2]],by="datehour")

clean_merge_4 <- merged_4[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_4["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_4)){
  
  
  bam_row = bam[which(bam$datehour == clean_merge_4$datehour[i]), ]
  
  clean_merge_4$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_4_with_BAM = clean_merge_4

for (k in length(out_of_range_rows):1){
  
  clean_merge_4_with_BAM = clean_merge_4_with_BAM[-out_of_range_rows[k],]
  
}


#### 5th Sensor
merged_5 <- merge( RF_training[[5]][[1]], RF_training[[5]][[2]],by="datehour")

clean_merge_5 <- merged_5[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_5["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_5)){
  
  
  bam_row = bam[which(bam$datehour == clean_merge_5$datehour[i]), ]
  
  clean_merge_5$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_5_with_BAM = clean_merge_5

for (k in length(out_of_range_rows):1){
  
  clean_merge_5_with_BAM = clean_merge_5_with_BAM[-out_of_range_rows[k],]
  
}



#### 6th Sensor
merged_6 <- merge( RF_training[[6]][[1]], RF_training[[6]][[2]],by="datehour")

clean_merge_6 <- merged_6[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_6["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_6)){
  
  
  bam_row = bam[which(bam$datehour == clean_merge_6$datehour[i]), ]
  
  clean_merge_6$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_6_with_BAM = clean_merge_6

for (k in length(out_of_range_rows):1){
  
  clean_merge_6_with_BAM = clean_merge_6_with_BAM[-out_of_range_rows[k],]
  
}

clean_merge_6_with_BAM = clean_merge_6_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]



###Comp 7
merged_7 <- merge( RF_training[[7]][[1]], RF_training[[7]][[2]],by="datehour")

clean_merge_7 <- merged_7[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_7["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_7)){
  
  bam_row = bam[which(bam$datehour == clean_merge_7$datehour[i]), ]
  
  clean_merge_7$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_7_with_BAM = clean_merge_7

for (k in length(out_of_range_rows):1){
  
  clean_merge_7_with_BAM = clean_merge_7_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 8
merged_8 <- merge( RF_training[[8]][[1]], RF_training[[8]][[2]],by="datehour")

clean_merge_8 <- merged_8[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_8["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_8)){
  
  bam_row = bam[which(bam$datehour == clean_merge_8$datehour[i]), ]
  
  clean_merge_8$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_8_with_BAM = clean_merge_8

for (k in length(out_of_range_rows):1){
  
  clean_merge_8_with_BAM = clean_merge_8_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 9
merged_9 <- merge( RF_training[[9]][[1]], RF_training[[9]][[2]],by="datehour")

clean_merge_9 <- merged_9[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_9["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_9)){
  
  bam_row = bam[which(bam$datehour == clean_merge_9$datehour[i]), ]
  
  clean_merge_9$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_9_with_BAM = clean_merge_9

for (k in length(out_of_range_rows):1){
  
  clean_merge_9_with_BAM = clean_merge_9_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 10
merged_10 <- merge( RF_training[[10]][[1]], RF_training[[10]][[2]],by="datehour")

clean_merge_10 <- merged_10[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_10["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_10)){
  
  bam_row = bam[which(bam$datehour == clean_merge_10$datehour[i]), ]
  
  clean_merge_10$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_10_with_BAM = clean_merge_10

for (k in length(out_of_range_rows):1){
  
  clean_merge_10_with_BAM = clean_merge_10_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 11
merged_11 <- merge( RF_training[[11]][[1]], RF_training[[11]][[2]],by="datehour")

clean_merge_11 <- merged_11[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_11["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_11)){
  
  bam_row = bam[which(bam$datehour == clean_merge_11$datehour[i]), ]
  
  clean_merge_11$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_11_with_BAM = clean_merge_11

for (k in length(out_of_range_rows):1){
  
  clean_merge_11_with_BAM = clean_merge_11_with_BAM[-out_of_range_rows[k],]
  
}



###Comp 12
merged_12 <- merge( RF_training[[12]][[1]], RF_training[[12]][[2]],by="datehour")

clean_merge_12 <- merged_12[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_12["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_12)){
  
  bam_row = bam[which(bam$datehour == clean_merge_12$datehour[i]), ]
  
  clean_merge_12$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_12_with_BAM = clean_merge_12

# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_12_with_BAM = clean_merge_12_with_BAM[-out_of_range_rows[k],]
#   
# }


###Comp 13
merged_13 <- merge( RF_training[[13]][[1]], RF_training[[13]][[2]],by="datehour")

clean_merge_13 <- merged_13[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_13["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_13)){
  
  bam_row = bam[which(bam$datehour == clean_merge_13$datehour[i]), ]
  
  clean_merge_13$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_13_with_BAM = clean_merge_13

for (k in length(out_of_range_rows):1){
  
  clean_merge_13_with_BAM = clean_merge_13_with_BAM[-out_of_range_rows[k],]
  
}


###Comp 14
merged_14 <- merge( RF_training[[14]][[1]], RF_training[[14]][[2]],by="datehour")

clean_merge_14 <- merged_14[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_14["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_14)){
  
  bam_row = bam[which(bam$datehour == clean_merge_14$datehour[i]), ]
  
  clean_merge_14$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_14_with_BAM = clean_merge_14

for (k in length(out_of_range_rows):1){
  
  clean_merge_14_with_BAM = clean_merge_14_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 15
merged_15 <- merge( RF_training[[15]][[1]], RF_training[[15]][[2]],by="datehour")

clean_merge_15 <- merged_15[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_15["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_15)){
  
  bam_row = bam[which(bam$datehour == clean_merge_15$datehour[i]), ]
  
  clean_merge_15$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_15_with_BAM = clean_merge_15

for (k in length(out_of_range_rows):1){
  
  clean_merge_15_with_BAM = clean_merge_15_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 16
merged_16 <- merge( RF_training[[16]][[1]], RF_training[[16]][[2]],by="datehour")

clean_merge_16 <- merged_16[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_16["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_16)){
  
  bam_row = bam[which(bam$datehour == clean_merge_16$datehour[i]), ]
  
  clean_merge_16$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_16_with_BAM = clean_merge_16

for (k in length(out_of_range_rows):1){
  
  clean_merge_16_with_BAM = clean_merge_16_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 17
merged_17 <- merge( RF_training[[17]][[1]], RF_training[[17]][[2]],by="datehour")

clean_merge_17 <- merged_17[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_17["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_17)){
  
  bam_row = bam[which(bam$datehour == clean_merge_17$datehour[i]), ]
  
  clean_merge_17$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_17_with_BAM = clean_merge_17

for (k in length(out_of_range_rows):1){
  
  clean_merge_17_with_BAM = clean_merge_17_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 18
merged_18 <- merge( RF_training[[18]][[1]], RF_training[[18]][[2]],by="datehour")

clean_merge_18 <- merged_18[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_18["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_18)){
  
  bam_row = bam[which(bam$datehour == clean_merge_18$datehour[i]), ]
  
  clean_merge_18$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_18_with_BAM = clean_merge_18

for (k in length(out_of_range_rows):1){
  
  clean_merge_18_with_BAM = clean_merge_18_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 19
merged_19 <- merge( RF_training[[19]][[1]], RF_training[[19]][[2]],by="datehour")

clean_merge_19 <- merged_19[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_19["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_19)){
  
  bam_row = bam[which(bam$datehour == clean_merge_19$datehour[i]), ]
  
  clean_merge_19$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_19_with_BAM = clean_merge_19

for (k in length(out_of_range_rows):1){
  
  clean_merge_19_with_BAM = clean_merge_19_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 20
merged_20 <- merge( RF_training[[20]][[1]], RF_training[[20]][[2]],by="datehour")

clean_merge_20 <- merged_20[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_20["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_20)){
  
  bam_row = bam[which(bam$datehour == clean_merge_20$datehour[i]), ]
  
  clean_merge_20$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_20_with_BAM = clean_merge_20

for (k in length(out_of_range_rows):1){
  
  clean_merge_20_with_BAM = clean_merge_20_with_BAM[-out_of_range_rows[k],]
  
}

###Comp 21
merged_21 <- merge( RF_training[[21]][[1]], RF_training[[21]][[2]],by="datehour")

clean_merge_21 <- merged_21[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]

clean_merge_21["ConcHR.ug.m3."] <- NA

out_of_range_rows = vector()

for(i in 1:nrow(clean_merge_21)){
  
  bam_row = bam[which(bam$datehour == clean_merge_21$datehour[i]), ]
  
  clean_merge_21$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
  
  if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
    
    out_of_range_rows = append(out_of_range_rows,i)
    
  }
  
}

clean_merge_21_with_BAM = clean_merge_21

for (k in length(out_of_range_rows):1){
  
  clean_merge_21_with_BAM = clean_merge_21_with_BAM[-out_of_range_rows[k],]
  
}

# ###Comp 22
# merged_22 <- merge( RF_training[[22]][[1]], RF_training[[22]][[2]],by="datehour")
# 
# clean_merge_22 <- merged_22[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]
# 
# clean_merge_22["ConcHR.ug.m3."] <- NA
# 
# out_of_range_rows = vector()
# 
# for(i in 1:nrow(clean_merge_22)){
#   
#   bam_row = bam[which(bam$datehour == clean_merge_22$datehour[i]), ]
#   
#   clean_merge_22$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
#   
#   if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
#     
#     out_of_range_rows = append(out_of_range_rows,i)
#     
#   }
#   
# }
# 
# clean_merge_22_with_BAM = clean_merge_22
# 
# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_22_with_BAM = clean_merge_22_with_BAM[-out_of_range_rows[k],]
#   
# }
# 
# ###Comp 23
# merged_23 <- merge( RF_training[[23]][[1]], RF_training[[23]][[2]],by="datehour")
# 
# clean_merge_23 <- merged_23[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]
# 
# clean_merge_23["ConcHR.ug.m3."] <- NA
# 
# out_of_range_rows = vector()
# 
# for(i in 1:nrow(clean_merge_23)){
#   
#   bam_row = bam[which(bam$datehour == clean_merge_23$datehour[i]), ]
#   
#   clean_merge_23$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
#   
#   if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
#     
#     out_of_range_rows = append(out_of_range_rows,i)
#     
#   }
#   
# }
# 
# clean_merge_23_with_BAM = clean_merge_23
# 
# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_23_with_BAM = clean_merge_23_with_BAM[-out_of_range_rows[k],]
#   
# }
# 
# ###Comp 24
# merged_24 <- merge( RF_training[[24]][[1]], RF_training[[24]][[2]],by="datehour")
# 
# clean_merge_24 <- merged_24[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]
# 
# clean_merge_24["ConcHR.ug.m3."] <- NA
# 
# out_of_range_rows = vector()
# 
# for(i in 1:nrow(clean_merge_24)){
#   
#   bam_row = bam[which(bam$datehour == clean_merge_24$datehour[i]), ]
#   
#   clean_merge_24$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
#   
#   if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
#     
#     out_of_range_rows = append(out_of_range_rows,i)
#     
#   }
#   
# }
# 
# clean_merge_24_with_BAM = clean_merge_24
# 
# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_24_with_BAM = clean_merge_24_with_BAM[-out_of_range_rows[k],]
#   
# }
# 
# ###Comp 25
# merged_25 <- merge( RF_training[[25]][[1]], RF_training[[25]][[2]],by="datehour")
# 
# clean_merge_25 <- merged_25[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]
# 
# clean_merge_25["ConcHR.ug.m3."] <- NA
# 
# out_of_range_rows = vector()
# 
# for(i in 1:nrow(clean_merge_25)){
#   
#   bam_row = bam[which(bam$datehour == clean_merge_25$datehour[i]), ]
#   
#   clean_merge_25$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
#   
#   if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
#     
#     out_of_range_rows = append(out_of_range_rows,i)
#     
#   }
#   
# }
# 
# clean_merge_25_with_BAM = clean_merge_25
# 
# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_25_with_BAM = clean_merge_25_with_BAM[-out_of_range_rows[k],]
#   
# }
# 
# ###Comp 26
# merged_26 <- merge( RF_training[[26]][[1]], RF_training[[26]][[2]],by="datehour")
# 
# clean_merge_26 <- merged_26[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]
# 
# clean_merge_26["ConcHR.ug.m3."] <- NA
# 
# out_of_range_rows = vector()
# 
# for(i in 1:nrow(clean_merge_26)){
#   
#   bam_row = bam[which(bam$datehour == clean_merge_26$datehour[i]), ]
#   
#   clean_merge_26$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
#   
#   if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
#     
#     out_of_range_rows = append(out_of_range_rows,i)
#     
#   }
#   
# }
# 
# clean_merge_26_with_BAM = clean_merge_26
# 
# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_26_with_BAM = clean_merge_26_with_BAM[-out_of_range_rows[k],]
#   
# }
# 
# ###Comp 27
# merged_27 <- merge( RF_training[[27]][[1]], RF_training[[27]][[2]],by="datehour")
# 
# clean_merge_27 <- merged_27[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]
# 
# clean_merge_27["ConcHR.ug.m3."] <- NA
# 
# out_of_range_rows = vector()
# 
# for(i in 1:nrow(clean_merge_27)){
#   
#   bam_row = bam[which(bam$datehour == clean_merge_27$datehour[i]), ]
#   
#   clean_merge_27$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
#   
#   if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
#     
#     out_of_range_rows = append(out_of_range_rows,i)
#     
#   }
#   
# }
# 
# clean_merge_27_with_BAM = clean_merge_27
# 
# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_27_with_BAM = clean_merge_27_with_BAM[-out_of_range_rows[k],]
#   
# }
# 
# ###Comp 28
# merged_28 <- merge( RF_training[[28]][[1]], RF_training[[28]][[2]],by="datehour")
# 
# clean_merge_28 <- merged_28[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]
# 
# clean_merge_28["ConcHR.ug.m3."] <- NA
# 
# out_of_range_rows = vector()
# 
# for(i in 1:nrow(clean_merge_28)){
#   
#   bam_row = bam[which(bam$datehour == clean_merge_28$datehour[i]), ]
#   
#   clean_merge_28$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
#   
#   if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
#     
#     out_of_range_rows = append(out_of_range_rows,i)
#     
#   }
#   
# }
# 
# clean_merge_28_with_BAM = clean_merge_28
# 
# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_28_with_BAM = clean_merge_28_with_BAM[-out_of_range_rows[k],]
#   
# }
# 
# ###Comp 29
# merged_29 <- merge( RF_training[[29]][[1]], RF_training[[29]][[2]],by="datehour")
# 
# clean_merge_29 <- merged_29[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]
# 
# clean_merge_29["ConcHR.ug.m3."] <- NA
# 
# out_of_range_rows = vector()
# 
# for(i in 1:nrow(clean_merge_29)){
#   
#   bam_row = bam[which(bam$datehour == clean_merge_29$datehour[i]), ]
#   
#   clean_merge_29$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
#   
#   if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
#     
#     out_of_range_rows = append(out_of_range_rows,i)
#     
#   }
#   
# }
# 
# clean_merge_29_with_BAM = clean_merge_29
# 
# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_29_with_BAM = clean_merge_29_with_BAM[-out_of_range_rows[k],]
#   
# }
# 
# ###Comp 30
# merged_30 <- merge( RF_training[[30]][[1]], RF_training[[30]][[2]],by="datehour")
# 
# clean_merge_30 <- merged_30[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]
# 
# clean_merge_30["ConcHR.ug.m3."] <- NA
# 
# out_of_range_rows = vector()
# 
# for(i in 1:nrow(clean_merge_30)){
#   
#   bam_row = bam[which(bam$datehour == clean_merge_30$datehour[i]), ]
#   
#   clean_merge_30$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
#   
#   if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
#     
#     out_of_range_rows = append(out_of_range_rows,i)
#     
#   }
#   
# }
# 
# clean_merge_30_with_BAM = clean_merge_30
# 
# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_30_with_BAM = clean_merge_30_with_BAM[-out_of_range_rows[k],]
#   
# }
# 
# ###Comp 31
# merged_31 <- merge( RF_training[[31]][[1]], RF_training[[31]][[2]],by="datehour")
# 
# clean_merge_31 <- merged_31[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]
# 
# clean_merge_31["ConcHR.ug.m3."] <- NA
# 
# out_of_range_rows = vector()
# 
# for(i in 1:nrow(clean_merge_31)){
#   
#   bam_row = bam[which(bam$datehour == clean_merge_31$datehour[i]), ]
#   
#   clean_merge_31$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
#   
#   if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
#     
#     out_of_range_rows = append(out_of_range_rows,i)
#     
#   }
#   
# }
# 
# clean_merge_31_with_BAM = clean_merge_31
# 
# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_31_with_BAM = clean_merge_31_with_BAM[-out_of_range_rows[k],]
#   
# }
# 
# ###Comp 32
# merged_32 <- merge( RF_training[[32]][[1]], RF_training[[32]][[2]],by="datehour")
# 
# clean_merge_32 <- merged_32[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]
# 
# clean_merge_32["ConcHR.ug.m3."] <- NA
# 
# out_of_range_rows = vector()
# 
# for(i in 1:nrow(clean_merge_32)){
#   
#   bam_row = bam[which(bam$datehour == clean_merge_32$datehour[i]), ]
#   
#   clean_merge_32$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
#   
#   if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
#     
#     out_of_range_rows = append(out_of_range_rows,i)
#     
#   }
#   
# }
# 
# clean_merge_32_with_BAM = clean_merge_32
# 
# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_32_with_BAM = clean_merge_32_with_BAM[-out_of_range_rows[k],]
#   
# }
# 
# ###Comp 33
# merged_33 <- merge( RF_training[[33]][[1]], RF_training[[33]][[2]],by="datehour")
# 
# clean_merge_33 <- merged_33[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]
# 
# clean_merge_33["ConcHR.ug.m3."] <- NA
# 
# out_of_range_rows = vector()
# 
# for(i in 1:nrow(clean_merge_33)){
#   
#   bam_row = bam[which(bam$datehour == clean_merge_33$datehour[i]), ]
#   
#   clean_merge_33$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
#   
#   if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
#     
#     out_of_range_rows = append(out_of_range_rows,i)
#     
#   }
#   
# }
# 
# clean_merge_33_with_BAM = clean_merge_33
# 
# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_33_with_BAM = clean_merge_33_with_BAM[-out_of_range_rows[k],]
#   
# }
# 
# ###Comp 34
# merged_34 <- merge( RF_training[[34]][[1]], RF_training[[34]][[2]],by="datehour")
# 
# clean_merge_34 <- merged_34[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]
# 
# clean_merge_34["ConcHR.ug.m3."] <- NA
# 
# out_of_range_rows = vector()
# 
# for(i in 1:nrow(clean_merge_34)){
#   
#   bam_row = bam[which(bam$datehour == clean_merge_34$datehour[i]), ]
#   
#   clean_merge_34$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
#   
#   if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
#     
#     out_of_range_rows = append(out_of_range_rows,i)
#     
#   }
#   
# }
# 
# clean_merge_34_with_BAM = clean_merge_34
# 
# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_34_with_BAM = clean_merge_34_with_BAM[-out_of_range_rows[k],]
#   
# }
# 
# ###Comp 35
# merged_35 <- merge( RF_training[[35]][[1]], RF_training[[35]][[2]],by="datehour")
# 
# clean_merge_35 <- merged_35[,c("datehour", "PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa")]
# 
# clean_merge_35["ConcHR.ug.m3."] <- NA
# 
# out_of_range_rows = vector()
# 
# for(i in 1:nrow(clean_merge_35)){
#   
#   bam_row = bam[which(bam$datehour == clean_merge_35$datehour[i]), ]
#   
#   clean_merge_35$`ConcHR.ug.m3.`[i] = bam_row$`ConcHR(ug/m3)`
#   
#   if(bam_row$`ConcHR(ug/m3)` == 99999.0 | bam_row$`ConcHR(ug/m3)` < 0 ){
#     
#     out_of_range_rows = append(out_of_range_rows,i)
#     
#   }
#   
# }
# 
# clean_merge_35_with_BAM = clean_merge_35
# 
# for (k in length(out_of_range_rows):1){
#   
#   clean_merge_35_with_BAM = clean_merge_35_with_BAM[-out_of_range_rows[k],]
#   
# }

###########################Merge All train data

c1 = clean_merge_1_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c2 = clean_merge_2_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c3 = clean_merge_3_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c4 = clean_merge_4_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c5 = clean_merge_5_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c6 = clean_merge_6_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c7 = clean_merge_7_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c8 = clean_merge_8_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c9 = clean_merge_9_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c10 = clean_merge_10_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c11 = clean_merge_11_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c12 = clean_merge_12_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c13 = clean_merge_13_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c14 = clean_merge_14_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c15 = clean_merge_15_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c16 = clean_merge_16_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c17 = clean_merge_17_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c18 = clean_merge_18_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c19 = clean_merge_19_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c20 = clean_merge_20_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

c21 = clean_merge_21_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

# c22 = clean_merge_22_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]

# c23 = clean_merge_23_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# c24 = clean_merge_24_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# c25 = clean_merge_25_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# c26 = clean_merge_26_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# c27 = clean_merge_27_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# c28 = clean_merge_28_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# c29 = clean_merge_29_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# c30 = clean_merge_30_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# c31 = clean_merge_31_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# c32 = clean_merge_32_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# c33 = clean_merge_33_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# c34 = clean_merge_34_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# c35 = clean_merge_35_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "Temperature_F", "Humidity_%", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 



# clean_merge_1_with_BAM = clean_merge_1_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# clean_merge_2_with_BAM = clean_merge_2_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# clean_merge_3_with_BAM = clean_merge_3_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# clean_merge_4_with_BAM = clean_merge_4_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]
# 
# clean_merge_5_with_BAM = clean_merge_5_with_BAM[,c("PM1.0_CF1_ug/m3.x", "PM2.5_CF1_ug/m3.x", "PM10.0_CF1_ug/m3.x", "PM2.5_ATM_ug/m3.x", "PM1.0_CF1_ug/m3.y", "PM2.5_CF1_ug/m3.y", "PM10.0_CF1_ug/m3.y", "PM2.5_ATM_ug/m3.y", "Pressure_hpa", "ConcHR.ug.m3.")]


#train_data <- rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)

# ,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21) #c22, c23,c24,c25,c26,c27, c28, c29,c30,c31,c32,c33,c34,c35)

train_data <- rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21)
names(train_data) <- c("PM1A", "PM25A", "PM10A", "PM25ATMA", "PATemperature","PAHumidity", "PM1B", "PM25B", "PM10B", "PM25ATMB", "PAPressure","ConcHRugm3")
View(train_data)

##Create randomForest
randomForestModel <- randomForest(ConcHRugm3 ~ ., data = train_data, ntree=500, importance=TRUE)

##Review attributes of model
print(randomForestModel)
attributes(randomForestModel)

##Predictions
p1 = predict(randomForestModel, newdata=test_data)

##Number of trees providing lowest error rate, 
which.min(randomForestModel$mse)

##which is 492 trees providing an Conc error of 8.728 off  
sqrt(randomForestModel$mse[which.min(randomForestModel$mse)])


sqrt(randomForestModel$mse[which.min(randomForestModel$mse)])

##RMSE Calculation
#rmse <- function(p1, test_data$ConcHRugm3){ 
#  res <- sqrt(mean((p1-test_data$ConcHRugm3)^2)) return(res) 
#}


##Accuracy
acc <- sum(diag(p1))/sum(p1)
print(acc)
#Which is 0.7333 pretty remarkable considering individual 
#trees accuracy, for which Forest Tree does. It increases 
#the model accuracy according with the test explanation. 

##Final Error Catagory Analysis (0-.5)
# Find errors in predictions
test_data$Absolute_Error = abs(p1 - test_data$ConcHRugm3)
test_data$RF_correction_Error_raw = (p1 / test_data$ConcHRugm3)
test_data$RF_correction_Error_cat = cut(test_data$RF_correction_Error_raw, c(-Inf, .05, .15, .25, .35, .45, .55, .65, .75, .85, 
                                                                             .95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, Inf), labels = c("<.05", ".05 - .15", 
                                                                                                                                                               ".15 - .25", ".25 - .35", ".35 - .45", ".45 - .55", ".55 - .65", ".65 - .75", ".75 - .85",
                                                                                                                                                               ".85 - .95", ".95 - 1.05", "1.05 - 1.15", "1.15 - 1.25", "1.25 - 1.35", "1.35 - 1.45", "1.45 - 1.55", "1.55 - 1.65",
                                                                                                                                                               "1.65 - 1.75", "1.75 - 1.85", "1.85 - 1.95", ">1.95"))

View(test_data)

Error_raw_A_B = c((test_data$PM25A  / test_data$ConcHRugm3), (test_data$PM25B  / test_data$ConcHRugm3))
Error_cat_A_B = cut(Error_raw_A_B, c(-Inf, .05, .15, .25, .35, .45, .55, .65, .75, .85, 
                                     .95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, Inf), labels = c("<.05", ".05 - .15", 
                                                                                                                       ".15 - .25", ".25 - .35", ".35 - .45", ".45 - .55", ".55 - .65", ".65 - .75", ".75 - .85",
                                                                                                                       ".85 - .95", ".95 - 1.05", "1.05 - 1.15", "1.15 - 1.25", "1.25 - 1.35", "1.35 - 1.45", "1.45 - 1.55", "1.55 - 1.65",
                                                                                                                       "1.65 - 1.75", "1.75 - 1.85", "1.85 - 1.95", ">1.95"))

View(test_data)

# Compare original accuracy to HMM accuracy by determing how many data points for each fall into the .95 - 1.05 error category
orig_accuracy = sum(Error_cat_A_B == ".95 - 1.05", na.rm = T)
RF_accuracy = sum(test_data$RF_correction_Error_cat == ".95 - 1.05", na.rm = T)
stats_error = test_data$Absolute_Error

print(paste("Original Accuracy (# of values falling in .95 - 1.05 category)", orig_accuracy/(nrow(test_data)*2), sep = " "))
print(paste("RF Accuracy (# of values falling in .95 - 1.05 category)", RF_accuracy/nrow(test_data), sep = " "))
print(paste("Median Absolute Error", median(test_data$Absolute_Error, na.rm = T), sep = " "))


#get Pearson Corr. BAM vs. A sensors
corrAvBAM <- vector()

for (sen_idx in 1:length(csvs)) {
  tempBam = bam[(bam$datehour %in% csvs[[sen_idx]][[1]]$datehour) & (bam$`ConcHR(ug/m3)` != 99999),]
  tempCsv = csvs[[sen_idx]][[1]][csvs[[sen_idx]][[1]]$datehour %in% tempBam$datehour,]
  print(paste(sen_idx, ": ", cor(tempCsv$`PM2.5_CF1_ug/m3`, tempBam$`ConcHR(ug/m3)`)))
  corrAvBAM[sen_idx] = cor(tempCsv$`PM2.5_CF1_ug/m3`, tempBam$`ConcHR(ug/m3)`)
  
}


#barplot Pearson Corr
barplot(corrAvBAM, main = "Pearson Correlation BAM vs PurpleAir A sensors",
        xlab = "PurpleAir sensors (index and lettering)",
        ylab = "Pearson Correlation",
        names.arg = c("01_A", "02_A", "03_A", "04_A", "05_A", "06_A", "07_A", "08_A", "09_A", "10_A", "11_A", "12_A", "13_A", "14_A", "15_A", "16_A", "17_A", "18_A", "19_A", "20_A", "21_A", "22_A", "23_A", "24_A", "25_A", "26_A", "27_A", "28_A", "29_A", "30_A", "31_A", "32_A", "33_A", "34_A", "35_A", "36_A"),
        col = "darkred",
        horiz = FALSE,
        width = 1,
        ylim = c(0.0 , 1.0)
)

#get Pearson Corr. BAM vs. uncorrected time segment
corrAvBAM <- vector()

for (sen_idx in 1:length(markov_testing)) {
  tempBam = bam[(bam$datehour %in% markov_testing[[sen_idx]]$datehour) & (bam$`ConcHR(ug/m3)` != 99999),]
  tempCsv = markov_testing[[sen_idx]][markov_testing[[sen_idx]]$datehour %in% tempBam$datehour,]
  print(paste(sen_idx, ": ", cor(tempCsv$`PM2.5_CF1_ug/m3`, tempBam$`ConcHR(ug/m3)`)))
  corrAvBAM[sen_idx] = cor(tempCsv$`PM2.5_CF1_ug/m3`, tempBam$`ConcHR(ug/m3)`)
  
}


#barplot Pearson Corr
barplot(corrAvBAM, main = "Pearson Correlation BAM vs Uncorrected PurpleAir Continuous Time Segments",
        xlab = "PurpleAir sensors (index)",
        ylab = "Pearson Correlation",
        col = "darkred",
        horiz = FALSE,
        width = 1,
        ylim = c(-0.5 , 0.5)
)

#get Pearson Corr. BAM vs. corrected time segment
corrAvBAM <- vector()

for (sen_idx in 1:length(markov_testing)) {
  tempBam = bam[(bam$datehour %in% markov_testing[[sen_idx]]$datehour) & (bam$`ConcHR(ug/m3)` != 99999),]
  tempCsv = markov_testing[[sen_idx]][markov_testing[[sen_idx]]$datehour %in% tempBam$datehour,]
  tempBam = tempBam[!is.na(tempCsv$HMM_correction),]
  tempCsv = tempCsv[!is.na(tempCsv$HMM_correction),]
  print(paste(sen_idx, ": ", cor(tempCsv$HMM_correction, tempBam$`ConcHR(ug/m3)`)))
  corrAvBAM[sen_idx] = cor(tempCsv$HMM_correction, tempBam$`ConcHR(ug/m3)`)
  
}


#barplot Pearson Corr
barplot(corrAvBAM, main = "Pearson Correlation BAM vs Corrected PurpleAir Continuous Time Segments",
        xlab = "PurpleAir sensors (index)",
        ylab = "Pearson Correlation",
        col = "darkred",
        horiz = FALSE,
        width = 1,
        ylim = c(-0.5 , 0.5)
)

##Figure7
y_axis = randomForestModel$mse
x_axis = seq(1,randomForestModel$ntree,1)
plot(x_axis, y_axis, main="RandomForest Model Error Rate: Error vs. Number of Trees", ylab="Error", xlab= "Number of Trees")

cor(randomForestModel$predicted, train_data$ConcHRugm3)
##Figure8 
plot(randomForestModel$predicted, train_data$ConcHRugm3, main="RandomForest Model: Corrected PM2.5 concentrations vs. BAM concentrations", ylab="Corrected BAM concentrations", xlab= "BAM concentrations") 
legend("topleft", title="Correlation",  c("0.504702"))
legend(cor(randomForestModel$predicted, train_data$ConcHRugm3), legend=c("Line 1", "Line 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8)

##Figure9
#plot(train_data$ConcHRugm3, main="RandomForest Model: Train_Data BAM concentrations vs. Number of BAM measurements", ylab="Train_Data BAM concentrations", xlab= "Number of BAM measurements")


##Figure10
varImpPlot(randomForestModel)

