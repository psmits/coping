read_one_stan_csv <- function(csvfile) {
  if (length(csvfile) != 1) 
    stop("'csvfile' must be of length 1")
  if (!file.exists(csvfile))
    stop("'csvfile' does not exist on the disk")

  mark <- 0L
  fields <- character()
  while(length(fields) == 0) {
    mark <- mark + 1L
    fields <- scan(csvfile, what = character(), sep = ",", 
                   comment.char = "#", nlines = mark, quiet = TRUE)
  }

  comments <- scan(csvfile, what = character(), sep = "\n", comment.char = "", 
                   nlines = mark + 2L, quiet = TRUE)
  comments <- gsub("#", "", comments, fixed = TRUE)
  comments <- gsub("(Default)", "", comments, fixed = TRUE)
  comments <- grep("=", comments, fixed = TRUE, value = TRUE)
  comments <- strsplit(comments, split = "=", fixed = TRUE)
  comments <- lapply(comments, FUN = trimws)
  comments <- sapply(comments, FUN = function(x) {
                     y <- x[2]
                     names(y) <- x[1]
                     return(y)
                   })

  method <- comments["algorithm"]
  if (method %in% c("meanfield", "fullrank")) {
    draws <- scan(csvfile, what = double(), sep = ",", comment.char = "", 
                  quiet = TRUE, skip = mark + 2L, 
                  nlines = mark + as.integer(comments["output_samples"]) + 3L)
    timings <- NULL
  }
  else { # sampling
    iter <- as.integer(comments["iter"])
    draws <- scan(csvfile, what = double(), sep = ",", comment.char = "", 
                  quiet = TRUE, skip = mark, nlines = mark + iter)
    timings <- scan(csvfile, what = character(), sep = "\n", comment.char = "", 
                    quiet = TRUE, skip = mark + iter)
  }
  draws <- matrix(draws, ncol = length(fields), byrow = TRUE)
  colnames(draws) <- fields
  draws <- as.data.frame(draws)

  attributes(draws)$comments <- comments
  attributes(draws)$timings <- timings
  return(draws)
}

