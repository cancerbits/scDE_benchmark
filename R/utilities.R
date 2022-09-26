# example of an R file with a helper function

# pretty print a time difference between two proc.time() calls
time_diff <- function(start_time, end_time = NULL) {
  if (is.null(end_time)) {
    end_time <- proc.time()
  }
  dt_cpu <- lubridate::make_difftime(num = sum(end_time[c('user.self', 'sys.self')] - start_time[c('user.self', 'sys.self')]))
  dt_elapsed <- lubridate::make_difftime(num = end_time['elapsed'] - start_time['elapsed'])
  
  sprintf('Elapsed time: %1.2f %s; CPU time: %1.2f %s', 
          dt_elapsed, attr(x = dt_elapsed, which = 'units'),
          dt_cpu, attr(x = dt_cpu, which = 'units'))
}
