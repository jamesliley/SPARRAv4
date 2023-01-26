#' get_episodes_within_time
#' This function returns a binary vector the same length as input times,
#' which is 1 iff time_x is up to max_days before a time_y for which condition is true
#' @param times an input vector of episode times in UNIX format (!!!), we expect it to be ordered
#' @param condition a binary vector that signifies times for which an event is true
#' @param max_days is the maximum time window (in days) considered between times such that an event at the later time affects the earlier time 
#' 
#' @return a binary vector the same length as input times, 
#' which is 1 iff time_x is up to max_days before a time_y for which condition is true
#' 
#' @export
get_episodes_within_time <- function(times, condition, max_days=365){

  # Convert days to seconds (so that it is consistent with UNIX time)
  max_seconds = max_days*24*60*60
  
  out = vector(mode="logical", length=length(times))
  i = length(times)
  cur_state = FALSE
  last_true_time = times[i]
  
  while (i>=1){
    # Set the output to cur_state (with certain conditions)
    if (cur_state == TRUE){
      #timediff = difftime(last_true_time, times[i], units = "days")
      timediff = last_true_time - times[i]
      # Check if enough time has passed since last_true_time to set cur_state to 0
      if (timediff > max_seconds){
        cur_state = FALSE
      }
      
      # Set the output
      if (timediff>0){
        out[i] = cur_state 
      } else {
        # Co-occurring events, can't influcence eachother
        # Set the output to the same as the other event
        out[i] = out[i+1]
      }
    } else {
      out[i] = FALSE
    }
    
    
    # Update cur_state
    if (condition[i]){
      cur_state = TRUE
      last_true_time = times[i]
    } 
    
    # Step to previous time
    i = i-1
    
  }
  
  out
}