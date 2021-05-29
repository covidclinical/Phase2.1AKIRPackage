#' Generates AKI KDIGO grades using previous serum creatinine values
#' @param x data.table contaning serum creatinine values, the baseline in past 90 days and past 48h
#' @noRd

aki_kdigo_grade <- function(x) {
  creat = as.numeric(x[4])
  baseline_90d = as.numeric(x[5])
  baseline_48h = as.numeric(x[6])
  grade = 0
  diff = creat - baseline_48h
  ratio = round(creat/baseline_90d,2)
  if(diff > 0.3 || ratio >= 1.5) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 4 || ratio >= 3) {
    grade = 3
  }
  grade
}

#' Generates AKI KDIGO grades using future serum creatinine values
#' @param x data.table contaning serum creatinine values, the baseline in the next 7 days and next 48h
#' @noRd

aki_kdigo_grade_retro <- function(x) {
  # Instead of using the pure KDIGO definition, this looks at future Cr values, and 
  # determines whether the current value, in retrospect, represents an episode of AKI
  creat = as.numeric(x[4])
  baseline_7d = as.numeric(x[7])
  baseline_48h = as.numeric(x[8])
  grade = 0
  diff = creat - baseline_48h
  ratio = round(creat/baseline_7d,2)
  if(diff > 0.3 || ratio >= 1.5) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 4 || ratio >= 3) {
    grade = 3
  }
  grade
}

#' Generates AKI KDIGO grades but on serum creatinine values 7 days later using previous serum creatinine values
#' @param x data.table contaning serum creatinine values, the baseline in past 90 days and past 48h
#' @noRd

akd_grade_7d <- function(x) {
  creat = as.numeric(x[4])
  baseline_90d = as.numeric(x[5])
  baseline_48h = as.numeric(x[6])
  baseline_7d_retro = as.numeric(x[7])
  baseline_48h_retro = as.numeric(x[8])
  baseline = min(baseline_90d,baseline_48h,baseline_7d_retro,baseline_48h_retro)
  cr_7d = as.numeric(x[9])
  grade = 0
  ratio = round(cr_7d/baseline,2)
  diff = cr_7d - baseline
  # We will code grade B/C as 0.5
  if(ratio > 1.25) {
    grade = 0.5
  }
  if(ratio >= 1.5 & ratio < 2) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 4 || ratio >= 3) {
    grade = 3
  }
  grade
}

#' Generates AKI KDIGO grades but on serum creatinine values 90 days later using previous serum creatinine values
#' @param x data.table contaning serum creatinine values, the baseline in past 90 days and past 48h
#' @noRd

akd_grade_90d <- function(x) {
  creat = as.numeric(x[4])
  baseline_90d = as.numeric(x[5])
  baseline_48h = as.numeric(x[6])
  baseline_7d_retro = as.numeric(x[7])
  baseline_48h_retro = as.numeric(x[8])
  baseline = min(baseline_90d,baseline_48h,baseline_7d_retro,baseline_48h_retro)
  cr_90d = as.numeric(x[10])
  grade = 0
  ratio = round(cr_90d/baseline,2)
  diff = cr_90d - baseline
  # We will code grade B/C as 0.5
  if(ratio > 1.25) {
    grade = 0.5
  }
  if(ratio >= 1.5 & ratio < 2) {
    grade = 1
  }
  if(ratio >= 2 & ratio < 3) {
    grade = 2
  }
  if(diff > 4 || ratio >= 3) {
    grade = 3
  }
  grade
}

#' Gives the position of the lowest serum creatinine in a specified window
#' @param cr Vector containing serum creatinine values
#' @param day Vector containing time points of respective serum creatinine values
#' @param lag Specified whether to look at retrospective or future data. Default = TRUE
#' @param gap Specifies rolling window length to use. Default = 7 days
#' @noRd

pos_min <- function(cr,day,lag=TRUE,gap=7) {
  len = length(cr)
  day_pos = day
  for(i in 1:len) {
    if(lag) {
      j = max(1,i-gap)
      pos = j-1+which.min(cr[j:i])
      day_pos[i] = day[pos]
    } else {
      j = min(i+gap,len)
      pos = i-1+which.min(cr[i:j])
      day_pos[i] = day[pos]
    }
  }
  day_pos
}

#' Boolean function which returns whether a certain value is minimum/maximum
#' @param x Vector containing serial creatinine values
#' @param partial Specifies whether to count the start/end of the array as a min/max value
#' @param decreasing Specified whether to find min (TRUE) or max (FALSE)
#' @noRd

which.peaks <- function(x,partial=TRUE,decreasing=FALSE) {
  if(decreasing) {
    if(partial) {
      which(diff(c(FALSE,diff(x) > 0,TRUE)) > 0)
    } else {
      which(diff(diff(x)>0)>0) + 1
    }
  } else {
    if(partial) {
      which(diff(c(TRUE,diff(x) >= 0,FALSE)) < 0)
    } else {
      which(diff(diff(x)>=0)<0) + 1
    }
  }
}

#' Returns the time point at which the normalised serum Cr falls below a certain value
#' @param ratio Vector containing normalised serum Cr values
#' @param time_from_peak Day of interest, expressed as time from peak
#' @param target Threshold of interest, default = 1.25
#' @noRd
get_day <- function(ratio,time_from_peak,target=1.25) {
  index = purrr::detect_index(ratio,function(x) x <= 1.25)
  if(index > 0) {
    day = time_from_peak[index]
  } else {
    day = time_from_peak[length(time_from_peak)]
  }
  day
}

#' Generates MELD score
#' @param bil serum bilirubin in mg/dL
#' @param inr INR
#' @param sCr serum creatinine in mg/dL
#' @noRd
meld_score <- function(bil,inr,sCr,Na = 137) {
  # Unable to implement corrections for sodium at present as not available in the 4CE labs
  # However, still hard-coded in, in case this changes in future
  bil_corr <- ifelse((bil < 1 | is.na(bil)),1,bil)
  inr_corr <- ifelse((inr < 1 | is.na(inr)),1,inr)
  sCr_corr <- ifelse((sCr < 1 | is.na (sCr)),1,min(4,sCr))
  na_corr <- ifelse(is.na(Na),137,Na)
  meld_init <- 3.78 * log(bil_corr) + 11.2 * log(inr_corr) + 9.57 * log(sCr_corr) + 6.43
  meld_corr <- meld_init + 1.32 * (137 - na_corr) - 0.033 * meld_init * (137 - na_corr)
  meld_corr <- round(meld_corr)
  meld_corr
}

#' Display function to show mean(SD) for continuous variables in the table1 demographics table
#' @noRd
my.render.cont <- function(x) {
  with(table1::stats.apply.rounding(table1::stats.default(x), digits=2),c("","Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}

#' Display function to show N (%) for continuous variables in the table1 demographics table
#' @noRd
my.render.cat <- function(x) {
  c("", sapply(table1::stats.default(x), function(y) with(y, sprintf("%d (%0.0f %%)", FREQ, PCT))))
}

#' Computes t-test statistics from summary statistics
#' @param m1 sample 1 mean
#' @param m2 sample 2 mean
#' @param s1 sample 1 SD
#' @param s2 sample 2 SD
#' @param n1 sample 1 size
#' @param n2 sample 2 size
#' @param m0 the null value for the difference in means to be tested for. Default is 0. 
#' @paramequal.variance whether or not to assume equal variance. Default is FALSE. 
#' @noRd
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}