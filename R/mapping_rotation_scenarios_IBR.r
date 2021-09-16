
R --vanilla

if (!requireNamespace("data.table")) install.packages("data.table")
if (!requireNamespace("mobr")) install.packages("mobr")

library(data.table)
library(mobr)

## home function
# iNEXT coverage
Chat.Ind <- function(x, m){
  x <- x[x>0]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  Sub <- function(m){
    #if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m < n) {
      xx <- x[(n-x)>=m]
      out <- 1-sum(xx / n * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }
    if(m == n) out <- 1-f1/n*A
    if(m > n) out <- 1-f1/n*A^(m-n+1)
    out
  }
  sapply(m, Sub)
}

smpl_space_y <- function(j, a, max_effort, samp_method){
                  x <- data.table(reshape2::dcast(a[year == j, .N, by = .(parcel_id, crop_name)], parcel_id ~ crop_name, fun.agg = function(x) sum(x), value.var = "N"))
                  y_rar <- as.numeric(rarefaction(as.matrix(x[, safety := 0][, -c("parcel_id")]), method = samp_method, effort = c(1:max_effort), extrapolate = TRUE, quiet_mode = TRUE))
              }

smpl_time_y <- function(j, a, max_effort, samp_method){
                  x <- data.table(reshape2::dcast(a[parcel_id == j, .N, by = .(year, crop_name)], year ~ crop_name, fun.agg = function(x) sum(x), value.var = "N"))
                  f_rar <- as.numeric(rarefaction(as.matrix(x[, safety := 0][, -c("year")]), method = samp_method, effort = c(1:max_effort), extrapolate = TRUE, quiet_mode = TRUE))
              }

cov_time_y <- function(j, a, max_effort){
                    x <- data.table(reshape2::dcast(a[parcel_id == j, .N, by = .(year, crop_name)], year ~ crop_name, fun.agg = function(x) sum(x), value.var = "N"))
                    f_rar <- as.numeric(Chat.Ind(colSums(as.matrix(x[, safety := 0][, -c("year")])), c(1:max_effort)))
              }

trans.matrix <- function(X, prob=T){
              tt <- table( c(X[,-ncol(X)]), c(X[,-1]) )
              if(prob) tt <- tt / rowSums(tt)
              tt
            }


## Scenarios (6)
## =============
## MONOCULTURE -space and time
## MONOCULTURE -time but diverse in space
## Synchronized rotation - Temporal synchrony
## Asynchrony in rotation
## Diverse rotations but synchronized
## Diverse rotation and asynchony

mono <- data.table(parcel_id = c(rep(c(1,2,3,4), each = 4)),
            year = c(rep(c(1,2,3,4),c(4))),
            crop_name = c(rep("a", 16)))
            
mono_div <- data.table(parcel_id = c(rep(c(1,2,3,4), each = 4)),
            year = c(rep(c(1,2,3,4),c(4))),
            crop_name = c(rep(c("a", "a", "b", "c"), each = 4)))

rot_sync <- data.table(parcel_id = c(rep(c(1,2,3,4), each = 4)),
            year = c(rep(c(1,2,3,4),c(4))),
            crop_name = c("a", "a", "b", "c",
                          "a", "a", "b", "c", 
                          "a", "a", "b", "c",
                          "a", "a", "b", "c"))

rot_desync <- data.table(parcel_id = c(rep(c(1,2,3,4), each = 4)),
            year = c(rep(c(1,2,3,4),c(4))),
            crop_name = c("a", "a", "b", "c",
                          "c", "a", "a", "b", 
                          "b", "c", "a", "a",
                          "a", "b", "c", "a"))

div_rot_sync <- data.table(parcel_id = c(rep(c(1,2,3,4), each = 4)),
            year = c(rep(c(1,2,3,4),c(4))),
            crop_name = c("a", "a", "b", "c",
                          "e", "e", "f", "g", 
                          "a", "a", "b", "c",
                          "e", "e", "f", "g"))

div_rot_desync <- data.table(parcel_id = c(rep(c(1,2,3,4), each = 4)),
            year = c(rep(c(1,2,3,4),c(4))),
            crop_name = c("a", "a", "b", "c",
                          "g", "e", "e", "f", 
                          "b", "c", "a", "a",
                          "f", "g", "e", "e"))

pdf_name <- c("mono.pdf",               ## MONOCULTURE -space and time
              "mono_div.pdf",           ## MONOCULTURE -time but diverse in space
              "rot_sync.pdf",           ## Synchronized rotation - Temporal synchrony
              "rot_desync.pdf",         ## Asynchrony in rotation
              "div_rot_sync.pdf",       ## Diverse rotations but synchronized
              "div_rot_desync.pdf")     ## Diverse rotation and asynchony

scenario_data <- list(mono, mono_div, rot_sync, rot_desync, div_rot_sync, div_rot_desync)

max_effort <- 16

for(i in seq_along(pdf_name)){

  a_t <- scenario_data[[i]]

  x <- data.table(reshape2::dcast(a_t[, .N, by = .(parcel_id, crop_name)], parcel_id ~ crop_name, fun.agg = function(x) sum(x), value.var = "N"))

  land_div <- as.numeric(rarefaction(as.matrix(x[, safety := 0][, -c("parcel_id")]), 
                                      method = 'IBR', 
                                      effort = c(1 : max_effort, seq(from = 2*max_effort, to = max_effort^2, by = max_effort)),
                                      extrapolate = TRUE, 
                                      quiet_mode = TRUE))

  spce_div <- rowMeans(sapply(1:4, FUN = smpl_space_y, a = a_t, samp_method = 'IBR', max_effort = max_effort))
  time_div <- rowMeans(sapply(1:4, FUN = smpl_time_y, a = a_t, samp_method = 'IBR', max_effort = max_effort))

pdf(paste0("output/scenarios/",pdf_name[i]))
  plot(land_div[1:20], type = 'l', ylim = c(0, 10))
  points(spce_div[1:8], pch= 19)
  points(spce_div[1:20], type = 'l', lty = 2)
  points(time_div[1:8], pch = 22)
  points(time_div[1:20], type = 'l', lty = 3)
dev.off()

}