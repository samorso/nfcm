# -----------
# Simulations
# -----------
# general setting

# simulation specifics
MC <- 1000 # number of simulations

spline_type <- as.character(Sys.getenv("SPLINE"))
assumed_model <- as.character(Sys.getenv("MODEL"))
n <- as.integer(Sys.getenv("N")) # sample size
param <- as.list(as.double(Sys.getenv(paste0("PARAM",1:5))))
names(param) <- as.character(Sys.getenv(paste0("NPARAM",1:5)))
d <- as.integer(Sys.getenv("D")) # number of dimension (min 2)
k <- 6

starting_value = matrix(nrow = MC, ncol = k * k)
spline_coef_G = matrix(nrow = MC, ncol = k * k)
spline_coef_H = matrix(nrow = MC, ncol = k * k)
oracle_G = matrix(nrow = MC, ncol = k * k)
oracle_H = matrix(nrow = MC, ncol = k * k)
convergence = rep(NA_integer_, MC)
iteration = rep(NA_integer_, MC)
time = matrix(nrow = MC, ncol = 5)
nll = matrix(nrow = MC, ncol = 3)
mise = matrix(nrow = MC, ncol = 5)
rmse = matrix(nrow = MC, ncol = 2)

##------------------ Slurm specs --------------
n_array <- 1000
ind <- matrix(seq_len(MC),nr=n_array,byr=T)

for(i in seq_len(n_array)){
  if(!file.exists(file=paste0("tmp/",assumed_model,"_id_",i,".rds"))) next
  load(file=paste0("tmp/",assumed_model,"_id_",i,".rds"))
  starting_value[ind[i,],] <- res$starting_value[ind[i,],]
  spline_coef_G[ind[i,],] <- res$spline_coef_G[ind[i,],]
  spline_coef_H[ind[i,],] <- res$spline_coef_H[ind[i,],]
  oracle_G[ind[i,],] <- res$oracle_G[ind[i,],]
  oracle_H[ind[i,],] <- res$oracle_H[ind[i,],]
  convergence[ind[i,]] <- res$convergence[ind[i,]]
  iteration[ind[i,]] <- res$iteration[ind[i,]]
  time[ind[i,],] <- res$time[ind[i,],]
  nll[ind[i,],] <- res$nll[ind[i,],]
  mise[ind[i,],] <- res$mise[ind[i,],]
  rmse[ind[i,],] <- res$rmse[ind[i,],]
}

res <- list(
  starting_value = starting_value,
  spline_coef_G = spline_coef_G,
  spline_coef_H = spline_coef_H,
  oracle_G = oracle_G,
  oracle_H = oracle_H,
  convergence = convergence,
  iteration = iteration,
  time = time,
  nll = nll,
  mise = mise,
  rmse = rmse
)

save(res, file=paste0("data/",assumed_model,"_",spline_type,"spline","_n_",n,"_d_",d,"_param_",param[[1]],".rds"))
