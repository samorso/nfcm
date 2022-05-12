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
k <- 7

starting_value = matrix(nrow = MC, ncol = k * k)
spline_coef = matrix(nrow = MC, ncol = k * k)
oracle = matrix(nrow = MC, ncol = k * k)
convergence = rep(NA_integer_, MC)
iteration = rep(NA_integer_, MC)
time = matrix(nrow = MC, ncol = 3)
nll = matrix(nrow = MC, ncol = 3)

##------------------ Slurm specs --------------
n_array <- 1000
ind <- matrix(seq_len(MC),nr=n_array,byr=T)

for(i in seq_len(n_array)){
  if(!file.exists(file=paste0("tmp/",assumed_model,"_id_",i,".rds"))) next
  load(file=paste0("tmp/",assumed_model,"_id_",i,".rds"))
  starting_value[ind[i,],] <- res$starting_value[ind[i,],]
  spline_coef[ind[i,],] <- res$spline_coef[ind[i,],]
  oracle[ind[i,],] <- res$oracle[ind[i,],]
  convergence[ind[i,]] <- res$convergence[ind[i,]]
  iteration[ind[i,]] <- res$iteration[ind[i,]]
  time[ind[i,],] <- res$time[ind[i,],]
  nll[ind[i,],] <- res$nll[ind[i,],]
}

res <- list(
  starting_value = starting_value,
  spline_coef = spline_coef,
  oracle = oracle,
  convergence = convergence,
  iteration = iteration,
  time = time,
  nll = nll
)

save(res, file=paste0("data/",assumed_model,"_",spline_type,"spline","_n_",n,"_d_",d,"_param_",param[[1]],".rds"))
