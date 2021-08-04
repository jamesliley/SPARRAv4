######################################################
## Code to generate constituent models for super-   ##
##  learner, given a fully constructed data matrix  ##
######################################################
##
## James Liley, 2020
##

## More models can be added. Specifications for a model are:
## Called fit_model.${model}. The SPARRAv4fit function will search the environment for any functions of this form
## Parameters X (training data), seed  (random seed), hyperparameter grid values or preset values, and return_hyp
## If return_hyp=T, should return a list of: 'fit' - a function taking an input X and making a prediction, hyp_options: list of options for hyperparameter values, hyp_performance: performance at each hyperparameter value
## If return_hyp=F, just return 'fit'
## All functions will be called with default values for any hyperparameters. To put in several values, define a general public function fit.${model} and private functions fit_model.${model}. See below
## See functions below for examples


######################################################
## SPARRA version 3                                 ##
######################################################

##' Fits a model in the same way as SPARRAv3. 
##' 
##' 
##' @name fit.v3
##' @param X data matrix, assumed to contain column 'target'
##' @param seed random seed; defaults to 214. Does not do anything in this case as model is deterministic.
##' @param return_hyp set to TRUE to additionally return an object detailing fit at each hyperparameter value. Returns 'null' here.
##' @param ... other parameters. For consistency only, ignored here.
##' @return a list of four objects. 1. pred: a string defining a function taking a data-frame argument X (including target) and a model which predicts outcome - note this is returned as a CHARACTER STRING, NOT A CLOSURE. 2. model: a model which can be saved; makes fit work sometimes. 3. hyp_options: list containing options for hyperparameters. 4. hyp_performance: performance at each of these hyperparameter options.
fit.v3=function(X,seed=214,return_hyp=F,...) {
  
  set.seed(seed)	
  
  func="function(X,model) X$v3score"

  if (return_hyp) 
    return(list(pred=func,model=NULL,hyp_options=NULL,hyp_performance=NULL))
  else return(eval(parse(text=func)))
}

fit_model.v3=function(...) fit.v3(...)






######################################################
## Tree models: random forest and XGboost           ##
######################################################

##' Fit random forest to data using Ranger package
##' 
##' @name fit.rf.ranger
##' @param X data matrix, assumed to contain column 'target'
##' @param num.trees number of trees for random forest
##' @param seed random seed; defaults to 214
##' @param return_hyp set to TRUE to additionally return an object detailing fit at each hyperparameter value. Returns 'null' here.
##' @param ... further parameters passed to function ranger
##' @return a list of four objects. 1. pred: a string defining a function taking a data-frame argument X (including target) and a model which predicts outcome - note this is returned as a CHARACTER STRING, NOT A CLOSURE. 2. model: a model which can be saved; makes fit work sometimes. 3. hyp_options: list containing options for hyperparameters. 4. hyp_performance: performance at each of these hyperparameter options.
fit.rf.ranger=function(X,num.trees=1000,seed=214,return_hyp=F,...) {

set.seed(seed)	

# Deal with target; needs to be quite specific
X$target=as.logical(X$target)  	
  
model <- ranger(target ~ .,
		X %>% mutate_at(vars(-target), replace_na, replace=-999) %>% 
      mutate_at(vars(target),as.numeric),
		num.trees = num.trees,
		probability = TRUE,
    max.depth=50,
	  ...)

func="function(X,model=NULL) predict(model,
	X %>% mutate_at(vars(-target), 
		replace_na, replace=-999)  %>%
    mutate_at(vars(target),as.numeric))$predictions[,2]"

if (return_hyp) 
  return(list(pred=func,model=model,hyp_options=NULL,hyp_performance=NULL))
else return(eval(parse(text=func)))
}


# The ranger model tends to take an impractically long time to fit.
# fit_model.rf.ranger=function(...) fit.rf.ranger(num.trees=15,...) 



##' Fit random forest to data using h2o package
##' 
##' @name fit.rf.h2o
##' @param X data matrix, assumed to contain column 'target'
##' @param num.trees number of trees for random forest
##' @param max.depth maximum depth of trees
##' @param seed random seed; defaults to 214
##' @param init set to TRUE to run h2o.removeAll() first
##' @param return_hyp set to TRUE to additionally return an object detailing fit at each hyperparameter value. Returns 'null' here.
##' @param ... further parameters passed to function h2o.randomForest
##' @return a list of four objects. 1. pred: a string defining a function taking a data-frame argument X (including target) and a model which predicts outcome - note this is returned as a CHARACTER STRING, NOT A CLOSURE. 2. model: a model which can be saved; makes fit work sometimes. 3. hyp_options: list containing options for hyperparameters. 4. hyp_performance: performance at each of these hyperparameter options.
fit.rf.h2o=function(X,num.trees=1000,seed=214,return_hyp=FALSE,init=FALSE,max.depth=20,...) {

if (init) h2o.removeAll()

tr <- as.h2o(X)

N=dim(X)[1]
fit_time= 1e6 #round(N*log(N)/6000)*((max.depth/20))

model <- h2o.randomForest(training_frame = tr,
  		y = "target",
			ntrees = num.trees,
      max_depth=max.depth,
			balance_classes = FALSE,
			max_runtime_secs = fit_time,
			seed = seed,
      ...)

func="function(X,model=NULL) as.data.frame(h2o.predict(model, as.h2o(X)))[,3]"

if (return_hyp) 
  return(list(pred=func,model=model,hyp_options=NULL,hyp_performance=NULL))
else return(eval(parse(text=func)))
}

## Internal function for using the above in super-learner
fit_model.rf1.h2o=function(...) fit.rf.h2o(max.depth=20,num.trees=500,...)
fit_model.rf2.h2o=function(...) fit.rf.h2o(max.depth=40,num.trees=50,...)


##' Fit boosted tree models to data using xgboost
##' 
##' @name fit.xgb
##' @param X data matrix, assumed to contain column 'target'
##' @param train_proportion proportion of samples to use for training set; defaults to 0.8
##' @param seed random seed; defaults to 215
##' @param max_depth maximum tree depth
##' @param eta scale for learning rate: scale contribution of new trees
##' @param gamma minimum loss reduction to add a leaf node
##' @param return_hyp set to TRUE to additionally return an object detailing fit at each hyperparameter value. Returns 'null' here.
##' @param ... further parameters passed to function xgb.train
##' @return a list of four objects. 1. pred: a string defining a function taking a data-frame argument X (including target) and a model which predicts outcome - note this is returned as a CHARACTER STRING, NOT A CLOSURE. 2. model: a model which can be saved; makes fit work sometimes. 3. hyp_options: list containing options for hyperparameters. 4. hyp_performance: performance at each of these hyperparameter options.
fit.xgb=function(X,train_proportion=0.8,seed=215,max_depth=4,eta=0.3,gamma=0,objective="binary:logistic",return_hyp=FALSE,...) {

set.seed(seed)	

# Training and hyperparameter (ntree) submatrices, converted to xgb form
ntot=dim(X)[1]; nfit=round((1-train_proportion)*ntot);
ind_fit=order(runif(ntot))[1:nfit]; ind_hyp=setdiff(1:ntot,ind_fit)

ix=as.numeric(X %>% pull(target))
offset=min(ix,na.rm=T)

xgfit <- xgb.DMatrix(data = as.matrix(X[ind_fit,] %>% select(-target)), label = ix[ind_fit]-offset)
xghyp <- xgb.DMatrix(data = as.matrix(X[ind_hyp,] %>% select(-target)), label = ix[ind_hyp]-offset)

model <- xgb.train(params =
	list(objective = objective, 
			 max_depth = max_depth,
		   gamma=gamma,
		   eta=eta),
			 data = xgfit,
			 nrounds = 10000,
			 watchlist = list(train = xgfit, test = xghyp),
			 early_stopping_rounds = 25,
			 nthread = 20,
	     ...)


func="function(X,model=NULL) predict(model,
	xgb.DMatrix(data = as.matrix(X %>% select(-target))))"

if (return_hyp) 
  return(list(pred=func,model=model,hyp_options=NULL,hyp_performance=NULL))
else return(eval(parse(text=func)))
}

## Internal function for using the above in super-learner
fit_model.xgb1.xgb=function(...) fit.xgb(...,max_depth=4,eta=0.3,gamma=0)
fit_model.xgb2.xgb=function(...) fit.xgb(...,max_depth=8,eta=0.075,gamma=5)
fit_model.xgb3.xgb=function(...) fit.xgb(...,max_depth=3,eta=0.3,gamma=0)



	

##' Fit feed-forward neural network to data using h2o package, fitting hyperparameters if specified
##' If any of the parameters n_nodes, n_layer, l1_pen, l2_pen has length>1, then a hyperparameter fitting procedure is performed.
##' 
##' @name fit.ann.h2o
##' @param X data matrix, assumed to contain column 'target'
##' @param n_nodes number of nodes in each layer
##' @param n_layer number of layers
##' @param l1_pen L1-penalty to edge inclusion (add l1_pen*sum(|edge weight|) to objective)
##' @param l2_pen L1-penalty to edge inclusion (add l2_pen*sum(|edge weight|^2) to objective)
##' @param seed random seed; defaults to 214
##' @param init set to TRUE to run h2o.removeAll() first
##' @param return_hyp set to TRUE to additionally return an object detailing fit at each hyperparameter value
##' @param hyp_proportion proportion of training set used in hyperparameter optimisation
##' @param ... further parameters passed to function h2o.deeplearning
##' @return a list of four objects. 1. pred: a string defining a function taking a data-frame argument X (including target) and a model which predicts outcome - note this is returned as a CHARACTER STRING, NOT A CLOSURE. 2. model: a model which can be saved; makes fit work sometimes. 3. hyp_options: list containing options for hyperparameters. 4. hyp_performance: performance at each of these hyperparameter options.
fit.ann.h2o=function(X,
n_nodes=c(128,256), n_layer=1:2, 
input_dropout_ratio=0.2,
seed=214,init=FALSE,return_hyp=FALSE,hyp_proportion=0.8,...) {
	
set.seed(seed)
	
if (init) h2o.removeAll()
	
if (max(c(length(n_nodes),length(n_layer),length(input_dropout_ratio)))>1) hyper=T else hyper=F

	
if (hyper) { 	## Optimise hyperparameters
  	## List of hyperparameter options  
  	hyp_options=list() # list of options
  	for (i in 1:length(n_nodes)) 
  	  for (j in 1:length(n_layer)) 
  	    for (k in 1:length(input_dropout_ratio)) {
  		hyp_options[[1+length(hyp_options)]]=
  		  list(hidden=rep(n_nodes[i],n_layer[j]), input_dropout_ratio=input_dropout_ratio[k])
	  }
	
		# Training and hyperparameter (ntree) submatrices, converted to xgb form
		ntot=dim(X)[1]; nfit=round((1-hyp_proportion)*ntot);
		ind_fit=order(runif(ntot))[1:nfit]; ind_hyp=setdiff(1:ntot,ind_fit)
		
		# Parameter and hyparparameter fitting matrices
		fit <- as.h2o(X[ind_fit,] %>% mutate_at(vars(-target), replace_na, replace=-99))
		hyp <- as.h2o(X[ind_hyp,] %>% mutate_at(vars(-target), replace_na, replace=-99))  
		
		# Decide hyperparameters
		hyp_performance=lapply(hyp_options,hyp_evaluate,
				method=h2o.deeplearning,
				train=fit,test=hyp,metric="AUC",
				y = "target",
				activation = "RectifierWithDropout",
				balance_classes =FALSE,
				epochs = 5000,
				classification_stop = -1,
				train_samples_per_iteration = -1,
				max_runtime_secs = 8*3600,...)  # 3*3600
    
		print("Fitted deep learning hyperparameters")

		hyp_top=hyp_options[[which.max(unlist(hyp_performance))]]
    hyp_options_out=list()
    for (ix in 1:length(hyp_options))
      hyp_options_out[[ix]]=list(n_nodes=hyp_options[[ix]]$hidden[1],
        n_layer=length(hyp_options[[ix]]$hidden),
        input_dropout_ratio=hyp_options[[ix]]$input_dropout_ratio)
    h2o.rm(fit); h2o.rm(hyp); h2o.removeAll(); rm(list=c("fit","hyp"))
} else {
  hyp_options=NULL
  hyp_options_out=NULL
  hyp_performance=NULL
  hyp_top=list(hidden=rep(n_nodes,n_layer), l1=l1_pen, l2=l2_pen)
}

tr <- as.h2o(X %>% mutate_at(vars(-target), replace_na, replace=-99))

print("Fitting main model")


## Fit model for fold given optimal hyperparameters
model <- h2o.deeplearning(training_frame = tr,
	y = "target",
  hidden = hyp_top$hidden,
  activation = "RectifierWithDropout",
	balance_classes = FALSE,
	input_dropout_ratio = hyp_top$input_dropout_ratio,
	epochs = 5000,
	classification_stop = -1,
	train_samples_per_iteration = -1,
	max_runtime_secs = 8*3600,...)


func="function(X,model=NULL) as.data.frame(h2o.predict(model, as.h2o(X)))[,3]"

if (return_hyp) 
	return(list(pred=func,model=model,hyp_options=hyp_options_out,hyp_performance=hyp_performance))
else return(eval(parse(text=func)))
}

## Internal function for use in super_leaner
fit_model.ann.h2o=function(...) fit.ann.h2o(...)



####################################################################
## Support vector machine                                         ##
####################################################################

##' Fit SVM to data using h2o package (radial kernel), fitting hyperparameters if specified
##' If any of the parameters hyper_param,gamma has length>1, then a hyperparameter fitting procedure is performed.
##' 
##' @name fit.svm.h2o
##' @param X data matrix, assumed to contain column 'target'
##' @param hyper_param penalty parameter (leniency for objective function)
##' @param gamma kernel coefficient
##' @param seed random seed; defaults to 214
##' @param init set to TRUE to run h2o.removeAll() first
##' @param return_hyp set to TRUE to additionally return an object detailing fit at each hyperparameter value
##' @param hyp_proportion proportion of training set used in hyperparameter optimisation
##' @param ... further parameters passed to function h2o.psvm
##' @return a list of four objects. 1. pred: a string defining a function taking a data-frame argument X (including target) and a model which predicts outcome - note this is returned as a CHARACTER STRING, NOT A CLOSURE. 2. model: a model which can be saved; makes fit work sometimes. 3. hyp_options: list containing options for hyperparameters. 4. hyp_performance: performance at each of these hyperparameter options.
fit.svm.h2o=function(X,
	hyper_param=c(0.1,0.5,0.8,1,1.2,2,10), gamma=c(0.1,0.5,0.8,1,1.2,2,10),
	seed=214,init=FALSE,return_hyp=FALSE,hyp_proportion=0.8,...) {
	
	set.seed(seed)
	
	if (init) h2o.removeAll()
	
	if (max(c(length(hyper_param),length(gamma)))>1) hyper=T else hyper=F
	
	
	if (hyper) { 	## Optimise hyperparameters
		## List of hyperparameter options  
		hyp_options=list() # list of options
		for (i in 1:length(hyper_param)) for (j in 1:length(gamma)) {
			hyp_options[[1+length(hyp_options)]]=list(hyper_param=hyper_param[i],gamma=gamma[j])
		}
		
		# Training and hyperparameter (ntree) submatrices, converted to xgb form
		ntot=dim(X)[1]; nfit=round((1-hyp_proportion)*ntot);
		ind_fit=order(runif(ntot))[1:nfit]; ind_hyp=setdiff(1:ntot,ind_fit)
		
		# Parameter and hyparparameter fitting matrices
		fit <- as.h2o(X[ind_fit,] %>% mutate_at(vars(-target), replace_na, replace=-99))
		hyp <- as.h2o(X[ind_hyp,] %>% mutate_at(vars(-target), replace_na, replace=-99))  
		
		# Decide hyperparameters.
		hyp_performance=lapply(hyp_options,hyp_evaluate,
			method=h2o.psvm,
			train=fit,test=hyp,metric="AUC",
			y = "target",
			max_iterations=200,...)
		
		hyp_top=hyp_options[[which.max(unlist(hyp_performance))]]
	} else {
	  hyp_options=NULL
	  hyp_performance=NULL
	  hyp_top=list(hyper_param=hyper_param,gamma=gamma)
	}
	
	tr <- as.h2o(X %>% mutate_at(vars(-target), replace_na, replace=-99))
	
	## Fit model for fold given optimal hyperparameters
	model <- h2o.psvm(training_frame = fit,
		hyper_param=hyp_top$hyper_param,
		gamma=hyp_top$gamma,
		y = "target",
		max_iterations = 200,...)
	
	pred="function(X,model=NULL) as.data.frame(h2o.predict(model, as.h2o(X)))[,3]"
	
	if (return_hyp) 
		return(list(pred=func,model=model,hyp_options=hyp_options,hyp_performance=hyp_performance))
	else return(eval(parse(text=func)))
}

## Internal function for use in super-learner
## fit_model.svm.h2o=function(...) fit.svm.h2o(...) ## not currently working



	
####################################################################
## Naive Bayes                                                    ##
####################################################################

##' Fit naive Bayes model to data using h2o package, fitting hyperparameters if specified
##' If parameter laplace has length>1, then a hyperparameter fitting procedure is performed.
##' 
##' @name fit.nb.h2o
##' @param X data matrix, assumed to contain column 'target'
##' @param laplace Laplace hyperparameter values to consider
##' @param seed random seed; defaults to 214
##' @param init set to TRUE to run h2o.removeAll() first
##' @param return_hyp set to TRUE to additionally return an object detailing fit at each hyperparameter value
##' @param hyp_proportion proportion of training set used in hyperparameter optimisation
##' @param ... further parameters passed to function h2o.nb
##' @return a list of four objects. 1. pred: a string defining a function taking a data-frame argument X (including target) and a model which predicts outcome - note this is returned as a CHARACTER STRING, NOT A CLOSURE. 2. model: a model which can be saved; makes fit work sometimes. 3. hyp_options: list containing options for hyperparameters. 4. hyp_performance: performance at each of these hyperparameter options.
fit.nb.h2o=function(X,
	laplace=0:4,
	seed=214,init=FALSE,return_hyp=FALSE,hyp_proportion=0.8,...) {
	
	set.seed(seed)
	
	if (init) h2o.removeAll()
	
	if (max(c(length(laplace)))>1) hyper=T else hyper=F
	
	
	if (hyper) { 	## Optimise hyperparameters
		## List of hyperparameter options  
		hyp_options=list() # list of options
		for (i in 1:length(laplace)) {
			hyp_options[[1+length(hyp_options)]]=list(laplace=laplace[i])
		}
		
		# Training and hyperparameter (ntree) submatrices, converted to xgb form
		ntot=dim(X)[1]; nfit=round((1-hyp_proportion)*ntot);
		ind_fit=order(runif(ntot))[1:nfit]; ind_hyp=setdiff(1:ntot,ind_fit)
		
		# Parameter and hyparparameter fitting matrices
		fit <- as.h2o(X[ind_fit,] %>% mutate_at(vars(-target), replace_na, replace=-99))
		hyp <- as.h2o(X[ind_hyp,] %>% mutate_at(vars(-target), replace_na, replace=-99))  
		
		# Decide hyperparameters.
		hyp_performance=lapply(hyp_options,hyp_evaluate,
			method=h2o.naiveBayes,
			train=fit,test=hyp,metric="AUC",
			y = "target",
			max_runtime_secs=10*60,...)

		hyp_top=hyp_options[[which.max(unlist(hyp_performance))]]
	} else {
	  hyp_options=NULL
	  hyp_performance=NULL
	  hyp_top=list(laplace=laplace)
	}
	
	
	tr <- as.h2o(X %>% mutate_at(vars(-target), replace_na, replace=-99))
	
	## Fit model for fold given optimal hyperparameters
	model <- h2o.naiveBayes(training_frame = tr,
		y = "target",
		laplace=hyp_top$laplace,
		max_runtime_secs = 60,...)

	func="function(X,model=NULL) as.data.frame(h2o.predict(model, as.h2o(X)))[,3]"
	
	if (return_hyp) 
		return(list(pred=func,model=model,hyp_options=hyp_options,hyp_performance=hyp_performance))
	else return(eval(parse(text=func)))
}

## Internal function for super learner

fit_model.nb.h2o=function(...) fit.nb.h2o(...)


####################################################################
## Generalised linear model, data not partitioned first           ##
####################################################################

##' Fit generalised linear model to data using h2o package, fitting hyperparameters if specified
##' If either of lambda, alpha has length>1, then a hyperparameter fitting procedure is performed.
##' 
##' @name fit.glm.h2o
##' @param X data matrix, assumed to contain column 'target'
##' @param lambda Ratio of L1/L2 parameters (elastic net)
##' @param alpha Sum of L1/L2 parameters (elastic net)
##' @param seed random seed; defaults to 214
##' @param init set to TRUE to run h2o.removeAll() first
##' @param return_hyp set to TRUE to additionally return an object detailing fit at each hyperparameter value
##' @param hyp_proportion proportion of training set used in hyperparameter optimisation
##' @param ... further parameters passed to function h2o.nb
##' @return a list of four objects. 1. pred: a string defining a function taking a data-frame argument X (including target) and a model which predicts outcome - note this is returned as a CHARACTER STRING, NOT A CLOSURE. 2. model: a model which can be saved; makes fit work sometimes. 3. hyp_options: list containing options for hyperparameters. 4. hyp_performance: performance at each of these hyperparameter options.
fit.glm.h2o=function(X,
	lambda=c(0,10^-seq(5,1,length.out=5)), alpha=c(0,0.5,1),
	seed=214,init=FALSE,return_hyp=FALSE,hyp_proportion=0.8,...) {
	
	set.seed(seed)
	
	if (init) h2o.removeAll()
	
	if (max(c(length(lambda),length(alpha)))>1) hyper=T else hyper=F
	
	
	if (hyper) { 	## Optimise hyperparameters
		## List of hyperparameter options  
		hyp_options=list() # list of options
		for (i in 1:length(lambda)) for (j in 1:length(alpha)) {
			hyp_options[[1+length(hyp_options)]]=list(lambda=lambda[i],alpha=alpha[j])
		}
		

		# Training and hyperparameter (ntree) submatrices, converted to xgb form
		ntot=dim(X)[1]; nfit=round((1-hyp_proportion)*ntot);
		ind_fit=order(runif(ntot))[1:nfit]; ind_hyp=setdiff(1:ntot,ind_fit)
		
		# Parameter and hyparparameter fitting matrices
		fit <- as.h2o(X[ind_fit,] %>% mutate_at(vars(-target), replace_na, replace=-99))
		hyp <- as.h2o(X[ind_hyp,] %>% mutate_at(vars(-target), replace_na, replace=-99))  
		
		# Decide hyperparameters.
		hyp_performance=lapply(hyp_options,hyp_evaluate,
			method=h2o.glm,
			train=fit,test=hyp,metric="AUC",
			y = "target",
			remove_collinear_columns=TRUE,
			intercept=TRUE,
			family="binomial",...)
		
		hyp_top=hyp_options[[which.max(unlist(hyp_performance))]]
	} else {
	  hyp_options=NULL
	  hyp_performance=NULL
	  hyp_top=list(lambda=lambda,alpha=alpha)
	}
	
	tr <- as.h2o(X %>% mutate_at(vars(-target), replace_na, replace=-99))
	
	## Fit model for fold given optimal hyperparameters
	model <- h2o.glm(training_frame=tr,
		y = "target",
		remove_collinear_columns=TRUE,
		intercept=TRUE,
		lambda=hyp_top$lambda,
		alpha=hyp_top$alpha,
		family="binomial",...)

	func="function(X,model=NULL) as.data.frame(h2o.predict(model, as.h2o(X)))[,3]"
	
	if (return_hyp) 
		return(list(pred=func,model=model,hyp_options=hyp_options,hyp_performance=hyp_performance))
	else return(eval(parse(text=func)))
}

## Internal function for super-learner
fit_model.glm.h2o=function(...) fit.glm.h2o(...)