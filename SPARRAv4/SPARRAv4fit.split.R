## SPARRA v4 R package

# Preliminaries to calling SPARRAv4.fit
# Clean and load raw data
# Fit topic model in fold-dependent way
# This is a version of the package which can be run across different machines with different packages


##' SPARRAv4.fit
##' Function to fit a SPARRAv4 - type model to a data matrix
##' 
##' Will save the following files. Determines a name based on input X if model_id is not set.
##' ${model_id}.RDS: a list containing the following
##'  1. train,test:         Indices of raining and test sets
##'  2. mod.train.{model}:  All models trained to training set, including any hyperparameter values and performance at those values
##'  3. pred.train.{model}: Predictions on test set from all models in 2.
##'  4. mod.full{model}:    All models trained to full dataset, including any hyperparameter values and performance at those values
##'  5. slfit:              Fitted super-learner model: linear combiantion of predictions in (3) which best fit target
##'  6. ctransform:         Fitted model with calibrating function: monotonic function applied to output of (5) to optimise calibration
##' ${model_id}.pred.rds
##'  1. function.{model}:   Predictor derived from model fitted to full dataset, saved as STRING describing function not function (closure) itself. Note this is only part of mod.full.{model} above. In general, these will need the appropriate model loaded in order to work.
##'  2. slfit:              Fitted super-learner model: as above
##'  3. ctransform:         Fitted model with calibrating function: as above
##' ${model_id}.{model}
##'  For each constituent model of the super-learner, this is a saved version of it, which must be loaded in order to use the predictor function.
##' 
##' SPARRAv4.fit returns a function which can be applied to a new dataset and will generate the predictors. 
##' We recommend using the separate function predict.sparra() which uses the saved models, and will work even if R is restarted between running this function and making predictions.
##' 
##' @name SPARRAv4.fit
##' @param X data matrix including topic probabilities, assumed to contain a column 'target'
##' @param options list of models and parameters for models. Each list item is of the form model=list(option1=value1,option2=value2...). The value model should correspond to a function called fit_model.[model]. If options=NULL, the workspace is searched for functions of this type, and each is used with default options. See constituent_models.R for a specification of the function definition.
##' @param train list of indices of training subset of X
##' @param seed random seed
##' @param super_learner_loss passed to cv.glmnet; loss metric for super-learner. One of 'deviance', 'class', 'auc', 'mse'. Defaults to AUC.
##' @param model_id save model as ${model_id}.RData; can be fit on different machines. Set as NULL to not save
##' @param verbose print progress
##' @param linux TEMPORARY set to TRUE if on linux machine. 
##' @return object of SPARRAv4 class
SPARRAv4.fit.split = function(X,options=NULL,
	train=NULL,seed=220,super_learner_loss="auc",model_id=NULL,verbose=TRUE, linux=FALSE) {
	
	### Preliminaries
	if (verbose) print("Starting")

  ### In order to work, this function needs to save various elements. This path indicates where to do it.
  ### If this is not set, then chooses a file named deterministically according to X. 
  # if (is.null(model_id)) model_id=paste0("mod",as.character(as.numeric(Sys.time())*10000))
  if (is.null(model_id)) {
    set.seed(X)
    model_id=paste0("SPARRAv4mod",sample(1e8,1))
  }
  
  
  ## Short auxiliary functions
  # Saves a model
  savemodel=function(model,model_type,file) {
    if (grepl("h2o",model_type)) {
      fx=h2o.saveModel(model,path=dirname(file),force=TRUE)
      file.rename(fx,file)
    }
    if (!grepl("h2o",model_type)) saveRDS(model,file=file)
  }
  # Predicts on a large matrix
  longpred=function(model,Xpred) {
    npred=dim(Xpred)[1]     # Following block of code to prevent memory problems
    if (npred> 2e6) {
      fullpred=c()
      nblock=floor(dim(Xpred)[1]/(2e6))
      for (ii in 1:(1+nblock)) {
        print(paste0("Computing predictions for block ",ii," of ",1+nblock))
        xsubpred=Xpred[(1+(ii-1)*(2e6)):min(npred,ii*(2e6)),]
        modelfit=eval(parse(text=model$pred)) # The fitting function is stored as a string, to avoid saving a huge closure
        subpred=modelfit(xsubpred,model=model$model)
        fullpred=c(fullpred,subpred)
        rm(list=c("xsubpred","subpred"))
      }
    } else {
    	modelfit=eval(parse(text=model$pred))
    	fullpred=modelfit(Xpred,model=model$model)
    }
    return(fullpred)
  }	
  # Given a list of names of variables, turn them into a list and save as RDS
  gsave=function(xnames,file) {
    if (length(xnames)>0) {
    xlist=list()
    for (i in 1:length(xnames)) xlist[[i]]=get(xnames[i])
    names(xlist)=xnames
    saveRDS(xlist,file)
    }
  }
  # Save two important lists
  xsave=function() {
    gsave(save_list,id_file)
	  gsave(short_save_list,short_file)
  }
  
  
  # Random seed. Note this has to be after the code above determining model_id.
	set.seed(seed)	


	# Search workspace for models
	if (is.null(options)) {
	  models=gsub("fit_model.","",grep("fit_model.*",lsf.str(pos=1),val=T))
	  options=list(); for (i in 1:length(models)) options[[i]]=list()
	  names(options)=models
	} else {
	  models=names(options)
	}	

	### Names of files to save. id_file is the location where everything relevant to this function. short_file is a list of only the objects needed to use the derived predictor.
	id_file=paste0(model_id,".RDS")
	short_file=paste0(model_id,".pred.RDS")

	### Load if already exists
	if (file.exists(id_file)) {
	  #  full save list
	  saved=readRDS(id_file)
	  save_list=names(saved)
	  for (i in 1:length(save_list)) assign(save_list[i],saved[[i]])
	  rm(list="saved")
	  #  short save list
	  short_saved=readRDS(short_file)
	  short_save_list=names(short_saved)
	  for (i in 1:length(short_save_list)) assign(short_save_list[i],short_saved[[i]])
	  rm(list="short_saved")
	} else {
	  save_list=c()
	  short_save_list=c() 
	}
	# Throughout:
	#  save_list is a comprehensive list of all variables used; 
	#  short_save_list is a minimal list for using the derived predictor function.
	if (verbose) print("Loaded existing models if already saved")
	
	
	
	### Training and testing sets.  Generate if not already generated.
	if (!all(c("train","test") %in% save_list)) { # ie, if already generated
	if (is.null(train))	train=which(runif(dim(X)[1])>0.5)
	test=setdiff(1:dim(X)[1],train)
	save_list=unique(c(save_list,"train","test"))
	xsave()
	}
	
	
	# Does NOT save Xtrain or Xtest as they are usually prohibitively large. Sets them here.
	Xtrain=X[train,]
	Xtest=X[test,]
	if (verbose) print("Determined test and training sets")
	
	
	### Initialise h2o
	if (linux) {   ############### TEMPORARY WHILE SPLITTING MACHINES ########################
	h2o.init(bind_to_localhost = FALSE, max_mem_size = "125G")
	}   ############### TEMPORARY WHILE SPLITTING MACHINES ########################
  if (verbose) print("Initialised h2o")
	
	
	### Fit constituent models to Xtrain, fitting hyperparameters where appropriate (1)
  for (i in 1:length(models)) {
  	############### TEMPORARY WHILE SPLITTING MACHINES ########################
  	if ((linux & grepl("h2o",models[i])) | (!linux & !grepl("h2o",models[i]))) {
  	############### TEMPORARY WHILE SPLITTING MACHINES ########################
  	
  	  
  	if (!(paste0("mod.full.",models[i]) %in% save_list)) { # Fit model if it has not already been saved

  	# Get relevant fitting function
   	modf=get(paste0("fit_model.",models[i]))
  	
  	# Model trained only to training set, with predictions on test set
    mod.train=do.call(modf,c(list(X = Xtrain,return_hyp = TRUE,seed=seed+i),options[[i]]))
    assign(paste0("pred.train.",models[i]),longpred(mod.train,Xtest)) # predictions on test set. Note that 'pred.train' refers to the fact that the model was fitted to the training data, although predictions are made on the test data.
    # savemodel(mod.train$model,models[i],paste0(model_id,".train.",models[i])) # In practice we don't usually need this object.
   	mod.train$model=NULL # We don't need the actual model anymore
    assign(paste0("mod.train.",models[i]),mod.train)
    rm(list=c("mod.train"))
    if (linux) h2o.removeAll()
    gc()
    print("Trained model to training set and predicted on test set")
    
    # Model trained to full dataset, with predictions on prediction set
    mod.full=do.call(modf,c(list(X = X,return_hyp = TRUE,seed=seed+i),options[[i]])) 
    savemodel(mod.full$model,models[i],paste0(model_id,".",models[i]))
   	mod.full$model=NULL # We don't need the actual model anymore
    assign(paste0("mod.full.",models[i]),mod.full)
   	assign(paste0("function.",models[i]),mod.full$pred) # predictor function
   	rm(list=c("mod.full"))
   	if (linux) h2o.removeAll()
   	gc()

   	# Save new variables
   	save_list=unique(c(save_list,
   	  paste0(c("mod.train.","pred.train."),models[i]),
   	  paste0(c("mod.full.","function."),models[i])))
   	short_save_list=unique(c(short_save_list,paste0("function.",models[i])))
   	xsave()
  	}
  	  
  	  
   	if (verbose) print(paste0("Fitted model ",models[i]))

   	}   	############### TEMPORARY WHILE SPLITTING MACHINES ########################
  }

	############### TEMPORARY WHILE SPLITTING MACHINES ########################
	if (!linux) {    
	############### TEMPORARY WHILE SPLITTING MACHINES ########################
		
	### Assemble matrix of predicted outcomes on test set
  pred.train=c()
  for (i in 1:length(models)) pred.train=cbind(pred.train,get(paste0("pred.train.",models[i])))
  colnames(pred.train)=models
  pred.train=as.data.frame(pred.train)

  ### Fit super learner
  if (!("slfit" %in% save_list)) {
	slfit=cv.glmnet(as.matrix(pred.train),
	  as.matrix(Xtest$target),
		family = "binomial",
		type.measure = super_learner_loss,
		weights= rep(1,dim(pred.train)[1])) # pred.train$glm.h2o^2)
	
	save_list=unique(c(save_list,"slfit"))
  short_save_list=unique(c(short_save_list,"slfit"))
  xsave()
  }
  
	if (verbose) print("Fitted super learner")
	
	### Calibration-forcing transform
	pred_fit=predict(slfit,as.matrix(pred.train),type="response",s="lambda.min")[,1]
	
	if (!("ctransform" %in% save_list)) {
	
	ctransform=forcecal(as.logical(Xtest$target),as.vector(pred_fit))
	save_list=unique(c(save_list,"ctransform"))
	short_save_list=unique(c(short_save_list,"ctransform"))
  xsave()
	}

	if (verbose) print("Fitted calibrating transform")

	### Predictor function
	predictor=function(Xp) {
		# Matrix of constituent predictors
		pred.full=c()
		for (i in 1:length(models)) {
			mod.pred=get(paste0("mod.full.",models[i]))
			pred.full=cbind(pred.full,longpred(mod.pred,Xp))
		}
		colnames(pred.full)=models
		
		sl.full=predict(slfit,as.matrix(pred.full),type="response",lambda=lambda.min)[,1]
		final.full=ctransform$transform(sl.full)
		return(list(constituents=pred.full,predicted=final.full))
	}
	
	### Return final object
  return(predictor)

	############### TEMPORARY WHILE SPLITTING MACHINES ########################
} else return(NULL) 
  ############### TEMPORARY WHILE SPLITTING MACHINES ########################
} 




##' Stable predictor for SPARRAv4.fit
##' 
##' @name predict.sparra
##' @param X data matrix including topic probabilities, to be predicted on. NOTE - to be agnostic to training data, the topic probabilities must be determined using a topic  model fit ONLY to the training data.
##' @param model_id file path corresponding to a previously fitted SPARRAv4 model
##' @param save_id save results in {save_id}.RDS. If NULL, determined deterministically according to X
##' @param verbose print progress
##' @param linux TEMPORARY set to TRUE if on linux machine. 
##' @param use_loaded this function requires loading several objects from file corresponding to the model. To speed things up when the same model is used to make several predictions (eg for Shapley value estimation) this parameter can carry a named list of values which avoid reloading all variables. It should contain everything in {model_id}.pred.RDS, and every model (machine-relevant) in {model_id}.{model}
##' @return object of SPARRAv4 class
predict.sparra = function(X,model_id,save_id=NULL,verbose=TRUE,linux=FALSE,use_loaded=NULL) {
  
  if (verbose) print("Starting")
  
  # Determine file to save to  
  if (is.null(save_id)) {
    set.seed(X)
    save_id=paste0(model_id,".xpred",sample(1e8,1))
  }
  sfile=paste0(save_id,".RDS")
  
  
  # Must run on Linux machine first, then on Windows.
  if (!linux & !file.exists(sfile)) {
    stop("Run on linux machine first")
  }
  if (!linux & file.exists(sfile)) {
    xlist=readRDS(sfile)
    for (i in 1:length(xlist)) assign(names(xlist)[i],xlist[[i]])
    save_list=names(xlist)
    rm(xlist)
    if (verbose) print("Loaded predictions from h2o models")
  } else save_list=c()

  
  ## Initialise h2o
	if (linux & is.null(use_loaded)) {
	h2o.init(bind_to_localhost = FALSE, max_mem_size = "125G")
	} 
  if (verbose) print("Initialised h2o")

  
  # Load model file
  if (is.null(use_loaded)) {
     xlist=readRDS(paste0(model_id,".pred.RDS"))
     nlist=names(xlist)
     for (i in 1:length(xlist)) assign(nlist[i],xlist[[i]])
     rm(list="xlist")
  } else {
    nlist=names(use_loaded)
    for (i in 1:length(nlist)) assign(nlist[i],use_loaded[[i]])
  }
  if (verbose) print("Loaded predictor functions")  

  # Auxiliary function
  loadmodel=function(model_type,file) {
    if (grepl("h2o",model_type)) {
      fit=h2o.loadModel(file)
      return(fit)
    }
    if (!grepl("h2o",model_type)) {
      fit=readRDS(file)
      return(fit)
    }
  }
  # Given a list of names of variables, turn them into a list and save as RDS
  gsave=function(xnames,file) {
    if (length(xnames)>0) {
      xlist=list()
      for (i in 1:length(xnames)) xlist[[i]]=get(xnames[i])
      names(xlist)=xnames
      saveRDS(xlist,file)
    }
  }

  if (linux) X=as.h2o(X)

  ### Get constituent model 
  models=gsub("function.","",grep("function.",nlist,value=TRUE))
  for (i in 1:length(models)) {
    if ((linux & grepl("h2o",models[i])) | (!linux & !grepl("h2o",models[i]))) {
      if (!exists(paste0("pred.",models[i]))) {
      if (is.null(use_loaded)) fit=loadmodel(models[i],paste0(model_id,".",models[i])) else fit=get(models[i])
      map=eval(parse(text=get(paste0("function.",models[i]))))
      
      assign(paste0("pred.",models[i]),map(X,model=fit))
      
      save_list=unique(c(save_list,paste0("pred.",models[i])))
      }
      if (verbose) print(paste0("Computed predictions for model ",models[i]))
    }   	############### TEMPORARY WHILE SPLITTING MACHINES ########################
  }
  gsave(save_list,sfile)

  if (linux & !is.null(use_loaded)) {
    hx=h2o.ls()
    h2o.rm(setdiff(hx[,1],hx[1:(1+length(grep("h2o",models))),1]))
  }
    
  if (!linux) {
  
  pred.full=c()
  for (i in 1:length(models)) {
    pred.full=cbind(pred.full,get(paste0("pred.",models[i])))
  }
  colnames(pred.full)=models
  
  sl.full=predict(slfit,as.matrix(pred.full),type="response",s="lambda.min")[,1]
  final.full=ctransform$transform(sl.full)
  save_list=unique(c(save_list,"sl.full","final.full"))
  gsave(save_list,sfile)
  
  return(final.full)

  } else return(NULL)
}
  