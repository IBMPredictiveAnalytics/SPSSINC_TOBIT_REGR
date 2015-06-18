#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 1989, 2014
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

helptext="The SPSSINC TOBIT REGR command requires the R Integration Plug-in
and the R AER package.

SPSSINC TOBIT REGR DEPENDENT=dependent variable 
ENTER=independent variables
[DISTRIBUTION={GAUSSIAN** }]
              {EXPONENTIAL}
              {LOGISTIC   }
              {LOGNORMAL  }
              {WEIBULL    }
              {LOGLOGISTIC}

[LOWERBOUND=number]
[UPPERBOUND=number]
[/OPTIONS [MISSING={LISTWISE**}] [EXECUTE={TRUE**}] ]
                   {FAIL      }           {FALSE }
[/SAVE [COEFSDATASET=datasetname] [PROGRAMFILE=filespec] ]

Split files and weight are not honored by this command.

SPSSINC TOBIT REGR /HELP prints this information and does nothing else.

Example:
SPSSINC TOBIT REGR DEPENDENT=mpg ENTER=engine weight
LOWERBOUND=0.

Execute the tobit function from the R AER package.
DEPENDENT and ENTER specify the dependent and independent
variable names.

Categorical independent variables are automatically converted
appropriately to factors.  A constant term is automatically included.

LOWERBOUND and/or UPPERBOUND specify the lower and/or
upper limits for the dependent variable.  If neither is used, there
are no limits applied, and the method amounts to ordinary regression.

By default the error distribution is assumed to be Gaussian, but
any of the other distributions listed can be used instead.

MISSING=LISTWISE causes listwise deletion of missing values.  
FAIL stops the procedure if missing values are encountered.

EXECUTE=FALSE runs the command syntax without running the tobit regression.  
This is mainly useful in combination with SAVE PROGRAMFILE.

/SAVE COEFSDATASET causes a dataset containing the coefficients and the 
error scale parameter to be created.  The dataset name must not already
exist.

PROGRAMFILE causes the R code that implements the tobit regression to be written
to the specified file. Since the tobit function has features not exposed in this 
extension command, the generated program can be a useful starting point for 
additional specifications.
"

tobitreg<-function(dep, enter, distribution="gaussian", missing="listwise", lowerbound=NULL, upperbound=NULL, 
          coefsdataset=NULL){
        
    domain<-"SPSSINC_TOBIT_REGR"
    setuplocalization(domain)
    
    tryCatch(library(AER), error=function(e){
        stop(gettextf("The R %s package is required but could not be loaded.","AER",domain=domain),call.=FALSE)
        }
    )

    if (identical(missing,"listwise")) {missing<-na.exclude} else {missing<-na.fail}

    allvars <- c(dep,enter)
    model <- paste(dep,"~",paste(enter,collapse="+"))

    dta<-spssdata.GetDataFromSPSS(allvars,missingValueToNA=TRUE,factorMode="labels")
        
    lbnd<-lowerbound
    if (is.null(lowerbound)){
        lowerbound<- -Inf
        lbnd<-gettext("None",domain=domain)
    }
    ubnd<-upperbound
    if (is.null(upperbound)){
        upperbound<-Inf
        ubnd<-gettext("None",domain=domain)
    }
    
    bounds = gettextf(" Lower bound: %s, Upper bound: %s",lbnd,ubnd,domain=domain)

    res <- tryCatch(
            summary(restobit <- eval(call("tobit",as.formula(model),left=lowerbound,right=upperbound,dist=distribution,data=dta,na.action=missing)),
            error=function(e) {return(c(gettext("ERROR:",domain=domain),e))})
           )

    if (!is.null(res$message)) {print(res$message)} else {

         wt = waldtest(restobit)
         miss<-ifelse(identical(missing,na.exclude),"na.exclude","na.fail")
         call<-paste("tobit(formula = ",model,", left = ",lowerbound,", right = ",upperbound,", dist = ",dQuote(distribution),
                ", data = dta, na.action = ",miss,")",sep="")
         caption = paste(bounds,"\n",call, 
            sprintf(paste("\n",gettext("Scale: ",domain=domain),"%.4f\n",gettext("Residual d.f.: ",domain=domain),"%s\n",
                gettext("Log likelihood: ",domain=domain),"%.3f",gettext("  D.f.: ",domain=domain),"%.0f\n",
                gettext("Wald statistic: ",domain=domain),"%.3f",gettext("  D.f.: ",domain=domain),"%.0f"),    
                res$scale, restobit$df.residual, res$loglik[[2]], res$df, res$wald, abs(wt$Df[[2]])))
                
        coeff<-res$coefficients[,1:4]
        for (i in 1:length(attributes(coeff)$dimnames[[1]])){
           attributes(coeff)$dimnames[[1]][[i]]=gettext(attributes(coeff)$dimnames[[1]][[i]],domain=domain)
        }

        collabels = c(gettext("Coefficient",domain=domain), gettext("Std. Error",domain=domain), 
                    gettext("z Value",domain=domain), gettext("Sig.",domain=domain))
                    
        StartProcedure(gettext("Tobit Regression",domain=domain),"SPSSINC TOBIT REGR")
        spsspivottable.Display(coeff, 
            title=gettext("Coefficients",domain=domain), templateName="SPSSINCTOBITREGR",
            caption=caption, collabels=collabels,
            isSplit=FALSE)
        spsspkg.EndProcedure()
        
        if (!is.null(coefsdataset)){
            dict<- spssdictionary.CreateSPSSDictionary(c("term", gettext("Variable or Factor Value",domain=domain), 100, "A100", "nominal"),
            c("coefficient", gettext("Estimated Coefficient",domain=domain), 0, "F10.3", "scale"))
            tryCatch({
                spssdictionary.SetDictionaryToSPSS(coefsdataset, dict)
                spssdata.SetDataToSPSS(coefsdataset, data.frame(row.names(res$coef), res$coef[,1]))
                },
                error=function(e) {print(e)
                cat(gettext("Failed to create coefficients dataset. Dataset name must not already exist: ",domain=domain),coefsdataset)
                }
            )  
        }
        
        spssdictionary.EndDataStep()
    }

    res <- tryCatch(rm(list=ls()),warning=function(e){return(NULL)})
}

StartProcedure<-function(procname, omsid){
if (as.integer(substr(spsspkg.GetSPSSVersion(),1, 2)) >= 19)
   spsspkg.StartProcedure(procname,omsid)
else
   spsspkg.StartProcedure(omsid)
}

caller<-function(dep, enter, distribution="gaussian", missing="listwise", lowerbound=NULL, upperbound=NULL,  
          coefsdataset=NULL, programfile=NULL, execute=TRUE){
 
    if(!is.null(programfile)){
        lines<-c("tobitreg<-",
            attr(tobitreg,"source"),
            paste("dep<-",dQuote(dep),sep=""),
            paste("enter<-",deparse(enter),sep=""),
            paste("distribution<-",dQuote(distribution),sep=""),
            paste("missing<-",dQuote(missing),sep=""))
        func<-"tobitreg(dep, enter, distribution, missing"
        if(!is.null(lowerbound)){
            func<-paste(func,", lowerbound=",lowerbound,sep="")
        }
        if(!is.null(upperbound)){
            func<-paste(func,", upperbound=",upperbound,sep="")
        }
        if(!is.null(coefsdataset)){
            func<-paste(func,", coefsdataset=",dQuote(coefsdataset),sep="")
        }
        func<-paste(func,")",sep="")
        lines<-c(lines,func)        
        f<-file(description=programfile,open="wt",encoding="UTF-8")
        writeLines(lines,con=f)
        close(f)
    }
    
    if (execute) tobitreg(dep, enter, distribution, missing, lowerbound, upperbound, coefsdataset)
    
}

setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 

Run<-function(args){
    cmdname = args[[1]]
    args <- args[[2]]
    oobj<-spsspkg.Syntax(templ=list(
                spsspkg.Template("DEPENDENT", subc="",  ktype="existingvarlist", var="dep", islist=FALSE),
                spsspkg.Template("ENTER", subc="",  ktype="existingvarlist", var="enter", islist=TRUE),
                spsspkg.Template("DISTRIBUTION", subc="",  ktype="str", var="distribution", islist=FALSE),
                spsspkg.Template("LOWERBOUND", subc="",  ktype="float", var="lowerbound", islist=FALSE),
                spsspkg.Template("UPPERBOUND", subc="",  ktype="float", var="upperbound", islist=FALSE),
                spsspkg.Template("MISSING", subc="OPTIONS",ktype="str", var="missing"),
                spsspkg.Template("COEFSDATASET", subc="SAVE", ktype="literal", var="coefsdataset"),
                spsspkg.Template("EXECUTE", subc="OPTIONS", ktype="bool", var="execute"),
                spsspkg.Template("PROGRAMFILE", subc="SAVE", ktype="literal", var="programfile")
                ))

    if ("HELP" %in% attr(args,"names"))
        #writeLines(helptext)
        helper(cmdname)
    else
        res <- spsspkg.processcmd(oobj,args,"caller")
}


helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}