# Portable core models. No patient data or machine-specific paths are embedded.
suppressPackageStartupMessages(library(lmerTest))

validate_data <- function(d) {
  columns <- c('PATNO','time','y','x','parent','y0','up0','age0','sex',
               'duration0','ledd0','genetic_group')
  stopifnot(all(columns %in% names(d)))
  numeric_columns <- setdiff(columns, c('PATNO','sex','genetic_group','parent'))
  stopifnot(all(vapply(d[numeric_columns], function(z) is.numeric(z) && all(is.finite(z)), logical(1))))
  stopifnot(nrow(d)>0, !anyNA(d[c('PATNO','sex','genetic_group')]),
            all(d$time>0), all(d$y>=0 & d$y<=52),
            !anyDuplicated(d[c('PATNO','time')]), all(table(d$PATNO)>=2))
  for (v in setdiff(columns,c('time','y'))) {
    stopifnot(all(vapply(split(d[[v]],d$PATNO),function(x) length(unique(x))==1,logical(1))))
  }
  d$PATNO <- factor(d$PATNO); d$sex <- factor(d$sex)
  d$genetic_group <- factor(d$genetic_group)
  stopifnot(nlevels(d$sex)>=2, nlevels(d$genetic_group)>=2)
  d$time_sq <- (d$time/5)^2
  d
}

model_formula <- function(parent=FALSE, structure='diagonal', exposure=TRUE) {
  stopifnot(structure %in% c('diagonal','correlated','quadratic'))
  terms <- c(if(exposure) 'x','y0','up0','age0','sex','duration0','ledd0','genetic_group',if(parent) 'parent')
  time <- if(structure=='quadratic') '(time+time_sq)' else 'time'
  random <- if(structure=='correlated') '(1+time|PATNO)' else '(1+time||PATNO)'
  as.formula(paste('y ~',time,'* (',paste(terms,collapse='+'),') +',random))
}

fit_model <- function(d, parent=FALSE, structure='diagonal', exposure=TRUE, reml=TRUE) {
  m <- lmer(model_formula(parent,structure,exposure),data=d,REML=reml,
            na.action=na.fail,
            control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000)))
  stopifnot(is.null(attr(getME(m,'X'),'col.dropped')),
            is.null(m@optinfo$conv$lme4$messages),!isSingular(m,tol=1e-4),
            all(is.finite(coef(summary(m)))))
  m
}

design <- function(m,d) {
  model.matrix(delete.response(terms(lme4::nobars(formula(m)))),d,
               contrasts.arg=attr(getME(m,'X'),'contrasts'))[,names(fixef(m)),drop=FALSE]
}

clinical_contrasts <- function(m,d,label) {
  people <- d[!duplicated(d$PATNO),,drop=FALSE]
  grid <- function(at_time, exposure_value=NULL) {
    g <- people
    g$time <- rep(at_time,nrow(g)); g$time_sq <- (g$time/5)^2
    if(!is.null(exposure_value)) g$x <- rep(exposure_value,nrow(g))
    g
  }
  rows <- list()
  for (h in c(1,3,5)) {
    L <- (colMeans(design(m,grid(h,.5)))-colMeans(design(m,grid(h,-.5))))-
         (colMeans(design(m,grid(0,.5)))-colMeans(design(m,grid(0,-.5))))
    # Independent analytic check catches the historical data.table name-shadow bug.
    expected <- setNames(rep(0,length(fixef(m))),names(fixef(m)))
    expected['time:x'] <- h
    if('time_sq:x' %in% names(expected)) expected['time_sq:x'] <- (h/5)^2
    stopifnot(max(abs(L-expected))<1e-10)
    cohort <- colMeans(design(m,grid(h)))-colMeans(design(m,grid(0)))
    for(kind in c('change_difference_per_1SD','cohort_standardized_change')) {
      vec <- if(kind=='change_difference_per_1SD') L else cohort
      a <- as.data.frame(contest1D(m,vec,confint=TRUE))
      a$model <- label; a$contrast <- kind; a$year <- h
      rows[[length(rows)+1]] <- a
    }
  }
  do.call(rbind,rows)
}

model_metrics <- function(m,d,label) {
  vc <- as.data.frame(VarCorr(m)); vr <- 0
  for(i in seq_len(nrow(vc))) {
    a <- vc[i,]; if(a$grp=='Residual') next
    z1 <- if(a$var1=='(Intercept)') rep(1,nrow(d)) else d[[a$var1]]
    if(is.na(a$var2)) vr <- vr+a$vcov*mean(z1^2) else {
      z2 <- if(a$var2=='(Intercept)') rep(1,nrow(d)) else d[[a$var2]]
      vr <- vr+2*a$vcov*mean(z1*z2)
    }
  }
  vf <- var(predict(m,re.form=NA)); ve <- sigma(m)^2
  data.frame(model=label,REML=isREML(m),people=nlevels(droplevels(d$PATNO)),
             observations=nrow(d),AIC=AIC(m),BIC=BIC(m),logLik=as.numeric(logLik(m)),
             marginal_R2=vf/(vf+vr+ve),conditional_R2=(vf+vr)/(vf+vr+ve))
}
