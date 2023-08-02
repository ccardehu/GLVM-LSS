# For PISA 2018 application
# ~~~~~~~~~~~~~~~~~~~~~~~~~

q1D <- function(x, quant, item, model){
  mu <- (model$b$mu[item,1] + model$b$mu[item,3]*x)
  sg <- exp(model$b$si[item,1] + model$b$si[item,3]*x)
  nu <- plogis(model$b$nu[item,1] + model$b$nu[item,3]*x)
  
  p2nuc <- function(p){
    g1max <- sqrt(2)*(4-pi)/(pi-2)^(3/2)
    nuc <- (p - g1max/(2*g1max))*(2*g1max)
    nuc
  }
  mc2md <- function(muc,sgc,nuc){
    b <- sqrt(2/pi)
    nud <- nc2nd(nuc)
    d <- nud/sqrt(1+nud^2)
    sgd <- sc2sd(sgc,nuc)
    mud <- muc - b*sgd*d
    mud
  }
  sc2sd <- function(sgc,nuc){
    b <- sqrt(2/pi)
    nud <- nc2nd(nuc)
    d <- nud/sqrt(1+nud^2)
    sgd <- sgc/sqrt(1-b^2*d^2)
    sgd
  }
  nc2nd <- function(nuc){
    b <- sqrt(2/pi)
    R <- ((2*abs(nuc))/(4-pi))^(1/3) * sign(nuc)
    nud <- R/sqrt(b^2-(1-b^2)*R^2)
    nud 
  }
  
  nuc = p2nuc(nu)
  nud = nc2nd(nuc = nuc)
  sgd = sc2sd(sgc = sg, nuc = nuc)
  mud = mc2md(muc = mu, sgc = sg, nuc = nuc)
  
  sapply(1:length(x), function(i) sn::qsn(quant, xi = mud[i], omega = sgd[i], alpha = nud[i]))
}

gg_Dist <- function(item){
  texttitle <- paste0("Item ",i, " $(log(t_",i,") )$")
  p1 <- ggplot(data, aes(y = get(paste0("T",item)))) + 
    # ylim(min(data[,paste0("T",item)]), max(data[,paste0("T",item)])) +
    scale_x_continuous(limit = c(-3,3), breaks = round(seq(-3,3,by = 1),1)) +
    theme_classic() +
    xlab(TeX(r"(Latent speed factor ($z_2$))")) + 
    ylab(TeX(r"(Response time (log-minutes))")) +
    ggtitle(TeX(texttitle)) +
    geom_function(fun = ~ mod_SN_comp$b$mu[paste0("T",item),1] + mod_SN_comp$b$mu[paste0("T",item),3]*.x, colour = "black", linewidth = 1) +
    geom_function(fun = q1D, args = list(quant = 0.025, item = paste0("T",item), model = mod_SN_comp), colour = "black", linetype="dotted", linewidth = 0.7) +
    geom_function(fun = q1D, args = list(quant = 0.10, item = paste0("T",item), model = mod_SN_comp), colour = "black", linetype="dotted", linewidth = 0.7) +
    geom_function(fun = q1D, args = list(quant = 0.25, item = paste0("T",item), model = mod_SN_comp), colour = "black", linetype="dotted", linewidth = 0.7) + 
    geom_function(fun = q1D, args = list(quant = 0.50, item = paste0("T",item), model = mod_SN_comp), colour = "black", linetype="longdash", linewidth = 0.7) + 
    geom_function(fun = q1D, args = list(quant = 0.75, item = paste0("T",item), model = mod_SN_comp), colour = "black", linetype="dotted", linewidth = 0.7) + 
    geom_function(fun = q1D, args = list(quant = 0.90, item = paste0("T",item), model = mod_SN_comp), colour = "black", linetype="dotted", linewidth = 0.7) + 
    geom_function(fun = q1D, args = list(quant = 0.975, item = paste0("T",item), model = mod_SN_comp), colour = "black", linetype="dotted", linewidth = 0.7) + 
    theme(plot.title = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0), size = 15),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 13),
          axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0), size = 13),
          axis.title =   element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.margin = margin(t = 15, r = 15, b = 15, l = 15)) + 
    annotation_custom(text_high,xmin=3, xmax=3, ymin= ylablim[item], ymax=ylablim[item]) +
    annotation_custom(text_low,xmin=-3, xmax=-3,ymin= ylablim[item], ymax= ylablim[item]) +
    coord_cartesian(clip="off")
  p1
}

gg_Prob <- function(item){
  texttitle <- paste0("Item ",i, " $(y_",i,")$")
  axistitle <- paste0("$P(y_",item," = 1 \\,|\\, z_1)$")
  p1 <- ggplot(data, aes(y = get(paste0("Y",item)))) + 
    ylim(0,1) +
    scale_x_continuous(limit = c(-3,3), breaks = round(seq(-3,3,by = 1),1)) +
    theme_classic() +
    xlab(TeX(r"(Latent ability factor ($z_1$))")) + 
    ylab(TeX(axistitle)) +
    ggtitle(TeX(texttitle)) +
    geom_function(fun = ~ plogis(mod_SN_comp$b$mu[paste0("Y",item),1] + mod_SN_comp$b$mu[paste0("Y",item),2]*.x), colour = "black", linewidth = 1) +
    theme(plot.title = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0), size = 15),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 13),
          axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0), size = 13),
          axis.title =   element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.margin = margin(t = 15, r = 15, b = 15, l = 15))
  p1
}

margPlot <- function(item,m1,m2,m3){
  
  p2nuc <- function(p){
    g1max <- sqrt(2)*(4-pi)/(pi-2)^(3/2)
    nuc <- (p - g1max/(2*g1max))*(2*g1max)
    nuc
  }
  mc2md <- function(muc,sgc,nuc){
    b <- sqrt(2/pi)
    nud <- nc2nd(nuc)
    d <- nud/sqrt(1+nud^2)
    sgd <- sc2sd(sgc,nuc)
    mud <- muc - b*sgd*d
    mud
  }
  sc2sd <- function(sgc,nuc){
    b <- sqrt(2/pi)
    nud <- nc2nd(nuc)
    d <- nud/sqrt(1+nud^2)
    sgd <- sgc/sqrt(1-b^2*d^2)
    sgd
  }
  nc2nd <- function(nuc){
    b <- sqrt(2/pi)
    R <- ((2*abs(nuc))/(4-pi))^(1/3) * sign(nuc)
    nud <- R/sqrt(b^2-(1-b^2)*R^2)
    nud 
  }
  
  mod1 <- m1
  form1 <- prep_form(mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, nu.eq = ~ Z1+Z2, ta.eq = NULL)
  ghQ1 <- prep_ghq(nQP = 25, form = form1, Rz = mod1$Rz)
  
  mod2 <- m2
  form2 <- prep_form(mu.eq = ~ Z1+Z2, sg.eq = ~ Z1+Z2, nu.eq = NULL, ta.eq = NULL)
  ghQ2 <- prep_ghq(nQP = 25, form = form2, Rz = mod2$Rz)
  
  mod3 <- m3
  form3 <- prep_form(mu.eq = ~ Z1+Z2, sg.eq = ~ 1, nu.eq = NULL, ta.eq = NULL)
  ghQ3 <- prep_ghq(nQP = 25, form = form3, Rz = mod3$Rz)
  
  muc = as.matrix(ghQ1$out$mu)%*%mod1$b$mu[item,]
  sgc = exp(as.matrix(ghQ1$out$si)%*%mod1$b$si[item,])
  nup = plogis(as.matrix(ghQ1$out$nu)%*%mod1$b$nu[item,])
  
  nuc = p2nuc(nup)
  nud = nc2nd(nuc = nuc)
  sgd = sc2sd(sgc = sgc, nuc = nuc)
  mud = mc2md(muc = muc, sgc = sgc, nuc = nuc)
  
  mu2 = as.matrix(ghQ2$out$mu)%*%mod2$b$mu[item,]
  sg2 = exp(as.matrix(ghQ2$out$si)%*%mod2$b$si[item,])
  
  mu3 = as.matrix(ghQ3$out$mu)%*%mod3$b$mu[item,]
  sg3 = exp(as.matrix(ghQ3$out$si)%*%mod3$b$si[item,])
  
  rtime <- seq(from = min(data[,1:9]), to = max(data[,1:9]), length.out = 100)
  
  dens1 <- rowSums(sapply(1:length(ghQ1$weights), function(i) sn::dsn(x = rtime, xi = mud[i], omega = sgd[i], alpha = nud[i])*ghQ1$weights[i]))
  dens2 <- rowSums(sapply(1:length(ghQ2$weights), function(i) dnorm(x = rtime, mean = mu2[i], sd = sg2[i])*ghQ2$weights[i]))
  dens3 <- rowSums(sapply(1:length(ghQ3$weights), function(i) dnorm(x = rtime, mean = mu3[i], sd = sg3[i])*ghQ3$weights[i]))
  
  return(data.frame(x = rtime, d1 = dens1, d2 = dens2, d3 = dens3))
  
}

# For ANES 2020 application 
# ~~~~~~~~~~~~~~~~~~~~~~~~~

q0De2 <- function(x, quant, item, model){
  mu <- plogis(model$b$mu[item,1] + model$b$mu[item,2]*x)
  sg <- plogis(model$b$si[item,1])
  shape1 = mu*(1-sg^2)/(sg^2)
  shape2 = (1-mu)*(1-sg^2)/(sg^2)
  return(qbeta(quant, shape1, shape2))
}

q1De2 <- function(x, quant, item, model){
  mu <- plogis(model$b$mu[item,1] + model$b$mu[item,2]*x)
  sg <- plogis(model$b$si[item,1] + model$b$si[item,2]*x)
  shape1 = mu*(1-sg^2)/(sg^2)
  shape2 = (1-mu)*(1-sg^2)/(sg^2)
  return(qbeta(quant, shape1, shape2))
}

gg_qq_empirical <- function(a, b, quantiles = seq(0, 1, 0.01)){
  a_lab <- deparse(substitute(a))
  if(missing(b)) {
    b <- rnorm(length(a), mean(a), sd(a))
    b_lab <- "normal distribution"
  }
  else b_lab <- deparse(substitute(b))
  
  ggplot(mapping = aes(x = quantile(a, quantiles, na.rm = T), 
                       y = quantile(scale(b), quantiles, na.rm = T))) + 
    geom_point() +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2)
}

gg_q_data <- function(data){
  ggplot(mapping = aes(data)) + 
  stat_ecdf(geom = "step", size = 1, na.rm = T, alpha = 1, colour = "black", pad = T) + 
  ylab(TeX(r"(Empirical percentiles (CDF))")) + 
  xlab(TeX(r"(Thermometer rating (scaled))"))
}

margPlotANES <- function(item,m1,m2,m3,m4){
  
  mod1 <- m1 # Normal, homoscedastic
  form1 <- prep_form(mu.eq = ~ Z1, sg.eq = ~ 1, nu.eq = NULL, ta.eq = NULL)
  ghQ1 <- prep_ghq(nQP = 100, form = form1, Rz = mod1$Rz)
  
  mod2 <- m2 # Normal, heteroscedastic
  form2 <- prep_form(mu.eq = ~ Z1, sg.eq = ~ Z1, nu.eq = NULL, ta.eq = NULL)
  ghQ2 <- prep_ghq(nQP = 100, form = form2, Rz = mod2$Rz)
  
  mod3 <- m3 # Beta, homoscedastic
  form3 <- prep_form(mu.eq = ~ Z1, sg.eq = ~ 1, nu.eq = NULL, ta.eq = NULL)
  ghQ3 <- prep_ghq(nQP = 100, form = form3, Rz = mod3$Rz)
  
  mod4 <- m4 # Beta, heteroscedastic
  form4 <- prep_form(mu.eq = ~ Z1, sg.eq = ~ Z1, nu.eq = NULL, ta.eq = NULL)
  ghQ4 <- prep_ghq(nQP = 100, form = form4, Rz = mod3$Rz)
  
  mu1 = as.matrix(ghQ1$out$mu)%*%mod1$b$mu[item,]
  sg1 = exp(as.matrix(ghQ1$out$si)%*%mod1$b$si[item,])
  
  mu2 = as.matrix(ghQ2$out$mu)%*%mod2$b$mu[item,]
  sg2 = exp(as.matrix(ghQ2$out$si)%*%mod2$b$si[item,])
  
  mu3 = plogis(as.matrix(ghQ3$out$mu)%*%mod3$b$mu[item,])
  sg3 = plogis(as.matrix(ghQ3$out$si)%*%mod3$b$si[item,])
  mu3d = mu3*(1-sg3^2)/sg3^2
  sg3d = (1-mu3)*(1-sg3^2)/sg3^2
  
  mu4 = plogis(as.matrix(ghQ4$out$mu)%*%mod4$b$mu[item,])
  sg4 = plogis(as.matrix(ghQ4$out$si)%*%mod4$b$si[item,])
  mu4d = mu4*(1-sg4^2)/sg4^2
  sg4d = (1-mu4)*(1-sg4^2)/sg4^2
  
  rtime <- seq(from = 0+0.001, to = 1-0.001, length.out = 100)
  
  dens1 <- rowSums(sapply(1:length(ghQ1$weights), function(i) dnorm(x = rtime, mean = mu1[i], sd = sg1[i])*ghQ1$weights[i]))
  dens2 <- rowSums(sapply(1:length(ghQ2$weights), function(i) dnorm(x = rtime, mean = mu2[i], sd = sg2[i])*ghQ2$weights[i]))
  dens3 <- rowSums(sapply(1:length(ghQ3$weights), function(i) dbeta(x = rtime, shape1 = mu3d[i], shape2 = sg3d[i])*ghQ3$weights[i]))
  dens4 <- rowSums(sapply(1:length(ghQ4$weights), function(i) dbeta(x = rtime, shape1 = mu4d[i], shape2 = sg4d[i])*ghQ4$weights[i]))
  
  return(data.frame(x = rtime, d1 = dens1, d2 = dens2, d3 = dens3, d4 = dens4))
  
}

gg_Dist2 <- function(item, mod){
  texttitle <- paste0("Item ",i, " $(log(t_",i,") )$")
  p1 <- ggplot(data, aes(y = get(paste0("T",item)))) + 
    # ylim(min(data[,paste0("T",item)]), max(data[,paste0("T",item)])) +
    scale_x_continuous(limit = c(-3,3), breaks = round(seq(-3,3,by = 1),1)) +
    theme_classic() +
    xlab(TeX(r"(Latent speed factor ($z_2$))")) + 
    ylab(TeX(r"(Response time (log-minutes))")) +
    ggtitle(TeX(texttitle)) +
    geom_function(fun = ~ mod$b$mu[paste0("T",item),1] + mod$b$mu[paste0("T",item),3]*.x, colour = "black", linewidth = 1) +
    geom_function(fun = q1D, args = list(quant = 0.025, item = paste0("T",item), model = mod), colour = "black", linetype="dotted", linewidth = 0.7) +
    geom_function(fun = q1D, args = list(quant = 0.10, item = paste0("T",item), model = mod), colour = "black", linetype="dotted", linewidth = 0.7) +
    geom_function(fun = q1D, args = list(quant = 0.25, item = paste0("T",item), model = mod), colour = "black", linetype="dotted", linewidth = 0.7) + 
    geom_function(fun = q1D, args = list(quant = 0.50, item = paste0("T",item), model = mod), colour = "black", linetype="longdash", linewidth = 0.7) + 
    geom_function(fun = q1D, args = list(quant = 0.75, item = paste0("T",item), model = mod), colour = "black", linetype="dotted", linewidth = 0.7) + 
    geom_function(fun = q1D, args = list(quant = 0.90, item = paste0("T",item), model = mod), colour = "black", linetype="dotted", linewidth = 0.7) + 
    geom_function(fun = q1D, args = list(quant = 0.975, item = paste0("T",item), model = mod), colour = "black", linetype="dotted", linewidth = 0.7) + 
    theme(plot.title = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0), size = 15),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 13),
          axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0), size = 13),
          axis.title =   element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.margin = margin(t = 15, r = 15, b = 15, l = 15)) + 
    annotation_custom(text_high,xmin=3, xmax=3, ymin= ylablim[item], ymax=ylablim[item]) +
    annotation_custom(text_low,xmin=-3, xmax=-3,ymin= ylablim[item], ymax= ylablim[item]) +
    coord_cartesian(clip="off")
  p1
}

gg_Prob2 <- function(item, mod){
  texttitle <- paste0("Item ",i, " $(y_",i,")$")
  axistitle <- paste0("$P(y_",item," = 1 \\,|\\, z_1)$")
  p1 <- ggplot(data, aes(y = get(paste0("Y",item)))) + 
    ylim(0,1) +
    scale_x_continuous(limit = c(-3,3), breaks = round(seq(-3,3,by = 1),1)) +
    theme_classic() +
    xlab(TeX(r"(Latent ability factor ($z_1$))")) + 
    ylab(TeX(axistitle)) +
    ggtitle(TeX(texttitle)) +
    geom_function(fun = ~ plogis(mod$b$mu[paste0("Y",item),1] + mod$b$mu[paste0("Y",item),2]*.x), colour = "black", linewidth = 1) +
    theme(plot.title = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0), size = 15),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 13),
          axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0), size = 13),
          axis.title =   element_text(margin = margin(t = 0, r = 0, b = 10, l = 0)),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.margin = margin(t = 15, r = 15, b = 15, l = 15))
  p1
}