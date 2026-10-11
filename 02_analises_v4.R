# =============================================================================
# 02_analises_v4.R
# Todas as análises da versão 4 do artigo: não atuação (Tabela 4), modelos
# multinível (Tabela 5), correção para múltiplos testes, robustez (Tabela A3),
# validação da preparação, exposição por período (Tabela A2), interações,
# Heckman e descritivos. As saídas vão para saida_analises_v4.txt.
# Pacotes: lme4, lmerTest, clubSandwich, sandwich, lmtest, glmmTMB.
# Leva alguns minutos por causa dos modelos binomiais com efeito por observação.
# =============================================================================
suppressMessages({library(lme4); library(lmerTest); library(clubSandwich); library(sandwich); library(lmtest); library(glmmTMB)})
base <- read.csv("base_suas_desastres_2025.csv", fileEncoding="UTF-8", stringsAsFactors=FALSE)
base$regiao <- factor(base$regiao, levels=c("Sudeste","Centro Oeste","Nordeste","Norte","Sul"))
base$porte <- factor(base$porte, levels=c("Pequeno I","Pequeno II","Médio","Grande","Metrópole"))
base$uf <- factor(base$uf)
z <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
base <- within(base, {
  vigil <- an_vigilancia_formal; diag <- an_diagnostico; sistema <- an_sistema_informacao
  estimativa <- an_info_federal_risco; dc_mapeam <- dc_mapeamento_2020
  capacit <- op_capacitacao_2025; dirig_sup <- op_dirigente_superior; secret <- op_secretaria_exclusiva
  rh10 <- z(log1p(rh_por_10mil_hab)); estat <- z(rh_prop_estatutarios); gasto <- z(log1p(op_gasto_as_pc))
  pmas <- op_pmas; compdec <- dc_compdec_2020; dc_anos <- dc_anos_gasto; cofin <- sis_cofin_estadual_calamidade
  consorcio <- op_consorcio; hab <- hab_prioriza_risco_2024
  conselho <- z(pol_conselho_indice); ouvid <- pol_ouvidoria; cib <- pol_acompanha_cib
  alin <- pol_alinhamento_governador; troca <- pol_troca_prefeito
  ameaca <- z(risco_ameaca_max); estresse <- z(risco_estresse_hidr)
  hid <- z(log1p(atl_rec_nhid_2020_2024)); seca <- z(log1p(atl_rec_nseca_2020_2024))
  hid_reg <- z(log1p(atl_nhid_2020_2024)); seca_reg <- z(log1p(atl_nseca_2020_2024))
  hid_1319 <- z(log1p(atl_rec_nhid_2013_2019)); seca_1319 <- z(log1p(atl_rec_nseca_2013_2019))
  hid_25 <- z(log1p(atl_rec_nhid_2025)); seca_25 <- z(log1p(atl_rec_nseca_2025))
  lista <- cc_lista; vuln <- z(soc_vulnerab_renda); logpop <- z(log(populacao_2022))
  idx_provisao4 <- idx_provisao
  idx_provisao <- rowMeans(cbind(acao_beneficio_eventual, acao_cestas_basicas, acao_donativos)) * 100
  dif <- idx_provisao - idx_estrategica
  PeqI <- as.numeric(porte == "Pequeno I")
})
CAP <- c("vigil","diag","sistema","estimativa","capacit","dirig_sup","secret","rh10","estat","gasto","pmas","consorcio","cofin",
         "conselho","ouvid","cib","compdec","dc_mapeam","dc_anos","hab","alin","troca")
CTX <- c("hid","seca","ameaca","estresse","lista","vuln","logpop")
EVT <- c("atuou_regular","evento_enchente_inundacao","evento_seca_estiagem")
A <- subset(base, atuou_2025 == 1)
cc <- function(d, v) { d <- d[complete.cases(d[, v]), ]; d$uf <- droplevels(d$uf); d }
st <- function(p) ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*","")))
sink("saida_analises_v4.txt")

## ---------- 1. Etapa 1: não atuação ----------
L <- subset(base, rec25_algum == 1); L$nao <- 1 - L$atuou_2025; L$desab <- as.numeric(L$rec25_desabrigados > 0)
L$q4 <- as.numeric(L$rec25_mes_primeiro >= 10)
cat("Municípios com reconhecimento 2025:", nrow(L), " % não atuou:", round(100*mean(L$nao),1), "\n")
print(round(tapply(L$nao, L$rec25_so_seca, mean)*100,1)); print(table(L$rec25_so_seca)); print(table(L$rec25_so_seca, L$desab))
cat("só seca no RS:", sum(L$rec25_so_seca==1 & L$uf=="RS"), " no NE:", sum(L$rec25_so_seca==1 & L$regiao=="Nordeste"), "\n")
XL <- c("rec25_so_seca", CAP, CTX, "q4")
logit_tab <- function(d, X, lab) {
  m <- glm(reformulate(c(X, "regiao"), "nao"), data=d, family=binomial)
  c2 <- coef_test(m, vcov="CR2", cluster=d$uf, test="Satterthwaite")
  V2 <- vcovCR(m, cluster=d$uf, type="CR2")
  se <- sqrt(diag(as.matrix(V2)))
  o <- data.frame(OR=exp(coef(m)), lo=exp(coef(m)-1.96*se), hi=exp(coef(m)+1.96*se), p=c2$p_Satt)
  cat("\n####", lab, " n=", nrow(d), " G=", nlevels(d$uf), " McFadden=", round(1-m$deviance/m$null.deviance,3), "\n")
  print(round(o,3)); invisible(list(m=m, o=o))
}
L1 <- cc(L, c("nao", XL, "regiao", "uf", "desab"))
r1 <- logit_tab(L1, XL, "LOGIT sem desabrigados")
r2 <- logit_tab(L1, c(XL, "desab"), "LOGIT com desabrigados")
# probabilidades preditas médias (AME) so_seca, Nordeste
m <- r1$m; d0 <- L1; d1 <- L1; d0$rec25_so_seca <- 0; d1$rec25_so_seca <- 1
cat("\nAME só seca (pp):", round(100*mean(predict(m,d1,type="response")-predict(m,d0,type="response")),1), "\n")
ne <- L1[L1$regiao=="Nordeste",]; ne0<-ne; ne1<-ne; ne0$rec25_so_seca<-0; ne1$rec25_so_seca<-1
cat("NE preditos médios: demais", round(100*mean(predict(m,ne0,type="response")),1), " só seca", round(100*mean(predict(m,ne1,type="response")),1), "\n")
# EF de UF
mfe <- glm(reformulate(c(XL, "uf"), "nao"), data=L1, family=binomial)
cfe <- coef_test(mfe, vcov="CR2", cluster=L1$uf, test="Satterthwaite")
cat("EF UF so_seca OR", round(exp(coef(mfe)["rec25_so_seca"]),2), " p", round(cfe$p_Satt[names(coef(mfe))=="rec25_so_seca"],3), "\n")
# vigilância × só seca
mv <- glm(reformulate(c(XL, "regiao", "vigil:rec25_so_seca"), "nao"), data=L1, family=binomial)
Vv <- as.matrix(vcovCR(mv, cluster=L1$uf, type="CR2")); b <- coef(mv)
k <- grep("vigil:rec25_so_seca|rec25_so_seca:vigil", names(b), value=TRUE); bs <- b["vigil"]+b[k]; ses <- sqrt(Vv["vigil","vigil"]+Vv[k,k]+2*Vv["vigil",k])
cv <- coef_test(mv, vcov="CR2", cluster=L1$uf, test="Satterthwaite")
cat("Vigilância: hidro OR", round(exp(b["vigil"]),2), " p", round(cv$p_Satt[names(b)=="vigil"],3), "| seca OR", round(exp(bs),2), " p(z)", round(2*pnorm(-abs(bs/ses)),3), "| interação OR", round(exp(b[k]),2), " p", round(cv$p_Satt[names(b)==k],3), "\n")
# histórico de secas reconhecidas 2013-19 × só seca
mh <- glm(reformulate(c(XL, "hid_1319", "seca_1319", "regiao", "seca_1319:rec25_so_seca"), "nao"), data=L1, family=binomial)
Vh <- as.matrix(vcovCR(mh, cluster=L1$uf, type="CR2")); b <- coef(mh); k <- "rec25_so_seca:seca_1319"; if(!k %in% names(b)) k <- "seca_1319:rec25_so_seca"
ch <- coef_test(mh, vcov="CR2", cluster=L1$uf, test="Satterthwaite")
bs <- b["seca_1319"]+b[k]; ses <- sqrt(Vh["seca_1319","seca_1319"]+Vh[k,k]+2*Vh["seca_1319",k])
cat("Histórico seca: interação OR", round(exp(b[k]),2), " p", round(ch$p_Satt[names(b)==k],3), "| só seca OR", round(exp(bs),2), " p(z)", round(2*pnorm(-abs(bs/ses)),3),
    "| demais OR", round(exp(b["seca_1319"]),2), " p", round(ch$p_Satt[names(b)=="seca_1319"],3), "| hid_1319 OR", round(exp(b["hid_1319"]),2), " p", round(ch$p_Satt[names(b)=="hid_1319"],3), "\n")
# robustez: registros em vez de reconhecidos
XLr <- c("rec25_so_seca", CAP, "hid_reg","seca_reg","ameaca","estresse","lista","vuln","logpop","q4")
L2 <- cc(L, c("nao", XLr, "regiao", "uf")); mr <- glm(reformulate(c(XLr,"regiao"),"nao"), data=L2, family=binomial)
cr <- coef_test(mr, vcov="CR2", cluster=L2$uf, test="Satterthwaite"); cat("Registros: so_seca OR", round(exp(coef(mr)["rec25_so_seca"]),2), " p", round(cr$p_Satt[2],3), "\n")
# sem RS
L3 <- droplevels(L1[L1$uf!="RS",]); m3 <- glm(reformulate(c(XL,"regiao"),"nao"), data=L3, family=binomial)
c3 <- coef_test(m3, vcov="CR2", cluster=L3$uf, test="Satterthwaite"); cat("Sem RS: so_seca OR", round(exp(coef(m3)["rec25_so_seca"]),2), " p", round(c3$p_Satt[2],3), "\n")

## ---------- 2. Etapa 2: modelos multinível ----------
mlm <- function(y, d, X) {
  d <- cc(d, c(y, X, "regiao", "uf"))
  m <- lmer(reformulate(c(X, "regiao", "(1|uf)"), y), data=d, REML=TRUE)
  m0 <- lmer(reformulate("(1|uf)", y), data=d, REML=TRUE); vc <- as.data.frame(VarCorr(m0))$vcov
  s <- summary(m)$coefficients
  list(m=m, n=nrow(d), icc=vc[1]/sum(vc), tab=data.frame(b=s[,1], se=s[,2], df=s[,"df"], p=s[,5]), d=d)
}
M <- list(prep = mlm("idx_preparacao", base, c(CAP, CTX)),
          prov = mlm("idx_provisao", A, c(CAP, CTX, EVT)),
          plan = mlm("idx_estrategica", A, c(CAP, CTX, EVT)),
          dif  = mlm("dif", A, c(CAP, CTX, EVT)))
for (k in names(M)) { cat("\n####", k, " n=", M[[k]]$n, " ICC nulo=", round(M[[k]]$icc,3), " var UF=", round(as.data.frame(VarCorr(M[[k]]$m))$vcov[1],2), "\n"); print(round(M[[k]]$tab,3)) }
# FDR sobre os coeficientes de CAP+CTX(+EVT) dos três desfechos substantivos
P <- do.call(rbind, lapply(c("prep","prov","plan"), function(k) { t <- M[[k]]$tab; v <- intersect(c(CAP,CTX), rownames(t)); data.frame(mod=k, var=v, p=t[v,"p"]) }))
P$q <- p.adjust(P$p, "BH"); cat("\nFDR (BH) q<0,05 entre", nrow(P), "testes:", sum(P$q<.05), " p<0,05:", sum(P$p<.05), "\n"); print(P[P$p<.05,], row.names=FALSE)
write.csv(P, "fdr.csv", row.names=FALSE)
saveRDS(lapply(M, function(x) list(tab=x$tab, n=x$n, icc=x$icc)), "mlm_v4.rds")
saveRDS(list(r1=r1$o, r2=r2$o), "logit_v4.rds")

## ---------- 3. Robustez ----------
KEY <- c("cofin","estimativa","conselho","ouvid","capacit","diag","vigil","cib","pmas","consorcio","gasto","estat","secret","hid","seca","lista","hab")
rob <- list()
for (k in c("prep","prov","plan")) {
  y <- c(prep="idx_preparacao", prov="idx_provisao", plan="idx_estrategica")[k]; dat <- if (k=="prep") base else A
  X <- if (k=="prep") c(CAP,CTX) else c(CAP,CTX,EVT)
  d <- cc(dat, c(y, X, "uf"))
  fe <- lm(reformulate(c(X,"uf"), y), data=d); cf <- coef_test(fe, vcov="CR2", cluster=d$uf, test="Satterthwaite")
  rr <- data.frame(EF_CR2=paste0(round(coef(fe)[KEY],1), st(cf$p_Satt[match(KEY,names(coef(fe)))])))
  dn <- droplevels(d[d$uf!="RS",]); ms <- lmer(reformulate(c(X,"regiao","(1|uf)"), y), data=dn)
  ss <- summary(ms)$coefficients; rr$semRS <- paste0(round(ss[KEY,1],1), st(ss[KEY,5]))
  if (k!="prep") { Xn <- c(CAP,CTX); mn <- lmer(reformulate(c(Xn,"regiao","(1|uf)"), y), data=cc(dat,c(y,Xn,"uf","regiao")))
    sn <- summary(mn)$coefficients; rr$semEVT <- paste0(round(sn[KEY,1],1), st(sn[KEY,5])) }
  Xr <- sub("^hid$","hid_reg", sub("^seca$","seca_reg", X)); mrg <- lmer(reformulate(c(Xr,"regiao","(1|uf)"), y), data=cc(dat,c(y,Xr,"uf","regiao")))
  sr <- summary(mrg)$coefficients; KR <- sub("^hid$","hid_reg", sub("^seca$","seca_reg", KEY)); rr$registros <- paste0(round(sr[KR,1],1), st(sr[KR,5]))
  # GLMM binomial com efeito em nível de observação (OLRE)
  nit <- if (k=="prep") 4 else if (k=="prov") 3 else 7
  d$succ <- round(d[[y]]/100*nit); d$fail <- nit - d$succ; d$obs <- factor(seq_len(nrow(d)))
  g <- glmer(as.formula(paste("cbind(succ,fail)~", paste(c(X,"regiao"), collapse="+"), "+(1|uf)+(1|obs)")), data=d, family=binomial,
             control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  sg <- summary(g)$coefficients; rr$GLMM_OR <- paste0(round(exp(sg[KEY,1]),2), st(sg[KEY,4]))
  rownames(rr) <- KEY; rob[[k]] <- rr; cat("\n#### ROBUSTEZ", k, "\n"); print(rr)
}
saveRDS(rob, "rob_v4.rds")

## ---------- 4. Preparação: validação ----------
nao <- subset(base, atuou_2025 == 0)
mn <- mlm("idx_preparacao", nao, c(CAP, CTX)); cat("\n#### PREPARAÇÃO entre não atuantes n=", mn$n, "\n"); print(round(mn$tab[KEY,],3))
cat("Correlações preparação × MUNIC 2020 (DC):\n")
for (v in c("dc_compdec_2020","dc_plano_contingencia_2020","dc_alerta_2020","dc_mapeamento_2020")) cat(v, round(cor(base$idx_preparacao, base[[v]], use="complete.obs"),3), "\n")
cat("Média preparação por plano de contingência DC 2020:", round(tapply(base$idx_preparacao, base$dc_plano_contingencia_2020, mean),1), "\n")

## ---------- 5. Exposição por período (reconhecidos) ----------
EXP <- c("hid_1319","hid","hid_25","seca_1319","seca","seca_25")
for (k in c("prep","prov","plan")) {
  y <- c(prep="idx_preparacao", prov="idx_provisao", plan="idx_estrategica")[k]; dat <- if (k=="prep") base else A
  X <- c(CAP, setdiff(CTX,c("hid","seca")), EXP, if (k!="prep") EVT)
  r <- mlm(y, dat, X); cat("\n#### PERÍODOS", k, "\n"); print(round(r$tab[EXP,],3))
}
cat("cor hid_1319 × hid:", round(cor(base$hid_1319, base$hid),2), " seca:", round(cor(base$seca_1319, base$seca),2), "\n")

## ---------- 6. Interações ----------
binint <- function(y, dat, X, a, bv, nit) {
  d <- cc(dat, c(y, X, a, bv, "uf")); d$succ <- round(d[[y]]/100*nit); d$fail <- nit-d$succ; d$int <- d[[a]]*d[[bv]]
  g <- glm(as.formula(paste("cbind(succ,fail)~int+", paste(unique(c(X,a,bv,"uf")), collapse="+"))), data=d, family=quasibinomial)
  cg <- coef_test(g, vcov="CR2", cluster=d$uf, test="Satterthwaite"); j <- which(names(coef(g))=="int")
  lin <- lm(as.formula(paste(y,"~int+", paste(unique(c(X,a,bv,"uf")), collapse="+"))), data=d); cl <- coef_test(lin, vcov="CR2", cluster=d$uf, test="Satterthwaite")
  cat(sprintf("%s %s×%s: linear b=%.2f p=%.3f | binomial OR=%.2f p=%.3f\n", y, a, bv, coef(lin)["int"], cl$p_Satt[2], exp(coef(g)[j]), cg$p_Satt[j]))
}
cat("\n#### INTERAÇÕES (EF de UF, CR2)\n")
for (cv in c("conselho","cofin","capacit","ouvid","diag","vigil")) binint("idx_preparacao", base, c(CAP,CTX), "hid", cv, 4)
for (cv in c("conselho","cofin","capacit","ouvid")) binint("idx_estrategica", A, c(CAP,CTX,EVT), "hid", cv, 7)
binint("idx_estrategica", A, c(CAP,CTX,EVT), "cofin", "lista", 7)
binint("idx_estrategica", A, c(CAP,CTX,EVT), "conselho", "cofin", 7)
binint("idx_preparacao", base, c(CAP,CTX,"PeqI"), "cofin", "PeqI", 4)
# cofinanciamento × lista: efeitos condicionais (multinível)
d <- cc(A, c("idx_estrategica", CAP, CTX, EVT, "regiao", "uf")); mi <- lmer(reformulate(c(CAP,CTX,EVT,"regiao","cofin:lista","(1|uf)"),"idx_estrategica"), data=d)
b <- fixef(mi); V <- as.matrix(vcov(mi)); k <- grep("cofin:lista", names(b), value=TRUE)
cat("cofin|lista=0:", round(b["cofin"],1), " cofin|lista=1:", round(b["cofin"]+b[k],1), " int:", round(b[k],1), " p", round(summary(mi)$coefficients[k,5],3), "\n")

## ---------- 7. Heckman (seleção para atuação) ----------
base$rec25_algum0 <- ifelse(is.na(base$rec25_algum),0,base$rec25_algum)
base$n_rec25 <- log1p(base$atl_rec_nhid_2025 + base$atl_rec_nseca_2025)
dS <- cc(base, c("atuou_2025", CAP, CTX, "rec25_algum0", "n_rec25", "uf"))
pr <- glm(reformulate(c(CAP, CTX, "rec25_algum0", "n_rec25", "uf"), "atuou_2025"), data=dS, family=binomial(link="probit"))
cat("\n#### HECKMAN: exclusão rec25 z=", round(summary(pr)$coefficients["rec25_algum0",3],1), " n_rec25 z=", round(summary(pr)$coefficients["n_rec25",3],1), "\n")
xb <- predict(pr, dS, type="link"); dS$imr <- dnorm(xb)/pnorm(xb)
for (y in c("idx_provisao","idx_estrategica")) {
  d <- cc(dS[dS$atuou_2025==1,], c(y, CAP, CTX, EVT, "imr", "uf"))
  m0 <- lm(reformulate(c(CAP,CTX,EVT,"uf"), y), data=d); m1 <- lm(reformulate(c(CAP,CTX,EVT,"imr","uf"), y), data=d)
  c1 <- coef_test(m1, vcov="CR2", cluster=d$uf, test="Satterthwaite")
  cat(y, " IMR b=", round(coef(m1)["imr"],2), " p=", round(c1$p_Satt[names(coef(m1))=="imr"],3), "\n")
  print(round(cbind(sem=coef(m0)[c("cofin","conselho","ouvid","capacit","estimativa")], com=coef(m1)[c("cofin","conselho","ouvid","capacit","estimativa")]),2))
}
## ---------- 8. Descritivos ----------
cat("\n#### DESCRITIVOS\n")
cat("prov3 média", round(mean(A$idx_provisao),1), " CV", round(sd(A$idx_provisao)/mean(A$idx_provisao),2), " | plan média", round(mean(A$idx_estrategica),1), " CV", round(sd(A$idx_estrategica)/mean(A$idx_estrategica),2), "\n")
print(round(rbind(prov=tapply(A$idx_provisao, A$porte, mean), plan=tapply(A$idx_estrategica, A$porte, mean), prep=tapply(base$idx_preparacao, base$porte, mean)),1))
print(round(rbind(prov=tapply(A$idx_provisao, A$regiao, mean), plan=tapply(A$idx_estrategica, A$regiao, mean), prep=tapply(base$idx_preparacao, base$regiao, mean)),1))
print(table(base$porte)); print(table(A$porte)); cat("cor prov4-prov3", round(cor(A$idx_provisao, A$idx_provisao4),2), "\n")
cat("prep sem conselho:", round(mean(base$idx_preparacao[is.na(base$pol_conselho_indice)]),1), " com:", round(mean(base$idx_preparacao[!is.na(base$pol_conselho_indice)]),1), "\n")
mp4 <- mlm("idx_provisao4", A, c(CAP, CTX, EVT)); cat("PROV4 sensibilidade\n"); print(round(mp4$tab[KEY,],3))
sink()
