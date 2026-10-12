# =============================================================================
# 01_reproduz_modelos.R  (versão final do artigo, v6)
# Reproduz as principais tabelas do artigo a partir da base unificada:
#   Tabela 3  - preparação e repertórios por porte e região
#   Tabela 4  - logit da não atuação entre municípios com desastre reconhecido em 2025
#   Tabela 5  - modelos multinível (preparação, provisão, planejamento e diferença)
#
# Rode antes o 00_prepara_base.R (gera o arquivo .rds).
# Pacotes (instale uma vez):
#   install.packages(c("lme4", "lmerTest", "clubSandwich"))
# =============================================================================

library(lme4)
library(lmerTest)       # graus de liberdade de Satterthwaite nos modelos multinível
library(clubSandwich)   # erros-padrão agrupados CR2 (correção para poucos grupos)

base <- readRDS("base_suas_desastres_2025.rds")

# -----------------------------------------------------------------------------
# 1. Variáveis dos modelos
# Contínuas são padronizadas (média 0, DP 1) na base inteira, antes de qualquer
# recorte. Os coeficientes delas são, portanto, "pontos por desvio-padrão".
# A exposição usa apenas desastres RECONHECIDOS pelo governo federal.
# -----------------------------------------------------------------------------
z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

base <- within(base, {
  # nível individual
  capacit    <- op_capacitacao_2025
  dirig_sup  <- op_dirigente_superior
  # analítica, organizacional
  vigil      <- an_vigilancia_formal
  diag       <- an_diagnostico
  sistema    <- an_sistema_informacao
  # operacional, organizacional
  secret     <- op_secretaria_exclusiva
  rh10       <- z(log1p(rh_por_10mil_hab))
  estat      <- z(rh_prop_estatutarios)
  gasto      <- z(log1p(op_gasto_as_pc))
  pmas       <- op_pmas
  consorcio  <- op_consorcio
  # política e relacional
  conselho   <- z(pol_conselho_indice)
  ouvid      <- pol_ouvidoria
  cib        <- pol_acompanha_cib
  # nível sistêmico (recursos vindos de outros entes)
  estimativa <- an_info_federal_risco
  cofin      <- sis_cofin_estadual_calamidade
  # contexto institucional
  compdec    <- dc_compdec_2020
  dc_mapeam  <- dc_mapeamento_2020
  dc_anos    <- dc_anos_gasto            # em anos, não padronizada
  hab        <- hab_prioriza_risco_2024
  alin       <- pol_alinhamento_governador
  troca      <- pol_troca_prefeito
  # contexto de risco e social
  hid        <- z(log1p(atl_rec_nhid_2020_2024))
  seca       <- z(log1p(atl_rec_nseca_2020_2024))
  ameaca     <- z(risco_ameaca_max)
  estresse   <- z(risco_estresse_hidr)
  lista      <- cc_lista
  vuln       <- z(soc_vulnerab_renda)
  logpop     <- z(log(populacao_2022))
  # diferença entre repertórios (testa igualdade dos coeficientes)
  dif        <- idx_provisao - idx_estrategica
})

CAP <- c("vigil","diag","sistema","estimativa","capacit","dirig_sup","secret","rh10","estat","gasto",
         "pmas","consorcio","cofin","conselho","ouvid","cib","compdec","dc_mapeam","dc_anos","hab","alin","troca")
CTX <- c("hid","seca","ameaca","estresse","lista","vuln","logpop")
EVT <- c("atuou_regular","evento_enchente_inundacao","evento_seca_estiagem")   # só para quem atuou
cc  <- function(d, v) { d <- d[complete.cases(d[, v]), ]; d$uf <- droplevels(d$uf); d }

# -----------------------------------------------------------------------------
# 2. Tabela 3: médias por porte e por região
# -----------------------------------------------------------------------------
A <- subset(base, atuou_2025 == 1)
for (g in c("porte", "regiao")) {
  cat("\nTabela 3 por", g, "\n")
  print(round(rbind(Preparacao  = tapply(base$idx_preparacao, base[[g]], mean),
                    Provisao    = tapply(A$idx_provisao, A[[g]], mean),
                    Planejamento= tapply(A$idx_estrategica, A[[g]], mean)), 1))
}

# -----------------------------------------------------------------------------
# 3. Tabela 4: não atuação entre municípios com desastre reconhecido em 2025
# -----------------------------------------------------------------------------
L <- subset(base, rec25_algum == 1)
L$nao   <- 1 - L$atuou_2025
L$desab <- as.numeric(L$rec25_desabrigados > 0)
L$q4    <- as.numeric(L$rec25_mes_primeiro >= 10)   # 1º reconhecimento no último trimestre
XL <- c("rec25_so_seca", CAP, CTX, "q4")
L <- cc(L, c("nao", XL, "desab", "regiao", "uf"))

tab_logit <- function(X) {
  m  <- glm(reformulate(c(X, "regiao"), "nao"), data = L, family = binomial)
  ct <- coef_test(m, vcov = "CR2", cluster = L$uf, test = "Satterthwaite")
  ci <- conf_int(m, vcov = "CR2", cluster = L$uf, test = "Satterthwaite", level = 0.95)
  # intervalos com os mesmos graus de liberdade de Satterthwaite usados nos p-valores
  data.frame(OR = round(exp(coef(m)), 2), IC_inf = round(exp(ci$CI_L), 2),
             IC_sup = round(exp(ci$CI_U), 2), p = round(ct$p_Satt, 3))
}
cat("\nTabela 4 (n =", nrow(L), ")\nModelo 1\n"); print(tab_logit(XL)[c("rec25_so_seca","diag","capacit","cofin"), ])
cat("Modelo 2 (com desabrigados)\n");          print(tab_logit(c(XL, "desab"))[c("rec25_so_seca","desab"), ])

# -----------------------------------------------------------------------------
# 4. Tabela 5: modelos multinível com intercepto aleatório por UF
# -----------------------------------------------------------------------------
roda <- function(y, dados, X) {
  d <- cc(dados, c(y, X, "regiao", "uf"))
  m <- lmer(reformulate(c(X, "regiao", "(1 | uf)"), y), data = d, REML = TRUE)
  s <- summary(m)$coefficients     # p-valores com Satterthwaite (lmerTest)
  cat(sprintf("\n%s: n = %d\n", y, nrow(d)))
  data.frame(coef = round(s[, 1], 1), ep = round(s[, 2], 1), p = round(s[, 5], 3))
}
t_prep <- roda("idx_preparacao",  base, c(CAP, CTX))
t_prov <- roda("idx_provisao",    A,    c(CAP, CTX, EVT))
t_plan <- roda("idx_estrategica", A,    c(CAP, CTX, EVT))
t_dif  <- roda("dif",             A,    c(CAP, CTX, EVT))
chave <- c("cofin","estimativa","conselho","ouvid","capacit","diag","vigil","cib","gasto","estat","hid","seca")
print(data.frame(Preparacao = t_prep[chave, "coef"], Provisao = t_prov[chave, "coef"],
            Planejamento = t_plan[chave, "coef"], Diferenca = t_dif[chave, "coef"],
            row.names = chave))
# Para ver um modelo completo: t_plan
