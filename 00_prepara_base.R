# =============================================================================
# 00_prepara_base.R
# Lê a base municipal unificada (CSV), define os tipos das variáveis,
# anexa os rótulos do dicionário e salva uma versão .rds pronta para uso.
#
# Como usar: coloque este script na mesma pasta dos arquivos
#   base_suas_desastres_2025.csv e dicionario_variaveis.csv
# e defina essa pasta como diretório de trabalho (setwd) ou abra um
# projeto do RStudio nela. Depois rode o script inteiro.
# Pacotes necessários: nenhum além do R base.
# =============================================================================

# 1. Leitura --------------------------------------------------------------------
base <- read.csv("base_suas_desastres_2025.csv",
                 fileEncoding = "UTF-8",
                 na.strings = "NA", stringsAsFactors = FALSE)
dic  <- read.csv("dicionario_variaveis.csv",
                 fileEncoding = "UTF-8",
                 stringsAsFactors = FALSE)

cat("Municípios:", nrow(base), " | Variáveis:", ncol(base), "\n")

# 2. Variáveis categóricas como fatores ---------------------------------------
# Região com o Sudeste como categoria de referência (como nos modelos do artigo)
base$regiao <- factor(base$regiao,
                      levels = c("Sudeste", "Centro Oeste", "Nordeste", "Norte", "Sul"))
# Porte em ordem crescente
base$porte <- factor(base$porte,
                     levels = c("Pequeno I", "Pequeno II", "Médio", "Grande", "Metrópole"))
base$uf <- factor(base$uf)

# 3. Rótulos (ficam visíveis no View() do RStudio e em pacotes como labelled) ---
for (v in names(base)) {
  r <- dic$descricao[dic$variavel == v]
  if (length(r) == 1) attr(base[[v]], "label") <- r
}

# 4. Salvar ----------------------------------------------------------------------
saveRDS(base, "base_suas_desastres_2025.rds")
cat("Base salva em base_suas_desastres_2025.rds\n")

# Dicas rápidas -----------------------------------------------------------------
# Ver o rótulo de uma variável:   attr(base$idx_preparacao, "label")
# Procurar variáveis por tema:    dic[grepl("conselho", dic$descricao), c("variavel","descricao")]
# Variáveis de um bloco:          dic$variavel[dic$bloco == "4. Analítica"]
