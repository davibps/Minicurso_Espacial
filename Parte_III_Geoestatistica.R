
# Pacotes -----------------------------------------------------------------

library(geoR)         # Análise geoestatística
library(geobr)        # Obtenção de shapefiles do IBGE
library(dplyr)        # Manipulação de dados
library(MASS)         # Função boxcox (transformação de Box-Cox)
library(sf)           # Manipulação de dados espaciais
library(ggplot2)      # Criação de gráficos

# Leitura dos dados -------------------------------------------------------

# Carrega os dados de precipitação a partir de um arquivo CSV
df_precipitacao <- read.csv("Dados/Geoestatistica/focos_qmd_inpe_2025-01-01_2025-01-31_36.470290.csv")

# Converte os dados para o formato 'geodata' do pacote geoR,
# especificando as colunas de coordenadas (colunas 11 e 10) e a coluna de dados (coluna 8)
tb1 <- as.geodata(df_precipitacao, coords.col = c(11,10), data.col = 8)
summary(tb1)        # Exibe um resumo dos dados
plot(tb1, low = TRUE)  # Plota os dados com a opção 'low' para visualização dos pontos de baixa densidade

# Ajuste de dados para evitar problemas com logaritmo zero
tb1$data <- tb1$data + 0.01  
boxcox(tb1)         # Aplica transformação Box-Cox (para verificar o melhor parâmetro de transformação)

# Transformação logarítmica dos dados de precipitação
df_precipitacao$Precipitacao <- log10(df_precipitacao$Precipitacao + 0.01)
# Reconverte os dados transformados para o objeto 'geodata'
tb1 <- as.geodata(df_precipitacao, coords.col = c(11,10), data.col = 8)
summary(tb1)        # Exibe resumo após a transformação

# Plota os dados novamente, com diferentes tendências
plot(tb1, low = TRUE)
plot(tb1, low = TRUE, trend = "1st")  # Considera tendência linear (1ª ordem)
plot(tb1, low = TRUE, trend = "2nd")  # Considera tendência quadrática (2ª ordem)


# Definição dos limites (borders) para a área de previsão -----------------

# Lê o shapefile do estado de Minas Gerais (ano 2020)
bor <- read_state(year = 2020) |> filter(name_state == "Minas Gerais")
# Transforma o polígono em um data frame com coordenadas (longitude e latitude)
bor_coords <- st_coordinates(bor)       # Extrai as coordenadas do polígono
bor_df <- as.data.frame(bor_coords[, 1:2]) # Seleciona apenas X (Longitude) e Y (Latitude)
colnames(bor_df) <- c("Longitude", "Latitude")

# Adiciona os limites (borders) ao objeto 'tb1'
tb1$borders <- with(bor_df, cbind(Longitude, Latitude))


# Cálculo do variograma ---------------------------------------------------


# Calcula o variograma empírico com tendência de 1ª ordem e distância máxima de 5
grama <- variog(tb1, max.dist = 5, trend = "1st")
# Gera envelope (simulação) do variograma para avaliação de incertezas com 100 simulações
env.var <- variog.mc.env(tb1, obj.v = grama, nsim = 100)
ef1 = eyefit(grama)  # Código comentado para ajuste visual do variograma
summary(ef1)
# Plota o variograma empírico juntamente com o envelope de simulação
plot(grama, env = env.var)


# Ajuste de modelos de variograma usando máxima verossimilhança -----------

# Valores iniciais e nugget fixo (sigma2, phi e tau2)
# sigma2 = 1.31; phi = 1.6; tau2 = 0.16

# Cria uma lista para armazenar os ajustes de diferentes modelos
tb1.ml0 <- list()
# Ajusta modelo com covariância 'matern'
tb1.ml0$fit1 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "matern", trend = "1st")
# Ajusta modelo com covariância 'exponential'
tb1.ml0$fit2 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "exponential", trend = "1st")
# Ajusta modelo com covariância 'gaussian'
tb1.ml0$fit3 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "gaussian", trend = "1st")
# Ajusta modelo com covariância 'spherical'
tb1.ml0$fit4 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "spherical", trend = "1st")
# Ajusta modelo com covariância 'circular'
tb1.ml0$fit5 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "circular", trend = "1st")
# Ajusta modelo com covariância 'cubic'
tb1.ml0$fit6 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "cubic", trend = "1st")
# Ajusta modelo com covariância 'wave'
tb1.ml0$fit7 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "wave", trend = "1st")
#tb1.ml0$fit8 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "power", trend = "1st") # Código comentado
# Ajusta modelo com covariância 'powered.exponential'
tb1.ml0$fit9 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "powered.exponential", trend = "1st")
# Ajusta modelo com covariância 'cauchy'
tb1.ml0$fit10 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "cauchy", trend = "1st")
#tb1.ml0$fit11 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "gencauchy", trend = "1st") # Código comentado
# Ajusta modelo com covariância 'gneiting'
tb1.ml0$fit12 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "gneiting", trend = "1st")
#tb1.ml0$fit13 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "gneiting.matern", trend = "1st") # Código comentado
# Ajusta modelo com covariância 'pure.nugget' (apenas nugget, sem estrutura espacial)
tb1.ml0$fit14 <- likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "pure.nugget", trend = "1st")

# Compara os modelos utilizando o critério de informação Akaike (AIC)
sapply(tb1.ml0, AIC)

# Exibe o resumo do modelo com covariância 'cauchy'
cauchymodel <- summary(likfit(tb1, ini = c(1.31, 1.6), nug = 0.16, cov.model = "cauchy", trend = "1st"))
cauchymodel

# Plota os pontos dos dados com divisão em quintis e legendas
points(tb1, pt.div = "quintile", xlab = "Coord X", 
       ylab = "Coord Y", x.leg = -15.8, y.leg = -48, trend = "1st")

#colocar legenda

# Cria uma grade de previsão dentro dos limites definidos -----------------

gr <- pred_grid(tb1$borders, by = .1)   # Cria uma grade com espaçamento de 0.1
gr0 <- gr[.geoR_inout(gr, tb1$borders), ] # Seleciona os pontos que estão dentro dos limites


# Krigagem ----------------------------------------------------------------

# Define o controle de krigagem usando o modelo 'cauchy' ajustado anteriormente
KC <- krige.control(type = "SK", obj.model = tb1.ml0$fit10)
# Realiza a krigagem simples nos pontos da grade
tb1.kc <- krige.conv(tb1, loc = gr0, krige = KC)


# Mapeia os resultados da krigagem ----------------------------------------

# Imagem do valor interpolado (utilizando o objeto 'tb1.kc')
image(tb1.kc, gr,
      tb1$borders,
      col = terrain.colors(200),
      x.leg = c(-45, -40.8), y.leg = c(-24, -23.5), cex = 0.7, trend = "1st")

# Imagem da variância da krigagem (incerteza)
image(tb1.kc, loc = gr,
      bor = tb1$borders, col = terrain.colors(21),
      val = sqrt(tb1.kc$krige.var),
      x.leg = c(-50.5, -40.8), y.leg = c(-24, -23.5), cex = 0.6)


# Saída condicional da krigagem -------------------------------------------


# Define controles de saída com limiar e quantis para análise posterior
OC <- output.control(thres = median(tb1$data),
                     quantile = c(0.1, 0.9))

# Realiza a krigagem com controle de saída
tb1.kc <- krige.conv(tb1, loc = gr0, krige = KC, out = OC)

# Imagem da probabilidade condicional (1 - probabilidade)
image(tb1.kc, loc = gr, bor = tb1$borders, col = terrain.colors(21), 
      val = 1 - tb1.kc$prob, x.leg = c(-50.5, -40.8), y.leg = c(-24, -23.5), cex = 0.7)

# Imagem do quantil de 90% da krigagem
image(tb1.kc, loc = gr, bor = tb1$borders, col = terrain.colors(21), 
      val = tb1.kc$quantile$q90, x.leg = c(-50.5, -40.8), y.leg = c(-24, -23.5), cex = 0.7)

# Imagem do quantil de 10% da krigagem
image(tb1.kc, loc = gr, bor = tb1$borders, col = terrain.colors(21), 
      val = tb1.kc$quantile$q10, x.leg = c(-50.5, -40.8), y.leg = c(-24, -23.5), cex = 0.7)


# Simulações  -------------------------------------------------------------

# Calcula a probabilidade de que as simulações ultrapassem a mediana dos dados
p091 <- apply(tb1.kc$simul, 2, function(x) mean(x > median(tb1$data)))
# Plota um histograma das probabilidades obtidas com a densidade sobreposta
hist(p091, prob = TRUE)
lines(density(p091))
