
# Pacotes -----------------------------------------------------------------

# Instala apenas os pacotes que não estão instalados
install.packages(setdiff(c("geobr", "dplyr", "sp", "sf", "spatstat", "ggspatial",
                           "ggrepel"), 
                         rownames(installed.packages())), dependencies = TRUE)

library(geobr)        # Obtenção de shapefiles do IBGE
library(dplyr)        # Manipulação de dados
library(sp)           # Análise de dados espaciais
library(sf)           # Manipulação de dados espaciais
library(spatstat)     # Análise de padrões pontuais
library(ggspatial)    # Extensão do ggplot2 para dados espaciais
library(ggrepel)      # # Evita sobreposição de rótulos em gráficos ggplot2

# Leitura dos dados -------------------------------------------------------

# Carrega os dados de queimadas a partir de um arquivo CSV
df_queimadas <- read.csv("Dados/Pontual/focos_qmd_inpe_2025-01-01_2025-01-31_36.470290.csv")

# Exibe um resumo estatístico dos dados
summary(df_queimadas)

# Exibe a quantidade total de focos de queimadas no dataset
nrow(df_queimadas)

# Malha estadual ----------------------------------------------------------

# Obtém a malha estadual de Minas Gerais (ano 2020)
MGbord <- read_state(year = 2020) |> filter(name_state == "Minas Gerais")

# Transforma para o sistema de coordenadas UTM (Zona 22S)
MGbord_projetada <- st_transform(MGbord, crs = 32722)

# Converte as coordenadas para quilômetros para melhor interpretação
MGbord_projetada_KM <- st_geometry(MGbord_projetada) / 1000

# Distribuição ------------------------------------------------------------

# Seleciona apenas as colunas de Longitude e Latitude dos focos de queimadas
xy_vetor <- df_queimadas[,c("Longitude","Latitude")]
head(xy_vetor)  # Exibe as primeiras linhas dos dados de coordenadas

# Define o sistema de projeção das coordenadas geográficas (WGS84)
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Cria um objeto SpatialPointsDataFrame com os pontos de queimadas
sp.points <- SpatialPointsDataFrame(coords=xy_vetor, data=df_queimadas, proj4string=proj)

# Obtém a caixa delimitadora (bounding box) dos pontos espaciais
bbox(sp.points)

# Converte a malha estadual para uma janela de observação espacial
window <- as.owin(MGbord_projetada_KM)

# Converte os pontos espaciais para formato sf (Simple Features)
sp.points_sf <- st_as_sf(sp.points)

# Transforma os pontos para o sistema de coordenadas UTM (Zona 22S)
sp.points_projetado <- st_transform(sp.points_sf, crs = 32722)

# Converte as coordenadas dos pontos para quilômetros
sp.points_projetado <- st_geometry(sp.points_projetado) / 1000

# Extrai as coordenadas dos pontos projetados
coords <- st_coordinates(sp.points_projetado)

# Cria um objeto de padrão de pontos (PPP - Planar Point Pattern)
MG.ppp <- ppp(x=coords[,1], y=coords[,2], window=window)

# Plota os pontos sobre o mapa da região
plot(MG.ppp, pch=16, cex=0.5, main=" ")


# Intensidade Kolmogorov-Smirnov -------------------------------------------

# Estimando a intensidade média dos pontos (lambda)
lamb <- summary(MG.ppp)$intensity
lamb  # Exibe a intensidade média

# Teste de Kolmogorov-Smirnov para avaliar se os pontos seguem um processo de Poisson homogêneo
KS <- cdf.test(MG.ppp, "x", test="ks")

# Plota o resultado do teste de Kolmogorov-Smirnov
plot(KS, main= "")

# Obtém e exibe o valor-p do teste
pval <- KS$p.value
pval
KS


# Mapas de Intensidade ----------------------------------------------------

# Estima a densidade dos pontos (suavização da intensidade espacial)
density <- density(MG.ppp)

# Plota o mapa de intensidade dos focos de queimadas
plot(density, xlab="Distância", col=terrain.colors(200),
     ylab="Distância (quilômetros)", use.marks=TRUE,
     main="")

# Adiciona os pontos ao mapa de intensidade
plot(MG.ppp, add=T, cex=0.1, use.marks=TRUE)

# Adiciona contornos ao mapa de intensidade
contour(density, main="Mapa de Contornos", add=TRUE)


# Função F ----------------------------------------------------------------

# Define um vetor de distâncias para análise das funções F, G e K (até 50 km)
r <- seq(0, 50, length=10000)

# Calcula envelopes da função F
F <- envelope(MG.ppp, Fest, nsim = 500, r=r, correction="border")

# Plota a função F observada e os envelopes teóricos
plot(F$r, F$theo, main = "", type="n",
     xlab="Distância (Quilometros)", ylab="F(x)")

# Adiciona linhas representando os diferentes valores da função F
lines(F$r, F$theo, lty=3, col="blue")  # F teórica
lines(F$r, F$hi, lty=2, col="red")    # Envelope superior
lines(F$r, F$lo, lty=2, col="red")    # Envelope inferior
lines(F$r, F$obs, lty=1, col="black") # F observada

# Adiciona legenda ao gráfico
legend("bottomright", legend=c("F teórica", "F superior", "F inferior", "F observada"),
       col=c("blue", "red", "red", "black"),
       lty=c(3, 2, 2, 1), cex = 0.7)



# Função G ----------------------------------------------------------------


# Calcula envelopes da função G
G <- envelope(MG.ppp, Gest, nsim = 500, r=r, correction="border")

# Plota a função G observada e os envelopes teóricos
plot(G$r, G$theo, main = "", type="n",
     xlab="Distância (Quilometros)", ylab="G(y)")

# Adiciona linhas representando os diferentes valores da função G
lines(G$r, G$theo, lty=3, col="blue")  # G teórica
lines(G$r, G$hi, lty=2, col="red")    # Envelope superior
lines(G$r, G$lo, lty=2, col="red")    # Envelope inferior
lines(G$r, G$obs, lty=1, col="black") # G observada

# Adiciona legenda ao gráfico
legend("bottomright", legend=c("G teórica", "G superior", "G inferior", "G observada"),
       col=c("blue", "red", "red", "black"),
       lty=c(3, 2, 2, 1), cex = 0.6)


# Função K ----------------------------------------------------------------

# Calcula envelopes da função K de Ripley para análise de agrupamento espacial
L <- envelope(MG.ppp, Lest, nsim = 500, r=r, correction="border")

# Plota a função K observada e os envelopes teóricos
plot(L$r, L$theo, main = "", type="n",
     xlab="Distância (Quilometros)", ylab="K(d)")

# Adiciona linhas representando os diferentes valores da função K
lines(L$r, L$theo, lty=3, col="blue")  # K teórica
lines(L$r, L$hi, lty=2, col="red")    # Envelope superior
lines(L$r, L$lo, lty=2, col="red")    # Envelope inferior
lines(L$r, L$obs, lty=1, col="black") # K observada

# Adiciona legenda ao gráfico
legend("bottomright", legend=c("K teórica", "K superior", "K inferior", "K observada"),
       col=c("blue", "red", "red", "black"),
       lty=c(3, 2, 2, 1), cex = 0.7)


# Municipios com mais focos -----------------------------------------------

# Fazer cada linha do banco de dados vale 1
df_queimadas$foco <- 1

# Selecionar os 5 municipios com mais focos
focos <- df_queimadas |> 
  group_by(Municipio) |> 
  summarise(focos = sum(foco), 
            lon = mean(Longitude),
            lat = mean(Latitude)) |> 
  arrange(desc(focos)) |> 
  slice_max(focos, n = 6)

head(focos)
# Criação do gráfico
ggplot(data = MGbord) +
  geom_sf(fill = "white", color = "gray", size = 0.2) +  # Mapa com bordas cinza escuro
  geom_point(data = focos, aes(x = lon, y = lat), color = "red", size = 4) +  # Pontos dos focos
  geom_text_repel(data = focos, aes(x = lon, y = lat, label = paste(Municipio, "\n", focos, "focos")), 
                  size = 2.5, nudge_y = 0.5) +  # Rótulos dos focos
  theme_minimal() +  # Tema minimalista
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.5),  # Grade do fundo
    panel.background = element_rect(fill = "white"),  # Fundo branco
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"  # Remover a legenda
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) +  # Escala gráfica
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering) +  # Seta de norte
  labs(
    title = "",  # Título do mapa
    x = "Longitude",  # Título do eixo x
    y = "Latitude"  # Título do eixo y
  )

####################################################################################
####################################################################################
####################################################################################

# Leitura dos dados -------------------------------------------------------

# Carrega os dados de acidentes a partir de um arquivo CSV
df_acidentes <- read.csv("Dados/Pontual/acidentes2024_todas_causas_tipos/Acidentes_MG_2024.csv")

# Exibe um resumo estatístico dos dados
summary(df_acidentes)

# corrigindo as coordenadas
df_acidentes$longitude[df_acidentes$longitude < -60] <- df_acidentes$longitude[df_acidentes$longitude < -60] / 1000
df_acidentes$latitude[df_acidentes$latitude < -30] <- df_acidentes$latitude[df_acidentes$latitude < -30] / 1000

summary(df_acidentes)

# Exibe a quantidade total de acidentes no dataset
nrow(df_acidentes)

# Malha estadual ----------------------------------------------------------

# Obtém a malha estadual de Minas Gerais (ano 2020)
MGbord <- read_state(year = 2020) |> filter(name_state == "Minas Gerais")

# Transforma para o sistema de coordenadas UTM (Zona 22S)
MGbord_projetada <- st_transform(MGbord, crs = 32722)

# Converte as coordenadas para quilômetros para melhor interpretação
MGbord_projetada_KM <- st_geometry(MGbord_projetada) / 1000

# Distribuição ------------------------------------------------------------

# Seleciona apenas as colunas de Longitude e Latitude dos focos de queimadas
xy_vetor <- df_acidentes[,c("longitude","latitude")]
head(xy_vetor)  # Exibe as primeiras linhas dos dados de coordenadas

# Define o sistema de projeção das coordenadas geográficas (WGS84)
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Cria um objeto SpatialPointsDataFrame com os pontos de queimadas
sp.points <- SpatialPointsDataFrame(coords=xy_vetor, data=df_acidentes, proj4string=proj)

# Obtém a caixa delimitadora (bounding box) dos pontos espaciais
bbox(sp.points)

# Converte a malha estadual para uma janela de observação espacial
window <- as.owin(MGbord_projetada_KM)

# Converte os pontos espaciais para formato sf (Simple Features)
sp.points_sf <- st_as_sf(sp.points)

# Transforma os pontos para o sistema de coordenadas UTM (Zona 22S)
sp.points_projetado <- st_transform(sp.points_sf, crs = 32722)

# Converte as coordenadas dos pontos para quilômetros
sp.points_projetado <- st_geometry(sp.points_projetado) / 1000

# Extrai as coordenadas dos pontos projetados
coords <- st_coordinates(sp.points_projetado)

# Cria um objeto de padrão de pontos (PPP - Planar Point Pattern)
MG.ppp <- ppp(x=coords[,1], y=coords[,2], window=window)

# Plota os pontos sobre o mapa da região
plot(MG.ppp, pch=16, cex=0.5, main=" ")


# Intensidade Kolmogorov-Smirnov -------------------------------------------

# Estimando a intensidade média dos pontos (lambda)
lamb <- summary(MG.ppp)$intensity
lamb  # Exibe a intensidade média

# Teste de Kolmogorov-Smirnov para avaliar se os pontos seguem um processo de Poisson homogêneo
KS <- cdf.test(MG.ppp, "x", test="ks")

# Plota o resultado do teste de Kolmogorov-Smirnov
plot(KS, main= "")

# Obtém e exibe o valor-p do teste
pval <- KS$p.value
pval
KS


# Mapas de Intensidade ----------------------------------------------------

# Estima a densidade dos pontos (suavização da intensidade espacial)
density <- density(MG.ppp)

# Plota o mapa de intensidade dos acidentes
plot(density, xlab="Distância", col=terrain.colors(200),
     ylab="Distância (quilômetros)", use.marks=TRUE,
     main="")

# Adiciona os pontos ao mapa de intensidade
plot(MG.ppp, add=T, cex=0.1, use.marks=TRUE)

# Adiciona contornos ao mapa de intensidade
contour(density, main="Mapa de Contornos", add=TRUE)


# Função F ----------------------------------------------------------------

# Define um vetor de distâncias para análise das funções F, G e K (até 5 km)
r <- seq(0, 5, length=10000)

# Calcula envelopes da função F
F <- envelope(MG.ppp, Fest, nsim = 5, r=r, correction="border")

# Plota a função F observada e os envelopes teóricos
plot(F$r, F$theo, main = "", type="n",
     xlab="Distância (Quilometros)", ylab="F(x)")

# Adiciona linhas representando os diferentes valores da função F
lines(F$r, F$theo, lty=3, col="blue")  # F teórica
lines(F$r, F$hi, lty=2, col="red")    # Envelope superior
lines(F$r, F$lo, lty=2, col="red")    # Envelope inferior
lines(F$r, F$obs, lty=1, col="black") # F observada

# Adiciona legenda ao gráfico
legend("right", legend=c("F teórica", "F superior", "F inferior", "F observada"),
       col=c("blue", "red", "red", "black"),
       lty=c(3, 2, 2, 1), cex = 0.8)



# Função G ----------------------------------------------------------------


# Calcula envelopes da função G
G <- envelope(MG.ppp, Gest, nsim = 5, r=r, correction="border")


# Plota a função G observada e os envelopes teóricos
plot(G$r, G$theo, main = "", type="n",
     xlab="Distância (Quilometros)", ylab="G(y)")

# Adiciona linhas representando os diferentes valores da função G
lines(G$r, G$theo, lty=3, col="blue")  # G teórica
lines(G$r, G$hi, lty=2, col="red")    # Envelope superior
lines(G$r, G$lo, lty=2, col="red")    # Envelope inferior
lines(G$r, G$obs, lty=1, col="black") # G observada

# Adiciona legenda ao gráfico
legend("bottomright", legend=c("G teórica", "G superior", "G inferior", "G observada"),
       col=c("blue", "red", "red", "black"),
       lty=c(3, 2, 2, 1), cex = 0.8)


# Função K ----------------------------------------------------------------

# Calcula envelopes da função K de Ripley para análise de agrupamento espacial
L <- envelope(MG.ppp, Lest, nsim = 5, r=r, correction="border")

# Plota a função K observada e os envelopes teóricos
plot(L$r, L$theo, main = "", type="n",
     xlab="Distância (Quilometros)", ylab="K(d)",ylim = c(0,35))

# Adiciona linhas representando os diferentes valores da função K
lines(L$r, L$theo, lty=3, col="blue")  # K teórica
lines(L$r, L$hi, lty=2, col="red")    # Envelope superior
lines(L$r, L$lo, lty=2, col="red")    # Envelope inferior
lines(L$r, L$obs, lty=1, col="black") # K observada

# Adiciona legenda ao gráfico
legend("right", legend=c("K teórica", "K superior", "K inferior", "K observada"),
       col=c("blue", "red", "red", "black"),
       lty=c(3, 2, 2, 1), cex = 0.8)


# Municipios com mais acidentes -----------------------------------------------

# Fazer cada linha do banco de dados vale 1
df_acidentes$acidente <- 1

# Selecionar os 5 municipios com mais acidentes
acidentes <- df_acidentes |> 
  group_by(municipio) |> 
  summarise(acidentes = sum(acidente), 
            lon = mean(longitude),
            lat = mean(latitude)) |> 
  arrange(desc(acidentes)) |> 
  slice_max(acidentes, n = 5)

head(acidentes)
# Criação do gráfico
ggplot(data = MGbord) +
  geom_sf(fill = "white", color = "gray", size = 0.2) +  # Mapa com bordas cinza escuro
  geom_point(data = acidentes, aes(x = lon, y = lat), color = "red", size = 4) +  # Pontos dos focos
  geom_text_repel(data = acidentes, aes(x = lon, y = lat, label = paste(municipio, "\n", acidentes, "acidentes")), 
                  size = 2.5, nudge_y = 0.5) +  # Rótulos dos acidentes
  theme_minimal() +  # Tema minimalista
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.5),  # Grade do fundo
    panel.background = element_rect(fill = "white"),  # Fundo branco
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"  # Remover a legenda
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) +  # Escala gráfica
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering) +  # Seta de norte
  labs(
    title = "",  # Título do mapa
    x = "Longitude",  # Título do eixo x
    y = "Latitude"  # Título do eixo y
  )
