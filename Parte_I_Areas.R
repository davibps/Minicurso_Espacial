
# Pacotes -----------------------------------------------------------------

# Instala apenas os pacotes que não estão instalados
install.packages(setdiff(c("tidyverse", "dplyr", "geobr", "sf", "readxl", 
                           "stringr", "ggspatial", "spdep", "sp"), 
                         rownames(installed.packages())), dependencies = TRUE)

# Carrega os pacotes
library(tidyverse) # Manipulação de dados
library(dplyr) # Manipulação de dados
library(geobr) # Shapefiles IBGE
library(sf) # Manipulação de dados espaciais
library(readxl) # Leitura de dados formato Excel
library(stringr) # Manipulação de dados qualitativos
library(ggspatial) # Extensão do ggplot2 para dados espaciais
library(spdep) #Manipulação de dados espaciais
library(sp) #Manipulação de dados espaciais


# Polígonos ---------------------------------------------------------------

# Obter a malha / shapefile do IBGE (Estados do Brasil)
shp <- geobr::read_municipality(year = 2022) # último ano disponível: 2022

# Visualização gráfica - Brasil
shp |> ggplot() + geom_sf()

# Obter a malha do Estado de Minas Gerais
shp_MG <- shp |> filter(name_state=="Minas Gerais")

# Visualização gráfica - Minas Gerais
shp_MG |> ggplot() + geom_sf()

# Salvar a malha na pasta do projeto
shp_MG |> st_write("Malhas/MG_Estado_2022.shp",delete_layer = TRUE) 

# Adicional: Bairros
read_neighborhood(year = 2022) |> 
  filter(name_muni=="Belo Horizonte") |> 
  ggplot() + geom_sf()

# Também é possível obter as malhas no site do IBGE:
# https://www.ibge.gov.br/geociencias/organizacao-do-territorio/malhas-territoriais/15774-malhas.html


# Leitura dos dados -------------------------------------------------------

df_dengue_MG <- read_excel("Dados/Areas/Dados_Dengue_Pop.xlsx")
head(df_dengue_MG) # Verificar estrutura da base
summary(df_dengue_MG) # Resumo descritivo dos dados

# Preparação dos dados Espaciais ------------------------------------------

ind <- match(shp_MG$name_muni, df_dengue_MG$Município); sum(is.na(ind))
shp_MG$name_muni[is.na(ind)] # Verificar os municípios divergentes

# Correção na base de dados df_dengue_MG
df_dengue_MG$Município <- str_replace(df_dengue_MG$Município, "Barão do Monte Alto", "Barão de Monte Alto")
df_dengue_MG$Município <- str_replace(df_dengue_MG$Município, "Olhos-d'Água", "Olhos-D'água")
df_dengue_MG$Município <- str_replace(df_dengue_MG$Município, "Pingo-d'Água", "Pingo-D'água")
df_dengue_MG$Município <- str_replace(df_dengue_MG$Município, "São João del Rei", "São João Del Rei")

ind <- match(shp_MG$name_muni, df_dengue_MG$Município); sum(is.na(ind))

# Reoordendando e ajustando a base
df_dengue_ind <- df_dengue_MG[ind,]
row.names(df_dengue_ind) <- shp_MG$name_muni

# Unir os dataframes
ShapeMG <- left_join(shp_MG, df_dengue_ind, by = c("name_muni" = "Município"))
head(ShapeMG)


# Mapa de Quantis ------------------------------------------------------

# Definir os intervalos
quantile(ShapeMG$Casos_mil_hab, seq(0.2,1,by=0.2)) # Calcula os quantis
intervalos <- c(-1, 24.424, 43.564, 70.510, 114.416, 319.930) # colocar os intervalos (no primeiro subtrair 1 no último adicionar 1)
legendas <- c("Entre 0 e 24,424", "Entre 24,424 e 43,564", "Entre 43,564 e 70,510", 
               "Entre 70,510 e 114,416","Entre 114,416 e 318,930")
colors <- c("Entre 0 e 24,424"="white", "Entre 24,424 e 43,564"="#FFC600", "Entre 43,564 e 70,510"="#FF8D00", 
             "Entre 70,510 e 114,416"="#FF3800", "Entre 114,416 e 318,930"="darkred")


# Adicionar a coluna de classe ao ShapeMG
ShapeMG$quantis <- cut(ShapeMG$Casos_mil_hab, breaks = intervalos, labels = legendas)

# Mapa
ggplot(data = ShapeMG) +
  geom_sf(aes(fill = quantis), color = "darkgray", size = 0.2) +
  scale_fill_manual(values = colors, name = "Taxa 1:1000") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "white", color = "black"),
    axis.title = element_blank()  # Remove títulos dos eixos
  ) +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering)


# Vizinhança --------------------------------------------------------------

ShapeMG.nb <- poly2nb(ShapeMG, queen = T)
vizinhanca <- nb2listw(ShapeMG.nb, style="C"); vizinhanca # Matriz de pesos simétrica

plot(ShapeMG$geom, col = "lightgray")
plot(ShapeMG.nb, ShapeMG$geom, add = TRUE, col = "blue", lwd = 0.5)


# Índice de Moran Global --------------------------------------------------

Mglobal <- moran.test(ShapeMG$Casos_mil_hab, listw=vizinhanca)
Mglobal


# Índice de Moran Local ---------------------------------------------------

ShapeMG.mloc <- localmoran(ShapeMG$Casos_mil_hab, listw=vizinhanca,
                            zero.policy=T, 
                            alternative = "two.sided")

head(ShapeMG.mloc)


# LISA Cluster ------------------------------------------------------------

# Calcular desvios
Sd_1 <- ShapeMG$Casos_mil_hab - mean(ShapeMG$Casos_mil_hab)
mI_1 <- ShapeMG.mloc[, 1]
# Determinar os quadrantes
quadrant <- vector(mode = "numeric", length = nrow(ShapeMG))
quadrant[Sd_1 >= 0 & mI_1 >= 0] <- 1
quadrant[Sd_1 <= 0 & mI_1 >= 0] <- 2
quadrant[Sd_1 >= 0 & mI_1 <= 0] <- 3
quadrant[Sd_1 <= 0 & mI_1 <= 0] <- 4
# Significância
signif <- 0.05
quadrant[ShapeMG.mloc[, 5] > signif] <- 5
# Adicionar quadrantes ao dataframe
ShapeMG$quadrant <- factor(quadrant, levels = 1:5, 
                            labels = c("alto-alto", "baixo-baixo", "alto-baixo", 
                                       "baixo-alto", "N.Sgf"))
# Definir as cores
colors <- c("alto-alto" = "red", "baixo-baixo" = "blue", "alto-baixo" = "lightpink", 
            "baixo-alto" = "skyblue2", "N.Sgf" = "white")
# Criar o gráfico com ggplot2
ggplot(data = ShapeMG) +
  geom_sf(aes(fill = quadrant), color = "gray", size = 0.2) +
  scale_fill_manual(values = colors, name = "Quadrantes de Moran") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "white", color = "black"),
    axis.title = element_blank()  # Remove títulos dos eixos
  ) +
  # Adicionar a camada de legendas manualmente
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering)

# LISA Significância ------------------------------------------------------

# Definir os intervalos de significância LISA e cores
intervalos <- c(0, 0.001, 0.01, 0.05, 1)
legendas <- c("0.1%", "1.0%", "5.0%", "N.Sgf")
colors <- c("0.1%"="red", "1.0%"="blue", "5.0%"="skyblue", "N.Sgf"="white")

# Adicionar a coluna de classe de significância ao ShapeMG
ShapeMG$signif_classe <- cut(ShapeMG.mloc[,5], breaks = intervalos, labels = legendas, right = TRUE)

# Criar o gráfico com ggplot2
ggplot(data = ShapeMG) +
  geom_sf(aes(fill = signif_classe), color = "darkgray", size = 0.2) +
  scale_fill_manual(values = colors, name = "Significância") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "white", color = "black"),
    axis.title = element_blank()  # Remove títulos dos eixos
  ) +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering)
