install.packages("vtree")

library(formattable)
library(dplyr)
library(tidyr)
library(data.table)
library(janitor)
library(vtree)

customGreen0 = "#DeF7E9"

customGreen = "#71CA97"

customRed = "#ff7f7f"

deepsky = "#00BFFF"
hotpink = "#FF69B4"

skyblue = "#87CEEB"
lightcoral = "#F08080"
tabh2<- tabyl(tab, SingleR.labels, seurat_clusters)

formattable(tabh2)


formattable(tabh2, 
            align =c("l","c","c","c","c", "c", "c", "c", "c", "c", "c", "c"), 
            list(`seurat_clusters` = formatter(
              "span", style = ~ style(color = "grey", font.weight = "bold")),
              `0`= color_tile(customGreen, customGreen0),
              `Activated 2`= color_tile(customGreen, customGreen0),
              `Basal 1`= color_tile(customGreen, customGreen0),
              `Basal 2`= color_tile(skyblue, lightcoral),
              `CD21low 1`= color_tile(customGreen, customGreen0),
              `CD21low 2`= color_tile(customGreen, customGreen0),
              `FCRL2/3high MBC`= color_tile(customGreen, customGreen0),
              `FCRL4+`= color_tile(customGreen, customGreen0),
              `HSP-response`= color_tile(customGreen, customGreen0),
              `IFN-response`= color_tile(customGreen, customGreen0),
              `preGC MBC`= color_tile(customGreen, customGreen0)
              
            ))
         

finaltab <-formattable(tabh, 
                 align =c("l","c","c","c","c", "c", "c", "c", "c", "c", "c", "c"), 
                 list(`seurat_clusters` = formatter(
                    "span", style = ~ style(color = "grey",font.weight = "bold")),
                        `Activated 1`= color_tile(skyblue, lightcoral),
                        `Activated 2`= color_tile(skyblue, lightcoral),
                         `Basal 1`= color_tile(skyblue, lightcoral),
                         `Basal 2`= color_tile(skyblue, lightcoral),
                         `CD21low 1`= color_tile(skyblue, lightcoral),
                         `CD21low 2`= color_tile(skyblue, lightcoral),
                         `FCRL2/3high MBC`= color_tile(skyblue, lightcoral),
                         `FCRL4+`= color_tile(skyblue, lightcoral),
                         `HSP-response`= color_tile(skyblue, lightcoral),
                         `IFN-response`= color_tile(skyblue, lightcoral),
                         `preGC MBC`= color_tile(skyblue, lightcoral)
                  
                          ))   

finaltab
?formattable

