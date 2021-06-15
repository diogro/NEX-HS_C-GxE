if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(superheat)){install.packages("superheat"); library(superheat)}
if(!require(mcclust)){install.packages("mcclust"); library(mcclust)}

head = read_csv("data/output/6128genes-non_hierarchical_sbm-fit_head_df.csv")
body = read_csv("data/output/6176genes-non_hierarchical_sbm-fit_body_df.csv")

head$tissue = "head"
body$tissue = "body"

head = dplyr::semi_join(head, body, by = "Gene")
body = dplyr::semi_join(body, head, by = "Gene")
head = dplyr::semi_join(head, body, by = "Gene")
body = dplyr::semi_join(body, head, by = "Gene")

table(head$Block) |> length()
table(body$Block) |> length()

body_expr = t(read.table("./data/VOOMCounts_CPM1_body_ctrl_Jun2.20.txt",row.names = 1))
body_expr = body_expr[,body$Gene]
body_expr = as_tibble(body_expr)

head_expr = t(read.table("./data/VOOMCounts_CPM1_head_ctrl_Jun2.20.txt",row.names = 1))
head_expr = head_expr[,head$Gene]
head_expr = as_tibble(head_expr)

cor_head = cor(head_expr)
cor_body = cor(body_expr)

superheat(cor_head, 
          membership.rows = head$Block, 
          membership.cols = head$Block, 
          smooth.heat = TRUE)
superheat(cor_body, 
          membership.rows = body$Block, 
          membership.cols = body$Block, 
          smooth.heat = TRUE)

vi.dist(head$Block, body$Block, parts = TRUE)
